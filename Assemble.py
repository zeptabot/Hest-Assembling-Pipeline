#!/usr/bin/env python3
import csv, glob, io, json, os, re, shutil, sys, tarfile, tempfile, zipfile

import numpy as np
import tifffile
import pandas as pd
import scanpy as sc
from huggingface_hub import hf_hub_download
from PIL import Image
from hest.HESTData import STHESTData
from hest.trident_compat import wsi_factory
from hest.utils import find_pixel_size_from_spot_coords, register_downscale_img, SpotPacking

# ── Args ───────────────────────────────────────────────────────────────────
if len(sys.argv) != 2:
    print("Usage: Assemble.py <path/to/cleaned_metadata.csv>")
    sys.exit(1)

CLEANED_META = sys.argv[1]

output_base = input("Where do you wish to store the data? Please input directory: ").strip()
DATA_DIR    = os.path.join(output_base, "data")
print(f"Notice: data will be stored under newly created folder: {DATA_DIR}")
print("This might take a long time. Converting...")
os.makedirs(DATA_DIR, exist_ok=True)

# ── Constants ──────────────────────────────────────────────────────────────
TECH_SUFFIX = {
    "ST":       "st",
    "Visium":   "visium",
    "VisiumHD": "visiumhd",
}

SPOT_PARAMS = {
    "ST":       {"spot_diameter": 100., "inter_spot_dist": 200.},
    "Visium":   {"spot_diameter": 55.,  "inter_spot_dist": 100.},
    "VisiumHD": {"spot_diameter": 2.,   "inter_spot_dist": 2.},
}

# Each entry: (target_pixel_size µm/px, output name)
# 256px patch at 0.5 µm/px = 128 µm footprint (20x equivalent)
# 256px patch at 1.0 µm/px = 256 µm footprint (10x equivalent)
# 256px patch at 2.0 µm/px = 512 µm footprint  (5x equivalent)
PATCH_SCALES = [
    (0.5, "patches_20x"),
    (1.0, "patches_10x"),
    (2.0, "patches_5x"),
]

# ── Metadata helpers ────────────────────────────────────────────────────────
def _first_title(raw):
    m = re.findall(r"Title \d+: (.+?)(?= Title \d+:|$)", raw)
    return m[0].strip() if m else raw.strip()

def _norm_species(s):
    return {"human": "Homo sapiens", "mouse": "Mus musculus"}.get(s.strip().lower(), s.strip().capitalize())

def _first_study_link(raw):
    pmids = [p.strip() for p in re.split(r"[\n,;]+", raw) if p.strip().isdigit()]
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmids[0]}/" if pmids else None

# ── Helpers ────────────────────────────────────────────────────────────────
def _load_stimage_data(slide_name, stimage_dir):
    img_array = np.array(Image.open(os.path.join(stimage_dir, f"{slide_name}.png")).convert("RGB"))
    counts    = pd.read_csv(os.path.join(stimage_dir, f"{slide_name}_count.csv"), index_col=0)
    coord     = pd.read_csv(os.path.join(stimage_dir, f"{slide_name}_coord.csv"), index_col=0)
    return img_array, counts, coord

def _tnbc_role(filename):
    """Classify a TNBC byArray filename into a required role, or None to skip."""
    if filename == "all.RData":
        return "counts"
    if "spot_data-all-" in filename and filename.endswith(".tsv"):
        return "coords"
    if "_HE-" in filename and filename.endswith(".jpg") and not filename.startswith("._"):
        return "he"
    return None

_TNBC_REQUIRED_ROLES = frozenset({"counts", "coords", "he"})

def _load_tnbc_data(array_id, data_path):
    try:
        import pyreadr
    except ImportError:
        raise RuntimeError("pip install pyreadr  (required to read TNBC .RData files)")

    array_dir = os.path.join(data_path, "byArray", array_id)
    files = os.listdir(array_dir)

    # H&E image — long dataset-specific name containing "_HE-"
    he_file = next(
        (f for f in files if "_HE-" in f and f.endswith(".jpg") and not f.startswith("._")),
        None
    )
    if he_file is None:
        raise RuntimeError(f"No H&E image found in {array_dir}. Files: {files}")
    img_array = np.array(Image.open(os.path.join(array_dir, he_file)).convert("RGB"))

    # Count matrix
    counts_obj = pyreadr.read_r(os.path.join(array_dir, "all.RData"))
    counts = list(counts_obj.values())[0]

    # Spot coordinates — prefer TSV (spot_data-all-*.tsv), fall back to RData
    tsv_file = next(
        (f for f in files if "spot_data-all-" in f and f.endswith(".tsv")), None
    )
    if tsv_file:
        spots = pd.read_csv(os.path.join(array_dir, tsv_file), sep="\t")
        # Build spot-ID index from grid x,y to match count matrix rownames
        spots.index = spots["x"].astype(int).astype(str) + "x" + spots["y"].astype(int).astype(str)
    else:
        rdata_file = next(
            (f for f in files if "allSpots" in f and f.endswith(".RData")), None
        )
        if rdata_file is None:
            raise RuntimeError(f"No spot coordinate file found in {array_dir}. Files: {files}")
        spots = list(pyreadr.read_r(os.path.join(array_dir, rdata_file)).values())[0]

    rename = {}
    for c in spots.columns:
        if c.lower() == "pixel_x":
            rename[c] = "xaxis"
        elif c.lower() == "pixel_y":
            rename[c] = "yaxis"
        elif c.lower() in ("x", "col") and "pixel_x" not in [cc.lower() for cc in spots.columns]:
            rename[c] = "xaxis"
        elif c.lower() in ("y", "row") and "pixel_y" not in [cc.lower() for cc in spots.columns]:
            rename[c] = "yaxis"
    coord = spots.rename(columns=rename)
    if "xaxis" not in coord.columns or "yaxis" not in coord.columns:
        raise RuntimeError(
            f"Spot file has no recognisable pixel_x/pixel_y columns. "
            f"Columns found: {list(spots.columns)}"
        )

    # Orient counts so rows = spots
    if not counts.index.isin(coord.index).any():
        counts = counts.T
    if not counts.index.isin(coord.index).any():
        raise RuntimeError(
            f"Cannot align counts to spot coordinates — no shared IDs.\n"
            f"  counts index (first 3): {list(counts.index[:3])}\n"
            f"  coord  index (first 3): {list(coord.index[:3])}"
        )

    return img_array, counts, coord

def _iter_tnbc_zenodo(record_id):
    """Stream byArray.tar from Zenodo, yielding (array_id, {filename: bytes}) one array at a time.
    Peak disk usage: zero (files buffered in memory, written to tmpdir by caller).
    """
    import requests
    url = f"https://zenodo.org/records/{record_id}/files/byArray.tar"
    print(f"  Connecting to {url} …", flush=True)
    with requests.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        resp.raw.decode_content = True
        buffers    = {}   # array_id -> {role -> (filename, bytes)}
        current    = None
        bytes_read = 0
        with tarfile.open(fileobj=resp.raw, mode="r|") as tar:
            for member in tar:
                if not member.isfile():
                    continue
                parts = member.name.replace("\\", "/").split("/")
                if len(parts) < 4:   # byArray/CN_N/subslide/filename
                    continue
                filename = parts[-1]
                array_id = f"{parts[-3]}_{parts[-2]}"  # e.g. CN1_C2
                role = _tnbc_role(filename)
                if role is None:
                    continue
                if array_id != current:
                    current = array_id
                    print(f"    {array_id}  ({bytes_read/1e6:.0f} MB streamed)…",
                          end="", flush=True)
                fobj = tar.extractfile(member)
                if fobj is None:
                    continue
                data = fobj.read()
                bytes_read += len(data)
                buffers.setdefault(array_id, {})[role] = (filename, data)
                if set(buffers[array_id].keys()) == _TNBC_REQUIRED_ROLES:
                    print(" ready", flush=True)
                    yield array_id, buffers.pop(array_id)

def _fix_tif_mpp(tif_path, pixel_size_um):
    """Fix HEST's tiff_save unit bug: pyvips writes mpp 10× too small.

    HEST computes xres as px/cm but pyvips interprets it as px/mm, then
    multiplies ×10 when encoding resunit=CM → stored mpp is 0.1× the truth.
    This corrects all IFD levels in-place without re-encoding any pixel data.
    """
    px_per_cm = 10000.0 / pixel_size_um          # correct value: px/cm
    rational = (round(px_per_cm * 1000), 1000)   # rational: num/den
    with open(tif_path, "r+b") as f:
        with tifffile.TiffFile(f) as tif:
            for page in tif.pages:
                for tag_name in ("XResolution", "YResolution"):
                    if tag_name in page.tags:
                        page.tags[tag_name].overwrite(rational)

# ── Download ───────────────────────────────────────────────────────────────
def download_slide(slide_name, technology):
    suffix      = TECH_SUFFIX.get(technology, technology.lower())
    stimage_dir = os.path.join(DATA_DIR, f"stimage_{suffix}", slide_name)
    hest_dir    = os.path.join(DATA_DIR, f"hest_{suffix}",    slide_name)
    os.makedirs(stimage_dir, exist_ok=True)
    os.makedirs(hest_dir,    exist_ok=True)

    for file_type, filename in {
        "coord":    f"{slide_name}_coord.csv",
        "gene_exp": f"{slide_name}_count.csv",
        "image":    f"{slide_name}.png",
    }.items():
        dest_path = os.path.join(stimage_dir, filename)
        if os.path.exists(dest_path):
            continue
        try:
            with tempfile.TemporaryDirectory(dir=DATA_DIR) as tmp:
                src = hf_hub_download(
                    repo_id="jiawennnn/STimage-1K4M",
                    filename=f"{technology}/{file_type}/{filename}",
                    repo_type="dataset",
                    local_dir=tmp,
                )
                shutil.move(src, dest_path)
        except Exception as e:
            print(f"    x {filename}: {e}")
            return None

    return {"stimage_dir": stimage_dir, "hest_dir": hest_dir, "technology": technology}

def _hest_folder(source, data_path, technology):
    if source in ("tnbc", "tnbc_zenodo"):
        return "hest_tnbc"
    if source == "visium_geo":
        return f"hest_{data_path.lower()}"            # hest_gse278936
    if source == "visium_ebi":
        return f"hest_{data_path.lower().replace('-', '_')}"  # hest_e_mtab_13530
    if source == "visium_zenodo":
        return "hest_visium_tls"
    return f"hest_{TECH_SUFFIX.get(technology, technology.lower())}"

def _load_spaceranger(slide_dir, img_path=None):
    """Load a Space Ranger output directory (h5 or MTX) → (img_array, counts, coord)."""
    # ── counts ────────────────────────────────────────────────────────────────
    h5_path = os.path.join(slide_dir, "filtered_feature_bc_matrix.h5")
    fbc_dir = os.path.join(slide_dir, "filtered_feature_bc_matrix")
    if os.path.exists(h5_path):
        ca = sc.read_10x_h5(h5_path)
    elif os.path.isdir(fbc_dir):
        ca = sc.read_10x_mtx(fbc_dir, var_names="gene_symbols", cache=False)
    else:
        raise RuntimeError(f"No count matrix found in {slide_dir}")
    X      = ca.X.toarray() if hasattr(ca.X, "toarray") else np.array(ca.X)
    counts = pd.DataFrame(X, index=ca.obs_names, columns=ca.var_names)

    # ── spatial directory ─────────────────────────────────────────────────────
    spatial_dir = os.path.join(slide_dir, "spatial")
    if not os.path.isdir(spatial_dir):
        spatial_dir = slide_dir

    # ── tissue positions ──────────────────────────────────────────────────────
    pos_file = None
    for name in ("tissue_positions.csv", "tissue_positions_list.csv"):
        p = os.path.join(spatial_dir, name)
        if os.path.exists(p):
            pos_file = p; break
    if pos_file is None:
        for root, _, fns in os.walk(slide_dir):
            for fn in fns:
                if fn in ("tissue_positions.csv", "tissue_positions_list.csv"):
                    pos_file = os.path.join(root, fn); break
            if pos_file:
                break
    if pos_file is None:
        raise RuntimeError(f"No tissue_positions CSV found under {slide_dir}")

    raw = pd.read_csv(pos_file, header=None, dtype=str)
    if raw.iloc[0, 0].strip().lower() == "barcode":
        raw.columns = [c.strip() for c in raw.iloc[0]]
        raw = raw.iloc[1:].reset_index(drop=True)
    else:
        raw.columns = ["barcode", "in_tissue", "array_row", "array_col",
                       "pxl_row_in_fullres", "pxl_col_in_fullres"]
    raw = raw[raw["in_tissue"].str.strip() == "1"].copy()
    raw.index = raw["barcode"].str.strip()

    # ── scalefactors ──────────────────────────────────────────────────────────
    sf_file = os.path.join(spatial_dir, "scalefactors_json.json")
    if not os.path.exists(sf_file):
        for root, _, fns in os.walk(slide_dir):
            if "scalefactors_json.json" in fns:
                sf_file = os.path.join(root, "scalefactors_json.json"); break
    with open(sf_file) as fj:
        sf = json.load(fj)
    r_fullres = float(sf["spot_diameter_fullres"]) / 2.0

    # ── image ─────────────────────────────────────────────────────────────────
    using_hires = False
    if img_path is None or not os.path.exists(img_path):
        tifs = glob.glob(os.path.join(slide_dir, "*.tif")) + \
               glob.glob(os.path.join(slide_dir, "*.tiff"))
        if tifs:
            img_path = tifs[0]
        else:
            hp = os.path.join(spatial_dir, "tissue_hires_image.png")
            if os.path.exists(hp):
                img_path = hp; using_hires = True
    if img_path is None or not os.path.exists(img_path):
        raise RuntimeError(f"No image found in {slide_dir}")
    img_array = np.array(Image.open(img_path).convert("RGB"))

    # ── align barcodes ────────────────────────────────────────────────────────
    common = counts.index.intersection(raw.index)
    counts = counts.loc[common]
    raw    = raw.loc[common]

    if using_hires:
        hires_scalef = float(sf.get("tissue_hires_scalef", 1.0))
        coord = pd.DataFrame({
            "xaxis": raw["pxl_col_in_fullres"].astype(float) * hires_scalef,
            "yaxis": raw["pxl_row_in_fullres"].astype(float) * hires_scalef,
            "r":     r_fullres * hires_scalef,
        }, index=raw.index)
    else:
        coord = pd.DataFrame({
            "xaxis": raw["pxl_col_in_fullres"].astype(float),
            "yaxis": raw["pxl_row_in_fullres"].astype(float),
            "r":     r_fullres,
        }, index=raw.index)

    return img_array, counts, coord

def _download_geo_raw(geo_acc, cache_dir):
    try:
        import requests
    except ImportError:
        raise RuntimeError("pip install requests")
    os.makedirs(cache_dir, exist_ok=True)
    tar_path = os.path.join(cache_dir, f"{geo_acc}_RAW.tar")
    if os.path.exists(tar_path):
        print(f"  Using cached {geo_acc}_RAW.tar")
        return tar_path
    n   = geo_acc[3:]  # strip "GSE"
    nnn = (n[:-3] + "nnn") if len(n) > 3 else (n + "nnn")
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{nnn}"
           f"/{geo_acc}/suppl/{geo_acc}_RAW.tar")
    print(f"  Downloading {url} …")
    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(tar_path, "wb") as fh:
            mb = 0
            for chunk in r.iter_content(chunk_size=8_388_608):
                fh.write(chunk); mb += len(chunk) / 1e6
                print(f"\r  {mb:.0f} MB…", end="", flush=True)
    print()
    return tar_path

def _extract_all_geo_gsm(tar_path, gsm_out_dirs):
    """Single-pass extraction of multiple GSM samples from one RAW.tar.

    GEO file naming: {GSM}_{SAMPLE_LABEL}_{filetype}[.gz]
    e.g. GSM8557976_BPH_1_matrix.mtx.gz
         GSM8557976_BPH_1_scalefactors_json.json.gz  (extra .gz added by GEO)
    We match by suffix so the sample label in the middle is irrelevant.
    Spatial files that GEO extra-gzipped are decompressed on the fly.
    """
    import gzip
    # (suffix_to_match, dest_subdir, dest_basename, decompress_gz)
    _RULES = [
        ("_matrix.mtx.gz",               "filtered_feature_bc_matrix", "matrix.mtx.gz",               False),
        ("_barcodes.tsv.gz",              "filtered_feature_bc_matrix", "barcodes.tsv.gz",             False),
        ("_features.tsv.gz",              "filtered_feature_bc_matrix", "features.tsv.gz",             False),
        ("_genes.tsv.gz",                 "filtered_feature_bc_matrix", "genes.tsv.gz",                False),
        ("_matrix.mtx",                   "filtered_feature_bc_matrix", "matrix.mtx",                  False),
        ("_barcodes.tsv",                 "filtered_feature_bc_matrix", "barcodes.tsv",                False),
        ("_features.tsv",                 "filtered_feature_bc_matrix", "features.tsv",                False),
        ("_genes.tsv",                    "filtered_feature_bc_matrix", "genes.tsv",                   False),
        ("_scalefactors_json.json.gz",    "spatial", "scalefactors_json.json",      True),
        ("_scalefactors_json.json",       "spatial", "scalefactors_json.json",      False),
        ("_tissue_positions_list.csv.gz", "spatial", "tissue_positions_list.csv",   True),
        ("_tissue_positions_list.csv",    "spatial", "tissue_positions_list.csv",   False),
        ("_tissue_positions.csv.gz",      "spatial", "tissue_positions.csv",        True),
        ("_tissue_positions.csv",         "spatial", "tissue_positions.csv",        False),
        ("_tissue_hires_image.png.gz",    "spatial", "tissue_hires_image.png",      True),
        ("_tissue_hires_image.png",       "spatial", "tissue_hires_image.png",      False),
        ("_tissue_lowres_image.png.gz",   "spatial", "tissue_lowres_image.png",     True),
        ("_tissue_lowres_image.png",      "spatial", "tissue_lowres_image.png",     False),
    ]
    for out_dir in gsm_out_dirs.values():
        os.makedirs(os.path.join(out_dir, "filtered_feature_bc_matrix"), exist_ok=True)
        os.makedirs(os.path.join(out_dir, "spatial"), exist_ok=True)
    prefix_map = {gsm_id + "_": (gsm_id, out_dir)
                  for gsm_id, out_dir in gsm_out_dirs.items()}
    print(f"  Extracting {len(gsm_out_dirs)} samples from {os.path.basename(tar_path)} …",
          flush=True)
    with tarfile.open(tar_path, "r") as tar:
        for member in tar.getmembers():
            fname = os.path.basename(member.name)
            matched = next(((gsm_id, od) for pfx, (gsm_id, od) in prefix_map.items()
                            if fname.startswith(pfx)), None)
            if matched is None:
                continue
            gsm_id, out_dir = matched
            rule = next(((sub, bn, dec) for sfx, sub, bn, dec in _RULES
                         if fname.endswith(sfx)), None)
            if rule is None:
                continue
            subdir, dest_basename, decompress = rule
            dest = os.path.join(out_dir, subdir, dest_basename)
            if os.path.exists(dest):
                continue
            fobj = tar.extractfile(member)
            if fobj is None:
                continue
            data = fobj.read()
            if decompress:
                data = gzip.decompress(data)
            with open(dest, "wb") as fh:
                fh.write(data)

def _download_tls_zip(record_url, cache_dir):
    try:
        import requests
    except ImportError:
        raise RuntimeError("pip install requests")
    import zipfile as _zf
    os.makedirs(cache_dir, exist_ok=True)
    zip_path  = os.path.join(cache_dir, "TLS_VISIUM_USZ.zip")
    record_id = record_url.rstrip("/").split("/")[-1]
    url       = f"https://zenodo.org/records/{record_id}/files/TLS_VISIUM_USZ.zip"

    # If file exists and is a valid zip, use it as-is
    if os.path.exists(zip_path):
        try:
            with _zf.ZipFile(zip_path): pass
            print("  Using cached TLS_VISIUM_USZ.zip")
            return zip_path
        except _zf.BadZipFile:
            pass  # partial — resume below

    MAX_RETRIES = 10
    for attempt in range(1, MAX_RETRIES + 1):
        existing = os.path.getsize(zip_path) if os.path.exists(zip_path) else 0
        headers  = {"Range": f"bytes={existing}-"} if existing else {}
        mode     = "ab" if existing else "wb"
        if existing:
            print(f"  Resuming TLS_VISIUM_USZ.zip from {existing/1e6:.0f} MB "
                  f"(attempt {attempt}/{MAX_RETRIES})…")
        else:
            print(f"  Downloading {url} (~2 GB) …")
        try:
            with requests.get(url, stream=True, headers=headers, timeout=600) as r:
                r.raise_for_status()
                total = existing + int(r.headers.get("Content-Length", 0))
                with open(zip_path, mode) as fh:
                    mb = existing / 1e6
                    for chunk in r.iter_content(chunk_size=8_388_608):
                        fh.write(chunk); mb += len(chunk) / 1e6
                        pct = f" ({mb/total*100:.0f}%)" if total else ""
                        print(f"\r  {mb:.0f} MB{pct}…", end="", flush=True)
            print()
            # Validate before returning
            try:
                with _zf.ZipFile(zip_path): pass
                return zip_path
            except _zf.BadZipFile:
                print(f"  ZIP invalid after download (attempt {attempt}); retrying…")
        except Exception as exc:
            print(f"\n  Download error (attempt {attempt}/{MAX_RETRIES}): {exc}")
            if attempt == MAX_RETRIES:
                raise RuntimeError(
                    f"Failed to download TLS zip after {MAX_RETRIES} attempts"
                ) from exc
            print("  Retrying…")

    raise RuntimeError("Failed to download TLS zip")

def _load_tls_sample(sample_name, zip_path):
    with zipfile.ZipFile(zip_path, "r") as zf:
        with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp:
            tmp.write(zf.read(f"tif_slides/{sample_name}.tif"))
            tif_tmp = tmp.name
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
            tmp.write(zf.read(f"h5ad_preprocessed/{sample_name}.h5ad"))
            h5ad_tmp = tmp.name
    try:
        img_array = np.array(Image.open(tif_tmp).convert("RGB"))
        adata     = sc.read(h5ad_tmp)
    finally:
        for p in (tif_tmp, h5ad_tmp):
            try: os.unlink(p)
            except: pass

    if hasattr(adata, "raw") and adata.raw is not None:
        X = adata.raw.X; var_names = adata.raw.var_names
    elif "counts" in adata.layers:
        X = adata.layers["counts"]; var_names = adata.var_names
    else:
        X = adata.X; var_names = adata.var_names
    if hasattr(X, "toarray"):
        X = X.toarray()
    counts = pd.DataFrame(X, index=adata.obs_names, columns=var_names)

    # Squidpy/Visium convention: spatial[:, 0]=x (col), spatial[:, 1]=y (row)
    sp = adata.obsm["spatial"]
    coord = pd.DataFrame({
        "xaxis": sp[:, 0].astype(float),
        "yaxis": sp[:, 1].astype(float),
        # no 'r' — convert_slide falls back to find_pixel_size_from_spot_coords
    }, index=adata.obs_names)

    return img_array, counts, coord

def _download_ebi_sample(sample_name, ebi_acc, cache_dir):
    """Download {sample}-filtered_feature_bc_matrix.h5 and {sample}-spatial.tar
    from EBI BioStudies (flat file layout, no per-sample subdirs in the archive).
    """
    try:
        import requests
    except ImportError:
        raise RuntimeError("pip install requests")
    sample_dir  = os.path.join(cache_dir, "ebi", ebi_acc, sample_name)
    h5_dest     = os.path.join(sample_dir, "filtered_feature_bc_matrix.h5")
    spatial_dir = os.path.join(sample_dir, "spatial")
    base_url    = f"https://www.ebi.ac.uk/biostudies/files/{ebi_acc}"

    if os.path.exists(h5_dest) and os.path.isdir(spatial_dir):
        print(f"    Using cached {sample_name}/")
        return sample_dir

    os.makedirs(sample_dir,  exist_ok=True)
    os.makedirs(spatial_dir, exist_ok=True)

    def _fetch(remote_name, dest_path, label):
        if os.path.exists(dest_path):
            return
        url = f"{base_url}/{remote_name}"
        print(f"    {label}…", end="", flush=True)
        with requests.get(url, stream=True, timeout=300) as fr:
            fr.raise_for_status()
            with open(dest_path, "wb") as fh:
                for chunk in fr.iter_content(chunk_size=8_388_608):
                    fh.write(chunk)
        print(" ok")

    # Count matrix (h5)
    _fetch(f"{sample_name}-filtered_feature_bc_matrix.h5", h5_dest, "h5")

    # Spatial tar → extract into spatial/
    tar_tmp = os.path.join(sample_dir, "_spatial.tar")
    _fetch(f"{sample_name}-spatial.tar", tar_tmp, "spatial.tar")
    if os.path.exists(tar_tmp):
        with tarfile.open(tar_tmp, "r:*") as tar:
            for member in tar.getmembers():
                fname = os.path.basename(member.name)
                if not fname or member.isdir():
                    continue
                dest = os.path.join(spatial_dir, fname)
                if os.path.exists(dest):
                    continue
                fobj = tar.extractfile(member)
                if fobj:
                    with open(dest, "wb") as fh:
                        fh.write(fobj.read())
        os.unlink(tar_tmp)

    return sample_dir

# ── Convert ────────────────────────────────────────────────────────────────
def convert_slide(slide_name, hest_dir, technology, row, img_array, counts, coord):
    common        = counts.index.intersection(coord.index)
    counts, coord = counts.loc[common], coord.loc[common]

    if technology == "ST":
        try:
            short_ids = [re.search(r"(\d+x\d+)$", idx).group(1) for idx in counts.index]
        except AttributeError:
            short_ids = list(counts.index.astype(str))
    else:
        prefix    = slide_name + "_"
        short_ids = [idx[len(prefix):] if idx.startswith(prefix) else idx for idx in counts.index]
    counts.index = coord.index = short_ids

    adata = sc.AnnData(counts)
    adata.obsm["spatial"] = coord[["xaxis", "yaxis"]].values

    my_df = pd.DataFrame(
        adata.obsm["spatial"], index=adata.obs_names,
        columns=["pxl_col_in_fullres", "pxl_row_in_fullres"]
    )
    if technology == "ST":
        try:
            my_df["array_row"] = [round(float(i.split("x")[0])) for i in my_df.index]
            my_df["array_col"] = [round(float(i.split("x")[1])) for i in my_df.index]
        except (ValueError, IndexError):
            _, y_inv = np.unique(my_df["pxl_row_in_fullres"].round().astype(int).values, return_inverse=True)
            _, x_inv = np.unique(my_df["pxl_col_in_fullres"].round().astype(int).values, return_inverse=True)
            my_df["array_row"] = y_inv
            my_df["array_col"] = x_inv
    else:
        _, y_inv = np.unique(my_df["pxl_row_in_fullres"].round().astype(int).values, return_inverse=True)
        _, x_inv = np.unique(my_df["pxl_col_in_fullres"].round().astype(int).values, return_inverse=True)
        my_df["array_row"] = y_inv
        my_df["array_col"] = x_inv

    adata.obs["array_row"]          = my_df["array_row"].values
    adata.obs["array_col"]          = my_df["array_col"].values
    adata.obs["pxl_col_in_fullres"] = my_df["pxl_col_in_fullres"].values
    adata.obs["pxl_row_in_fullres"] = my_df["pxl_row_in_fullres"].values
    adata.obs["in_tissue"]          = True

    if technology == "ST":
        pixel_size, spot_estimate_dist = find_pixel_size_from_spot_coords(
            my_df, inter_spot_dist=SPOT_PARAMS[technology]["inter_spot_dist"],
            packing=SpotPacking.GRID_PACKING
        )
    else:
        if "r" in coord.columns:
            # spot-radius column present (Space Ranger output) — most accurate
            pixel_size = SPOT_PARAMS[technology]["spot_diameter"] / 2.0 / coord["r"].median()
        else:
            # no radius known (pre-processed h5ad); estimate from spatial grid
            pixel_size, _ = find_pixel_size_from_spot_coords(
                my_df, inter_spot_dist=SPOT_PARAMS[technology]["inter_spot_dist"],
                packing=SpotPacking.GRID_PACKING
            )
        spot_estimate_dist = SPOT_PARAMS[technology]["inter_spot_dist"] / pixel_size

    # wsi_factory raises ValueError for mpp >= 2.4; STimage Visium PNGs are
    # ~2000px thumbnails (~3-9 um/px). Upscale to 1.0 um/px so HEST accepts them.
    if pixel_size >= 2.0:
        scale = pixel_size / 1.0
        h, w = img_array.shape[:2]
        img_array = np.array(Image.fromarray(img_array).resize(
            (round(w * scale), round(h * scale)), Image.LANCZOS
        ))
        my_df["pxl_col_in_fullres"] *= scale
        my_df["pxl_row_in_fullres"] *= scale
        adata.obs["pxl_col_in_fullres"] = my_df["pxl_col_in_fullres"].values
        adata.obs["pxl_row_in_fullres"] = my_df["pxl_row_in_fullres"].values
        adata.obsm["spatial"]           = my_df[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values
        pixel_size       /= scale
        spot_estimate_dist *= scale

    wsi = wsi_factory(img_array, mpp=pixel_size)
    register_downscale_img(adata, wsi, pixel_size,
                           spot_size=SPOT_PARAMS[technology]["spot_diameter"])

    if technology == "ST" and all(re.match(r"^\d+x\d+$", i) for i in list(adata.obs.index)[:5]):
        adata.obs.index = [
            i.split("x")[0].zfill(3) + "x" + i.split("x")[1].zfill(3)
            for i in adata.obs.index
        ]

    adata.var["mito"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
    adata.obs["n_counts"] = adata.obs["total_counts"]

    gse_m = re.match(r"(GSE\d+)_(GSM\d+)", slide_name)
    gene_num = row.get("gene_num", "").strip()

    meta = {
        # computed
        "pixel_size_um_embedded":  None,
        "pixel_size_um_estimated": pixel_size,
        "fullres_height":          img_array.shape[0],
        "fullres_width":           img_array.shape[1],
        "spots_under_tissue":      len(adata.obs),
        "spot_estimate_dist":      int(spot_estimate_dist),
        "spot_diameter":           SPOT_PARAMS[technology]["spot_diameter"],
        "inter_spot_dist":         SPOT_PARAMS[technology]["inter_spot_dist"],
        # from STimage metadata
        "id":                      slide_name,
        "image_filename":          f"{slide_name}.tif",
        "dataset_title":           _first_title(row["title"]),
        "organ":                   row["tissue"].strip() or None,
        "tissue":                  row["tissue"].strip() or None,
        "species":                 _norm_species(row["species"]),
        "st_technology":           row["tech"],
        "disease_state":           "Cancer" if str(row["involve_cancer"]).strip().lower() in ("true", "1", "yes") else "Normal",
        "nb_genes":                int(gene_num) if gene_num.isdigit() else None,
        "study_link":              _first_study_link(row["pmid"]),
        "download_page_link1":     f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_m.group(1)}" if gse_m else None,
        "subseries":               gse_m.group(2) if gse_m else slide_name,
        # not available in STimage
        "patient":                 None,
        "oncotree_code":           None,
        "data_publication_date":   None,
        "license":                 None,
        "preservation_method":     None,
        "magnification":           None,
        "treatment_comment":       None,
        "disease_comment":         None,
    }

    st = STHESTData(adata, img_array, pixel_size, meta)
    st.save(path=hest_dir, save_img=True, pyramidal=True, bigtiff=False, plot_pxl_size=True)
    _fix_tif_mpp(os.path.join(hest_dir, "aligned_fullres_HE.tif"), pixel_size)
    st.save_spatial_plot(save_path=hest_dir)
    st.segment_tissue(method="otsu")
    st.save_tissue_contours(hest_dir, "tissue_seg")
    st.save_tissue_vis(hest_dir, "tissue_seg")
    for target_pixel_size, name in PATCH_SCALES:
        st.dump_patches(patch_save_dir=hest_dir, name=name,
                        target_patch_size=256, target_pixel_size=target_pixel_size,
                        dump_visualization=True)

# ── Main loop ──────────────────────────────────────────────────────────────
with open(CLEANED_META, newline="") as f:
    slides = list(csv.DictReader(f))

# Separate sources that need batch/streaming downloads from per-slide sources
zenodo_groups        = {}   # data_path -> {array_id -> row}   (TNBC streaming)
geo_groups           = {}   # data_path -> {gsm_id   -> row}   (GEO batch extract)
visium_zenodo_groups = {}   # data_path -> {sample   -> row}   (TLS zip download)
other_slides         = []
for row in slides:
    src = row.get("source", "")
    if src == "tnbc_zenodo":
        zenodo_groups.setdefault(row["data_path"], {})[row["slide"]] = row
    elif src == "visium_geo":
        geo_groups.setdefault(row["data_path"], {})[row["slide"]] = row
    elif src == "visium_zenodo":
        visium_zenodo_groups.setdefault(row["data_path"], {})[row["slide"]] = row
    else:
        other_slides.append(row)

failed  = []
skipped = 0
done    = 0

# ── Non-streaming sources ───────────────────────────────────────────────────
for i, row in enumerate(other_slides):
    slide_name = row["slide"]
    technology = row["tech"]
    source     = row.get("source", "stimage")
    folder     = _hest_folder(source, row.get("data_path", ""), technology)
    hest_dir   = os.path.join(DATA_DIR, folder, slide_name)

    if os.path.exists(os.path.join(hest_dir, "aligned_adata.h5ad")):
        skipped += 1
        continue

    print(f"[{i+1}/{len(other_slides)}] {slide_name} ({technology}) [{source}]")
    try:
        if source == "stimage":
            result = download_slide(slide_name, technology)
            if result is None:
                raise RuntimeError("download failed")
            img_array, counts, coord = _load_stimage_data(slide_name, result["stimage_dir"])
        elif source == "tnbc":
            os.makedirs(hest_dir, exist_ok=True)
            img_array, counts, coord = _load_tnbc_data(slide_name, row["data_path"])
        elif source == "visium_ebi":
            ebi_acc    = row["data_path"]
            sample_dir = _download_ebi_sample(
                slide_name, ebi_acc, os.path.join(DATA_DIR, "raw_cache"))
            os.makedirs(hest_dir, exist_ok=True)
            img_array, counts, coord = _load_spaceranger(sample_dir)
        else:
            raise ValueError(f"Unknown source: {source!r}")
        convert_slide(slide_name, hest_dir, technology, row, img_array, counts, coord)
        done += 1
        print(f"  done")
    except Exception as e:
        print(f"  ERROR: {e}")
        failed.append({"slide": slide_name, "technology": technology, "error": str(e)})

# ── Zenodo streaming sources (one pass per record) ─────────────────────────
for data_path, array_rows in zenodo_groups.items():
    record_id  = data_path.rstrip("/").split("/")[-1]
    sample_row = next(iter(array_rows.values()))
    technology = sample_row["tech"]
    suffix     = TECH_SUFFIX.get(technology, technology.lower())
    folder     = "hest_tnbc"

    n_todo = sum(
        1 for aid in array_rows
        if not os.path.exists(os.path.join(DATA_DIR, folder, aid, "aligned_adata.h5ad"))
    )
    print(f"\n[tnbc_zenodo] record={record_id}  arrays={len(array_rows)}  to_do={n_todo}")
    if n_todo == 0:
        skipped += len(array_rows)
        continue

    for array_id, file_bufs in _iter_tnbc_zenodo(record_id):
        if array_id not in array_rows:
            continue
        row      = array_rows[array_id]
        hest_dir = os.path.join(DATA_DIR, folder, array_id)
        if os.path.exists(os.path.join(hest_dir, "aligned_adata.h5ad")):
            skipped += 1
            continue
        print(f"  {array_id}")
        try:
            with tempfile.TemporaryDirectory(dir=DATA_DIR) as tmpdir:
                arr_dir = os.path.join(tmpdir, "byArray", array_id)
                os.makedirs(arr_dir)
                for role, (fname, data) in file_bufs.items():
                    with open(os.path.join(arr_dir, fname), "wb") as fh:
                        fh.write(data)
                img_array, counts, coord = _load_tnbc_data(array_id, tmpdir)
            convert_slide(array_id, hest_dir, technology, row, img_array, counts, coord)
            done += 1
            print(f"    done")
        except Exception as e:
            print(f"    ERROR: {e}")
            failed.append({"slide": array_id, "technology": technology, "error": str(e)})

# ── GEO batch (one tar pass per accession) ────────────────────────────────
for geo_acc, gsm_rows in geo_groups.items():
    sample_row = next(iter(gsm_rows.values()))
    technology = sample_row["tech"]
    folder     = _hest_folder("visium_geo", geo_acc, technology)
    geo_cache  = os.path.join(DATA_DIR, "raw_cache", "geo")

    n_todo = sum(
        1 for gsm in gsm_rows
        if not os.path.exists(os.path.join(DATA_DIR, folder, gsm, "aligned_adata.h5ad"))
    )
    print(f"\n[visium_geo] acc={geo_acc}  samples={len(gsm_rows)}  to_do={n_todo}")
    if n_todo == 0:
        skipped += len(gsm_rows)
        continue

    try:
        tar_path = _download_geo_raw(geo_acc, geo_cache)
        gsm_out_dirs = {gsm: os.path.join(geo_cache, geo_acc, gsm)
                        for gsm in gsm_rows}
        _extract_all_geo_gsm(tar_path, gsm_out_dirs)
    except Exception as e:
        print(f"  ERROR downloading/extracting {geo_acc}: {e}")
        for gsm in gsm_rows:
            failed.append({"slide": gsm, "technology": technology, "error": str(e)})
        continue

    for slide_name, row in gsm_rows.items():
        hest_dir = os.path.join(DATA_DIR, folder, slide_name)
        if os.path.exists(os.path.join(hest_dir, "aligned_adata.h5ad")):
            skipped += 1; continue
        print(f"  {slide_name}")
        try:
            os.makedirs(hest_dir, exist_ok=True)
            img_array, counts, coord = _load_spaceranger(gsm_out_dirs[slide_name])
            convert_slide(slide_name, hest_dir, technology, row, img_array, counts, coord)
            done += 1; print(f"    done")
        except Exception as e:
            print(f"    ERROR: {e}")
            failed.append({"slide": slide_name, "technology": technology, "error": str(e)})

# ── TLS Zenodo batch (one zip download per record) ─────────────────────────
for record_url, sample_rows in visium_zenodo_groups.items():
    sample_row = next(iter(sample_rows.values()))
    technology = sample_row["tech"]
    folder     = "hest_visium_tls"
    tls_cache  = os.path.join(DATA_DIR, "raw_cache", "tls")

    n_todo = sum(
        1 for sn in sample_rows
        if not os.path.exists(os.path.join(DATA_DIR, folder, sn, "aligned_adata.h5ad"))
    )
    print(f"\n[visium_zenodo] record={record_url.split('/')[-1]}  "
          f"samples={len(sample_rows)}  to_do={n_todo}")
    if n_todo == 0:
        skipped += len(sample_rows)
        continue

    try:
        zip_path = _download_tls_zip(record_url, tls_cache)
    except Exception as e:
        print(f"  ERROR downloading TLS zip: {e}")
        for sn in sample_rows:
            failed.append({"slide": sn, "technology": technology, "error": str(e)})
        continue

    for slide_name, row in sample_rows.items():
        hest_dir = os.path.join(DATA_DIR, folder, slide_name)
        if os.path.exists(os.path.join(hest_dir, "aligned_adata.h5ad")):
            skipped += 1; continue
        print(f"  {slide_name}")
        try:
            os.makedirs(hest_dir, exist_ok=True)
            img_array, counts, coord = _load_tls_sample(slide_name, zip_path)
            convert_slide(slide_name, hest_dir, technology, row, img_array, counts, coord)
            done += 1; print(f"    done")
        except Exception as e:
            print(f"    ERROR: {e}")
            failed.append({"slide": slide_name, "technology": technology, "error": str(e)})

print(f"\nFinished: {done} converted, {skipped} skipped (already done), {len(failed)} failed")
if failed:
    print("\nFailed slides:")
    for f in failed:
        print(f"  {f['slide']} ({f['technology']}): {f['error']}")
