#!/usr/bin/env python3
# =============================================================================
# Assemble.py
# =============================================================================
# PURPOSE: Take the cleaned_metadata.csv produced by Clean.py and actually
# download + convert every slide into HEST format.
#
# WHAT IT DOES IN PLAIN ENGLISH:
#   For each tissue slide listed in the CSV:
#   1. Download the raw data (H&E image + gene expression counts + spot
#      coordinates) from whichever public repository it lives in.
#   2. Align the spot coordinates to the image so we know exactly where on
#      the tissue each sequenced spot sits.
#   3. Calculate quality metrics (total counts per spot, mitochondrial
#      percentage, etc.).
#   4. Save everything into the standardised HEST format:
#        - aligned_adata.h5ad     : gene expression + coordinates (scanpy format)
#        - aligned_fullres_HE.tif : whole-slide image as a pyramidal TIFF
#        - patches_20x/10x/5x.h5  : 256×256 pixel image crops around each spot
#        - tissue_seg_*.geojson   : tissue mask (which pixels are tissue vs. background)
#        - spatial_plots.png      : visualisation of spots overlaid on the image
#        - pixel_size_vis.png     : visualisation of the calibrated pixel size
#
# HOW DOWNLOADING WORKS PER DATASET:
#   STimage-1K4M  → downloads individual files from HuggingFace (one at a time)
#   TNBC Zenodo   → streams byArray.tar from Zenodo without saving the full tar
#   GEO (GSE*)    → downloads one big RAW.tar from NCBI FTP, extracts all samples
#   EBI (E-MTAB*) → downloads one .h5 + one spatial.tar per sample from EBI
#   TLS Zenodo    → downloads one big ZIP from Zenodo (with resume on failure)
# =============================================================================

# ── Standard library imports ──────────────────────────────────────────────────
# csv      : read the cleaned_metadata.csv produced by Clean.py
# glob     : find files matching a wildcard pattern (e.g. "*.tif")
# io       : work with in-memory byte streams
# json     : read JSON files (scalefactors_json.json from Space Ranger)
# os       : interact with the filesystem
# re       : regular expressions for pattern matching
# shutil   : high-level file operations (move, copy, delete directory trees)
# sys      : read command-line arguments, exit the program
# tarfile  : open .tar archive files
# tempfile : create temporary files and directories
# zipfile  : open .zip archive files (used for TLS dataset)
import csv, glob, io, json, os, re, shutil, sys, tarfile, tempfile, zipfile

# ── Third-party library imports ───────────────────────────────────────────────
# numpy (np)         : fast numerical arrays — we use it to manipulate image
#                      pixel arrays and spot coordinate matrices
# tifffile           : read/write TIFF image files at the byte level, used to
#                      fix metadata (resolution tags) in the output TIFF
# pandas (pd)        : data tables — we store counts and coordinates as DataFrames
# scanpy (sc)        : the standard Python toolkit for single-cell / spatial
#                      transcriptomics; handles the .h5ad file format
# huggingface_hub    : download files from HuggingFace datasets (STimage lives there)
# PIL.Image          : open and resize images (JPEG, PNG, TIFF)
# hest.HESTData      : the HEST library's main data container — wraps scanpy AnnData
#                      and adds methods like save(), segment_tissue(), dump_patches()
# hest.trident_compat: wsi_factory() creates a "whole slide image" object from a
#                      numpy array so HEST can process it
# hest.utils         : helper functions — find_pixel_size_from_spot_coords() estimates
#                      the image resolution from the spacing between spots;
#                      register_downscale_img() embeds a thumbnail in the h5ad;
#                      SpotPacking tells the estimator whether spots are in a square
#                      grid (ST) or hexagonal grid (Visium)
import numpy as np
import tifffile
import pandas as pd
import scanpy as sc
from huggingface_hub import hf_hub_download
from PIL import Image
from hest.HESTData import STHESTData
from hest.trident_compat import wsi_factory
from hest.utils import find_pixel_size_from_spot_coords, register_downscale_img, SpotPacking

# =============================================================================
# STEP 1 — Read the command-line argument
# =============================================================================
if len(sys.argv) != 2:
    print("Usage: Assemble.py <path/to/cleaned_metadata.csv>")
    sys.exit(1)

# Path to the CSV that Clean.py produced.
CLEANED_META = sys.argv[1]

# Ask the user where to store output.  Same directory as Clean.py used.
output_base = input("Where do you wish to store the data? Please input directory: ").strip()
# All converted HEST folders go inside a "data/" subfolder.
DATA_DIR = os.path.join(output_base, "data")
print(f"Notice: data will be stored under newly created folder: {DATA_DIR}")
print("This might take a long time. Converting...")
os.makedirs(DATA_DIR, exist_ok=True)

# =============================================================================
# STEP 2 — Constants
# =============================================================================

# Maps the technology name to a short suffix used in the output folder name.
# e.g. technology="Visium" → output folder named "hest_visium"
TECH_SUFFIX = {
    "ST":       "st",
    "Visium":   "visium",
    "VisiumHD": "visiumhd",
}

# Physical dimensions of capture spots, in micrometres (µm), per technology.
# These are fixed hardware specifications — not estimated from data.
#
#   ST (Spatial Transcriptomics, older technology):
#     - Spots are 100 µm in diameter
#     - Spots are 200 µm apart (centre-to-centre)
#     - Arranged in a regular square grid
#
#   Visium (10x Genomics, current standard):
#     - Spots are 55 µm in diameter
#     - Spots are 100 µm apart (centre-to-centre)
#     - Arranged in a hexagonal grid
#
#   VisiumHD (10x Genomics, newest high-definition):
#     - 2 µm bins (essentially continuous coverage, no gaps)
SPOT_PARAMS = {
    "ST":       {"spot_diameter": 100., "inter_spot_dist": 200.},
    "Visium":   {"spot_diameter": 55.,  "inter_spot_dist": 100.},
    "VisiumHD": {"spot_diameter": 2.,   "inter_spot_dist": 2.},
}

# The three image resolutions at which we extract 256×256 pixel patches.
# A "patch" is a square crop of the H&E image centred on one sequenced spot.
# We extract the same spot at three different physical scales so a downstream
# model can choose the level of detail it needs:
#
#   0.5 µm/px  →  each pixel represents 0.5 µm of real tissue
#               →  a 256×256 patch covers 128×128 µm  ≈ 20x microscope objective
#               →  fine detail: you can see individual cell nuclei clearly
#
#   1.0 µm/px  →  each pixel = 1 µm
#               →  a 256×256 patch covers 256×256 µm  ≈ 10x objective
#               →  medium detail: cell clusters and tissue architecture visible
#
#   2.0 µm/px  →  each pixel = 2 µm
#               →  a 256×256 patch covers 512×512 µm  ≈ 5x objective
#               →  low detail: large-scale tissue structure, neighbourhood context
#
# To extract a patch at 0.5 µm/px from an image whose native resolution is
# e.g. 0.25 µm/px, HEST downsamples the image by 2× before cropping.
# The output patch is always 256×256 pixels regardless of the scale.
PATCH_SCALES = [
    (0.5, "patches_20x"),   # (target µm/px, output folder name)
    (1.0, "patches_10x"),
    (2.0, "patches_5x"),
]

# =============================================================================
# STEP 3 — Small metadata helper functions
# =============================================================================

def _first_title(raw):
    """Extract the first paper title from STimage's combined title string.

    STimage sometimes packs multiple papers into one cell like:
      "Title 1: Foo paper. Title 2: Bar paper."
    We extract just the first one.  If the string doesn't follow that pattern,
    we return it as-is.
    """
    m = re.findall(r"Title \d+: (.+?)(?= Title \d+:|$)", raw)
    return m[0].strip() if m else raw.strip()

def _norm_species(s):
    """Convert informal species names to their formal scientific names.

    "human" → "Homo sapiens"
    "mouse" → "Mus musculus"
    Anything else → returned capitalised but unchanged.
    """
    return {"human": "Homo sapiens", "mouse": "Mus musculus"}.get(
        s.strip().lower(), s.strip().capitalize()
    )

def _first_study_link(raw):
    """Turn a PubMed ID (or comma-separated list of IDs) into a PubMed URL.

    "39550375"  →  "https://pubmed.ncbi.nlm.nih.gov/39550375/"
    If no valid integer ID is found, returns None.
    """
    pmids = [p.strip() for p in re.split(r"[\n,;]+", raw) if p.strip().isdigit()]
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmids[0]}/" if pmids else None

# =============================================================================
# STEP 4 — STimage data loader
# =============================================================================

def _load_stimage_data(slide_name, stimage_dir):
    """Read the three files that make up one STimage slide from a local directory.

    STimage stores each slide as three separate files:
      {id}.png           — the H&E image (full resolution, RGB)
      {id}_count.csv     — gene expression matrix (rows=spots, cols=genes)
      {id}_coord.csv     — spot coordinates (xaxis, yaxis, r columns)

    Returns a tuple: (img_array, counts, coord)
      img_array : numpy array of shape (height, width, 3) — the image pixels
      counts    : pandas DataFrame (spots × genes) — RNA counts
      coord     : pandas DataFrame with columns xaxis, yaxis, r (spot radius)
    """
    img_array = np.array(Image.open(os.path.join(stimage_dir, f"{slide_name}.png")).convert("RGB"))
    counts    = pd.read_csv(os.path.join(stimage_dir, f"{slide_name}_count.csv"), index_col=0)
    coord     = pd.read_csv(os.path.join(stimage_dir, f"{slide_name}_coord.csv"), index_col=0)
    return img_array, counts, coord

# =============================================================================
# STEP 5 — TNBC data helpers
# =============================================================================

def _tnbc_role(filename):
    """Classify a filename from the TNBC byArray directory into a data role.

    The TNBC dataset stores three types of file per array:
      all.RData                      → the gene expression count matrix (R format)
      spot_data-all-{n}.tsv         → spot pixel coordinates (tab-separated)
      {long_name}_HE-{id}.jpg       → the H&E image

    Returns one of "counts", "coords", "he", or None (if the file is irrelevant).
    We use this to know which file serves which purpose when streaming from Zenodo.
    """
    if filename == "all.RData":
        return "counts"
    if "spot_data-all-" in filename and filename.endswith(".tsv"):
        return "coords"
    if "_HE-" in filename and filename.endswith(".jpg") and not filename.startswith("._"):
        # Filenames starting with "._" are macOS metadata artifacts — skip them.
        return "he"
    return None

# The set of roles every array must have before we can convert it.
_TNBC_REQUIRED_ROLES = frozenset({"counts", "coords", "he"})

def _load_tnbc_data(array_id, data_path):
    """Load H&E image, counts, and spot coordinates for one TNBC array.

    The TNBC dataset (whether downloaded from Zenodo or local) is organised as:
      {data_path}/byArray/{array_id}/
        all.RData                   ← gene expression (R binary format)
        spot_data-all-{n}.tsv       ← spot coordinates
        {long_name}_HE-{id}.jpg     ← H&E image

    Returns: (img_array, counts, coord)
    """
    try:
        import pyreadr   # reads R binary data files (.RData, .RDS) into pandas
    except ImportError:
        raise RuntimeError("pip install pyreadr  (required to read TNBC .RData files)")

    array_dir = os.path.join(data_path, "byArray", array_id)
    files = os.listdir(array_dir)

    # --- H&E image ---
    # The image filename is long and dataset-specific, but always contains "_HE-".
    he_file = next(
        (f for f in files if "_HE-" in f and f.endswith(".jpg") and not f.startswith("._")),
        None
    )
    if he_file is None:
        raise RuntimeError(f"No H&E image found in {array_dir}. Files: {files}")
    img_array = np.array(Image.open(os.path.join(array_dir, he_file)).convert("RGB"))

    # --- Count matrix ---
    # all.RData contains a single R data frame with spots as rows and genes as columns.
    counts_obj = pyreadr.read_r(os.path.join(array_dir, "all.RData"))
    counts = list(counts_obj.values())[0]   # get the first (only) object from the R file

    # --- Spot coordinates ---
    # Prefer the TSV file; fall back to an RData file if TSV doesn't exist.
    tsv_file = next(
        (f for f in files if "spot_data-all-" in f and f.endswith(".tsv")), None
    )
    if tsv_file:
        spots = pd.read_csv(os.path.join(array_dir, tsv_file), sep="\t")
        # The TSV has columns "x" and "y" representing the spot's grid position
        # (row × col in the capture array).  We build the spot ID as "ROWxCOL"
        # so it matches the row labels in the count matrix (e.g. "10x13").
        spots.index = spots["x"].astype(int).astype(str) + "x" + spots["y"].astype(int).astype(str)
    else:
        rdata_file = next(
            (f for f in files if "allSpots" in f and f.endswith(".RData")), None
        )
        if rdata_file is None:
            raise RuntimeError(f"No spot coordinate file found in {array_dir}. Files: {files}")
        spots = list(pyreadr.read_r(os.path.join(array_dir, rdata_file)).values())[0]

    # Rename pixel coordinate columns to the standard names our pipeline expects:
    #   "pixel_x" or "x" or "col"  →  "xaxis"
    #   "pixel_y" or "y" or "row"  →  "yaxis"
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

    # The count matrix might have spots as columns and genes as rows (transposed).
    # We need spots as ROWS.  Try the original orientation first; if spot IDs don't
    # match, transpose and try again.
    if not counts.index.isin(coord.index).any():
        counts = counts.T
    if not counts.index.isin(coord.index).any():
        raise RuntimeError(
            f"Cannot align counts to spot coordinates — no shared IDs.\n"
            f"  counts index (first 3): {list(counts.index[:3])}\n"
            f"  coord  index (first 3): {list(coord.index[:3])}"
        )

    return img_array, counts, coord

# =============================================================================
# STEP 6 — TNBC Zenodo streaming downloader
# =============================================================================

def _iter_tnbc_zenodo(record_id):
    """Stream byArray.tar from Zenodo and yield one completed array at a time.

    WHY STREAMING?
    The byArray.tar file is very large (~50 GB).  If we downloaded it all to
    disk first we'd need 50 GB of free space.  Instead, we open the tar file
    directly over HTTP and read it entry by entry.  As soon as we have all three
    required files for one array (H&E image, counts, spot coords), we yield that
    array's data and discard it from memory.  Peak memory use: one array at a time.

    Yields: (array_id, file_dict)
      array_id  : e.g. "CN1_C2"
      file_dict : {"counts": (filename, bytes), "coords": ..., "he": ...}
    """
    import requests
    url = f"https://zenodo.org/records/{record_id}/files/byArray.tar"
    print(f"  Connecting to {url} …", flush=True)

    with requests.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        resp.raw.decode_content = True   # decompress gzip on the fly if needed

        buffers    = {}    # accumulates files for each array as we see them
        current    = None  # the array ID we're currently printing progress for
        bytes_read = 0

        # mode="r|" means "streaming tar" — we can only read forward, not seek.
        with tarfile.open(fileobj=resp.raw, mode="r|") as tar:
            for member in tar:
                if not member.isfile():
                    continue   # skip directory entries

                # The tar path looks like: byArray/CN_N/C2/filename
                # We combine the last two path components to form the array ID.
                parts = member.name.replace("\\", "/").split("/")
                if len(parts) < 4:
                    continue   # skip top-level files we don't need
                filename = parts[-1]
                array_id = f"{parts[-3]}_{parts[-2]}"   # e.g. "CN1_C2"

                # Check whether this file is one we care about.
                role = _tnbc_role(filename)
                if role is None:
                    continue

                if array_id != current:
                    current = array_id
                    print(f"    {array_id}  ({bytes_read/1e6:.0f} MB streamed)…",
                          end="", flush=True)

                # Read the file's bytes into memory.
                fobj = tar.extractfile(member)
                if fobj is None:
                    continue
                data = fobj.read()
                bytes_read += len(data)

                # Store under the array's buffer dict.
                buffers.setdefault(array_id, {})[role] = (filename, data)

                # Once we have all three required files for this array, yield it
                # and remove it from the buffer to free memory.
                if set(buffers[array_id].keys()) == _TNBC_REQUIRED_ROLES:
                    print(" ready", flush=True)
                    yield array_id, buffers.pop(array_id)

# =============================================================================
# STEP 7 — TIFF metadata fix
# =============================================================================

def _fix_tif_mpp(tif_path, pixel_size_um):
    """Fix a bug in HEST's TIFF writer where the resolution tag is 10× too small.

    BACKGROUND:
    TIFF files store image resolution as a "pixels per centimetre" fraction in
    their metadata tags (XResolution, YResolution).  HEST uses a library called
    pyvips to write the TIFF.  Due to a unit mismatch in pyvips, it writes the
    resolution as if the unit is mm even though the tag says cm — resulting in a
    value that is 10× smaller than the correct one.

    WHAT WE DO:
    After HEST saves the TIFF, we open it at the raw byte level with tifffile and
    overwrite just the resolution tags with the correct value.  We do NOT re-encode
    any pixel data (that would take a long time and degrade image quality).

    pixel_size_um : the true pixel size in micrometres (µm per pixel), e.g. 0.5
    """
    # Convert µm/px to the standard TIFF unit: pixels per centimetre.
    # 1 cm = 10,000 µm, so px/cm = 10,000 / (µm/px).
    px_per_cm = 10000.0 / pixel_size_um
    # TIFF stores rationals as (numerator, denominator).  We multiply by 1000
    # to preserve three decimal places of precision.
    rational = (round(px_per_cm * 1000), 1000)

    with open(tif_path, "r+b") as f:   # "r+b" = read + write, binary mode
        with tifffile.TiffFile(f) as tif:
            # A pyramidal TIFF has multiple "pages" (resolution levels).
            # We must fix the tag on every page so all levels are correct.
            for page in tif.pages:
                for tag_name in ("XResolution", "YResolution"):
                    if tag_name in page.tags:
                        page.tags[tag_name].overwrite(rational)

# =============================================================================
# STEP 8 — STimage downloader (HuggingFace)
# =============================================================================

def download_slide(slide_name, technology):
    """Download one STimage slide from the HuggingFace dataset repository.

    The STimage-1K4M dataset is hosted on HuggingFace at:
      jiawennnn/STimage-1K4M
    Each slide has three files stored under technology-named subfolders:
      {technology}/image/{slide_name}.png
      {technology}/coord/{slide_name}_coord.csv
      {technology}/gene_exp/{slide_name}_count.csv

    Returns a dict with paths to the downloaded files, or None on failure.
    """
    suffix      = TECH_SUFFIX.get(technology, technology.lower())
    # stimage_dir : where we cache the raw downloaded files
    stimage_dir = os.path.join(DATA_DIR, f"stimage_{suffix}", slide_name)
    # hest_dir    : where the final converted HEST output will go
    hest_dir    = os.path.join(DATA_DIR, f"hest_{suffix}", slide_name)
    os.makedirs(stimage_dir, exist_ok=True)
    os.makedirs(hest_dir,    exist_ok=True)

    for file_type, filename in {
        "coord":    f"{slide_name}_coord.csv",
        "gene_exp": f"{slide_name}_count.csv",
        "image":    f"{slide_name}.png",
    }.items():
        dest_path = os.path.join(stimage_dir, filename)
        if os.path.exists(dest_path):
            continue   # already downloaded on a previous run, skip
        try:
            # hf_hub_download() fetches a single file from a HuggingFace
            # repository and saves it to a local directory.
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
            return None   # signal failure to the caller

    return {"stimage_dir": stimage_dir, "hest_dir": hest_dir, "technology": technology}

# =============================================================================
# STEP 9 — Output folder name helper
# =============================================================================

def _hest_folder(source, data_path, technology):
    """Return the name of the top-level output folder for a given data source.

    We use different folder names per dataset so all outputs can coexist under
    the same DATA_DIR without overwriting each other.

    Examples:
      source="tnbc"         → "hest_tnbc"
      source="visium_geo"   → "hest_gse278936"
      source="visium_ebi"   → "hest_e_mtab_13530"
      source="visium_zenodo"→ "hest_visium_tls"
      source="stimage"      → "hest_visium" or "hest_st" (depends on technology)
    """
    if source in ("tnbc", "tnbc_zenodo"):
        return "hest_tnbc"
    if source == "visium_geo":
        # Use the GEO accession as the suffix, lowercased.
        return f"hest_{data_path.lower()}"
    if source == "visium_ebi":
        # Use the EBI accession, lowercased, with hyphens replaced by underscores.
        return f"hest_{data_path.lower().replace('-', '_')}"
    if source == "visium_zenodo":
        return "hest_visium_tls"
    # Default for STimage: use the technology name.
    return f"hest_{TECH_SUFFIX.get(technology, technology.lower())}"

# =============================================================================
# STEP 10 — Space Ranger output loader
# =============================================================================

def _load_spaceranger(slide_dir, img_path=None):
    """Load data from a 10x Space Ranger output directory.

    Space Ranger is 10x Genomics' official pipeline for processing Visium data.
    It produces a standardised folder structure that GEO, EBI, and the TLS
    dataset all use.  This function reads:
      - the gene expression count matrix (either .h5 or MTX format)
      - the spot pixel coordinates from tissue_positions[_list].csv
      - the spot diameter from scalefactors_json.json
      - the H&E image (from img_path, or a .tif in the folder, or the hires PNG)

    Returns: (img_array, counts, coord)
      coord has columns: xaxis (pixel x), yaxis (pixel y), r (spot radius in px)
    """
    # ── Count matrix ─────────────────────────────────────────────────────────
    h5_path = os.path.join(slide_dir, "filtered_feature_bc_matrix.h5")
    fbc_dir = os.path.join(slide_dir, "filtered_feature_bc_matrix")

    if os.path.exists(h5_path):
        # .h5 format: a single HDF5 file containing the full sparse matrix.
        ca = sc.read_10x_h5(h5_path)
    elif os.path.isdir(fbc_dir):
        # MTX format: three separate files (matrix.mtx.gz, barcodes.tsv.gz,
        # features.tsv.gz) in a subdirectory.  This is the older Space Ranger format.
        ca = sc.read_10x_mtx(fbc_dir, var_names="gene_symbols", cache=False)
    else:
        raise RuntimeError(f"No count matrix found in {slide_dir}")

    # Convert the sparse count matrix to a dense numpy array, then to a DataFrame.
    # Rows = spot barcodes (e.g. "ACGT...TTAG-1"), columns = gene names.
    X      = ca.X.toarray() if hasattr(ca.X, "toarray") else np.array(ca.X)
    counts = pd.DataFrame(X, index=ca.obs_names, columns=ca.var_names)

    # ── Spatial directory ─────────────────────────────────────────────────────
    # Space Ranger puts coordinate files in a "spatial/" subfolder.
    # If that subfolder doesn't exist, fall back to the slide directory itself.
    spatial_dir = os.path.join(slide_dir, "spatial")
    if not os.path.isdir(spatial_dir):
        spatial_dir = slide_dir

    # ── Tissue positions ──────────────────────────────────────────────────────
    # tissue_positions.csv (newer Space Ranger) or tissue_positions_list.csv (older)
    # contains one row per spot with columns:
    #   barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
    # "in_tissue" = 1 means the spot overlaps tissue (as determined by Space Ranger).
    pos_file = None
    for name in ("tissue_positions.csv", "tissue_positions_list.csv"):
        p = os.path.join(spatial_dir, name)
        if os.path.exists(p):
            pos_file = p; break

    # If not found in the expected location, search the whole slide directory.
    if pos_file is None:
        for root, _, fns in os.walk(slide_dir):
            for fn in fns:
                if fn in ("tissue_positions.csv", "tissue_positions_list.csv"):
                    pos_file = os.path.join(root, fn); break
            if pos_file:
                break
    if pos_file is None:
        raise RuntimeError(f"No tissue_positions CSV found under {slide_dir}")

    # Read the CSV.  Newer versions have a header row; older ones don't.
    raw = pd.read_csv(pos_file, header=None, dtype=str)
    if raw.iloc[0, 0].strip().lower() == "barcode":
        # Header row present — use it as column names and drop the first data row.
        raw.columns = [c.strip() for c in raw.iloc[0]]
        raw = raw.iloc[1:].reset_index(drop=True)
    else:
        # No header — assign column names manually.
        raw.columns = ["barcode", "in_tissue", "array_row", "array_col",
                       "pxl_row_in_fullres", "pxl_col_in_fullres"]

    # Keep only spots that are under tissue (in_tissue == "1").
    raw = raw[raw["in_tissue"].str.strip() == "1"].copy()
    raw.index = raw["barcode"].str.strip()   # use barcode as the row index

    # ── Scale factors ─────────────────────────────────────────────────────────
    # scalefactors_json.json stores numbers that relate the full-resolution image
    # to the downscaled images Space Ranger outputs.
    # We need "spot_diameter_fullres" — the diameter of a Visium spot in pixels,
    # measured in the full-resolution image.  Half of that is the spot radius (r).
    sf_file = os.path.join(spatial_dir, "scalefactors_json.json")
    if not os.path.exists(sf_file):
        for root, _, fns in os.walk(slide_dir):
            if "scalefactors_json.json" in fns:
                sf_file = os.path.join(root, "scalefactors_json.json"); break
    with open(sf_file) as fj:
        sf = json.load(fj)
    # spot radius in pixels at full resolution
    r_fullres = float(sf["spot_diameter_fullres"]) / 2.0

    # ── Image ─────────────────────────────────────────────────────────────────
    # Priority order for finding the image:
    #   1. The img_path argument (caller explicitly supplied a path)
    #   2. Any .tif or .tiff file in the slide directory
    #   3. The tissue_hires_image.png that Space Ranger generates (lower resolution)
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

    # Load the image into a numpy array of shape (height, width, 3).
    img_array = np.array(Image.open(img_path).convert("RGB"))

    # ── Align barcodes ────────────────────────────────────────────────────────
    # The count matrix and the tissue positions file might not have exactly the
    # same set of barcodes (Space Ranger sometimes has minor mismatches).
    # We keep only the barcodes present in both.
    common = counts.index.intersection(raw.index)
    counts = counts.loc[common]
    raw    = raw.loc[common]

    # Build the coordinate DataFrame.
    # If we're using the hires PNG instead of the full-res image, the pixel
    # coordinates in tissue_positions.csv refer to the FULL-res image and need
    # to be scaled down by tissue_hires_scalef (typically 0.2 = 5× downscale).
    if using_hires:
        hires_scalef = float(sf.get("tissue_hires_scalef", 1.0))
        coord = pd.DataFrame({
            "xaxis": raw["pxl_col_in_fullres"].astype(float) * hires_scalef,
            "yaxis": raw["pxl_row_in_fullres"].astype(float) * hires_scalef,
            "r":     r_fullres * hires_scalef,
        }, index=raw.index)
    else:
        # Full-res image: coordinates are already in the correct pixel space.
        coord = pd.DataFrame({
            "xaxis": raw["pxl_col_in_fullres"].astype(float),
            "yaxis": raw["pxl_row_in_fullres"].astype(float),
            "r":     r_fullres,   # spot radius in full-res pixels
        }, index=raw.index)

    return img_array, counts, coord

# =============================================================================
# STEP 11 — GEO downloader
# =============================================================================

def _download_geo_raw(geo_acc, cache_dir):
    """Download the _RAW.tar archive for a GEO series from the NCBI FTP server.

    NCBI packages all supplementary files for a GEO series into one tar archive
    named {accession}_RAW.tar and stores it on their FTP server at:
      ftp.ncbi.nlm.nih.gov/geo/series/GSE{nnn}/{acc}/suppl/{acc}_RAW.tar

    The {nnn} in the path is the accession number with the last 3 digits replaced
    by "nnn" — a NCBI convention for grouping accessions in directories.
    E.g. GSE278936 → GSE278nnn

    The file is cached locally so repeated runs don't re-download it.
    Returns the local path to the downloaded tar file.
    """
    try:
        import requests
    except ImportError:
        raise RuntimeError("pip install requests")

    os.makedirs(cache_dir, exist_ok=True)
    tar_path = os.path.join(cache_dir, f"{geo_acc}_RAW.tar")

    if os.path.exists(tar_path):
        print(f"  Using cached {geo_acc}_RAW.tar")
        return tar_path

    # Build the NCBI FTP URL.  Strip "GSE" prefix to get the numeric part.
    n   = geo_acc[3:]
    # Replace the last 3 digits with "nnn" (NCBI directory grouping convention).
    nnn = (n[:-3] + "nnn") if len(n) > 3 else (n + "nnn")
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{nnn}"
           f"/{geo_acc}/suppl/{geo_acc}_RAW.tar")

    print(f"  Downloading {url} …")
    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(tar_path, "wb") as fh:
            mb = 0
            # Download in 8 MB chunks and print running total.
            for chunk in r.iter_content(chunk_size=8_388_608):
                fh.write(chunk); mb += len(chunk) / 1e6
                print(f"\r  {mb:.0f} MB…", end="", flush=True)
    print()
    return tar_path

# =============================================================================
# STEP 12 — GEO extractor
# =============================================================================

def _extract_all_geo_gsm(tar_path, gsm_out_dirs):
    """Extract all GSM samples from one GEO RAW.tar in a single pass.

    WHY SINGLE PASS?
    The RAW.tar can be several gigabytes.  If we looped over it once per sample
    we'd read the whole file N times.  Instead, we read it exactly once and
    route each file to the correct sample directory based on the GSM prefix.

    GEO FILE NAMING CONVENTION:
    Files inside the tar are named: {GSM_ID}_{SAMPLE_LABEL}_{filetype}[.gz]
    e.g.:  GSM8557976_BPH_1_matrix.mtx.gz
           GSM8557976_BPH_1_scalefactors_json.json.gz

    Note: GEO sometimes adds an extra .gz to spatial files that weren't
    originally gzip-compressed (e.g. scalefactors_json.json becomes
    scalefactors_json.json.gz).  We detect and decompress these on the fly.

    _RULES maps filename suffixes to destination paths.  Each rule is:
      (suffix_to_match, destination_subfolder, output_filename, decompress_gz)
    """
    import gzip

    # Each rule: if the filename ends with the given suffix, copy the file to
    # the given subfolder with the given output name.  If decompress=True,
    # run gzip decompression before writing.
    _RULES = [
        # Count matrix files (MTX format) — keep gzip compression as-is
        ("_matrix.mtx.gz",               "filtered_feature_bc_matrix", "matrix.mtx.gz",             False),
        ("_barcodes.tsv.gz",              "filtered_feature_bc_matrix", "barcodes.tsv.gz",           False),
        ("_features.tsv.gz",              "filtered_feature_bc_matrix", "features.tsv.gz",           False),
        ("_genes.tsv.gz",                 "filtered_feature_bc_matrix", "genes.tsv.gz",              False),
        ("_matrix.mtx",                   "filtered_feature_bc_matrix", "matrix.mtx",                False),
        ("_barcodes.tsv",                 "filtered_feature_bc_matrix", "barcodes.tsv",              False),
        ("_features.tsv",                 "filtered_feature_bc_matrix", "features.tsv",              False),
        ("_genes.tsv",                    "filtered_feature_bc_matrix", "genes.tsv",                 False),
        # Spatial files — GEO sometimes double-gzips these, so decompress them
        ("_scalefactors_json.json.gz",    "spatial", "scalefactors_json.json",     True),
        ("_scalefactors_json.json",       "spatial", "scalefactors_json.json",     False),
        ("_tissue_positions_list.csv.gz", "spatial", "tissue_positions_list.csv",  True),
        ("_tissue_positions_list.csv",    "spatial", "tissue_positions_list.csv",  False),
        ("_tissue_positions.csv.gz",      "spatial", "tissue_positions.csv",       True),
        ("_tissue_positions.csv",         "spatial", "tissue_positions.csv",       False),
        ("_tissue_hires_image.png.gz",    "spatial", "tissue_hires_image.png",     True),
        ("_tissue_hires_image.png",       "spatial", "tissue_hires_image.png",     False),
        ("_tissue_lowres_image.png.gz",   "spatial", "tissue_lowres_image.png",    True),
        ("_tissue_lowres_image.png",      "spatial", "tissue_lowres_image.png",    False),
    ]

    # Create destination subdirectories for every sample in advance.
    for out_dir in gsm_out_dirs.values():
        os.makedirs(os.path.join(out_dir, "filtered_feature_bc_matrix"), exist_ok=True)
        os.makedirs(os.path.join(out_dir, "spatial"), exist_ok=True)

    # Build a lookup: "GSM8557976_" → ("GSM8557976", "/path/to/output/GSM8557976")
    # We use the prefix (GSM ID + underscore) to quickly identify which sample
    # each file belongs to.
    prefix_map = {gsm_id + "_": (gsm_id, out_dir)
                  for gsm_id, out_dir in gsm_out_dirs.items()}

    print(f"  Extracting {len(gsm_out_dirs)} samples from {os.path.basename(tar_path)} …",
          flush=True)

    with tarfile.open(tar_path, "r") as tar:
        for member in tar.getmembers():
            fname = os.path.basename(member.name)

            # Identify which GSM sample this file belongs to.
            matched = next(((gsm_id, od) for pfx, (gsm_id, od) in prefix_map.items()
                            if fname.startswith(pfx)), None)
            if matched is None:
                continue   # file doesn't belong to any of our samples
            gsm_id, out_dir = matched

            # Find the matching rule for this filename's suffix.
            rule = next(((sub, bn, dec) for sfx, sub, bn, dec in _RULES
                         if fname.endswith(sfx)), None)
            if rule is None:
                continue   # unrecognised file type — skip
            subdir, dest_basename, decompress = rule

            dest = os.path.join(out_dir, subdir, dest_basename)
            if os.path.exists(dest):
                continue   # already extracted on a previous run

            fobj = tar.extractfile(member)
            if fobj is None:
                continue
            data = fobj.read()

            # Decompress if GEO added an extra .gz wrapper.
            if decompress:
                data = gzip.decompress(data)

            with open(dest, "wb") as fh:
                fh.write(data)

# =============================================================================
# STEP 13 — TLS zip downloader (with resume + retry)
# =============================================================================

def _download_tls_zip(record_url, cache_dir):
    """Download TLS_VISIUM_USZ.zip from Zenodo, with automatic resume and retry.

    This is a ~2 GB file.  On slow or unreliable connections (common on HPC
    clusters with shared internet) the download may be interrupted.

    HOW RESUME WORKS:
    HTTP supports a "Range" request header that tells the server to start
    sending from a specific byte offset.  If we already have a partial file,
    we tell the server "start from byte N" and append to the existing file.

    HOW VALIDATION WORKS:
    A ZIP file has an "end of central directory" record at the very end.
    If the download was cut short, the ZIP is missing this record.
    zipfile.ZipFile() raises BadZipFile if it's missing — we catch that to
    detect partial downloads and resume from where we left off.

    HOW RETRY WORKS:
    If the connection drops mid-download, we catch the exception, wait, and
    loop back to try again.  The partial file is kept on disk so the next
    attempt resumes (not restarts).  We try up to 10 times before giving up.
    """
    try:
        import requests
    except ImportError:
        raise RuntimeError("pip install requests")
    import zipfile as _zf

    os.makedirs(cache_dir, exist_ok=True)
    zip_path  = os.path.join(cache_dir, "TLS_VISIUM_USZ.zip")
    record_id = record_url.rstrip("/").split("/")[-1]
    url       = f"https://zenodo.org/records/{record_id}/files/TLS_VISIUM_USZ.zip"

    # Fast path: if we already have a complete, valid ZIP, use it immediately.
    if os.path.exists(zip_path):
        try:
            with _zf.ZipFile(zip_path): pass   # just opens and validates
            print("  Using cached TLS_VISIUM_USZ.zip")
            return zip_path
        except _zf.BadZipFile:
            pass   # file is incomplete — fall through to resume logic below

    MAX_RETRIES = 10
    for attempt in range(1, MAX_RETRIES + 1):
        # How many bytes do we already have on disk?
        existing = os.path.getsize(zip_path) if os.path.exists(zip_path) else 0
        # Tell the server to start from that byte offset (resume).
        headers  = {"Range": f"bytes={existing}-"} if existing else {}
        # "ab" = append bytes (resume); "wb" = write bytes from scratch.
        mode     = "ab" if existing else "wb"

        if existing:
            print(f"  Resuming TLS_VISIUM_USZ.zip from {existing/1e6:.0f} MB "
                  f"(attempt {attempt}/{MAX_RETRIES})…")
        else:
            print(f"  Downloading {url} (~2 GB) …")

        try:
            with requests.get(url, stream=True, headers=headers, timeout=600) as r:
                r.raise_for_status()
                # Content-Length tells us how many bytes remain to download.
                # Adding "existing" gives the total expected file size.
                total = existing + int(r.headers.get("Content-Length", 0))
                with open(zip_path, mode) as fh:
                    mb = existing / 1e6
                    for chunk in r.iter_content(chunk_size=8_388_608):
                        fh.write(chunk); mb += len(chunk) / 1e6
                        pct = f" ({mb/total*100:.0f}%)" if total else ""
                        print(f"\r  {mb:.0f} MB{pct}…", end="", flush=True)
            print()

            # Verify the download completed correctly.
            try:
                with _zf.ZipFile(zip_path): pass
                return zip_path   # success
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

# =============================================================================
# STEP 14 — TLS sample loader
# =============================================================================

def _load_tls_sample(sample_name, zip_path):
    """Extract one TLS sample from the ZIP and load it via _load_spaceranger.

    The TLS ZIP file (TLS_VISIUM_USZ.zip) contains all 8 samples in one archive
    with this internal folder structure:
      TLS_VISIUM_USZ/
        10x_Visium/{sample}/          ← full Space Ranger output
          filtered_feature_bc_matrix/ ← MTX count matrix
          spatial/                    ← tissue_positions.csv, scalefactors, images
        tif_slides/{sample}.tif       ← full-resolution H&E image

    Because _load_spaceranger() needs files on disk (not inside a ZIP), we:
      1. Extract the Space Ranger folder to a temporary directory
      2. Extract the TIF to a temporary file
      3. Call _load_spaceranger() on those temp files
      4. Delete the temp files before returning (to save disk space)
    """
    sr_prefix = f"TLS_VISIUM_USZ/10x_Visium/{sample_name}/"   # Space Ranger folder
    tif_entry = f"TLS_VISIUM_USZ/tif_slides/{sample_name}.tif" # full-res image

    tmp_dir = tempfile.mkdtemp()   # create a unique temporary directory
    tif_tmp = None
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            # Extract all Space Ranger files for this sample, stripping the
            # ZIP-internal prefix so the result looks like a plain Space Ranger folder.
            for member in zf.infolist():
                name = member.filename
                if name.startswith(sr_prefix) and not name.endswith("/"):
                    rel  = name[len(sr_prefix):]   # path relative to this sample
                    dest = os.path.join(tmp_dir, rel)
                    os.makedirs(os.path.dirname(dest), exist_ok=True)
                    with zf.open(member) as src, open(dest, "wb") as dst:
                        dst.write(src.read())

            # Extract the TIF to a named temp file (tifffile needs a real path).
            with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp:
                tmp.write(zf.read(tif_entry))
                tif_tmp = tmp.name

        # Load the data using the same Space Ranger loader as GEO and EBI.
        return _load_spaceranger(tmp_dir, img_path=tif_tmp)

    finally:
        # Always clean up temp files, even if an exception was raised.
        shutil.rmtree(tmp_dir, ignore_errors=True)
        if tif_tmp:
            try: os.unlink(tif_tmp)
            except: pass

# =============================================================================
# STEP 15 — EBI per-sample downloader
# =============================================================================

def _download_ebi_sample(sample_name, ebi_acc, cache_dir):
    """Download and unpack one sample from EBI BioStudies / ArrayExpress.

    EBI stores each sample as exactly two files at a flat URL path:
      {base_url}/{sample_name}-filtered_feature_bc_matrix.h5   ← count matrix
      {base_url}/{sample_name}-spatial.tar                      ← spatial files

    The spatial.tar contains: scalefactors_json.json,
    tissue_positions.csv, tissue_hires_image.png, tissue_lowres_image.png.

    We download both files, extract the tar into a "spatial/" subfolder, and
    return the local directory path (which _load_spaceranger can then read).
    Files are cached — if already downloaded, we skip the download.
    """
    try:
        import requests
    except ImportError:
        raise RuntimeError("pip install requests")

    sample_dir  = os.path.join(cache_dir, "ebi", ebi_acc, sample_name)
    h5_dest     = os.path.join(sample_dir, "filtered_feature_bc_matrix.h5")
    spatial_dir = os.path.join(sample_dir, "spatial")
    base_url    = f"https://www.ebi.ac.uk/biostudies/files/{ebi_acc}"

    # Fast path: both files already on disk.
    if os.path.exists(h5_dest) and os.path.isdir(spatial_dir):
        print(f"    Using cached {sample_name}/")
        return sample_dir

    os.makedirs(sample_dir,  exist_ok=True)
    os.makedirs(spatial_dir, exist_ok=True)

    def _fetch(remote_name, dest_path, label):
        """Download one file from EBI if it isn't already on disk."""
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

    # Download the count matrix HDF5 file.
    _fetch(f"{sample_name}-filtered_feature_bc_matrix.h5", h5_dest, "h5")

    # Download the spatial tar and extract it.
    tar_tmp = os.path.join(sample_dir, "_spatial.tar")
    _fetch(f"{sample_name}-spatial.tar", tar_tmp, "spatial.tar")
    if os.path.exists(tar_tmp):
        with tarfile.open(tar_tmp, "r:*") as tar:
            for member in tar.getmembers():
                fname = os.path.basename(member.name)
                if not fname or member.isdir():
                    continue   # skip directory entries
                dest = os.path.join(spatial_dir, fname)
                if os.path.exists(dest):
                    continue
                fobj = tar.extractfile(member)
                if fobj:
                    with open(dest, "wb") as fh:
                        fh.write(fobj.read())
        os.unlink(tar_tmp)   # delete the tar after extracting (saves disk space)

    return sample_dir

# =============================================================================
# STEP 16 — Core conversion function
# =============================================================================

def convert_slide(slide_name, hest_dir, technology, row, img_array, counts, coord):
    """Convert one tissue slide into HEST format and write all output files.

    This is the heart of the pipeline.  It takes the raw data (image + counts
    + coordinates) and produces the full standardised HEST output folder.

    Arguments:
      slide_name  : unique ID for this slide (e.g. "GSM8557976" or "KC1")
      hest_dir    : output directory to write into
      technology  : "ST" or "Visium"
      row         : the CSV row dict with metadata (tissue, species, pmid, etc.)
      img_array   : numpy array (H×W×3) — the H&E image
      counts      : DataFrame (spots × genes) — raw RNA counts
      coord       : DataFrame with xaxis, yaxis[, r] columns — spot pixel positions
    """

    # ── 1. Align counts to coordinates ───────────────────────────────────────
    # Keep only spots that appear in BOTH the count matrix and coordinate table.
    # Discards any spot that has counts but no recorded position, or vice versa.
    common        = counts.index.intersection(coord.index)
    counts, coord = counts.loc[common], coord.loc[common]

    # ── 2. Shorten spot IDs ───────────────────────────────────────────────────
    # ST slides: spot IDs look like "SAMPLE_10x13" — we strip the prefix to get
    # just the grid position "10x13".
    # Visium slides: spot IDs look like "GSM8557976_ACGT...TTAG-1" — we strip
    # the slide-name prefix to get just the barcode "ACGT...TTAG-1".
    if technology == "ST":
        try:
            short_ids = [re.search(r"(\d+x\d+)$", idx).group(1) for idx in counts.index]
        except AttributeError:
            short_ids = list(counts.index.astype(str))
    else:
        prefix    = slide_name + "_"
        short_ids = [idx[len(prefix):] if idx.startswith(prefix) else idx for idx in counts.index]
    counts.index = coord.index = short_ids

    # ── 3. Create the AnnData object ──────────────────────────────────────────
    # AnnData (Annotated Data) is the standard container used by scanpy and HEST.
    # It holds the count matrix plus any per-spot and per-gene metadata.
    adata = sc.AnnData(counts)

    # Store pixel coordinates as a 2D array in obsm["spatial"].
    # HEST and scanpy/squidpy both look for spatial coordinates here.
    # Column 0 = x (horizontal), column 1 = y (vertical).
    adata.obsm["spatial"] = coord[["xaxis", "yaxis"]].values

    # ── 4. Build the obs (per-spot) metadata table ───────────────────────────
    my_df = pd.DataFrame(
        adata.obsm["spatial"], index=adata.obs_names,
        columns=["pxl_col_in_fullres", "pxl_row_in_fullres"]
    )

    # Compute array_row and array_col — the integer grid position of each spot.
    # For ST slides, the grid position is encoded in the spot ID ("ROWxCOL").
    # For Visium, we rank-order the unique pixel positions to recover the grid.
    if technology == "ST":
        try:
            my_df["array_row"] = [round(float(i.split("x")[0])) for i in my_df.index]
            my_df["array_col"] = [round(float(i.split("x")[1])) for i in my_df.index]
        except (ValueError, IndexError):
            # Fallback: derive grid indices from pixel positions.
            _, y_inv = np.unique(my_df["pxl_row_in_fullres"].round().astype(int).values, return_inverse=True)
            _, x_inv = np.unique(my_df["pxl_col_in_fullres"].round().astype(int).values, return_inverse=True)
            my_df["array_row"] = y_inv
            my_df["array_col"] = x_inv
    else:
        _, y_inv = np.unique(my_df["pxl_row_in_fullres"].round().astype(int).values, return_inverse=True)
        _, x_inv = np.unique(my_df["pxl_col_in_fullres"].round().astype(int).values, return_inverse=True)
        my_df["array_row"] = y_inv
        my_df["array_col"] = x_inv

    # Copy all per-spot metadata into adata.obs.
    adata.obs["array_row"]          = my_df["array_row"].values
    adata.obs["array_col"]          = my_df["array_col"].values
    adata.obs["pxl_col_in_fullres"] = my_df["pxl_col_in_fullres"].values
    adata.obs["pxl_row_in_fullres"] = my_df["pxl_row_in_fullres"].values
    adata.obs["in_tissue"]          = True   # all spots in our data are on-tissue

    # ── 5. Calculate pixel size (µm per pixel) ────────────────────────────────
    # This is the most critical calibration step.  We need to know how large
    # each pixel is in real physical units (micrometres) so that:
    #   a) The TIFF file metadata is correct (needed for pathology software)
    #   b) We can extract patches at the correct physical magnification
    #
    # METHOD A (most accurate): use the spot radius from Space Ranger.
    # Space Ranger reports spot_diameter_fullres in its scalefactors_json.json.
    # Since we know the true physical diameter of a Visium spot is exactly 55 µm,
    # pixel_size = (55 µm / 2) / r_pixels = 27.5 / r
    #
    # METHOD B (fallback): estimate from the spacing between spots.
    # We know adjacent Visium spots are 100 µm apart.  By measuring how many
    # pixels apart the nearest neighbours are, we can compute µm/px.
    if technology == "ST":
        pixel_size, spot_estimate_dist = find_pixel_size_from_spot_coords(
            my_df, inter_spot_dist=SPOT_PARAMS[technology]["inter_spot_dist"],
            packing=SpotPacking.GRID_PACKING   # ST uses a square grid
        )
    else:
        if "r" in coord.columns:
            # Method A: spot radius is known from Space Ranger output.
            # spot_diameter (µm) / 2 / r_pixels = µm per pixel
            pixel_size = SPOT_PARAMS[technology]["spot_diameter"] / 2.0 / coord["r"].median()
        else:
            # Method B: estimate from spot spacing (less accurate but always works).
            pixel_size, _ = find_pixel_size_from_spot_coords(
                my_df, inter_spot_dist=SPOT_PARAMS[technology]["inter_spot_dist"],
                packing=SpotPacking.GRID_PACKING
            )
        # Convert inter-spot distance from µm to pixels.
        spot_estimate_dist = SPOT_PARAMS[technology]["inter_spot_dist"] / pixel_size

    # ── 6. Upsample very low-resolution images ───────────────────────────────
    # STimage Visium slides are stored as ~2000-pixel PNG thumbnails, giving a
    # resolution of roughly 3–9 µm/px.  HEST's internal library (wsi_factory)
    # rejects images with pixel_size >= 2.4 µm/px.
    # Solution: enlarge the image so that 1 pixel = 1 µm, then rescale the
    # pixel coordinates accordingly.
    if pixel_size >= 2.0:
        scale = pixel_size / 1.0   # how much to enlarge (e.g. 4.0 → enlarge 4×)
        h, w = img_array.shape[:2]
        # Resize the image using LANCZOS (high-quality downsampling algorithm).
        img_array = np.array(Image.fromarray(img_array).resize(
            (round(w * scale), round(h * scale)), Image.LANCZOS
        ))
        # Rescale the pixel coordinates to match the larger image.
        my_df["pxl_col_in_fullres"] *= scale
        my_df["pxl_row_in_fullres"] *= scale
        adata.obs["pxl_col_in_fullres"] = my_df["pxl_col_in_fullres"].values
        adata.obs["pxl_row_in_fullres"] = my_df["pxl_row_in_fullres"].values
        adata.obsm["spatial"]           = my_df[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values
        # After enlarging, each pixel is now 1 µm instead of the original size.
        pixel_size       /= scale
        spot_estimate_dist *= scale

    # ── 7. Register image with HEST ───────────────────────────────────────────
    # wsi_factory() wraps the numpy image array in HEST's WSI object.
    # register_downscale_img() creates a downscaled thumbnail and embeds it
    # inside adata.uns["spatial"]["ST"]["images"]["downscaled_fullres"] —
    # this is required by the HEST format specification.
    wsi = wsi_factory(img_array, mpp=pixel_size)   # mpp = microns per pixel
    register_downscale_img(adata, wsi, pixel_size,
                           spot_size=SPOT_PARAMS[technology]["spot_diameter"])

    # ── 8. Pad ST spot indices to zero-padded strings ─────────────────────────
    # For ST slides, normalise index format from "10x13" to "010x013" so that
    # alphabetical sorting matches spatial ordering (10x2 < 10x13 when zero-padded).
    if technology == "ST" and all(re.match(r"^\d+x\d+$", i) for i in list(adata.obs.index)[:5]):
        adata.obs.index = [
            i.split("x")[0].zfill(3) + "x" + i.split("x")[1].zfill(3)
            for i in adata.obs.index
        ]

    # ── 9. Calculate quality control (QC) metrics ────────────────────────────
    # Mark mitochondrial genes.  In human/mouse, mitochondrial genes start with
    # "MT-".  High mitochondrial percentage often indicates a dying or damaged cell.
    adata.var["mito"] = adata.var_names.str.startswith("MT-")

    # scanpy's calculate_qc_metrics() computes:
    #   obs columns: n_genes_by_counts, total_counts, pct_counts_mito,
    #                log1p_* versions, pct_counts_in_top_N_genes
    #   var columns: n_cells_by_counts, mean_counts, pct_dropout_by_counts,
    #                total_counts, log1p_* versions
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)

    # HEST expects "n_counts" as an alias for "total_counts".
    adata.obs["n_counts"] = adata.obs["total_counts"]

    # ── 10. Build the slide-level metadata dictionary ─────────────────────────
    # HEST stores per-slide metadata as a JSON-like dict inside the STHESTData
    # object.  We populate it from both the computed values above and the row
    # from cleaned_metadata.csv.

    # If the slide name looks like "GSE278936_GSM8557976", extract both parts.
    gse_m    = re.match(r"(GSE\d+)_(GSM\d+)", slide_name)
    gene_num = row.get("gene_num", "").strip()

    meta = {
        # --- Computed from the image and spot layout ---
        "pixel_size_um_embedded":  None,          # not embedded in original source
        "pixel_size_um_estimated": pixel_size,    # our calibration result
        "fullres_height":          img_array.shape[0],
        "fullres_width":           img_array.shape[1],
        "spots_under_tissue":      len(adata.obs),
        "spot_estimate_dist":      int(spot_estimate_dist),
        "spot_diameter":           SPOT_PARAMS[technology]["spot_diameter"],
        "inter_spot_dist":         SPOT_PARAMS[technology]["inter_spot_dist"],
        # --- From the metadata CSV (originally from the paper) ---
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
        # For GEO slides, link back to the series page on NCBI.
        "download_page_link1":     f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_m.group(1)}" if gse_m else None,
        "subseries":               gse_m.group(2) if gse_m else slide_name,
        # Fields not available in our sources (would need manual curation).
        "patient":                 None,
        "oncotree_code":           None,
        "data_publication_date":   None,
        "license":                 None,
        "preservation_method":     None,
        "magnification":           None,
        "treatment_comment":       None,
        "disease_comment":         None,
    }

    # ── 11. Save to HEST format ───────────────────────────────────────────────
    # STHESTData wraps our data + metadata in HEST's container class.
    st = STHESTData(adata, img_array, pixel_size, meta)

    # save() writes:
    #   aligned_adata.h5ad      — the scanpy data object
    #   aligned_fullres_HE.tif  — the H&E image as a pyramidal TIFF
    #                             (pyramidal = multiple resolutions in one file,
    #                              like a Google Maps tile pyramid, for fast
    #                              viewing at any zoom level)
    #   downscaled_fullres.jpeg — a small thumbnail of the full slide
    #   metrics.json            — pixel size and spot count metadata
    st.save(path=hest_dir, save_img=True, pyramidal=True, bigtiff=False, plot_pxl_size=True)

    # Fix the pixel size metadata that HEST's TIFF writer gets wrong (see Step 7).
    _fix_tif_mpp(os.path.join(hest_dir, "aligned_fullres_HE.tif"), pixel_size)

    # spatial_plot: saves spatial_plots.png — the H&E image with coloured dots
    # showing where each sequenced spot sits.
    st.save_spatial_plot(save_path=hest_dir)

    # Tissue segmentation: classify every pixel as "tissue" or "background"
    # using Otsu thresholding (a standard image processing algorithm that
    # automatically finds the brightness cut-off between tissue and white space).
    st.segment_tissue(method="otsu")

    # Save the segmentation as:
    #   tissue_seg_contours.geojson — polygon outlines of tissue regions
    #   tissue_seg_vis.jpg          — visualisation overlaid on the downscaled image
    st.save_tissue_contours(hest_dir, "tissue_seg")
    st.save_tissue_vis(hest_dir, "tissue_seg")

    # Extract 256×256 pixel patches centred on each spot, at three magnifications.
    # See PATCH_SCALES at the top for explanation of the three resolutions.
    # dump_visualization=True also saves a *_patch_vis.png overview image.
    for target_pixel_size, name in PATCH_SCALES:
        st.dump_patches(patch_save_dir=hest_dir, name=name,
                        target_patch_size=256, target_pixel_size=target_pixel_size,
                        dump_visualization=True)

# =============================================================================
# STEP 17 — Main loop: read the CSV and dispatch each slide
# =============================================================================

# Read all rows from the cleaned_metadata.csv that Clean.py produced.
with open(CLEANED_META, newline="") as f:
    slides = list(csv.DictReader(f))

# =============================================================================
# STEP 18 — Group slides by download strategy
# =============================================================================
# Different datasets require fundamentally different download strategies.
# We separate them up front so each group can use the optimal approach:
#
#   zenodo_groups        : TNBC from Zenodo → stream one big tar, never save it
#   geo_groups           : GEO datasets → download one RAW.tar, extract all samples
#   visium_zenodo_groups : TLS from Zenodo → download one ZIP, load samples from it
#   other_slides         : STimage (HuggingFace) + TNBC local + EBI → one at a time

zenodo_groups        = {}   # Zenodo URL → {array_id: row}
geo_groups           = {}   # GEO accession → {GSM_id: row}
visium_zenodo_groups = {}   # Zenodo URL → {sample_name: row}
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

# Counters for the final summary line.
failed  = []
skipped = 0
done    = 0

# =============================================================================
# STEP 19 — Process non-streaming slides (STimage, TNBC local, EBI)
# =============================================================================
# These are handled one at a time.  Each slide's data is downloaded, loaded,
# converted, and written before moving to the next.

for i, row in enumerate(other_slides):
    slide_name = row["slide"]
    technology = row["tech"]
    source     = row.get("source", "stimage")
    folder     = _hest_folder(source, row.get("data_path", ""), technology)
    hest_dir   = os.path.join(DATA_DIR, folder, slide_name)

    # If this slide was already converted in a previous run, skip it.
    # This means interrupted runs can be safely restarted without re-doing work.
    if os.path.exists(os.path.join(hest_dir, "aligned_adata.h5ad")):
        skipped += 1
        continue

    print(f"[{i+1}/{len(other_slides)}] {slide_name} ({technology}) [{source}]")
    try:
        if source == "stimage":
            # Download from HuggingFace, then load the three STimage files.
            result = download_slide(slide_name, technology)
            if result is None:
                raise RuntimeError("download failed")
            img_array, counts, coord = _load_stimage_data(slide_name, result["stimage_dir"])

        elif source == "tnbc":
            # Data is already on disk locally — just read it.
            os.makedirs(hest_dir, exist_ok=True)
            img_array, counts, coord = _load_tnbc_data(slide_name, row["data_path"])

        elif source == "visium_ebi":
            # Download the two EBI files for this sample, then load via Space Ranger.
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

# =============================================================================
# STEP 20 — Process TNBC Zenodo (streaming)
# =============================================================================
# We stream byArray.tar from Zenodo without saving it to disk.
# The generator _iter_tnbc_zenodo() yields one array at a time as it streams.

for data_path, array_rows in zenodo_groups.items():
    record_id  = data_path.rstrip("/").split("/")[-1]
    sample_row = next(iter(array_rows.values()))
    technology = sample_row["tech"]
    folder     = "hest_tnbc"

    # Count how many arrays still need converting (not already done).
    n_todo = sum(
        1 for aid in array_rows
        if not os.path.exists(os.path.join(DATA_DIR, folder, aid, "aligned_adata.h5ad"))
    )
    print(f"\n[tnbc_zenodo] record={record_id}  arrays={len(array_rows)}  to_do={n_todo}")
    if n_todo == 0:
        skipped += len(array_rows); continue

    for array_id, file_bufs in _iter_tnbc_zenodo(record_id):
        if array_id not in array_rows:
            continue   # this array isn't in our CSV (shouldn't happen, but safe)
        row      = array_rows[array_id]
        hest_dir = os.path.join(DATA_DIR, folder, array_id)
        if os.path.exists(os.path.join(hest_dir, "aligned_adata.h5ad")):
            skipped += 1; continue

        print(f"  {array_id}")
        try:
            # Write the streamed files to a temporary directory, load them,
            # then delete the temp directory immediately after loading.
            with tempfile.TemporaryDirectory(dir=DATA_DIR) as tmpdir:
                arr_dir = os.path.join(tmpdir, "byArray", array_id)
                os.makedirs(arr_dir)
                for role, (fname, data) in file_bufs.items():
                    with open(os.path.join(arr_dir, fname), "wb") as fh:
                        fh.write(data)
                img_array, counts, coord = _load_tnbc_data(array_id, tmpdir)
            convert_slide(array_id, hest_dir, technology, row, img_array, counts, coord)
            done += 1; print(f"    done")
        except Exception as e:
            print(f"    ERROR: {e}")
            failed.append({"slide": array_id, "technology": technology, "error": str(e)})

# =============================================================================
# STEP 21 — Process GEO datasets (batch tar extraction)
# =============================================================================
# Download the whole RAW.tar once, extract all samples in a single pass,
# then convert each sample.

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
        skipped += len(gsm_rows); continue

    try:
        # Download the RAW.tar (or use the cached copy).
        tar_path = _download_geo_raw(geo_acc, geo_cache)
        # Create one output directory per sample.
        gsm_out_dirs = {gsm: os.path.join(geo_cache, geo_acc, gsm)
                        for gsm in gsm_rows}
        # Extract all samples from the tar in one pass.
        _extract_all_geo_gsm(tar_path, gsm_out_dirs)
    except Exception as e:
        print(f"  ERROR downloading/extracting {geo_acc}: {e}")
        for gsm in gsm_rows:
            failed.append({"slide": gsm, "technology": technology, "error": str(e)})
        continue

    # Now convert each extracted sample.
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

# =============================================================================
# STEP 22 — Process TLS Zenodo (zip download + per-sample extraction)
# =============================================================================
# Download the single 2 GB ZIP once (with resume on failure), then extract
# and convert each of the 8 samples from it.

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
        skipped += len(sample_rows); continue

    try:
        # Download the ZIP (resumes automatically if partially downloaded).
        zip_path = _download_tls_zip(record_url, tls_cache)
    except Exception as e:
        print(f"  ERROR downloading TLS zip: {e}")
        for sn in sample_rows:
            failed.append({"slide": sn, "technology": technology, "error": str(e)})
        continue

    # Extract and convert each sample from the ZIP.
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

# =============================================================================
# STEP 23 — Final summary
# =============================================================================
print(f"\nFinished: {done} converted, {skipped} skipped (already done), {len(failed)} failed")
if failed:
    print("\nFailed slides:")
    for f in failed:
        print(f"  {f['slide']} ({f['technology']}): {f['error']}")
