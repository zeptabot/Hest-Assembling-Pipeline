#!/usr/bin/env python3
import csv, os, re, shutil, sys, tempfile

import numpy as np
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

# ── Metadata helpers ────────────────────────────────────────────────────────
def _first_title(raw):
    m = re.findall(r"Title \d+: (.+?)(?= Title \d+:|$)", raw)
    return m[0].strip() if m else raw.strip()

def _norm_species(s):
    return {"human": "Homo sapiens", "mouse": "Mus musculus"}.get(s.strip().lower(), s.strip().capitalize())

def _first_study_link(raw):
    pmids = [p.strip() for p in re.split(r"[\n,;]+", raw) if p.strip().isdigit()]
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmids[0]}/" if pmids else None

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
            with tempfile.TemporaryDirectory() as tmp:
                src = hf_hub_download(
                    repo_id="jiawennnn/STimage-1K4M",
                    filename=f"{technology}/{file_type}/{filename}",
                    repo_type="dataset",
                    local_dir=tmp,
                )
                shutil.copy(src, dest_path)
        except Exception as e:
            print(f"    x {filename}: {e}")
            return None

    return {"stimage_dir": stimage_dir, "hest_dir": hest_dir, "technology": technology}

# ── Convert ────────────────────────────────────────────────────────────────
def convert_slide(slide_name, stimage_dir, hest_dir, technology, row):
    img_array = np.array(Image.open(os.path.join(stimage_dir, f"{slide_name}.png")).convert("RGB"))
    counts    = pd.read_csv(os.path.join(stimage_dir, f"{slide_name}_count.csv"), index_col=0)
    coord     = pd.read_csv(os.path.join(stimage_dir, f"{slide_name}_coord.csv"), index_col=0)

    common       = counts.index.intersection(coord.index)
    counts, coord = counts.loc[common], coord.loc[common]

    short_ids    = [re.search(r"(\d+x\d+)$", idx).group(1) for idx in counts.index]
    counts.index = coord.index = short_ids

    adata = sc.AnnData(counts)
    adata.obsm["spatial"] = coord[["xaxis", "yaxis"]].values

    my_df = pd.DataFrame(
        adata.obsm["spatial"], index=adata.obs_names,
        columns=["pxl_col_in_fullres", "pxl_row_in_fullres"]
    )
    my_df["array_row"] = [round(float(i.split("x")[0])) for i in my_df.index]
    my_df["array_col"] = [round(float(i.split("x")[1])) for i in my_df.index]

    adata.obs["array_row"]          = my_df["array_row"].values
    adata.obs["array_col"]          = my_df["array_col"].values
    adata.obs["pxl_col_in_fullres"] = my_df["pxl_col_in_fullres"].values
    adata.obs["pxl_row_in_fullres"] = my_df["pxl_row_in_fullres"].values
    adata.obs["in_tissue"]          = True

    pixel_size, spot_estimate_dist = find_pixel_size_from_spot_coords(
        my_df, inter_spot_dist=SPOT_PARAMS[technology]["inter_spot_dist"],
        packing=SpotPacking.GRID_PACKING
    )

    wsi = wsi_factory(img_array, mpp=pixel_size)
    register_downscale_img(adata, wsi, pixel_size,
                           spot_size=SPOT_PARAMS[technology]["spot_diameter"])

    adata.obs.index = [
        i.split("x")[0].zfill(3) + "x" + i.split("x")[1].zfill(3)
        for i in adata.obs.index
    ]

    adata.var["mito"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)

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
    st.save_spatial_plot(save_path=hest_dir)
    st.segment_tissue(method="otsu")
    st.save_tissue_contours(hest_dir, "tissue_seg")
    st.save_tissue_vis(hest_dir, "tissue_seg")
    st.dump_patches(patch_save_dir=hest_dir, name="patches",
                    target_patch_size=256, target_pixel_size=0.5,
                    dump_visualization=True)

# ── Main loop ──────────────────────────────────────────────────────────────
with open(CLEANED_META, newline="") as f:
    slides = list(csv.DictReader(f))

failed  = []
skipped = 0
done    = 0

for i, row in enumerate(slides):
    slide_name = row["slide"]
    technology = row["tech"]
    suffix     = TECH_SUFFIX.get(technology, technology.lower())
    hest_dir   = os.path.join(DATA_DIR, f"hest_{suffix}", slide_name)

    if os.path.exists(os.path.join(hest_dir, "aligned_adata.h5ad")):
        skipped += 1
        continue

    print(f"[{i+1}/{len(slides)}] {slide_name} ({technology})")
    try:
        result = download_slide(slide_name, technology)
        if result is None:
            raise RuntimeError("download failed")
        convert_slide(slide_name, result["stimage_dir"], result["hest_dir"], technology, row)
        done += 1
        print(f"  done")
    except Exception as e:
        print(f"  ERROR: {e}")
        failed.append({"slide": slide_name, "technology": technology, "error": str(e)})

print(f"\nFinished: {done} converted, {skipped} skipped (already done), {len(failed)} failed")
if failed:
    print("\nFailed slides:")
    for f in failed:
        print(f"  {f['slide']} ({f['technology']}): {f['error']}")
