# HEST Assembling Pipeline

Converts spatial transcriptomics (ST) datasets from public repositories into [HEST](https://github.com/mahmoodlab/HEST) format for model training. Two scripts handle all datasets: `Clean.py` produces a metadata CSV, `Assemble.py` downloads raw data and writes HEST-format output.

![VisualONE](https://raw.githubusercontent.com/zeptabot/Hest-Assembling-Pipeline/main/VisualONE.png)

---

## Setup

Requires the `hest` conda environment:

```bash
conda activate hest
```

Both scripts prompt for an output directory at runtime. All output lands under `<output_dir>/data/` and metadata under `<output_dir>/meta/`.

---

## Usage

Every dataset follows the same two-step pattern:

```
python Clean.py <source>
# prompt: output directory
python Assemble.py <path/to/meta/cleaned_metadata.csv>
# prompt: output directory
```

`Clean.py` accepts six source types:

| Source type | Argument |
|---|---|
| STimage-1K4M CSV (local or URL) | path or `https://` URL ending in `.csv` |
| TNBC local directory | path to directory with `byArray/` |
| TNBC Zenodo | `https://zenodo.org/records/14204217` |
| TLS Visium Zenodo | `https://zenodo.org/records/14620362` |
| GEO accession | e.g. `GSE278936` |
| EBI accession | e.g. `E-MTAB-13530` |

Run each Clean + Assemble pair to completion before starting the next — the `cleaned_metadata.csv` is overwritten by each `Clean.py` run.

---

## Supported Datasets

### 1. STimage-1K4M

~934 Visium/ST slides not already in HEST. `Clean.py` filters out HEST duplicates and flags ambiguous overlaps for manual review.

```bash
python Clean.py https://raw.githubusercontent.com/JiawenChenn/STimage-1K4M/main/meta/meta_all_gene.csv
# prompt → /path/to/output
python Assemble.py /path/to/output/meta/cleaned_metadata.csv
# prompt → /path/to/output
```

`Clean.py` classifies each of the 1,149 STimage slides:
- **No HEST title match** → safe to convert (946 slides)
- **Title matches + subseries in slide name** → confirmed duplicate, silently dropped (123)
- **Title matches but subseries not found** → ambiguous, written to `ambiguous_metadata.csv` for manual review (80)

After cleaning, you are prompted:
```
Merge remaining ambiguous rows into cleaned_metadata.csv? [Y/N]:
```
- `Y` — appends surviving ambiguous rows
- `N` — leaves `cleaned_metadata.csv` unchanged

Output folder: `hest_{suffix}` (e.g. `hest_visium`, `hest_st`)

---

### 2. TNBC (Triple-Negative Breast Cancer)

ST-technology slides from [Zenodo record 14204217](https://zenodo.org/records/14204217).

```bash
python Clean.py https://zenodo.org/records/14204217
# prompt → /path/to/output
python Assemble.py /path/to/output/meta/cleaned_metadata.csv
# prompt → /path/to/output
```

Or, if the Zenodo data is already downloaded locally:

```bash
python Clean.py /path/to/local/tnbc_dir
# prompt → /path/to/output
python Assemble.py /path/to/output/meta/cleaned_metadata.csv
# prompt → /path/to/output
```

Output folder: `hest_tnbc`

---

### 3. GSE278936 — Prostate Cancer (GEO Visium)

10x Visium from [GSE278936](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278936) (52 samples). `Clean.py` fetches the sample list from NCBI GEO; `Assemble.py` downloads the `_RAW.tar` and extracts each GSM.

```bash
python Clean.py GSE278936
# prompt → /path/to/output
python Assemble.py /path/to/output/meta/cleaned_metadata.csv
# prompt → /path/to/output
```

Output folder: `hest_gse278936`

---

### 4. E-MTAB-13530 — Lung NSCLC (EBI Visium)

10x Visium from [E-MTAB-13530](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530) (40 samples). `Clean.py` enumerates samples from the BioStudies API; `Assemble.py` downloads each sample's `.h5` and `spatial.tar` individually.

```bash
python Clean.py E-MTAB-13530
# prompt → /path/to/output
python Assemble.py /path/to/output/meta/cleaned_metadata.csv
# prompt → /path/to/output
```

Output folder: `hest_e_mtab_13530`

---

### 5. TLS Visium — Kidney & Lung (Zenodo)

10x Visium from [Zenodo record 14620362](https://zenodo.org/records/14620362) (8 samples: KC1–KC3 kidney, LC1–LC5 lung). A single 2 GB ZIP is downloaded once and cached; subsequent runs use the cache. The download is resume-capable and retries automatically on connection drops (up to 10 attempts).

```bash
python Clean.py https://zenodo.org/records/14620362
# prompt → /path/to/output
python Assemble.py /path/to/output/meta/cleaned_metadata.csv
# prompt → /path/to/output
```

Output folder: `hest_visium_tls`

---

## Output Format (per sample)

Each converted sample produces a flat directory:

```
{sample_id}/
├── aligned_adata.h5ad          # ST expression + spatial coords
├── aligned_fullres_HE.tif      # Pyramidal TIFF whole slide image
├── downscaled_fullres.jpeg     # Thumbnail
├── metrics.json                # Pixel size and spot metrics
├── patches_20x.h5              # 256×256 patches at 0.5 µm/px
├── patches_10x.h5              # 256×256 patches at 1.0 µm/px
├── patches_5x.h5               # 256×256 patches at 2.0 µm/px
├── patches_20x_patch_vis.png   # Patch extraction visualization
├── patches_10x_patch_vis.png
├── patches_5x_patch_vis.png
├── pixel_size_vis.png          # Pixel size calibration visualization
├── spatial_plots.png           # ST spots overlaid on WSI
├── tissue_seg_contours.geojson # Tissue segmentation mask
└── tissue_seg_vis.jpg          # Tissue mask visualization
```

### `aligned_adata.h5ad` fields

**`adata.obs`** (per spot):

| Field | Description |
|---|---|
| `in_tissue` | Whether spot is within tissue |
| `pxl_col_in_fullres` | Spot x-coordinate in full-res image |
| `pxl_row_in_fullres` | Spot y-coordinate in full-res image |
| `array_col` | Spot column in array grid |
| `array_row` | Spot row in array grid |
| `n_counts` | Total UMI counts |
| `n_genes_by_counts` | Number of detected genes |
| `log1p_n_genes_by_counts` | log1p of above |
| `total_counts` | Total counts |
| `log1p_total_counts` | log1p of above |
| `pct_counts_in_top_50_genes` | % counts in top 50 genes |
| `pct_counts_in_top_100_genes` | % counts in top 100 genes |
| `pct_counts_in_top_200_genes` | % counts in top 200 genes |
| `pct_counts_in_top_500_genes` | % counts in top 500 genes |
| `total_counts_mito` | Mitochondrial counts |
| `log1p_total_counts_mito` | log1p of above |
| `pct_counts_mito` | % mitochondrial counts |

**`adata.var`** (per gene):

| Field | Description |
|---|---|
| `n_cells_by_counts` | Number of spots with non-zero counts |
| `mean_counts` | Mean counts across spots |
| `log1p_mean_counts` | log1p of above |
| `pct_dropout_by_counts` | % spots with zero counts |
| `total_counts` | Total counts across all spots |
| `log1p_total_counts` | log1p of above |
| `mito` | Whether gene is mitochondrial |

**`adata.uns['spatial']['ST']['images']['downscaled_fullres']`** — downscaled WSI array

**`adata.obsm['spatial']`** — `(N, 2)` array of pixel coordinates `[x, y]`

---

## STimage-1K4M Input Format

For reference, the STimage-1K4M dataset uses this layout (handled internally by `Assemble.py`):

```
image/
  {id}.png                  # Full-resolution H&E slide
coord/
  {id}_coord.csv            # xaxis, yaxis, r (spot radius in pixels)
gene_exp/
  {id}_count.csv            # Spots × genes count matrix
```

---

## Notes

- **Resumable downloads**: large files (GEO RAW tar, TLS zip) are cached under `<output_dir>/data/raw_cache/`. If a run is interrupted, restarting picks up from where it left off.
- **Skip completed**: `Assemble.py` skips any sample whose output directory already exists, so interrupted runs are safe to restart.
- **HPC use**: no manual pre-downloading needed. All downloads including the 2 GB TLS zip are handled automatically with retry logic.
