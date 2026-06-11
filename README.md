# HEST Sample Assmeble Pipeline.

``` bash
(base) bradzap@Brads-MacBook-Air ~ % conda activate hest
(hest) bradzap@Brads-MacBook-Air ~ % "/Users/bradzap/Developer/GitHub/Hest Assembling Pipeline/Clean.py" /Users/bradzap/Developer/GitHub/STimage-1K4M/meta/meta_all_gene.csv
Where do you wish to store cleaned metadata? Please input directory: /Users/bradzap/Developer
Notice: metadata will be stored under newly created folder: /Users/bradzap/Developer/meta
Cleaning…
Total STimage slides : 1149
Confirmed duplicates : 123  (skipped)
Ambiguous            : 80   → ambiguous_metadata.csv
Safe to convert      : 946  → cleaned_metadata.csv

Ambiguous slides are in: /Users/bradzap/Developer/meta/ambiguous_metadata.csv
Please open it, review each slide, and delete any rows that are
confirmed duplicates of HEST samples. Save the file when done.

Merge remaining ambiguous rows into cleaned_metadata.csv? [Y/N]: N
Ambiguous rows not added. cleaned_metadata.csv unchanged.
(hest) bradzap@Brads-MacBook-Air ~ % python3 "/Users/bradzap/Developer/GitHub/Hest Assembling Pipeline/Assemble.py" /Users/bradzap/Developer/meta/cleaned_metadata.csv
Where do you wish to store the data? Please input directory: /Users/bradzap/Developer
Notice: data will be stored under newly created folder: /Users/bradzap/Developer/data
This might take a long time. Converting…
```

## Target Format
- **wsis/**: H&E-stained whole slide images in pyramidal Generic TIFF (or pyramidal Generic BigTIFF if >4.1GB)
- **st/**: Spatial transcriptomics expressions in a scanpy `.h5ad` object
    - **Observations** (`st.adata.obs`)
        - `in_tissue`: Indicator if the observation is within the tissue (`in_tissue` comes from the initial Visium/Xenium run and might not be accurate, prefer the segmentation obtained by `st.segment_tissue()` instead)
        - `pxl_col_in_fullres`: Pixel column position of the patch/spot centroid in the full resolution image
        - `pxl_row_in_fullres`: Pixel row position of the patch/spot centroid in the full resolution image
        - `array_col`: Patch/spot column position in the array
        - `array_row`: Patch/spot row position in the array
        - `n_counts`: Number of counts for each observation
        - `n_genes_by_counts`: Number of genes detected by counts in each observation
        - `log1p_n_genes_by_counts`: Log-transformed number of genes detected by counts
        - `total_counts`: Total counts per observation
        - `log1p_total_counts`: Log-transformed total counts
        - `pct_counts_in_top_50_genes`: Percentage of counts in the top 50 genes
        - `pct_counts_in_top_100_genes`: Percentage of counts in the top 100 genes
        - `pct_counts_in_top_200_genes`: Percentage of counts in the top 200 genes
        - `pct_counts_in_top_500_genes`: Percentage of counts in the top 500 genes
        - `total_counts_mito`: Total mitochondrial counts per observation *(may not be accurate)*
        - `log1p_total_counts_mito`: Log-transformed total mitochondrial counts *(may not be accurate)*
        - `pct_counts_mito`: Percentage of counts that are mitochondrial *(may not be accurate)*
    - **Variables** (`st.adata.var`)
        - `n_cells_by_counts`: Number of cells detected by counts for each variable
        - `mean_counts`: Mean counts per variable
        - `log1p_mean_counts`: Log-transformed mean counts
        - `pct_dropout_by_counts`: Percentage of dropout events by counts
        - `total_counts`: Total counts per variable
        - `log1p_total_counts`: Log-transformed total counts
        - `mito`: Indicator if the gene is mitochondrial
    - **Unstructured** (`st.adata.uns`)
        - `spatial`: Contains a downscaled version of the full resolution image in `st.adata.uns['spatial']['ST']['images']['downscaled_fullres']`
    - **Observation-wise Multidimensional** (`st.adata.obsm`)
        - `spatial`: Pixel coordinates of spots/patches centroids on the full resolution image (first column = x axis, second column = y axis)
- **metadata/**: Metadata
- **spatial_plots/**: Overlay of the WSI with the ST spots
- **thumbnails/**: Downscaled version of the WSI
- **tissue_seg/**: Tissue segmentation masks
    - `{id}.geojson`: Tissue segmentation mask
    - `{id}_vis.jpg`: Visualization of the tissue mask on the downscaled WSI
- **pixel_size_vis/**: Visualization of the pixel size
- **patches/**: 256×256 H&E patches (0.5µm/px) extracted around ST spots in a `.h5` object optimized for deep learning; each patch is matched to the corresponding ST profile (see **st/**) with a barcode
- **patches_vis/**: Visualization of the mask and patches on a downscaled WSI
- **transcripts/**: Individual transcripts aligned to H&E for Xenium samples; read with `pandas.read_parquet`; aligned coordinates in pixels are in columns `['he_x', 'he_y']`
- **cellvit_seg/**: CellViT nuclei segmentation
- **xenium_seg/**: Xenium segmentation on DAPI aligned to H&E

## Input Format (Current Version)
- **image**: H&E-stained histopathology image
    - `{id}.png`: Full resolution tissue slide image
- **coord/**: Spot coordinates
    - `{id}_coord.csv`: Pixel positions of each capture spot on the full resolution image
        - `yaxis`: Pixel row position of the spot centroid (vertical)
        - `xaxis`: Pixel column position of the spot centroid (horizontal)
        - `r`: Radius of the spot in pixels
- **gene_exp/**: Gene expression count matrix
    - `{id}_count.csv`: RNA counts per spot per gene
        - Rows: spot IDs in format `{sample}_{row}x{col}` (e.g. `ST_A1_10x13`)
        - Columns: gene names (e.g. `SAMD11`, `NOC2L`, ...) — up to ~15,000 genes
        - Values: integer RNA read counts (sparse — mostly zeros)

## Visual Representations
![VisualONE](https://raw.githubusercontent.com/zeptabot/Hest-Assembling-Pipeline/main/VisualONE.png)

# Clean.py

- Set directory where you want to store the cleaned metadata
- **Inputs**
  - `STimage-1K4M/meta/meta_all_gene.csv` — 1,149 slides
  - `HEST/assets/HEST_v1_1_0.csv` — all HEST samples
- **Build HEST lookup**
  - Key: `dataset_title` (lowercased)
  - Value: list of `subseries` labels for that study
- **For each STimage slide → 3-way classification**
  - Extract study title(s) from the `title` field
    - Handles multi-paper cells: `"Title 1: Foo. Title 2: Bar."`
  - **No title match** → `cleaned` (946) — not in HEST, safe to convert
  - **Title matches + subseries is substring of slide name** → `duplicates` (123) — confirmed in HEST, silently dropped
  - **Title matches but subseries not found in slide name** → `ambiguous` (80) — same study, different naming convention, needs human review
- **Outputs**
  - `cleaned_metadata.csv` — 946 rows, same columns as input
  - `ambiguous_metadata.csv` — 80 rows + two extra debug columns:
    - `_matched_hest_title` — which HEST study matched
    - `_hest_subseries` — HEST's subseries labels for that study
- Tells user to open `ambiguous_metadata.csv` and delete confirmed duplicates
- Prompts `[Y/N]`
  - **Y** — appends surviving ambiguous rows into `cleaned_metadata.csv`
  - **N** — leaves `cleaned_metadata.csv` unchanged
- Final `cleaned_metadata.csv` = the definitive list of STimage slides to convert