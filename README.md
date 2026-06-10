# HEST Sample Assmeble Pipeline.

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
![VisualTWO](https://raw.githubusercontent.com/zeptabot/Hest-Assembling-Pipeline/main/VisualTWO.png)

```
╔══════════════════════════════════════════════════════════════════════════════════╗
║                         STimage → HEST TRANSFORMATION                          ║
╚══════════════════════════════════════════════════════════════════════════════════╝

┌─────────────────────────────────────────────────────────────────────────────────┐
│ INPUTS (STimage)                                                                │
├──────────────────┬──────────────────────────┬──────────────────────────────────┤
│ {id}.png         │ {id}_coord.csv           │ {id}_count.csv                   │
│                  │                          │                                  │
│ H&E image        │ spot_id | xaxis | yaxis  │ spot_id | GENE1 | GENE2 | ...    │
│                  │ 10x13   | 1923  | 2622   │ 10x13   |   0   |   1   | ...    │
│                  │ 10x14   | 1926  | 2844   │ 10x14   |   0   |   0   | ...    │
└──────────────────┴──────────────────────────┴──────────────────────────────────┘
         │                      │                           │
         │           coord.rename(xaxis→X, yaxis→Y)        │
         │                      │                           │
         └──────────────────────┴───────────────────────────┘
                                │
                                ▼
                     STReader().read(img_path,
                       raw_counts_path, spot_coord_path)
                                │
                                ▼
          ┌─────────────────────────────────────────────────┐
          │  HESTData object (st)                           │
          │                                                 │
          │  adata.X  ←  count.csv entire matrix           │
          │                                                 │
          │  adata.obs  (one row per spot):                 │
          │   ├─ pxl_col_in_fullres  ←  coord xaxis        │
          │   ├─ pxl_row_in_fullres  ←  coord yaxis        │
          │   ├─ array_row           ←  "10" from 10x13    │
          │   ├─ array_col           ←  "13" from 10x13    │
          │   ├─ in_tissue           ←  hardcode True      │
          │   │                                             │
          │   │  [sc.pp.calculate_qc_metrics() — AUTO]     │
          │   ├─ n_counts            ←  row.sum()          │
          │   ├─ total_counts        ←  row.sum()          │
          │   ├─ log1p_total_counts  ←  log(total+1)       │
          │   ├─ n_genes_by_counts   ←  (row>0).sum()      │
          │   ├─ log1p_n_genes_by_counts ← log(n_genes+1)  │
          │   ├─ pct_counts_in_top_50_genes                 │
          │   ├─ pct_counts_in_top_100_genes  ← sort row   │
          │   ├─ pct_counts_in_top_200_genes    desc, sum  │
          │   ├─ pct_counts_in_top_500_genes    top N /    │
          │   ├─ total_counts_mito   ← row[MT-*].sum()     │
          │   ├─ log1p_total_counts_mito                    │
          │   └─ pct_counts_mito     ← mito/total×100      │
          │                                                 │
          │  adata.var  (one row per gene):                 │
          │   ├─ mito               ← name.startswith(MT-) │
          │   │                                             │
          │   │  [sc.pp.calculate_qc_metrics() — AUTO]     │
          │   ├─ total_counts       ←  col.sum()           │
          │   ├─ log1p_total_counts ←  log(total+1)        │
          │   ├─ mean_counts        ←  col.mean()          │
          │   ├─ log1p_mean_counts  ←  log(mean+1)         │
          │   ├─ n_cells_by_counts  ←  (col>0).sum()       │
          │   └─ pct_dropout_by_counts ← (col==0)/n×100    │
          │                                                 │
          │  adata.uns['spatial']['ST']['images']           │
          │   └─ downscaled_fullres ← png downscaled       │
          │      [via register_downscale_img() — AUTO]      │
          │                                                 │
          │  adata.obsm['spatial']                          │
          │   └─ [[xaxis, yaxis], ...]  ←  coord CSV       │
          └──────────────────────┬──────────────────────────┘
                                 │
       ┌─────────────────────────┼──────────────────────────────────────┐
       │                         │                                      │
       ▼                         ▼                                      ▼
st.save(                 st.save_spatial_plot()           st.segment_tissue()
  save_img=True,                 │                        st.save_tissue_seg_pkl()
  plot_pxl_size=True)            │                                      │
       │                         │                                      │
       ▼                         ▼                                      ▼
┌─────────────────┐   ┌──────────────────┐               ┌─────────────────────┐
│ wsis/           │   │ spatial_plots/   │               │ tissue_seg/         │
│  .tif           │   │  spots overlaid  │               │  {id}.geojson       │
│                 │   │  on downscaled   │               │  {id}_vis.jpg       │
│ st/             │   │  WSI             │               └─────────────────────┘
│  .h5ad          │   └──────────────────┘
│                 │
│ metadata/       │              ▼
│  metrics.json   │   st.dump_patches(
│                 │     target_patch_size=256,
│ thumbnails/     │     target_pixel_size=0.5)
│  .jpeg          │              │
│                 │     ┌────────┴────────┐
│ pixel_size_vis/ │     ▼                ▼
│  .png           │  ┌──────────┐  ┌─────────────┐
└─────────────────┘  │ patches/ │  │ patches_vis/│
                     │  .h5     │  │  _patch_    │
                     │          │  │  vis.png    │
                     └──────────┘  └─────────────┘

┌──────────────────────────────────────────────┐
│ N/A for STimage (Xenium only — skip)         │
│  transcripts/   cellvit_seg/   xenium_seg/   │
└──────────────────────────────────────────────┘
```
