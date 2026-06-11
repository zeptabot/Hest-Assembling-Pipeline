# Example run — STimage → HEST pipeline

A real end-to-end session: `Clean.py` deduplicates STimage against HEST,
then `Assemble.py` converts the safe slides.

```shell
(base) bradzap@Brads-MacBook-Air ~ % conda activate hest
(hest) bradzap@Brads-MacBook-Air ~ % python3 "/Users/bradzap/Developer/GitHub/Hest Assembling Pipeline/Clean.py" \
    /Users/bradzap/Developer/GitHub/STimage-1K4M/meta/meta_all_gene.csv
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

(hest) bradzap@Brads-MacBook-Air ~ % python3 "/Users/bradzap/Developer/GitHub/Hest Assembling Pipeline/Assemble.py" \
    /Users/bradzap/Developer/meta/cleaned_metadata.csv
Where do you wish to store the data? Please input directory: /Users/bradzap/Developer
Notice: data will be stored under newly created folder: /Users/bradzap/Developer/data
This might take a long time. Converting…

[1/946] GSE144239_GSM4284316 (ST)
GSE144239_GSM4284316_coord.csv: 39.4kB [00:00, 23.6MB/s]
ST/gene_exp/GSE144239_GSM4284316_count.c(…): 100%|█████████| 45.9M/45.9M [00:01<00:00, 45.5MB/s]
ST/image/GSE144239_GSM4284316.png: 100%|███████████████████| 11.8M/11.8M [00:00<00:00, 12.5MB/s]
objc[27700]: Class GNotificationCenterDelegate is implemented in both
  …libopenslide.1.dylib and …libgio-2.0.0.dylib. One of the duplicates must be removed or renamed.
saving to pyramidal tiff... can be slow
…/hest/HESTData.py:1657: FutureWarning: Use `squidpy.pl.spatial_scatter` instead.
  fig = sc.pl.spatial(adata, show=False, img_key="downscaled_fullres", color=[key], …)
…/trident/segmentation_models/model_zoo/otsu.py:49: FutureWarning: Parameter `min_size` is
  deprecated since 0.26.0; use `max_size` instead.
  otsu_masking = sk_morphology.remove_small_objects(otsu_masking, 60)
…/trident/segmentation_models/model_zoo/otsu.py:56: FutureWarning: Parameter `area_threshold`
  is deprecated since 0.26.0; use `max_size` instead.
  otsu_masking = sk_morphology.remove_small_holes(otsu_masking, 5000)
  ✓ done

[2/946] GSE144239_GSM4284317 (ST)
GSE144239_GSM4284317_coord.csv: 38.1kB [00:00, 89.0MB/s]
ST/gene_exp/GSE144239_GSM4284317_count.c(…): 100%|█████████| 45.1M/45.1M [00:00<00:00, 55.9MB/s]
ST/image/GSE144239_GSM4284317.png: 100%|███████████████████| 11.8M/11.8M [00:00<00:00, 14.1MB/s]
^CTraceback (most recent call last):
  File "…/Hest Assembling Pipeline/Assemble.py", line 159, in <module>
    convert_slide(slide_name, result["stimage_dir"], result["hest_dir"], technology)
  File "…/Hest Assembling Pipeline/Assemble.py", line 72, in convert_slide
    img_array = np.array(Image.open(os.path.join(stimage_dir, f"{slide_name}.png")).convert("RGB"))
  File "…/PIL/Image.py", line 1070, in convert
    self.load()
  File "…/PIL/ImageFile.py", line 412, in load
    n, err_code = decoder.decode(b)
KeyboardInterrupt
(hest) bradzap@Brads-MacBook-Air ~ %
```

> **Note:** The `KeyboardInterrupt` at slide 2 is a manual **Ctrl-C**, not an error.
> Re-running `Assemble.py` resumes automatically: slide 1 is skipped (already has
> `aligned_adata.h5ad`) and slide 2 restarts from its cached download.
