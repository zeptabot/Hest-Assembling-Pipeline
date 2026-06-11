# Example run — STimage → HEST pipeline

A real end-to-end session: `Clean.py` deduplicates STimage against HEST, then `Assemble.py` converts the safe slides. (GitHub renders the colors below.)

```ansi
[1;32m(base) bradzap@Brads-MacBook-Air ~ % [0m[1;36mconda activate hest[0m
[1;32m(hest) bradzap@Brads-MacBook-Air ~ % [0m[1;36mpython3 "/Users/bradzap/Developer/GitHub/Hest Assembling Pipeline/Clean.py" /Users/bradzap/Developer/GitHub/STimage-1K4M/meta/meta_all_gene.csv[0m
[0;37mWhere do you wish to store cleaned metadata? Please input directory: [0m[1;33m/Users/bradzap/Developer[0m
[0;37mNotice: metadata will be stored under newly created folder: /Users/bradzap/Developer/meta[0m
[0;37mCleaning…[0m
[0;37mTotal STimage slides : 1149[0m
[0;37mConfirmed duplicates : 123  (skipped)[0m
[0;37mAmbiguous            : 80   → ambiguous_metadata.csv[0m
[0;37mSafe to convert      : 946  → cleaned_metadata.csv[0m
[0;37m[0m
[0;37mAmbiguous slides are in: /Users/bradzap/Developer/meta/ambiguous_metadata.csv[0m
[0;37mPlease open it, review each slide, and delete any rows that are[0m
[0;37mconfirmed duplicates of HEST samples. Save the file when done.[0m
[0;37m[0m
[0;37mMerge remaining ambiguous rows into cleaned_metadata.csv? [Y/N]: [0m[1;33mN[0m
[0;37mAmbiguous rows not added. cleaned_metadata.csv unchanged.[0m
[1;32m(hest) bradzap@Brads-MacBook-Air ~ % [0m[1;36mpython3 "/Users/bradzap/Developer/GitHub/Hest Assembling Pipeline/Assemble.py" /Users/bradzap/Developer/meta/cleaned_metadata.csv[0m
[0;37mWhere do you wish to store the data? Please input directory: [0m[1;33m/Users/bradzap/Developer[0m
[0;37mNotice: data will be stored under newly created folder: /Users/bradzap/Developer/data[0m
[0;37mThis might take a long time. Converting…[0m
[0;37m[1/946] GSE144239_GSM4284316 (ST)[0m
[2;37mGSE144239_GSM4284316_coord.csv: 39.4kB [00:00, 23.6MB/s][0m
[2;37mST/gene_exp/GSE144239_GSM4284316_count.c(…): 100%|█████████| 45.9M/45.9M [00:01<00:00, 45.5MB/s][0m
[2;37mST/image/GSE144239_GSM4284316.png: 100%|███████████████████| 11.8M/11.8M [00:00<00:00, 12.5MB/s][0m
[0;33mobjc[27700]: Class GNotificationCenterDelegate is implemented in both …libopenslide.1.dylib and …libgio-2.0.0.dylib. One of the duplicates must be removed or renamed.[0m
[0;37msaving to pyramidal tiff... can be slow[0m
[0;33m…/hest/HESTData.py:1657: FutureWarning: Use `squidpy.pl.spatial_scatter` instead.[0m
[0;33m  fig = sc.pl.spatial(adata, show=False, img_key="downscaled_fullres", color=[key], …)[0m
[0;33m…/trident/segmentation_models/model_zoo/otsu.py:49: FutureWarning: Parameter `min_size` is deprecated since 0.26.0; use `max_size` instead.[0m
[0;33m  otsu_masking = sk_morphology.remove_small_objects(otsu_masking, 60)[0m
[0;33m…/trident/segmentation_models/model_zoo/otsu.py:56: FutureWarning: Parameter `area_threshold` is deprecated since 0.26.0; use `max_size` instead.[0m
[0;33m  otsu_masking = sk_morphology.remove_small_holes(otsu_masking, 5000)[0m
[1;32m  ✓ done[0m
[0;37m[2/946] GSE144239_GSM4284317 (ST)[0m
[2;37mGSE144239_GSM4284317_coord.csv: 38.1kB [00:00, 89.0MB/s][0m
[2;37mST/gene_exp/GSE144239_GSM4284317_count.c(…): 100%|█████████| 45.1M/45.1M [00:00<00:00, 55.9MB/s][0m
[2;37mST/image/GSE144239_GSM4284317.png: 100%|███████████████████| 11.8M/11.8M [00:00<00:00, 14.1MB/s][0m
[1;31m^CTraceback (most recent call last):[0m
[1;31m  File "…/Hest Assembling Pipeline/Assemble.py", line 159, in <module>[0m
[1;31m    convert_slide(slide_name, result["stimage_dir"], result["hest_dir"], technology)[0m
[1;31m  File "…/Hest Assembling Pipeline/Assemble.py", line 72, in convert_slide[0m
[1;31m    img_array = np.array(Image.open(os.path.join(stimage_dir, f"{slide_name}.png")).convert("RGB"))[0m
[1;31m  File "…/PIL/Image.py", line 1070, in convert[0m
[1;31m    self.load()[0m
[1;31m  File "…/PIL/ImageFile.py", line 412, in load[0m
[1;31m    n, err_code = decoder.decode(b)[0m
[1;31mKeyboardInterrupt[0m
[1;32m(hest) bradzap@Brads-MacBook-Air ~ %[0m
```

**Legend** — 🟢 prompt · 🔵 typed command · 🟡 your input / warnings · ⚪ output · 🔴 traceback

> The `KeyboardInterrupt` at slide 2 is a manual **Ctrl-C**, not an error. Re-running `Assemble.py` resumes: slide 1 is skipped (already has `aligned_adata.h5ad`) and slide 2 restarts from its cached download.
