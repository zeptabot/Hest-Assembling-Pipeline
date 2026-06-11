#!/usr/bin/env python3
import csv, os, re, sys

# ── Per-machine constant ───────────────────────────────────────────────────
HEST_META = "/Users/bradzap/Developer/GitHub/HEST/assets/HEST_v1_1_0.csv"

# ── Args ───────────────────────────────────────────────────────────────────
if len(sys.argv) != 2:
    print("Usage: Clean.py <path/to/meta_all_gene.csv>")
    sys.exit(1)

STIMAGE_META = sys.argv[1]

output_base = input("Where do you wish to store cleaned metadata? Please input directory: ").strip()
OUTPUT_DIR  = os.path.join(output_base, "meta")
print(f"Notice: metadata will be stored under newly created folder: {OUTPUT_DIR}")
os.makedirs(OUTPUT_DIR, exist_ok=True)

CLEANED_PATH   = os.path.join(OUTPUT_DIR, "cleaned_metadata.csv")
AMBIGUOUS_PATH = os.path.join(OUTPUT_DIR, "ambiguous_metadata.csv")

# ── If outputs already exist, skip cleaning and go straight to merge ───────
if os.path.exists(CLEANED_PATH) and os.path.exists(AMBIGUOUS_PATH):
    print(f"\nExisting output detected — skipping cleaning step.")
else:
    print("Cleaning…")

    with open(HEST_META, encoding="utf-8-sig") as f:
        hest = list(csv.DictReader(f))
    with open(STIMAGE_META) as f:
        stimage = list(csv.DictReader(f))

    hest_title_to_subseries: dict[str, list[str]] = {}
    for row in hest:
        title = row["dataset_title"].strip().lower()
        sub   = row["subseries"].strip()
        hest_title_to_subseries.setdefault(title, [])
        if sub:
            hest_title_to_subseries[title].append(sub)

    def extract_titles(raw: str) -> list[str]:
        titles = re.findall(r"Title \d+: (.+?)(?= Title \d+:|$)", raw)
        return [t.strip() for t in titles] if titles else [raw.strip()]

    cleaned, ambiguous, duplicates = [], [], []

    for row in stimage:
        slide  = row["slide"]
        titles = extract_titles(row["title"])

        matched_title = None
        for t in titles:
            if t.lower() in hest_title_to_subseries:
                matched_title = t.lower()
                break

        if matched_title is None:
            cleaned.append(row)
            continue

        subseries_list = hest_title_to_subseries[matched_title]
        confirmed = any(sub and sub in slide for sub in subseries_list)

        if confirmed:
            duplicates.append(row)
        else:
            row["_matched_hest_title"] = matched_title
            row["_hest_subseries"]     = " | ".join(subseries_list[:10])
            ambiguous.append(row)

    fieldnames = list(stimage[0].keys())

    with open(CLEANED_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(cleaned)

    amb_fields = fieldnames + ["_matched_hest_title", "_hest_subseries"]
    with open(AMBIGUOUS_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=amb_fields, extrasaction="ignore")
        w.writeheader(); w.writerows(ambiguous)

    print(f"Total STimage slides : {len(stimage)}")
    print(f"Confirmed duplicates : {len(duplicates)}  (skipped)")
    print(f"Ambiguous            : {len(ambiguous)}   → ambiguous_metadata.csv")
    print(f"Safe to convert      : {len(cleaned)}  → cleaned_metadata.csv")

# ── Manual review gate ─────────────────────────────────────────────────────
print(f"\nAmbiguous slides are in: {AMBIGUOUS_PATH}")
print("Please open it, review each slide, and delete any rows that are")
print("confirmed duplicates of HEST samples. Save the file when done.")
print()

while True:
    answer = input("Merge remaining ambiguous rows into cleaned_metadata.csv? [Y/N]: ").strip().upper()
    if answer in ("Y", "N"):
        break
    print("Please enter Y or N.")

if answer == "Y":
    with open(AMBIGUOUS_PATH, newline="") as f:
        amb_rows = list(csv.DictReader(f))
    with open(CLEANED_PATH, newline="") as f:
        cleaned_fields = csv.DictReader(f).fieldnames
    with open(CLEANED_PATH, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cleaned_fields, extrasaction="ignore")
        w.writerows(amb_rows)
    total = sum(1 for _ in open(CLEANED_PATH)) - 1
    print(f"Appended {len(amb_rows)} rows. Total safe-to-convert slides: {total}")
else:
    print("Ambiguous rows not added. cleaned_metadata.csv unchanged.")
