#!/usr/bin/env python3
import csv, os, re, sys, tarfile, tempfile

# ── Per-machine constant ───────────────────────────────────────────────────
HEST_META = "/Users/bradzap/Developer/GitHub/HEST/assets/HEST_v1_1_0.csv"

# ── Args ───────────────────────────────────────────────────────────────────
if len(sys.argv) != 2:
    print("Usage: Clean.py <path/to/meta_all_gene.csv>          (STimage-1K4M)")
    print("       Clean.py <path/to/local_zenodo_dir>           (TNBC local)")
    print("       Clean.py https://zenodo.org/records/14204217  (TNBC Zenodo)")
    print("       Clean.py https://zenodo.org/records/14620362  (TLS Visium Zenodo)")
    print("       Clean.py GSE278936                            (GEO Visium)")
    print("       Clean.py E-MTAB-13530                         (EBI Visium)")
    sys.exit(1)

ARG = sys.argv[1]

ZENODO_PREFIX = "https://zenodo.org/records/"
TNBC_RECORD   = "14204217"
TLS_RECORD    = "14620362"
GEO_RE        = re.compile(r'^GSE\d+$',       re.IGNORECASE)
EBI_RE        = re.compile(r'^E-[A-Z]+-\d+$', re.IGNORECASE)

is_zenodo = ARG.startswith(ZENODO_PREFIX)
is_geo    = bool(GEO_RE.match(ARG))
is_ebi    = bool(EBI_RE.match(ARG))
is_remote = is_zenodo or is_geo or is_ebi

if not is_remote and not os.path.exists(ARG):
    print(f"ERROR: path does not exist: {ARG}")
    sys.exit(1)

output_base = input("Where do you wish to store cleaned metadata? Please input directory: ").strip()
OUTPUT_DIR  = os.path.join(output_base, "meta")
print(f"Notice: metadata will be stored under newly created folder: {OUTPUT_DIR}")
os.makedirs(OUTPUT_DIR, exist_ok=True)

CLEANED_PATH   = os.path.join(OUTPUT_DIR, "cleaned_metadata.csv")
AMBIGUOUS_PATH = os.path.join(OUTPUT_DIR, "ambiguous_metadata.csv")

# Shared fieldnames for all non-STimage sources
FIELDS = ["slide", "tech", "title", "tissue", "species",
          "involve_cancer", "pmid", "gene_num", "source", "data_path"]

def _write_empty_ambiguous(fieldnames):
    with open(AMBIGUOUS_PATH, "w", newline="") as f:
        csv.DictWriter(f, fieldnames=fieldnames).writeheader()

# ── Known metadata per accession ──────────────────────────────────────────
_GEO_META = {
    "GSE278936": {
        "title":          ("Single cell and spatial transcriptomics highlight the interaction of "
                           "club-like cells with immunosuppressive myeloid cells in prostate cancer"),
        "tissue":         "Prostate",
        "species":        "human",
        "involve_cancer": "True",
        "pmid":           "39550375",
        "gene_num":       "",
        "tech":           "Visium",
    },
}

_EBI_META = {
    "E-MTAB-13530": {
        "title":          "Single cell and spatial transcriptomics analysis of non-small cell lung cancer",
        "tissue":         "Lung",
        "species":        "human",
        "involve_cancer": "True",
        "pmid":           "38782905",
        "gene_num":       "",
        "tech":           "Visium",
    },
}

_TNBC_META = {
    "title":          ("Spatial transcriptomics reveals substantial heterogeneity in "
                       "triple-negative breast cancer with potential clinical implications"),
    "tissue":         "Breast",
    "species":        "human",
    "involve_cancer": "True",
    "pmid":           "39562547",
    "gene_num":       "",
    "tech":           "ST",
}

_TLS_META_BASE = {
    "title":          ("DeepSpot: Leveraging Spatial Context for Enhanced Spatial "
                       "Transcriptomics Prediction from H&E Images"),
    "species":        "human",
    "involve_cancer": "True",
    "pmid":           "",
    "gene_num":       "",
    "tech":           "Visium",
}
_TLS_SAMPLES = ["KC1", "KC2", "KC3", "LC1", "LC2", "LC3", "LC4", "LC5"]

# ── Zenodo mode ────────────────────────────────────────────────────────────
if is_zenodo:
    try:
        import requests
    except ImportError:
        print("ERROR: pip install requests"); sys.exit(1)

    record_id = ARG.rstrip("/").split("/")[-1]

    # ── TNBC (record 14204217) ─────────────────────────────────────────────
    if record_id == TNBC_RECORD:
        try:
            import pyreadr
        except ImportError:
            print("ERROR: pip install pyreadr"); sys.exit(1)

        clinical_url = f"https://zenodo.org/records/{record_id}/files/Clinical.tar"
        print(f"Fetching array IDs from Clinical.tar (Zenodo record {record_id})…")

        ids_bytes = None
        with requests.get(clinical_url, stream=True, timeout=60) as resp:
            resp.raise_for_status()
            resp.raw.decode_content = True
            with tarfile.open(fileobj=resp.raw, mode="r|") as tar:
                for m in tar:
                    if m.name.endswith("ids.RDS") and m.isfile():
                        ids_bytes = tar.extractfile(m).read()
                        break

        if ids_bytes is None:
            print("ERROR: ids.RDS not found in Clinical.tar"); sys.exit(1)

        with tempfile.NamedTemporaryFile(suffix=".RDS", delete=False) as tmp:
            tmp.write(ids_bytes); tmp_path = tmp.name
        try:
            ids_df = list(pyreadr.read_r(tmp_path).values())[0]
        finally:
            os.unlink(tmp_path)

        array_ids = list(ids_df.index)
        with open(CLEANED_PATH, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            w.writeheader()
            for array_id in array_ids:
                w.writerow({**_TNBC_META, "slide": array_id,
                            "source": "tnbc_zenodo", "data_path": ARG})

        _write_empty_ambiguous(FIELDS)
        print(f"TNBC Zenodo: wrote {len(array_ids)} arrays → {CLEANED_PATH}")
        print(f"Ready: python Assemble.py {CLEANED_PATH}")
        sys.exit(0)

    # ── TLS Visium (record 14620362) ───────────────────────────────────────
    elif record_id == TLS_RECORD:
        with open(CLEANED_PATH, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            w.writeheader()
            for sn in _TLS_SAMPLES:
                tissue = "Kidney" if sn.startswith("KC") else "Lung"
                w.writerow({**_TLS_META_BASE, "slide": sn, "tissue": tissue,
                            "source": "visium_zenodo", "data_path": ARG})

        _write_empty_ambiguous(FIELDS)
        print(f"TLS Zenodo: wrote {len(_TLS_SAMPLES)} samples → {CLEANED_PATH}")
        print(f"Ready: python Assemble.py {CLEANED_PATH}")
        sys.exit(0)

    else:
        print(f"ERROR: Unsupported Zenodo record: {record_id}")
        print(f"  Supported: {TNBC_RECORD} (TNBC), {TLS_RECORD} (TLS Visium)")
        sys.exit(1)

# ── GEO Visium mode ────────────────────────────────────────────────────────
if is_geo:
    try:
        import requests
    except ImportError:
        print("ERROR: pip install requests"); sys.exit(1)

    geo_acc = ARG.upper()
    print(f"Fetching sample list from GEO {geo_acc}…")

    r = requests.get(
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
        params={"acc": geo_acc, "targ": "gsm", "view": "brief", "form": "text"},
        timeout=60)
    r.raise_for_status()
    gsm_ids = re.findall(r'^\^SAMPLE = (GSM\d+)', r.text, re.MULTILINE)

    if not gsm_ids:
        print(f"ERROR: No GSM samples found for {geo_acc}"); sys.exit(1)

    meta = _GEO_META.get(geo_acc, {
        "title": geo_acc, "tissue": "", "species": "human",
        "involve_cancer": "True", "pmid": "", "gene_num": "", "tech": "Visium",
    })

    with open(CLEANED_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for gsm_id in gsm_ids:
            w.writerow({**meta, "slide": gsm_id,
                        "source": "visium_geo", "data_path": geo_acc})

    _write_empty_ambiguous(FIELDS)
    print(f"GEO {geo_acc}: wrote {len(gsm_ids)} samples → {CLEANED_PATH}")
    print(f"Ready: python Assemble.py {CLEANED_PATH}")
    sys.exit(0)

# ── EBI Visium mode ────────────────────────────────────────────────────────
if is_ebi:
    try:
        import requests
    except ImportError:
        print("ERROR: pip install requests"); sys.exit(1)

    ebi_acc = ARG.upper()
    print(f"Fetching sample list from EBI {ebi_acc}…")

    r = requests.get(
        f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{ebi_acc}",
        timeout=60)
    r.raise_for_status()

    def _collect_paths(obj):
        paths = []
        if isinstance(obj, dict):
            if obj.get("type") == "file" and "path" in obj:
                paths.append(obj["path"])
            for v in obj.values():
                paths.extend(_collect_paths(v))
        elif isinstance(obj, list):
            for item in obj:
                paths.extend(_collect_paths(item))
        return paths

    all_paths = _collect_paths(r.json())
    suffix = "-filtered_feature_bc_matrix.h5"
    sample_names = list(dict.fromkeys(
        p[:-len(suffix)] for p in all_paths if p.endswith(suffix)
    ))

    if not sample_names:
        print(f"ERROR: Could not enumerate samples from {ebi_acc}")
        print(f"  Paths found: {all_paths[:5]}")
        sys.exit(1)

    meta = _EBI_META.get(ebi_acc, {
        "title": ebi_acc, "tissue": "", "species": "human",
        "involve_cancer": "True", "pmid": "", "gene_num": "", "tech": "Visium",
    })

    with open(CLEANED_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for sn in sample_names:
            w.writerow({**meta, "slide": sn,
                        "source": "visium_ebi", "data_path": ebi_acc})

    _write_empty_ambiguous(FIELDS)
    print(f"EBI {ebi_acc}: wrote {len(sample_names)} samples → {CLEANED_PATH}")
    print(f"Ready: python Assemble.py {CLEANED_PATH}")
    sys.exit(0)

# ── TNBC local mode: argument is a directory ───────────────────────────────
if os.path.isdir(ARG):
    by_array = os.path.join(ARG, "byArray")
    if not os.path.isdir(by_array):
        print(f"ERROR: {ARG} is a directory but has no byArray/ subdirectory."); sys.exit(1)

    if os.path.exists(CLEANED_PATH):
        print(f"\nExisting output detected — skipping TNBC metadata generation.")
    else:
        arrays = sorted(
            d for d in os.listdir(by_array)
            if os.path.isdir(os.path.join(by_array, d))
        )
        with open(CLEANED_PATH, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            w.writeheader()
            for array_id in arrays:
                w.writerow({**_TNBC_META, "slide": array_id,
                            "source": "tnbc", "data_path": ARG})
        print(f"TNBC: wrote {len(arrays)} arrays → {CLEANED_PATH}")

    _write_empty_ambiguous(FIELDS)
    print(f"Ready: python Assemble.py {CLEANED_PATH}")
    sys.exit(0)

# ── STimage mode (CSV file) ────────────────────────────────────────────────
STIMAGE_META = ARG

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

    cleaned, ambiguous, duplicates, non_cancer = [], [], [], []

    for row in stimage:
        if str(row.get("involve_cancer", "")).strip().lower() not in ("true", "1", "yes"):
            non_cancer.append(row)
            continue

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
    print(f"Non-cancer           : {len(non_cancer)}  (skipped)")
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
