#!/usr/bin/env python3
# =============================================================================
# Clean.py
# =============================================================================
# PURPOSE: Produce a clean list of slides to convert, before any data is
# downloaded.  Every supported dataset goes through this script first.
#
# WHAT IT DOES IN PLAIN ENGLISH:
#   1. You tell it what dataset you want (a URL, an accession code, or a file).
#   2. It figures out exactly which individual tissue slides are in that dataset
#      by either reading a local file, or calling a public web API.
#   3. It writes those slide names + study metadata (tissue type, species,
#      cancer yes/no, PubMed ID, technology) into a single CSV called
#      cleaned_metadata.csv.
#   4. For the STimage dataset only, it also cross-checks against HEST to
#      remove slides that already exist in HEST (so we don't create duplicates).
#   5. Assemble.py then reads that CSV and does the actual downloading and
#      conversion.
#
# USAGE:
#   python Clean.py <source>
#
# WHERE <source> IS ONE OF:
#   A URL ending in .csv          → STimage-1K4M slide list from GitHub
#   A local path to a .csv file   → same, from disk
#   A local path to a directory   → TNBC dataset already downloaded locally
#   https://zenodo.org/records/14204217  → TNBC dataset on Zenodo
#   https://zenodo.org/records/14620362  → TLS kidney/lung dataset on Zenodo
#   GSE278936                     → prostate cancer dataset on NCBI GEO
#   E-MTAB-13530                  → lung cancer dataset on EBI BioStudies
# =============================================================================

# ── Standard library imports ──────────────────────────────────────────────────
# csv      : read and write comma-separated value files
# os       : interact with the filesystem (check paths, create folders, etc.)
# re       : regular expressions — pattern matching on text strings
# sys      : access command-line arguments and exit the program
# tarfile  : open .tar archive files (used for TNBC Clinical.tar on Zenodo)
# tempfile : create temporary files that are automatically deleted when done
import csv, os, re, sys, tarfile, tempfile

# ── Path to the HEST master catalogue ────────────────────────────────────────
# This CSV lists every slide already published in the HEST dataset.
# We use it in STimage mode to avoid creating duplicates.
# This path must point to a locally cloned copy of the HEST repository.
HEST_META = "/Users/bradzap/Developer/GitHub/HEST/assets/HEST_v1_1_0.csv"

# =============================================================================
# STEP 1 — Read the command-line argument
# =============================================================================
# The user must pass exactly one argument: the data source.
# If they pass nothing, or too many arguments, print instructions and quit.
if len(sys.argv) != 2:
    print("Usage: Clean.py <path/to/meta_all_gene.csv>          (STimage-1K4M)")
    print("       Clean.py <path/to/local_zenodo_dir>           (TNBC local)")
    print("       Clean.py https://zenodo.org/records/14204217  (TNBC Zenodo)")
    print("       Clean.py https://zenodo.org/records/14620362  (TLS Visium Zenodo)")
    print("       Clean.py GSE278936                            (GEO Visium)")
    print("       Clean.py E-MTAB-13530                         (EBI Visium)")
    sys.exit(1)

# Store the argument in a variable called ARG.
ARG = sys.argv[1]

# =============================================================================
# STEP 2 — Detect which type of source was provided
# =============================================================================

# The Zenodo URL prefix all our Zenodo records start with.
ZENODO_PREFIX = "https://zenodo.org/records/"

# The Zenodo record numbers for the two datasets we support via Zenodo.
TNBC_RECORD = "14204217"   # Triple-negative breast cancer (ST technology)
TLS_RECORD  = "14620362"   # Tumour-lymphocyte structures, kidney + lung (Visium)

# Regular expressions that match accession-code formats:
#   GEO accessions look like:  GSE278936
#   EBI accessions look like:  E-MTAB-13530
GEO_RE = re.compile(r'^GSE\d+$',       re.IGNORECASE)
EBI_RE = re.compile(r'^E-[A-Z]+-\d+$', re.IGNORECASE)

# Detect which category the argument falls into.
is_zenodo  = ARG.startswith(ZENODO_PREFIX)
is_geo     = bool(GEO_RE.match(ARG))
is_ebi     = bool(EBI_RE.match(ARG))
# A remote CSV is any http/https URL that ends in ".csv" (e.g. the raw GitHub link
# for STimage's meta_all_gene.csv).
is_csv_url = (ARG.startswith("http://") or ARG.startswith("https://")) and ARG.endswith(".csv")
# "is_remote" is True if we need to download anything at all (vs. a local path).
is_remote  = is_zenodo or is_geo or is_ebi or is_csv_url

# If the argument is not a recognised remote source AND the path doesn't exist
# on disk, we cannot proceed — tell the user and exit.
if not is_remote and not os.path.exists(ARG):
    print(f"ERROR: path does not exist: {ARG}")
    sys.exit(1)

# =============================================================================
# STEP 3 — If the source is a remote CSV URL, download it to a temporary file
# =============================================================================
# This handles the case where the user passes the raw GitHub URL for
# meta_all_gene.csv instead of a local copy.  We download it to a temp file
# so the rest of the STimage code can read it exactly like a local file.
if is_csv_url:
    try:
        import requests as _req
    except ImportError:
        print("ERROR: pip install requests"); sys.exit(1)
    print(f"Downloading {ARG} …")
    _r = _req.get(ARG, timeout=120)
    _r.raise_for_status()   # crash loudly if the server returns an error
    import tempfile as _tmp
    # NamedTemporaryFile creates a file that persists on disk until we delete it
    # (delete=False means it won't vanish when we close it).
    _tf = _tmp.NamedTemporaryFile(mode="wb", suffix=".csv", delete=False)
    _tf.write(_r.content); _tf.close()
    # Replace ARG with the path to the temp file so the rest of the code is
    # unaware anything was downloaded.
    ARG = _tf.name

# =============================================================================
# STEP 4 — Ask the user where to save the output
# =============================================================================
# We don't hardcode output paths so the same script works on any machine or HPC.
output_base = input("Where do you wish to store cleaned metadata? Please input directory: ").strip()
# We always write into a "meta/" subfolder to keep things tidy.
OUTPUT_DIR  = os.path.join(output_base, "meta")
print(f"Notice: metadata will be stored under newly created folder: {OUTPUT_DIR}")
# exist_ok=True means this doesn't fail if the folder already exists.
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Full file paths for the two output CSVs.
CLEANED_PATH   = os.path.join(OUTPUT_DIR, "cleaned_metadata.csv")   # safe to convert
AMBIGUOUS_PATH = os.path.join(OUTPUT_DIR, "ambiguous_metadata.csv") # needs human review

# =============================================================================
# STEP 5 — Define the column names for the output CSV
# =============================================================================
# Every row in cleaned_metadata.csv represents one tissue slide.
# These are the ten columns we write for every non-STimage source:
#   slide        : unique ID for this slide / sample
#   tech         : sequencing technology (ST or Visium)
#   title        : full paper title (so Assemble.py can populate HEST metadata)
#   tissue       : organ of origin (e.g. "Prostate", "Lung")
#   species      : "human" or "mouse"
#   involve_cancer : "True" or "False"
#   pmid         : PubMed ID of the paper (links to the publication)
#   gene_num     : number of genes in the panel (blank if unknown)
#   source       : internal tag telling Assemble.py how to download this slide
#   data_path    : the dataset-level identifier (e.g. the GEO accession or Zenodo URL)
FIELDS = ["slide", "tech", "title", "tissue", "species",
          "involve_cancer", "pmid", "gene_num", "source", "data_path"]

# Small helper: write an empty ambiguous file (just the header, no data rows).
# This is needed because Assemble.py expects the file to always exist.
def _write_empty_ambiguous(fieldnames):
    with open(AMBIGUOUS_PATH, "w", newline="") as f:
        csv.DictWriter(f, fieldnames=fieldnames).writeheader()

# =============================================================================
# STEP 6 — Hardcoded metadata dictionaries
# =============================================================================
# For GEO, EBI, TNBC, and TLS we know in advance (from reading the papers)
# what tissue type, species, and cancer status each dataset has.
# Rather than scraping this from the web every time, we store it here as
# Python dictionaries.  Each entry maps an accession code to its metadata.
#
# WHERE DID THESE VALUES COME FROM?
#   - tissue / species / involve_cancer: stated in the abstract of each paper
#   - title : copied verbatim from the paper title
#   - pmid  : the PubMed ID found on pubmed.ncbi.nlm.nih.gov
#   - tech  : stated in the methods section (all three use 10x Visium)

# --- GEO datasets ---
_GEO_META = {
    "GSE278936": {
        # Paper: "Single cell and spatial transcriptomics highlight the interaction of
        # club-like cells with immunosuppressive myeloid cells in prostate cancer"
        # Source: https://pubmed.ncbi.nlm.nih.gov/39550375/
        "title":          ("Single cell and spatial transcriptomics highlight the interaction of "
                           "club-like cells with immunosuppressive myeloid cells in prostate cancer"),
        "tissue":         "Prostate",
        "species":        "human",
        "involve_cancer": "True",
        "pmid":           "39550375",
        "gene_num":       "",    # not specified in the paper
        "tech":           "Visium",
    },
}

# --- EBI datasets ---
_EBI_META = {
    "E-MTAB-13530": {
        # Paper: "Single cell and spatial transcriptomics analysis of non-small cell lung cancer"
        # Source: https://pubmed.ncbi.nlm.nih.gov/38782905/
        "title":          "Single cell and spatial transcriptomics analysis of non-small cell lung cancer",
        "tissue":         "Lung",
        "species":        "human",
        "involve_cancer": "True",
        "pmid":           "38782905",
        "gene_num":       "",
        "tech":           "Visium",
    },
}

# --- TNBC (Triple-Negative Breast Cancer) ---
# Paper: "Spatial transcriptomics reveals substantial heterogeneity in
# triple-negative breast cancer with potential clinical implications"
# Source: https://pubmed.ncbi.nlm.nih.gov/39562547/
# Note: this dataset uses the older "ST" technology (not Visium).
_TNBC_META = {
    "title":          ("Spatial transcriptomics reveals substantial heterogeneity in "
                       "triple-negative breast cancer with potential clinical implications"),
    "tissue":         "Breast",
    "species":        "human",
    "involve_cancer": "True",
    "pmid":           "39562547",
    "gene_num":       "",
    "tech":           "ST",    # Spatial Transcriptomics (older technology, square grid spots)
}

# --- TLS Visium (Tumour-Lymphocyte Structures, kidney + lung) ---
# Paper: "DeepSpot: Leveraging Spatial Context for Enhanced Spatial
# Transcriptomics Prediction from H&E Images"
# Note: tissue type differs per sample — KC = kidney cancer, LC = lung cancer.
#       We assign tissue per-sample in the loop below, not here.
_TLS_META_BASE = {
    "title":          ("DeepSpot: Leveraging Spatial Context for Enhanced Spatial "
                       "Transcriptomics Prediction from H&E Images"),
    "species":        "human",
    "involve_cancer": "True",
    "pmid":           "",       # preprint, no PubMed ID yet
    "gene_num":       "",
    "tech":           "Visium",
}
# The 8 sample names are hardcoded because the Zenodo archive is a fixed,
# known set of files.  KC = Kidney Cancer, LC = Lung Cancer.
_TLS_SAMPLES = ["KC1", "KC2", "KC3", "LC1", "LC2", "LC3", "LC4", "LC5"]

# =============================================================================
# BRANCH A — Zenodo URL (TNBC or TLS)
# =============================================================================
# If the user passed a zenodo.org URL, figure out which record it is and
# handle it accordingly.
if is_zenodo:
    try:
        import requests
    except ImportError:
        print("ERROR: pip install requests"); sys.exit(1)

    # Extract the numeric record ID from the end of the URL.
    # e.g. "https://zenodo.org/records/14204217" → "14204217"
    record_id = ARG.rstrip("/").split("/")[-1]

    # ── TNBC (Zenodo record 14204217) ────────────────────────────────────────
    if record_id == TNBC_RECORD:
        # The TNBC Zenodo archive contains a file called Clinical.tar.
        # Inside that tar is an R data file called ids.RDS which lists the
        # unique array (slide) IDs.  We download Clinical.tar, stream it,
        # extract ids.RDS, and read it with pyreadr (an R-to-Python converter).
        # This avoids downloading the full ~50 GB dataset just to get the names.
        try:
            import pyreadr
        except ImportError:
            print("ERROR: pip install pyreadr"); sys.exit(1)

        clinical_url = f"https://zenodo.org/records/{record_id}/files/Clinical.tar"
        print(f"Fetching array IDs from Clinical.tar (Zenodo record {record_id})…")

        ids_bytes = None
        # stream=True means we receive the file in chunks rather than all at once —
        # important because Clinical.tar is large.
        with requests.get(clinical_url, stream=True, timeout=60) as resp:
            resp.raise_for_status()
            resp.raw.decode_content = True
            # Open the tar archive directly from the network stream (mode="r|"
            # means "streaming, uncompressed").  We scan every file inside it
            # until we find ids.RDS, read its bytes, then stop.
            with tarfile.open(fileobj=resp.raw, mode="r|") as tar:
                for m in tar:
                    if m.name.endswith("ids.RDS") and m.isfile():
                        ids_bytes = tar.extractfile(m).read()
                        break

        if ids_bytes is None:
            print("ERROR: ids.RDS not found in Clinical.tar"); sys.exit(1)

        # Write ids.RDS to a temporary file on disk so pyreadr can read it
        # (pyreadr needs a real file path, not just bytes in memory).
        with tempfile.NamedTemporaryFile(suffix=".RDS", delete=False) as tmp:
            tmp.write(ids_bytes); tmp_path = tmp.name
        try:
            # pyreadr.read_r() reads an R .RDS or .RData file into a Python dict.
            # The dict values are pandas DataFrames; we take the first (and only) one.
            ids_df = list(pyreadr.read_r(tmp_path).values())[0]
        finally:
            os.unlink(tmp_path)   # delete the temp file immediately after reading

        # The index of the DataFrame contains the array IDs (e.g. "CN1_C1").
        array_ids = list(ids_df.index)

        # Write one row per array ID into cleaned_metadata.csv.
        # {**_TNBC_META, "slide": array_id, ...} merges the shared metadata
        # dict with the per-slide fields.
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

    # ── TLS Visium (Zenodo record 14620362) ──────────────────────────────────
    elif record_id == TLS_RECORD:
        # For TLS we already know all 8 sample names, so no web API call needed.
        # We just write one row per sample, assigning tissue based on the prefix:
        #   KC* → Kidney Cancer
        #   LC* → Lung Cancer
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

# =============================================================================
# BRANCH B — GEO accession (e.g. GSE278936)
# =============================================================================
# NCBI GEO (Gene Expression Omnibus) is the main US repository for genomics
# data.  A "GSE" (GEO Series) groups many individual "GSM" (GEO Sample) entries.
# We call the GEO web API to get the list of all GSM sample IDs in the series.
if is_geo:
    try:
        import requests
    except ImportError:
        print("ERROR: pip install requests"); sys.exit(1)

    geo_acc = ARG.upper()
    print(f"Fetching sample list from GEO {geo_acc}…")

    # The GEO text API returns a plain-text summary of the series.
    # Parameters:
    #   acc=GSE278936   → the series we want
    #   targ=gsm        → include all GSM samples
    #   view=brief      → minimal output (faster)
    #   form=text       → return plain text, not HTML
    r = requests.get(
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
        params={"acc": geo_acc, "targ": "gsm", "view": "brief", "form": "text"},
        timeout=60)
    r.raise_for_status()

    # The text response contains lines like:
    #   ^SAMPLE = GSM8557976
    # We extract every GSM number using a regular expression.
    # The leading ^ is literal (part of the NCBI format, not a regex anchor),
    # so we write \^ to match it literally.
    gsm_ids = re.findall(r'^\^SAMPLE = (GSM\d+)', r.text, re.MULTILINE)

    if not gsm_ids:
        print(f"ERROR: No GSM samples found for {geo_acc}"); sys.exit(1)

    # Look up the hardcoded metadata for this accession, or use a generic
    # fallback if we haven't seen it before.
    meta = _GEO_META.get(geo_acc, {
        "title": geo_acc, "tissue": "", "species": "human",
        "involve_cancer": "True", "pmid": "", "gene_num": "", "tech": "Visium",
    })

    # Write one row per GSM sample ID.
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

# =============================================================================
# BRANCH C — EBI accession (e.g. E-MTAB-13530)
# =============================================================================
# EBI BioStudies / ArrayExpress is the main European repository for genomics
# data.  We call the BioStudies JSON API to get the file listing for the study,
# then extract sample names from it.
if is_ebi:
    try:
        import requests
    except ImportError:
        print("ERROR: pip install requests"); sys.exit(1)

    ebi_acc = ARG.upper()
    print(f"Fetching sample list from EBI {ebi_acc}…")

    # The BioStudies API returns a JSON document describing the study,
    # including a nested list of every file it contains.
    r = requests.get(
        f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{ebi_acc}",
        timeout=60)
    r.raise_for_status()

    # _collect_paths() recursively walks the JSON tree and collects every
    # "path" field it finds (these are relative file paths within the study).
    def _collect_paths(obj):
        paths = []
        if isinstance(obj, dict):
            # If this dict node represents a file, grab its path.
            if obj.get("type") == "file" and "path" in obj:
                paths.append(obj["path"])
            # Then recurse into all child values.
            for v in obj.values():
                paths.extend(_collect_paths(v))
        elif isinstance(obj, list):
            for item in obj:
                paths.extend(_collect_paths(item))
        return paths

    all_paths = _collect_paths(r.json())

    # Every sample in E-MTAB-13530 has a file called:
    #   {sample_name}-filtered_feature_bc_matrix.h5
    # We find all such files and strip the suffix to recover the sample name.
    # dict.fromkeys() preserves order and removes duplicates.
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

    # Write one row per sample name.
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

# =============================================================================
# BRANCH D — Local TNBC directory
# =============================================================================
# If the user downloaded the TNBC Zenodo archive manually and unpacked it,
# they can pass the local directory path.  We scan the byArray/ subfolder
# to get the list of array IDs from the directory names.
if os.path.isdir(ARG):
    by_array = os.path.join(ARG, "byArray")
    if not os.path.isdir(by_array):
        print(f"ERROR: {ARG} is a directory but has no byArray/ subdirectory."); sys.exit(1)

    if os.path.exists(CLEANED_PATH):
        # If we already ran this before, don't overwrite the output.
        print(f"\nExisting output detected — skipping TNBC metadata generation.")
    else:
        # List every subdirectory inside byArray/ — each is one tissue array.
        arrays = sorted(
            d for d in os.listdir(by_array)
            if os.path.isdir(os.path.join(by_array, d))
        )
        with open(CLEANED_PATH, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            w.writeheader()
            for array_id in arrays:
                # source="tnbc" (not "tnbc_zenodo") tells Assemble.py to read
                # from the local directory instead of downloading from Zenodo.
                w.writerow({**_TNBC_META, "slide": array_id,
                            "source": "tnbc", "data_path": ARG})
        print(f"TNBC: wrote {len(arrays)} arrays → {CLEANED_PATH}")

    _write_empty_ambiguous(FIELDS)
    print(f"Ready: python Assemble.py {CLEANED_PATH}")
    sys.exit(0)

# =============================================================================
# BRANCH E — STimage-1K4M CSV
# =============================================================================
# This is the largest and most complex branch.  The STimage-1K4M dataset
# contains 1,149 tissue slides from many different published studies.
# Some of those studies are already in HEST, so we must cross-check and
# remove duplicates before converting.
#
# The input CSV (meta_all_gene.csv) has one row per slide with columns:
#   slide, tech, title, tissue, species, involve_cancer, pmid, gene_num, ...
#
# We classify each slide into one of three groups:
#   1. CLEANED  — the study is NOT in HEST at all → safe to convert
#   2. DUPLICATE — the study IS in HEST and we can confirm this specific
#                  slide is in HEST → skip it
#   3. AMBIGUOUS — the study IS in HEST but we can't confirm whether this
#                  specific slide is included → flag for human review

STIMAGE_META = ARG
# Remember whether ARG was a downloaded temp file (we'll delete it at the end).
_stimage_tmp = ARG if is_csv_url else None

if os.path.exists(CLEANED_PATH) and os.path.exists(AMBIGUOUS_PATH):
    # Both output files already exist — skip the cleaning step and go
    # straight to the human review gate below.
    print(f"\nExisting output detected — skipping cleaning step.")
else:
    print("Cleaning…")

    # Read the HEST master catalogue into memory.
    # This file maps study titles to subseries labels (short identifiers HEST
    # uses to name individual slides within a multi-sample study).
    with open(HEST_META, encoding="utf-8-sig") as f:
        hest = list(csv.DictReader(f))

    # Read the STimage slide list.
    with open(STIMAGE_META) as f:
        stimage = list(csv.DictReader(f))

    # Build a lookup table: lowercase study title → list of subseries labels.
    # Example: "visium of breast cancer" → ["CN1", "CN2", "CN3"]
    # We lowercase titles so comparison is case-insensitive.
    hest_title_to_subseries: dict[str, list[str]] = {}
    for row in hest:
        title = row["dataset_title"].strip().lower()
        sub   = row["subseries"].strip()
        hest_title_to_subseries.setdefault(title, [])
        if sub:
            hest_title_to_subseries[title].append(sub)

    # Some STimage rows reference multiple papers in a single "title" field,
    # formatted like: "Title 1: First paper. Title 2: Second paper."
    # This function extracts each individual title from that combined string.
    def extract_titles(raw: str) -> list[str]:
        titles = re.findall(r"Title \d+: (.+?)(?= Title \d+:|$)", raw)
        return [t.strip() for t in titles] if titles else [raw.strip()]

    cleaned, ambiguous, duplicates, non_cancer = [], [], [], []

    for row in stimage:
        # FILTER 1: Skip non-cancer slides.
        # We only want slides where involve_cancer is True/1/yes.
        if str(row.get("involve_cancer", "")).strip().lower() not in ("true", "1", "yes"):
            non_cancer.append(row)
            continue

        slide  = row["slide"]
        titles = extract_titles(row["title"])

        # FILTER 2: Check if this study's title appears in the HEST catalogue.
        # We check all titles in the field (in case there are multiple).
        matched_title = None
        for t in titles:
            if t.lower() in hest_title_to_subseries:
                matched_title = t.lower()
                break

        # If the study is NOT in HEST at all → it's new, safe to convert.
        if matched_title is None:
            cleaned.append(row)
            continue

        # The study IS in HEST.  Now check whether this specific slide is
        # already in HEST by looking for a subseries label inside the slide ID.
        # Example: HEST says subseries for this study are ["CN1", "CN2"].
        #          If the slide ID contains "CN1" → confirmed duplicate.
        subseries_list = hest_title_to_subseries[matched_title]
        confirmed = any(sub and sub in slide for sub in subseries_list)

        if confirmed:
            # This exact slide is already in HEST — drop it silently.
            duplicates.append(row)
        else:
            # The study is in HEST but we can't match this slide to a known
            # subseries.  Might be a different naming convention.
            # Flag it for human review with extra debug columns.
            row["_matched_hest_title"] = matched_title
            row["_hest_subseries"]     = " | ".join(subseries_list[:10])
            ambiguous.append(row)

    # Get the column names from the first row of the input.
    fieldnames = list(stimage[0].keys())

    # Write the safe-to-convert slides.
    with open(CLEANED_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(cleaned)

    # Write the ambiguous slides (with two extra debug columns appended).
    amb_fields = fieldnames + ["_matched_hest_title", "_hest_subseries"]
    with open(AMBIGUOUS_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=amb_fields, extrasaction="ignore")
        w.writeheader(); w.writerows(ambiguous)

    print(f"Total STimage slides : {len(stimage)}")
    print(f"Non-cancer           : {len(non_cancer)}  (skipped)")
    print(f"Confirmed duplicates : {len(duplicates)}  (skipped)")
    print(f"Ambiguous            : {len(ambiguous)}   → ambiguous_metadata.csv")
    print(f"Safe to convert      : {len(cleaned)}  → cleaned_metadata.csv")

# =============================================================================
# STEP 7 — Human review gate (STimage only)
# =============================================================================
# The ambiguous slides might or might not be duplicates.  A human needs to
# open ambiguous_metadata.csv, look at each row, and delete the ones that
# really are in HEST.  Then we merge the survivors into cleaned_metadata.csv.
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
    # Read the (possibly edited) ambiguous file and append it to cleaned.
    with open(AMBIGUOUS_PATH, newline="") as f:
        amb_rows = list(csv.DictReader(f))
    with open(CLEANED_PATH, newline="") as f:
        cleaned_fields = csv.DictReader(f).fieldnames
    with open(CLEANED_PATH, "a", newline="") as f:
        # extrasaction="ignore" drops the two debug columns we added earlier.
        w = csv.DictWriter(f, fieldnames=cleaned_fields, extrasaction="ignore")
        w.writerows(amb_rows)
    total = sum(1 for _ in open(CLEANED_PATH)) - 1   # subtract 1 for header
    print(f"Appended {len(amb_rows)} rows. Total safe-to-convert slides: {total}")
else:
    print("Ambiguous rows not added. cleaned_metadata.csv unchanged.")

# Clean up the temporary file if we downloaded the CSV from a URL.
if _stimage_tmp and os.path.exists(_stimage_tmp):
    os.unlink(_stimage_tmp)
