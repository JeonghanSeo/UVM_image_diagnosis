"""
TCGA-UVM RNA-seq (STAR-Counts) GDC 다운로드 스크립트
사용법:
  1. manifest 생성:  python download_rnaseq.py --manifest
  2. 다운로드 실행:  python download_rnaseq.py --download
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import requests
import json
import os
import argparse

BASE       = "https://api.gdc.cancer.gov"
OUT_DIR    = r"C:\Users\RokitGenomics\.openclaw\workspace\UVM_image_diagnosis\01.Raw_Data\03.Genomics\01.RNAseq"
MANIFEST   = os.path.join(OUT_DIR, "gdc_manifest_rnaseq.txt")
GDC_CLIENT = r"C:\Users\RokitGenomics\.openclaw\workspace\UVM_image_diagnosis\99.Scripts\gdc-client\gdc-client.exe"

os.makedirs(OUT_DIR, exist_ok=True)

# ── GDC 필터 ─────────────────────────────────────────────────────────
FILTERS = {
    "op": "and",
    "content": [
        {"op": "=", "content": {"field": "cases.project.project_id",   "value": "TCGA-UVM"}},
        {"op": "=", "content": {"field": "data_category",              "value": "Transcriptome Profiling"}},
        {"op": "=", "content": {"field": "data_type",                  "value": "Gene Expression Quantification"}},
        {"op": "=", "content": {"field": "experimental_strategy",      "value": "RNA-Seq"}},
        {"op": "=", "content": {"field": "analysis.workflow_type",     "value": "STAR - Counts"}},
    ]
}

FIELDS = ",".join([
    "file_id", "file_name", "md5sum", "file_size", "state",
    "cases.submitter_id", "cases.case_id",
])


def get_file_list():
    params = {
        "filters": json.dumps(FILTERS),
        "fields":  FIELDS,
        "size":    200,
        "format":  "json",
    }
    resp = requests.get(f"{BASE}/files", params=params, timeout=30)
    resp.raise_for_status()
    hits = resp.json()["data"]["hits"]
    print(f"Found {len(hits)} files")
    return hits


def make_manifest(hits):
    with open(MANIFEST, "w", encoding="utf-8") as f:
        f.write("id\tfilename\tmd5\tsize\tstate\n")
        for h in hits:
            f.write(f"{h['file_id']}\t{h['file_name']}\t{h.get('md5sum','')}\t{h.get('file_size',0)}\treleased\n")
    print(f"Manifest saved: {MANIFEST}")
    print(f"Total files: {len(hits)}")

    total_gb = sum(h.get('file_size', 0) for h in hits) / 1e9
    print(f"Total size: {total_gb:.2f} GB")

    # 파일 목록 미리보기
    print("\nSample file list (first 5):")
    for h in hits[:5]:
        case_id = h['cases'][0]['submitter_id'] if h.get('cases') else 'unknown'
        size_mb = h.get('file_size', 0) / 1e6
        print(f"  {case_id}  {h['file_name']}  ({size_mb:.1f} MB)")


def download(n_threads=4):
    if not os.path.exists(MANIFEST):
        print("Manifest not found. Run with --manifest first.")
        return
    if not os.path.exists(GDC_CLIENT):
        print(f"gdc-client not found: {GDC_CLIENT}")
        return

    print(f"Downloading to: {OUT_DIR}")
    print(f"Threads: {n_threads}")
    print("Running gdc-client...")
    os.system(
        f'"{GDC_CLIENT}" download '
        f'--manifest "{MANIFEST}" '
        f'--dir "{OUT_DIR}" '
        f'--n-processes {n_threads}'
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", action="store_true", help="manifest 파일 생성")
    parser.add_argument("--download", action="store_true", help="다운로드 실행")
    parser.add_argument("--threads",  type=int, default=4)
    args = parser.parse_args()

    if args.manifest or (not args.manifest and not args.download):
        hits = get_file_list()
        make_manifest(hits)

    if args.download:
        download(args.threads)
