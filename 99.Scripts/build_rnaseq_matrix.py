"""
80개 STAR-Counts TSV -> 환자 x 유전자 행렬 생성
출력:
  rnaseq_tpm_matrix.csv        (80 x ~19k, protein_coding, TPM)
  rnaseq_counts_matrix.csv     (80 x ~19k, protein_coding, raw counts)
  rnaseq_tpm_top5000.csv       (80 x 5000, 분산 상위 유전자, TPM)
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import pandas as pd
import numpy as np
import os, glob, requests, json

RNASEQ_DIR = r"C:\Users\RokitGenomics\.openclaw\workspace\UVM_image_diagnosis\01.Raw_Data\03.Genomics\01.RNAseq"
OUT_DIR    = r"C:\Users\RokitGenomics\.openclaw\workspace\UVM_image_diagnosis\01.Raw_Data\03.Genomics"

# ── Step 1: file_id → patient_id 매핑 ───────────────────────────────
print("Step 1: GDC API에서 file→patient 매핑 가져오기...")
resp = requests.get(
    "https://api.gdc.cancer.gov/files",
    params={
        "filters": json.dumps({
            "op": "and", "content": [
                {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-UVM"}},
                {"op": "=", "content": {"field": "data_type",                "value": "Gene Expression Quantification"}},
                {"op": "=", "content": {"field": "analysis.workflow_type",   "value": "STAR - Counts"}},
            ]
        }),
        "fields": "file_id,file_name,cases.submitter_id",
        "size": 200,
    }, timeout=30
)
hits = resp.json()["data"]["hits"]
# file_name의 첫 부분 (UUID)이 폴더명
file2patient = {h["file_id"]: h["cases"][0]["submitter_id"] for h in hits}
fname2patient = {h["file_name"]: h["cases"][0]["submitter_id"] for h in hits}
print(f"  매핑 완료: {len(file2patient)}개")

# ── Step 2: TSV 파일 목록 ────────────────────────────────────────────
tsv_files = glob.glob(os.path.join(RNASEQ_DIR, "**", "*.tsv"), recursive=True)
tsv_files = [f for f in tsv_files if "gene_counts" in f]
print(f"\nStep 2: TSV 파일 {len(tsv_files)}개 발견")

# ── Step 3: 전체 행렬 구축 ───────────────────────────────────────────
print("\nStep 3: 행렬 구축 중...")
tpm_dict    = {}
counts_dict = {}

for i, fpath in enumerate(tsv_files):
    fname = os.path.basename(fpath)
    patient_id = fname2patient.get(fname)
    if patient_id is None:
        # 폴더명(UUID)으로 매핑 시도
        folder_uuid = os.path.basename(os.path.dirname(fpath))
        patient_id = file2patient.get(folder_uuid, f"UNKNOWN_{i}")

    df = pd.read_csv(fpath, sep='\t', comment='#', low_memory=False)
    # ENSG 행만 & protein_coding 유전자만
    df = df[df['gene_id'].astype(str).str.startswith('ENSG')]
    df = df[df['gene_type'] == 'protein_coding'].copy()
    # gene_id에서 버전 제거 (ENSG00000000003.15 -> ENSG00000000003)
    df['gene_id_clean'] = df['gene_id'].str.split('.').str[0]
    df = df.set_index('gene_id_clean')

    tpm_dict[patient_id]    = df['tpm_unstranded']
    counts_dict[patient_id] = df['unstranded']

    if (i+1) % 20 == 0:
        print(f"  {i+1}/{len(tsv_files)} 처리 완료...")

print(f"  전체 {len(tsv_files)}개 처리 완료!")

# ── Step 4: DataFrame 변환 (환자 x 유전자) ───────────────────────────
print("\nStep 4: 행렬 변환 중...")
tpm_df    = pd.DataFrame(tpm_dict).T      # 환자 x 유전자
counts_df = pd.DataFrame(counts_dict).T

tpm_df.index.name    = 'patient_id'
counts_df.index.name = 'patient_id'
tpm_df    = tpm_df.sort_index()
counts_df = counts_df.sort_index()

print(f"  TPM 행렬: {tpm_df.shape}  (환자 x 유전자)")
print(f"  Counts 행렬: {counts_df.shape}")

# gene_name 매핑 추가 (컬럼 이름을 gene_name으로 변경)
first_file = pd.read_csv(tsv_files[0], sep='\t', comment='#', low_memory=False)
first_file = first_file[first_file['gene_id'].astype(str).str.startswith('ENSG')]
first_file = first_file[first_file['gene_type'] == 'protein_coding']
first_file['gene_id_clean'] = first_file['gene_id'].str.split('.').str[0]
id2name = dict(zip(first_file['gene_id_clean'], first_file['gene_name']))

# ── Step 5: 분산 기준 상위 5000 유전자 선택 ─────────────────────────
print("\nStep 5: 분산 상위 유전자 선택...")
# 최소 50% 환자에서 TPM > 0인 유전자만
expressed = (tpm_df > 0).mean(axis=0) >= 0.5
tpm_filtered = tpm_df.loc[:, expressed]
print(f"  발현된 유전자 (>=50% 환자): {expressed.sum()}개")

gene_var = tpm_filtered.var(axis=0)
top5000_genes = gene_var.nlargest(5000).index
tpm_top5000 = tpm_filtered[top5000_genes].copy()
print(f"  분산 상위 5000개 선택 완료")

# ── Step 6: 저장 ─────────────────────────────────────────────────────
print("\nStep 6: CSV 저장 중...")

tpm_df.to_csv(os.path.join(OUT_DIR, "rnaseq_tpm_matrix.csv"), encoding='utf-8-sig')
print(f"  rnaseq_tpm_matrix.csv 저장 완료  ({tpm_df.shape[0]}x{tpm_df.shape[1]})")

counts_df.to_csv(os.path.join(OUT_DIR, "rnaseq_counts_matrix.csv"), encoding='utf-8-sig')
print(f"  rnaseq_counts_matrix.csv 저장 완료")

tpm_top5000.to_csv(os.path.join(OUT_DIR, "rnaseq_tpm_top5000.csv"), encoding='utf-8-sig')
print(f"  rnaseq_tpm_top5000.csv 저장 완료  ({tpm_top5000.shape[0]}x{tpm_top5000.shape[1]})")

# ── Step 7: 간단 요약 ────────────────────────────────────────────────
print("\n" + "="*50)
print("완료! 요약")
print("="*50)
print(f"  환자 수:        {tpm_df.shape[0]}명")
print(f"  전체 유전자:    {tpm_df.shape[1]}개 (protein_coding)")
print(f"  발현 유전자:    {expressed.sum()}개 (>=50% 환자)")
print(f"  Top 5000 유전자: {tpm_top5000.shape[1]}개")
print(f"\n  TPM 값 범위 (전체): {tpm_df.values.min():.2f} ~ {tpm_df.values.max():.2f}")
print(f"  중앙값: {np.median(tpm_df.values[tpm_df.values > 0]):.2f}")

# 상위 발현 유전자 (전체 평균 기준)
top_expressed = tpm_df.mean(axis=0).nlargest(10)
print(f"\n  평균 발현 상위 10개 유전자:")
for gid, val in top_expressed.items():
    gname = id2name.get(gid, gid)
    print(f"    {gname:15s}  {val:.1f} TPM")
