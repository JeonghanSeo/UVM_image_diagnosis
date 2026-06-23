import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
import warnings
warnings.filterwarnings('ignore')

# ── paths ──────────────────────────────────────────────────────────
DATA  = r'C:\Users\RokitGenomics\.openclaw\workspace\UVM_image_diagnosis\01.Raw_Data\02.Clinical\master_clinical_table.csv'
FIGS  = r'C:\Users\RokitGenomics\.openclaw\workspace\UVM_image_diagnosis\07.Figures'

# ── load ────────────────────────────────────────────────────────────
df = pd.read_csv(DATA)
df['scna_cluster'] = df['scna_cluster'].astype(int)
df['OS_months']  = df['OS_time_days']  / 30.44
df['PFI_months'] = df['PFI_time_days'] / 30.44

CLUSTER_COLORS = {1: '#2196F3', 2: '#4CAF50', 3: '#FF9800', 4: '#F44336'}
CLUSTER_LABELS = {
    1: 'Cluster 1\n(D3, EIF1AX)',
    2: 'Cluster 2\n(D3, SF3B1)',
    3: 'Cluster 3\n(M3, BAP1)',
    4: 'Cluster 4\n(M3, BAP1+iso8q)',
}

print(f'Loaded: {df.shape[0]} patients x {df.shape[1]} columns')
print(f'SCNA clusters: {dict(df["scna_cluster"].value_counts().sort_index())}')

# ════════════════════════════════════════════════════════════════════
# Fig 1: SCNA Cluster distribution (pie + bar)
# ════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('TCGA-UVM: SCNA Cluster Distribution (n=80)', fontsize=14, fontweight='bold')

counts = df['scna_cluster'].value_counts().sort_index()
colors = [CLUSTER_COLORS[c] for c in counts.index]
labels_pie = [f'Cluster {c}\n(n={n})' for c, n in zip(counts.index, counts.values)]

axes[0].pie(counts.values, labels=labels_pie, colors=colors,
            autopct='%1.1f%%', startangle=90,
            textprops={'fontsize': 10})
axes[0].set_title('Proportion', fontsize=12)

bars = axes[1].bar([f'C{c}' for c in counts.index], counts.values, color=colors, edgecolor='white', linewidth=1.5)
for bar, val in zip(bars, counts.values):
    axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                 str(val), ha='center', fontsize=11, fontweight='bold')
axes[1].set_xlabel('SCNA Cluster', fontsize=11)
axes[1].set_ylabel('Number of Patients', fontsize=11)
axes[1].set_title('Count', fontsize=12)
axes[1].set_ylim(0, 28)
axes[1].set_xticks(range(4))
axes[1].set_xticklabels([f'C{c}\n({["D3,EIF1AX","D3,SF3B1","M3,BAP1","M3,BAP1+8q"][i]})' for i, c in enumerate(counts.index)], fontsize=8)

plt.tight_layout()
plt.savefig(f'{FIGS}\\fig1_cluster_distribution.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: fig1_cluster_distribution.png')

# ════════════════════════════════════════════════════════════════════
# Fig 2: Kaplan-Meier OS by SCNA Cluster
# ════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('TCGA-UVM: Kaplan-Meier Survival by SCNA Cluster', fontsize=14, fontweight='bold')

for ax, (time_col, event_col, title) in zip(axes, [
    ('OS_months',  'OS',  'Overall Survival (OS)'),
    ('PFI_months', 'PFI', 'Progression-Free Interval (PFI)'),
]):
    kmf = KaplanMeierFitter()
    for cluster in [1, 2, 3, 4]:
        mask = df['scna_cluster'] == cluster
        T = df.loc[mask, time_col].dropna()
        E = df.loc[mask, event_col].dropna()
        idx = T.index.intersection(E.index)
        kmf.fit(T[idx], E[idx], label=f'Cluster {cluster} (n={mask.sum()})')
        kmf.plot_survival_function(ax=ax, ci_show=True, color=CLUSTER_COLORS[cluster], linewidth=2)

    # Log-rank test (4-group)
    valid = df[[time_col, event_col, 'scna_cluster']].dropna()
    res = multivariate_logrank_test(
        valid[time_col], valid['scna_cluster'], valid[event_col]
    )
    ax.set_title(f'{title}\nLog-rank p = {res.p_value:.4f}', fontsize=12)
    ax.set_xlabel('Time (months)', fontsize=11)
    ax.set_ylabel('Survival Probability', fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=9, loc='lower left')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{FIGS}\\fig2_kaplan_meier.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: fig2_kaplan_meier.png')

# ════════════════════════════════════════════════════════════════════
# Fig 3: Clinical features by SCNA Cluster
# ════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle('TCGA-UVM: Clinical Features by SCNA Cluster', fontsize=14, fontweight='bold')

# 3-1: Age distribution
ax = axes[0, 0]
data_by_cluster = [df[df['scna_cluster'] == c]['age_at_diagnosis'].dropna().values for c in [1,2,3,4]]
bp = ax.boxplot(data_by_cluster, patch_artist=True,
                medianprops=dict(color='black', linewidth=2))
for patch, c in zip(bp['boxes'], [1,2,3,4]):
    patch.set_facecolor(CLUSTER_COLORS[c])
    patch.set_alpha(0.7)
ax.set_xticklabels([f'C{c}' for c in [1,2,3,4]])
ax.set_xlabel('SCNA Cluster')
ax.set_ylabel('Age at Diagnosis (years)')
ax.set_title('Age Distribution')
ax.grid(True, alpha=0.3, axis='y')

# 3-2: Gender distribution
ax = axes[0, 1]
gender_ct = pd.crosstab(df['scna_cluster'], df['gender'])
gender_pct = gender_ct.div(gender_ct.sum(axis=1), axis=0) * 100
bottom = np.zeros(4)
for gender, color in [('male', '#5C9BD6'), ('female', '#E8A0A0')]:
    if gender in gender_pct.columns:
        vals = gender_pct[gender].values
        ax.bar([f'C{c}' for c in [1,2,3,4]], vals, bottom=bottom,
               label=gender.capitalize(), color=color, alpha=0.8)
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v > 5:
                ax.text(i, b + v/2, f'{v:.0f}%', ha='center', va='center', fontsize=9, color='white', fontweight='bold')
        bottom += vals
ax.set_xlabel('SCNA Cluster')
ax.set_ylabel('Percentage (%)')
ax.set_title('Gender Distribution')
ax.legend(loc='upper right', fontsize=9)
ax.set_ylim(0, 110)

# 3-3: AJCC Stage distribution
ax = axes[1, 0]
stage_map = {}
for s in df['ajcc_stage'].dropna().unique():
    s_str = str(s)
    if 'IV' in s_str: stage_map[s] = 'Stage IV'
    elif 'IIIC' in s_str: stage_map[s] = 'Stage IIIC'
    elif 'IIIB' in s_str: stage_map[s] = 'Stage IIIB'
    elif 'IIIA' in s_str: stage_map[s] = 'Stage IIIA'
    elif 'IIC' in s_str: stage_map[s] = 'Stage IIC'
    elif 'IIB' in s_str: stage_map[s] = 'Stage IIB'
    elif 'IIA' in s_str: stage_map[s] = 'Stage IIA'
    elif 'I' in s_str: stage_map[s] = 'Stage I'
    else: stage_map[s] = s_str
df['ajcc_stage_grouped'] = df['ajcc_stage'].map(stage_map)

stage_ct = pd.crosstab(df['scna_cluster'], df['ajcc_stage_grouped'])
stage_colors = ['#B39DDB', '#7986CB', '#42A5F5', '#26A69A', '#66BB6A', '#FFA726', '#EF5350']
stage_ct.plot(kind='bar', ax=ax, color=stage_colors[:len(stage_ct.columns)], alpha=0.8, edgecolor='white')
ax.set_xticklabels([f'C{c}' for c in stage_ct.index], rotation=0)
ax.set_xlabel('SCNA Cluster')
ax.set_ylabel('Number of Patients')
ax.set_title('AJCC Stage Distribution')
ax.legend(fontsize=7, loc='upper right', ncol=2)
ax.grid(True, alpha=0.3, axis='y')

# 3-4: OS event rate by cluster
ax = axes[1, 1]
os_rates = df.groupby('scna_cluster')['OS'].mean() * 100
bar_colors = [CLUSTER_COLORS[c] for c in os_rates.index]
bars = ax.bar([f'C{c}' for c in os_rates.index], os_rates.values, color=bar_colors, alpha=0.8, edgecolor='white', linewidth=1.5)
for bar, val in zip(bars, os_rates.values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f'{val:.1f}%', ha='center', fontsize=10, fontweight='bold')
ax.set_xlabel('SCNA Cluster')
ax.set_ylabel('Death Rate (%)')
ax.set_title('OS Event Rate by Cluster')
ax.set_ylim(0, 80)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(f'{FIGS}\\fig3_clinical_by_cluster.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: fig3_clinical_by_cluster.png')

# ════════════════════════════════════════════════════════════════════
# Fig 4: Mutation landscape
# ════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('TCGA-UVM: Mutation Landscape by SCNA Cluster', fontsize=14, fontweight='bold')

mut_cols = {
    'BAP1_mutation_present': 'BAP1',
    'EIF1AX_mutation_present': 'EIF1AX',
    'SF3B1_mutation_present': 'SF3B1',
    'GNAQ_mutation_present': 'GNAQ',
    'GNA11_mutation_present': 'GNA11',
}
available_mut = {k: v for k, v in mut_cols.items() if k in df.columns}

# 4-1: Mutation rate heatmap
mut_rates = {}
for cluster in [1, 2, 3, 4]:
    mask = df['scna_cluster'] == cluster
    rates = {}
    for col, label in available_mut.items():
        rates[label] = df.loc[mask, col].mean() * 100
    mut_rates[f'C{cluster}'] = rates

mut_df = pd.DataFrame(mut_rates).T
sns.heatmap(mut_df, ax=axes[0], annot=True, fmt='.0f', cmap='RdYlGn',
            linewidths=0.5, linecolor='white',
            cbar_kws={'label': 'Mutation Rate (%)'})
axes[0].set_title('Mutation Rate (%) per Cluster', fontsize=12)
axes[0].set_ylabel('SCNA Cluster')

# 4-2: Stacked bar - mutation frequency
width = 0.6
x = np.arange(4)
cluster_labels_short = ['C1\n(EIF1AX)', 'C2\n(SF3B1)', 'C3\n(BAP1)', 'C4\n(BAP1+8q)']
mut_colors = {'BAP1': '#E53935', 'EIF1AX': '#1E88E5', 'SF3B1': '#43A047',
              'GNAQ': '#FB8C00', 'GNA11': '#8E24AA'}

bottom = np.zeros(4)
for gene in ['BAP1', 'EIF1AX', 'SF3B1', 'GNAQ', 'GNA11']:
    if gene in mut_df.columns:
        vals = [mut_df.loc[f'C{c}', gene] for c in [1,2,3,4]]
        axes[1].bar(x, vals, width, bottom=bottom,
                    label=gene, color=mut_colors.get(gene, 'gray'), alpha=0.85)
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v > 8:
                axes[1].text(x[i], b + v/2, f'{v:.0f}%', ha='center', va='center',
                             fontsize=8, color='white', fontweight='bold')
        bottom += np.array(vals)

axes[1].set_xticks(x)
axes[1].set_xticklabels(cluster_labels_short, fontsize=9)
axes[1].set_ylabel('Mutation Rate (%)')
axes[1].set_title('Stacked Mutation Rate by Cluster', fontsize=12)
axes[1].legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0, fontsize=9)
axes[1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(f'{FIGS}\\fig4_mutation_landscape.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: fig4_mutation_landscape.png')

# ════════════════════════════════════════════════════════════════════
# Fig 5: Survival summary table + D3 vs M3 KM
# ════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('TCGA-UVM: Disomy 3 (D3) vs Monosomy 3 (M3) Survival', fontsize=14, fontweight='bold')

for ax, (time_col, event_col, title) in zip(axes, [
    ('OS_months',  'OS',  'Overall Survival (OS)'),
    ('PFI_months', 'PFI', 'Progression-Free Interval (PFI)'),
]):
    kmf = KaplanMeierFitter()
    for status, color, ls in [('D3', '#2196F3', '-'), ('M3', '#F44336', '--')]:
        mask = df['chromosome3_status'] == status
        T = df.loc[mask, time_col].dropna()
        E = df.loc[mask, event_col].dropna()
        idx = T.index.intersection(E.index)
        n = mask.sum()
        kmf.fit(T[idx], E[idx], label=f'{status} (n={n})')
        kmf.plot_survival_function(ax=ax, ci_show=True, color=color, linestyle=ls, linewidth=2.5)

    # Log-rank D3 vs M3
    d3 = df[df['chromosome3_status'] == 'D3'][[time_col, event_col]].dropna()
    m3 = df[df['chromosome3_status'] == 'M3'][[time_col, event_col]].dropna()
    res = logrank_test(
        d3[time_col], m3[time_col],
        d3[event_col], m3[event_col]
    )
    ax.set_title(f'{title}\nLog-rank p = {res.p_value:.4f}', fontsize=12)
    ax.set_xlabel('Time (months)', fontsize=11)
    ax.set_ylabel('Survival Probability', fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=10, loc='lower left')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{FIGS}\\fig5_d3_vs_m3_survival.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: fig5_d3_vs_m3_survival.png')

# ════════════════════════════════════════════════════════════════════
# Summary statistics
# ════════════════════════════════════════════════════════════════════
print()
print('='*60)
print('SUMMARY STATISTICS')
print('='*60)
print(f'\nTotal patients: {len(df)}')
print(f'  - D3 (Cluster 1+2): {(df["chromosome3_status"]=="D3").sum()} ({(df["chromosome3_status"]=="D3").mean()*100:.1f}%)')
print(f'  - M3 (Cluster 3+4): {(df["chromosome3_status"]=="M3").sum()} ({(df["chromosome3_status"]=="M3").mean()*100:.1f}%)')
print(f'\nAge: mean={df["age_at_diagnosis"].mean():.1f}, std={df["age_at_diagnosis"].std():.1f}, range={df["age_at_diagnosis"].min():.0f}-{df["age_at_diagnosis"].max():.0f}')
print(f'\nGender: {dict(df["gender"].value_counts())}')
print(f'\nOS events: {int(df["OS"].sum())} deaths ({df["OS"].mean()*100:.1f}%)')
print(f'OS time: median={df["OS_months"].median():.1f} months')
print(f'\nPFI events: {int(df["PFI"].sum())} ({df["PFI"].mean()*100:.1f}%)')
print(f'\nAJCC Stage:\n{df["ajcc_stage"].value_counts().to_string()}')
print()

# ════════════════════════════════════════════════════════════════════
# Fig 0: Cluster Summary Table
# ════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(14, 5))
fig.patch.set_facecolor('#F8F9FA')
ax.set_facecolor('#F8F9FA')
ax.axis('off')
fig.suptitle('TCGA-UVM: SCNA Cluster Summary', fontsize=15, fontweight='bold', y=0.98)

# Compute stats per cluster from data
cluster_stats = []
for c in [1, 2, 3, 4]:
    sub = df[df['scna_cluster'] == c]
    os_med = sub['OS_months'].median()
    pfi_med = sub['PFI_months'].median()
    death_rate = sub['OS'].mean() * 100
    cluster_stats.append({
        'cluster': c,
        'n': len(sub),
        'os_med': os_med,
        'pfi_med': pfi_med,
        'death_rate': death_rate,
    })

col_labels = [
    'Cluster', 'Chr3\nStatus', 'Key\nMutation', 'Old\nClass',
    'n', 'Death\nRate', 'Median\nOS', 'Median\nPFI', 'Prognosis'
]

rows = [
    ['Cluster 1', 'Disomy 3\n(D3)', 'EIF1AX\n(point mut)',  'Class 1A',
     f'{cluster_stats[0]["n"]}',
     f'{cluster_stats[0]["death_rate"]:.0f}%',
     f'{cluster_stats[0]["os_med"]:.1f} mo',
     f'{cluster_stats[0]["pfi_med"]:.1f} mo',
     'Best'],
    ['Cluster 2', 'Disomy 3\n(D3)', 'SF3B1\n(splice mut)', 'Class 1B',
     f'{cluster_stats[1]["n"]}',
     f'{cluster_stats[1]["death_rate"]:.0f}%',
     f'{cluster_stats[1]["os_med"]:.1f} mo',
     f'{cluster_stats[1]["pfi_med"]:.1f} mo',
     'Very Good'],
    ['Cluster 3', 'Monosomy 3\n(M3)', 'BAP1\n(del/mut)',   'Class 2',
     f'{cluster_stats[2]["n"]}',
     f'{cluster_stats[2]["death_rate"]:.0f}%',
     f'{cluster_stats[2]["os_med"]:.1f} mo',
     f'{cluster_stats[2]["pfi_med"]:.1f} mo',
     'Poor'],
    ['Cluster 4', 'Monosomy 3\n(M3)', 'BAP1 +\niso8q(MYC)', 'Class 2',
     f'{cluster_stats[3]["n"]}',
     f'{cluster_stats[3]["death_rate"]:.0f}%',
     f'{cluster_stats[3]["os_med"]:.1f} mo',
     f'{cluster_stats[3]["pfi_med"]:.1f} mo',
     'Worst'],
]

row_colors_map = ['#BBDEFB', '#C8E6C9', '#FFE0B2', '#FFCDD2']
prog_colors = {
    'Best':      '#1565C0',
    'Very Good': '#2E7D32',
    'Poor':      '#E65100',
    'Worst':     '#B71C1C',
}

table = ax.table(
    cellText=rows,
    colLabels=col_labels,
    cellLoc='center',
    loc='center',
    bbox=[0, 0.02, 1, 0.92],
)
table.auto_set_font_size(False)
table.set_fontsize(9.5)

# Header style
for col_idx in range(len(col_labels)):
    cell = table[0, col_idx]
    cell.set_facecolor('#37474F')
    cell.set_text_props(color='white', fontweight='bold')
    cell.set_height(0.22)

# Row styles
for row_idx, (row_data, row_bg) in enumerate(zip(rows, row_colors_map), start=1):
    for col_idx in range(len(col_labels)):
        cell = table[row_idx, col_idx]
        cell.set_height(0.19)
        # Cluster name column
        if col_idx == 0:
            cell.set_facecolor(row_bg)
            cell.set_text_props(fontweight='bold')
        # Prognosis column (last)
        elif col_idx == len(col_labels) - 1:
            prog = row_data[-1]
            cell.set_facecolor(prog_colors[prog])
            cell.set_text_props(color='white', fontweight='bold')
        # Chr3 status
        elif col_idx == 1:
            color = '#BBDEFB' if 'Disomy' in row_data[1] else '#FFCDD2'
            cell.set_facecolor(color)
        else:
            cell.set_facecolor('#FAFAFA')

# Add footnote
fig.text(0.5, 0.01,
         'SCNA = Somatic Copy Number Alteration  |  iso8q = Isochromosome 8q (MYC amplification)  |  '
         'OS = Overall Survival  |  PFI = Progression-Free Interval  |  Source: Robertson et al. Cancer Cell 2017',
         ha='center', fontsize=7.5, color='#666666', style='italic')

plt.savefig(f'{FIGS}\\fig0_cluster_summary_table.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: fig0_cluster_summary_table.png')

print('All figures saved to 07.Figures/')
