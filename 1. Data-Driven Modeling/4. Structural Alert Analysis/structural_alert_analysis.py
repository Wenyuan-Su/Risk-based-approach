# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 16:40:29 2026

@author: Zz
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from collections import Counter
import warnings

warnings.filterwarnings("ignore")

print("=== 启动数据驱动的致癌骨架提取 (引入 Butina 聚类降重) ===")

# =======================================================
# 1. 加载共识数据集
# =======================================================
df = pd.read_csv("all_compounds_9_datasets_activity.csv")
df = df.dropna(subset=['canon_smiles', 'Final_Label']).reset_index(drop=True)
df['Final_Label'] = df['Final_Label'].astype(int)

active_smiles = df[df['Final_Label'] == 1]['canon_smiles'].tolist()
inactive_smiles = df[df['Final_Label'] == 0]['canon_smiles'].tolist()

n_active = len(active_smiles)
n_inactive = len(inactive_smiles)
print(f"1. 数据加载完毕: 阳性致癌物 {n_active} 个, 阴性安全物 {n_inactive} 个。")

# =======================================================
# 2. 剥离侧链，提取 Bemis-Murcko 骨架
# =======================================================
print("2. 正在自动剥离侧链，提取分子的核心骨架...")

def get_scaffold(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        try:
            core = MurckoScaffold.GetScaffoldForMol(mol)
            core_smi = Chem.MolToSmiles(core)
            return core_smi if core_smi != '' else None
        except:
            return None
    return None

active_scaffolds = [get_scaffold(smi) for smi in active_smiles]
inactive_scaffolds = [get_scaffold(smi) for smi in inactive_smiles]

active_counts = Counter([s for s in active_scaffolds if s])
inactive_counts = Counter([s for s in inactive_scaffolds if s])

# 获取所有至少出现过 3 次的有效骨架（初步剔除极度罕见的长尾噪音）
all_unique_scaffolds = set(active_counts.keys()).union(set(inactive_counts.keys()))
valid_scaffolds = [s for s in all_unique_scaffolds if (active_counts.get(s, 0) + inactive_counts.get(s, 0)) >= 3]

print(f"   初步提取到 {len(valid_scaffolds)} 种有效骨架，准备开始聚类...")

# =======================================================
# 3. 核心升级：利用 Butina 算法进行骨架聚类
# =======================================================
print("3. 正在运行 Butina 算法对相似骨架进行聚类合并...")

mols = [Chem.MolFromSmiles(smi) for smi in valid_scaffolds]
# 计算 Morgan 指纹用于聚类 (半径=2, 1024位)
fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(m, 2, 1024) for m in mols]

# 计算距离矩阵 (1 - Tanimoto相似度)
dists = []
nfps = len(fps)
for i in range(1, nfps):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
    dists.extend([1 - x for x in sims])

# 运行 Butina 聚类
# cutoff = 0.35 意味着两张骨架指纹的相似度 >= 0.65 就会被强行合并到一起
# 你可以调大这个值 (比如 0.4) 来让合并更加激进，多环芳烃会被合得更狠
cluster_cutoff = 0.5
clusters = Butina.ClusterData(dists, nfps, cluster_cutoff, isDistData=True)

print(f"   聚类完成！原本 {nfps} 种骨架被合并为了 {len(clusters)} 个超级类簇！")

# =======================================================
# 4. 汇总超级类簇的致癌富集得分
# =======================================================
print("4. 正在汇总各聚类簇的致癌富集得分...")

clustered_results = []

for cluster in clusters:
    center_idx = cluster[0] # Butina 簇的第一个元素是中心代表骨架
    rep_scaffold_smi = valid_scaffolds[center_idx]
    
    # 将该簇内所有变体骨架的频次相加
    total_act = sum(active_counts.get(valid_scaffolds[idx], 0) for idx in cluster)
    total_inact = sum(inactive_counts.get(valid_scaffolds[idx], 0) for idx in cluster)
    
    # 重新计算该超级类簇的整体发生率
    f_act = total_act / n_active if n_active > 0 else 0
    f_inact = total_inact / n_inactive if n_inactive > 0 else 0
    
    enrichment_score = f_act - f_inact
    
    clustered_results.append({
        'Rep_Scaffold_SMILES': rep_scaffold_smi,
        'Cluster_Size': len(cluster), # 包含了多少个子骨架
        'Total_in_Active': total_act,
        'Total_in_Inactive': total_inact,
        'Freq_in_Active': round(f_act, 4),
        'Enrichment_Score': round(enrichment_score, 4)
    })

df_clustered = pd.DataFrame(clustered_results)
df_clustered = df_clustered.sort_values(by='Enrichment_Score', ascending=False).reset_index(drop=True)

# 导出聚类后的骨架 CSV
csv_filename = "consensus_Clustered_Scaffolds_Alerts.csv"
df_clustered.to_csv(csv_filename, index=False)
print(f"✅ 统计完成！超级类簇报告已导出至: {csv_filename}")
print("\n[Top 5 最危险的超级致癌骨架类簇]:")
print(df_clustered[['Rep_Scaffold_SMILES', 'Cluster_Size', 'Enrichment_Score']].head(5))

# =======================================================
# 5. 可视化：提取排名第一的代表性骨架，并在原分子上高亮
# =======================================================
print("\n5. 正在生成最危险代表性骨架的高亮可视化图片...")

if not df_clustered.empty:
    top_1_scaffold_smi = df_clustered.iloc[0]['Rep_Scaffold_SMILES']
    patt_top_1 = Chem.MolFromSmarts(top_1_scaffold_smi)
    
    # 找出包含这个代表性骨架的真实致癌分子
    target_mols_df = df[(df['Final_Label'] == 1) & 
                        (df['canon_smiles'].apply(lambda x: Chem.MolFromSmiles(x) is not None and Chem.MolFromSmiles(x).HasSubstructMatch(patt_top_1)))]
    
    sample_smiles = target_mols_df['canon_smiles'].head(6).tolist()
    
    plot_mols = []
    plot_highlights = []
    
    for smi in sample_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            match_indices = mol.GetSubstructMatch(patt_top_1)
            if match_indices:
                plot_mols.append(mol)
                plot_highlights.append(match_indices)
    
    if plot_mols:
        img = Draw.MolsToGridImage(
            plot_mols, 
            molsPerRow=3, 
            subImgSize=(350, 300), 
            highlightAtomLists=plot_highlights,
            legends=["Clustered Core Scaffold"] * len(plot_mols)
        )
        
        img_filename = "Highlight_Top1_Clustered_Scaffold.png"
        img.save(img_filename)
        print(f"✅ 核心骨架高亮图片已成功导出为: {img_filename}")
        
        # 在 Jupyter Notebook 中可以直接调用：
        # display(img)
else:
    print("未提取到满足条件的骨架。")