# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 16:40:29 2026

@author: Zz
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Draw import rdMolDraw2D # 引入底层高清绘图引擎
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from collections import Counter
from PIL import Image # 用于处理 TIFF 和 DPI
import io
import warnings
import os

warnings.filterwarnings("ignore")

print("=== 启动数据驱动的致癌骨架提取 (Butina 聚类 + 最小骨架提取 + 出版级加粗优化) ===")

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

all_unique_scaffolds = set(active_counts.keys()).union(set(inactive_counts.keys()))
valid_scaffolds = [s for s in all_unique_scaffolds if (active_counts.get(s, 0) + inactive_counts.get(s, 0)) >= 3]
print(f"   初步提取到 {len(valid_scaffolds)} 种有效骨架，准备开始聚类...")

# =======================================================
# 3. 利用 Butina 算法进行骨架聚类
# =======================================================
print("3. 正在运行 Butina 算法对相似骨架进行聚类合并...")

mols = [Chem.MolFromSmiles(smi) for smi in valid_scaffolds]
fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(m, 2, 1024) for m in mols]

dists = []
nfps = len(fps)
for i in range(1, nfps):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
    dists.extend([1 - x for x in sims])

cluster_cutoff = 0.5
clusters = Butina.ClusterData(dists, nfps, cluster_cutoff, isDistData=True)
print(f"   聚类完成！原本 {nfps} 种骨架被合并为了 {len(clusters)} 个超级类簇！")

# =======================================================
# 4. 汇总超级类簇的致癌富集得分，并提取最小核心骨架
# =======================================================
print("4. 正在汇总各聚类簇的致癌富集得分，并提取最小核心骨架...")

clustered_results = []
for cluster in clusters:
    scaffolds_in_cluster = [valid_scaffolds[idx] for idx in cluster]
    
    # 强制选取簇内重原子数最少的骨架作为代表 (Minimal Toxicophore)
    rep_scaffold_smi = min(
        scaffolds_in_cluster, 
        key=lambda smi: Chem.MolFromSmiles(smi).GetNumHeavyAtoms() if Chem.MolFromSmiles(smi) else float('inf')
    )
    
    total_act = sum(active_counts.get(valid_scaffolds[idx], 0) for idx in cluster)
    total_inact = sum(inactive_counts.get(valid_scaffolds[idx], 0) for idx in cluster)
    
    f_act = total_act / n_active if n_active > 0 else 0
    f_inact = total_inact / n_inactive if n_inactive > 0 else 0
    
    enrichment_score = f_act - f_inact
    
    clustered_results.append({
        'Rep_Scaffold_SMILES': rep_scaffold_smi,
        'Cluster_Size': len(cluster),
        'Total_in_Active': total_act,
        'Total_in_Inactive': total_inact,
        'Freq_in_Active': round(f_act, 4),
        'Enrichment_Score': round(enrichment_score, 4)
    })

df_clustered = pd.DataFrame(clustered_results)
df_clustered = df_clustered.sort_values(by='Enrichment_Score', ascending=False).reset_index(drop=True)

csv_filename = "consensus_Clustered_Minimal_Scaffolds_Alerts.csv"
df_clustered.to_csv(csv_filename, index=False)
print(f"✅ 统计完成！最小超级类簇报告已导出至: {csv_filename}")


# =======================================================
# 核心外挂：定义底层加粗字体绘图函数 (完美高亮修复版)
# =======================================================
def draw_bold_molecule_to_tiff(mol, match_atoms, match_bonds, tiff_filename):
    """使用 rdMolDraw2D 底层引擎：实现加粗化学键与放大原子的 300 DPI TIFF 导出"""
    # 1. 建立高清画布
    drawer = rdMolDraw2D.MolDraw2DCairo(1000, 1000)
    opts = drawer.drawOptions()
    
    # --- 【你需要的最强视觉控制面板】 ---
    opts.bondLineWidth = 4               # 加粗基础化学键 (调到4已经非常清晰)
    opts.highlightBondWidthMultiplier = 3  # 让红色的高亮线是黑色键的3倍粗！彻底包裹住化学键！
    opts.minFontSize = 45                # 强制放大所有的原子字母(如 N, O, Cl)字号！
    # -----------------------------------
    
    # 2. 避免 C++ 底层报错，强制把提取出来的 tuple 索引转为 list
    match_atoms_list = list(match_atoms) if match_atoms else []
    match_bonds_list = list(match_bonds) if match_bonds else []
    
    # 3. 开始画图
    drawer.DrawMolecule(mol, highlightAtoms=match_atoms_list, highlightBonds=match_bonds_list)
    drawer.FinishDrawing()
    
    # 4. 抓取图片数据并转为 300 DPI 无损 TIFF
    png_data = drawer.GetDrawingText()
    img = Image.open(io.BytesIO(png_data))
    img.save(tiff_filename, format="TIFF", dpi=(300, 300))


# =======================================================
# 5. 出版级可视化：分离导出 Top 10 SAs 高清 TIFF (300 DPI)
# =======================================================
print("\n5. 正在为 Top 10 最小 SAs 生成独立的 300 DPI 高清 TIFF 图片，并提取匹配分子 SMILES...")

output_dir = "SA_HighRes_Images"
os.makedirs(output_dir, exist_ok=True)

top_10_scaffolds = df_clustered.head(10)
active_mols_smiles = df[df['Final_Label'] == 1]['canon_smiles'].tolist()

# >>> 新增：用来收集最后画图用的 10 个分子 SMILES <<<
highlighted_smiles_records = [] 

for idx, row in top_10_scaffolds.iterrows():
    scaffold_smi = row['Rep_Scaffold_SMILES']
    enrichment = row['Enrichment_Score']
    patt = Chem.MolFromSmarts(scaffold_smi)
    
    if not patt:
        continue
        
    matched_mol = None
    matched_smi = None  # 用于存储当前匹配到的完整 SMILES
    match_atoms = None
    match_bonds = []
    
    # 寻找匹配分子
    for smi in active_mols_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol and mol.HasSubstructMatch(patt):
            matched_mol = mol
            matched_smi = smi  # 锁定这个成功匹配的 SMILES
            
            # 获取匹配到的原子索引
            match_atoms = mol.GetSubstructMatch(patt)
            
            # 获取匹配到的化学键索引
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() in match_atoms and bond.GetEndAtomIdx() in match_atoms:
                    match_bonds.append(bond.GetIdx())
            break
            
    if matched_mol:
        tiff_filename = os.path.join(output_dir, f"Rank_{idx+1}_Score_{enrichment:.2f}.tiff")
        
        # >>> 调用自定义的高清加粗出图函数 <<<
        draw_bold_molecule_to_tiff(matched_mol, match_atoms, match_bonds, tiff_filename)
        
        # 在控制台醒目地输出这个分子的 SMILES
        print(f"   -> [Rank {idx+1}] 已生成图片: {tiff_filename}")
        print(f"      📍 匹配分子的SMILES: {matched_smi}")
        
        # 把信息存入列表，留作后用
        highlighted_smiles_records.append({
            'Rank': idx + 1,
            'Scaffold_SMARTS': scaffold_smi,
            'Matched_Molecule_SMILES': matched_smi
        })
        
    else:
        print(f"   ⚠️ 警告: Rank {idx+1} 骨架未找到匹配的可视化分子。")

# >>> 循环结束后，将这 10 个 SMILES 统一存入一个 CSV 文件 <<<
if highlighted_smiles_records:
    df_matched = pd.DataFrame(highlighted_smiles_records)
    smiles_output_file = os.path.join(output_dir, "Top10_Highlighted_Molecules.csv")
    df_matched.to_csv(smiles_output_file, index=False)
    print(f"\n✅ Top 10 画图所用分子的 SMILES 已汇总导出至: {smiles_output_file}")
    
# =======================================================
# 6. 单分子测试仪 (输入 SMILES 自动匹配高危 SAs 并出图)
# =======================================================
print("\n=== 6. 专属测试台：目标分子高亮分析 ===")

def analyze_custom_molecule(target_smiles, df_scaffolds, output_name="Custom_Target.tiff"):
    mol = Chem.MolFromSmiles(target_smiles)
    if not mol:
        print("❌ 无效的 SMILES，无法解析。")
        return
    
    print(f"正在扫描分子: {target_smiles}")
    hit_flag = False
    
    for idx, row in df_scaffolds.iterrows():
        patt = Chem.MolFromSmarts(row['Rep_Scaffold_SMILES'])
        
        if mol.HasSubstructMatch(patt):
            match_atoms = mol.GetSubstructMatch(patt)
            match_bonds = []
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() in match_atoms and bond.GetEndAtomIdx() in match_atoms:
                    match_bonds.append(bond.GetIdx())
            
            # >>> 同样调用自定义的高清加粗出图函数 <<<
            draw_bold_molecule_to_tiff(mol, match_atoms, match_bonds, output_name)
            
            print(f"🔥 警报命中！该分子含有排名第 {idx+1} 的核心致癌骨架！")
            print(f"   致癌富集得分: {row['Enrichment_Score']:.4f}")
            print(f"✅ 高清分析图已保存为: {output_name}")
            
            hit_flag = True
            break 
            
    if not hit_flag:
        print("🟢 扫描通过：该分子未命中数据集中提取的任何致癌核心骨架。")

# --- 调用单分子测试仪演示 ---
print("\n[演示调用单分子测试仪]")
test_smiles = "C[N+](=C1C=CC(=C(c2ccc(cc2)N(C)C)c2ccc(cc2)N(C)C)C=C1)C" 
analyze_custom_molecule(test_smiles, df_clustered, output_name="My_Pollutant_Alert.tiff")