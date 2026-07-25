"""
光敏剂分子 MaxMin 采样多维度科研评估工具 (v11.0)
===================================================
更新说明：
    - 标题修正：第一列 [采样算法]:[投影方法]；第二列 [投影方法]:[来源分布]。
    - 增加 MDS 支持：使 MaxMin 采样在视觉上呈现更真实的均匀感。
    - 四图联动：空间覆盖、来源溯源、NND分布、采样强度比。
"""

import pandas as pd
import numpy as np
import argparse
import os
import warnings
from tqdm import tqdm

# 环境配置与警告屏蔽
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

from sklearn.manifold import TSNE, MDS
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances

import matplotlib
matplotlib.use('Agg')  # 适配服务器环境
import matplotlib.pyplot as plt
import seaborn as sns

def load_pictos_labels(pictos_file):
    """
    从pictos.txt文件加载pic名称到标签的映射
    
    Args:
        pictos_file: pictos.txt文件路径
    
    Returns:
        字典，键为pic名称（如"pic3"），值为标签（如"s9v1"）
    """
    labels = {}
    if not os.path.exists(pictos_file):
        return labels
    
    try:
        with open(pictos_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if ':' in line:
                    parts = line.split(':', 1)
                    pic_name = parts[0].strip()
                    label = parts[1].strip()
                    labels[pic_name] = label
                else:
                    parts = line.split()
                    if len(parts) >= 2:
                        pic_name = parts[0]
                        label = parts[1]
                        labels[pic_name] = label
    except Exception as e:
        print(f"警告: 无法读取pictos文件: {e}")
    
    return labels

def extract_pic_name(filename):
    """
    从文件名中提取pic名称（如 pic3inpdftv22.csv -> pic3）
    """
    basename = os.path.basename(filename)
    # 移除扩展名
    name_without_ext = os.path.splitext(basename)[0]
    # 尝试匹配 pic + 数字
    import re
    match = re.match(r'(pic\d+)', name_without_ext, re.IGNORECASE)
    if match:
        return match.group(1).lower()
    return None

def parse_args():
    parser = argparse.ArgumentParser(description="MaxMin 采样科研评估工具 v11.0")
    parser.add_argument('input', nargs='+', help='输入CSV文件')
    parser.add_argument('--n_pick', '-n', type=int, default=10, help='采样数量')
    parser.add_argument('--col', '-c', default='SMILES', help='SMILES列名')
    parser.add_argument('--output', '-o', default='maxmin_samples.csv', help='输出文件名')
    parser.add_argument('--plot-dir', default='results/figures/chemical_space', help='图和投影坐标输出目录')
    parser.add_argument('--pictos', default='config/pictos.txt', help='批次到展示标签的映射文件')
    parser.add_argument('--method', choices=['umap', 'tsne', 'mds'], default='tsne', help='可视化投影算法')
    parser.add_argument('--pre_sample', type=int, default=50000, help='MaxMin预选池大小')
    parser.add_argument('--radius', type=int, default=2, help='Morgan 指纹半径')
    parser.add_argument('--fp-size', type=int, default=2048, help='Morgan 指纹位数')
    parser.add_argument('--seed', type=int, default=42, help='采样和降维随机种子')
    parser.add_argument('--umap-neighbors', type=int, default=30, help='UMAP 邻居数')
    parser.add_argument('--umap-min-dist', type=float, default=0.1, help='UMAP 最小距离')
    parser.add_argument('--tsne-perplexity', type=float, default=30.0, help='t-SNE perplexity')
    parser.add_argument('--first-as-bg', action='store_true', help='将第一个输入文件的点显示为灰色背景')
    args = parser.parse_args()
    if args.n_pick < 1:
        parser.error("--n_pick must be at least 1")
    if args.pre_sample < 1:
        parser.error("--pre_sample must be at least 1")
    if args.radius < 1 or args.fp_size < 1:
        parser.error("--radius and --fp-size must be positive")
    if args.umap_neighbors < 2 or args.tsne_perplexity <= 0:
        parser.error("UMAP neighbors must be >= 2 and t-SNE perplexity must be > 0")
    return args

def load_data(files, smiles_col, pictos_labels=None):
    all_data = []
    first_file_name = None
    file_to_label = {}  # 文件名到标签的映射
    print("\n--- Step 1: 加载数据与溯源 ---")
    for idx, f in enumerate(files):
        if not os.path.exists(f): continue
        df_temp = pd.read_csv(f)
        target = next((c for c in [smiles_col, 'SMILES', 'smiles', 'Smiles'] if c in df_temp.columns), None)
        if target:
            source_name = os.path.basename(f)
            if idx == 0:
                first_file_name = source_name
            
            # 尝试从文件名提取pic名称并映射到标签
            pic_name = extract_pic_name(f)
            display_label = None
            if pictos_labels and pic_name:
                display_label = pictos_labels.get(pic_name)
                if display_label:
                    # 统一格式为大写S
                    if display_label.lower().startswith('s'):
                        display_label = f"S{display_label[1:].upper()}"
                    file_to_label[source_name] = display_label
                    print(f"   {source_name} -> {pic_name} -> {display_label}")
            
            subset = df_temp[[target]].copy()
            subset.columns = ['smiles']
            subset['source_file'] = source_name
            all_data.append(subset)
    full_df = pd.concat(all_data).drop_duplicates(subset=['smiles']).reset_index(drop=True)
    full_counts = full_df['source_file'].value_counts().sort_index()
    print(f"✅ 数据合并完成。总分子数: {len(full_df)}")
    if first_file_name:
        print(f"   第一个文件: {first_file_name}")
    return full_df, full_counts, first_file_name, file_to_label

def calculate_nnd(sampled_fps):
    nnd_list = []
    # 如果没有采样值或只有一个采样值，返回空列表
    if len(sampled_fps) == 0:
        return nnd_list
    if len(sampled_fps) == 1:
        # 只有一个点时，无法计算NND，返回一个默认值（例如1.0，表示完全不同的距离）
        return [1.0]
    
    for i in range(len(sampled_fps)):
        others = sampled_fps[:i] + sampled_fps[i+1:]
        # 如果others为空（理论上不应该发生，但为了安全）
        if len(others) == 0:
            nnd_list.append(1.0)
            continue
        sims = DataStructs.BulkTanimotoSimilarity(sampled_fps[i], others)
        # 如果sims为空，跳过或使用默认值
        if len(sims) == 0:
            nnd_list.append(1.0)
        else:
            nnd_list.append(1 - max(sims))
    return nnd_list

def main():
    args = parse_args()
    os.makedirs(args.plot_dir, exist_ok=True)
    
    # 加载pictos映射
    pictos_labels = load_pictos_labels(args.pictos)
    if pictos_labels:
        print(f"已加载 {len(pictos_labels)} 个pic标签映射")
    
    df, full_counts, first_file_name, file_to_label = load_data(args.input, args.col, pictos_labels)
    smiles_list = df['smiles'].tolist()

    # 1. 指纹生成
    gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=args.radius, fpSize=args.fp_size
    )
    fps, valid_idx = [], []
    for i, sm in enumerate(tqdm(smiles_list, desc="Fingerprinting")):
        mol = Chem.MolFromSmiles(str(sm))
        if mol:
            fps.append(gen.GetFingerprint(mol))
            valid_idx.append(i)

    # 2. 采样逻辑
    if not fps:
        raise SystemExit("没有可解析的 SMILES，无法采样。")
    n_pick = min(args.n_pick, len(fps))
    np.random.seed(args.seed)
    
    # MaxMin 采样
    pre_indices = np.random.choice(range(len(fps)), min(len(fps), args.pre_sample), replace=False)
    subset_fps = [fps[i] for i in pre_indices]
    picker = MaxMinPicker()
    maxmin_sub_idx = picker.LazyBitVectorPick(
        subset_fps, len(subset_fps), n_pick, seed=args.seed
    )
    maxmin_indices = [pre_indices[i] for i in maxmin_sub_idx]

    # 随机采样
    random_indices = np.random.choice(range(len(fps)), n_pick, replace=False)

    # 保存 MaxMin 选中的分子，供后续计算或人工复核
    selected_df = df.iloc[[valid_idx[i] for i in maxmin_indices]].copy()
    selected_df = selected_df.rename(columns={"smiles": args.col})
    output_parent = os.path.dirname(args.output)
    if output_parent:
        os.makedirs(output_parent, exist_ok=True)
    selected_df.to_csv(args.output, index=False, encoding="utf-8")
    print(f"✅ MaxMin 样本已保存至: {args.output} (共 {len(selected_df)} 个分子)")

    # 3. 统计比例
    full_probs = full_counts / full_counts.sum()
    source_list = full_counts.index.tolist()

    def get_stats(indices):
        counts = df.iloc[[valid_idx[i] for i in indices]]['source_file'].value_counts().reindex(source_list, fill_value=0)
        return counts / counts.sum(), (counts / full_counts) * 100

    maxmin_probs, maxmin_rate = get_stats(maxmin_indices)
    random_probs, random_rate = get_stats(random_indices)

    # 4. 投影降维
    vis_method = args.method.upper()
    print(f"\n--- Step 4: 执行 {vis_method} 空间投影 ---")
    
    bg_limit = 5000 if args.method == 'mds' else 12000
    bg_size = min(bg_limit, len(fps))
    bg_idx = np.random.choice(range(len(fps)), bg_size, replace=False)
    
    # 组合点
    combined_fps_list = [fps[i] for i in list(bg_idx) + list(maxmin_indices) + list(random_indices)]
    
    if args.method == 'umap':
        try:
            import umap
        except ModuleNotFoundError as exc:
            raise SystemExit(
                "UMAP 方法需要 umap-learn；请先按 environment.yml 安装环境，"
                "或改用 --method tsne/mds。"
            ) from exc
        reducer = umap.UMAP(
            metric='jaccard',
            n_neighbors=min(args.umap_neighbors, max(2, len(combined_fps_list) - 1)),
            min_dist=args.umap_min_dist,
            random_state=args.seed,
        )
        embedding = reducer.fit_transform(np.array(combined_fps_list))
    elif args.method == 'tsne':
        fingerprint_array = np.array(combined_fps_list)
        pca_components = min(50, fingerprint_array.shape[0], fingerprint_array.shape[1])
        pca_data = PCA(n_components=pca_components).fit_transform(fingerprint_array)
        perplexity = min(args.tsne_perplexity, max(1.0, len(pca_data) - 1.0))
        reducer = TSNE(
            n_components=2,
            perplexity=perplexity,
            random_state=args.seed,
            init='pca',
        )
        embedding = reducer.fit_transform(pca_data)
    elif args.method == 'mds':
        bool_fps = np.array([np.array(fp) for fp in combined_fps_list], dtype=bool)
        dist_matrix = pairwise_distances(bool_fps, metric='jaccard', n_jobs=-1)
        reducer = MDS(
            n_components=2,
            dissimilarity='precomputed',
            random_state=args.seed,
            n_jobs=-1,
        )
        embedding = reducer.fit_transform(dist_matrix)

    # 坐标切片
    bg_embed = embedding[:bg_size]
    maxmin_embed = embedding[bg_size : bg_size + n_pick]
    random_embed = embedding[bg_size + n_pick :]
    bg_sources = df.iloc[[valid_idx[i] for i in bg_idx]]['source_file'].values

    # 保存所有点的坐标和标签到CSV文件
    coord_data = []
    
    # 收集所有采样点的索引（用于标记）
    sampled_indices_set = set()
    if len(maxmin_indices) > 0:
        sampled_indices_set.update(maxmin_indices)
    if len(random_indices) > 0:
        sampled_indices_set.update(random_indices)
    
    # 保存所有背景点
    for i in range(len(bg_embed)):
        original_idx = bg_idx[i]
        source_file = bg_sources[i]
        pic_name = extract_pic_name(source_file)
        label = file_to_label.get(source_file, None)
        is_sampled = original_idx in sampled_indices_set
        sampling_method = ''
        if is_sampled:
            if original_idx in maxmin_indices:
                sampling_method = 'MaxMin'
            elif original_idx in random_indices:
                sampling_method = 'Random'
        
        coord_data.append({
            'x': bg_embed[i, 0],
            'y': bg_embed[i, 1],
            'source_file': source_file,
            'pic_name': pic_name if pic_name else '',
            'label': label if label else '',
            'is_sampled': is_sampled,
            'sampling_method': sampling_method
        })
    
    # 保存MaxMin采样点（如果不在背景点中）
    if len(maxmin_embed) > 0 and len(maxmin_indices) > 0:
        maxmin_sources = df.iloc[[valid_idx[i] for i in maxmin_indices]]['source_file'].values
        for i in range(len(maxmin_embed)):
            original_idx = maxmin_indices[i]
            # 如果这个点已经在背景点中保存过了，跳过
            if original_idx in bg_idx:
                continue
            source_file = maxmin_sources[i]
            pic_name = extract_pic_name(source_file)
            label = file_to_label.get(source_file, None)
            coord_data.append({
                'x': maxmin_embed[i, 0],
                'y': maxmin_embed[i, 1],
                'source_file': source_file,
                'pic_name': pic_name if pic_name else '',
                'label': label if label else '',
                'is_sampled': True,
                'sampling_method': 'MaxMin'
            })
    
    # 保存Random采样点（如果不在背景点中）
    if len(random_embed) > 0 and len(random_indices) > 0:
        random_sources = df.iloc[[valid_idx[i] for i in random_indices]]['source_file'].values
        for i in range(len(random_embed)):
            original_idx = random_indices[i]
            # 如果这个点已经在背景点中保存过了，跳过
            if original_idx in bg_idx:
                continue
            source_file = random_sources[i]
            pic_name = extract_pic_name(source_file)
            label = file_to_label.get(source_file, None)
            coord_data.append({
                'x': random_embed[i, 0],
                'y': random_embed[i, 1],
                'source_file': source_file,
                'pic_name': pic_name if pic_name else '',
                'label': label if label else '',
                'is_sampled': True,
                'sampling_method': 'Random'
            })
    
    # 保存为CSV
    coord_df = pd.DataFrame(coord_data)
    output_stem = os.path.splitext(os.path.basename(args.output))[0]
    coord_output_path = os.path.join(
        args.plot_dir, f"{output_stem}_{args.method}_coordinates.csv"
    )
    coord_df.to_csv(coord_output_path, index=False, encoding='utf-8')
    print(f"✅ 坐标数据已保存至: {coord_output_path} (共 {len(coord_data)} 个点)")

    # 如果启用 first-as-bg，识别第一个文件的背景点
    bg_first_file_mask = None
    if args.first_as_bg and first_file_name:
        bg_first_file_mask = bg_sources == first_file_name
        print(f"   将第一个文件 '{first_file_name}' 的点显示为灰色背景")

    # 5. 绘图矩阵
    fig, axes = plt.subplots(2, 4, figsize=(32, 14), dpi=300)
    
    def plot_row(row_idx, embed_data, sampling_name, color, fps_indices, probs, rate, file_to_label):
        # A. 第一列：[算法]:[方法]
        if args.first_as_bg and bg_first_file_mask is not None:
            # 分别绘制第一个文件（灰色）和其他文件（浅灰色）
            first_bg_mask = bg_first_file_mask
            other_bg_mask = ~first_bg_mask
            if np.any(first_bg_mask):
                axes[row_idx, 0].scatter(bg_embed[first_bg_mask, 0], bg_embed[first_bg_mask, 1], 
                                        c='#95a5a6', s=5, alpha=0.6, label='First file (bg)')
            if np.any(other_bg_mask):
                axes[row_idx, 0].scatter(bg_embed[other_bg_mask, 0], bg_embed[other_bg_mask, 1], 
                                        c='#ecf0f1', s=5, alpha=0.4)
        else:
            axes[row_idx, 0].scatter(bg_embed[:, 0], bg_embed[:, 1], c='#ecf0f1', s=5, alpha=0.4)
        axes[row_idx, 0].scatter(embed_data[:, 0], embed_data[:, 1], c=color, s=25, edgecolors='white', linewidths=0.3)
        axes[row_idx, 0].set_title(f"{sampling_name}: {vis_method}", fontsize=16, fontweight='bold')
        
        # B. 第二列：[方法]: Source Distribution
        # 这一列背景是投影方法生成的来源云图
        if args.first_as_bg and bg_first_file_mask is not None:
            # 为第一个文件使用灰色，其他文件使用彩色
            first_bg_mask = bg_first_file_mask
            other_bg_mask = ~first_bg_mask
            if np.any(first_bg_mask):
                axes[row_idx, 1].scatter(bg_embed[first_bg_mask, 0], bg_embed[first_bg_mask, 1], 
                                        c='#ecf0f1', s=10, alpha=0.8, label='First file (bg)')
            if np.any(other_bg_mask):
                # 只显示有映射的文件在图例中
                other_sources = bg_sources[other_bg_mask]
                other_embed = bg_embed[other_bg_mask]
                # 过滤出有标签映射的文件
                mapped_mask = np.array([file_to_label.get(src, None) is not None for src in other_sources])
                # 绘制所有点（包括没有映射的），但不显示在图例中
                unmapped_mask = ~mapped_mask
                if np.any(unmapped_mask):
                    # 没有映射的点用浅灰色绘制，不显示在图例中
                    axes[row_idx, 1].scatter(other_embed[unmapped_mask, 0], other_embed[unmapped_mask, 1], 
                                            c='#ecf0f1', s=10, alpha=0.3)
                if np.any(mapped_mask):
                    mapped_sources = other_sources[mapped_mask]
                    mapped_embed = other_embed[mapped_mask]
                    # 将文件名替换为标签
                    mapped_labels = [file_to_label.get(src, src) for src in mapped_sources]
                    sns.scatterplot(x=mapped_embed[:, 0], y=mapped_embed[:, 1], 
                                  hue=mapped_labels, s=10, alpha=0.5, palette='husl', 
                                  ax=axes[row_idx, 1], legend=(row_idx==0))
                    
                    # 计算每个标签的中心位置并添加文本标注
                    unique_labels = list(set(mapped_labels))
                    for label in unique_labels:
                        label_mask = np.array([l == label for l in mapped_labels])
                        label_points = mapped_embed[label_mask]
                        if len(label_points) > 0:
                            # 计算中心位置（使用中位数，更稳健）
                            center_x = np.median(label_points[:, 0])
                            center_y = np.median(label_points[:, 1])
                            # 添加文本标注（无背景框和边框）
                            axes[row_idx, 1].text(center_x, center_y, label, 
                                                fontsize=12, fontweight='bold',
                                                ha='center', va='center')
        else:
            # 只显示有映射的文件在图例中
            mapped_mask = np.array([file_to_label.get(src, None) is not None for src in bg_sources])
            # 绘制所有点（包括没有映射的），但不显示在图例中
            unmapped_mask = ~mapped_mask
            if np.any(unmapped_mask):
                # 没有映射的点用浅灰色绘制，不显示在图例中
                axes[row_idx, 1].scatter(bg_embed[unmapped_mask, 0], bg_embed[unmapped_mask, 1], 
                                        c='#ecf0f1', s=10, alpha=0.3)
            if np.any(mapped_mask):
                mapped_sources = bg_sources[mapped_mask]
                mapped_embed = bg_embed[mapped_mask]
                # 将文件名替换为标签
                mapped_labels = [file_to_label.get(src, src) for src in mapped_sources]
                sns.scatterplot(x=mapped_embed[:, 0], y=mapped_embed[:, 1], 
                              hue=mapped_labels, s=10, alpha=0.5, palette='husl', 
                              ax=axes[row_idx, 1], legend=(row_idx==0))
                
                # 计算每个标签的中心位置并添加文本标注
                unique_labels = list(set(mapped_labels))
                for label in unique_labels:
                    label_mask = np.array([l == label for l in mapped_labels])
                    label_points = mapped_embed[label_mask]
                    if len(label_points) > 0:
                        # 计算中心位置（使用中位数，更稳健）
                        center_x = np.median(label_points[:, 0])
                        center_y = np.median(label_points[:, 1])
                        # 添加文本标注（无背景框和边框）
                        axes[row_idx, 1].text(center_x, center_y, label, 
                                            fontsize=12, fontweight='bold',
                                            ha='center', va='center')
        axes[row_idx, 1].scatter(embed_data[:, 0], embed_data[:, 1], facecolors='none', edgecolors='black', s=30, linewidths=0.8)
        axes[row_idx, 1].set_title(f"{vis_method}: Source Distribution", fontsize=16, fontweight='bold')
        if row_idx == 0: axes[row_idx, 1].legend(title="Files", bbox_to_anchor=(1.05, 1), loc='upper left')

        # C. 第三列：NND
        nnd_vals = calculate_nnd([fps[i] for i in fps_indices])
        if len(nnd_vals) > 0:
            sns.histplot(nnd_vals, bins=30, kde=True, ax=axes[row_idx, 2], color=color)
            mean_nnd = np.mean(nnd_vals)
            axes[row_idx, 2].set_title(f"NND Distribution (Mean={mean_nnd:.2f})", fontsize=16, fontweight='bold')
        else:
            axes[row_idx, 2].text(0.5, 0.5, 'No NND data', ha='center', va='center', transform=axes[row_idx, 2].transAxes)
            axes[row_idx, 2].set_title(f"NND Distribution (No data)", fontsize=16, fontweight='bold')

        # D. 第四列：比例与强度 (双轴)
        ax_l = axes[row_idx, 3]
        ax_r = ax_l.twinx()
        x = np.arange(len(source_list))
        width = 0.3
        ax_l.bar(x - width/2, full_probs.values, width, label='Full Library Prop.', color='#bdc3c7', alpha=0.7)
        ax_l.bar(x + width/2, probs.values, width, label='Sampled Prop.', color=color, alpha=0.8)
        ax_r.plot(x, rate.values, color='#2c3e50', marker='o', linestyle='--', linewidth=2, label='Internal Sampling Rate (%)')
        ax_l.set_xticks(x)
        ax_l.set_xticklabels(source_list, rotation=45, ha='right', fontsize=9)
        ax_l.set_title(f"Ratio & Sampling Intensity", fontsize=16, fontweight='bold')
        
        lines, labels = ax_l.get_legend_handles_labels()
        lines2, labels2 = ax_r.get_legend_handles_labels()
        ax_l.legend(lines + lines2, labels + labels2, loc='upper right', fontsize='x-small')

    plot_row(0, maxmin_embed, "MaxMin", "#e74c3c", maxmin_indices, maxmin_probs, maxmin_rate, file_to_label)
    plot_row(1, random_embed, "Random", "#2ecc71", random_indices, random_probs, random_rate, file_to_label)

    plt.tight_layout()
    img_path = os.path.join(
        args.plot_dir, f"{output_stem}_{args.method}_v11_final.png"
    )
    plt.savefig(img_path, bbox_inches='tight')
    print(f"\n📊 终极科研评估图已保存至: {img_path}")

if __name__ == "__main__":
    main()
