import os
import gc
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')

# 设置路径
workDir = "/mnt/public/zhengyeyang/CRC_project"
figurePath = os.path.join(workDir, "CRC_Figures")
save_path = os.path.join(workDir,  "integration_results")
subcluster_path = os.path.join(figurePath, "subclustering")
marker_table_path = os.path.join(workDir, "Marker_Genes")
os.makedirs(subcluster_path, exist_ok=True)
os.makedirs(marker_table_path, exist_ok=True)

# 加载已注释主要细胞类型的数据
print("加载数据...")
final_integrated = sc.read_h5ad(f"{save_path}/major_anno_all.h5ad")
print(f"数据加载完成，总细胞数: {final_integrated.shape[0]}")

# 设置期刊风格绘图参数
plt.rcParams.update({
    'font.size': 8,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.linewidth': 0.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'legend.frameon': False,
    'legend.fontsize': 7,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'figure.dpi': 300
})

# 科研级气泡图函数 - 优化版本
def create_research_bubble_plot(bubble_data, cell_type):
    """创建符合科研出版质量的气泡图 - 优化版本"""
    print("  🎨 绘制优化科研级气泡图...")
    
    # 获取唯一的cluster和gene
    unique_clusters = sorted(bubble_data['Cluster'].unique(), key=lambda x: int(x))
    unique_genes = bubble_data['Gene'].unique()
    
    # 设置图形大小 - 动态调整宽度以适应基因数量
    n_genes = len(unique_genes)
    base_width = 8  # 增加基础宽度以适应更多基因
    width = max(base_width, n_genes * 0.8)  # 每个基因0.8英寸
    height = max(5, len(unique_clusters) * 0.6)  # 增加高度
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(width, height), dpi=300)
    
    # 创建白色到紫色渐变颜色映射
    colors = ["#FFFFFF", "#F2F0F7", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D"]
    cmap = LinearSegmentedColormap.from_list("white_purple", colors)
    
    # 创建映射
    gene_to_idx = {gene: i for i, gene in enumerate(unique_genes)}
    cluster_to_idx = {cluster: i for i, cluster in enumerate(unique_clusters)}
    
    # 计算最大表达值用于颜色标准化
    if bubble_data['Mean_Expression'].max() > 0:
        max_exp = bubble_data['Mean_Expression'].max()
    else:
        max_exp = 1  # 防止除零错误
    
    # 绘制气泡
    for idx, row in bubble_data.iterrows():
        x = gene_to_idx[row['Gene']]  # 基因在x轴
        y = cluster_to_idx[row['Cluster']]  # 亚群在y轴
        
        # 增大气泡大小差异
        size = row['Pct_Expressed'] * 1.5  # 线性缩放，使大小差异更明显
        
        # 如果表达量为0，使用白色
        if row['Mean_Expression'] == 0:
            color = "#FFFFFF"
        else:
            color_val = row['Mean_Expression'] / max_exp
            color = cmap(color_val)
        
        # 设置边缘颜色
        edgecolor = 'black' if row['Pct_Expressed'] > 0 else 'lightgray'
        linewidth = 0.5 if row['Pct_Expressed'] > 0 else 0.3
        
        ax.scatter(
            x, y, 
            s=size, 
            c=[color],
            edgecolor=edgecolor,
            linewidth=linewidth,
            alpha=0.9,
            zorder=2
        )
    
    # 添加标签 - 优化显示（使用Arial细字体）
    plt.xticks(
        range(len(unique_genes)), 
        unique_genes, 
        rotation=45,  # 45度旋转避免重叠
        ha='right',   # 右对齐
        fontsize=10,  # 增大字体
        rotation_mode='anchor',  # 更好的锚点旋转
        fontname='Arial',
        fontweight='light'  # 细字体
    )
    
    plt.yticks(
        range(len(unique_clusters)), 
        [f'C{i}' for i in unique_clusters], 
        fontsize=11,  # 增大字体
        fontname='Arial',
        fontweight='light'  # 细字体
    )
    
    # 设置坐标轴标签（使用Arial字体）
    plt.xlabel('Marker Genes', fontsize=12, fontweight='bold', labelpad=12, fontname='Arial')
    plt.ylabel('Cell Subclusters', fontsize=12, fontweight='bold', labelpad=12, fontname='Arial')
    
    # 添加标题（使用Arial字体）
    plt.title(f'{cell_type} Subclusters - Marker Gene Expression', 
             fontsize=14, pad=15, fontweight='bold', fontname='Arial')
    
    # 调整边界 - 增加左右边距
    plt.xlim(-0.8, len(unique_genes)-0.2)
    plt.ylim(-0.5, len(unique_clusters)-0.5)
    
    # 增加基因名之间的间距
    plt.subplots_adjust(bottom=0.25)  # 为基因名留出更多空间
    
    # ===== 图例1: 气泡大小 (%表达) =====
    # 放置在右侧上方
    sizes = [20, 40, 60, 80]  # 表达百分比
    size_labels = [f'{s}%' for s in sizes]
    size_scales = [s * 1.5 for s in sizes]  # 与实际计算匹配
    
    # 创建图例（使用Arial细字体）
    size_legend = ax.legend(
        handles=[
            plt.scatter([], [], s=s, c='gray', alpha=0.8, edgecolor='black', label=label)
            for s, label in zip(size_scales, size_labels)
        ],
        title='% Expressed',
        title_fontsize=10,  # 增大字体
        fontsize=9,        # 增大字体
        loc='upper left',
        bbox_to_anchor=(1.02, 0.95),  # 放在右上角
        frameon=False,
        handletextpad=2.0,  # 增加图例项之间的水平间距
        labelspacing=1.8,    # 增加图例项之间的垂直间距
        prop={'family': 'Arial', 'weight': 'light'}  # 使用Arial细字体
    )
    ax.add_artist(size_legend)
    
    # ===== 图例2: 颜色条 (平均表达) =====
    # 放置在右侧下方
    if max_exp > 0:
        norm = mpl.colors.Normalize(vmin=0, vmax=max_exp)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        
        # 创建颜色条轴 - 放在右侧下方，与大小图例留出空间
        cax = fig.add_axes([0.92, 0.25, 0.02, 0.3])  # [x, y, width, height]
        cbar = plt.colorbar(sm, cax=cax)
        cbar.set_label('Mean Expression', fontsize=10, labelpad=5, fontname='Arial', fontweight='light')
        cbar.ax.tick_params(labelsize=9)
        # 设置刻度标签使用Arial细字体
        for label in cbar.ax.get_yticklabels():
            label.set_family('Arial')
            label.set_weight('light')
    
    # 调整布局
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # 为右侧图例留出空间
    
    # 保存图像
    save_file = os.path.join(subcluster_path, f'research_markers_{cell_type}.pdf')
    plt.savefig(save_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 优化气泡图保存至: {save_file}")
    return save_file

# 保存标记基因表格
def save_marker_genes_table(adata_sub, markers_dict, cell_type):
    """保存标记基因信息到CSV文件"""
    print("  💾 保存标记基因表格...")
    
    # 创建数据框存储所有标记基因
    all_markers = []
    
    # 收集每个亚群的标记基因信息
    for cluster, genes in markers_dict.items():
        # 尝试获取该亚群的完整标记基因数据
        try:
            markers_df = sc.get.rank_genes_groups_df(adata_sub, group=cluster)
            # 只保留当前选择的基因
            markers_df = markers_df[markers_df['names'].isin(genes)]
            # 添加亚群信息
            markers_df['Cluster'] = cluster
            # 添加到总表
            all_markers.append(markers_df)
        except:
            # 如果无法获取完整数据，创建简化版
            cluster_markers = pd.DataFrame({
                'names': genes,
                'Cluster': cluster,
                'logfoldchanges': [np.nan] * len(genes),
                'pvals_adj': [np.nan] * len(genes)
            })
            all_markers.append(cluster_markers)
    
    # 合并所有亚群的数据
    if all_markers:
        combined_df = pd.concat(all_markers)
        
        # 重新排序列
        combined_df = combined_df[['Cluster', 'names', 'logfoldchanges', 'pvals_adj']]
        combined_df.columns = ['Cluster', 'Gene', 'log2FC', 'Adjusted_p_value']
        
        # 计算表达百分比和平均表达量
        combined_df['Pct_Expressed'] = 0.0
        combined_df['Mean_Expression'] = 0.0
        
        for idx, row in combined_df.iterrows():
            gene = row['Gene']
            cluster = row['Cluster']
            
            # 获取目标亚群细胞
            cluster_mask = adata_sub.obs['subcluster'] == cluster
            adata_cluster = adata_sub[cluster_mask]
            
            # 如果基因不存在，保留为0
            if gene not in adata_cluster.var_names:
                combined_df.at[idx, 'Pct_Expressed'] = 0
                combined_df.at[idx, 'Mean_Expression'] = 0
                continue
                
            adata_gene = adata_cluster[:, gene]
            
            # 计算表达比例
            if isinstance(adata_gene.X, np.ndarray):
                expressed = (adata_gene.X > 0).sum()
            else:
                # 处理稀疏矩阵
                expressed = (adata_gene.X.toarray() > 0).sum()
            
            if len(adata_cluster) > 0:
                pct = expressed / len(adata_cluster) * 100
            else:
                pct = 0
            
            # 计算平均表达量 (非零细胞)
            if isinstance(adata_gene.X, np.ndarray):
                non_zero_exp = adata_gene.X[adata_gene.X > 0]
            else:
                dense_data = adata_gene.X.toarray().flatten()
                non_zero_exp = dense_data[dense_data > 0]
            
            if len(non_zero_exp) > 0:
                mean_exp = non_zero_exp.mean()
            else:
                mean_exp = 0
            
            combined_df.at[idx, 'Pct_Expressed'] = pct
            combined_df.at[idx, 'Mean_Expression'] = mean_exp
        
        # 保存到CSV
        save_file = os.path.join(marker_table_path, f'marker_genes_{cell_type}.csv')
        combined_df.to_csv(save_file, index=False)
        print(f"  ✅ 标记基因表格保存至: {save_file}")
        return combined_df
    else:
        print("  ⚠️ 没有标记基因数据可保存")
        return None

# 计算基因在特定亚群中的表达量
def calculate_gene_expression(adata_sub, gene, cluster):
    """计算基因在特定亚群中的表达百分比和平均表达量"""
    cluster_mask = adata_sub.obs['subcluster'] == cluster
    adata_cluster = adata_sub[cluster_mask]
    
    if gene not in adata_cluster.var_names:
        return 0.0, 0.0
    
    adata_gene = adata_cluster[:, gene]
    
    # 计算表达比例
    if isinstance(adata_gene.X, np.ndarray):
        expressed = (adata_gene.X > 0).sum()
    else:
        expressed = (adata_gene.X.toarray() > 0).sum()
    
    n_cells = len(adata_cluster)
    pct = (expressed / n_cells * 100) if n_cells > 0 else 0.0
    
    # 计算平均表达量 (非零细胞)
    if isinstance(adata_gene.X, np.ndarray):
        non_zero_exp = adata_gene.X[adata_gene.X > 0]
    else:
        dense_data = adata_gene.X.toarray().flatten()
        non_zero_exp = dense_data[dense_data > 0]
    
    mean_exp = non_zero_exp.mean() if len(non_zero_exp) > 0 else 0.0
    
    return pct, mean_exp

# 子聚类和可视化函数 - 使用8个标志基因
def subcluster_celltype(adata, cell_type, n_markers=8):  # 使用8个标志基因
    """对特定细胞类型进行子聚类并可视化标志基因 - 优化版本"""
    print(f"\n🔍 开始处理 {cell_type} 细胞...")
    
    # 扩展的已知生物学标志基因字典（每种细胞类型12个基因）
    known_markers = {
        'T': ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'FOXP3', 'CCR7', 'GZMB', 'PDCD1', 'CTLA4', 'IL7R', 'SELL'],
        'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'CDH1', 'KRT19', 'MUC1', 'CLDN4', 'CD24', 'CD44', 'MUC2', 'LGR5', 'OLFM4'],
        'B': ['CD19', 'CD79A', 'MS4A1', 'CD27', 'CD38', 'CD1D', 'CD22', 'CD40', 'CD80', 'CD86', 'CD69', 'CD83'],
        'Plasma': ['CD38', 'SDC1', 'JCHAIN', 'MZB1', 'PRDM1', 'XBP1', 'IGHA1', 'IGHG1', 'IGHG2', 'IGKC', 'IGLC2', 'TNFRSF17'],
        'Myeloid': ['CD14', 'CD68', 'ITGAM', 'CD163', 'CD1C', 'FCER1A', 'CD86', 'HLA-DRA', 'CD80', 'CD209', 'CLEC10A', 'CD206'],
        'Stromal': ['COL1A1', 'COL1A2', 'ACTA2', 'FAP', 'PDGFRA', 'PDGFRB', 'DCN', 'LUM', 'THY1', 'VIM', 'FN1', 'TAGLN'],
        'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'CLDN5', 'FLT1', 'KDR', 'CD34', 'ICAM1', 'VCAM1', 'SELE', 'CD93', 'ESAM'],
        'Mast': ['KIT', 'TPSAB1', 'CPA3', 'FCER1A', 'TPSB2', 'CMA1', 'HDC', 'MS4A2', 'ENPP3', 'CD63', 'CD69', 'CD9']
    }
    
    # 提取目标细胞类型
    cell_mask = adata.obs['Major_type'] == cell_type
    adata_sub = adata[cell_mask].copy()
    
    # 检查细胞数量
    n_cells = adata_sub.shape[0]
    min_cells = 100  # 最小细胞数要求
    if n_cells < min_cells:
        print(f"  ⚠️ 细胞数量不足 ({n_cells} < {min_cells}), 跳过 {cell_type}")
        return None, None
    
    print(f"  📊 细胞数量: {n_cells}")
    
    # 降采样（如果细胞数量过多）
    target_cells = 10000
    if n_cells > target_cells:
        print(f"  ⬇️ 降采样到 {target_cells} 个细胞...")
        # 随机降采样
        sampled_indices = np.random.choice(adata_sub.shape[0], size=target_cells, replace=False)
        adata_sub = adata_sub[sampled_indices].copy()
        n_cells = target_cells
    
    # 使用已有的Harmony嵌入
    if 'X_pca_harmony' in adata_sub.obsm:
        print("  🔄 使用已有的Harmony嵌入")
        adata_sub.obsm['X_pca'] = adata_sub.obsm['X_pca_harmony']
    else:
        print("  ⚠️ 未找到Harmony嵌入，使用原始PCA")
        if 'X_pca' not in adata_sub.obsm:
            # 创建伪PCA结果
            adata_sub.obsm['X_pca'] = np.zeros((adata_sub.shape[0], 2))
    
    # 构建邻接图
    print("  🕸️ 构建邻接图...")
    n_neighbors = min(15, n_cells-1)
    n_pcs = min(30, adata_sub.obsm['X_pca'].shape[1])
    
    try:
        sc.pp.neighbors(adata_sub, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep='X_pca')
    except Exception as e:
        print(f"  ⚠️ 邻接图构建失败: {str(e)}")
        print("  ⚠️ 使用简单聚类")
        # 创建伪聚类结果
        adata_sub.obs['subcluster'] = '0'
        return adata_sub, None
    
    # Leiden聚类
    print("  🧬 进行Leiden聚类...")
    try:
        sc.tl.leiden(adata_sub, resolution=0.8, key_added='subcluster')
    except Exception as e:
        print(f"  ⚠️ Leiden聚类失败: {str(e)}")
        print("  ⚠️ 使用简单聚类")
        adata_sub.obs['subcluster'] = '0'
    
    # 计算UMAP
    print("  🗺️ 计算UMAP...")
    try:
        sc.tl.umap(adata_sub)
    except Exception as e:
        print(f"  ⚠️ UMAP计算失败: {str(e)}")
    
    # 获取该细胞类型的已知标志基因
    celltype_markers = known_markers.get(cell_type, [])
    print(f"  🔬 使用已知生物学标志基因: {', '.join(celltype_markers)}")
    
    # 执行差异表达分析
    print("  🔬 计算亚群差异表达基因...")
    try:
        if len(adata_sub.obs['subcluster'].unique()) > 1:
            sc.tl.rank_genes_groups(
                adata_sub, 
                groupby='subcluster',
                method='wilcoxon',
                use_raw=False
            )
    except Exception as e:
        print(f"  ⚠️ 差异表达分析失败: {str(e)}")
    
    # 存储每个亚群的标志基因
    markers_dict = {}
    
    # 为每个亚群选择表达量最高的8个生物学标志基因
    for cluster in adata_sub.obs['subcluster'].unique():
        # 1. 计算每个生物学标志基因在该亚群中的表达得分 (pct * mean_exp)
        gene_scores = []
        for gene in celltype_markers:
            if gene in adata_sub.var_names:
                pct, mean_exp = calculate_gene_expression(adata_sub, gene, cluster)
                score = pct * mean_exp  # 综合表达得分
                gene_scores.append((gene, score))
        
        # 2. 按得分降序排序，选择前8个
        gene_scores.sort(key=lambda x: x[1], reverse=True)
        selected_genes = [gene for gene, score in gene_scores[:n_markers]]
        
        # 3. 如果生物学标志基因不足8个，添加差异表达基因
        if len(selected_genes) < n_markers:
            try:
                markers_df = sc.get.rank_genes_groups_df(adata_sub, group=cluster)
                if markers_df is not None and not markers_df.empty:
                    # 过滤显著基因 (p<0.05)
                    markers_df = markers_df[markers_df['pvals_adj'] < 0.05]
                    # 选择logfoldchange最大的基因
                    markers_df = markers_df.sort_values('logfoldchanges', ascending=False)
                    
                    # 添加差异表达基因（排除已选择的生物学知识基因）
                    for _, row in markers_df.iterrows():
                        gene = row['names']
                        if gene not in selected_genes and gene not in celltype_markers:
                            selected_genes.append(gene)
                            if len(selected_genes) >= n_markers:
                                break
            except Exception as e:
                print(f"    ⚠️ 获取差异表达基因失败: {str(e)}")
        
        # 4. 确保最终有8个基因
        markers_dict[cluster] = selected_genes[:n_markers]
        print(f"    ✅ 亚群 {cluster} - 最终标志基因: {', '.join(selected_genes[:n_markers])}")
    
    # 保存标记基因表格
    marker_table = save_marker_genes_table(adata_sub, markers_dict, cell_type)
    
    # 准备气泡图数据 - 确保包含所有标志基因
    print("  📊 准备气泡图数据...")
    all_genes = []
    cluster_list = []
    
    # 确保每个亚群包含所有标志基因
    for cluster, genes in markers_dict.items():
        # 确保每个亚群都有完整的标志基因列表
        cluster_genes = genes.copy()
        if len(cluster_genes) < n_markers:
            # 如果基因数量不足，用占位符填充
            cluster_genes += [''] * (n_markers - len(cluster_genes))
        
        all_genes.extend(cluster_genes)
        cluster_list.extend([cluster] * len(cluster_genes))
    
    if not all_genes:
        print("  ⚠️ 没有可用的标志基因，跳过气泡图")
        return adata_sub, marker_table
    
    # 创建数据框
    bubble_data = pd.DataFrame({
        'Gene': all_genes,
        'Cluster': cluster_list
    })
    
    # 计算表达比例和平均表达量
    bubble_data['Pct_Expressed'] = 0.0
    bubble_data['Mean_Expression'] = 0.0
    
    for idx, row in bubble_data.iterrows():
        gene = row['Gene']
        cluster = row['Cluster']
        
        # 跳过空基因
        if gene == '':
            continue
            
        # 获取目标亚群细胞
        cluster_mask = adata_sub.obs['subcluster'] == cluster
        adata_cluster = adata_sub[cluster_mask]
        
        # 如果基因不存在，保留为0
        if gene not in adata_cluster.var_names:
            bubble_data.at[idx, 'Pct_Expressed'] = 0
            bubble_data.at[idx, 'Mean_Expression'] = 0
            continue
            
        adata_gene = adata_cluster[:, gene]
        
        # 计算表达比例
        if isinstance(adata_gene.X, np.ndarray):
            expressed = (adata_gene.X > 0).sum()
        else:
            expressed = (adata_gene.X.toarray() > 0).sum()
        
        if len(adata_cluster) > 0:
            pct = expressed / len(adata_cluster) * 100
        else:
            pct = 0
        
        # 计算平均表达量 (非零细胞)
        if isinstance(adata_gene.X, np.ndarray):
            non_zero_exp = adata_gene.X[adata_gene.X > 0]
        else:
            dense_data = adata_gene.X.toarray().flatten()
            non_zero_exp = dense_data[dense_data > 0]
        
        if len(non_zero_exp) > 0:
            mean_exp = non_zero_exp.mean()
        else:
            mean_exp = 0
        
        bubble_data.at[idx, 'Pct_Expressed'] = pct
        bubble_data.at[idx, 'Mean_Expression'] = mean_exp
    
    # 创建科研级气泡图 - 使用优化版本
    create_research_bubble_plot(bubble_data, cell_type)
    
    # 创建UMAP可视化
    if 'X_umap' in adata_sub.obsm:
        plt.figure(figsize=(7, 6))  # 增大图形尺寸
        sc.pl.umap(
            adata_sub, 
            color='subcluster',
            title=f'{cell_type} Subclusters',
            legend_loc='on data',
            frameon=False,
            show=False,
            size=40 if n_cells < 5000 else 20,  # 增大点的大小
            palette='tab20'
        )
        plt.tight_layout()
        save_umap_file = os.path.join(subcluster_path, f'subcluster_umap_{cell_type}.pdf')
        plt.savefig(save_umap_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✅ UMAP图保存至: {save_umap_file}")
    
    # 保存子聚类结果
    save_data_file = os.path.join(save_path, f"subclustered_{cell_type}.h5ad")
    adata_sub.write(save_data_file)
    print(f"  ✅ 数据保存至: {save_data_file}")
    
    return adata_sub, marker_table

# 主要细胞类型列表
cell_types_to_analyze = ['T', 'Epithelial', 'B', 'Plasma', 'Myeloid', 'Stromal', 'Endothelial', 'Mast']

# 创建主报告
report_df = pd.DataFrame(columns=['Cell Type', 'Number of Subclusters', 'Markers CSV', 'Bubble Plot', 'UMAP'])

# 执行子聚类和可视化 - 每个亚群使用8个标志基因
for cell_type in cell_types_to_analyze:
    try:
        print(f"\n{'='*50}")
        print(f"开始处理: {cell_type}")
        print(f"{'='*50}")
        result, marker_table = subcluster_celltype(final_integrated, cell_type, n_markers=8)
        
        if result is None:
            print(f"  ⚠️ {cell_type} 处理未完成")
            continue
        
        # 添加到报告
        n_clusters = len(result.obs['subcluster'].unique())
        csv_path = os.path.join(marker_table_path, f'marker_genes_{cell_type}.csv')
        bubble_path = os.path.join(subcluster_path, f'research_markers_{cell_type}.pdf')
        umap_path = os.path.join(subcluster_path, f'subcluster_umap_{cell_type}.pdf') if 'X_umap' in result.obsm else "N/A"
        
        report_df = report_df.append({
            'Cell Type': cell_type,
            'Number of Subclusters': n_clusters,
            'Markers CSV': csv_path,
            'Bubble Plot': bubble_path,
            'UMAP': umap_path
        }, ignore_index=True)
        
    except Exception as e:
        print(f"⚠️ 处理 {cell_type} 时出错: {str(e)}")
        continue

# 保存总报告
report_file = os.path.join(workDir, "subcluster_analysis_report.csv")
report_df.to_csv(report_file, index=False)
print(f"\n✅ 分析报告保存至: {report_file}")

print("\n🎉 所有细胞类型的小群聚类和可视化完成！")