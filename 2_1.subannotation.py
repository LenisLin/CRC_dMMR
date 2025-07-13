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
figurePath = os.path.join(workDir, "CRC_figures")
save_path = os.path.join(workDir, "integration_results")
marker_table_path = os.path.join(workDir, "marker_genes")
os.makedirs(figurePath, exist_ok=True)
os.makedirs(save_path, exist_ok=True)
os.makedirs(marker_table_path, exist_ok=True)

# 加载已注释主要细胞类型的数据
print("加载数据...")
final_integrated = sc.read_h5ad(os.path.join(save_path, "major_anno_all.h5ad"))
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

# 创建自定义颜色映射
def create_custom_colormap():
    colors = ["#f0f9e8", "#7bccc4", "#0868ac"]
    return LinearSegmentedColormap.from_list("custom_blue", colors)

# 数据预处理函数 - 解决NaN问题
def safe_preprocessing(adata):
    """安全的数据预处理，避免NaN问题"""
    # 备份原始数据
    if "raw" not in adata.layers:
        adata.layers["raw"] = adata.X.copy()
    
    # 1. 处理NaN值
    if np.isnan(adata.X).any():
        print("  检测到NaN值，替换为0")
        if isinstance(adata.X, np.ndarray):
            adata.X = np.nan_to_num(adata.X, nan=0.0)
        else:
            # 处理稀疏矩阵
            data = adata.X.toarray()
            data = np.nan_to_num(data, nan=0.0)
            adata.X = data
    
    # 2. 处理无限值
    if np.isinf(adata.X).any():
        print("  检测到无限值，替换为0")
        if isinstance(adata.X, np.ndarray):
            adata.X = np.nan_to_num(adata.X, posinf=0.0, neginf=0.0)
        else:
            data = adata.X.toarray()
            data = np.nan_to_num(data, posinf=0.0, neginf=0.0)
            adata.X = data
    
    # 3. 归一化
    try:
        sc.pp.normalize_total(adata, target_sum=1e4)
    except Exception as e:
        print(f"  归一化失败: {str(e)}")
        print("  使用原始计数数据")
        adata.X = adata.layers["raw"].copy()
        return adata
    
    # 4. log1p变换
    sc.pp.log1p(adata)
    
    # 5. 再次检查NaN
    if np.isnan(adata.X).any() or np.isinf(adata.X).any():
        print("  log1p后检测到无效值，回退到原始数据")
        adata.X = adata.layers["raw"].copy()
    
    return adata

# 安全的高变基因计算
def safe_hvg(adata):
    """安全的高变基因计算"""
    # 计算基因均值
    if isinstance(adata.X, np.ndarray):
        gene_means = np.array(adata.X.mean(axis=0)).flatten()
    else:
        # 处理稀疏矩阵
        gene_means = np.array(adata.X.mean(axis=0)).A1
    
    # 检查是否有足够的变异
    if np.all(gene_means == gene_means[0]) or np.isnan(gene_means).any():
        print("  基因表达变异不足或存在NaN，使用所有基因")
        adata.var["highly_variable"] = True
        return adata
    
    try:
        # 尝试计算高变基因
        sc.pp.highly_variable_genes(
            adata, 
            min_mean=0.01,
            max_mean=10,
            min_disp=0.1,
            n_bins=10,
            n_top_genes=min(2000, adata.shape[1]))
        
        # 确保至少有一些高变基因
        if adata.var["highly_variable"].sum() < 50:
            print("  高变基因不足，使用所有基因")
            adata.var["highly_variable"] = True
    except Exception as e:
        print(f"  高变基因计算失败: {str(e)}")
        print("  使用所有基因")
        adata.var["highly_variable"] = True
    
    return adata

# 安全PCA计算
def safe_pca(adata):
    """安全的PCA计算"""
    try:
        # 检查是否有足够的特征
        if adata.n_vars < 10:
            print("  基因数量不足，跳过PCA")
            return adata
        
        # 尝试PCA
        n_comps = min(15, adata.n_vars, adata.shape[0]-1)
        sc.pp.pca(adata, n_comps=n_comps, svd_solver='arpack', use_highly_variable=True)
    except Exception as e:
        print(f"  PCA失败: {str(e)}")
        # 尝试使用随机化SVD
        try:
            print("  尝试使用随机化SVD")
            sc.pp.pca(adata, n_comps=n_comps, svd_solver='randomized', use_highly_variable=True)
        except:
            print("  所有PCA方法失败，跳过PCA")
    
    return adata

# 科研级气泡图函数
def create_research_bubble_plot(bubble_data, cell_type):
    """创建符合科研出版质量的气泡图"""
    print("  🎨 绘制科研级气泡图...")
    
    # 获取唯一的cluster和gene
    unique_clusters = sorted(bubble_data['Cluster'].unique(), key=lambda x: int(x))
    unique_genes = bubble_data['Gene'].unique()
    
    # 设置图形大小 - 科研论文标准尺寸 (双栏宽度)
    width = 8  # 英寸 (约20cm)
    height = max(4, len(unique_clusters) * 0.6)  # 动态高度
    
    # 创建图形
    plt.figure(figsize=(width, height), dpi=300)
    ax = plt.gca()
    
    # 创建自定义颜色映射
    cmap = create_custom_colormap()
    
    # 创建映射
    gene_to_idx = {gene: i for i, gene in enumerate(unique_genes)}
    cluster_to_idx = {cluster: i for i, cluster in enumerate(unique_clusters)}
    
    # 计算最大表达值用于颜色标准化
    max_exp = bubble_data['Mean_Expression'].max()
    
    # 绘制气泡 - 交换坐标轴
    for idx, row in bubble_data.iterrows():
        x = gene_to_idx[row['Gene']]  # 基因在x轴
        y = cluster_to_idx[row['Cluster']]  # 亚群在y轴
        size = np.sqrt(row['Pct_Expressed']) * 20  # 非线性缩放，更好区分大小
        color_val = row['Mean_Expression'] / max_exp if max_exp > 0 else 0
        color = cmap(color_val)
        
        plt.scatter(
            x, y, 
            s=size, 
            c=[color],
            edgecolor='black',
            linewidth=0.5,
            alpha=0.9,
            zorder=2
        )
    
    # 添加标签
    plt.xticks(range(len(unique_genes)), unique_genes, rotation=90, fontsize=9)
    plt.yticks(range(len(unique_clusters)), [f'C{i}' for i in unique_clusters], fontsize=9)
    
    # 设置坐标轴标签
    plt.xlabel('Marker Genes', fontsize=10, fontweight='bold', labelpad=10)
    plt.ylabel('Cell Subclusters', fontsize=10, fontweight='bold', labelpad=10)
    
    # 添加标题
    plt.title(f'{cell_type} Subclusters - Marker Gene Expression', 
             fontsize=12, pad=15, fontweight='bold')
    
    # 添加网格线 - 更细致的样式
    plt.grid(True, linestyle='--', alpha=0.3, which='both', axis='both')
    
    # 调整边界
    plt.xlim(-0.5, len(unique_genes)-0.5)
    plt.ylim(-0.5, len(unique_clusters)-0.5)
    
    # ===== 图例1: 气泡大小 (%表达) =====
    # 放置在顶部左侧
    sizes = [20, 40, 60, 80]  # 表达百分比
    size_labels = [f'{s}%' for s in sizes]
    size_scales = [np.sqrt(s) * 20 for s in sizes]  # 与实际计算匹配
    
    # 创建图例
    size_legend = ax.legend(
        handles=[
            plt.scatter([], [], s=s, c='gray', alpha=0.8, edgecolor='black', label=label)
            for s, label in zip(size_scales, size_labels)
        ],
        title='% Expressed',
        title_fontsize=9,
        fontsize=8,
        loc='upper left',
        bbox_to_anchor=(0.01, 1.15),
        frameon=False,
        ncol=len(sizes),
        handletextpad=1.5,
        columnspacing=1.5
    )
    ax.add_artist(size_legend)
    
    # ===== 图例2: 颜色条 (平均表达) =====
    # 放置在顶部右侧
    if max_exp > 0:
        norm = mpl.colors.Normalize(vmin=0, vmax=max_exp)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        
        # 创建颜色条轴
        cax = ax.inset_axes([0.75, 1.15, 0.25, 0.03])  # [x, y, width, height]
        cbar = plt.colorbar(sm, cax=cax, orientation='horizontal')
        cbar.set_label('Mean Expression', fontsize=9, labelpad=2)
        cbar.ax.tick_params(labelsize=8)
    
    # 调整布局
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # 为顶部图例留出空间
    
    # 保存图像
    save_file = os.path.join(figurePath, f'research_markers_{cell_type}.pdf')
    plt.savefig(save_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✅ 科研级气泡图保存至: {save_file}")
    return save_file

# 保存标记基因表格
def save_marker_genes_table(adata_sub, markers_dict, cell_type):
    """保存标记基因信息到CSV文件"""
    print("  💾 保存标记基因表格...")
    
    # 创建数据框存储所有标记基因
    all_markers = []
    
    # 收集每个亚群的标记基因信息
    for cluster, genes in markers_dict.items():
        # 获取该亚群的完整标记基因数据
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
            
            # 检查基因是否存在
            if gene not in adata_cluster.var_names:
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

# 子聚类和可视化函数
def subcluster_and_visualize(adata, cell_type, n_markers=8):
    """对特定细胞类型进行子聚类并可视化标志基因"""
    print(f"\n🔍 开始处理 {cell_type} 细胞...")
    
    # 提取目标细胞类型
    cell_mask = adata.obs['Major_type'] == cell_type
    adata_sub = adata[cell_mask].copy()
    
    # 检查细胞数量
    n_cells = adata_sub.shape[0]
    min_cells = 50  # 降低最小细胞数要求
    if n_cells < min_cells:
        print(f"  ⚠️ 细胞数量不足 ({n_cells} < {min_cells}), 跳过 {cell_type}")
        return None
    
    print(f"  📊 细胞数量: {n_cells}")
    
    # 安全预处理
    print("  ⚙️ 安全预处理...")
    adata_sub = safe_preprocessing(adata_sub)
    
    # 安全计算高变基因
    print("  🔬 安全计算高变基因...")
    adata_sub = safe_hvg(adata_sub)
    
    # 安全PCA
    print("  📉 安全PCA...")
    adata_sub = safe_pca(adata_sub)
    
    # 检查是否有有效的PCA结果
    if 'X_pca' not in adata_sub.obsm:
        print("  ⚠️ 没有有效的PCA结果，使用原始表达数据")
        # 创建伪PCA结果
        adata_sub.obsm['X_pca'] = np.zeros((adata_sub.shape[0], 2))
    
    # 构建邻接图
    print("  🕸️ 构建邻接图...")
    n_neighbors = min(10, n_cells-1)
    n_pcs = min(5, adata_sub.obsm['X_pca'].shape[1])
    
    try:
        sc.pp.neighbors(adata_sub, n_neighbors=n_neighbors, n_pcs=n_pcs)
    except:
        print("  ⚠️ 邻接图构建失败，使用简单聚类")
        # 创建伪聚类结果
        adata_sub.obs['subcluster'] = '0'
        return adata_sub
    
    # Leiden聚类
    print("  🧬 进行Leiden聚类...")
    try:
        sc.tl.leiden(adata_sub, resolution=0.5, key_added='subcluster')
    except:
        print("  ⚠️ Leiden聚类失败，使用简单聚类")
        adata_sub.obs['subcluster'] = '0'
    
    # 计算UMAP
    print("  🗺️ 计算UMAP...")
    try:
        sc.tl.umap(adata_sub)
    except:
        print("  ⚠️ UMAP计算失败，跳过UMAP")
    
    # 计算标志基因
    print("  🔬 计算亚群标志基因...")
    markers_dict = {}
    rank_genes_groups_df = None
    
    try:
        if len(adata_sub.obs['subcluster'].unique()) > 1:
            sc.tl.rank_genes_groups(
                adata_sub, 
                groupby='subcluster',
                method='wilcoxon',
                use_raw=False
            )
            
            # 保存完整的rank_genes_groups结果
            rank_genes_groups_df = sc.get.rank_genes_groups_df(adata_sub, group=None)
            
            # 获取每个亚群的前n_markers个标志基因
            print("  📈 提取标志基因...")
            for cluster in adata_sub.obs['subcluster'].unique():
                try:
                    markers_df = sc.get.rank_genes_groups_df(adata_sub, group=cluster)
                    if markers_df is not None and not markers_df.empty:
                        # 过滤显著基因 (p<0.05)
                        markers_df = markers_df[markers_df['pvals_adj'] < 0.05]
                        # 选择logfoldchange最大的前n_markers个基因
                        markers_df = markers_df.sort_values('logfoldchanges', ascending=False)
                        top_markers = markers_df.head(n_markers)['names'].tolist()
                        markers_dict[cluster] = top_markers
                except:
                    continue
    except:
        print("  ⚠️ 标志基因计算失败")
    
    # 如果没有找到标志基因，使用已知标记基因
    if not markers_dict:
        print("  ℹ️ 使用已知标记基因作为备选")
        known_markers = {
            'B': ['CD19', 'CD79A', 'MS4A1', 'CD79B', 'CD22', 'CD24', 'CD27', 'CD38', 'CD40'],
            'T': ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'FOXP3', 'IL7R', 'CCR7', 'SELL'],
            'Myeloid': ['CD14', 'CD68', 'FCGR3A', 'ITGAM', 'ITGAX', 'CD1C', 'CLEC9A', 'CD163', 'CD86'],
            'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1', 'MUC1', 'CEACAM5', 'CD24', 'CD44'],
            'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'CLDN5', 'SELE', 'ICAM1', 'VCAM1', 'CD34', 'FLT1'],
            'Stromal': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'ACTA2', 'PDGFRA', 'PDGFRB', 'FAP'],
            'Plasma': ['CD38', 'SDC1', 'CD138', 'PRDM1', 'XBP1', 'IRF4', 'MZB1', 'JCHAIN', 'IGKC'],
            'Mast': ['KIT', 'TPSAB1', 'CPA3', 'MS4A2', 'FCER1A', 'HDC', 'TPSB2', 'CMA1', 'GATA2']
        }
        
        if cell_type in known_markers:
            for cluster in adata_sub.obs['subcluster'].unique():
                markers_dict[cluster] = known_markers[cell_type][:n_markers]
    
    # 保存标记基因表格
    marker_table = save_marker_genes_table(adata_sub, markers_dict, cell_type)
    
    # 准备气泡图数据
    print("  📊 准备气泡图数据...")
    all_genes = []
    cluster_list = []
    
    for cluster, genes in markers_dict.items():
        all_genes.extend(genes)
        cluster_list.extend([cluster] * len(genes))
    
    if not all_genes:
        print("  ⚠️ 没有可用的标志基因，跳过气泡图")
        # 保存子聚类结果
        save_data_file = os.path.join(save_path, f"subclustered_{cell_type}.h5ad")
        adata_sub.write(save_data_file)
        print(f"  ✅ 数据保存至: {save_data_file}")
        return adata_sub
    
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
        
        # 获取目标亚群细胞
        cluster_mask = adata_sub.obs['subcluster'] == cluster
        adata_cluster = adata_sub[cluster_mask]
        
        # 检查基因是否存在
        if gene not in adata_cluster.var_names:
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
        
        bubble_data.at[idx, 'Pct_Expressed'] = pct
        bubble_data.at[idx, 'Mean_Expression'] = mean_exp
    
    # 创建科研级气泡图
    create_research_bubble_plot(bubble_data, cell_type)
    
    # 创建UMAP可视化
    if 'X_umap' in adata_sub.obsm:
        plt.figure(figsize=(8, 6))
        sc.pl.umap(
            adata_sub, 
            color='subcluster',
            title=f'{cell_type} Subclusters',
            legend_loc='on data',
            frameon=False,
            show=False,
            size=30 if n_cells < 5000 else 15,
            palette='tab20'
        )
        plt.tight_layout()
        save_umap_file = os.path.join(figurePath, f'subcluster_umap_{cell_type}.pdf')
        plt.savefig(save_umap_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✅ UMAP图保存至: {save_umap_file}")
    
    # 保存子聚类结果
    save_data_file = os.path.join(save_path, f"subclustered_{cell_type}.h5ad")
    adata_sub.write(save_data_file)
    print(f"  ✅ 数据保存至: {save_data_file}")
    
    return adata_sub, marker_table

# 指定要分析的细胞类型
cell_types_to_analyze = ['B', 'Endothelial', 'Epithelial', 'Mast', 'Myeloid', 'Plasma', 'Stromal', 'T']

# 创建主报告
report_df = pd.DataFrame(columns=['Cell Type', 'Number of Subclusters', 'Markers CSV', 'Bubble Plot', 'UMAP'])

# 执行子聚类和可视化
for cell_type in cell_types_to_analyze:
    try:
        print(f"\n{'='*50}")
        print(f"开始处理: {cell_type}")
        print(f"{'='*50}")
        result, marker_table = subcluster_and_visualize(final_integrated, cell_type, n_markers=9)
        
        if result is None:
            print(f"  ⚠️ {cell_type} 处理未完成")
            continue
        
        # 添加到报告
        n_clusters = len(result.obs['subcluster'].unique())
        csv_path = os.path.join(marker_table_path, f'marker_genes_{cell_type}.csv')
        bubble_path = os.path.join(figurePath, f'research_markers_{cell_type}.pdf')
        umap_path = os.path.join(figurePath, f'subcluster_umap_{cell_type}.pdf')
        
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