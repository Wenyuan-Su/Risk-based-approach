# 1. 安装并加载必要的包
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

# 2. 读取数据
# 假设文件在当前目录下
data <- read.csv("bioprofile_matrix.csv", row.names = 1, check.names = FALSE)

# 清洗列名 (去掉 .0)
colnames(data) <- as.integer(as.numeric(colnames(data)))

# 3. 设置颜色 (根据您的要求)
# -1: 蓝色 (Inactive)
#  0: 白色/浅灰 (Inconclusive) - 用白色可以让有颜色的点更突出
#  1: 红色 (Active)
my_colors <- c("#3B4CC0", "#DDDCDC", "#B40426") 
my_breaks <- c(-1.5, -0.5, 0.5, 1.5) # 严格划分区间

# 4. 【核心步骤】预先计算聚类，为了做“稀疏标签”
# 使用曼哈顿距离 (Manhattan) 和 Ward 算法，适合这种离散数据
dist_rows <- dist(data, method = "euclidean")
dist_cols <- dist(t(data), method = "euclidean")
hc_rows <- hclust(dist_rows, method = "complete")
hc_cols <- hclust(dist_cols, method = "complete")

# 5. 【核心步骤】生成稀疏标签
# 逻辑：根据聚类后的顺序(hc_rows$order)，每隔 N 个显示原始名称，其余设为空

# --- 设置 Y 轴 (Compound ID) 的稀疏度 ---
# round(nrow(data)/40) 意思是总共大概显示 40 个标签，您可以手动改比如 20
step_row <- max(1, round(nrow(data) / 40)) 
labels_row_sparse <- rep("", nrow(data)) # 先全部初始化为空

# 只在特定位置填回标签
# hc_rows$order[i] 表示热图中第 i 行对应原始数据的哪一行
for (i in 1:nrow(data)) {
  if (i %% step_row == 0) {
    original_index <- hc_rows$order[i]
    labels_row_sparse[original_index] <- rownames(data)[original_index]
  }
}

# --- 设置 X 轴 (Assay ID) 的稀疏度 ---
step_col <- max(1, round(ncol(data) / 30)) # 总共大概显示 30 个标签
labels_col_sparse <- rep("", ncol(data))

for (i in 1:ncol(data)) {
  if (i %% step_col == 0) {
    original_index <- hc_cols$order[i]
    labels_col_sparse[original_index] <- colnames(data)[original_index]
  }
}

# 6. 绘图
png(filename = "pubchem_final_heatmap.png", width = 12, height = 10, units = "in", res = 300)

pheatmap(data,
         color = my_colors,
         breaks = my_breaks,
         
         # 使用我们预先计算好的聚类对象
         cluster_rows = hc_rows,
         cluster_cols = hc_cols,
         
         # 使用我们要生成的稀疏标签
         labels_row = labels_row_sparse,
         labels_col = labels_col_sparse,
         
         # 视觉微调
         border_color = NA,      # 去掉网格线
         fontsize_row = 8,       # 标签字体大小
         fontsize_col = 8,
         angle_col = 45,         # 列名倾斜
         
         # 图例设置
         legend_breaks = c(-1, 0, 1),
         legend_labels = c("Inactive (-1)", "Inconclusive (0)", "Active (1)"),
         
         main = "Bioactivity Heatmap (Sparse Labels & Clustered)"
)

dev.off()
print("图片已保存为: pubchem_final_heatmap.png")