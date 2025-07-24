# 2020/04/21

# Single-cell RNA-seq analysis - QC
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(MAST)
library(RColorBrewer)

options (future.globals.maxSize = 4000 * 1024^5)
dir_prefix <- "E:/fkh-and-hongyuan-single-cell/hd5_files-hy/"
dir_suffix <- "/raw_feature_bc_matrix.h5"


for (file in c("IMR904", "NC3_1", "SA5", "SA7")){
  seurat_data <- Read10X_h5(paste0(dir_prefix, file, dir_suffix))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 200, 
                                   project = file)
  assign(file, seurat_obj)
}

rm(seurat_data, seurat_obj)
# 创建一个合并的Seurat对象
merged_seurat <- merge(x = IMR904, 
                       y = c(NC3_1, SA5, SA7), 
                       add.cell.id = c("IMR904", "NC3_1", "SA5", "SA7"))

rm(IMR904, NC3_1, SA5, SA7)

# 将每个细胞每个UMI的基因数目添加到元数据中
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# 计算线粒体比率
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# 创建元数据数据框
metadata <- merged_seurat@meta.data

# 为元数据添加细胞ID
metadata$cells <- rownames(metadata)

# 重命名列
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# 创建样本列
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^IMR"))] <- "CTRL1"
metadata$sample[which(str_detect(metadata$cells, "^NC"))] <- "CTRL2"
metadata$sample[which(str_detect(metadata$cells, "^SA5"))] <- "MDD1"
metadata$sample[which(str_detect(metadata$cells, "^SA7"))] <- "MDD2"

# 将元数据添加回Seurat对象中
merged_seurat@meta.data <- metadata

merged_seurat <- readRDS("new_merged_seurat.rds")
metadata <- filtered_seurat@meta.data

# 可视化每个样本的细胞计数
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# 可视化每个细胞的UMI/转录本数目
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# 通过频数图可视化每个细胞检测出的基因数分布
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# 通过箱线图可视化每个细胞检测到的基因的分布
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# 可视化检测到的基因数和UMI数之间的关系，并且观察是否存在大量低数目的基因数/UMI数的细胞
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# 可视化每个细胞检测到的线粒体基因表达分布
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# 通过可视化每一个UMI检测到的基因数来可视化基因表达的整体复杂性
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

VlnPlot(merged_seurat, 
        features = c("nGene", "nUMI", "mitoRatio"), 
        ncol = 3)

## Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & (nUMI < 35000) &
                            (nGene >= 500) & (nGene < 7500) &
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.10))

table(filtered_seurat$sample)

# 可视化筛选后每个样本的细胞计数
filtered_seurat@meta.data %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# 提取计数
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# 根据在每个细胞的计数是否大于0为每个基因输出一个逻辑向量
nonzero <- counts > 0

# 将所有TRUE值相加，如果每个基因的TRUE值超过10个，则返回TRUE。
keep_genes <- Matrix::rowSums(nonzero) >= 10

# 仅保留那些在10个以上细胞中表达的基因
filtered_counts <- counts[keep_genes, ]

# 重新赋值给经过过滤的Seurat对象
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# 标准化计数
seurat_phase <- NormalizeData(filtered_seurat)

# 导入细胞周期标记物
load("Cell_cycle&annotation/cycle.rda")

# 给细胞的细胞周期评分
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# 观察细胞周期评分和分配给每个细胞的周期 
View(head(seurat_phase@meta.data))
write.csv(head(seurat_phase@meta.data),"data/cell_cycle_phase.csv")

# 识别最大变异基因
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# 缩放计数
seurat_phase <- ScaleData(seurat_phase)

# 运行PCA
seurat_phase <- RunPCA(seurat_phase)

# 绘制以细胞周期着色的PCA图
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

options(future.globals.maxSize = 4000 * 1024^5)
options()

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

split_seurat
split_seurat <- split_seurat[c("CTRL1", "CTRL2", "MDD1", "MDD2")]

library(future)
plan()
# set multiprocess
plan("multiprocess", workers = 4)

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio", "S.Score", "G2M.Score"))
}

rm(counts, filtered_counts, filtered_seurat, nonzero, merged_seurat,seurat_phase)
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
plan("multiprocess", workers = 1)
saveRDS(integ_anchors, "new_integ_anchors.rds")
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
saveRDS(seurat_integrated, "data/seurat_after_SCT_integration.rds")
seurat_integrated <- readRDS("data/seurat_after_SCT_integration.rds")
# 运行PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# 绘制PCA
PCAPlot(seurat_integrated,group.by="Phase",split.by="Phase")  

# 运行UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# 绘制UMAP                             
DimPlot(seurat_integrated, reduction = "umap", group.by = "sample",pt.size = 0.1) 

# 运行TSNE
seurat_integrated <- RunTSNE(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# 绘制TSNE                             
DimPlot(seurat_integrated, reduction = "tsne", group.by = "sample",pt.size = 0.1) 

# 确定k-近邻图
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# 确定聚类的不同分辨率                               
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

# 将分辨率设为0.8
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# 绘制UMAP图
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

# 从seurat对象中提取身份和样本信息，以确定每个类群中每个样本的细胞数
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample) %>%
  tidyr::spread(ident, n)

# 查看表格
View(n_cells)

Cluster_sample_composition <- as.data.frame(table(seurat_integrated$sample, 
                                                  Idents(seurat_integrated)))
colnames(Cluster_sample_composition) <- c("Sample", "Cluster", "Num")
colpalette <- c("#D99F6C", "#F2EA79", "#F2AE2E", "#8C4E03")
ggplot(Cluster_sample_composition, aes(fill=Sample, y=Num, x=Cluster)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values=colpalette) +
  theme(
    # Hide panel borders and remove grid lines
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

# 选择RNA计数插槽作为默认分析
DefaultAssay(seurat_integrated) <- "RNA"

# 导入基因注释
annotations <- read.csv("data/annotation.csv")

# 创建可以获得每个类群保留标记物的函数
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

levels(Idents(seurat_integrated))

# 使用创建的函数迭代每个类群
conserved_markers_total <- map_dfr(c(0:25), get_conserved)

write.csv(conserved_markers_total, "data/conserved_markers_total.csv")
seurat_integrated <- readRDS("data/seurat_after_clustering.rds")

seurat_integrated <- subset(seurat_integrated_backup, subset=integrated_snn_res.0.6!=0)
Idents(seurat_integrated) <- "integrated_snn_res.0.6"
seurat_integrated <- RenameIdents(seurat_integrated, 
                                  "1" = "aRG",
                                  "2" = "EN",
                                  "3" = "IN",
                                  "4" = "EN",
                                  "5" = "CH",
                                  "6" = "EN",
                                  "7" = "EN",
                                  "8" = "CH",
                                  "9" = "VP",
                                  "10" = "dRG",
                                  "11" = "IPC",
                                  "12" = "oRG",
                                  "13" = "EN",
                                  "14" = "aRG",
                                  "15" = "ET",
                                  "16" = "dRG",
                                  "17" = "oRG",
                                  "18" = "ET",
                                  "19" = "CH",
                                  "20" = "Others",
                                  "21" = "VP")
DimPlot(seurat_integrated, reduction = "umap", label = T)
seurat_integrated$new_celltypes <- Idents(seurat_integrated)

list <- c("#F57F73", "#418C41", "#A52929", "#2A276A", "#B74B9B", "#6DCDDD", "#ED1E24", "#8150A0", "#4D95CF", "#CCCCCC")
DimPlot(seurat_integrated, label = F, cols = list, label.box = F) + NoAxes()

Idents(seurat_integrated) <- "new_celltypes"
Progenitor_cells <- subset(seurat_integrated, subset=new_celltypes%in%c("dRG", "oRG", "aRG", "CH", "VP", "IPC"))
Cluster_sample_composition <- as.data.frame(table(Progenitor_cells$group, Progenitor_cells$Phase))
colpalette <- c("#FC7F01","#FCBD6F","#349E2A")
colnames(Cluster_sample_composition) <- c("Group", "Cluster", "Num")
Cluster_sample_composition %>% 
  ggplot(aes(fill=Cluster, y=Num, x=Group)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values=colpalette) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
table(Progenitor_cells$group, Progenitor_cells$Phase)

Progenitor_DEGs <- FindMarkers(Progenitor_cells, 
                               ident.1 = "MDD", 
                               ident.2 = "CTRL", 
                               group.by = "group", 
                               logfc.threshold=0.25,
                               min.pct=0.1, 
                               min.diff.pct=0.1)


Cluster_sample_composition <- as.data.frame(table(seurat_integrated$group, seurat_integrated$new_celltypes))
colnames(Cluster_sample_composition) <- c("Group", "Cluster", "Num")
colpalette <- c("#FC7F01","#FCBD6F","#349E2A", "#B1DC89")
colpalette <- c("#F57F73", "#418C41", "#A52929", "#2A276A", "#B74B9B", "#6DCDDD", "#ED1E24", "#8150A0", "#4D95CF", "#CCCCCC")
Cluster_sample_composition %>% 
  ggplot(aes(fill=Cluster, y=Num, x=Group)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values=colpalette) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

Control_compositon <- Cluster_sample_composition %>%
  filter(Group=="CTRL")
MDD_compositon <- Cluster_sample_composition %>%
  filter(Group=="MDD")

seurat_integrated_cds <- as.cell_data_set(seurat_integrated)
seurat_integrated_cds <- cluster_cells(cds = seurat_integrated_cds, reduction_method = "UMAP")
seurat_integrated_cds <- learn_graph(seurat_integrated_cds, use_partition = TRUE)
plot_cells(seurat_integrated_cds,
           color_cells_by = "new_celltypes",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

metadata <- seurat_integrated@meta.data
CH_cells <- metadata %>%
  filter(new_celltypes=="CH")
seurat_integrated_cds <- order_cells(seurat_integrated_cds, reduction_method = "UMAP", root_cells = CH_cells$cells)

plot_cells(seurat_integrated_cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")

seurat_integrated <- AddMetaData(
  object = seurat_integrated,
  metadata = as.numeric(seurat_integrated_cds@principal_graph_aux@listData$UMAP$pseudotime),
  col.name = "monocle3_pseudotime"
)

# 将元数据添加回Seurat对象中
metadata <- seurat_integrated@meta.data
library(ggridges)
# 可视化每个样本的细胞计数
metadata %>% 
  ggplot(aes(y=sample, x=monocle3_pseudotime, fill=sample)) + 
  theme_classic()

# 通过频数图可视化每个细胞检测出的基因数分布
metadata %>% 
  filter(monocle3_pseudotime > 0.0001) %>% 
  ggplot(aes(color=sample, x=monocle3_pseudotime, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()



metadata %>% 
  filter(monocle3_pseudotime > 0.0001) %>% 
  ggplot(aes(y=new_celltypes, x=monocle3_pseudotime, fill=new_celltypes)) + 
  geom_density_ridges(scale=1) + 
  theme_classic() +
  scale_x_log10()

library(DAseq)

python2use <- "/home/hongyuan/anaconda3/bin/python3.8"
GPU <- 1

head(X.label.info)

labels_res <- X.label.info[X.label.info$condition == "R", "label"]
labels_nonres <- X.label.info[X.label.info$condition == "NR", "label"]
da_cells <- getDAcells(
  X = X.melanoma,
  cell.labels = X.label.melanoma,
  labels.1 = labels_res,
  labels.2 = labels_nonres,
  k.vector = seq(50, 500, 50),
  plot.embedding = X.2d.melanoma
)

head(X.melanoma)
head(X.label.melanoma)
head(X.2d.melanoma)

# 获取细胞编号
cell_ids <- row.names(seurat_integrated@meta.data)
# 获取细胞来源
cell_source <- factor(seurat_integrated$group)
# 根据细胞来源换分细胞，然后抽样
cell_ids_list <- lapply(split(cell_ids, cell_source), function(x){ sample(x, 6000) })
# 提取细胞编号
cell_id <- unlist(cell_ids_list)
#seurat_integrated <- subset(seurat_integrated_backup,subset=cells%in%cell_id)
table(seurat_integrated$sample)
X.myPC <- as.data.frame(Embeddings(seurat_integrated, reduction = "pca")[, 1:10])
X.myUMAP <- as.data.frame(Embeddings(seurat_integrated, reduction = "umap")[, 1:2])
X.mylabel <- seurat_integrated_backup@meta.data[rownames(X.myPC), "sample"]
labels_mdd <- c("MDD1", "MDD2")
labels_ctrl <- c("CTRL1", "CTRL2")
my_da_cells <- getDAcells(
  X = X.myPC,
  cell.labels = X.mylabel,
  labels.1 = labels_mdd,
  labels.2 = labels_ctrl,
  k.vector = seq(50, 500, 50),
  plot.embedding = X.myUMAP
)
seurat_integrated@meta.data[rownames(X.myPC), "sample"]
my_da_cells$pred.plot

my_da_cells <- updateDAcells(
  X = my_da_cells, pred.thres = c(-0.5,0.5),
  plot.embedding = X.myUMAP
)

my_da_regions <- getDAregion(
  X = X.myPC,
  da.cells = my_da_cells,
  cell.labels = X.mylabel,
  labels.1 = labels_mdd,
  labels.2 = labels_ctrl,
  resolution = 0.02,
  plot.embedding = X.myUMAP,
)

my_da_regions$da.region.plot

my_da_regions$da.region.label
head(seurat_integrated@assays$RNA)
seurat_integrated@assays$RNA
X.data.melanoma <- read.table(
  "./GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
  sep = "\t", header = F, row.names = 1, stringsAsFactors = F, skip = 2
)
HY.data <- as.matrix(seurat_integrated@assays$RNA@data)
X.data.melanoma <- as.matrix(X.data.melanoma[,-16292])

STG_markers <- STGmarkerFinder(
  X = HY.data,
  da.regions = my_da_regions,
  lambda = 1.5, n.runs = 5, return.model = T,
  python.use = python2use, GPU = GPU
)


da4_genes <- head(STG_markers$da.markers[["4"]], n=5)
da4_genes <- da4_genes$gene   



da_4 <- STG_markers$da.markers[["4"]]
"BMP7" %in% da_4$gene
seurat_integrated$da_regions <- my_da_regions$da.region.label
seurat_integrated$da_scores <- my_da_cells$da.pred

DimPlot(seurat_integrated_backup, reduction = "umap", group.by = "da_regions",
        cols = c("#BEBEBD", "#FF776E", "#A2A600", "#00BF7B", "#FF6BF7"),
        pt.size = 1)


metadata <- seurat_integrated@meta.data

FeaturePlot(seurat_integrated_backup, 
            reduction = "umap", 
            features = "da_scores",
            cols = colorRampPalette(colors = c("blue","white","red"))(100))
library(tidyverse)
da_score_results_regions <- metadata %>%
  group_by(da_regions) %>%
  summarise(mean_score=mean(da_scores))

colnames(da_score_results_regions) <- c("Clusters", "mean_score")

da_score_results_celltype <- metadata %>%
  group_by(new_celltypes) %>%
  summarise(mean_score=mean(da_scores))  

colnames(da_score_results_celltype) <- c("Clusters", "mean_score")

total_da_scores <- rbind(da_score_results_regions, da_score_results_celltype)
total_da_scores <- as.data.frame(total_da_scores[c(-1),])
total_da_scores$Clusters <- NULL
rownames(total_da_scores) <- total_da_scores$Clusters
library(pheatmap)
da_score_results_celltype <- as.data.frame(da_score_results_celltype)
rownames(da_score_results_celltype) <- da_score_results_celltype$new_celltypes
da_score_results_celltype$new_celltypes <- NULL
pheatmap(da_score_results_celltype,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         border_color = "white", display_numbers = T)

a <- DotPlot(seurat_integrated, 
             features = c("ID1", "ID2", "ID3", "BMP7", "HES5"),
             group.by = "da_regions")
data <- a$data

library(RColorBrewer)
pal <- rev(brewer.pal(11, 'RdYlBu'))[2:10]

newcol <- colorRampPalette(pal)
ncols <- 50
pal <- newcol(ncols)
data %>% 
  ggplot(aes(y=id, x=features.plot)) +
  geom_point(aes(fill=avg.exp.scaled, size = pct.exp), 
             pch=21,colour="black", stroke=1) +  
  
  scale_fill_gradientn(colours = pal) +   
  scale_size(range = c(2, 10), name="pct.exp")  +
  theme( panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         axis.line = element_line(colour = "black"),
         axis.text.x = element_text(angle = 90)
  )
STG_markers$da.markers[["4"]]
DimPlot(seurat_integrated_backup, reduction = "umap", group.by = "group")

library(RColorBrewer)
pal <- rev(brewer.pal(11, 'RdYlBu'))[2:10]

newcol <- colorRampPalette(pal)
ncols <- 50
pal <- newcol(ncols)

control_cells <- row.names(subset(pData(seurat_integrated_cds),
                                  group == "CTRL"))
MDD_cells <- row.names(subset(pData(seurat_integrated_cds),
                              group == "MDD"))
control_cds <- seurat_integrated_cds[,control_cells]
MDD_cds <- seurat_integrated_cds[,MDD_cells]

genes_selected <- c("MKI67","CDC20","CDK1", "TOP2A", "TACC3","AURKA")
pt.matrix <- exprs(control_cds)[match(genes_selected,rownames(rowData(control_cds))),order(pseudotime(control_cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes_selected;
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
#K means with 6 groups
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  #km = 6,
  row_title_rot                = 0,
  cluster_rows                 = F,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
print(htkm)
#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
print(htkm)

library(hdWGCNA)

seurat_obj  <- SetupForWGCNA(
  seurat_integrated_backup,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "new_WGCNA" # the name of the hdWGCNA experiment
)

str(seurat_obj@misc)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("new_celltypes", "group"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'new_celltypes' # set the Idents of the metacell seurat object
)
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("aRG", "EN", "IN", "CH", "VP", "dRG", "IPC", "oRG", "ET"),
  group.by='new_celltypes'
)
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

library("patchwork")
# assemble with patchwork

wrap_plots(plot_list, ncol=2)

seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'NPC' # name of the topoligical overlap matrix written to disk
)


seurat_integrated$group <- factor(x = seurat_integrated$group, levels = c('MDD', 'CTRL'))
