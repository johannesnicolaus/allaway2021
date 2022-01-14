library(tidyverse)
library(Seurat)
source("scripts/R/functions/general_functions.R")


# preprocessing of data ---------------------------------------------------

dir.create("results/plots/preprocessing")

# read data
data <- Read10X("results/data/cellranger_rna/neg_pos/outs/count/filtered_feature_bc_matrix/")
data_seurat <- CreateSeuratObject(counts = data, project = "allaway2021", min.cells = 3, min.features = 200)

# read GTF
gtf <- read_gtf("data/reference/gencode.vM25.annotation.gtf.gz")

# calculate percentage of mitochondrial genes
data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^mt-")

# visualize data
VlnPlot(data_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave("results/plots/preprocessing/vln_preprocess.svg", width = 7, height = 7)

# filter based on features and percentage of MT
data_seurat <- subset(data_seurat, subset = nFeature_RNA > 2500 & 
                        nFeature_RNA < 5000 &
                        percent.mt < 3.5 &
                        percent.mt > 1.5 &
                        nCount_RNA < 20000 &
                        nCount_RNA > 2000)


# cell cycle scoring ------------------------------------------------------
# Normalize data and identify highly variable features
data_seurat <- NormalizeData(data_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
data_seurat <- FindVariableFeatures(data_seurat, selection.method = "vst", nfeatures = 2000)
data_seurat <- ScaleData(data_seurat, features = rownames(data_seurat))

# run PCA and visualize
data_seurat <- RunPCA(data_seurat, features = VariableFeatures(object = data_seurat))
DimPlot(data_seurat, reduction = "pca")

# read cell cycle genes
cell_cycle_genes <- read_csv(url("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv"))

# add gene names
cell_cycle_markers <- left_join(cell_cycle_genes, 
                                gtf %>% select(gene_id, gene_name), 
                                by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

# Perform cell cycle scoring
data_seurat <- CellCycleScoring(data_seurat,
                                g2m.features = g2m_genes,
                                s.features = s_genes)

# Perform PCA and color by cell cycle phase
data_seurat <- RunPCA(data_seurat)

# Visualize the PCA, grouping by cell cycle phase
DimPlot(data_seurat,
        reduction = "pca",
        group.by= "Phase")       

ggsave("results/plots/preprocessing/pca_cell_cycle.svg", height = 7, width = 7)


# perform data normalization using baynorm ---------------------------------------------------------
# filter features
counts_filtered <- data_seurat@assays$RNA@counts

bay_out <- bayNorm::bayNorm(data_seurat@assays$RNA@counts, mean_version = TRUE, UMI_sffl = 10)

bay_out_counts <- bay_out$Bay_out %>% as.data.frame() %>% rownames_to_column("Geneid") %>% as_tibble()

# save normalized data as csv
dir.create("results/data/rna_preprocessing")
write_csv(bay_out_counts, "results/data/rna_preprocessing/baynorm_normalized.csv")

# create baynorm seurat object 
baynorm_seurat <- CreateSeuratObject(counts = bay_out$Bay_out, assay = "bayNorm")

# combine assays ----------------------------------------------------------
# store normalized seurat data in misc
data_seurat@misc[["seurat_data"]] <- as.matrix(x = data_seurat@misc)

# change RNA matrix to baynorm normalized data
data_seurat@assays$RNA@data <- baynorm_seurat@assays$bayNorm@data
data_seurat@assays$RNA@counts <- baynorm_seurat@assays$bayNorm@data

# data normalization ------------------------------------------------------
data_seurat <- NormalizeData(data_seurat, verbose = TRUE)
data_seurat <- FindVariableFeatures(data_seurat, selection.method = "vst", nfeatures = 2000)

# scale data based on mitoratio and cell cycle scoring
data_seurat <- ScaleData(data_seurat, vars.to.regress = c("S.Score", "G2M.Score"))

# check PCA after regression
data_seurat <- RunPCA(data_seurat)
DimPlot(data_seurat,
        reduction = "pca",
        group.by= "Phase")  
ggsave("results/plots/preprocessing/pca_cell_cycle.svg", height = 7, width = 7)

# perform clustering ------------------------------------------------------

# check variable features
LabelPoints(VariableFeaturePlot(data_seurat), points = head(VariableFeatures(data_seurat), 20), repel = T)

dir.create("results/plots/rna_clustering")
ggsave("results/plots/rna_clustering/variable_features.svg", width = 7, height = 7)

# check number of PC to use for clustering
ElbowPlot(data_seurat)


# use PC 1-15 for clustering
data_seurat <- FindNeighbors(data_seurat, dims = 1:15)
data_seurat <- FindClusters(data_seurat, resolution = 0.05)

# run UMAP
data_seurat <- RunUMAP(data_seurat, dims = 1:15)

DimPlot(data_seurat, reduction = "umap",
        group.by= "Phase")
ggsave("results/plots/rna_clustering/umap_cellcycle.svg", height = 7, width = 7)

DimPlot(data_seurat, reduction = "umap")
ggsave("results/plots/rna_clustering/umap.svg", height = 7, width = 7)

# calculate DEGs
marker_genes <- FindAllMarkers(data_seurat, logfc.threshold = 0.1) %>%
  as_tibble() %>% filter(cluster == 1) %>% filter(p_val_adj < 0.05)

# create featureplot
FeaturePlot(data_seurat, features = marker_genes %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 12) %>% pull(gene))




