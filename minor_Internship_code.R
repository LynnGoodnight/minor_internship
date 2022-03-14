###################################  Loading single cell RNA-seq results of #########
# Feature-Barcode Matrices format  into R
library(Matrix)
matrix_dir = "D:/minor_intership_data/minor Internship in house scRNA-seq Human_lymph_node_stromal_cells_data/Human_lymph_node_stromal_cells_data/minusDataset/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

################################### install Seurat#########
# Enter commands in R (or R studio, if installed)
install.packages('Seurat')



library(dplyr)
library(Seurat)
library(patchwork)

################################### Seurat - Guided Clustering Tutorial#########
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "D:/minor_intership_data/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

################################### import raw and processed&annotated in-house #########
# lymph nodes sc-RNA seq data to seurat object
import_in_house_sc <- function()
{
  library(dplyr)
  library(Seurat)
  library(patchwork)
  
  # in_house_minus.data <- Read10X(data.dir = "D:/minor_intership_data/minor Internship in house scRNA-seq Human_lymph_node_stromal_cells_data/Human_lymph_node_stromal_cells_data/minusDataset")
  # in_house_minus <- CreateSeuratObject(counts = in_house_minus.data, project = "in_house_minus", min.cells = 3, min.features = 200)
  # in_house_minus
  # 
  # in_house_plus.data <- Read10X(data.dir = "D:/minor_intership_data/minor Internship in house scRNA-seq Human_lymph_node_stromal_cells_data/Human_lymph_node_stromal_cells_data/plusDataset")
  # in_house_plus <- CreateSeuratObject(counts = in_house_plus.data, project = "in_house_plus", min.cells = 3, min.features = 200)
  # in_house_plus
  
  
  in_house_sc <- readRDS("D:/minor_intership_data/5_sc.data.rds")
  head(in_house_sc@meta.data$blue.main, 2)
  return(in_house_sc)
}
in_house_sc = import_in_house_sc()
# Show QC metrics for the first 5 cells
head(in_house_sc@meta.data, 5)
table(reference@cell_types) 
in_house_sc_sub <- subset(in_house_sc, subset = blue.main == 'Adipocytes' | 
                            blue.main == 'Endothelial cells' |
                            blue.main == 'Astrocytes'| 
                            blue.main == 'B-cells' |
                            blue.main == 'Fibroblasts' |
                            blue.main == 'Myocytes' |
                            blue.main == 'Skeletal muscle'  )

################################### import cell2location and Tabula h5ad to seurat object --------

# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
# devtools::install_github('satijalab/seurat-data')

library(Seurat)
library(SeuratData)
library(SeuratDisk)


# url <- "https://seurat.nygenome.org/pbmc3k_final.h5ad"
# curl::curl_download(url, basename(url))
# Convert("pbmc3k_final.h5ad", dest = "h5seurat", overwrite = TRUE)
# pbmc3k <- LoadH5Seurat("pbmc3k_final.h5seurat")
# pbmc3k
# #> An object of class Seurat 
# #> 13714 features across 2638 samples within 1 assay 
# #> Active assay: RNA (13714 features, 0 variable features)
# #>  2 dimensional reductions calculated: pca, umap

import_cell2location_sc <- function()
{
  Convert("D:/minor_intership_data/sc.h5ad", dest = "h5seurat", overwrite = TRUE)
  cell2location_sc <- LoadH5Seurat("D:/minor_intership_data/sc.h5seurat",assays  = "RNA")
  cell2location_sc
  head(cell2location_sc@meta.data, 2)
  saveRDS(cell2location_sc,'cell2location_sc.rds')
  return(cell2location_sc)
}
cell2location_sc = import_cell2location_sc()
rds_file='cell2location_sc.rds'
cell2location_sc <- readRDS(rds_file)


import_Tabula_sc <- function()
{
  Convert("D:/minor_intership_data/TS_Lymph_Node.h5ad", dest = "h5seurat", overwrite = TRUE)
  Tabula_sc <- LoadH5Seurat("D:/minor_intership_data/TS_Lymph_Node.h5seurat",assays  = "RNA")
  Tabula_sc
  head(Tabula_sc@meta.data, 5)
  saveRDS(Tabula_sc,'Tabula_ref.rds')
  return(Tabula_sc)
}
Tabula_sc = import_Tabula_sc()
rds_file='Tabula_ref.rds'
Tabula_sc <- readRDS(rds_file)

# head(Matrix::rowSums(cell2location_sc@assays$RNA))
# head(Matrix::rowSums(in_house_sc@assays$RNA))
# head(Matrix::rowSums(Tabula_sc@assays$RNA))
# 
# # check how many features are zero in all cells
# sum(Matrix::rowSums(Tabula_sc@assays$RNA)==0)

################################### integrate two scRNA-seq data sets: in-house + Tabula_sc datasets --------

library(Seurat)
library(SeuratData)
library(patchwork)

# install dataset
InstallData("ifnb")

# # load dataset
# LoadData("ifnb")
# 
# # split the dataset into a list of two seurat objects (stim and CTRL)
# ifnb.list <- SplitObject(ifnb, split.by = "stim")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = list(Tabula_sc=Tabula_sc,in_house_sc=in_house_sc), FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)


# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(immune.combined, reduction = "umap", split.by = "stim")
# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
head(nk.markers)

DietSeurat

###################################  RCTD tutorail  spatial transcriptomics vignette--------
# install.packages("devtools")
# devtools::install_github("dmcable/spacexr", build_vignettes = TRUE)
library(spacexr)
library(Matrix)

### Load in/preprocess your data, this might vary based on your file type
refdir <- system.file("extdata",'Reference/Vignette',package = 'spacexr') # directory for the reference
counts <- read.csv(file.path(refdir,"dge.csv")) # load in counts matrix
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
meta_data <- read.csv(file.path(refdir,"meta_data.csv")) # load in meta_data (barcodes, clusters, and nUMI)
cell_types <- meta_data$cluster; names(cell_types) <- meta_data$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- meta_data$nUMI; names(nUMI) <- meta_data$barcode # create nUMI named list

### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)
#> Warning in Reference(counts, cell_types, nUMI): Reference: nUMI does not match
#> colSums of counts. If this is unintended, please correct this discrepancy. If
#> this is intended, there is no problem.

## Examine reference object (optional)
print(dim(reference@counts)) #observe Digital Gene Expression matrix
#> [1] 384 475
table(reference@cell_types) #number of occurences for each cell type
#> 
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 
#> 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25

## Save RDS object (optional)
saveRDS(reference, file.path(refdir,'SCRef.rds'))

datadir <- system.file("extdata",'SpatialRNA/Vignette',package = 'spacexr') # directory for sample Slide-seq dataset
counts <- read.csv(file.path(datadir,"MappedDGEForR.csv")) # load in counts matrix
coords <- read.csv(file.path(datadir,"BeadLocationsForR.csv"))
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI

### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)

## Examine SpatialRNA object (optional)
print(dim(puck@counts)) # observe Digital Gene Expression matrix
hist(log(puck@nUMI,2)) # histogram of log_2 nUMI
print(head(puck@coords)) # start of coordinate data.frame
barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 

# This list can be restricted if you want to crop the puck e.g. 
# puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
# on the plot:
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI') 
myRCTD <- create.RCTD(puck, reference, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

#########################import Spatial transcriptomics dataset: 10X Visium Human Lymph Node and #########
# scRNA-seq dataset: Tabula Sapiens into RCTD
# 


init_require_packages <- function(){
  library(Seurat)
  library("readxl")
  library(SeuratData)
  library(patchwork)
  library(dplyr)
  library(RCTD)
  library(Matrix)
  library(ggplot2)
  library(SeuratDisk)
  library(spacexr)
  library(Matrix)
}

init_require_packages()

in_house_creat_ref_RCTD <- function(in_house_sc) {
  # load in counts matrix
  counts = as.matrix(x = GetAssayData(object = in_house_sc, assay = "RNA", slot = "counts"))
  # downsampled = SampleUMI(data = counts, max.umi = 1000, upsample = TRUE, verbose = TRUE)
  
  #rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
  meta_data <- in_house_sc@meta.data[,c('blue.main','nCount_RNA')]# load in meta_data (barcodes, clusters, and nUMI)
  cell_types <- meta_data$blue.main; names(cell_types) <- rownames(meta_data)#meta_data$barcode # create cell_types named list
  cell_types <- as.factor(cell_types) # convert to factor data type
  nUMI <- meta_data$nCount_RNA; names(nUMI) <- rownames(meta_data)#meta_data$barcode # create nUMI named list
  ### Create the Reference object
  reference <- Reference(counts, cell_types, nUMI)
  
  ## Examine reference object (optional)
  print(dim(reference@counts)) 
  # gene number; cell number
  table(reference@cell_types) 
  return(reference)
}

# reference = creat_ref_RCTD(in_house_sc)
# saveRDS(reference,'in_house_ref.rds')
# rds_file='in_house_ref.rds'
# reference <- readRDS(rds_file)

qc_filter_sc <- function(Tabula_sc){
  colnames(x = Tabula_sc[[]])
  Tabula_sc[['nCount_RNA']] = Tabula_sc$n_counts_UMIs
  
  # head(Tabula_sc@meta.data)
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  Tabula_sc[["percent.mt"]] <- PercentageFeatureSet(Tabula_sc, pattern = "^MT-")
  # # Visualize QC metrics as a violin plot
  # VlnPlot(Tabula_sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  Tabula_sc
  # filter cells based on QC metrics
  Tabula_sc <- subset(Tabula_sc, subset = n_genes > 200 & n_genes < 2500 & percent.mt < 5 )
  return(Tabula_sc)
  
}
# i filter cells that have gene number over 2,500 or less than 200
# i filter cells that have >5% mitochondrial counts
# these qc requirement comes from Seurat tutorial, and only retain 9532/53275 cells
Tabula_sc_filtered = qc_filter_sc(Tabula_sc)

subsample_cells <- function(Tabula_sc)
{
  Idents(Tabula_sc) <- 'free_annotation'
  table(Idents(Tabula_sc))
  unique(Idents(Tabula_sc))
  
  
  cell.list <- WhichCells(Tabula_sc, idents = unique(Idents(Tabula_sc)), 
                          downsample = 1000)
  Tabula_sc.downsampled <- Tabula_sc[, cell.list]
  table(Tabula_sc.downsampled$free_annotation)
  # saveRDS(Tabula_sc.downsampled,'Tabula_sc.downsampled.rds')
  # rds_file='Tabula_sc.downsampled.rds'
  # Tabula_sc.downsampled <- readRDS(rds_file)
  return(Tabula_sc.downsampled)
  
}
# subsample fixed numbers of cells from differently sized cell types in a seurat object.
Tabula_sc.downsampled = subsample_cells(Tabula_sc)

# Tabula_sc_filtered <- FindVariableFeatures(Tabula_sc_filtered, selection.method = "vst", nfeatures = 2000)

# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(Tabula_sc_filtered), 10)
# 
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(Tabula_sc_filtered)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# unsure whether a global-scaling normalization method ??LogNormalize?? is necesssary for FindVariableFeatures

############################################## import sc dataset to RCTD####
#find cell types that have less than 25 cells/cell type, delete these cell types and cells belong to them
delete_rare_cells <- function(Tabula_sc_filtered)
{
  colnames(x = Tabula_sc_filtered[[]])
  Idents(Tabula_sc_filtered) <- 'free_annotation'
  
  rare_cell_types = names(table(Idents(Tabula_sc_filtered))[!table(Idents(Tabula_sc_filtered)) > 25])
  unique(Idents(Tabula_sc_filtered))
  # delete cells belong to rare cell types
  for (char in rare_cell_types)
  {
    # remian cells whose cell type is not char
    Tabula_sc_filtered <- subset(Tabula_sc_filtered, subset = free_annotation != char)  
  }  
  
  table(Idents(Tabula_sc_filtered))
  
  table(Tabula_sc_filtered[['free_annotation']])
  return(Tabula_sc_filtered)
}
Tabula_sc_filtered = delete_rare_cells(Tabula_sc_filtered)
Create_Reference_object <- function(Tabula_sc_filtered)
{
  # load in counts matrix
  counts = as.matrix(Tabula_sc_filtered@assays$RNA@counts)
  
  # counts = as.matrix(x = GetAssayData(object = Tabula_sc.downsampled, assay = "RNA", slot = "counts"))
  # downsampled = SampleUMI(data = counts, max.umi = 1000, upsample = TRUE, verbose = TRUE)
  colnames(x = Tabula_sc_filtered[[]])
  #rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
  meta_data <- Tabula_sc_filtered@meta.data[,c("free_annotation",'nCount_RNA')]# load in meta_data (barcodes, clusters, and nUMI)
  cell_types <- meta_data$free_annotation; names(cell_types) <- rownames(meta_data)#meta_data$barcode # create cell_types named list
  cell_types <- as.factor(cell_types) # convert to factor data type
  cell_types <- droplevels(cell_types)
  nUMI <- meta_data$nCount_RNA; names(nUMI) <- rownames(meta_data)#meta_data$barcode # create nUMI named list
  ### Create the Reference object
  reference <- Reference(counts, cell_types, nUMI)
  
  ## Examine reference object (optional)
  print(dim(reference@counts)) 
  # gene number; cell number
  table(reference@cell_types) 
  return(reference)
}
reference = Create_Reference_object(Tabula_sc_filtered)

############################################## import spatial dataset to Seurat object####
setwd("C:/Users/12895/Documents/cell2location_Visium_RCTD")
rds_file='cell2location_sc_filtered.downsampled.rds'
cell2location_sc_filtered.downsampled <- readRDS(rds_file)



Create_Reference_object <- function(Tabula_sc_filtered)
{
  
  # load in counts matrix
  counts = as.matrix(Tabula_sc_filtered@assays$RNA@counts)
  
  # counts = as.matrix(x = GetAssayData(object = Tabula_sc.downsampled, assay = "RNA", slot = "counts"))
  # downsampled = SampleUMI(data = counts, max.umi = 1000, upsample = TRUE, verbose = TRUE)
  colnames(x = Tabula_sc_filtered[[]])
  #rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
  meta_data <- Tabula_sc_filtered@meta.data[,c("CellType",'nCount_RNA')]# load in meta_data (barcodes, clusters, and nUMI)
  cell_types <- meta_data$CellType
  cell_types = str_replace_all( cell_types,"/", "_") 
  names(cell_types) <- rownames(meta_data)#meta_data$barcode # create cell_types named list
  cell_types <- as.factor(cell_types) # convert to factor data type
  cell_types <- droplevels(cell_types)
  nUMI <- meta_data$nCount_RNA; names(nUMI) <- rownames(meta_data)#meta_data$barcode # create nUMI named list
  ### Create the Reference object
  reference <- Reference(counts, cell_types, nUMI)
  ## Examine reference object (optional)
  print(dim(reference@counts)) 
  # gene number; cell number
  table(reference@cell_types) 
  return(reference)
}



reference = Create_Reference_object(cell2location_sc_filtered.downsampled)











creat_spatial_RCTD <- function(data.dir){
  spatial = Load10X_Spatial( data.dir , filename = "V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5" )
  #https://github.com/dmcable/RCTD/issues/26
  # coords <- spatial@images$slice1@image
  coords <-GetTissueCoordinates(spatial, scale = NULL) 
  # pulling unscaled tissue coordinates
  #spatialcounts <- as.matrix(GetAssayData(spatial, assay = "Spatial", slot = "counts"))
  spatialcounts <- as.matrix(spatial@assays$Spatial@counts)
  ### Create SpatialRNA object
  puck <- SpatialRNA(coords, spatialcounts)
  saveRDS(puck,'spatial.rds')
  return(puck)
}
 
data.dir = 'D:/minor_intership_data/Visium_Spatial_Lymph_Node_Data/data'
puck = creat_spatial_RCTD(data.dir)
rds_file='spatial.rds'
puck <- readRDS(rds_file)

myRCTD <- create.RCTD(puck, reference, max_cores = 4)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full') #doublet_mode = 'doublet'
saveRDS(myRCTD,'myRCTD.rds')
rds_file='myRCTD.rds'
myRCTD  <- readRDS(rds_file)
# utils::methods(class = 'Assay')
show_results <- function(myRCTD)
{
  results <- myRCTD@results
  # normalize the cell type proportions to sum to 1.
  norm_weights = normalize_weights(results$weights) 
  head(norm_weights)
  
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- myRCTD@spatialRNA
  resultsdir <- 'cell2location_GeoMx_RCTD' ## you may change this to a more accessible directory on your computer.
  dir.create(resultsdir)
  #> Warning in dir.create(resultsdir): 'RCTD_Plots' already exists
  # make the plots 
  # Plots the confident weights for each cell type as in full_mode (saved as 
  # 'results/cell_type_weights_unthreshold.pdf')
  plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
  
  # Plots all weights for each cell type as in full_mode. (saved as 
  # 'results/cell_type_weights.pdf')
  plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) +geom_point(size = 10)

  
  
  
  # Plots the number of  pixels of each cell type where they confidently located in under 'full_mode'. (saved as 
  # 'results/cell_type_occur.pdf')
  plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
  return(norm_weights)
  
}
norm_weights_object = show_results(myRCTD)



creat_spatial_RCTD_GeoMx <- function()
{
  coords <- read_excel("D:/minor_intership_data/GeoMx WTA Human Lymph Node Data_count_results/Export4_NormalizationQ3.xlsx ",sheet = "SegmentProperties") 
  coords <- as.data.frame(coords[,c("SegmentDisplayName","ROICoordinateX","ROICoordinateY")]) 
  rownames(coords) = coords[,"SegmentDisplayName"]
  coords$SegmentDisplayName <- NULL
  counts <- read_excel("D:/minor_intership_data/GeoMx WTA Human Lymph Node Data_count_results/Export4_NormalizationQ3.xlsx ",sheet = "TargetCountMatrix") 
  counts <- as.data.frame(counts)
  rownames(counts) = counts[,1]
  counts$TargetName <- NULL
  puck <- SpatialRNA(coords, counts,require_int = FALSE)
  saveRDS(puck,'spatial_GeoMx.rds')
  return(puck)
}
puck = creat_spatial_RCTD_GeoMx()
rds_file='spatial_GeoMx.rds'
puck <- readRDS(rds_file)
myRCTD <- create.RCTD(puck, reference, max_cores = 4)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi') #doublet_mode = 'doublet'
saveRDS(myRCTD,'myRCTD3.rds')
rds_file='myRCTD3.rds'
myRCTD  <- readRDS(rds_file)

norm_weights_object = show_results(myRCTD)



# crutinize which cell types are delete after filtering, compare the cell types composition before and after qc filtering, see if this changed
#####
cell_type_table = as.data.frame(table(Tabula_sc[['free_annotation']]))
colnames(cell_type_table) = c("cell_type","before")
cell_type_table$after = as.data.frame(table(Tabula_sc_filtered[['free_annotation']]))[,2]
cell_type_table[,-1] = apply(cell_type_table[,-1],2,function(x) x/sum(x))
library(reshape2)
cell_type_table <- melt(cell_type_table, id.vars = "cell_type",
                        variable.name = "state",  value.name = "proportion")

ggplot(data = cell_type_table, mapping = aes(x = cell_type, 
                                             y = proportion, 
                                             colour = state,group = state )) + 
geom_point() +geom_line(aes(color=state))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

#####

# calculate percent of mitochondria gene
# Tabula_sc <- LoadH5Seurat(paste0(Data, 'Raw/TS_Lymph_Node.h5seurat'), assays = 'RNA')
mt.genes <- grep(pattern = '^MT-', x = rownames(x = Tabula_sc@assays$RNA),value = TRUE)

percent.mt <- Matrix::colSums(Tabula_sc@assays$RNA[mt.genes, ]) / Matrix::colSums(Tabula_sc@assays$RNA)

Tabula_sc <- AddMetaData(object = Tabula_sc,metadata = percent.mt,
                         col.name = 'percent.mt')


# Analysis, visualization, and integration tutorial of spatial datasets with Seurat####

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)



SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, 
                                                          idents = c(2, 1, 4, 3,5, 6)),
               facet.highlight = TRUE, ncol = 3)

cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2





allen_reference <- readRDS("allen_cortex.rds")
# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells, this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 300, verbose = FALSE)
%>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

# now remove additional cells, use SpatialDimPlots to visualize what to remove
# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
#change variable into mydata########################################################
plot1 <- VlnPlot(visium_spatial, features = "nCount_Spatial", pt.size = 0.1) + 
  NoLegend()
plot2 <- SpatialFeaturePlot(visium_spatial, features = "nCount_Spatial") +
  theme(legend.position = "right")
wrap_plots(plot1, plot2)

head(cell2location_sc_filtered.downsampled@assays$RNA@meta.features)
mt.genes <- grep(pattern = 'ENSG00000198804', 
                 x = cell2location_sc_filtered.downsampled@assays$RNA@meta.features[,1],value = TRUE)



######################################################### SPOTlight ####
# Installe SPOTlight
# version
install.packages("installr")
library(installr)
updateR()
BiocManager::install(version =
                       '3.14')
# install.packages("BiocManager")
# BiocManager::install("SPOTlight")# Or the devel version
# BiocManager::install("SPOTlight", version = "devel")

install.packages("devtools")
library(devtools)
install_github("https://github.com/MarcElosua/SPOTlight")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SpatialExperiment")
BiocManager::install("scater")
BiocManager::install("scran")
remove.packages("glue")
install.packages("glue")
devtools::install_github("fmicompbio/TabulaMurisSenisData")
BiocManager::install("TENxVisiumData")
devtools::install_github("kassambara/ggcorrplot")


#plot bar#####
cell_type_table = as.data.frame(table(sc_dataset_filtered$free_annotation))
colnames(cell_type_table) = c("cell_type","frequency")
library(ggplot2)
# Basic barplot
p<-ggplot(data=cell_type_table, aes(x=cell_type, y=frequency,fill=cell_type)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
p
browseVignettes("SpatialExperiment")
#####

dir <- system.file(
  file.path("extdata", "10xVisium"),
  package = "SpatialExperiment")

sample_ids <- c("section1", "section2")
samples <- file.path(dir, sample_ids)
spatial_dir = 'D:/minor_intership_data/Visium_Spatial_Lymph_Node_Data/data'



spatial <- read10xVisium(spatial_dir,
                         type = "HDF5", data = "raw",
                        images = "lowres",load = FALSE)




#spotlight tutorial older version######
install.packages("gt")
library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)


path_to_data <- system.file(package = "SPOTlight")
cortex_sc <- readRDS("cell2location_sc_dataset.downsampled.rds")


if (! "stxBrain" %in% SeuratData::AvailableData()[, "Dataset"]) {
  # If dataset not downloaded proceed to download it
  SeuratData::InstallData("stxBrain")
}

# Load data
anterior <- SeuratData::LoadData("stxBrain", type = "anterior1")



set.seed(123)
cortex_sc <- Seurat::SCTransform(cortex_sc, verbose = FALSE) %>%
  Seurat::RunPCA(., verbose = FALSE) %>%
  Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)

Seurat::DimPlot(cortex_sc,
                group.by = "CellType",
                label = TRUE) + Seurat::NoLegend()


cortex_sc@meta.data %>%
  dplyr::count(CellType) %>%
  gt::gt(.[-1, ]) %>%
  gt::tab_header(
    title = "Cell types present in the reference dataset",
  ) %>%
  gt::cols_label(
    subclass = gt::html("Cell Type")
  )



Seurat::Idents(object = cortex_sc) <- cortex_sc@meta.data$CellType
cluster_markers_all <- Seurat::FindAllMarkers(object = cortex_sc, 
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)

saveRDS(object = cluster_markers_all,
        "markers_sc.RDS")



set.seed(123)

spotlight_ls <- spotlight_deconvolution(
  se_sc = cortex_sc,
  counts_spatial = anterior@assays$Spatial@counts,
  clust_vr = "CellType", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 3000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

saveRDS(object = spotlight_ls, file = here::here("inst/spotlight_ls.rds"))

spotlight_ls <- readRDS(file = here::here("inst/spotlight_ls.rds"))

nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]



#! /usr/bin/Rscript
# screen -S spotlight_cell2location34
# Rscript run_spotlight.R &>nohup.out&

# .libPaths()
# 
# 
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("SpatialExperiment")

if(file.exists("Visium_spatial.rds")){
  
  rds_file='Visium_spatial.rds'
  spatial <- readRDS(rds_file)
}else{
  # data.dir = 'D:/minor_intership_data/Visium_Spatial_Lymph_Node_Data/data'
  # 
  # spatial = Load10X_Spatial( data.dir , filename =
  #                              "V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5" )
  # saveRDS(spatial,"Visium_spatial.rds")
}




###############################
#' @rdname plotTopicProfiles
#' @name plotTopicProfiles
#' @title Plot NMF topic profiles
#'
#' @description This function takes in the fitted NMF model and returns the
#'   topic profiles learned for each cell \code{facet = FALSE} or cell type
#'   \code{facet = TRUE}. Ideal training will return all the cell from the same
#'   cell type to share a unique topic profile.
#'
#' @param x \code{\link{NMFfit}} object
#' @param y vector of group labels. Should be of length \code{ncol(coef(x))}.
#' @param facet logical indicating whether to stratify by group.
#'   If \code{FALSE} (default), weights will be the median across cells
#'   for each group (point = topic weight for a given cell type).
#'   If \code{TRUE}, cell-specific weights will be shown
#'   (point = topic weight of a given cell).
#' @param min_prop scalar in [0,1]. When \code{facet = TRUE},
#'   only cells with a weight > \code{min_prop} will be included.
#' @param ncol integer scalar specifying the number of facet columns.
#' @param ... additional parameters
#' 
#' @return \code{ggplot} object
#'
#' @author Marc Elosua Bayes & Helena L Crowell
#'
#' @examples
#' library(ggplot2)
#' x <- mockSC()
#' y <- mockSP(x)
#' z <- getMGS(x)
#' 
#' res <- SPOTlight(x, y,
#'     groups = x$type,
#'     mgs = z,
#'     group_id = "type",
#'     verbose = FALSE)
#'
#' plotTopicProfiles(res[[3]], x$type, facet = TRUE)
#' plotTopicProfiles(res[[3]], x$type, facet = FALSE)
NULL

# try to convert anything to character
# (e.g., factor or numeric labels as input)
#' @rdname plotTopicProfiles
#' @export
setMethod(
  "plotTopicProfiles", c("NMFfit", "ANY"),
  function(x, y, ...) plotTopicProfiles(x, as.character(y), ...)
)

#' @rdname plotTopicProfiles
#' @importFrom NMF coef
#' @importFrom stats aggregate
#' @import ggplot2
#' @export
setMethod(
  "plotTopicProfiles",
  c("NMFfit", "character"),
  function(x, y,
           facet = FALSE,
           min_prop = 0.1,
           ncol = NULL) {
    
    # check validity of input arguments
    stopifnot(
      methods::is(x, "NMF"),
      length(y) == ncol(coef(x)),
      setequal(rownames(coef(x)), y),
      is.logical(facet), length(facet) == 1,
      is.numeric(min_prop), length(min_prop) == 1)
    # x = mod
    # y = cell_type
    # min_prop = 0.3 
    facet = TRUE
    y <- as.character(y)
    stopifnot(y %in% colnames(basis(x)))

    # get group proportions
    mat <- prop.table(t(coef(x)), 1)
    facet = TRUE
    if (facet) {
      # stretch for plotting
      df <- data.frame(
        id = seq_len(nrow(mat)),
        weight = c(mat),
        group = rep(y, ncol(mat)),#21088 cells, repeat 34 times
        topic = rep(seq_len(ncol(mat)), each = nrow(mat)))
      #34 topics, repeat 21088 times

      # drop cells with 'weight < min_prop'
      df <- df[df$weight >= 0.01, ]
      
      # set aesthetics
      x <- "id"
      f <- facet_wrap(~group, ncol = ncol, scales = "free_x")
    } else {
      # get topic medians
      df <- aggregate(mat, list(y), sum)
      group = df[,1]
      df <- df[,-1]
      
      # stretch for plotting
      df <- data.frame(
        weight = unlist(df),
        group = rep(group, ncol(df) ),
        topic = rep(seq_len(nrow(df)), each = nrow(df))
      )
      
      # set aesthetics
      x <- "group"
      f <- NULL
    }
    
    
    # fix topic order
    df$topic <- factor(df$topic, seq_along(unique(y)))
    
    # render plot
    ggplot(df, aes_string(x, "topic",
                           size = "weight")) +
     geom_point()  +
      guides(col = guide_legend(override.aes = list(size = 2))) +
      scale_size_continuous(range = c(0, 3)) +
      scale_color_continuous(low = "lightgrey", high = "blue") +
      xlab(if (facet) x) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1))
  }
)
# location = "server"
# # location = "my_computer"
# if(location == "my_computer"){
#   setwd("C:/Users/12895/Documents")
# 
# }else if(location == "server"){
#   setwd("/data/home/lyang")
#   resultsdir <- 'Visium_RCTD' 
#   dir.create(resultsdir)
#   setwd("/data/home/lyang/Visium_RCTD")
# }
##########################