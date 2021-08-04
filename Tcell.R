# integration Tcells
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
library(patchwork)
library(stringr)
library(future)
library(viridis)
library(rcartocolor)

setwd("/beegfs/scratch/ric.cosr/ric.alfano/AlfanoM_904_infertilita_epigenetics/Tcell")
outdir = "SeuratOutputs/"
dir.create(outdir, recursive = T)

# load Matrix samples
path_file = "./GSE126030_RAW/"
fl = list.files(path_file, pattern=glob2rx("*MOD.txt"))
input_datasets = str_sub(fl, 1, 10)
input_datasets
metadata = data.frame(Dataset = input_datasets,
                      Tissue  = c("Lung", "Lung", "BoneMarrow", "BoneMarrow", "LymphNode", "LymphNode",
                                  "Lung", "Lung", "BoneMarrow", "BoneMarrow", "LymphNode", "LymphNode",
                                  "Blood", "Blood", "Blood", "Blood"),
                      Stim    = c("Resting", "Activated","Resting", "Activated","Resting", "Activated",
                                  "Resting", "Activated","Resting", "Activated","Resting", "Activated",
                                   "Activated","Resting", "Activated", "Resting"),
                      Donor   = c("Donor1", "Donor1", "Donor1", "Donor1", "Donor1","Donor1",
                                  "Donor2", "Donor2", "Donor2", "Donor2", "Donor2","Donor2",
                                  "BloodDonorA", "BloodDonorA", "BloodDonorB","BloodDonorB"))
row.names(metadata) <- metadata$Dataset
metadata


first_time = FALSE
if (first_time) {
  Tcell.list = list()
  
  start_time <- Sys.time()
  for (i in 1:length(fl)) {
    dataset = str_sub(fl[i], 1, 10)
    print(dataset)
    # read file
    matrix <- read.table(paste(path_file,fl[i], sep=''),
                         sep="\t", header = TRUE, quote ="", 
                         row.names = 1)
    # remove gene column
    matrix_s <- as.sparse(matrix[,2:dim(matrix)[2]])
    # from ensembl to gene name
    row.names(matrix_s) <-  matrix$Gene
    # seurat object
    tcell <- CreateSeuratObject(counts       = matrix_s, 
                                project      = dataset, 
                                min.cells    = 5, 
                                min.features = 200)
    # add metadata
    tcell$Tissue = metadata[dataset, "Tissue"]
    tcell$Stim = metadata[dataset, "Stim"]
    tcell$Donor = metadata[dataset, "Donor"]
    # save object on a list
    Tcell.list[[i]] = tcell
  }
  end_time <- Sys.time()
  end_time - start_time
  
  Tcell.list
  names(Tcell.list) = metadata$Dataset
  save(Tcell.list, file = paste(outdir, "Tcell.list", sep=''))
} else {load( paste(outdir, "Tcell.list", sep=''))}


# ---------- ALL T-cell ------------------------------
TCELL = merge(x = Tcell.list[[1]], 
              y = Tcell.list[2:length(fl)],
              add.cell.ids = metadata$Dataset, 
              project = "GSE126030")
# MT%
TCELL[["percent.mt"]] <- PercentageFeatureSet(TCELL, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(TCELL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Normalize data
TCELL <- NormalizeData(TCELL, normalization.method = "LogNormalize", 
                       scale.factor = 10000)
# Find Variable Features
TCELL <- FindVariableFeatures(TCELL, selection.method = "vst", nfeatures = 2000)

# Scale Data
all.genes <- rownames(TCELL)
options(future.globals.maxSize = 3145728000) 
plan(strategy = "multicore", workers = 6)
TCELL <- ScaleData(TCELL, features = all.genes)

# PCA
TCELL <- RunPCA(TCELL, features = VariableFeatures(object = TCELL))
ElbowPlot(TCELL)

# cluster the cells 
nPC = 20; res = 0.2
TCELL <- FindNeighbors(TCELL, dims = 1:nPC)
TCELL <- FindClusters(TCELL, resolution = res)

# UMAP
TCELL <- RunUMAP(TCELL, dims = 1:nPC)
DimPlot(TCELL, reduction = "umap", group.by = "Donor")
DimPlot(TCELL, reduction = "umap", group.by = "Tissue")
DimPlot(TCELL, reduction = "umap")

save(TCELL, file = paste(outdir, "TCELL_GSE126030", sep=''))

TCELL.markers <- FindAllMarkers(TCELL, only.pos = TRUE, min.pct = 0.25, 
                                logfc.threshold = 0.25)
top20 = TCELL.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
top20[top20$cluster == '6',]
write.csv(TCELL.markers, file = paste(outdir, "ClusterMarker_allGSE.csv", sep=''), quote = F, row.names = F)
HM_top10 = DoHeatmap(TCELL, features = top10$gene) + NoLegend()
HM_top10

TCELL.markers[grep("PDZD8", TCELL.markers$gene),]

# ---------------- DONOR1 Tissue only ---------------
dat_subset = metadata[metadata$Donor == "Donor1", "Dataset"]
Tcell.list_Donor1 = Tcell.list[dat_subset]

TCELL_D1 = merge(x = Tcell.list_Donor1[[1]], 
                 y = Tcell.list_Donor1[2:length(dat_subset)],
                 add.cell.ids = dat_subset, 
                 project = "GSE126030")
TCELL_D1
# MT%
TCELL_D1[["percent.mt"]] <- PercentageFeatureSet(TCELL_D1, pattern = "^MT-")
#Normalize data
TCELL_D1 <- NormalizeData(TCELL_D1, normalization.method = "LogNormalize", scale.factor = 10000)
# Find Variable Features
TCELL_D1 <- FindVariableFeatures(TCELL_D1, selection.method = "vst", nfeatures = 2000)
# Scale Data
all.genes <- rownames(TCELL_D1)
TCELL_D1 <- ScaleData(TCELL_D1, features = all.genes)
# PCA
TCELL_D1 <- RunPCA(TCELL_D1, features = VariableFeatures(object = TCELL_D1))
ElbowPlot(TCELL_D1)
# cluster the cells 
nPC = 10; res = 0.3
TCELL_D1 <- FindNeighbors(TCELL_D1, dims = 1:nPC)
TCELL_D1 <- FindClusters(TCELL_D1, resolution = res)
# UMAP
TCELL_D1 <- RunUMAP(TCELL_D1, dims = 1:nPC)
u = DimPlot(TCELL_D1, reduction = "umap")
t = DimPlot(TCELL_D1, reduction = "umap", group.by = "Tissue")
s = DimPlot(TCELL_D1, reduction = "umap", group.by = "Stim")
TCELL_D1$Tissue.Stim = paste(TCELL_D1$Tissue, TCELL_D1$Stim, sep='_')
t.s = DimPlot(TCELL_D1, reduction = "umap", group.by = "Tissue.Stim")
TCELL_D1[["CD4_CD8"]] <- log2((FetchData(object = TCELL_D1, vars = "CD4")+.01)/(FetchData(object = TCELL_D1, vars = "CD8A")+.01))
cd = FeaturePlot(TCELL_D1, features = "CD4_CD8", min.cutoff = -4, max.cutoff = 4, cols = c("#00CC00", "#FF00FF"))

(u+t.s) / (s+cd)
ggsave(paste(outdir, "DimPLot_Donor1.png", sep=''), width = 10, height = 8)

# find markers for every cluster compared to all remaining cells, report only the positive ones
TCELL_D1.markers <- FindAllMarkers(TCELL_D1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 = TCELL_D1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10
write.csv(TCELL_D1.markers, file = paste(outdir, "ClusterMarker_D1.csv", sep=''), quote = F, row.names = F)
HM_top10 = DoHeatmap(TCELL_D1, features = top10$gene) + NoLegend()
HM_top10
ggsave(paste(outdir, "HM_Donor1.png", sep=''), HM_top10, width = 15, height = 20)

save(TCELL_D1, file = paste(outdir, "TCELL_D1", sep=''))

#myTCL = get(load("OSR_TCL"))

# ---------------- DONOR2 Tissue only ---------------
dat_subset = metadata[metadata$Donor == "Donor2", "Dataset"]
Tcell.list_Donor2 = Tcell.list[dat_subset]

TCELL_D2 = merge(x = Tcell.list_Donor2[[1]], 
                 y = Tcell.list_Donor2[2:length(dat_subset)],
                 add.cell.ids = dat_subset, 
                 project = "GSE126030")
TCELL_D2
# MT%
TCELL_D2[["percent.mt"]] <- PercentageFeatureSet(TCELL_D2, pattern = "^MT-")
#Normalize data
TCELL_D2 <- NormalizeData(TCELL_D2, normalization.method = "LogNormalize", scale.factor = 10000)
# Find Variable Features
TCELL_D2 <- FindVariableFeatures(TCELL_D2, selection.method = "vst", nfeatures = 2000)
# Scale Data
all.genes <- rownames(TCELL_D2)
TCELL_D2 <- ScaleData(TCELL_D2, features = all.genes)
# PCA
TCELL_D2 <- RunPCA(TCELL_D2, features = VariableFeatures(object = TCELL_D2))
ElbowPlot(TCELL_D2)
# cluster the cells 
nPC = 10; res = 0.3
TCELL_D2 <- FindNeighbors(TCELL_D2, dims = 1:nPC)
TCELL_D2 <- FindClusters(TCELL_D2, resolution = res)
# UMAP
TCELL_D2 <- RunUMAP(TCELL_D2, dims = 1:nPC)
u = DimPlot(TCELL_D2, reduction = "umap")
t = DimPlot(TCELL_D2, reduction = "umap", group.by = "Tissue")
s = DimPlot(TCELL_D2, reduction = "umap", group.by = "Stim")
TCELL_D2$Tissue.Stim = paste(TCELL_D2$Tissue, TCELL_D2$Stim, sep='_')
t.s = DimPlot(TCELL_D2, reduction = "umap", group.by = "Tissue.Stim")
TCELL_D2[["CD4_CD8"]] <- log2((FetchData(object = TCELL_D2, vars = "CD4")+.01)/(FetchData(object = TCELL_D2, vars = "CD8A")+.01))
cd = FeaturePlot(TCELL_D2, features = "CD4_CD8", min.cutoff = -4, max.cutoff = 4, cols = c("#00CC00", "#FF00FF"))

(u+t.s) / (s+cd)
ggsave(paste(outdir, "DimPLot_Donor2.png", sep=''), width = 10, height = 8)

# find markers for every cluster compared to all remaining cells, report only the positive ones
TCELL_D2.markers <- FindAllMarkers(TCELL_D2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 = TCELL_D2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10
write.csv(TCELL_D2.markers, file = paste(outdir, "ClusterMarker_D2.csv", sep=''), quote = F, row.names = F)
HM_top10 = DoHeatmap(TCELL_D2, features = top10$gene) + NoLegend()
HM_top10
ggsave(paste(outdir, "HM_Donor2.png", sep=''), HM_top10, width = 15, height = 20)
saveRDS(TCELL_D2, file = paste(outdir, "TCELL_D2.rds", sep=''))

# ----------------------- Project OSR cells -----------------------------
myTCL = get(load("TCL_OSR"))
myTCL$Donor = "iGCA"
myTCL$Tissue = "Testis"
myTCL$Stim = "NA"

dat_subset = metadata[metadata$Donor == "Donor1", "Dataset"]
Tcell.list_Donor1 = Tcell.list[dat_subset]

TCELL_D1_myTCL = merge(x = myTCL, 
                       y = Tcell.list_Donor1,
                       add.cell.ids = c("iNOA_TCL", dat_subset), 
                       project = "GSE126030")
TCELL_D1_myTCL
# MT%
TCELL_D1_myTCL[["percent.mt"]] <- PercentageFeatureSet(TCELL_D1_myTCL, pattern = "^MT-")
#Normalize data
TCELL_D1_myTCL <- NormalizeData(TCELL_D1_myTCL, normalization.method = "LogNormalize", scale.factor = 10000)
# Find Variable Features
TCELL_D1_myTCL <- FindVariableFeatures(TCELL_D1_myTCL, selection.method = "vst", nfeatures = 2000)
# Scale Data
all.genes <- rownames(TCELL_D1_myTCL)
TCELL_D1_myTCL <- ScaleData(TCELL_D1_myTCL, features = all.genes)
# PCA
TCELL_D1_myTCL <- RunPCA(TCELL_D1_myTCL, features = VariableFeatures(object = TCELL_D1_myTCL))
ElbowPlot(TCELL_D1_myTCL)
# cluster the cells 
nPC = 10; res = 0.3
TCELL_D1_myTCL <- FindNeighbors(TCELL_D1_myTCL, dims = 1:nPC)
TCELL_D1_myTCL <- FindClusters(TCELL_D1_myTCL, resolution = res)
# UMAP
TCELL_D1_myTCL <- RunUMAP(TCELL_D1_myTCL, dims = 1:nPC)
u = DimPlot(TCELL_D1_myTCL, reduction = "umap")
t = DimPlot(TCELL_D1_myTCL, reduction = "umap", group.by = "Tissue")
s = DimPlot(TCELL_D1_myTCL, reduction = "umap", group.by = "Stim")
TCELL_D1_myTCL$Tissue.Stim = paste(TCELL_D1_myTCL$Tissue, TCELL_D1_myTCL$Stim, sep='_')
t.s = DimPlot(TCELL_D1_myTCL, reduction = "umap", group.by = "Tissue.Stim")
TCELL_D1_myTCL[["CD4_CD8"]] <- log2((FetchData(object = TCELL_D1_myTCL, vars = "CD4")+.01)/(FetchData(object = TCELL_D1_myTCL, vars = "CD8A")+.01))
cd = FeaturePlot(TCELL_D1_myTCL, features = "CD4_CD8", min.cutoff = -4, max.cutoff = 4, cols = c("#00CC00", "#FF00FF"))

(u+t.s) / (s+cd)
ggsave(paste(outdir, "DimPLot_Donor1_myTCL.png", sep=''), width = 10, height = 8)

# find markers for every cluster compared to all remaining cells, report only the positive ones
TCELL_D1_myTCL.markers <- FindAllMarkers(TCELL_D1_myTCL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 = TCELL_D1_myTCL.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10
write.csv(TCELL_D1_myTCL.markers, file = paste(outdir, "ClusterMarker_D1.csv", sep=''), quote = F, row.names = F)
HM_top10 = DoHeatmap(TCELL_D1_myTCL, features = top10$gene) + NoLegend()
HM_top10
ggsave(paste(outdir, "HM_Donor1.png", sep=''), HM_top10, width = 15, height = 20)
save(TCELL_D1_myTCL, file = paste(outdir, "TCELL_D1_myTCL", sep=''))

#myTCL = get(load("OSR_TCL"))


# --------------- integrate my Tcell ---------------------
load(paste(outdir, "TCELL_D1.rds", sep=''))

samples.list = list(TCELL_D1, myTCL)

# normalize and identify variable features for each dataset independently
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = samples.list)
tcell.anchors <- FindIntegrationAnchors(object.list = samples.list, anchor.features = features,
                                        k.filter = 90)
tcell.combined <- IntegrateData(anchorset = tcell.anchors)

DefaultAssay(tcell.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
tcell.combined <- ScaleData(tcell.combined, verbose = FALSE)
tcell.combined <- RunPCA(tcell.combined, npcs = 30, verbose = FALSE)
tcell.combined <- RunUMAP(tcell.combined, reduction = "pca", dims = 1:30)
tcell.combined <- FindNeighbors(tcell.combined, reduction = "pca", dims = 1:30)
tcell.combined <- FindClusters(tcell.combined, resolution = 0.5)

u = DimPlot(tcell.combined, reduction = "umap")
t = DimPlot(tcell.combined, reduction = "umap", group.by = "Tissue", order = TRUE, pt.size = 2)
s = DimPlot(tcell.combined, reduction = "umap", group.by = "Stim", order = TRUE, pt.size = 2)
tcell.combined$Tissue.Stim = paste(tcell.combined$Tissue, tcell.combined$Stim, sep='_')
DefaultAssay(tcell.combined) <- "RNA"
t.s = DimPlot(tcell.combined, reduction = "umap", group.by = "Tissue.Stim", order = TRUE, pt.size = 2)
tcell.combined[["CD4_CD8"]] <- log2((FetchData(object = tcell.combined, vars = "CD4")+.01)/(FetchData(object = tcell.combined, vars = "CD8A")+.01))
cd = FeaturePlot(tcell.combined, features = "CD4_CD8", min.cutoff = -4, max.cutoff = 4, cols = c("#00CC00", "#FF00FF"))

(u+t.s) / (s+cd)
ggsave(paste(outdir, "DimPLot_Donor1_myTCL.png", sep=''), width = 10, height = 8)



# --------------- Integration ----------------
head(myTCL[[]])
myTCL
TCELL = merge(x = myTCL, 
              y = Tcell.list,
              add.cell.ids = c("OSR", metadata$Dataset), 
              project = "Tcell")
# MT%
TCELL[["percent.mt"]] <- PercentageFeatureSet(TCELL, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(TCELL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(TCELL, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TCELL, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# object is subsetted already

#Normalize data
TCELL <- NormalizeData(TCELL, normalization.method = "LogNormalize", 
                       scale.factor = 10000)

# Find Variable Features
TCELL <- FindVariableFeatures(TCELL, selection.method = "vst", nfeatures = 2000)

# Scale Data
all.genes <- rownames(TCELL)
options(future.globals.maxSize = 3145728000) 
plan(strategy = "multicore", workers = 4)
TCELL <- ScaleData(TCELL)

# PCA
TCELL <- RunPCA(TCELL, features = VariableFeatures(object = TCELL))
DimPlot(TCELL, reduction = "pca")
ElbowPlot(TCELL)

# cluster the cells 
nPC = 20; res = 0.1
TCELL <- FindNeighbors(TCELL, dims = 1:nPC)
TCELL <- FindClusters(TCELL, resolution = res)




TCELL <- RunUMAP(TCELL, dims = 1:nPC)
DimPlot(TCELL, reduction = "umap", group.by = "Tissue")
u = DimPlot(TCELL, reduction = "umap", cols = col)
t = DimPlot(TCELL, reduction = "umap", group.by = "Tissue")
s = DimPlot(TCELL, reduction = "umap", group.by = "Stim")
TCELL$Tissue.Stim = paste(TCELL$Tissue, TCELL$Stim, sep='_')
t.s = DimPlot(TCELL, reduction = "umap", group.by = "Tissue.Stim")
TCELL[["CD4_CD8"]] <- log2((FetchData(object = TCELL, vars = "CD4")+.001)/(FetchData(object = TCELL, vars = "CD8A")+.001))
cd = FeaturePlot(TCELL, features = "CD4_CD8", min.cutoff = -4, max.cutoff = 4, cols = c("#00CC00", "#FF00FF"))

(u+t.s) / (s+cd)
ggsave(paste(outdir, "DimPLot_ALL.png", sep=''), width = 10, height = 8)
TCELL.markers <- FindAllMarkers(TCELL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 = TCELL.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10
write.csv(TCELL.markers, file = paste(outdir, "ClusterMarker_D2.csv", sep=''), quote = F, row.names = F)
HM_top10 = DoHeatmap(TCELL, features = top10$gene) + NoLegend()
HM_top10
ggsave(paste(outdir, "HM_Tcell_all_osr.png", sep=''), HM_top10, width = 15, height = 20)
save(TCELL, file = paste(outdir, "TCELL", sep=''))


# --------------- label Transfer my Tcell ----------------
load(paste(outdir, "TCELL_new", sep=''))
new.labels = c("CD4 act 1", 
               "CD8 EM/TRM rest", 
               "CD4 rest 1", 
               "CD8 EM/TRM act", 
               "CD4 Treg",
               "CD8 TRM rest",
               "CD8 TRM act",
               "CD4 act 2",
               "CD4 rest 2",
               "CD8 TEMRA")

names(new.labels) <- levels(TCELL_new)
TCELL_new <- RenameIdents(TCELL_new, new.labels)
DimPlot(TCELL_new)


nColor <- length(new.labels)
col = (carto_pal(nColor, "Bold")); names(col) <- new.labels
DimPlot(TCELL_new, label = T)
u = DimPlot(TCELL_new, reduction = "umap", label = F, cols = col)
colT = rev((carto_pal(length(table(TCELL_new$Tissue)), "Geyser")))
names(colT) <- names(table(TCELL_new$Tissue))
t = DimPlot(TCELL_new, reduction = "umap", group.by = "Tissue", cols = colT)
s = DimPlot(TCELL_new, reduction = "umap", group.by = "Stim", 
            cols = c("#FF0033", "#0066FF"))
TCELL_new[["CD4_CD8"]] <- log2((FetchData(object = TCELL_new, vars = "CD4")+.05)/(FetchData(object = TCELL_new, vars = "CD8A")+.05))
cd = FeaturePlot(TCELL_new, features = "CD4_CD8", min.cutoff = -3, max.cutoff = 3, cols = c("#00CC00", "#FF00FF"))
cd
(u + t)/(s + cd)

ggsave(paste(outdir, "GSE126030_reanalysis.png", sep =''), width = 15, height = 10)
pdf(paste(outdir, "GSE126030_reanalysis.pdf", sep =''), width = 15, height = 10)
print((u + t)/(s + cd))
dev.off()

TCELL_new

save(TCELL_new, file= paste(outdir, "TCELL_new", sep=''))

# ------ predicition ---- 

TCELL_new$cell.ont <- Idents(TCELL_new)

anchors <- FindTransferAnchors(reference=TCELL_new, query=myTCL, k.filter = 90,dims=1:30)
predictions <- TransferData(anchorset=anchors, 
                            refdata=TCELL_new$cell.ont, dims=1:30)
predictions_Tissue <- TransferData(anchorset=anchors, 
                                   refdata=TCELL_new$Tissue, dims=1:30)
predictions_Stim <- TransferData(anchorset=anchors, 
                                 refdata=TCELL_new$Stim, dims=1:30)
table(predictions_Stim$predicted.id)

myTCL <- AddMetaData(myTCL, metadata=predictions)
myTCL$pred.Tissue = predictions_Tissue$predicted.id

table(myTCL$predicted.id)

DimPlot(myTCL, group.by = "predicted.id")
pt = length(table(predictions$predicted.id))

RedoUMAP = TRUE

if (!RedoUMAP) {
  p = DimPlot(myTCL, group.by = "predicted.id", pt.size = 4, cols = col) + 
    ggtitle("Projection of Tcell identity in Szabo et al.") + 
    theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=12, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 14),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
          axis.title.y = element_text(face = "bold", color = "black", size = 14),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          #legend.position=c(0.7,0.8),
          legend.position = c(1, 0),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.box.background = element_rect(color="black", fill = "#FFFFFF", size=1),
          legend.margin = margin(6, 6, 6, 6),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) 
  print(p) 
  umapr = "zoomUMAP"
  } else {
  myTCL2 = ScaleData(myTCL)
  myTCL2 = RunUMAP(myTCL2, cells.use = NULL, dims = 1:10, reduction.use = "pca")
  DimPlot(myTCL2)
  p = DimPlot(myTCL2, group.by = "predicted.id", pt.size = 6, cols = col) + 
    ggtitle("Projection of Tcell identity in Szabo et al.") + 
    theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=12, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 14),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
          axis.title.y = element_text(face = "bold", color = "black", size = 14),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = "top",
          #legend.position=c(0.7,0.8),
          #legend.position = c(1, 0),
          #legend.justification = c("right", "bottom"),
          #legend.box.just = "right",
          #legend.box.background = element_rect(color="black", fill = "#FFFFFF", size=1),
          #legend.margin = margin(6, 6, 6, 6),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) 
  p 
  umapr = "redoUMAP"
  
  fp_blend <- FeaturePlot(myTCL2, 
                          reduction = "umap", 
                          pt.size = 4,
                          features = c('CD8A','CD69'),
                          label = F, 
                          repel = T,
                          order= TRUE,
                          blend.threshold = 0.5,
                          #slot = 'counts',
                          blend = T,       
                          max.cutoff = 'q85',
                          cols = c("lightgrey", "red","green")) + theme(legend.position = 'right')
  
  pdf(paste(outdir, "Figure2_Tcell_CD8_CD69_blend_tcellonly_",umapr,".pdf",sep =''), width = 12, height = 4)
  print(fp_blend)
  dev.off()
  
  fp_blend <- FeaturePlot(myTCL2, 
                          reduction = "umap", 
                          pt.size = 4,
                          features = c('CCL4','CCL5'),
                          label = F, 
                          repel = T,
                          order= TRUE,
                          blend.threshold = 0.5,
                          #slot = 'counts',
                          blend = T,       
                          max.cutoff = 'q85',
                          cols = c("lightgrey", "red","green")) + theme(legend.position = 'right')
  
  pdf(paste(outdir, "Figure2_Tcell_CCL4_CCL5_blend_tcellonly_",umapr,".pdf",sep =''), width = 12, height = 4)
  print(fp_blend)
  dev.off()
  
  fp_blend <- FeaturePlot(myTCL2, 
                          reduction = "umap", 
                          pt.size = 4,
                          features = c('GZMM','GZMK'),
                          label = F, 
                          repel = T,
                          order= T,
                          blend.threshold = 0.5,
                          #slot = 'counts',
                          blend = T,       
                          max.cutoff = 'q85',
                          cols = c("lightgrey", "red","green")) + theme(legend.position = 'right')
  
  pdf(paste(outdir, "Figure2_Tcell_GZMM_GZMK_blend_tcellonly_",umapr,".pdf",sep =''),width = 12, height = 4)
  print(fp_blend)
  dev.off()
}

A = as.data.frame(table(myTCL$predicted.id, 
                        myTCL$celltype.cond)[,"TCL_iGCA"])
A2 = as.data.frame(prop.table(table(myTCL$predicted.id, 
                                    myTCL$celltype.cond))[,"TCL_iGCA"])
A$N = A2$`prop.table(table(myTCL$predicted.id, myTCL$celltype.cond))[, "TCL_iGCA"]`
A$prediction = row.names(A)
A$celltype = "TCL_iGCA"

colnames(A)  <- c("N", "perc", "prediction","celltype")
A
N <- ggplot(data=A, aes(x=celltype, y=N, fill=prediction)) +
  geom_bar(stat="identity",  color="black", position="fill") + 
  scale_fill_manual(values = col) +
  theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=12, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
        axis.title.y = element_text(face = "bold", color = "black", size = 14),
        legend.title = element_text(size = 0),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "", y = "fraction of cells")
N

A = as.data.frame(table(myTCL$pred.Tissue, 
                        myTCL$celltype)[,"TCL"])
A2 = as.data.frame(prop.table(table(myTCL$pred.Tissue, 
                                    myTCL$celltype))[,"TCL"])
A$N = A2$`prop.table(table(myTCL$pred.Tissue, myTCL$celltype))[, "TCL"]`
A$prediction = row.names(A)
A$celltype = "TCL_iGCA"

colnames(A)  <- c("N", "perc", "predictionT","celltype")

NT <- ggplot(data=A, aes(x=celltype, y=N, fill=predictionT)) +
  geom_bar(stat="identity",  color="black", position="fill") + 
  scale_fill_manual(values = colT) +
  theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=12, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
        axis.title.y = element_text(face = "bold", color = "black", size = 14),
        legend.title = element_text(size = 0),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "", y = "fraction of cells")
NT

layout <- "
AAAAAABB
AAAAAACC
"

wrap_plots(A = p, B = N, C=NT, design = layout)
ggsave(paste(outdir, "Tcell_R_mapping_",umapr,".png", sep=''), 
       wrap_plots(A = p, B = N, C=NT, design = layout), 
       width = 12, height = 6)
pdf(paste(outdir, "Tcell_mapping_",umapr,".pdf", sep=''), width = 12, height = 6)
print(wrap_plots(A = p, B = N, C=NT, design = layout))
dev.off()

# ------------ Shami --------
Shami_TCL = get(load("TCL_Shami"))
Shami_TCL
anchors_Shami <- FindTransferAnchors(reference=TCELL_new, 
                                     query=Shami_TCL, k.filter = 60,dims=1:30)
predictions_Shami <- TransferData(anchorset=anchors_Shami, 
                            refdata=TCELL_new$cell.ont, dims=1:30)
predictions_Tissue <- TransferData(anchorset=anchors_Shami, 
                                   refdata=TCELL_new$Tissue, dims=1:30)
predictions_Stim <- TransferData(anchorset=anchors_Shami, 
                                 refdata=TCELL_new$Stim, dims=1:30)
table(predictions_Shami$predicted.id)
table(predictions_Stim$predicted.id)
table(predictions_Tissue$predicted.id)

Shami_TCL <- AddMetaData(Shami_TCL, metadata=predictions_Shami)
Shami_TCL$pred.Tissue = predictions_Tissue$predicted.id

table(Shami_TCL$predicted.id)

merge()


# --------------- Integration with RPCA ---------------

# normalize and identify variable features for each dataset independently
load(paste(outdir, "TCELL_new", sep=''))
myTCL = get(load("TCL_OSR"))
myTCL$Donor = "iGCA"
myTCL$Tissue = "Testis"
myTCL$Stim = "NA"



FeaturePlot(myTCL, features = "CD8A")
myTCL
subset(myTCL, subset = CD8A>0)

allTCL.list = c(SplitObject(TCELL_new, split.by = "orig.ident"), list(myTCL))

allTCL.list <- lapply(X = allTCL.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = allTCL.list)

allTCL.list <- lapply(X = allTCL.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
options(future.globals.maxSize = 10145728000) 
TCL.anchors <- FindIntegrationAnchors(object.list = allTCL.list, 
                                         anchor.features = features, 
                                         reduction = "rpca", 
                                         k.anchor = 90)
TCL.combined <- IntegrateData(anchorset = TCL.anchors)
DefaultAssay(TCL.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
TCL.combined <- ScaleData(TCL.combined, verbose = FALSE)
TCL.combined <- RunPCA(TCL.combined, npcs = 30, verbose = FALSE)
TCL.combined <- RunUMAP(TCL.combined, reduction = "pca", dims = 1:30)
TCL.combined <- FindNeighbors(TCL.combined, reduction = "pca", dims = 1:30)
TCL.combined <- FindClusters(TCL.combined, resolution = 0.2)

p1 <- DimPlot(TCL.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(TCL.combined, reduction = "umap", 
              group.by = "seurat_annotations", label = TRUE, 
              repel = TRUE)
p1 + p2
