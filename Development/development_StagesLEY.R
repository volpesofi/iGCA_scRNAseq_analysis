#' developement script
#' 
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
library(patchwork)
library(stringr)
library(future)
library(viridis)
library(rcartocolor)
library("openxlsx")
library(ggpubr)

setwd("/beegfs/scratch/ric.cosr/ric.alfano/AlfanoM_904_infertilita_epigenetics/Development/")
indir = "/beegfs/scratch/ric.cosr/ric.alfano/AlfanoM_904_infertilita_epigenetics/SeuratObjects/"
outdir = "Outputs/"
dir.create(outdir, recursive = T)

######## my functions #######
source("/beegfs/scratch/ric.cosr/ric.alfano/AlfanoM_904_infertilita_epigenetics/Testis_scRNA/mySeuratfunctions.R")
source('/beegfs/scratch/ric.cosr/ric.alfano/AlfanoM_904_infertilita_epigenetics/Testis_scRNA/utility_functions.R')
source('/beegfs/scratch/ric.cosr/ric.alfano/AlfanoM_904_infertilita_epigenetics/Testis_scRNA/integration_function.R')


#----------- mycolor -----------
colors_Alf <- data.frame(celltype = c("SSCs","Differentiated S'gonia",
                                      "Early primary S'cytes",
                                      "Late primary S'cyte",
                                      "Round S'tidis",
                                      "Elong. S'tidis S1",
                                      "Elong. S'tidis S2",
                                      "Sperm",
                                      "Sperm1",
                                      "Sperm2",
                                      "Primary S'cytes",
                                      "LEY",
                                      "MYD",
                                      "SRT",
                                      "MCR",
                                      "TCL",
                                      "END",
                                      "END2",
                                      "STRO",
                                      "UND",
                                      "UND2",
                                      "CTLonly"), 
                         cols = c(brewer.pal(9, 'Oranges')[2:9],
                                  '#932F04',
                                  '#7F2704',
                                  '#FD9E54',
                                  '#E41A1C', #L
                                  '#3399FF', #M
                                  '#FFCC33', #S
                                  '#9900FF', #M
                                  '#FF33CC', #Tcell
                                  '#4DAF4A',
                                  '#4DAF9A',
                                  'cadetblue2',
                                  '#999999',
                                  "#CCCCCC",
                                  'grey'))
row.names(colors_Alf) <- colors_Alf$celltype

#----------- Load Guo list -----------
load( paste(indir, "GuoDevelop.list", sep=''))
names(GuoDevelop.list) = c("1yo", "7yo","11yo","13yo","14yo", "25yo")

#----------- Guo list integration -----------
GuoDevelop.list
dat = 'GuoDevelop_Integration'
dir = paste("/beegfs/scratch/ric.cosr/ric.alfano/AlfanoM_904_infertilita_epigenetics/Development/",
            dat,'/', sep ='')

min_nFeature_RNA = 500 
max_nFeature_RNA = 6500 
max_percent_MT = 20   
res = 0.1; nPC = 20

for (i in 1:length(GuoDevelop.list)) {
  GuoDevelop.list[[i]][["percent.mt"]] <- PercentageFeatureSet(GuoDevelop.list[[i]], 
                                                            pattern = "^MT-")
  GuoDevelop.list[[i]] <-  subset(GuoDevelop.list[[i]], 
                               subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)
  GuoDevelop.list[[i]] <- NormalizeData(GuoDevelop.list[[i]], verbose = FALSE)
  GuoDevelop.list[[i]] <- FindVariableFeatures(GuoDevelop.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

GuoDevelop.anchors <- FindIntegrationAnchors(object.list = GuoDevelop.list, dims = 1:30)
GuoDevelop.integrated <- IntegrateData(anchorset = GuoDevelop.anchors, dims = 1:30)

DefaultAssay(GuoDevelop.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering

GuoDevelop.integrated <- ScaleData(GuoDevelop.integrated, verbose = FALSE)
GuoDevelop.integrated <- RunPCA(GuoDevelop.integrated, npcs = 30, verbose = FALSE)
ElbowPlot(GuoDevelop.integrated)
GuoDevelop.integrated <- RunUMAP(GuoDevelop.integrated, reduction = "pca", dims = 1:30)
GuoDevelop.integrated <- FindNeighbors(GuoDevelop.integrated, dims = 1:nPC)
GuoDevelop.integrated <- FindClusters(GuoDevelop.integrated, resolution = res)
p1 <- DimPlot(GuoDevelop.integrated, reduction = "umap", group.by = "age")
p2 <- DimPlot(GuoDevelop.integrated, reduction = "umap", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 | p2

DefaultAssay(GuoDevelop.integrated) <- 'RNA'
pL <- Plot_sign(GuoDevelop.integrated,
                signature= c('CFD','DLK1','LUM'), 
                operator = mean, title = 'LEY')
pL
pM <- Plot_sign(GuoDevelop.integrated,
                signature= c('ACTA2','MYH11','DES'), 
                operator = mean, title = 'MYD')
pM
pS <- Plot_sign(GuoDevelop.integrated,
                signature= c('FATE1','CITED1','SOX9', 'AMH'), 
                operator = mean, title = 'SRT')
pS
pMa <- Plot_sign(GuoDevelop.integrated,
                 signature= c('CD14','CD74','HLA-DRA'), 
                 operator = mean, title = 'MCR')
pMa
pE <- Plot_sign(GuoDevelop.integrated,
                signature= c('VWF',"EGFL7", 'PRSS23'), 
                operator = mean, title = 'END')
pE
pT <- Plot_sign(GuoDevelop.integrated,
                signature= c('GZMA','CD8A','CCL5'), 
                operator = mean, title = 'TCL')
pT
pSTRO <- Plot_sign(GuoDevelop.integrated,
                   signature= c('RGS5','TPM2','IGFBP5'), 
                   operator = mean, title = 'STRO')
pSTRO
layout <- '
AABBCC##
AABBCC##
AABBCCDD
EEHHGGDD
EEHHGG##
EEHHGG##
'
png(paste(outdir, "MG_integration_",dat,"_iNOA_res",res,".png", sep=''),  width=1400, height=600)
print(wrap_plots(A = pL, B = pM, C = pS, D = pSTRO, E = pE, H = pMa, G =pT, design = layout))
dev.off()

new_cluster_id = c("SRT", "END", "MYD", "LEY", "MCR", "END")
names(new_cluster_id) <- levels(GuoDevelop.integrated)
GuoDevelop.integrated <- RenameIdents(GuoDevelop.integrated, new_cluster_id)
col = as.character(colors_Alf[levels(GuoDevelop.integrated),]$cols)

p1 = DimPlot(GuoDevelop.integrated, reduction = "umap", label = TRUE, pt.size = 1, order = TRUE) +
  scale_color_manual(values = (col)) +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=22, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")

GuoDevelop.integrated$age = factor(GuoDevelop.integrated$age, levels = c("1yo", "7yo","11yo","13yo","14yo", "25yo"))

p2 = DimPlot(GuoDevelop.integrated, reduction = "umap", group.by = "age", pt.size = 1, order = TRUE) +
  scale_color_manual(values = viridis(6)) +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=22, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")

p1 | p2
ggsave(paste(outdir, "UMAP_celltype_age.png", sep=''), p1 | p2, width = 12, height = 5)

GuoDevelop.integrated$celltype =  Idents(GuoDevelop.integrated)
GuoDevelop.integrated$donor = NULL
head(GuoDevelop.integrated[[]])

# -------------- select LEY + MYD -------------- 
Idents(GuoDevelop.integrated) <- "celltype"
DefaultAssay(GuoDevelop.integrated) <- "RNA"
LM = subset(GuoDevelop.integrated, idents = c("LEY", "MYD"))
LM.list = SplitObject(LM, split.by = "age")

# -------------- integration LEY + MYD -------------- 
for (i in 1:length(LM.list)) {
  LM.list[[i]] <- NormalizeData(LM.list[[i]], verbose = FALSE)
  LM.list[[i]] <- FindVariableFeatures(LM.list[[i]], selection.method = "vst", 
                                               nfeatures = 2000, verbose = FALSE)
  }
LM.anchors <- FindIntegrationAnchors(object.list = LM.list, dims = 1:20, 
                                     k.filter = min(200, sapply(LM.list, ncol)),
                                     k.score = 10)
LM.integrated <- IntegrateData(anchorset = LM.anchors, dims = 1:20)
DefaultAssay(LM.integrated) <- "integrated"
LM.integrated <- ScaleData(LM.integrated, verbose = FALSE)
LM.integrated <- RunPCA(LM.integrated, npcs = 30, verbose = FALSE)
ElbowPlot(LM.integrated)
LM.integrated <- RunUMAP(LM.integrated, reduction = "pca", dims = 1:20)
nPC = 20
LM.integrated <- FindNeighbors(LM.integrated, dims = 1:nPC)
LM.integrated <- FindClusters(LM.integrated, resolution = res)

p1 = DimPlot(LM.integrated, reduction = "umap", label = TRUE, pt.size = 1, order = TRUE) +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=22, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")

LM.integrated$age = factor(LM.integrated$age, levels = c("1yo", "7yo","11yo","13yo","14yo", "25yo"))

p2 = DimPlot(LM.integrated, reduction = "umap", group.by = "age", pt.size = 1, order = TRUE) +
  scale_color_manual(values = viridis(6)) +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=22, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")

p1
p2

DimPlot(LM.integrated, split.by = "celltype", group.by = "age", reduction = "pca")
LM.integrated[[]]
DefaultAssay(LM.integrated) <- 'RNA'
pL <- Plot_sign(LM.integrated,
                signature= c('CFD','DLK1','LUM'), 
                operator = mean, title = 'LEY')
pL
pM <- Plot_sign(LM.integrated,
                signature= c('ACTA2','MYH11','DES'), 
                operator = mean, title = 'MYD')
pM


#dev_int <- get(load( paste(indir, "int_dev.rds", sep='')))
#DefaultAssay(dev_int) <- "RNA"
#all.dev.list = SplitObject(dev_int, split.by = "source")

load(paste(indir, "cleanAdultIntegration", sep=''))
DefaultAssay(cleanAdultIntegration) <- "RNA"
Idents(cleanAdultIntegration) <- "celltype"

adult.list = SplitObject(cleanAdultIntegration, split.by = "source")

subset_iNOA3_LM = subset(adult.list$`iNOA donor 3`, idents = c("LEY", "MYD"))
subset_iNOA2_LM = subset(adult.list$`iNOA donor 2`, idents = c("LEY", "MYD"))
subset_iNOA1_LM = subset(adult.list$`iNOA donor 1`, idents = c("LEY", "MYD"))
subset_MW2_LM = subset(adult.list$GSE124263_adult_Donor2, idents = c("LEY", "MYD"))
subset_GD3_LM = subset(adult.list$GSE112013_adult_Donor3, idents = c("LEY", "MYD"))
subset_GD2_LM = subset(adult.list$GSE112013_adult_Donor2, idents = c("LEY", "MYD"))

DefaultAssay(LM.integrated) <- "integrated"

anchors <- FindTransferAnchors(reference=LM.integrated, query=subset_iNOA3_LM, dims=1:30)
predictionsAge <- TransferData(anchorset=anchors, 
                               refdata=LM.integrated$age, dims=1:30)
table(predictionsAge$predicted.id)

anchors <- FindTransferAnchors(reference=LM.integrated, query=subset_iNOA2_LM, dims=1:30)
predictionsAge <- TransferData(anchorset=anchors, 
                               refdata=LM.integrated$age, dims=1:30)
table(predictionsAge$predicted.id)

anchors <- FindTransferAnchors(reference=LM.integrated, query=subset_iNOA1_LM, dims=1:30)
predictionsAge <- TransferData(anchorset=anchors, 
                               refdata=LM.integrated$age, dims=1:30)
table(predictionsAge$predicted.id)
  

anchors <- FindTransferAnchors(reference=LM.integrated, query=subset_MW2_LM, dims=1:30)
predictionsAge <- TransferData(anchorset=anchors, 
                               refdata=LM.integrated$age, dims=1:30)
table(predictionsAge$predicted.id)

anchors <- FindTransferAnchors(reference=LM.integrated, query=subset_GD3_LM, dims=1:30, k.filter = 150)
predictionsAge <- TransferData(anchorset=anchors, 
                               refdata=LM.integrated$age, dims=1:30)
table(predictionsAge$predicted.id)

anchors <- FindTransferAnchors(reference=LM.integrated, query=subset_GD2_LM, dims=1:30, k.filter = 150)
predictionsAge <- TransferData(anchorset=anchors, 
                               refdata=LM.integrated$age, dims=1:30)
table(predictionsAge$predicted.id)


# ----------- merge Guo --------------
GuoDevelop.list
GuoMerge = merge(x = GuoDevelop.list[[1]], 
                 y = GuoDevelop.list[2:length(names(GuoDevelop.list))],
                 add.cell.ids = names(GuoDevelop.list), 
                 project = "Guodev")
#Normalize data
GuoMerge <- NormalizeData(GuoMerge, normalization.method = "LogNormalize", 
                       scale.factor = 10000)
# Find Variable Features
GuoMerge <- FindVariableFeatures(GuoMerge, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GuoMerge)
options(future.globals.maxSize = 3145728000) 
plan(strategy = "multicore", workers = 6)
GuoMerge <- ScaleData(GuoMerge, features = all.genes)
# PCA
GuoMerge <- RunPCA(GuoMerge, features = VariableFeatures(object = GuoMerge), npcs = 50)
ElbowPlot(GuoMerge, ndims = 40)
# cluster the cells 
nPC = 20; res = 0.2
GuoMerge <- FindNeighbors(GuoMerge, dims = 1:nPC)
GuoMerge <- FindClusters(GuoMerge, resolution = res)
# UMAP
GuoMerge <- RunUMAP(GuoMerge, dims = 1:nPC)
GuoMerge <- RunTSNE(GuoMerge, dims = 1:nPC)

DimPlot(GuoMerge, reduction = "umap", label = T)

FeaturePlot(GuoMerge, "DLK1", reduction = "umap", label = T, order = T)
FeaturePlot(GuoMerge, "MYH11", reduction = "umap", label = T, order = T)

DimPlot(GuoMerge, reduction = "tsne", group.by = "age")
DimPlot(GuoMerge, reduction = "umap", label = T)

new_cluster_id = c("END", "SRT", "END", "SRT", "MYD", "LEY", "MCR", "LEY", "STRO", "SRT", "STRO", "SRT")
names(new_cluster_id) <- levels(GuoMerge)
GuoMerge <- RenameIdents(GuoMerge, new_cluster_id)

col = as.character(colors_Alf[levels(GuoMerge),]$cols)
p1 = DimPlot(GuoMerge, reduction = "umap", label = TRUE, pt.size = 1, order = TRUE) +
  scale_color_manual(values = (col)) +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2") + NoLegend()
p1
ggsave(paste(outdir, "UMAP_merge_celltype.png", sep=''), p1, width = 7, height = 6)
pdf(paste(outdir, "UMAP_merge_celltype.pdf", sep=''),  width = 7, height = 6)
print(p1)
dev.off()
col2 = c("#CCCCCC" ,   "#CCCCCC" ,   "#3399FF",    "#E41A1C" ,   "#CCCCCC"  ,  "#CCCCCC")
p1 = DimPlot(GuoMerge, reduction = "umap", label = TRUE, pt.size = 1, order = TRUE) +
  scale_color_manual(values = (col2)) +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2") + NoLegend()
ggsave(paste(outdir, "UMAP_merge_celltype_2.png", sep=''), p1, width = 6, height = 6)
pdf(paste(outdir, "UMAP_merge_celltype_2.pdf", sep=''),  width = 6, height = 6)
print(p1)
dev.off()
GuoMerge$celltype = Idents(GuoMerge)

# ------------- subset MYD LEY merge --------------
LM_MERGE = subset(GuoMerge, idents = c("MYD","LEY"))
#DimPlot(LM_MERGE , reduction = "umap", group.by = "celltype")
old = as.character(LM_MERGE$seurat_clusters)
old[old == "4"] = "MYD"
old[old %in% c("5", "7")] = "LEY"
table(old)
LM_MERGE$celltype = as.factor(old)

LM_MERGE
LM_MERGE <- NormalizeData(LM_MERGE, normalization.method = "LogNormalize", 
                          scale.factor = 10000)
LM_MERGE <- FindVariableFeatures(LM_MERGE, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(LM_MERGE)
LM_MERGE <- ScaleData(LM_MERGE, features = all.genes)
LM_MERGE <- RunPCA(LM_MERGE, features = VariableFeatures(object = LM_MERGE), npcs = 50)
ElbowPlot(LM_MERGE, ndims = 50)
# cluster the cells 
nPC = 30; res = 0.2
LM_MERGE <- RunUMAP(LM_MERGE, dims = 1:nPC)
LM_MERGE <- RunTSNE(LM_MERGE, dims = 1:nPC)
LM_MERGE <- FindNeighbors(LM_MERGE, dims = 1:nPC)
LM_MERGE <- FindClusters(LM_MERGE, resolution = res)
LM_MERGE$age = factor(LM_MERGE$age, levels = c("1yo", "7yo","11yo","13yo","14yo", "25yo"))

p1 = DimPlot(LM_MERGE, reduction = "umap", group.by = "celltype", pt.size = 1, order = TRUE) +
  scale_color_manual(values = c("#3399FF", "#E41A1C")) +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2") 
p1

p2 = DimPlot(LM_MERGE, reduction = "umap", group.by = "age", pt.size = 1, order = TRUE) +
  scale_color_manual(values = viridis(length(levels(LM_MERGE$age)), option = "C")) +
  #scale_color_manual(values = carto_pal(length(levels(LM_MERGE$age)), "TealRose")) +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2") 
p2

p3 = DimPlot(LM_MERGE, reduction = "umap", pt.size = 1, order = TRUE) +
  #scale_color_manual(values = viridis(length(levels(LM_MERGE$age)), option = "C")) +
  #scale_color_manual(values = carto_pal(length(levels(LM_MERGE$age)), "TealRose")) +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2") 
p3

p1 | p2
ggsave(paste(outdir, "UMAP_LM_MERGE_celltype.png", sep=''), p1 | p2, width = 12, height = 6)
pdf(paste(outdir, "UMAP_LM_MERGE_age.pdf", sep=''),  width = 12, height = 6)
print(p1 | p2)
dev.off()

# -------- LEY Stages ----------
new.label = c("MYD", "Stage_a", "Stage_b", "Stage_c")
names(new.label) <- levels(LM_MERGE)
LM_MERGE <- RenameIdents(LM_MERGE, new.label)
DimPlot(LM_MERGE)
cm_Stage_c = FindMarkers(LM_MERGE, ident.1 = "Stage_c", ident.2 = c("Stage_b", "Stage_a"), only.pos = T,
                              test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
cm_Stage_c = cm_Stage_c[order(cm_Stage_c$avg_logFC, decreasing = T),]
head(cm_Stage_c, n = 20)

signature_Stage_c = row.names(cm_Stage_c[which(cm_Stage_c$p_val_adj < 0.05 & cm_Stage_c$avg_logFC > 0.5),])[1:23]
signature_Stage_c = signature_Stage_c[-c(1:3)]
signature_Stage_c 

cm_Stage_a = FindMarkers(LM_MERGE, ident.1 = "Stage_a", ident.2 = c("Stage_b", "Stage_c"), only.pos = T,
                         test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
cm_Stage_a = cm_Stage_a[order(cm_Stage_a$avg_logFC, decreasing = T),]
head(cm_Stage_a, n = 20)
signature_Stage_a = row.names(cm_Stage_a[which(cm_Stage_a$p_val_adj < 0.05 & cm_Stage_a$avg_logFC > 0.5),])[1:20]
signature_Stage_a

cm_Stage_b = FindMarkers(LM_MERGE, ident.1 = "Stage_b", ident.2 = c("Stage_a", "Stage_c"), only.pos = T,
                         test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
cm_Stage_b = cm_Stage_b[order(cm_Stage_b$avg_logFC, decreasing = T),]
head(cm_Stage_b, n = 20)
signature_Stage_b = row.names(cm_Stage_b[which(cm_Stage_b$p_val_adj < 0.05 & cm_Stage_b$avg_logFC > 0.5),])[1:20]
length(signature_Stage_b)
signature_Stage_b

Stages_list = list(Stage_A =  cm_Stage_a,
                   Stage_B = cm_Stage_b,
                   Stage_C = cm_Stage_c) 
write.xlsx(Stages_list, file = paste(outdir, "StagesLEY_MarkerGenes.xlsx", sep =''), row.names = T, asTable = T)
MG = c(row.names(cm_Stage_a[which(cm_Stage_a$p_val_adj < 0.05 & cm_Stage_a$avg_logFC > 0.25),]),
       row.names(cm_Stage_b[which(cm_Stage_b$p_val_adj < 0.05 & cm_Stage_b$avg_logFC > 0.25),]),
       row.names(cm_Stage_c[which(cm_Stage_c$p_val_adj < 0.05 & cm_Stage_c$avg_logFC > 0.25),]))
MGl = list( A = row.names(cm_Stage_a[which(cm_Stage_a$p_val_adj < 0.05 & cm_Stage_a$avg_logFC > 0.25),]),
       B = row.names(cm_Stage_b[which(cm_Stage_b$p_val_adj < 0.05 & cm_Stage_b$avg_logFC > 0.25),]),
       C = row.names(cm_Stage_c[which(cm_Stage_c$p_val_adj < 0.05 & cm_Stage_c$avg_logFC > 0.25),]))
lengths(MGl)

hm = DoHeatmap(subset(LM_MERGE, idents = c("Stage_a", "Stage_b","Stage_c")),
               features = MG, group.colors = viridis(3, option = "C"),
               disp.min = -1.5, disp.max = 1.5, draw.lines = TRUE) + 
  theme(axis.text.y = element_text(size = 0)) +
  scale_fill_carto_c(name = "Expression: ",
                     type = "diverging", palette = "Earth", direction = -1) 
hm
ggsave(filename = paste(outdir, "Heatmap_StagesLEY.png", sep =''), hm, height = 10, width = 9)
pdf(file = paste(outdir, "Heatmap_StagesLEY.pdf", sep =''), height = 10, width = 9)
print(hm)
dev.off()

signature_Stage_a
VlnPlot(LM_MERGE, c("JUNB","HES1" ,"WFDC1","EGR1", "DLK1"), 
        idents = c("Stage_a", "Stage_b","Stage_c"), pt.size = 0, ncol = 5, cols = viridis(3, option = "C"))  
ggsave(paste(outdir, "Signature_StageA_violin.png", sep=''), width = 12, height = 3)

signature_Stage_b
VlnPlot(LM_MERGE, c("IGFBP3", "CAVIN2","C7", "CLU", "COL6A3"), 
        idents = c("Stage_a", "Stage_b","Stage_c"), pt.size = 0, ncol = 5, cols = viridis(3, option = "C"))  
ggsave(paste(outdir, "Signature_StageB_violin.png", sep=''), width = 12, height = 3)

signature_Stage_c
VlnPlot(LM_MERGE, c("CFD", "APOD","PI16", "IFI27", "GPX3"), 
        idents = c("Stage_a", "Stage_b","Stage_c"), pt.size = 0, ncol = 5, cols = viridis(3, option = "C"))  
ggsave(paste(outdir, "Signature_StageC_violin.png", sep=''), width = 12, height = 3)


Plot_sign(Seraut.object = LM_MERGE, signature = signature_Stage_a, title = " Signature Stage a", operator = sum) +
  theme(plot.title = element_text(color="black", size=18, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=0, face="italic"))
Plot_sign(Seraut.object = LM_MERGE, signature = signature_Stage_b, title = " Signature Stage b", operator = sum) +
  theme(plot.title = element_text(color="black", size=18, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=0, face="italic"))
Plot_sign(Seraut.object = LM_MERGE, signature = signature_Stage_c, title = " Signature Stage c", operator = sum) +
  theme(plot.title = element_text(color="black", size=18, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=0, face="italic"))

signature_Stage_a
signature_Stage_b
signature_Stage_c

# -------- Stages Signature genes ----------
Seurat.object = cleanAdultIntegration
Seurat.object$condition = plyr::revalue(Seurat.object$condition, c("iNOA"="iGCA", "CTL"="CTL"))
Seurat.object$condition = factor(Seurat.object$condition, levels = c("iGCA", "CTL"))
VlnPlot(Seurat.object , idents = "LEY", assay = "RNA", features = c("HES1","EGR1", "JUNB", "FOS","DLK1","MEG3"),
        group.by = "condition", split.plot = F, pt.size = 0, cols = c('orange','dodgerblue2'), ncol = 6)
ggsave(paste(outdir, "Signature_StageA_iNOACTL_violin.png", sep=''), width = 18, height = 3)
pdf(paste(outdir, "Signature_StageA_iNOACTL_violin.pdf", sep=''), width = 18, height = 3)
print(VlnPlot(Seurat.object , idents = "LEY", assay = "RNA", features = c("HES1","EGR1", "JUNB", "FOS","DLK1","MEG3"),
              group.by = "condition", split.plot = F, pt.size = 0, cols = c('orange','dodgerblue2'), ncol = 6))
dev.off()

VlnPlot(Seurat.object , idents = "LEY", assay = "RNA", features = c("IGFBP3", "S100A6","C7","PLPP3", "COL6A3", 
                                                                    "SERPINF1"),
        group.by = "condition", split.plot = F, pt.size = 0, cols = c('orange','dodgerblue2'), ncol = 6)
ggsave(paste(outdir, "Signature_StageB_iNOACTL_violin.png", sep=''), width = 18, height = 3)
pdf(paste(outdir, "Signature_StageB_iNOACTL_violin.pdf", sep=''), width = 18, height = 3)
print(VlnPlot(Seurat.object , idents = "LEY", assay = "RNA", features = c("IGFBP3", "S100A6","C7","PLPP3", "COL6A3", 
                                                                          "SERPINF1"),
              group.by = "condition", split.plot = F, pt.size = 0, cols = c('orange','dodgerblue2'), ncol = 6))
dev.off()


VlnPlot(Seurat.object , idents = "LEY", assay = "RNA", features = c("PTGDS", "TIMP3", "IFI27", "CCN5", "GPX3", "MT1X"),
        group.by = "condition", split.plot = F, pt.size = 0, cols = c('orange','dodgerblue2'), ncol = 6)
ggsave(paste(outdir, "Signature_StageC_iNOACTL_violin.png", sep=''), width = 18, height = 3)
pdf(paste(outdir, "Signature_StageC_iNOACTL_violin.pdf", sep=''), width = 18, height = 3)
print(VlnPlot(Seurat.object , idents = "LEY", assay = "RNA", features = c("PTGDS", "TIMP3", "IFI27", "CCN5", "GPX3", "MT1X"),
              group.by = "condition", split.plot = F, pt.size = 0, cols = c('orange','dodgerblue2'), ncol = 6))
dev.off()

VlnPlot(cleanAdultIntegration, idents = "LEY", assay = "RNA", features = "IGFBP3", 
        group.by = "condition", split.plot = F, pt.size = 0)
VlnPlot(cleanAdultIntegration, idents = "LEY", assay = "RNA", features = "C7", 
        group.by = "condition", split.plot = F, pt.size = 0)
VlnPlot(cleanAdultIntegration, idents = "LEY", assay = "RNA", features = "JUNB", 
        group.by = "condition", split.plot = F, pt.size = 0)
VlnPlot(cleanAdultIntegration, idents = "LEY", assay = "RNA", features = "DLK1", 
        group.by = "condition", split.plot = F, pt.size = 0)

# -------- Stages Signature VlnPlot ----------
DefaultAssay(cleanAdultIntegration) <- "RNA"

x <- cleanAdultIntegration
signature = signature_Stage_a
length(signature_Stage_a)
head(cm_Stage_a, n = 20)
DefaultAssay(x) <- "RNA"
x[["SignatureA"]] <- apply(FetchData(object = x, 
                           vars = signature),
                 1,
                 sum)
VlnPlot(x,
        features = "SignatureA",
        group.by = 'condition', 
        pt.size = 0,
        idents = "LEY") +
  ggtitle("Stage a") +
  ylab("Sum of Expression") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=18, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) 

ct = "LEY"
so = subset(x, idents = ct) 
so[["SignatureA"]] <- FetchData(object = so, vars = "SignatureA")
df= so[[]]
df$condition = plyr::revalue(df$condition, c("iNOA"="iGCA_LEY", "CTL"="CTL_LEY"))
df$condition = factor(df$condition, levels = c("iGCA_LEY", "CTL_LEY"))
my_comparisons <- list(c("iGCA_LEY", "CTL_LEY"))

v <- ggviolin(df, x = "condition", y = "SignatureA", 
              fill = "condition", color = "black", width = 1,
              #add = c('jitter'), add.params = list(size =.5, shape = 1),
              #              add.params = list(size = .5, shape = 20),
              alpha = 0.8,
              #facet.by = 'celltype', 
              short.panel.labs = TRUE,
              palette = c('orange','dodgerblue2'), 
              error.plot = "crossbar", draw_quantiles = c(0.5), trim = T,
              #add = c('mean','mean_se'),
              title = "Signature Stage A", 
              panel.labs.background = list(color = "white", fill = "white", size = 0.5),
              panel.labs.font = list(color = "black", face = "bold", size = 16))
v = ggpar(v, ylab = 'Expression level', xlab = '', ylim = c(0, 60))
v = v + stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6) +
  theme(plot.title = element_text(color="black", size=18, face="bold.italic", hjust = 0.5),
        plot.title.position =  'panel',
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 14, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=14),
        axis.title.y = element_text(face = "bold", color = "black", size = 16),
        legend.text = element_text(color = "black", size = 12),
        legend.position = 'right',
        panel.spacing = unit('-0.1', "lines"))
v
v_signatureA = v

signature = signature_Stage_b
signature
DefaultAssay(x) <- "RNA"
x[["SignatureB"]] <- apply(FetchData(object = x, 
                                    vars = signature),
                          1,
                          sum)
VlnPlot(x,
        features = "SignatureB",
        group.by = 'condition', 
        pt.size = 0,
        idents = "LEY") +
  ggtitle("Stage b") +
  ylab("Sum of Expression") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=18, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) 

so = subset(x, idents = ct)
so[["SignatureB"]] <- FetchData(object = so, vars = "SignatureB")
df= so[[]]
df$condition = plyr::revalue(df$condition, c("iNOA"="iGCA_LEY", "CTL"="CTL_LEY"))
df$condition = factor(df$condition, levels = c("iGCA_LEY", "CTL_LEY"))
my_comparisons <- list(c("iGCA_LEY", "CTL_LEY"))

v <- ggviolin(df, x = "condition", y = "SignatureB", 
              fill = "condition", color = "black", width = 1,
              #add = c('jitter'), add.params = list(size =.5, shape = 1),
              #              add.params = list(size = .5, shape = 20),
              alpha = 0.8,
              #facet.by = 'celltype', 
              short.panel.labs = TRUE,
              palette = c('orange','dodgerblue2'), 
              error.plot = "crossbar", draw_quantiles = c(0.5), trim = T,
              #add = c('mean','mean_se'),
              title = "Signature Stage B", 
              panel.labs.background = list(color = "white", fill = "white", size = 0.5),
              panel.labs.font = list(color = "black", face = "bold", size = 16))
v = ggpar(v, ylab = 'Expression level', xlab = '', ylim = c(0, 60))
v = v + stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6) +
  theme(plot.title = element_text(color="black", size=18, face="bold.italic", hjust = 0.5),
        plot.title.position =  'panel',
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 14, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=14),
        axis.title.y = element_text(face = "bold", color = "black", size = 16),
        legend.text = element_text(color = "black", size = 12),
        legend.position = 'right',
        panel.spacing = unit('-0.1', "lines"))
v
v_signatureB = v

signature = signature_Stage_c
DefaultAssay(x) <- "RNA"
x[["SignatureC"]] <- apply(FetchData(object = x, 
                                    vars = signature),
                          1,
                         sum)
VlnPlot(x,
        features = "SignatureC",
        group.by = 'condition', 
        pt.size = 0,
        idents = "LEY") +
  ggtitle("Stage c") +
  ylab("Sum of Expression") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=18, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) 

so = subset(x, idents = ct) 
so[["SignatureC"]] <- FetchData(object = so, vars = "SignatureC")
df= so[[]]
df$condition = plyr::revalue(df$condition, c("iNOA"="iGCA_LEY", "CTL"="CTL_LEY"))
df$condition = factor(df$condition, levels = c("iGCA_LEY", "CTL_LEY"))
my_comparisons <- list(c("iGCA_LEY", "CTL_LEY"))

v <- ggviolin(df, x = "condition", y = "SignatureC", 
              fill = "condition", color = "black", width = 1,
              #add = c('jitter'), add.params = list(size =.5, shape = 1),
              #              add.params = list(size = .5, shape = 20),
              alpha = 0.8,
              #facet.by = 'celltype', 
              short.panel.labs = TRUE,
              palette = c('orange','dodgerblue2'), 
              error.plot = "crossbar", draw_quantiles = c(0.5), trim = T,
              #add = c('mean','mean_se'),
              title = "Signature Stage C", 
              panel.labs.background = list(color = "white", fill = "white", size = 0.5),
              panel.labs.font = list(color = "black", face = "bold", size = 16))
v = ggpar(v, ylab = 'Expression level', xlab = '', ylim = c(0, 60))
v = v + stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6) +
  theme(plot.title = element_text(color="black", size=18, face="bold.italic", hjust = 0.5),
        plot.title.position =  'panel',
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 14, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=14),
        axis.title.y = element_text(face = "bold", color = "black", size = 16),
        legend.text = element_text(color = "black", size = 12),
        legend.position = 'right',
        panel.spacing = unit('-0.1', "lines"))
v
v_signatureC = v
(v_signatureA + v_signatureB + v_signatureC) + plot_layout(guides = "collect") 
ggsave(paste(outdir, "Signature_Stages_iGCA_CTL.png", sep =''), width = 14, height = 5)
pdf(paste(outdir, "Signature_Stages_iGCA_CTL.pdf", sep =''), width = 14, height = 5)
print((v_signatureA + v_signatureB + v_signatureC) + plot_layout(guides = "collect"))
dev.off()



VlnPlot(x, y.max = 50,
        features = c("SignatureA","SignatureB","SignatureC"),
        group.by = 'condition', 
        pt.size = 0, 
        idents = "LEY")

  #geom_boxplot(width = 0.1, outlier.size=1) +
  #scale_fill_manual(values = viridis(8))
  #scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))
ggsave(paste(dpaper, "VlnPlot_",ls,".png", sep = ''), width = 12, height = 4)

anchors <- FindTransferAnchors(reference=LM_MERGE_prova, query=subset_GD3_LM, dims=1:30, k.filter = 150)
predictionsAge <- TransferData(anchorset=anchors, 
                               refdata=LM_MERGE_prova$age, dims=1:30)
table(predictionsAge$predicted.id)

anchors <- FindTransferAnchors(reference=LM_MERGE_prova, query=subset_iNOA3_LM, dims=1:30)
predictionsAge <- TransferData(anchorset=anchors, 
                               refdata=LM_MERGE_prova$age, dims=1:30)
table(predictionsAge$predicted.id)

anchors <- FindTransferAnchors(reference=LM_MERGE_prova, query=subset_iNOA2_LM, dims=1:30)
predictionsAge <- TransferData(anchorset=anchors, 
                               refdata=LM_MERGE_prova$age, dims=1:30)
table(predictionsAge$predicted.id)

anchors <- FindTransferAnchors(reference=LM_MERGE_prova, query=subset_iNOA1_LM, dims=1:30)
predictionsAge <- TransferData(anchorset=anchors, 
                               refdata=LM_MERGE_prova$age, dims=1:30)
table(predictionsAge$predicted.id)

anchors <- FindTransferAnchors(reference=LM_MERGE_prova, query=subset_MW2_LM, dims=1:30)
predictionsAge <- TransferData(anchorset=anchors, 
                               refdata=LM_MERGE_prova$age, dims=1:30)
table(predictionsAge$predicted.id)

anchors <- FindTransferAnchors(reference=LM_MERGE_prova, query=, dims=1:30)
predictionsAge <- TransferData(anchorset=anchors, 
                               refdata=LM_MERGE_prova$age, dims=1:30)
table(predictionsAge$predicted.id)

save.image(file = paste(outdir,"LeyDev.RData", sep=''))
