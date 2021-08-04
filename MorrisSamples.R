######## Morris Syndrome ###########

#### load library ####
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(repr))
suppressMessages(library(cowplot))
suppressMessages(library("viridis"))
suppressMessages(library('pals'))
suppressMessages(library(RColorBrewer))
suppressMessages(library(wesanderson))
suppressMessages(library("stringr"))
suppressMessages(library("magick"))
suppressMessages(library('pdftools'))
suppressMessages(library(rtracklayer))


Plot_sign_MG <- function(Seraut.object, signature, operator = sum, title = '') {
  x <- Seraut.object
  DefaultAssay(x) <- "RNA"
  x[["Sign_exp"]] <- apply(FetchData(object = x, 
                                     vars = signature),
                           1,
                           operator)
  FP <- FeaturePlot(x, reduction = "umap", 
                    features = 'Sign_exp', 
                    label = T, 
                    order=T,
                    repel = T, max.cutoff = 'q99',
                    cols = c("lightgrey", "red")) + labs(title = title)
  #cols = as.vector(coolwarm(5))) + 
  return(FP)
}
leydig_signature = c('CFD','DLK1','LUM','CALB2')
myoid_signature = c('ACTA2','MYH11','DES', 'MYL9')
sertoli_signature = c('FATE1','CITED1','BEX2','AMH')
macrophage_signature = c('CD14','CD74','HLA-DRA','HLA-DRB1')
endothelial_signature = c('VWF', "EGFL7",'CD34', 'PRSS23', "RBP7")
t_cell_signature =  c('GZMA','GZMK','TRAC','CD3E')

######## input ########
dataset="Morris"
counts_file = '/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/7209.counts.matrix.tsv.gz'

parameters = data.frame(Morris   = c(500, 6000,  20,5000))
row.names(parameters) <- c('min_nFeature_RNA','max_nFeature_RNA','max_percent_MT','N_hv_features')
parameters 
# constant
regresscellcycle = TRUE
nPC = 20
res = 1
thresh.use = 0.25
min.pct = 0.25
min.diff.pct = -Inf
testMG="wilcox"

matrix_sparse = as.sparse(read.table(counts_file, 
                                     sep="\t", 
                                     header = T, 
                                     quote = "", 
                                     row.names = 1))
dim(matrix_sparse)
testis <- CreateSeuratObject(counts       = matrix_sparse, 
                             project      = dataset, 
                             min.cells    = 5, 
                             min.features = 200)
testis
testis$dataset = dataset
testis$donor = dataset
###### QC plots ######
testis[["percent.mt"]] <- PercentageFeatureSet(testis, pattern = "^MT-")
VlnPlot(testis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = .1, group.by = 'dataset')
VlnPlot(testis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 1, pt.size = .1, group.by = 'donor')
plot1 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        cols = NULL, group.by = 'dataset') + NoLegend()
plot2 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                        group.by = 'dataset') + NoLegend()
plot1 / plot2

### filtering #####
testis <- subset(testis, 
                 subset = nFeature_RNA > parameters['min_nFeature_RNA',dataset] & 
                   nFeature_RNA < parameters['max_nFeature_RNA',dataset] & 
                   percent.mt < parameters['max_percent_MT',dataset])
testis
VlnPlot(testis, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
plot1 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)
plot2 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 | plot2
# normalize data
testis <- NormalizeData(testis, 
                        normalization.method = "LogNormalize", 
                        scale.factor = 10000)
###### variable feature ######
testis <- FindVariableFeatures(testis, 
                               selection.method = "vst", 
                               nfeatures = parameters['N_hv_features',dataset])
top20_2 <- head(VariableFeatures(testis), 20) # Identify the 20 most highly variableGenes
plot1 <- VariableFeaturePlot(testis)
LabelPoints(plot = plot1, points = top20_2, repel = TRUE)

####### scale data ######
# scale data and regress out %MT and nUMI
testis <- ScaleData(testis, 
                    vars.to.regress = c("percent.mt", "nFeature_RNA"),
                    features = VariableFeatures(testis))



####### pca ######
testis <- RunPCA(testis, 
                 features = VariableFeatures(object = testis))
# Examine and visualize PCA results a few different ways
print(testis[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(testis, dims = 1:3, reduction = "pca")
DimPlot(testis, 
        reduction = "pca")
#DimHeatmap(testis, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(testis, dims = 1:15, cells = 500, balanced = TRUE)

####### Cell cycle ######
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
testis <- CellCycleScoring(testis, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
testis <- RunPCA(testis, features = VariableFeatures(testis))
DimPlot(testis)
# regress out cell cycle
if (regresscellcycle) {
  testis <- ScaleData(testis, 
                      vars.to.regress = c("S.Score", "G2M.Score"), 
                      features = VariableFeatures(testis)) 
  testis <- RunPCA(testis, features = VariableFeatures(testis))
  DimPlot(testis)
}

###### nPC ######
ElbowPlot(testis, ndims = 50)

####### clustering and dimReduction ######
#clustering
res = 0.5
testis<- FindNeighbors(testis, dims = 1:nPC)
testis <- FindClusters(testis, resolution = res)
# dimentional reduction 
testis <- RunUMAP(testis, dims = 1:nPC)
#testis <- RunTSNE(testis, dims = 1:nPC)
# Dimplot    
DimPlot(testis, reduction = "umap", label = T, order = T, pt.size = 2,) + 
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=22, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2") 


# quality check
pmt = FeaturePlot(testis, reduction = "umap", pt.size = 2, features = 'percent.mt', label = T)
pRNA = FeaturePlot(testis, reduction = "umap", pt.size = 2, features = 'nFeature_RNA', label = T)
pmt + pRNA
# assign plot
FeaturePlot(testis, reduction = "umap", pt.size = 2, features = 'HMGA1', label = T)
FeaturePlot(testis, reduction = "umap", pt.size = 2, features = c('TBX3', 'HOXA3', 'DDX4'), label = T, order = T)
FeaturePlot(testis, reduction = "umap", pt.size = 2, features = c( 'HMGA1','FGFR3', 'UTF1', 'DDX4'), label = T, order = T, ncol = 3)
Seurat.object = testis
DefaultAssay(Seurat.object) <- 'RNA'
pL <- Plot_sign_MG(Seurat.object,
                   signature= leydig_signature, 
                   operator = mean, title = 'LEY')
pM <- Plot_sign_MG(Seurat.object,
                   signature= myoid_signature, 
                   operator = mean, title = 'MYD')
pS <- Plot_sign_MG(Seurat.object,
                   signature= sertoli_signature, 
                   operator = mean, title = 'SRT')
pMa <- Plot_sign_MG(Seurat.object,
                    signature= macrophage_signature, 
                    operator = mean, title = 'MCR')
pE <- Plot_sign_MG(Seurat.object,
                   signature= endothelial_signature, 
                   operator = mean, title = 'END')
pT <- Plot_sign_MG(Seurat.object,
                   signature= t_cell_signature, 
                   operator = mean, title = 'TCL')
pSTRO <- Plot_sign_MG(Seurat.object,
                      signature= c('RGS5','TPM2'), 
                      operator = mean, title = 'STRO')

(pL / pM / pS ) | (pMa / pT / pE) | pSTRO



DimPlot(testis)

####### MG #######
filename_xls = 'Morris_markergenes.xlsx'
test.use = testMG
cluster.markers = FindAllMarkers(testis, 
                                 thresh.use = thresh.use, 
                                 test.use=test.use, 
                                 min.pct=min.pct, 
                                 min.diff.pct=min.diff.pct, 
                                 only.pos=TRUE)

write.xlsx(cluster.markers,
           file= filename_xls, 
           row.names = T,
           asTable = T)

top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
hm = DoHeatmap(testis, features = top10$gene, angle = 90) + NoLegend()
hm



####### save Seurat obj ########
assign(paste('testis_',dataset, sep = ''),testis)
save(testis, file = paste('testis_',dataset, sep = ''))

