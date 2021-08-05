# libraries
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

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
setwd('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/')

##### my functions ####
source("mySeuratfunctions.R")
source('./utility_functions.R')
source('./integration_function.R')

leydig_signature = c('CFD','DLK1','LUM','CALB2')
myoid_signature = c('ACTA2','MYH11','DES', 'MYL9')
sertoli_signature = c('FATE1','CITED1','BEX2','AMH')
macrophage_signature = c('CD14','CD74','HLA-DRA','HLA-DRB1')
endothelial_signature = c('VWF', "EGFL7",'CD34', 'PRSS23', "RBP7")
t_cell_signature =  c('GZMA','GZMK','TRAC','CD3E')

somatic_clu_list = list()

##### dataset info and parameters #######
# Guo et al - infant testis
# import data from litterature
dataset_l = c("GSE120506", "GSE142585", "GSE109037", "GSE124263", 
              "GSE134144","GSE134144_t", "GSE125372", "GSE106487")
dataset = dataset_l[8]
is10x_l = list()
is10x_l[["GSE120506"]]   = FALSE
is10x_l[["GSE142585"]]   = FALSE
is10x_l[["GSE109037"]]   = TRUE
is10x_l[["GSE124263"]]   = TRUE
is10x_l[["GSE134144"]]   = FALSE
is10x_l[["GSE134144_t"]] = FALSE
is10x_l[["GSE125372"]]   = FALSE
is10x_l[["GSE106487"]]   = FALSE

counts_file_list = list()
counts_file_list[["GSE120506"]] = 'datasets/GSE120506_infant_combined_UMI.txt'
counts_file_list[["GSE142585"]] = 'datasets/GSE142585_MergedHumanTestis4_CountMatrix.txt'
counts_file_list[["GSE109037"]] = 'datasets/GSE109037_RAW/'
counts_file_list[["GSE124263"]] = 'datasets/GSE124263_RAW/'
counts_file_list[["GSE134144"]] = 'datasets/GSE134144_Pubertal_combined_UMI.txt'
counts_file_list[["GSE134144_t"]] = 'datasets/GSE134144_Transgender_combined_UMI.txt'
counts_file_list[["GSE125372"]] = c('datasets/GSE125372_RAW/GSM3572760_TESE1_1_gene_expression.tsv.gz',
                                    'datasets/GSE125372_RAW/GSM3572761_TESE1_2_gene_expression.tsv.gz',
                                    'datasets/GSE125372_RAW/GSM3572762_TESE2_1_gene_expression.tsv.gz',
                                    'datasets/GSE125372_RAW/GSM3572763_TESE2_2_gene_expression.tsv.gz')
counts_file_list[["GSE106487"]] =  'datasets/GSE106487_RAW/'

# parameters
# sample specific
parameters = data.frame(GSE120506   = c(500, 4000,  10,5000), 
                        GSE142585   = c(500, 6000,  10, 2000),
                        GSE109037   = c(500, 10000, 20, 2000),
                        GSE124263   = c(500, 6000,  20, 2000),
                        GSE134144   = c(500, 6000,  20, 2000),
                        GSE134144_t = c(500, 6000,  20, 2000),
                        GSE125372   = c(500, 6000,  20, 2000),
                        GSE106487   = c(500, 15000, 20, 2000))
row.names(parameters) <- c('min_nFeature_RNA','max_nFeature_RNA','max_percent_MT','N_hv_features')
# constant
regresscellcycle = T
nPC = 20
res = .5
thresh.use = 0.25
min.pct = 0.25
min.diff.pct = -Inf
testMG="wilcox"
filename_xls <- paste('marker_genes_',dataset ,'_res',res,'_nPC',nPC,'.xlsx',sep='')

is10x = is10x_l[[dataset]]

if (dataset != 'GSE125372') {
  if (is10x) {
    if (dataset == 'GSE109037') {
      matrix_sparse <- Read10X(data.dir = c(paste(counts_file_list[[dataset]],'AdultHuman_17-3', sep=''),
                                            paste(counts_file_list[[dataset]],'AdultHuman_17-4', sep=''),
                                            paste(counts_file_list[[dataset]],'AdultHuman_17-5', sep='')))
    } else if (dataset == 'GSE124263'){
      matrix_sparse <- Read10X(data.dir = c(paste(counts_file_list[[dataset]],'neonatal_day2', sep=''),
                                            paste(counts_file_list[[dataset]],'neonatal_day7', sep=''),
                                            paste(counts_file_list[[dataset]],'adult_Donor1', sep=''),
                                            paste(counts_file_list[[dataset]],'adult_Donor2', sep='')))
    }
  } else {
    if (dataset == 'GSE106487') {
      # Smart-seq2
      fl = list.files(counts_file_list[[dataset]], pattern=glob2rx("*_gene_expression_TPM.txt.gz"))
      count_list = list()
      for (file in fl) {
        count_list[[file]] = read.table(paste(counts_file_list[[dataset]],file, sep =''), 
                                        header = 1, check.names = F)
      }
      counts = Reduce(function(...) merge(..., by = c('Gene'), sort = F), count_list)
      counts_m = counts[,2:dim(counts)[2]]
      row.names(counts_m) = counts$Gene
      matrix_sparse = as.sparse(counts_m)
      rm(list = c("count_list","counts", "counts_m"))
    } else {
    matrix_sparse = as.sparse(read.table(counts_file_list[[dataset]], 
                                         sep="\t", 
                                         header = T, 
                                         quote = "", 
                                         row.names = 1))
    }
  }
  
  
  testis <- CreateSeuratObject(counts       = matrix_sparse, 
                               project      = dataset, 
                               min.cells    = 5, 
                               min.features = 200)
} else {
  for (i in 1:length(counts_file_list[[dataset]])) {
    print(counts_file_list[[dataset]][i])
    m = as.sparse(read.table(counts_file_list[[dataset]][i], 
                             sep="\t", 
                             header = T, 
                             quote = "", 
                             row.names = 1))
    t <- CreateSeuratObject(counts = m, 
                            project = dataset, 
                            min.cells = 5, 
                            min.features = 200)
    t$donor = str_sub(counts_file_list[[dataset]][i],35,41)
    assign(paste('t',i,sep='_'),t)
  }
  test = merge(x = t_1, y = t_2, add.cell.ids = c('TESE1_1', 'TESE1_2'))
  testi = merge(x =test, y = t_3, add.cell.ids = c('', 'TESE2_1'))
  testis = merge(x =testi,y = t_4, add.cell.ids = c('', 'TESE2_2'))
  rm(list = c("t_1","t_2", "t_3","t_4"))
}

testis
testis$dataset = dataset

if (dataset == 'GSE124263'){
  barcodes_label = data.frame(barcode = row.names(testis[[]]), l = str_sub(row.names(testis[[]]),1,2))
  barcodes_label$donor = "neonatal_day2"
  barcodes_label[grep('2_',barcodes_label$l),]$donor = "neonatal_day7"
  barcodes_label[grep('3_',barcodes_label$l),]$donor = "adult_Donor1"
  barcodes_label[grep('4_',barcodes_label$l),]$donor = "adult_Donor2"
  barcodes_label$dev = "adult"
  barcodes_label[grep('neonatal',barcodes_label$donor),]$dev = "neonatal"
  testis$donor = barcodes_label$donor
  testis$dev = barcodes_label$dev
  
}

if (dataset == 'GSE134144'){
  barcodes_label = data.frame(barcode = row.names(testis[[]]), 
                              donor = str_sub(row.names(testis[[]]),1,3))
  testis$donor =  barcodes_label$donor
}

if (dataset == 'GSE134144_t'){
  barcodes_label = data.frame(barcode = row.names(testis[[]]), 
                              donor = str_sub(row.names(testis[[]]),1,6))
  testis$donor =  barcodes_label$donor
}

if (dataset == 'GSE106487'){
  barcodes_label = data.frame(barcode = row.names(testis[[]]), 
                              donor = str_sub(row.names(testis[[]]),1,4))
  testis$donor = barcodes_label$donor
}


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
testis
###### filter ######
if (dataset != "GSE106487") {
  testis <- subset(testis, 
                   subset = nFeature_RNA > parameters['min_nFeature_RNA',dataset] & 
                     nFeature_RNA < parameters['max_nFeature_RNA',dataset] & 
                     percent.mt < parameters['max_percent_MT',dataset])
  testis
}

VlnPlot(testis, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
plot1 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)
plot2 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 | plot2

###### normalise data ######
testis_bk = testis
#testis = testis_bk

if (dataset == "GSE106487") {
  # log of TPM
  prior.count = 1
  testis <- SetAssayData(object = testis, 
                         slot = "data", 
                         new.data = as.sparse(log(testis@assays$RNA@counts+prior.count)))
} else {
  # normalize data
  testis <- NormalizeData(testis, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)
}
###### variable feature ######
testis <- FindVariableFeatures(testis, 
                               selection.method = "vst", 
                               nfeatures = parameters['N_hv_features',dataset])
top20_2 <- head(VariableFeatures(testis), 20) # Identify the 20 most highly variableGenes
plot1 <- VariableFeaturePlot(testis)
LabelPoints(plot = plot1, points = top20_2, repel = TRUE)

####### scale data ######
# scale data and regress out %MT and nUMI

#all.genes <- rownames(testis)
if (dataset == "GSE106487") {
  testis <- ScaleData(testis, 
                      vars.to.regress = c("nFeature_RNA"),
                      features = VariableFeatures(testis))
} else {
  testis <- ScaleData(testis, 
                      vars.to.regress = c("percent.mt", "nFeature_RNA"),
                      features = VariableFeatures(testis))
}

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
testis<- FindNeighbors(testis, dims = 1:nPC)
testis <- FindClusters(testis, resolution = res)
# dimentional reduction 
testis <- RunUMAP(testis, dims = 1:nPC)
testis <- RunTSNE(testis, dims = 1:nPC)
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

DimPlot(testis, reduction = "umap", group.by = 'orig.ident', label = F, order = T, pt.size = .5) + 
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

####### MG #######
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

######## save Seurat obj ########
assign(paste('testis_',dataset, sep = ''),testis)
save(testis, file = paste('testis_',dataset, sep = ''))


####### somatic #########
somatic_clu_list[[dataset]] =  c('0','1','2','3','4')
testis_somatic <- subset(testis, idents = somatic_clu_list[[dataset]] )
testis_somatic
table(testis_somatic$seurat_clusters)
assign(paste('testis_somatic_',dataset, sep = ''),testis_somatic)
save(testis_somatic, file = paste('testis_somatic_',dataset, sep = ''))
