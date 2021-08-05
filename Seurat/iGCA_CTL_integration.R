#' adults integration
#' CTL and iGCA somatic testis cell integration
#' performed with Seurat 

######## libraries ##########
library(plyr)
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
library(enrichR)
library(ggpubr)

###### load data ########
datasets_dir = ('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/')
OBJ_dir = '/Users/tascini.annasofia/Dropbox (HSR Global)/WORKSPACE/Alfano/AlfanoM_904_infertilita_epigenetics/7_bioinfo/SeuratObjects/'
SO_dir = '/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/Seurat_objects/'
load(paste(OBJ_dir,"Somatic_integration_beforeplots.RData", sep=''))

##### my functions ####
source("mySeuratfunctions.R")
source('./utility_functions.R')
source('./integration_function.R')

####### mycolor ####
colors_cellTypes <- data.frame(celltype = c("SSCs","Differentiated S'gonia",
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
row.names(colors_cellTypes) <- colors_cellTypes$celltype


# Download gff files from enseml 
# for the different datasets employed in this analysis
# hg19 v75
#curl ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz -o Homo_sapiens.GRCh37.75.gtf.g
#gunzip Homo_sapiens.GRCh37.75.gtf.gz
# hg38 v97
#curl ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz -o Homo_sapiens.GRCh38.97.gtf.gz
#gunzip Homo_sapiens.GRCh38.97.gtf.gz
# Cellranger References - 3.0.0 (February 7, 2018)
#curl ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz -o Homo_sapiens.GRCh37.87.gtf.gz
#gunzip Homo_sapiens.GRCh37.87.gtf.gz
# Cellranger References - 2.1.0 (February 7, 2018)
#curl ftp://ftp.ensembl.org/pub/grch37/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh37.82.gtf.gz -o Homo_sapiens.GRCh37.82.gtf.gz
#gunzip Homo_sapiens.GRCh37.82.gtf.gz

###### my Gene #######
# gene OSR
gtf_OSR = rtracklayer::import.gff('Homo_sapiens.GRCh38.97.gtf')
gtf_OSR_df = as.data.frame(gtf_OSR)[,c("gene_id", "gene_name")]
gtf_OSR_dictionary = gtf_OSR_df %>% distinct(gene_id, gene_name, .keep_all = TRUE)

# import count matrixes of INOA patients
counts_file = paste('iNOA1.counts.matrix.tsv.gz')
UMImatrix_sparse_iNOA1 = Matrix(as.matrix(read.table(counts_file, 
                                                   sep="\t", 
                                                   header = T, 
                                                   quote = "", 
                                                   row.names = 1)), 
                              sparse=TRUE)
counts_file = paste('iNOA2.counts.matrix.tsv.gz')
UMImatrix_sparse_iNOA2 = Matrix(as.matrix(read.table(counts_file, 
                                                   sep="\t", 
                                                   header = T, 
                                                   quote = "", 
                                                   row.names = 1)), 
                              sparse=TRUE)

counts_file = paste('iNOA3.counts.matrix.tsv.gz')
UMImatrix_sparse_iNOA3 = Matrix(as.matrix(read.table(counts_file, 
                                                   sep="\t", 
                                                   header = T, 
                                                   quote = "", 
                                                   row.names = 1)), 
                              sparse=TRUE)

OSR_Genes = unique(c(row.names(UMImatrix_sparse_iNOA1), 
                         row.names(UMImatrix_sparse_iNOA2), 
                         row.names(UMImatrix_sparse_iNOA3)))

# check correspondence gff, gene name
# wrong_genes_OSR = setdiff(OSR_Genes, gtf_OSR_dictionary$gene_name)
# if (length(wrong_genes_OSR) != 0) {print("Check OSR genes")} 

#------------- GSE124263 ---------------------------------------------------
# 2 samples from Sohni et al. (https://doi.org/10.1016/j.celrep.2019.01.045) 
# Miles Wilkinson group.
# Genome: hg19 CellRanger v2.2.0
# Cellranger References - 2.1.0 (February 7, 2018)

FirstTime_GSE124263 = TRUE
if (FirstTime_GSE124263) {
  dataset = 'GSE124263'
  # Gene Conversion
  counts_file_list = list()
  #GSE124263_RAW folder was download from GEO - GSE124263
  counts_file_list[["GSE124263"]] = 'datasets/GSE124263_RAW/'
  start_time <- Sys.time()
  matrix_sparse <- Read10X(data.dir = c(paste(counts_file_list[[dataset]],'neonatal_day2', sep=''),
                                        paste(counts_file_list[[dataset]],'neonatal_day7', sep=''),
                                        paste(counts_file_list[[dataset]],'adult_Donor1', sep=''),
                                        paste(counts_file_list[[dataset]],'adult_Donor2', sep='')),
                           gene.column = 1)
  end_time <- Sys.time()
  # print time
  end_time - start_time
  d_Genes = row.names(matrix_sparse)
  
  # GRCh37  Cell Ranger software (v2.2.0)
  gtf = rtracklayer::import.gff('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/datasets/Homo_sapiens.GRCh37.82.gtf')
  gtf_df = as.data.frame(gtf)[,c("gene_id", "gene_name")]
  # gene dataset
  gtf_dictionary = gtf_df %>% distinct(gene_id, gene_name, .keep_all = TRUE)
  d_Genes_name = gtf_dictionary[gtf_dictionary$gene_id %in% d_Genes,]$gene_name
  row.names(matrix_sparse) <- d_Genes_name 
  # check correspondence gff, gene name
  # wrong_genes = setdiff(d_Genes, gtf_dictionary$gene_id)
  # if (length(wrong_genes) != 0) {print("Check dataset genes")} 
  
  #conversion
  start_time <- Sys.time()
  conv = sapply(d_Genes_name, Gene_conversion, gtf_dictionary = gtf_dictionary, gtf_OSR_dictionary = gtf_OSR_dictionary)
  end_time <- Sys.time()
  end_time - start_time
  gene = data.frame(old_name = as.character(names(conv)), new_name = as.character(conv), stringsAsFactors = FALSE)
  for (inx in as.numeric(row.names(gene[grep(';',conv),]))) {
    if (length(gtf_OSR_dictionary[gtf_OSR_dictionary$gene_id %in% d_Genes[inx], "gene_name"]) > 0) {
      gene[inx, "new_name"] = gtf_OSR_dictionary[gtf_OSR_dictionary$gene_id %in% d_Genes[inx], "gene_name"]
    } else {
      print(inx)
      gene[inx, "new_name"] = gene[inx, "old_name"]}
  }
  # check intersection
  p1 = venn::venn(list(GSE124263 = gene$old_name, OSR=OSR_Genes))
  p2 = venn::venn(list(GSE124263 = gene$new_name, OSR=OSR_Genes))
  #
  matrix_sparse_c = matrix_sparse
  row.names(matrix_sparse_c) = gene$new_name
  # create seurat object
  testis_GSE124263.raw <- CreateSeuratObject(counts       = matrix_sparse_c, 
                                             project      = dataset, 
                                             min.cells    = 5, 
                                             min.features = 200)
  testis_GSE124263.raw
  testis = testis_GSE124263.raw
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
  testis_GSE124263.raw = testis
  
  # samples was previusly analysed to detect somatic cells
  # [see code GEOdatasets/testis_dataset_all.R]
  # load previously analised sample for selection somatic cells
  load(paste(datasets_dir, 'testis_',dataset, sep = ''))
  testis_somatic <- subset(testis, idents = c('0','3','4','5','8','11','15'))
  assign(paste('testis_somatic_',dataset, sep = ''), testis_somatic)
  Idents(testis_somatic) <- "donor"
  testis_somatic_adult <- subset(testis_somatic, idents = c('adult_Donor1', 'adult_Donor2'))
  Idents(testis_somatic_adult) <- "seurat_clusters"
  assign(paste('testis_somatic_adult_',dataset, sep = ''), testis_somatic_adult)
  DimPlot(testis_somatic_adult_GSE124263, label = TRUE)
  
  # cell to select
  whichcell = row.names(testis_somatic_adult_GSE124263[[]])
  testis_GSE124263.raw.somatic.adult = subset(testis_GSE124263.raw, cells = whichcell)
  testis_GSE124263.raw.somatic.adult$development = 'adult'
  testis_GSE124263.raw.somatic.adult$GEO = 'GSE124263'
  testis_GSE124263.raw.somatic.adult$source <- paste(testis_GSE124263.raw.somatic.adult$GEO,
                                                     testis_GSE124263.raw.somatic.adult$donor, sep = "_")
  head(testis_GSE124263.raw.somatic.adult[[]])
  testis_GSE124263.raw.somatic.adult <- RenameCells(testis_GSE124263.raw.somatic.adult, add.cell.id = 'GSE124263_')
  save(testis_GSE124263.raw.somatic.adult, file = paste(SO_dir, "testis_GSE124263.raw.somatic.adult", sep =''))
  
} else { 
  load(paste(SO_dir, "testis_GSE124263.raw.somatic.adult", sep =''))
}


#------------- GSE112013 ---------------------------------------------------
# Dataset of adult testis from 3 healthy donors 
# [Guo et al.](https://doi.org/10.1038/s41422-018-0099-2).
FirstTime_GSE112013 = FALSE
if (FirstTime_GSE112013) {
  dataset = 'GSE112013'
  # downloaded from GEO
  counts_file_list[['GSE112013']] = 'GSE112013_HealthyAdults_Combined_UMI_table.txt'
  start_time <- Sys.time()
  matrix_sparse <- as.sparse(read.table(counts_file_list[[dataset]], 
                                        sep="\t", 
                                        header = T, 
                                        quote = "", 
                                        row.names = 1))
  end_time <- Sys.time()
  end_time - start_time
  d_Genes = row.names(matrix_sparse)
  # GRCh38  Cell Ranger software (v3.0.0)
  gtf = rtracklayer::import.gff('Homo_sapiens.GRCh38.84.gtf')
  gtf_df = as.data.frame(gtf)[,c("gene_id", "gene_name")]
  # gene dataset
  gtf_dictionary = gtf_df %>% distinct(gene_id, gene_name, .keep_all = TRUE)
  
  # check correspondence gff, gene name
  wrong_genes = setdiff(d_Genes, gtf_dictionary$gene_name)
  if (length(wrong_genes) != 0) {print("Check dataset genes")} 
  # some of the genes do not correspond cause they have .1 at the end of the name
  d_Genes[d_Genes %in% wrong_genes] = gsub(wrong_genes, pattern = "\\.1", replacement = '')
  wrong_genes2 = setdiff(d_Genes, gtf_dictionary$gene_name)
  if (length(wrong_genes2) != 0) {print("Check dataset genes again")}
  
  #conversion
  start_time <- Sys.time()
  conv = sapply(d_Genes, Gene_conversion, gtf_dictionary = gtf_dictionary, gtf_OSR_dictionary = gtf_OSR_dictionary)
  end_time <- Sys.time()
  end_time - start_time
  gene = data.frame(old_name = as.character(names(conv)), new_name = as.character(conv), stringsAsFactors = FALSE)
  for (inx in as.numeric(row.names(gene[grep(';',conv),]))) {
    if (length(gtf_OSR_dictionary[gtf_OSR_dictionary$gene_id %in% d_Genes[inx], "gene_name"]) > 0) {
      gene[inx, "new_name"] = gtf_OSR_dictionary[gtf_OSR_dictionary$gene_id %in% d_Genes[inx], "gene_name"]
    } else {
      print(inx)
      gene[inx, "new_name"] = gene[inx, "old_name"]}
  }
  gene[grep(';',conv),]
  p1 = venn::venn(list(GSE112013 = gene$old_name, OSR=OSR_Genes))
  p2 = venn::venn(list(GSE112013 = gene$new_name, OSR=OSR_Genes))
  
  matrix_sparse_c = matrix_sparse
  row.names(matrix_sparse_c) = gene$new_name
  testis_GSE112013.raw <- CreateSeuratObject(counts       = matrix_sparse_c, 
                                             project      = dataset, 
                                             min.cells    = 5, 
                                             min.features = 200)
  testis_GSE112013.raw
  #select somatic cells
  # samples was previusly analysed to detect somatic cells
  # [see code GEOdatasets/testis_dataset_all.R]
  # load previously analised sample for selection somatic cells
  load(paste(datasets_dir, 'pz.GSE112013.somatic', sep=''))
  whichcell = row.names(pz.GSE112013.somatic[[]])
  testis_GSE112013.raw.somatic = subset(testis_GSE112013.raw, cells = whichcell)
  
  testis_GSE112013.raw.somatic$development = 'adult'
  testis_GSE112013.raw.somatic$condition = 'CTL'
  testis_GSE112013.raw.somatic$GEO = 'GSE112013'
  save(testis_GSE112013.raw.somatic, file = paste(SO_dir, "testis_GSE112013.raw.somatic", sep =''))
  
} else {
  load(paste(SO_dir, "testis_GSE112013.raw.somatic", sep =''))
}

# ------------ integration list -----------------------
#' extract Donors
#' GSE124263
testis_GSE124263.raw.somatic.adult$condition = "CTL" 
Idents(testis_GSE124263.raw.somatic.adult) <- "donor"
testis_GSE124263.raw.somatic_D1 = subset(testis_GSE124263.raw.somatic.adult, idents = "adult_Donor1")
Idents(testis_GSE124263.raw.somatic_D1) <- "orig.ident"
testis_GSE124263.raw.somatic_D2 = subset(testis_GSE124263.raw.somatic.adult, idents = "adult_Donor2")
Idents(testis_GSE124263.raw.somatic_D2) <- "orig.ident"
head(testis_GSE124263.raw.somatic.adult[[]])

#' GSE112013
D1_cells = row.names(testis_GSE112013.raw.somatic[[]])[grep("Donor1", row.names(testis_GSE112013.raw.somatic[[]]))]
D2_cells = row.names(testis_GSE112013.raw.somatic[[]])[grep("Donor2", row.names(testis_GSE112013.raw.somatic[[]]))]
D3_cells = row.names(testis_GSE112013.raw.somatic[[]])[grep("Donor3", row.names(testis_GSE112013.raw.somatic[[]]))]
testis_GSE112013_D1.somatic = subset(testis_GSE112013.raw.somatic, cells = D1_cells); 
testis_GSE112013_D1.somatic$source = "GSE112013_adult_Donor1"
testis_GSE112013_D2.somatic = subset(testis_GSE112013.raw.somatic, cells = D2_cells); 
testis_GSE112013_D2.somatic$source = "GSE112013_adult_Donor2"
testis_GSE112013_D3.somatic = subset(testis_GSE112013.raw.somatic, cells = D3_cells); 
testis_GSE112013_D3.somatic$source = "GSE112013_adult_Donor3"

#' iNOA
# load previously analised sample for selection somatic cells
load(paste(datasets_dir, 'iNOA_donor1', sep=''))
iNOA_donor1 <- RenameCells(iNOA_donor1, add.cell.id = 'iNOA1')
iNOA_donor1$development = 'iNOA'

load(paste(datasets_dir, 'iNOA_donor2', sep=''))
iNOA_donor2 <- RenameCells(iNOA_donor2, add.cell.id = 'iNOA2')
iNOA_donor2$development = 'iNOA'

load(paste(datasets_dir, 'iNOA_donor3', sep=''))
iNOA_donor3 <- RenameCells(iNOA_donor3, add.cell.id = 'iNOA3')
iNOA_donor3$development = 'iNOA'

#' list
adult.somatic.list = list(iNOA_donor1, iNOA_donor2, iNOA_donor3,
                          testis_GSE124263.raw.somatic_D1, testis_GSE124263.raw.somatic_D2,
                          testis_GSE112013_D1.somatic, testis_GSE112013_D2.somatic, testis_GSE112013_D3.somatic)

names(adult.somatic.list) = c('iNOA1', 'iNOA2', 'INOA3', 
                        'GSE124263_adult_Donor1', 'GSE124263_adult_Donor2',
                        'GSE112013_adult_Donor1', 'GSE112013_adult_Donor2','GSE112013_adult_Donor3')

save(adult.somatic.list, file = paste(SO_dir, "adult.somatic.list", sep=''))

# -------------- integration ----------------------

options(future.globals.maxSize = 3145728000) 
dat = "cleanAdultIntegration"
dir = paste('ric.cosr/AlfanoM_904_infertilita_epigenetics/',dat,'/', sep ='')

start_time <- Sys.time()
adult.somatic.list
int_obj = integration_workflow(somatic.list = adult.somatic.list, 
                               dataset_iNOA_dir = dir, 
                               dataset = dat)
end_time <- Sys.time()
end_time - start_time
DimPlot(int_obj, label = T, order = T)

# ----------- Immuno cell assignment check ----------------
TC = subset(int_obj, idents = '6')
table(TC$source)
DimPlot(MCR, split.by = 'condition', group.by = 'source')

DefaultAssay(int_obj) <- "RNA" 

# machrophage cells
MCR = subset(int_obj, idents = c('3', '6'))
matrixTcells = GetAssayData(MCR)
# cell identity prediction from marker gene expression 
MACRO.pred = read.table(file = "testis_MCR_prediction.tsv", quote = "")
WhichTcells = row.names(MACRO.pred[MACRO.pred$labels %in% c("T_cells"),])
DimPlot(int_obj, cells.highlight = WhichTcells, label = T)
DimPlot(subset(int_obj, cells = WhichTcells))

metadata = int_obj[[]]
metadata$newIDENT = as.character(metadata$seurat_clusters)
metadata$newIDENT[metadata$newIDENT == '6'] <- '3'
metadata[WhichTcells_f,]$newIDENT = 6

int_obj$Tcell_reassigned = factor(metadata$newIDENT)
table(int_obj$Tcell_reassigned)
DimPlot(int_obj, group.by = "Tcell_reassigned", order = T)
Idents(int_obj) <- "Tcell_reassigned"
int_obj$seurat_clusters <- int_obj$Tcell_reassigned

# ----------- integration post processing --------------
DefaultAssay(int_obj) <- "integrated"
table(Idents(int_obj))

cell_type <- c( "LEY",
                "MYD", 
                "END",
                "MCR",
                "STRO",
                "SRT",
                "TCL",
                "UND")

final_int_obj = integration_postprocessing(SO = int_obj,
                                           new_cluster_id = cell_type, 
                                           dataset = dat, 
                                           dataset_iNOA_dir = dir,
                                           res = 0.2, 
                                           c1 = "iNOA", c2 = "CTL",
                                           colors_cellTypes = colors_cellTypes)

# ------------------- MG ------------------------
Idents(final_int_obj) <- "celltype"
DimPlot(final_int_obj, split.by = 'condition')
DefaultAssay(final_int_obj) <- "integrated"
tc = FindMarkers(final_int_obj, ident.1 = 'TCL')
tc  %>% top_n(n = 10, wt = avg_logFC)
cluster.markers = FindAllMarkers(final_int_obj)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
filename =  paste("integration_",dataset,"_res",res,"_celltype.xlsx", sep = '')
write.xlsx(cluster.markers,
           file= paste(dir, filename, sep =''), 
           row.names = T,
           asTable = T)
HM =  DoHeatmap(final_int_obj, features = top10$gene, angle = 90) + NoLegend() + 
  scale_fill_gradientn(colours = coolwarm(200)) 
HM
png(paste(dir,"HM_celltype_integration_",dat,"_res",res,".png",sep=''),  width = 1400, height = 1600)
plot(HM)
dev.off()

# ------------------- DGE ------------------------
dir_DGE = paste(dir,'DGE/', sep='')
dir.create(dir_DGE)
SO = final_int_obj
DefaultAssay(SO) <- "RNA"
# relabel cells to perform DGE
SO$celltype.cond <- paste(Idents(SO), SO$condition, sep = "_")
SO$celltype <- Idents(SO)
Idents(SO) <- "celltype.cond"
new_cluster_id = cell_type
TCL_DGE = FALSE
c1 = 'iNOA'; c2='CTL'
if (!TCL_DGE) {new_cluster_id = new_cluster_id[new_cluster_id != 'TCL']}
DGEresults = lapply(new_cluster_id, 
                    myDGE, 
                    Seurat.object = SO, 
                    outdir = dir_DGE,
                    c1 = c1,
                    c2 = c2)
names(DGEresults) <- new_cluster_id
filename = paste(dir, "DGEall_integration_",dat,"_res",res,".xlsx", sep = '')
write.xlsx(x = DGEresults, file = filename, asTable = T, row.names = T)

# --------- Save object -------------------
cleanAdultIntegration = final_int_obj
save(cleanAdultIntegration, file = paste(SO_dir, dat, sep=''))