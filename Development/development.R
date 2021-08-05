#' development integration

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
datasets_dir = ('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/')
OBJ_dir = '/Users/tascini.annasofia/Dropbox (HSR Global)/WORKSPACE/Alfano/AlfanoM_904_infertilita_epigenetics/7_bioinfo/SeuratObjects/'
SO_dir = '/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/Seurat_objects/'
#load(paste(OBJ_dir,"Somatic_integration_beforeplots.RData", sep=''))

######## my functions #######
source("mySeuratfunctions.R")
source('./utility_functions.R')
source('./integration_function.R')


####### mycolor ####
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

#------------- my Gene -------------
# gene OSR
gtf_OSR = rtracklayer::import.gff('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/datasets/Homo_sapiens.GRCh38.97.gtf')
gtf_OSR_df = as.data.frame(gtf_OSR)[,c("gene_id", "gene_name")]
gtf_OSR_dictionary = gtf_OSR_df %>% distinct(gene_id, gene_name, .keep_all = TRUE)
wrong_genes_OSR = setdiff(Tascini_Genes, gtf_OSR_dictionary$gene_name)
if (length(wrong_genes_OSR) != 0) {print("Check OSR genes")} 


#------------- GSE124263 ---------------------------------------------------

# MK infant testis (2 samples)
# Genome: hg19 CellRanger v2.2.0
# Cellranger References - 2.1.0 (February 7, 2018)
FT_GSE124263 = FALSE
if (FT_GSE124263) {
  dataset = 'GSE124263'
  load(paste(datasets_dir, 'testis_',dataset, sep = ''))
  testis_somatic <- subset(testis, idents = c('0','3','4','5','8','11','15'))
  assign(paste('testis_somatic_',dataset, sep = ''), testis_somatic)
  Idents(testis_somatic) <- "donor"
  testis_somatic_infant <- subset(testis_somatic, idents = c('neonatal_day2', 'neonatal_day7'))
  Idents(testis_somatic_infant) <- "seurat_clusters"
  assign(paste('testis_somatic_infant_',dataset, sep = ''), testis_somatic_infant)
  DimPlot(testis_somatic_infant_GSE124263, label = TRUE)
  #dataset = "GSE124263"
  # Gene Conversion
  counts_file_list = list()
  counts_file_list[["GSE124263"]] = '/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/datasets/GSE124263_RAW/'
  start_time <- Sys.time()
  matrix_sparse <- Read10X(data.dir = c(paste(counts_file_list[[dataset]],'neonatal_day2', sep=''),
                                        paste(counts_file_list[[dataset]],'neonatal_day7', sep=''),
                                        paste(counts_file_list[[dataset]],'adult_Donor1', sep=''),
                                        paste(counts_file_list[[dataset]],'adult_Donor2', sep='')),
                           gene.column = 1)
  end_time <- Sys.time()
  end_time - start_time
  d_Genes = row.names(matrix_sparse)
  
  # GRCh37  Cell Ranger software (v2.2.0)
  gtf = rtracklayer::import.gff('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/datasets/Homo_sapiens.GRCh37.82.gtf')
  gtf_df = as.data.frame(gtf)[,c("gene_id", "gene_name")]
  # gene dataset
  gtf_dictionary = gtf_df %>% distinct(gene_id, gene_name, .keep_all = TRUE)
  d_Genes_name = gtf_dictionary[gtf_dictionary$gene_id %in% d_Genes,]$gene_name
  row.names(matrix_sparse) <- d_Genes_name 
  wrong_genes = setdiff(d_Genes, gtf_dictionary$gene_id)
  if (length(wrong_genes) != 0) {print("Check dataset genes")} 
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
  p1 = venn::venn(list(GSE124263 = gene$old_name, OSR=Tascini_Genes))
  p2 = venn::venn(list(GSE124263 = gene$new_name, OSR=Tascini_Genes))
  #
  matrix_sparse_c = matrix_sparse
  row.names(matrix_sparse_c) = gene$new_name
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
  #cell to select
  whichcell = row.names(testis_somatic_infant_GSE124263[[]])
  testis_GSE124263.raw.somatic.infant = subset(testis_GSE124263.raw, cells = whichcell)
  testis_GSE124263.raw.somatic.infant$development = 'neonatal'
  testis_GSE124263.raw.somatic.infant$GEO = 'GSE124263'
  testis_GSE124263.raw.somatic.infant$source <- paste(testis_GSE124263.raw.somatic.infant$GEO, 
                                                      testis_GSE124263.raw.somatic.infant$donor, sep = "_")
  head(testis_GSE124263.raw.somatic.infant[[]])
  testis_GSE124263.raw.somatic.infant <- RenameCells(testis_GSE124263.raw.somatic.infant, add.cell.id = 'GSE124263_')
  save(testis_GSE124263.raw.somatic.infant, file = paste(SO_dir, "testis_GSE124263.raw.somatic.infant", sep =''))

} else { 
  load(paste(SO_dir, "testis_GSE124263.raw.somatic.infant", sep =''))
}

#------------- GSE120506 ----------------------------------------------------
# Guo infant testis (1 sample)
# Cellranger References - 3.0.0 (February 7, 2018)
FT_GSE120506 = FALSE
if (FT_GSE120506) {
  dataset = 'GSE120506'
  load(paste(datasets_dir, 'testis_',dataset, sep = ''))
  testis_somatic <- testis # all cells
  assign(paste('testis_somatic_',dataset, sep = ''), testis_somatic)
  DimPlot(testis_somatic_GSE120506, label = TRUE)
  counts_file_list[["GSE120506"]] = '/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/datasets/GSE120506_infant_combined_UMI.txt'
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
  gtf = rtracklayer::import.gff('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/datasets/Homo_sapiens.GRCh38.84.gtf')
  gtf_df = as.data.frame(gtf)[,c("gene_id", "gene_name")]
  # gene dataset
  gtf_dictionary = gtf_df %>% distinct(gene_id, gene_name, .keep_all = TRUE)
  wrong_genes = setdiff(d_Genes, gtf_dictionary$gene_name)
  if (length(wrong_genes) != 0) {print("Check dataset genes")} 
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
  p1 = venn::venn(list(GSE120506 = gene$old_name, OSR=Tascini_Genes))
  p2 = venn::venn(list(GSE120506 = gene$new_name, OSR=Tascini_Genes))
  
  matrix_sparse_c = matrix_sparse
  row.names(matrix_sparse_c) = gene$new_name
  testis_GSE120506.raw.somatic <- CreateSeuratObject(counts       = matrix_sparse_c, 
                                                     project      = dataset, 
                                                     min.cells    = 5, 
                                                     min.features = 200)
  testis_GSE120506.raw.somatic
  testis = testis_GSE120506.raw.somatic
  testis_GSE120506.raw.somatic$development = 'infant'
  testis_GSE120506.raw.somatic$GEO = 'GSE120506'
  testis_GSE120506.raw.somatic$source = 'GSE120506_infant'
  save(testis_GSE120506.raw.somatic, file = paste(SO_dir, "testis_GSE120506.raw.somatic", sep =''))
} else { 
  load(paste(SO_dir, "testis_GSE120506.raw.somatic", sep =''))
  testis_GSE120506.raw.somatic
  }

#------------- GSE134144 -----------------------------------------------------------------
# Pubertal testis - Guo et al (from juvenile males of 7, 11, 13, and 14 years old)
FT_GSE134144 = FALSE
if (FT_GSE134144) {
  dataset = 'GSE134144'
  load(paste(datasets_dir, 'testis_',dataset, sep = ''))
  if (dataset == 'GSE134144'){
    barcodes_label = data.frame(barcode = row.names(testis[[]]), 
                                donor = str_sub(row.names(testis[[]]),1,3))
    testis$donor =  barcodes_label$donor
  }
  testis_somatic <- subset(testis, idents = c('0','1','2','3','4','6','7','11'))
  assign(paste('testis_somatic_',dataset, sep = ''), testis_somatic)
  
  # Gene Conversion
  counts_file_list[["GSE134144"]] = '/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/datasets/GSE134144_Pubertal_combined_UMI.txt'
  start_time <- Sys.time()
  matrix_sparse <- as.sparse(read.table(counts_file_list[[dataset]], 
                                        sep="\t", 
                                        header = T, 
                                        quote = "", 
                                        row.names = 1))
  end_time <- Sys.time()
  end_time - start_time
  d_Genes = row.names(matrix_sparse)
  
  # GRCh38  Cell Ranger software (v2.2.2)
  gtf <- rtracklayer::import.gff('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/datasets/Homo_sapiens.GRCh38.84.gtf')
  gtf_df = as.data.frame(gtf)[,c("gene_id", "gene_name")]
  gtf_dictionary = gtf_df %>% distinct(gene_id, gene_name, .keep_all = TRUE)
  wrong_genes = setdiff(d_Genes, gtf_dictionary$gene_name)
  if (length(wrong_genes) != 0) {print("Check dataset genes")}
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
  
  p1 = venn::venn(list(GSE134144 = gene$old_name, OSR=Tascini_Genes))
  p2 = venn::venn(list(GSE134144 = gene$new_name, OSR=Tascini_Genes))
  #
  matrix_sparse_c = matrix_sparse
  row.names(matrix_sparse_c) = gene$new_name
  testis_GSE134144.raw <- CreateSeuratObject(counts       = matrix_sparse_c, 
                                             project      = dataset, 
                                             min.cells    = 5, 
                                             min.features = 200)
  testis_GSE134144.raw
  testis = testis_GSE134144.raw
  if (dataset == 'GSE134144'){
    barcodes_label = data.frame(barcode = row.names(testis[[]]), 
                                donor = str_sub(row.names(testis[[]]),1,3))
    testis$donor =  barcodes_label$donor
  }
  testis_GSE134144.raw = testis
  whichcell = row.names(testis_somatic_GSE134144[[]])
  testis_GSE134144.raw.somatic= subset(testis_GSE134144.raw, cells = whichcell)
  testis_GSE134144.raw.somatic$development = 'prepubertal'
  testis_GSE134144.raw.somatic$GEO = 'GSE134144'
  testis_GSE134144.raw.somatic$source <- paste(testis_GSE134144.raw.somatic$GEO, 
                                               testis_GSE134144.raw.somatic$donor, sep = "_")
  head(testis_GSE134144.raw.somatic[[]])
  testis_GSE134144.raw.somatic <- RenameCells(testis_GSE134144.raw.somatic, add.cell.id = 'GSE134144_')
  save(testis_GSE134144.raw.somatic, file = paste(SO_dir, "testis_GSE134144.raw.somatic", sep =''))
} else {
  load(paste(SO_dir, "testis_GSE134144.raw.somatic", sep =''))
  testis_GSE134144.raw.somatic
}

#----------- Load adult testis ---------------------------------------------------
load(paste(SO_dir, "adult.somatic.list", sep=''))
adult.somatic.list

#-------------young testis list ---------------------------------------
# Separate donors
# GSE124263
SO = testis_GSE124263.raw.somatic.infant
Idents(SO) <- "donor"
testis_GSE124263.raw.somatic.infant_day2 = subset(SO, idents = "neonatal_day2")
Idents(testis_GSE124263.raw.somatic.infant_day2) <- "orig.ident"
testis_GSE124263.raw.somatic.infant_day7 = subset(SO, idents = "neonatal_day7")
Idents(testis_GSE124263.raw.somatic.infant_day7) <- "orig.ident"
# GSE134144
SO = testis_GSE134144.raw.somatic
Idents(SO) <- "donor"
table(SO$donor)
testis_GSE134144.raw.somatic_y13 = subset(SO, idents = "Y13")
Idents(testis_GSE134144.raw.somatic_y13) <- "orig.ident"
testis_GSE134144.raw.somatic_y11 = subset(SO, idents = "Y11")
Idents(testis_GSE134144.raw.somatic_y11) <- "orig.ident"
testis_GSE134144.raw.somatic_y14 = subset(SO, idents = "Y14")
Idents(testis_GSE134144.raw.somatic_y14) <- "orig.ident"
testis_GSE134144.raw.somatic_y7 = subset(SO, idents = "Y7.")
Idents(testis_GSE134144.raw.somatic_y7) <- "orig.ident"

## Young list
young.somatic.list = list(
  testis_GSE124263.raw.somatic.infant_day2,
  testis_GSE124263.raw.somatic.infant_day7,
  testis_GSE120506.raw.somatic,
  testis_GSE134144.raw.somatic_y7,
  testis_GSE134144.raw.somatic_y11,
  testis_GSE134144.raw.somatic_y13,
  testis_GSE134144.raw.somatic_y14
  )
names(young.somatic.list) <- c("GSE12426_infant_day2", "GSE12426_infant_day7",
                               "GSE120506_infant",
                               "GSE134144_y7", "GSE134144_y11", "GSE134144_y13", "GSE134144_y14")

save(young.somatic.list, file = paste(SO_dir, "young.somatic.list", sep=''))

for (i in names(young.somatic.list)) {
  young.somatic.list[[i]]$condition = "CTL"
}


######## morris list ########
load(paste(SO_dir, "iNOA_Morris_v2.integrated", sep=''))
iNOA_Morris.integrated = final_int_obj
list_Morris_iNOA = SplitObject(iNOA_Morris.integrated, 
                               split.by = "source")

head(list_Morris_iNOA$`Morris donor`[[]])
morris = list_Morris_iNOA$`Morris donor`
morris$condition <- NULL
morris$condition <- 'Morris'
morris$development <- 'Morris'
morris_list = list(morris)
names(morris_list) = "Morris_Donor"

#------------- development all list -------------
somatic.list = c(adult.somatic.list,young.somatic.list)
save(somatic.list, file = paste(SO_dir, "somatic.list_all.rds", sep=''))

#------------- integration ------------
options(future.globals.maxSize = 3145728000) 
dat = 'Final_Development_2'
dir = paste('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/',dat,'/', sep ='')

start_time <- Sys.time()
int_obj = integration_workflow(somatic.list = somatic.list, 
                               dataset_iNOA_dir = dir, 
                               dataset = dat)
end_time <- Sys.time()
end_time - start_time

int_obj
DimPlot(int_obj, label = T, order = T)
new_cluster_id = c('LEY', 'MYD', 'SRT', 'END', 'MCR', 'STRO', "END",'TCL')
names(new_cluster_id) <- levels(int_obj)
dev_int <- RenameIdents(int_obj, new_cluster_id)
dev_int$celltype <- Idents(dev_int)

col = as.character(colors_Alf[levels(dev_int),]$cols)
umap = DimPlot(dev_int, reduction = "umap", label = TRUE, pt.size = 1, order = TRUE) +
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
ggsave(filename = paste(dir, 'UMAP.png', sep=''), umap,
       width = 15, height = 12,
       units = 'cm')
umap

save(dev_int, file = paste(SO_dir, "int_dev.rds", sep =''))

load(paste(SO_dir, "int_dev.rds", sep =''))

ND = length(table(dev_int$source))

dev_fix = dev_int@meta.data[,c("development", "condition")]
table(dev_fix$development)
table(dev_fix$condition)

dev_fix$development[dev_fix$development == "infant"] = "neonatal"
dev_fix$development[dev_fix$condition == "iNOA"] = "iGCA"
dev_int$development = dev_fix$development
dev_int$development = factor(dev_int$development, levels = c('iGCA','neonatal','prepubertal', 'adult'))

final_int_obj = dev_int
A <- as.data.frame(prop.table(table(final_int_obj$celltype, 
                                    final_int_obj$development),margin = 2))
A2 <- as.data.frame(table(final_int_obj$celltype, 
                          final_int_obj$development),margin = 2)
A$N <- A2$Freq

col = as.character(colors_Alf[levels(final_int_obj$celltype),]$cols)
N <- ggplot(data=A, aes(x=Var2, y=N, fill=Var1)) +
  geom_bar(stat="identity",  color="black", position="fill") + 
  scale_fill_manual(values = col) +
  theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=12, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
        axis.title.y = element_text(face = "bold", color = "black", size = 14),
        legend.title = element_text(size = 0),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "left",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "", y = "fraction of cells")
N
A <- as.data.frame(prop.table(table(final_int_obj$celltype, 
                                    final_int_obj$source),margin = 2))
A2 <- as.data.frame(table(final_int_obj$celltype, 
                          final_int_obj$source),margin = 2)
A$N <- A2$Freq

col = as.character(colors_Alf[levels(final_int_obj$celltype),]$cols)
N2 <- ggplot(data=A, aes(x=Var2, y=N, fill=Var1)) +
  geom_bar(stat="identity",  color="black", position="fill") + 
  scale_fill_manual(values = col) +
  theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=8, hjust = .5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
        axis.title.y = element_text(face = "bold", color = "black", size = 14),
        legend.title = element_text(size = 0),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "left",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "", y = "fraction of cells")
N + N2 + plot_layout(guides = "collect")

pdf(paste(dir,'N_iNOAvsCTL.pdf',sep=''),width=18, height=6)
print(N + N2 + plot_layout(guides = "collect"))
dev.off()


#--------- DLK1 ------------
Seurat.object = final_int_obj
Seurat.object$source = factor(Seurat.object$source, 
                                levels = c("iNOA donor 1",
                                           "iNOA donor 2",
                                           "iNOA donor 3",
                                           "GSE124263_neonatal_day2",
                                           "GSE124263_neonatal_day7",
                                           "GSE120506_infant",
                                           "GSE134144_Y7.",
                                           "GSE134144_Y11",
                                           "GSE134144_Y13",
                                           "GSE134144_Y14",
                                           "GSE112013_adult_Donor1",
                                           "GSE112013_adult_Donor2",
                                           "GSE112013_adult_Donor3",
                                           "GSE124263_adult_Donor1",
                                           "GSE124263_adult_Donor2"))
DefaultAssay(Seurat.object) <- "RNA"
Idents(Seurat.object) <- "celltype"
Seurat.object_s <- subset(Seurat.object,  idents = c("LEY", "MYD"))
vp_D <- VlnPlot(Seurat.object_s,
                feature = 'DLK1',
                split.plot = F,
                pt.size = 0,
                split.by = 'development') +
  ylab("Expression Level") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  #geom_boxplot(width = 0.1, outlier.size=1) +
  scale_fill_manual(values = c('orange','#FFFFFF','#99FFFF', 'dodgerblue2'))
vp_D
ggsave(filename = paste(dir, "DLK1_bydevelp.png", sep =''), vp_D, width = 10, height = 5)
table(Seurat.object$source)

vp_D <- VlnPlot(Seurat.object_s,
                feature = 'DLK1',
                split.plot = F,
                pt.size = 0,
                split.by = 'source') +
  ylab("Expression Level") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 5),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  #geom_boxplot(width = 0.1, outlier.size=1) +
  scale_fill_manual(values = c('orange','orange','orange',
                               '#FFFFFF','#FFFFFF','#FFFFFF',
                               '#99FFFF','#99FFFF','#99FFFF','#99FFFF',
                               'dodgerblue2','dodgerblue2','dodgerblue2','dodgerblue2','dodgerblue2'))
vp_D
ggsave(filename = paste(dir, "DLK1_bysource.png", sep =''), vp_D, width = 15, height = 5)

#--------- INSL3 exp cells ------------
Seurat.object_s <- subset(Seurat.object,  idents = c("LEY"), subset = INSL3 >0)
vp_D <- VlnPlot(Seurat.object_s,
                feature = 'INSL3',
                split.plot = F,
                pt.size = 0,
                split.by = 'development') +
  ylab("Expression Level") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  #geom_boxplot(width = 0.1, outlier.size=1) +
  scale_fill_manual(values = c('orange','#FFFFFF','#99FFFF', 'dodgerblue2'))
vp_D
ggsave(filename = paste(dir, "INSL3_bydevelp.png", sep =''), vp_D, width = 6, height = 5)

vp_D <- VlnPlot(Seurat.object_s,
                feature = 'INSL3',
                split.plot = F,
                pt.size = 0,
                split.by = 'source') +
  ylab("Expression Level") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 5),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  #geom_boxplot(width = 0.1, outlier.size=1) +
  scale_fill_manual(values = c('orange','orange','orange',
                               '#FFFFFF','#FFFFFF','#FFFFFF',
                               '#99FFFF','#99FFFF','#99FFFF','#99FFFF',
                               'dodgerblue2','dodgerblue2','dodgerblue2','dodgerblue2','dodgerblue2'))
vp_D
ggsave(filename = paste(dir, "INSL3_bysource.png", sep =''), vp_D, width = 7, height = 5)


#--------- INSL3 allp cells ------------
Seurat.object_s <- subset(Seurat.object,  idents = c("LEY"))
vp_D <- VlnPlot(Seurat.object_s,
                feature = 'INSL3',
                split.plot = F,
                pt.size = 0,
                split.by = 'development') +
  ylab("Expression Level") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  #geom_boxplot(width = 0.1, outlier.size=1) +
  scale_fill_manual(values = c('orange','#FFFFFF','#99FFFF', 'dodgerblue2'))
vp_D
ggsave(filename = paste(dir, "INSL3_bydevelp_allCell.png", sep =''), vp_D, width = 6, height = 5)

vp_D <- VlnPlot(Seurat.object_s,
                feature = 'INSL3',
                split.plot = F,
                pt.size = 0,
                split.by = 'source') +
  ylab("Expression Level") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 5),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  #geom_boxplot(width = 0.1, outlier.size=1) +
    scale_fill_manual(values = c('orange','orange','orange',                                '#FFFFFF','#FFFFFF','#FFFFFF',                                '#99FFFF','#99FFFF','#99FFFF','#99FFFF',                                'dodgerblue2','dodgerblue2','dodgerblue2','dodgerblue2','dodgerblue2'))
vp_D
ggsave(filename = paste(dir, "INSL3_bysource_allCell.png", sep =''), vp_D, width = 7, height = 5)

#-------- IGF2 ---------------------
gene = c("IGF2", "NOTCH2", "HSD17B3")
for (g in gene) {
  Seurat.object = final_int_obj
  Seurat.object$source = factor(Seurat.object$source, 
                              levels = c("iNOA donor 1",
                                         "iNOA donor 2",
                                         "iNOA donor 3",
                                         "GSE124263_neonatal_day2",
                                         "GSE124263_neonatal_day7",
                                         "GSE120506_infant",
                                         "GSE134144_Y7.",
                                         "GSE134144_Y11",
                                         "GSE134144_Y13",
                                         "GSE134144_Y14",
                                         "GSE112013_adult_Donor1",
                                         "GSE112013_adult_Donor2",
                                         "GSE112013_adult_Donor3",
                                         "GSE124263_adult_Donor1",
                                         "GSE124263_adult_Donor2"))
DefaultAssay(Seurat.object) <- "RNA"
Idents(Seurat.object) <- "celltype"
Seurat.object_s <- subset(Seurat.object,  idents = c("LEY", "MYD"))
vp_D <- VlnPlot(Seurat.object_s,
                feature = g,
                split.plot = F,
                pt.size = 0,
                split.by = 'development') +
  ylab("Expression Level") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  #geom_boxplot(width = 0.1, outlier.size=1) +
  scale_fill_manual(values = c('orange','#FFFFFF','#99FFFF', 'dodgerblue2'))
print(vp_D)
ggsave(filename = paste(dir, g,"_bydevelp.png", sep =''), vp_D, width = 10, height = 5)
table(Seurat.object$source)

vp_D <- VlnPlot(Seurat.object_s,
                feature = g,
                split.plot = F,
                pt.size = 0,
                split.by = 'source') +
  ylab("Expression Level") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 5),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  #geom_boxplot(width = 0.1, outlier.size=1) +
    scale_fill_manual(values = c('orange','orange','orange','#FFFFFF','#FFFFFF','#FFFFFF','#99FFFF','#99FFFF','#99FFFF','#99FFFF','dodgerblue2','dodgerblue2','dodgerblue2','dodgerblue2','dodgerblue2'))
print(vp_D)
ggsave(filename = paste(dir, g,"_bysource.png", sep =''), vp_D, width = 15, height = 5)
}

table(Idents(Seurat.object))
FeaturePlot(subset(Seurat.object, idents = c("LEY")), 
            features = "HSD17B3", split.by = 'develo', order = T)
FeaturePlot(Seurat.object, 
            features = "HSD17B3", split.by = 'develo', order = T)



# ------------ Plot with statistics ---------------
Seurat.object = dev_int
Seurat.object$celltype.dev <- paste(Idents(Seurat.object), Seurat.object$development, sep = "_")
table(Seurat.object$celltype.dev)
c1 = "iGCA"; c2 = "neonatal"; 
features = c("DLK1")
Idents(Seurat.object) <- "celltype.dev"
dge_value_for_figure = data.frame()
for (ct in c("LEY", "MYD")) {
for (c2 in c("neonatal", "prepubertal", "adult")) { 
  value <- FindMarkers(Seurat.object,
                       ident.1 = paste(ct,"_",c1, sep=''),
                       ident.2 = paste(ct,"_",c2, sep=''), 
                       verbose = FALSE,
                       feature = features,
                       logfc.threshold=0,
                       min.pct = 0)
  value$celltype = ct
  value$comparison = paste(c1, "vs",c2, sep="_")
  value$gene = features
  dge_value_for_figure = rbind(dge_value_for_figure, value)
}
}
dge_value_for_figure = dge_value_for_figure[order(dge_value_for_figure$gene),
                                            c("gene", "celltype","comparison", "avg_logFC", "p_val", "p_val_adj", "pct.1","pct.2")]
dge_value_for_figure
write.xlsx(dge_value_for_figure, file = paste(outdir, "Figure5_DLK1_Devstatistics.xlsx", sep=''), 
           row.names = FALSE, asTable = T)


features = c("INSL3")
Idents(Seurat.object) <- "celltype.dev"
DefaultAssay(Seurat.object) <- "RNA"
dge_value_for_figure = data.frame()
for (ct in c("LEY")) {
  for (c2 in c("neonatal", "prepubertal", "adult")) { 
    value <- FindMarkers(subset(Seurat.object, subset = INSL3 >0),
                         ident.1 = paste(ct,"_",c1, sep=''),
                         ident.2 = paste(ct,"_",c2, sep=''), 
                         verbose = FALSE,
                         feature = features,
                         logfc.threshold=0,
                         min.pct = 0)
    value$celltype = ct
    value$comparison = paste(c1, "vs",c2, sep="_")
    value$gene = features
    dge_value_for_figure = rbind(dge_value_for_figure, value)
  }
}

dge_value_for_figure = dge_value_for_figure[order(dge_value_for_figure$gene),
                                            c("gene", "celltype","comparison", "avg_logFC", "p_val", "p_val_adj", "pct.1","pct.2")]
dge_value_for_figure
write.xlsx(dge_value_for_figure, file = paste(outdir, "Figure5_INL3_Devstatistics.xlsx", sep=''), 
           row.names = FALSE, asTable = T)




Seurat.object = dev_int 
DefaultAssay(Seurat.object) <- "RNA"
Idents(Seurat.object) <- "celltype"  
for (ct in c('LEY', 'MYD')) {
  so = subset(Seurat.object, idents = ct) 
  gene = "DLK1"
  so[[gene]] <- FetchData(object = so, vars = gene)
  df= so[[]]
  my_comparisons <- list( c("iGCA", "neonatal"), c("iGCA", "prepubertal"), c("iGCA", "adult") )
  #df$condition = plyr::revalue(df$condition, c("iNOA"="iGCA", "CTL"="OA"))
  v <- ggviolin(df, x = "development", y = gene, 
                fill = "development", color = "black", width = 1,
                #add = c('jitter'), add.params = list(size =.5, shape = 1),
                #              add.params = list(size = .5, shape = 20),
                alpha = 0.8,
  #facet.by = 'celltype', 
  short.panel.labs = TRUE,
  palette = c('orange','#FFFFFF','#99FFFF', 'dodgerblue2'), 
  error.plot = "crossbar", draw_quantiles = c(0.5), trim = T,
  #add = c('mean','mean_se'),
  title = ct, 
  panel.labs.background = list(color = "white", fill = "white", size = 0.5),
  panel.labs.font = list(color = "black", face = "bold", size = 16))
  v = ggpar(v, ylab = 'Expression level', xlab = '')
  v = v + stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 6) +
    theme(plot.title = element_text(color="black", size=22, face="bold.italic", hjust = 0.5),
          plot.title.position =  'panel',
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 14, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 14),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=14),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(color = "black", size = 12),
          legend.position = 'right',
          panel.spacing = unit('-0.1', "lines"))
  assign(paste('v',gene,ct,sep = '_'),v)
  #pdf(paste(dir, 'ImprintedGenes/VlnPlot_r/','vlnPlot_',gene,'.pdf', sep=''), width = 5, height = 4.5)
  print(v)
  #dev.off()
}

v_DLK1_LEY
pdf(paste(outdir, "DLK1_dev_withStatistic.pdf", sep=''), width = 15, height = 6)
print(v_DLK1_LEY + v_DLK1_MYD + plot_layout(guides = "collect"))
dev.off()

pdf(paste(outdir, "DLK1_dev_withStatistic_LEY.pdf", sep=''), width = 8, height = 6)
print(v_DLK1_LEY)
dev.off()
pdf(paste(outdir, "DLK1_dev_withStatistic_MYD.pdf", sep=''), width = 8, height = 6)
print(v_DLK1_MYD)
dev.off()

fig4F = v_DLK1_LEY$data[,c("development", "celltype", "DLK1")]
fig4F$cellID = row.names(fig4F)
df_dir = "/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/Dataframes_Figures/"
write.xlsx(fig4F[,colnames(fig4F)[c(4,1:3)]],paste(df_dir,"Figure4F_DLK1_L.xlsx",sep=''),row.names = F)

fig4F = v_DLK1_MYD$data[,c("development", "celltype", "DLK1")]
fig4F$cellID = row.names(fig4F)
df_dir = "/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/Dataframes_Figures/"
write.xlsx(fig4F[,colnames(fig4F)[c(4,1:3)]],paste(df_dir,"Figure4F_DLK1_M.xlsx",sep=''),row.names = F)


for (ct in c('LEY', 'MYD')) {
  so = subset(Seurat.object, idents = ct, subset = INSL3 >0) 
  gene = "INSL3"
  so[[gene]] <- FetchData(object = so, vars = gene)
  df= so[[]]
  my_comparisons <- list( c("iGCA", "neonatal"), c("iGCA", "prepubertal"), c("iGCA", "adult") )
  #df$condition = plyr::revalue(df$condition, c("iNOA"="iGCA", "CTL"="OA"))
  v <- ggviolin(df, x = "development", y = gene, 
                fill = "development", color = "black", width = 1,
                #add = c('jitter'), add.params = list(size =.5, shape = 1),
                #              add.params = list(size = .5, shape = 20),
                alpha = 0.8,
                #facet.by = 'celltype', 
                short.panel.labs = TRUE,
                palette = c('orange','#FFFFFF','#99FFFF', 'dodgerblue2'), 
                error.plot = "crossbar", draw_quantiles = c(0.5), 
                trim = T,
                title = ct, 
                panel.labs.background = list(color = "white", fill = "white", size = 0.5),
                panel.labs.font = list(color = "black", face = "bold", size = 16))
  v = ggpar(v, ylab = 'Expression level', xlab = '')
  v = v + stat_compare_means(comparisons = my_comparisons, 
                             label = "p.signif", 
                             size = 6) +
    theme(plot.title = element_text(color="black", size=22, face="bold.italic", hjust = 0.5),
          plot.title.position =  'panel',
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 14, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 14),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=14),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(color = "black", size = 12),
          legend.position = 'right',
          panel.spacing = unit('-0.1', "lines"))
  assign(paste('v',gene,ct,sep = '_'),v)
  #pdf(paste(dir, 'ImprintedGenes/VlnPlot_r/','vlnPlot_',gene,'.pdf', sep=''), width = 5, height = 4.5)
  print(v)
  #dev.off()
}

head(v_INSL3_LEY$data[,c("development", "celltype", "INSL3")])
fig4F = v_INSL3_LEY$data[,c("development", "celltype", "INSL3")]
fig4F$cellID = row.names(fig4F)
df_dir = "/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/Dataframes_Figures/"
write.xlsx(fig4F[,colnames(fig4F)[c(4,1:3)]],paste(df_dir,"Figure4F_INSL3.xlsx",sep=''),row.names = F)


pdf(paste(outdir, gene,"_dev_withStatistic_LEY.pdf", sep=''), width = 8, height = 6)
print(v_INSL3_LEY)
dev.off()

path_file = paste(outdir,"DevelopmentGene/statistics/",sep='')
fl = list.files(path_file, pattern=glob2rx("*.pdf"))
fl
#fl = paste("DevelopmentGene/statistcs/","Figure1G.umap.pdf", sep='')
for (file in fl){
  im = image_read_pdf(paste(path_file,file,sep=''),density = 200)
  image_write(im, path = paste(path_file,file,'.tiff', sep=''), format = "tiff")
}




#-------- 144 -------------#
## Young list
testis_GSE120506.raw.somatic$age = "1yo"
testis_GSE134144.raw.somatic_y7$age = "7yo"
testis_GSE134144.raw.somatic_y11$age = "11yo"
testis_GSE134144.raw.somatic_y13$age = "13yo"
testis_GSE134144.raw.somatic_y14$age = "14yo"
testis_GSE112013.raw.somatic_y25 = adult.somatic.list$GSE112013_adult_Donor3
testis_GSE112013.raw.somatic_y25$age = "25yo"
testis_GSE112013.raw.somatic_y25$development = "adult"

GuoDevelop.list = list(
  testis_GSE120506.raw.somatic,
  testis_GSE134144.raw.somatic_y7,
  testis_GSE134144.raw.somatic_y11,
  testis_GSE134144.raw.somatic_y13,
  testis_GSE134144.raw.somatic_y14,
  testis_GSE112013.raw.somatic_y25
)

save(GuoDevelop.list , file = paste(SO_dir, "GuoDevelop.list", sep =''))


