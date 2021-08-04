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
library(rtracklayer)


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
                                     "STRO",
                                     "UND",
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
                                'cadetblue2',
                                '#999999',
                                'grey'))
rownames(colors_Alf) <- colors_Alf$celltype


outdir = '/Users/tascini.annasofia/Dropbox (HSR Global)/Alfano_904_paperDraft/Final_updateiNOAvsCTL/'
# functions
source('/lustre2/scratch/bioinfotree/common/bioinfotree/prj/AlfanoM_904_infertilita_epigenetics/local/src/mySeuratfunctions.R')
local = T
if (local) {
  data.path = '/Users/tascini.annasofia/Documents/AlfanoM/'
} else {
  data.path = "/Users/tascini.annasofia/ctgb_cluster_root/lustre2/scratch/bioinfotree/common/bioinfotree/prj/AlfanoM_904_infertilita_epigenetics/dataset/20200110/seraut/"
}


############ convert Gene symbols #################

counts_file=('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/datasets/GSE112013_HealthyAdults_Combined_UMI_table.txt')
UMImatrix_sparse_Guo = Matrix(as.matrix(read.table(counts_file, 
                                               sep="\t", 
                                               header = T, 
                                               quote = "", 
                                               row.names = 1)), 
                          sparse=TRUE)
Guo_Genes = row.names(UMImatrix_sparse_Guo)


counts_file = paste('/Users/tascini.annasofia/Documents/AlfanoM/198.counts.matrix.tsv.gz')
UMImatrix_sparse_198 = Matrix(as.matrix(read.table(counts_file, 
                                               sep="\t", 
                                               header = T, 
                                               quote = "", 
                                               row.names = 1)), 
                          sparse=TRUE)
UMImatrix_sparse_198 = UMImatrix_sparse_198[which(row.names(UMImatrix_sparse_198)!= 'gene'),]

counts_file = paste('/Users/tascini.annasofia/Documents/AlfanoM/205.counts.matrix.tsv.gz')
UMImatrix_sparse_205 = Matrix(as.matrix(read.table(counts_file, 
                                                   sep="\t", 
                                                   header = T, 
                                                   quote = "", 
                                                   row.names = 1)), 
                              sparse=TRUE)

counts_file = paste('/Users/tascini.annasofia/Documents/AlfanoM/204sx.counts.matrix.tsv.gz')
UMImatrix_sparse_204 = Matrix(as.matrix(read.table(counts_file, 
                                                   sep="\t", 
                                                   header = T, 
                                                   quote = "", 
                                                   row.names = 1)), 
                              sparse=TRUE)
Tascini_Genes = unique(c(row.names(UMImatrix_sparse_198), 
                       row.names(UMImatrix_sparse_204), 
                       row.names(UMImatrix_sparse_205)))

venn::venn(list(guo = Guo_Genes, tasci=Tascini_Genes))

# import Guoetal cellranger gft
gtf_Guo <- rtracklayer::import.gff('/Users/tascini.annasofia/Downloads/genes.gtf')
gtf_Guo_df = as.data.frame(gtf_Guo)[,c("gene_id", "gene_name")]
gtf_Guo_df_f = gtf_Guo_df %>% distinct(gene_id, gene_name, .keep_all = TRUE)

#fix problem in Guo genes
wrong_genes = setdiff(Guo_Genes, gtf_Guo_df_f$gene_name)
real_name = unique(unlist(strsplit(wrong_genes, split = '.', fixed = TRUE)))
real_name = real_name[real_name != 1]

M = UMImatrix_sparse_Guo
for (i in 1:length(wrong_genes)) {
  if (real_name[i] %in% Guo_Genes) {
    M[real_name[i],] = colSums(M[c(wrong_genes[i], real_name[i]),], 
                               na.rm = FALSE, dims = 1, sparseResult = T)
  }
}
M = M[!(row.names(M) %in% wrong_genes),]
UMImatrix_sparse_Guo = M
library(biomaRt)
# conversion
ensembl97 = useMart("ensembl",
                    dataset="hsapiens_gene_ensembl", 
                    host = 'http://jul2019.archive.ensembl.org')

listFilters(ensembl97)
dictionary_Tascini = getBM(attributes= c('ensembl_gene_id', 'external_gene_name'),
                           filter = 'external_gene_name',
                           values = Tascini_Genes,
                           mart = ensembl97)

mapTASCINI = FALSE
if (mapTASCINI) {
# gene conversion function
dictionary = data.frame()
for (i in 1:length(Tascini_Genes)) {
  gene.query = Tascini_Genes[i]
  #print(gene.query)
  gene.ensembl = dictionary_Tascini[which(dictionary_Tascini$external_gene_name == gene.query),]$ensembl_gene_id
  gene.conversion = gtf_Guo_df_f[gtf_Guo_df_f$gene_id %in% gene.ensembl,]$gene_name
  #print(length(gene.conversion))
  if (length(gene.conversion) == 0) {
    dictionary.entry = data.frame(gene_Tascini = gene.query, gene_Guo = 'NA')
  } else if (length(gene.conversion) == 1) {
    dictionary.entry = data.frame(gene_Tascini = gene.query, gene_Guo = gene.conversion)
  } else {
    dictionary.entry = data.frame(gene_Tascini = gene.query,
                                  gene_Guo = paste(gene.conversion, collapse = ';'))
  }
  dictionary = rbind(dictionary, dictionary.entry)
}

dictionary_f = dictionary[which(dictionary$gene_Guo != 'NA'),]
nID = dictionary_f[which(as.character(dictionary_f$gene_Guo) != as.character(dictionary_f$gene_Tascini)),]


#lapply(Tascini_Genes, function(x) {
# gene.ensembl = dictionary_Tascini[which(dictionary_Tascini$external_gene_name == x),]$ensembl_gene_id
#gene.conversion = gtf_Guo_df_f[gtf_Guo_df_f$gene_id %in% gene.ensembl,]$gene_name
#paste(gene.conversion, collapse = ';')
#if (length(gene.conversion) == 0) {
# dictionary.entry = data.frame(gene_Tascini = x, gene_Guo = 'NA')
#} else if (length(gene.conversion) == 1) {
# dictionary.entry = data.frame(gene_Tascini = x, gene_Guo = gene.conversion)
#} else {
# dictionary.entry = data.frame(gene_Tascini = x,
#                    gene_Guo = paste(gene.conversion, collapse = ';'))
#}
#})


GeneGuo = row.names(UMImatrix_sparse_Guo)
.l = intersect(GeneGuo, nID$gene_Guo)
nID_use = nID[nID$gene_Guo %in% .l,]

GeneGuo_m = sapply(GeneGuo, function(x) {
  if (x %in% nID_use$gene_Guo) {
    as.character(nID_use[nID_use$gene_Guo == x,]$gene_Tascini)
  } else {
    as.character(x)
  }
})

GeneGuo_m = as.character(GeneGuo_m)
check = data.frame(new = GeneGuo_m, ex = GeneGuo)

venn(list(GUO = GeneGuo_m, Tasci = Tascini_Genes))
row.names(M) <- GeneGuo_m 

summary(UMImatrix_sparse_Guo['SELK',])
summary(M['SELENO',])

row.names(UMImatrix_sparse_Guo) <- GeneGuo_m 
}

mapGUO = TRUE
if (mapGUO) {
dictionary2 = data.frame()
for (i in 1:length(Guo_Genes)) {
  print(i)
  gene.query = Guo_Genes[i]
  gene.ensembl =  gtf_Guo_df_f[which(gtf_Guo_df_f$gene_name == gene.query), ]$gene_id
  gene.conversion = dictionary_Tascini[dictionary_Tascini$ensembl_gene_id %in% gene.ensembl,]$external_gene_name
  print(length(gene.conversion))
  if (length(gene.conversion) == 0) {
    dictionary.entry = data.frame(gene_Guo = gene.query, gene_Tascini = 'NA', N = length(gene.conversion))
  } else if (length(gene.conversion) == 1) {
    dictionary.entry = data.frame(gene_Guo = gene.query, gene_Tascini = gene.conversion, N = length(gene.conversion))
  } else {
    dictionary.entry = data.frame(gene_Guo = gene.query,
                                  gene_Tascini = paste(gene.conversion, collapse = ';'), 
                                  N = length(gene.conversion))
  }
  dictionary2 = rbind(dictionary2, dictionary.entry)
}

gene_Guo2 = dictionary2[dictionary2$N == 2,]
dictionary2_upd = dictionary2
dictionary2_upd$gene_Tascini = as.character(dictionary2_upd$gene_Tascini)
dictionary2_upd$gene_Tascini[dictionary2$N == 2] <- 
  as.character(dictionary2_upd$gene_Guo[dictionary2$N == 2])

dictionary2_upd_f = dictionary2_upd[which(dictionary2_upd$gene_Tascini != 'NA'),]
nID_2 = dictionary2_upd_f[which(as.character(dictionary2_upd_f$gene_Tascini) != as.character(dictionary2_upd_f$gene_Guo)),]


GeneGuo = row.names(UMImatrix_sparse_Guo)
.l = intersect(GeneGuo, nID_2$gene_Guo)
nID_use = nID_2[nID_2$gene_Guo %in% .l,]

GeneGuo_m = sapply(GeneGuo, function(x) {
  if (x %in% nID_use$gene_Guo) {
    as.character(nID_use[nID_use$gene_Guo == x,]$gene_Tascini)
  } else {
    as.character(x)
  }
})

GeneGuo_m = as.character(GeneGuo_m)
check = data.frame(new = GeneGuo_m, ex = GeneGuo)

venn(list(GUO = GeneGuo_m, Tasci = Tascini_Genes))
row.names(M) <- GeneGuo_m 

summary(UMImatrix_sparse_Guo['SELK',])
summary(M['SELENOK',])
summary(UMImatrix_sparse_Guo['SELM',])
summary(M['SELENOM',])

a = as.data.frame(table(GeneGuo_m))
d= as.character(a[a$Freq>1,]$GeneGuo_m)
check[check$new %in% d,]

row.names(UMImatrix_sparse_Guo) <- GeneGuo_m

}



pz.Guoetal =  CreateSeuratObject(counts = UMImatrix_sparse_Guo, 
                                 project = 'Guoetal')

#rescue = read.csv('/Users/tascini.annasofia/Downloads/hgnc-symbol-check-4.csv', skip = 1)
#table(rescue$Match.type)
#rescue = rescue[which(rescue$Match.type!= 'Unmatched'),]
#duplicate_check = as.data.frame(table(as.character(rescue$Input)), stringsAsFactors =F)
#repeated = duplicate_check$Var1[which(duplicate_check$Freq >1)]
#rescue = rescue[!(Guo_only$Input %in% rescue),]
#newID_Guo =rescue[which(as.character(rescue$Input) != as.character(rescue$Approved.symbol)),]
#new_names = intersect(rescue$Approved.symbol, Tascini_Genes)
#dict_rescue =  rescue[rescue$Approved.symbol %in% new_names,c("Input","Approved.symbol")]


other_conversion = FALSE
if (other_conversion) {
########### gene Symbols ##########
Gene_lit = row.names(pz.lit.somatic.best)
Gene_D1  = row.names(pz.204.raw)
Gene_D2  = row.names(pz.205.raw)
Gene_D3  = row.names(pz.198.raw)

library(venn)
venn(list(Guo=Gene_lit, D1 = Gene_D1, D2 = Gene_D2, D3 = Gene_D3))
list_intersection_2 = attr(x = b, "intersections")

venn(list(Guo=Gene_lit, Tascini = unique(c(Gene_D1, Gene_D2, Gene_D3))), intersections =F)
Guo_only = as.character(setdiff(Gene_lit, unique(c(Gene_D1, Gene_D2, Gene_D3))))

for (i in names(list_intersection)) {
  write.table(list_intersection[[i]], paste('Genes_',i,'.txt', sep=''), quote = F, row.names = F, col.names = F)
}

ensembl97 = useMart("ensembl",
                    dataset="hsapiens_gene_ensembl", 
                    host = 'http://jul2019.archive.ensembl.org')




GeneTascini = unique(c(Gene_D1, Gene_D2, Gene_D3))

dictionary_Tascini = getBM(attributes= c('ensembl_gene_id', 'external_gene_name'),
                           filter = 'external_gene_name',
                           #values = GeneTascini,
                           mart = ensembl97)



test =  gtf_Guo_df_f[gtf_Guo_df_f$gene_id %in% gene.ensembl,]
unique(gtf_Guo_df[gtf_Guo_df$gene_id %in% gene.ensembl,]$gene_name)


venn(list(Guo = dictionary_Guo$ensembl_gene_id, Tasci = dictionary_Tascini$ensembl_gene_id))
venn(list(Guo = dictionary_Guo$hgnc_symbol, Tasci = dictionary_Tascini$hgnc_symbol))


Guo_only = read.csv('/Users/tascini.annasofia/Downloads/hgnc-symbol-check.csv', skip = 1)
table(Guo_only$Match.type)
Guo_only = Guo_only[which(Guo_only$Match.type!= 'Unmatched'),]
duplicate_check = as.data.frame(table(as.character(Guo_only$Input)), stringsAsFactors =F)
repeated = duplicate_check$Var1[which(duplicate_check$Freq >1)]
Guo_only = Guo_only[!(Guo_only$Input %in% repeated),]

newID_Guo = Guo_only[which(as.character(Guo_only$Input) != as.character(Guo_only$Approved.symbol)),]
new_names = intersect(newID_Guo$Approved.symbol, list_intersection$Tascini)
dict_Guo =  Guo_only[Guo_only$Approved.symbol %in% new_names,c("Input","Approved.symbol")]


Guo_matrix = as.matrix(pz.Guoetal@assays$RNA@counts)

GeneGuo = row.names(Guo_matrix)
GeneGuo_m = replace(GeneGuo, which(GeneGuo %in% as.character(dict_Guo$Input)), as.character(dict_Guo$Approved.symbol))

venn(list(GUO = GeneGuo_m, Tasci = GeneTascini))
row.names(Guo_matrix) <- GeneGuo_m 
Guo_matrix_sparce = Matrix(Guo_matrix, sparse = T)

pz.Guoetal <- CreateSeuratObject(counts = Guo_matrix_sparce, 
                         project = 'Guoetal')

}


########### integration list ################
load(paste(data.path, 'pz.lit.somatic.best', sep=''))
load(paste(data.path, 'pz.198.raw', sep=''))
load(paste(data.path,'pz.205.raw', sep =''))
load(paste(data.path,'pz.204sx.raw.polished', sep =''))

# select somatic cells
whichcell = row.names(pz.lit.somatic.best[[]])
pz.Guoetal.somatic = subset(pz.Guoetal, cells = whichcell)

# obj name updating
pz.204.raw <- pz.204sx.raw.polished

# metadata samples
pz.198.raw$sample <- 'pz198'
pz.204.raw$sample <- 'pz204'
pz.205.raw$sample <- 'pz205'
pz.Guoetal.somatic$sample <- 'Guoetal'

pz.198.raw <- RenameCells(pz.198.raw, add.cell.id = 'pz198_')
pz.204.raw <- RenameCells(pz.204.raw, add.cell.id = 'pz204_')
pz.205.raw <- RenameCells(pz.205.raw, add.cell.id = 'pz205_')
pz.Guoetal.somatic <- RenameCells(pz.Guoetal.somatic, add.cell.id = 'Guoetal_')


# metadata condition
pz.198.raw$condition <- 'iNOA'
pz.204.raw$condition <- 'iNOA'
pz.205.raw$condition <- 'iNOA'
pz.Guoetal.somatic$condition<- 'CTL'

# list creation
somatic.list <- list()
somatic.list <- list(pz.198.raw, pz.204.raw, pz.205.raw, pz.Guoetal.somatic)
names(somatic.list) <- c('pz198', 'pz204','pz205',"Guoetal")
print("input_samples:")
print(somatic.list)


############  integration ############ 
#Use the same filtering for all the patients

min_nFeature_RNA = 500
max_nFeature_RNA = 6500
max_percent_MT = 20

for (i in 1:length(somatic.list)) {
    
    somatic.list[[i]][["percent.mt"]] <- PercentageFeatureSet(somatic.list[[i]], 
                                                              pattern = "^MT-")
    somatic.list[[i]] <-  subset(somatic.list[[i]],
                                     subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)
    somatic.list[[i]] <- NormalizeData(somatic.list[[i]], 
                                           verbose = FALSE)
    somatic.list[[i]] <- FindVariableFeatures(somatic.list[[i]], 
                                                  selection.method = 'vst',
                                                  nfeatures = 2000,
                                                  verbose = FALSE)
}
print("samples after QC filters")
print(somatic.list)

# integration
somatic.anchors <- FindIntegrationAnchors(object.list = somatic.list, dims = 1:30, verbose = FALSE)
somatic.integrated <- IntegrateData(anchorset = somatic.anchors, dims = 1:30, verbose = FALSE)
# switch to integrated assay. 
# The variable features of this assay are automatically set during IntegrateData

DefaultAssay(somatic.integrated) <- "RNA"
all.genes <- row.names(somatic.integrated)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
somatic.integrated <- CellCycleScoring(somatic.integrated, 
                                       s.features = s.genes, 
                                       g2m.features = g2m.genes, 
                                       set.ident = TRUE,
                                       verbose=F)
DefaultAssay(somatic.integrated) <- "integrated"

# SCALE and regress data
somatic.integrated <- ScaleData(somatic.integrated, 
                                    verbose = FALSE, 
                                    vars.to.regress = c("percent.mt", "nFeature_RNA"), 
                                    features = all.genes)
# Regress for Cell Cycle
somatic.integrated <- ScaleData(somatic.integrated, 
                                verbose = FALSE,
                                vars.to.regress = c("S.Score", "G2M.Score"), 
                                features = all.genes) 

somatic.integrated <- RunPCA(somatic.integrated, npcs = 30, verbose = FALSE)
ElbowPlot(somatic.integrated)

somatic.integrated <- RunUMAP(somatic.integrated, reduction = "pca", dims = 1:20, verbose = FALSE)

DimPlot(somatic.integrated, reduction = "umap", pt.size = 2)
DimPlot(somatic.integrated, reduction = "umap", pt.size = 3, split.by = 'sample', group.by = 'sample')


res = 0.2
nPC = 20
somatic.integrated <- FindNeighbors(somatic.integrated, dims = 1:nPC, verbose = FALSE)
somatic.integrated <- FindClusters(somatic.integrated, resolution = res, verbose = FALSE)
somatic.integrated <- RunUMAP(somatic.integrated, dims = 1:nPC, verbose = FALSE)
somatic.integrated <- RunTSNE(somatic.integrated, dims = 1:nPC, verbose = FALSE)

p1 = DimPlot(somatic.integrated, reduction = "umap", label = T, pt.size = 2) + 
  theme(plot.title = element_text(color="blue", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = 'dodgerblue4', size=22, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "dodgerblue2", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'dodgerblue4', size=22),
        axis.title.y = element_text(face = "bold", color = "dodgerblue2", size = 24),
        legend.text = element_text(face = "bold", color = "dodgerblue2", size = 22),
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(subtitle = paste('Res = ', res,', nPC = ',nPC, sep = ''), x = "UMAP 1", y = "UMAP 2") 
p2 = DimPlot(somatic.integrated, reduction = "umap", group.by = 'condition',label = F, pt.size = 2) + 
  theme(plot.title = element_text(color="blue", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = 'dodgerblue4', size=22, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "dodgerblue2", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'dodgerblue4', size=22),
        axis.title.y = element_text(face = "bold", color = "dodgerblue2", size = 24),
        legend.text = element_text(face = "bold", color = "dodgerblue2", size = 22),
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(subtitle = paste('Res = ', res,', nPC = ',nPC, sep = ''), x = "UMAP 1", y = "UMAP 2") 

p1 | p2


DimPlot(somatic.integrated, reduction = "umap", label = T, pt.size = 2, split.by = 'sample') + 
  theme(plot.title = element_text(color="blue", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = 'dodgerblue4', size=22, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "dodgerblue2", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'dodgerblue4', size=22),
        axis.title.y = element_text(face = "bold", color = "dodgerblue2", size = 24),
        legend.text = element_text(face = "bold", color = "dodgerblue2", size = 22),
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(subtitle = paste('Res = ', 0.2,', nPC = ',20, sep = ''), x = "UMAP 1", y = "UMAP 2") 

DimPlot(somatic.integrated, split.by = 'condition', label = T)

######### Tcell REassignemnt ###############
obj = somatic.integrated
#obj = somatic.integrated.new2
for (g in  c('CD3E','TRAC', 'GZMK')) {
  v = VlnPlot(subset(obj, idents = c('2','6')), features =g, 
        split.by = 'condition', assay = 'RNA', split.plot = F, pt.size = 0.5) + 
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = c(0.5, 0.92),
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  #geom_boxplot(width = 0.1, outlier.size=1) +
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTL'))
  assign(paste('v_', g, sep =''),v)
}

v_CD3E | v_TRAC | v_GZMK

for (g in  c('HLA-DRA','CD74', 'HLA-DPA1')) {
  v = VlnPlot(subset(obj, idents = c('2','6')), features =g, 
              split.by = 'condition', assay = 'RNA', pt.size = .5) + 
    theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 45, face = "bold", color = "black", size = 16, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.5, 0.92),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    #geom_boxplot(width = 0.1, outlier.size=1) +
    scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTL'))
  assign(paste('v_', g, sep =''),v)
}

`v_HLA-DRA` | v_CD74 | `v_HLA-DPA1`

pdf(paste(outdir,'Figure_S15.pdf', sep=''), width = 18, height = 12)
print((v_CD3E | v_TRAC | v_GZMK) / (`v_HLA-DRA` | v_CD74 | `v_HLA-DPA1`))
dev.off()

VlnPlot(somatic.integrated, features = c('HLA-DRA','CD74', 'HLA-DPA1'), split.by = 'condition', assay = 'RNA')

#table(somatic.integrated@meta.data[somatic.integrated@meta.data$integrated_snn_res.0.1=='6',]$condition)

# assign 14 wrongly assigned TCL (cluster 6) to macrophages (cluster 2) 
cell2rename <- rownames(somatic.integrated@meta.data[somatic.integrated@meta.data$condition == 'CTL' & 
                                                       somatic.integrated@meta.data$seurat_clusters == '6',])

somatic.integrated.mod <- SetIdent(object = somatic.integrated, cells = cell2rename, value = '2')

seeMG = F
if (seeMG) {
  mg = FindAllMarkers(somatic.integrated.mod)
  cluster.markers = mg
  top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  DoHeatmap(somatic.integrated.mod, features = top10$gene, angle = 90) + NoLegend()
}

DimPlot(somatic.integrated.mod, label = T)

#### clusterID
#Assign cluster to cell types 
new.cluster.ids <- c( "MCR",
                      "LEY", 
                     "MYD",
                     "END", 
                     "UND",
                     "SRT",
                     "TCL",
                     "STRO",
                     "STRO",'END')

names(new.cluster.ids) <- levels(somatic.integrated.mod)
somatic.integrated.new <- RenameIdents(somatic.integrated.mod, new.cluster.ids)
DimPlot(somatic.integrated.new, reduction = "umap", label = TRUE, pt.size = 2)

somatic.integrated.new@meta.data$cell_type=Idents(somatic.integrated.new)
somatic.integrated.new$cell_type <- somatic.integrated.new@meta.data$cell_type

names(new.cluster.ids) <- levels(somatic.integrated.mod)
somatic.integrated.new2 <- RenameIdents(somatic.integrated.mod, new.cluster.ids)
DimPlot(somatic.integrated.new2, reduction = "umap", label = TRUE, pt.size = 2)

somatic.integrated.new2@meta.data$cell_type=Idents(somatic.integrated.new2)
somatic.integrated.new2$cell_type <- somatic.integrated.new2@meta.data$cell_type

somatic.integrated.new$condition <- factor(x = somatic.integrated.new$condition, levels = c("iNOA", "CTL"))
somatic.integrated.new2$condition <- factor(x = somatic.integrated.new2$condition, levels = c("iNOA", "CTL"))

order.2.plot <- c("LEY",
                 "MYD",
                 "SRT",
                 "MCR",
                 "TCL",
                 "END",
                 "STRO",
                 "UND")
levels(somatic.integrated.new) <- order.2.plot 
levels(somatic.integrated.new2) <- order.2.plot
DimPlot(somatic.integrated.new, label = T, split.by = 'condition' )

######### save image ##########
save.image("Somatic_integration_beforeplots.RData")

OBJ_dir = '/Users/tascini.annasofia/Dropbox (HSR Global)/WORKSPACE/Alfano/904_infertilita_epigenetics/7_bioinfo/SeuratObjects/'
dir.create(OBJ_dir, recursive = T)

DimPlot(pz.literature.new)
SO_Guoetal <-  pz.literature.new

head(x = SO_Guoetal[[]])
# Create a new metadata column called "original" with the contents of "orig.ident"
SO_Guoetal$orig.ident <- NULL
SO_Guoetal$old.ident <- NULL
SO_Guoetal$source <- "GSE112013"
SO_Guoetal$old.ident  <- "Guoetal"
SO_Guoetal$orig.ident  <- "Guoetal"
SO_Guoetal$condition <- "CTL"
# Quickly preview the metadata again
head(x = SO_Guoetal[[]])

#
pz.179.update
OA_donor <- pz.179.update
head(x = OA_donor[[]])
OA_donor$source <- "OA donor"
OA_donor$orig.ident <- NULL
OA_donor$old.ident <- NULL
OA_donor$orig.ident <- "OA pz.179"
OA_donor$old.ident <- "OA"
OA_donor$condition <- "OA"
head(x = OA_donor[[]])

#
pz.198
iNOA_donor3 <- pz.198
head(x = iNOA_donor3[[]])
iNOA_donor3$source <- "iNOA donor 3"
iNOA_donor3$orig.ident <- NULL
iNOA_donor3$orig.ident <- "iNOA pz.198"
iNOA_donor3$condition <- "iNOA"
head(x = iNOA_donor3[[]])

#
pz.204
iNOA_donor1 <- pz.204
head(x = iNOA_donor1[[]])
iNOA_donor1$source <- "iNOA donor 1"
iNOA_donor1$orig.ident <- NULL
iNOA_donor1$orig.ident <- "iNOA pz.204"
iNOA_donor1$condition <- "iNOA"
head(x = iNOA_donor1[[]])

#
pz.205
iNOA_donor2 <- pz.205
head(x = iNOA_donor2[[]])
iNOA_donor2$source <- "iNOA donor 2"
iNOA_donor2$orig.ident <- NULL
iNOA_donor2$orig.ident <- "iNOA pz.205"
iNOA_donor2$condition <- "iNOA"
head(x = iNOA_donor2[[]])

somatic.integrated <- somatic.integrated.new
head(somatic.integrated[[]])
somatic.integrated$orig.ident <- NULL
somatic.integrated$old.ident <- NULL
head(x = somatic.integrated[[]])


save(SO_Guoetal, file = paste(OBJ_dir,'SO_Guoetal', sep =''))
save(OA_donor, file = paste(OBJ_dir,'OA_donor', sep =''))
save(iNOA_donor3, file = paste(OBJ_dir,'iNOA_donor3', sep =''))
save(iNOA_donor2, file = paste(OBJ_dir,'iNOA_donor2', sep =''))
save(iNOA_donor1, file = paste(OBJ_dir,'iNOA_donor1', sep =''))
save(somatic.integrated, file = paste(OBJ_dir,'somatic.integrated', sep =''))

save.image(paste(OBJ_dir,"Somatic_integration_beforeplots.RData", sep=''))

######  plots ######rm#
# HM
DefaultAssay(somatic.integrated.new2) <- 'integrated'
HM_azo_10 <- SaveMarkers(Seurat.object = somatic.integrated.new2, 
                      filename_xlsx = paste(outdir,'somatic.intergrated.MK_logFCpos.xlsx',sep=''),
                      LogFC.onlypos = TRUE)



Seurat.object = somatic.integrated.new2
Seurat.object.downsample <- subset(Seurat.object, downsample = 100)

col = c('#E41A1C', #L
        '#3399FF', #M
        '#FFCC33', #S
        '#9900FF', #M
        '#FF33CC', #Tcell
        '#4DAF4A', #E
        'cadetblue2', #U
        'grey83', #U
        '#999999') #U
#Seurat.object.downsample
HM <- DoHeatmap(Seurat.object, 
                features = HM_azo_10$gene, 
                group.colors = col,
                disp.min = -2,
                disp.max = 2,
                angle = 90) +
scale_fill_gradientn(colours = coolwarm(200)) 

HM

pdf(paste(outdir,"somatic.HM.UPD_allcell.pdf",sep=''),  width=14, height=14)
plot(HM)
dev.off()

# UMAP for figure
col = as.character(colors_Alf[levels(somatic.integrated.new2),]$cols)
p.inf <- DimPlot(somatic.integrated.new2, 
                 reduction = "umap", 
                 label = T, 
                 order = T,
                 pt.size = 2, 
                 split.by = 'condition') +
scale_color_manual(values = (col)) +
theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
      axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=22, hjust =1), 
      axis.title.x = element_text(face = "bold", color = "black", size = 24),
      axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
      axis.title.y = element_text(face = "bold", color = "black", size = 24),
      legend.text = element_text(face = "bold", color = "black", size = 12),
      legend.position="top",
      panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
labs(x = "UMAP 1", y = "UMAP 2") + NoLegend()

p.inf

pdf(paste(outdir,"INF.umap.pdf",sep=''),  width=14, height=6)
p.inf
dev.off()


Plot_sign <- function(Seraut.object, signature, operator = sum, title = '') {
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
                      repel = T,
                      cols = c("lightgrey", "red")) +
                      #cols = as.vector(coolwarm(5))) +
    theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=9, face="italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = 'black', size=12, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 14),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=12),
          axis.title.y = element_text(face = "bold", color = "black", size = 14),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(title = title, 
         subtitle = paste('MarkerGenes: ',toString(signature), sep=''), 
         x = "UMAP 1", y = "UMAP 2") 
    return(FP)
    }


Seurat.object = somatic.integrated.new

leydig_signature = c('CFD','DLK1','LUM','CALB2')
myoid_signature = c('ACTA2','MYH11','DES', 'MYL9')
sertoli_signature = c('FATE1','CITED1','SOX9','AMH','CLDN11')
macrophage_signature = c('CD14','CD74','HLA-DRA','HLA-DRB1')
endothelial_signature = c('VWF', "EGFL7",'CD34', 'PRSS23', "RBP7")
t_cell_signature =  c('GZMA','GZMK','CD2','CCL5','NKG7')

DefaultAssay(Seurat.object) <- 'RNA'
pL <- Plot_sign(Seurat.object,
          signature= c('CFD','DLK1','LUM'), 
          operator = mean, title = 'LEY')
pM <- Plot_sign(Seurat.object,
          signature= c('ACTA2','MYH11','DES'), 
          operator = mean, title = 'MYD')
pS <- Plot_sign(Seurat.object,
          signature= c('FATE1','CITED1','SOX9'), 
          operator = mean, title = 'SRT')
pMa <- Plot_sign(Seurat.object,
          signature= c('CD14','CD74','HLA-DRA'), 
          operator = mean, title = 'MCR')
pE <- Plot_sign(Seurat.object,
          signature= c('VWF',"EGFL7", 'PRSS23'), 
          operator = mean, title = 'END')
pT <- Plot_sign(Seurat.object,
          signature= c('GZMA','CD2','CCL5'), 
          operator = mean, title = 'TCL')
pSTRO <- Plot_sign(Seurat.object,
          signature= c('RGS5','TPM2','IGFBP5'), 
          operator = mean, title = 'STRO')
layout <- '
AABBCC##
AABBCC##
AABBCCDD
EEFFGGDD
EEFFGG##
EEFFGG##
'
 
wrap_plots(A = pL, B = pM, C = pS, G = pT, E = pE, F = pMa, D =pSTRO, design = layout)

pdf(paste(outdir,"MG_somaticintegration.pdf",sep=''),  width=14, height=6)
wrap_plots(A = pL, B = pM, C = pS, D = pSTRO, E = pE, F = pMa, G =pT, design = layout)
dev.off()


library(magick)

layout <- '
UUUUUUAABBCC##
UUUUUUAABBCC##
UUUUUUAABBCCDD
UUUUUUEEFFGGDD
UUUUUUEEFFGG##
UUUUUUEEFFGG##
'
wrap_plots(U = p.inf, 
           A = pL, B = pM, C = pS, D = pSTRO, E = pE, F = pMa, G =pT, design = layout) 
#+  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 8))

layout <- '
UAABBCC##
UAABBCC##
UAABBCCDD
UEEFFGGDD
UEEFFGG##
UEEFFGG##
'
  
#+plot_annotation(title = 'The surprising story about mtcars')
wrap_plots(U = wrap_elements(grid::textGrob('Integration \n somatic \n  cells \n iNOA + CTL')), 
           A = pL, B = pM, C = pS, D = pT, E = pE, F = pMa, G =pSTRO, design = layout)  

pdf(paste(outdir,"UMAP_MG_somaticintegratino.pdf",sep=''),  width=20, height=8)
wrap_plots(U = p.inf, 
           A = pL, B = pM, C = pS, D = pSTRO, E = pE, F = pMa, G =pT, design = layout) 
dev.off()

im = image_read_pdf(paste(outdir,"Figure2_panel_E.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure2/Figure2_panel_E.tiff", sep=''), format = "tiff")

im = image_read_pdf(paste(outdir,"Figure2/Figure2_panel_B.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure2/Figure2_panel_B.tiff", sep=''), format = "tiff")
im = image_read_pdf(paste(outdir,"Figure2/Figure2_panel_D.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure2/Figure2_panel_D.tiff", sep=''), format = "tiff")

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
                    split.by = 'condition',
                    cols = c("lightgrey", "red"))
  #cols = as.vector(coolwarm(5))) + 
  return(FP)
}

Seurat.object = somatic.integrated.new2

leydig_signature = c('CFD','DLK1','LUM','CALB2')
myoid_signature = c('ACTA2','MYH11','DES', 'MYL9')
sertoli_signature = c('FATE1','CITED1','SOX9','AMH')
macrophage_signature = c('CD14','CD74','HLA-DRA','HLA-DRB1')
endothelial_signature = c('VWF', "EGFL7",'CD34', 'PRSS23', "RBP7")
t_cell_signature =  c('GZMA','GZMK','TRAC','CD3E')

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


(pL / pM / pS ) | (pMa / pT / pE)

pdf(paste(outdir,"somatic_MG.pdf",sep=''),  width=14, height=9)
(pL / pM / pS ) | (pMa / pT / pE)
dev.off()

im = image_read_pdf(paste(outdir,"somatic_MG.pdf",sep=''),density = 200)
image_write(im, path = paste(outdir,"somatic_MG.tiff", sep=''), format = "tiff")

pdf(paste(outdir,"somatic_MG_stro.pdf",sep=''),  width=7, height=3)
pSTRO
dev.off()
im = image_read_pdf(paste(outdir,"somatic_MG_stro.pdf",sep=''),density = 200)
image_write(im, path = paste(outdir,"somatic_MG_stro.tiff", sep=''), format = "tiff")

########## DGE #############
# relabel cells to perform DGE
DefaultAssay(somatic.integrated.new) <- "RNA"
somatic.integrated.new$celltype.cond <- paste(Idents(somatic.integrated.new), 
                                                  somatic.integrated.new$condition, sep = "_")
somatic.integrated.new$celltype <- Idents(somatic.integrated.new)

Idents(somatic.integrated.new) <- "celltype.cond"

myDGE(somatic.integrated.new,'LEY', outdir = outdir)
myDGE(somatic.integrated.new,'MYD', outdir = outdir)
myDGE(somatic.integrated.new,'SRT', outdir = outdir)
myDGE(somatic.integrated.new,'MCR', outdir = outdir)
myDGE(somatic.integrated.new,'END', outdir = outdir)
myDGE(somatic.integrated.new,'STRO', outdir = outdir)
myDGE(somatic.integrated.new,'UND', outdir = outdir)



suppressMessages(library(enrichR))

dir_DGE = paste(outdir,'excel_tables/DGE_tables/', sep='')
#'/Users/tascini.annasofia/Dropbox (HSR Global)/Alfano_904_paperDraft/xlsx.table/somatics_DGE/'
dir_enrichR = paste(outdir,'enrichR/', sep='')
dir.create(dir_enrichR)

#'/Users/tascini.annasofia/Dropbox (HSR Global)/Alfano_904_paperDraft/xlsx.table/somatics_DGE_enrichR/'

databases <- listEnrichrDbs()
# databases to make the enrichment of
enrich.databases <- c("GO_Biological_Process_2018",
                      "GO_Cellular_Component_2018",
                      "GO_Molecular_Function_2018",
                      "Reactome_2016",
                      "KEGG_2016",
                      "WikiPathways_2016",
                      "BioCarta_2016",
                      "BioPlanet_2019")


cell_types= c("LEY", "MYD", "SRT", "MCR", "END", "STRO", "UND")
DGE_gene_up_list <- list()
DGE_gene_dw_list <- list()

for (celltype in cell_types) {
  print(celltype)
  DGE_file = paste(dir_DGE,"DGE_", celltype, ".iNOA.vs.CTRL.xlsx", sep= '')
  A <- read.xlsx(DGE_file)
  
  up.genes <- A$row.names[A$avg_logFC>0 & A$p_val_adj <=0.05]
  down.genes <- A$row.names[A$avg_logFC<0 & A$p_val_adj <=0.05]
  both.genes <- A$row.names[A$p_val_adj <=0.05]
  
  DGE_gene_up_list[[celltype]] <- up.genes
  DGE_gene_dw_list[[celltype]] <- down.genes
  
  enrichr.list <- list()
  enrichr.list <- lapply(list(up.genes,down.genes,both.genes),function(x) {
    enrichR::enrichr(genes = x, databases = enrich.databases)
  })
  
  names(enrichr.list) <-  c("up","down","both")
  
  for (i in names(enrichr.list)){
    for (dat in names(enrichr.list[[i]])) {
      enrichr.list[[i]][[dat]] = enrichr.list[[i]][[dat]][,
                                                          c('Term', 
                                                            'Overlap', 
                                                            'P.value', 
                                                            'Adjusted.P.value', 
                                                            'Odds.Ratio', 
                                                            'Combined.Score', 
                                                            'Genes')]
    }
  }
  for (i in 1:length(enrichr.list)) {
    filename = paste(dir_enrichR,'enrichR_DGE_',celltype,'_',names(enrichr.list)[i],'.xlsx', sep='')
    write.xlsx(x = enrichr.list[[i]], file = filename, asTable = T)
  }
}


library(venn)
venn::venn(DGE_gene_up_list, zcolor = rainbow(7))
venn::venn(DGE_gene_dw_list, zcolor = rainbow(7))

#Reduce(intersect, DGE_gene_up_list)
#Reduce(intersect, DGE_gene_dw_list)



PD="/Users/tascini.annasofia/Dropbox (HSR Global)/Alfano_904_paperDraft/"
DGEdir=paste(PD,"xlsx.table/somatics_DGE/", sep='')
AGEdir=paste(PD,"41467_2015_BFncomms9570/", sep='')
# age related gene from reference 41467_2015_BFncomms9570
allAgeGenesFile = paste(AGEdir,'41467_2015_BFncomms9570_MOESM436_ESM.xlsx',sep='')
allAgeGenesTable = read.xlsx(allAgeGenesFile, startRow = 3)
allAgeGenesID_Direction = allAgeGenesTable[1:1497,c('NEW-Gene-ID','Direction')]
allAgeGenesID_up = allAgeGenesID_Direction$`NEW-Gene-ID`[allAgeGenesID_Direction$Direction=='+']
allAgeGenesID_dw = allAgeGenesID_Direction$`NEW-Gene-ID`[allAgeGenesID_Direction$Direction=='-']
# consider also Zscore?

AgeGenesClustersFile = paste(AGEdir,'41467_2015_BFncomms9570_MOESM438_ESM.xlsx',sep='')
AgeGenesClustersTable_dw = read.xlsx(AgeGenesClustersFile, startRow = 4, sheet = 'A')
row.names(AgeGenesClustersTable_dw) <- AgeGenesClustersTable_dw$gene.name
AgeGenesClustersTable_up = read.xlsx(AgeGenesClustersFile, startRow = 2, sheet = 'N')
row.names(AgeGenesClustersTable_up) <- AgeGenesClustersTable_up$gene.name

cell_types= c("LEY", "MYD", "SRT", "MCR", "TCL","END", "STRO", "UND")
#cell_types = cell_types[4]
DGE_gene_list <- list()
for (celltype in cell_types) {
    #print(celltype)
    DGE_file = paste(DGEdir,"DGE_", celltype, ".iNOA.vs.CTRL.xlsx", sep= '')
    A <- read.xlsx(DGE_file)
    DGE_gene_list_tmp <- list()
    up.genes <- A$row.names[A$avg_logFC>0 & A$p_val_adj <0.05]
    down.genes <- A$row.names[A$avg_logFC<0 & A$p_val_adj <0.05]
    both.genes <- A$row.names[A$p_val_adj <0.05]    
    DGE_gene_list_tmp <- lapply(list(up.genes,down.genes,both.genes),function(x) {x})
    names(DGE_gene_list_tmp) <-  c("up","down","both")
    
    int_both = intersect(c(allAgeGenesID_dw,allAgeGenesID_up), DGE_gene_list_tmp[['both']])
    int_DW = intersect(allAgeGenesID_dw, DGE_gene_list_tmp[['down']])
    int_W_DW = intersect(allAgeGenesID_up, DGE_gene_list_tmp[['down']])
    int_UP = intersect(allAgeGenesID_up, DGE_gene_list_tmp[['up']])
    int_W_UP = intersect(allAgeGenesID_dw, DGE_gene_list_tmp[['up']])
    #print(c(length(int_DW), length(int_W_DW), length(DGE_gene_list_tmp[['down']]), length(allAgeGenesID_dw)))
    #print(c(length(int_UP), length(int_W_UP), length(DGE_gene_list_tmp[['up']]), length(allAgeGenesID_up)))
    #print(c(length(int_both),length(DGE_gene_list_tmp[['both']]), length(allAgeGenesID_up)+length(allAgeGenesID_dw)))
    
    DGE_gene_list[[celltype]] <- DGE_gene_list_tmp
    
    }






AgeGenesClustersTable_dw$gene.name[AgeGenesClustersTable_dw$cluster==3]

table(factor(Idents(somatic.integrated.new)))
#Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("pDC", "Eryth", "Mk", "DC", 
 #   "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))


#### Age Dotplot ####
Seurat_object <- somatic.integrated.new
order.2.plot = c('LEY_iNOA', 'LEY_CTL', 'MYD_iNOA','MYD_CTL','SRT_iNOA','SRT_CTL',
              'END_iNOA','END_CTL','MCR_iNOA','MCR_CTL','TCL_iNOA',
             'STRO_iNOA','STRO_CTL','UND_iNOA','UND_CTL')

levels(Seurat_object) <- order.2.plot

DefaultAssay(Seurat_object) <- "RNA"
markers.to.plot <- AgeGenesClustersTable_dw$gene.name[AgeGenesClustersTable_dw$cluster==3]
dp_RPL = DotPlot(Seurat_object, 
        col.min = -1,
        col.max = 1,
        features = rev(markers.to.plot), 
        cols = c("lightblue", "red"), 
        dot.scale = 8) +
        #split.by = "condition") + 
RotatedAxis() +
#coord_flip() +
theme(axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=20, vjust =.5, hjust =1, 
                                 family = 'mono'), 
      axis.title.x = element_text(face = "bold", color = "black", size = 24),
      axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=24),
      axis.title.y = element_text(face = "bold", color = "black", size = 24),
      legend.text = element_text(face = "bold", color = "black", size = 16),
      legend.position="right",
      panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
labs(x = "Genes", y = "Cell type")

dp_RPL
jpeg(paste(outdir,'DotPlot_RPLgenes.jpg',sep=''),width=1400, height=700, unit='px')
plot(dp_RPL)
dev.off()

pdf(paste(outdir,'DotPlot_RPLgenes.pdf',sep=''),width=14, height=7)
plot(dp_RPL)
dev.off()

im = image_read_pdf(paste(outdir,"DotPlot_RPLgenes.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure7_DotPlot_RPLgene.tiff", sep=''), format = "tiff")

#### ViolinPLots ####

legendiNOA='Germ Cell aplasia'
legendCTRL='Normal Spermatogenesis'

options(repr.plot.width=18, repr.plot.height=5)
Seurat.object = somatic.integrated.new2
# insulin-like
#genes = c('AR','IGF1','IGF1R','IGF2','IGF2R')
#genes = c('COL1A1','COL1A2','COL4A1','COL4A2','COL4A3','COL4A4','COL4A5','COL4A6')
genes = c('DLK1','NOTCH2','NOTCH3','HSD17B3')
#genes = c('MAGEB4', 'NANOS1', 'NR0B1','NR5A1','SOHLH1','SYCE1','TAF4B','TEX11','TEX15','WT1','ZMYND15')
#N = 1
for (N in 1:length(genes)) {
gene=genes[N]
DefaultAssay(Seurat.object) <- "RNA"
#Seurat.object <- subset(Seurat.object,  idents = c("LEY", "MYD","SRT","END"))
vp <- VlnPlot(Seurat.object,
        feature = gene,
        slot = "counts", 
        log = TRUE,
        pt.size = 0,
        split.plot = T,
       split.by = 'condition') +
ylab("Log Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.85, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c(legendiNOA,legendCTRL)) 
 
jpeg(paste(outdir,'VlnPlot.allCell.',gene,'.jpg',sep=''),width=1500, height=500, unit='px')
plot(vp)
dev.off()
    }

for (N in 1:length(genes)) {
gene=genes[N]
DefaultAssay(Seurat.object) <- "RNA"
Seurat.object <- subset(Seurat.object,  idents = c("LEY", "MYD","SRT"))
Seurat.object <- subset(Seurat.object,  idents = c("LEY", "MYD","SRT"))#,"END"))
vp <- VlnPlot(Seurat.object,
        feature = gene,
        slot = "counts", 
        log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.85, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c(legendiNOA,legendCTRL)) 
 
jpeg(paste(outdir,'VlnPlot.LMSE.',gene,'.jpg',sep=''),width=1200, height=550, unit='px')
plot(vp)
dev.off()
    }

paste(outdir,'VlnPlot.allCell.',gene,'.jpg')

# Receptor-ligand pairs

Seurat.object = somatic.integrated.new2
DefaultAssay(Seurat.object) <- "RNA"
#Seurat.object_s <- subset(Seurat.object,  idents = c("LEY", "MYD","SRT"))
Seurat.object_s <- somatic.integrated.new2
VlnPlot(Seurat.object,
        feature = c('DLK1'),
        slot = "counts", 
        log = TRUE,
        pt.size = 0,
        split.plot = T,
       split.by = 'condition') +
#ylab("Log Expression Level") +
#xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.75, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTRL')) 

#stat_summary(fun.data=c("mean_sdl"),  fun.args = list(mult=1), geom="pointrange", color = c("red"))

VlnPlot(Seurat.object,
        feature = c('NOTCH2'),
        slot = "counts", 
        log = TRUE,
        pt.size = 0,
       split.by = 'condition',
       split.plot = T,) +
ylab("Log Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5),
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.75, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTRL')) 



Seurat.object = somatic.integrated.new2
DefaultAssay(Seurat.object) <- "RNA"

#genes = c('AR','IGF1','IGF1R','IGF2','IGF2R')
#genes = c('COL1A1','COL1A2','COL4A1','COL4A2','COL4A3','COL4A4','COL4A5','COL4A6')
genes = c('DLK1','NOTCH2','NOTCH3','HSD17B3')

#genes = c('AR','IGF1','IGF1R','IGF2','IGF2R',
 #         'COL1A1','COL1A2','COL4A1','COL4A2','COL4A3','COL4A4','COL4A5','COL4A6',
  #        'DLK1','NOTCH2','NOTCH3','HSD17B3')
#genes = c('MAGEB4', 'NANOS1', 'NR0B1','NR5A1','SOHLH1','SYCE1','TAF4B','TEX11','TEX15','WT1','ZMYND15')
#N = 1
for (N in 1:length(genes)) {
    gene = genes[N]
    fp <- FeaturePlot(Seurat.object, 
            reduction = "umap", 
            pt.size = 1,
            features = gene,
            max.cutoff = 'q99',
            label = T, 
            order=T,
            cols = c("lightgrey", "red"),
            split.by='condition') + theme(legend.position = 'right')

    jpeg(paste(outdir,'FeaturePlot.',gene,'.jpg',sep=''),width=1000, height=500, unit='px')
    plot(fp)
    dev.off()
}



# Receptor-ligand pairs
options(repr.plot.width=11, repr.plot.height=5)
Seurat.object = somatic.integrated.new2
DefaultAssay(Seurat.object) <- "RNA"
#Seurat.object <- subset(Seurat.object,  idents = c("LEY", "MYD","END","SRT"))

VlnPlot(Seurat.object,
        feature = c('FCGR2A'),
        slot = "counts", 
        log = TRUE,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.75, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c(legendiNOA,legendCTRL)) 

#stat_summary(fun.data=c("mean_sdl"),  fun.args = list(mult=1), geom="pointrange", color = c("red"))

VlnPlot(Seurat.object,
        feature = c('IGF2'),
        slot = "counts", 
        log = TRUE,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5),
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.75, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c(legendiNOA,legendCTRL)) 



VlnPlot(Seurat.object,
        feature = c('IGF1R'),
        slot = "counts", 
        log = TRUE,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5),
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.75, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c(legendiNOA,legendCTRL)) 



VlnPlot(Seurat.object,
        feature = c('IGF2R'),
        slot = "counts", 
        log = TRUE,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5),
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.75, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c(legendiNOA,legendCTRL)) 




FeaturePlot(Seurat.object, 
            reduction = "umap", 
            pt.size = 2,
            features = c('TCF21'), 
            label = T, 
            order= T,
            split.by = 'condition',
            max.cutoff = 'q99',
            cols = c("lightgrey", "red")) + theme(legend.position = 'right')

FeaturePlot(Seurat.object, 
            reduction = "umap", 
            pt.size = 2,
            features = c('ACE2'), 
            label = T, 
            order= T,
            max.cutoff = 'q99',
            cols = c("lightgrey", "red"),
            split.by='condition') + theme(legend.position = 'right')


##### DLK1 NOTCH2 #####
#DLK1
Seurat.object = somatic.integrated.new2
DefaultAssay(Seurat.object) <- "RNA"
Seurat.object_s <- subset(Seurat.object,  idents = c("LEY", "MYD","SRT"))
#Seurat.object_s = Seurat.object
vp_D <- VlnPlot(Seurat.object_s,
        feature = 'DLK1',
        #slot = "counts", 
        #log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.85, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))

fp_D <- FeaturePlot(Seurat.object, 
            reduction = "umap", 
            pt.size = 1,
            features = c('DLK1'), 
            label = T, 
            order= T,
            max.cutoff = 'q99',
            cols = c("lightgrey", "red"),
            split.by='condition') + theme(legend.position = 'right')

vp_N <- VlnPlot(Seurat.object_s,
        feature = 'NOTCH2',
        #slot = "counts", 
        #log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.85, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))

Seurat.object$condition = revalue(Seurat.object$condition, c("iNOA"="iGCA", "CTL"="CTL"))
fp_N <- FeaturePlot(Seurat.object, 
            reduction = "umap", 
            pt.size = 1,
            features = c('NOTCH2'), 
            label = T, 
            order= T,
            max.cutoff = 'q99',
            cols = c("lightgrey", "red"),
            split.by='condition') + theme(legend.position = 'right')


fp_blend <- FeaturePlot(Seurat.object, 
            reduction = "umap", 
            pt.size = 1,
            features = c('DLK1','NOTCH2'),
            label = T, 
            order= T,
             blend = T,       
            max.cutoff = 'q99',
            cols = c("lightgrey", "red","green"),
            split.by='condition') + theme(legend.position = 'right')


layout <- '
AACCCC
BBCCCC
' 


pdf(paste(outdir,'Figure3_DLK1_NOTCH2_blend.pdf',sep=''),width=19, height=7)
wrap_plots(A = vp_D, B = vp_N, C = fp_blend, design = layout)
dev.off()

im = image_read_pdf(paste(outdir,"Figure3_DLK1_NOTCH2_blend.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure3_DLK1_NOTCH2_blend.tiff", sep=''), format = "tiff")
im


pdf(paste(outdir,'NewLabel/','Figure_DLK1.pdf',sep=''),width=6, height=3.5)
plot(vp_D)
dev.off()
pdf(paste(outdir,'NewLabel/','Figure_NOTCH2.pdf',sep=''),width=6, height=3.5)
plot(vp_N)
dev.off()
pdf(paste(outdir,'NewLabel/','Figure_DLK1-NOTCH2_blend.pdf',sep=''),width=12, height=7)
plot(fp_blend)
dev.off()


######## IGF1 IGF2##########
Seurat.object = somatic.integrated.new2
DefaultAssay(Seurat.object) <- "RNA"
Seurat.object_s <- subset(Seurat.object, idents = c("LEY", "MYD","SRT","END","STRO"))

vp_IG1 <- VlnPlot(Seurat.object_s,
        feature = c('IGF1'),
        #slot = "counts", 
        #log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.85, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))

vp_IG2 <- VlnPlot(Seurat.object_s,
        feature = c('IGF2'),
        #slot = "counts", 
        #log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.85, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))

Seurat.object$condition = revalue(Seurat.object$condition, c("iNOA"="iGCA", "CTL"="CTL"))
fp_bl <- FeaturePlot(Seurat.object, 
            reduction = "umap", 
            pt.size = 1,
            features = c('IGF1','IGF2'), 
            label = T, 
            order= T,
             blend = T,       
            max.cutoff = 'q99',
            #cols = c("lightgrey", "red"),
            split.by='condition') + theme(legend.position = 'right')

fp_bl

layout <- '
AACCCC
BBCCCC
' 
pdf(paste(outdir,'Figure3_IGFs_blend.pdf',sep=''),width=19, height=7)
wrap_plots(A = vp_IG1, B = vp_IG2, C = fp_bl, design = layout)
dev.off()

im = image_read_pdf(paste(outdir,"Figure3_IGFs_blend.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure3_IGFs_blend.tiff", sep=''), format = "tiff")
im


pdf(paste(outdir,'NewLabel/','Figure_IGF1.pdf',sep=''),width=6, height=3.5)
plot(vp_IG1)
dev.off()
pdf(paste(outdir,'NewLabel/','Figure_IGF2.pdf',sep=''),width=6, height=3.5)
plot(vp_IG2)
dev.off()
pdf(paste(outdir,'NewLabel/','Figure_IGF1-IGF2_blend.pdf',sep=''),width=12, height=7)
plot(fp_bl)
dev.off()

#### COL1A1 COL1A2 ######
Seurat.object = somatic.integrated.new2
DefaultAssay(Seurat.object) <- "RNA"
Seurat.object_s <- subset(Seurat.object,  idents = c("LEY", "MYD","STRO"))

vp_IG1 <- VlnPlot(Seurat.object_s,
        feature = c('COL1A1'),
        #slot = "counts", 
        #log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.85, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))

vp_IG2 <- VlnPlot(Seurat.object_s,
        feature = c('COL1A2'),
        #slot = "counts", 
        #log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.85, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))

Seurat.object$condition = revalue(Seurat.object$condition, c("iNOA"="iGCA", "CTL"="CTL"))
fp_bl <- FeaturePlot(Seurat.object, 
            reduction = "umap", 
            pt.size = 1,
            features = c('COL1A1','COL1A2'), 
            label = T, 
            order= T,
             blend = T,       
            max.cutoff = 'q99',
            #cols = c("lightgrey", "red"),
            split.by='condition') + theme(legend.position = 'right')


layout <- '
AACCCC
BBCCCC
' 
pdf(paste(outdir,'Figure5_COL1s_blend.pdf',sep=''),width=19, height=7)
wrap_plots(A = vp_IG1, B = vp_IG2, C = fp_bl, design = layout)
dev.off()

im = image_read_pdf(paste(outdir,"Figure5_COL1s_blend.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure5_COL1s_blend.tiff", sep=''), format = "tiff")
im


pdf(paste(outdir,'NewLabel/','Figure_COL1A1.pdf',sep=''),width=6, height=3.5)
plot(vp_IG1)
dev.off()
pdf(paste(outdir,'NewLabel/','Figure_COL1A2.pdf',sep=''),width=6, height=3.5)
plot(vp_IG2)
dev.off()
pdf(paste(outdir,'NewLabel/','Figure_COL1A1-COL1A2_blend.pdf',sep=''),width=12, height=7)
plot(fp_bl)
dev.off()

######## HSD17B3 ########
Seurat.object = somatic.integrated.new2
DefaultAssay(Seurat.object) <- "RNA"
Seurat.object_s <- subset(Seurat.object,  idents = c("LEY"))

vp_HD <- VlnPlot(Seurat.object_s,
        feature = c('HSD17B3'),
        slot = "counts", 
        log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.65, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTL'))

Seurat.object$condition = revalue(Seurat.object$condition, c("iNOA"="iGCA", "CTL"="CTL"))
fp_bl <- FeaturePlot(Seurat.object, 
            reduction = "umap", 
            pt.size = 1,
            features = c('HSD17B3'), 
            label = T, 
            order= T,
            blend = F,       
            max.cutoff = 'q99',
            cols = c("lightgrey", "red"),
            split.by='condition') + theme(legend.position = 'right')



layout <- '
ACCCC
ACCCC
' 

pdf(paste(outdir,'Figure_HSD17B3.pdf',sep=''),width=14, height=4)
wrap_plots(A = vp_HD, C = fp_bl, design = layout)
dev.off()

im = image_read_pdf(paste(outdir,'Figure_HSD17B3.pdf',sep=''),density = 140)
image_write(im, path = paste(outdir,'Figure_HSD17B3.tiff', sep=''), format = "tiff")
im

pdf(paste(outdir,'Figure_HSD17B3.pdf',sep=''),width=9.5, height=3.5)
plot(fp_bl)
dev.off()

im = image_read_pdf(paste(outdir,'Figure_HSD17B3.pdf',sep=''),density = 140)
image_write(im, path = paste(outdir,'Figure_HSD17B3.tiff', sep=''), format = "tiff")
im



#### CCL4 CCL5 ########
Seurat.object = somatic.integrated.new2
DefaultAssay(Seurat.object) <- "RNA"
Seurat.object_s <- subset(Seurat.object,  idents = c("TCL"))

vp_CCL4 <- VlnPlot(Seurat.object_s,
                   feature = c('CCL4'),
                   #slot = "counts", 
                   #log = TRUE,
                   split.plot = F,
                   pt.size = 0) +
  #       split.by = 'condition') +
  ylab("Expression Level") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = c(0.75, 0.9),
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  #geom_boxplot(width = 0.1, outlier.size=1) +
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))

vp_CCL5 <- VlnPlot(Seurat.object_s,
                   feature = c('CCL5'),
                   #slot = "counts", 
                   #log = TRUE,
                   split.plot =F,
                   pt.size = 0) +
  #split.by = 'condition') +
  ylab("Expression Level") +
  xlab("") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = c(0.75, 0.9),
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  #geom_boxplot(width = 0.1, outlier.size=1) +
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))

vp_CCL4 | vp_CCL5

pdf(paste(outdir,'Figure_CCL.pdf',sep=''),width=8, height=4)
plot(vp_CCL4 | vp_CCL5)
dev.off()

im = image_read_pdf(paste(outdir,'Figure_CCL.pdf',sep=''),density = 140)
image_write(im, path = paste(outdir,'Figure_CCL.tiff', sep=''), format = "tiff")


vp_IG1 <- VlnPlot(Seurat.object_s,
        feature = c('IGF1R'),
        slot = "counts", 
        log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.7, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTL'))

vp_IG2 <- VlnPlot(Seurat.object_s,
        feature = c('IGF2R'),
        slot = "counts", 
        log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.7, 0.9),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTL'))

fp_bl <- FeaturePlot(Seurat.object, 
            reduction = "umap", 
            pt.size = 1,
            features = c('IGF1R','IGF2R'), 
            label = T, 
            order= T,
             blend = T,       
            max.cutoff = 'q99',
            #cols = c("lightgrey", "red"),
            split.by='condition') + theme(legend.position = 'right')


pdf(paste(outdir,'Figure3_IGF2R.pdf',sep=''),width=6, height=4)
(vp_IG1) 
dev.off()
pdf(paste(outdir,'Figure3_IGF1R.pdf',sep=''),width=6, height=4)
(vp_IG2) 
dev.off()

im = image_read_pdf(paste(outdir,"Figure3_IGF2R.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure3_IGF2R.tiff", sep=''), format = "tiff")
im

im = image_read_pdf(paste(outdir,"Figure3_IGF1R.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure3_IGF1R.tiff", sep=''), format = "tiff")
im


layout <- '
AACCCC
BBCCCC
' 
pdf(paste(outdir,'Figure3_IGFRs_blend.pdf',sep=''),width=19, height=7)
wrap_plots(A = vp_IG1, B = vp_IG2, C = fp_bl, design = layout)
dev.off()

im = image_read_pdf(paste(outdir,"Figure3_IGFRs_blend.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure3_IGFRs_blend.tiff", sep=''), format = "tiff")


fp_bl

########### N ################
A <- as.data.frame(prop.table(table(Idents(somatic.integrated.new2), 
                                    somatic.integrated.new2$condition),margin = 2))
A2 <- as.data.frame(table(Idents(somatic.integrated.new2), 
                                    somatic.integrated.new2$condition),margin = 2)
A$N <- A2$Freq

N <- ggplot(data=A, aes(x=Var1, y=Freq*100, fill=Var2)) +
geom_bar(stat="identity",  color="black", position=position_dodge()) + 
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTRL')) +
theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
      axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5), 
      axis.title.x = element_text(face = "bold", color = "black", size = 24),
      axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
      axis.title.y = element_text(face = "bold", color = "black", size = 24),
      legend.title = element_text(size = 0),
      legend.text = element_text(face = "bold", color = "black", size = 12),
      legend.position = c(0.9, 0.8),
      panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
labs(x = "CellType", y = "%")

pdf(paste(outdir,'N_iNOAvsCTL.pdf',sep=''),width=10, height=4)
N
dev.off()

im = image_read_pdf(paste(outdir,'N_iNOAvsCTL.pdf',sep=''),density = 200)
image_write(im, path = paste(outdir,'N_iNOAvsCTL.pdf.tiff', sep=''), format = "tiff")
im

#geom_text(size = 7, position = position_stack(vjust = 0.5)) 

###### ACE2 #############
ct='LEY_healthy'
Seurat.object = somatic.integrated.new
Seurat.object <- subset(Seurat.object,  idents = ct)
a = summary(Seurat.object@assays$RNA@counts['ACE2',])
a
mean(Seurat.object@assays$RNA@counts['ACE2',])


###### sunset GEnes ############
CELL_INSL3 <- WhichCells(somatic.integrated.new, slot = 'counts', expression = INSL3 > 0)
somatic.integrated.new.INSL3 <- subset(somatic.integrated.new, cells = CELL_INSL3)
somatic.integrated.new
somatic.integrated.new.INSL3

CELL_ACE2 <- WhichCells(somatic.integrated.new, slot = 'counts', expression = ACE2> 0)
somatic.integrated.new.ACE2 <- subset(somatic.integrated.new, cells = CELL_ACE2)
somatic.integrated.new
somatic.integrated.new.ACE2

CELL_HSD17B3 <- WhichCells(somatic.integrated.new, slot = 'counts', expression = HSD17B3> 0)
somatic.integrated.new.HSD17B3 <- subset(somatic.integrated.new, cells = CELL_HSD17B3)
somatic.integrated.new
somatic.integrated.new.HSD17B3

table(Idents(somatic.integrated.new))


########### Custom VlnPLot ##############
legendiNOA='iGCA'
legendCTRL='CTL'

Seraut.object = somatic.integrated.new 
#.INSL3
#Seraut.object = somatic.integrated.new.ACE2
Seraut.object = somatic.integrated.new.INSL3
#Seraut.object = somatic.integrated.new.HSD17B3
cell_type.2.use = 'LEY'
cell_type = c(paste(cell_type.2.use,'_iNOA',sep=''),
              paste(cell_type.2.use,'_CTL',sep=''))
features = c('INSL3')
#features = c('NOTCH2')
#features = c('ACE2')
#features = c('HSD17B3')
#features = c('VAV3',"CAV1",'CALM1')
#features = c('COL1A1','COL1A2', 'COL4A1')
#features = c('COL1A1','COL1A2')
#features = c('COL4A1','COL4A2')#,'COL4A3','COL4A4')
#features = c('COL1A1','COL1A2','COL4A1')

    
    Cell.subset <- subset(Seraut.object, idents = cell_type , invert=F)
    #subset.counts <- Cell.subset@assays$RNA@counts
    subset.counts <- Cell.subset@assays$RNA@data
    gene.counts.dataframe = data.frame(gene = character(),
                                       cell_barcode = character(),
                                       cell_type = character(),
                                       counts = integer(),
                                       stringsAsFactors=FALSE)
    for (gene in features) {
        gene.counts.tmp = subset.counts[gene,]
        gene.counts.dataframe.entry = data.frame(gene = rep(gene, length(gene.counts.tmp)),
                                                 cell_barcode =  names(gene.counts.tmp),
                                                 condition = Seraut.object[[]][names(gene.counts.tmp),
                                                                               'celltype.cond'],
                                                 counts = (gene.counts.tmp))
        gene.counts.dataframe <- rbind(gene.counts.dataframe, gene.counts.dataframe.entry)
        }
gene.counts.dataframe$condition <- factor(x = gene.counts.dataframe$condition, 
                                          levels = c(paste(cell_type.2.use,"_iNOA",sep=''),
                                                     paste(cell_type.2.use,"_CTL",sep='')))
CVP = ggplot(gene.counts.dataframe, aes(x=gene, y = counts, fill=condition)) +
            geom_split_violin() + 
            #ylim(0,6) +
 #scale_y_log10() +
        scale_x_discrete('LEY') +
            #scale_y_continuous(trans='log') +
            scale_fill_manual(values = c('purple', 'orange')) +
            #geom_jitter(shape=16, size = 0.3, position=position_jitter(0.4)) +  
theme(plot.title = element_text(color="black", size=22, face="bold.italic", hjust = 0.5),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
      axis.text.x = element_blank(),
#          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5),
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.85, 0.87),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#+
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c(legendiNOA,legendCTRL)) +
labs(title = cell_type.2.use , 
     x = "Genes", 
     y = "Expression Level") 
CVP = CVP +labs(title= 'INSL3',
          x = "LEY", 
          y = "Expression Level") 
CVP
assign('INSL3_allcell',CVP)  
assign('INSL3_INSL3cell',CVP)

INSL3_allcell + plot_spacer() + INSL3_INSL3cell
p1 =INSL3_allcell; p2 =INSL3_INSL3cell;
(p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
(p2 + theme(plot.margin = unit(c(0,0,0,30), "pt")))
pdf(paste(outdir,'Figure_INSL3.pdf',sep=''),width=8.5, height=3.5)
plot((p1 + theme(plot.margin = unit(c(0,40,0,15), "pt"))) +
       (p2 + theme(plot.margin = unit(c(0,15,0,40), "pt"))))
dev.off()

im = image_read_pdf(paste(outdir,'Figure_INSL3.pdf',sep=''),density = 140)
image_write(im, path = paste(outdir,'Figure_INSL3.tiff', sep=''), format = "tiff")
im


SO = somatic.integrated.new2
DefaultAssay(SO) <- 'RNA'
c1 = VlnPlot(subset(SO, idents = 'LEY' , invert=F),
        features = c('INSL3'), ncol = 1,
        slot = "counts", 
        log = F,
        split.plot = F,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expr Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.58, 0.82),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTL'))
c1

SO = somatic.integrated.new2
DefaultAssay(SO) <- 'RNA'
VlnPlot(subset(SO, idents = 'TCL' , invert=F),
        features = c('CD3E'), ncol = 1,
        slot = "counts", 
        log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition')

options(repr.plot.width=10, repr.plot.height=3)
SO = somatic.integrated.new2
DefaultAssay(SO) <- 'RNA'
c1 = VlnPlot(subset(SO, idents = 'END' , invert=F),
        features = c('COL4A1'), ncol = 1,
        slot = "counts", 
        log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expr Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.58, 0.82),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTL'))
c2 = VlnPlot(subset(SO, idents = 'END' , invert=F),
        features = c('COL4A2'), ncol = 1,
        slot = "counts", 
        log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expr Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.58, 0.82),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTL'))
c3 = VlnPlot(subset(SO, idents = 'END' , invert=F),
        features = c('COL4A3'), ncol = 1,
        slot = "counts", 
        log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expr Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.58, 0.82),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTL'))
c4 = VlnPlot(subset(SO, idents = 'END' , invert=F),
        features = c('COL4A4'), ncol = 1,
        slot = "counts", 
        log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expr Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.58, 0.82),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTL'))
c1 | c2 | c3 | c4

pdf(paste(outdir,"Figure5_COL4_",cell_type.2.use,"_2.pdf",sep=''),  width=10, height=3)
c1 | c2 | c3 | c4
dev.off()
im = image_read_pdf(paste(outdir,"Figure5_COL4_",cell_type.2.use,"_2.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure5/Figure5_COL4_",cell_type.2.use,"_2.tiff", sep=''), format = "tiff")
im

save.image(paste(outdir,"integration.RData",sep=''))

c1 = VlnPlot(subset(SO, idents = 'SRT' , invert=F),
        features = c('COL1A1'), ncol = 1,
        slot = "counts", 
        log = TRUE,
        split.plot = T,
        pt.size = 0,
       split.by = 'condition') +
ylab("Log Expr Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = c(0.58, 0.82),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iNOA','CTL'))
c1 

pdf(paste(outdir,"Figure5_COL4_",'SRT',".pdf",sep=''),  width=3, height=3)
c1 
dev.off()
im = image_read_pdf(paste(outdir,"Figure5_COL4_",'SRT',".pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure5/Figure5_COL4_",'SRT',".tiff", sep=''), format = "tiff")
im



pdf(paste(outdir,"Figure3_INSL3_",cell_type.2.use,"_allcells.pdf",sep=''),  width=4.4, height=6)
CVP
dev.off()
im = image_read_pdf(paste(outdir,"Figure3_INSL3_",cell_type.2.use,"_allcells.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure3_INSL3_",cell_type.2.use,"_allcells.tiff", sep=''), format = "tiff")
im

mt = FindMarkers(Seraut.object,
            ident.1 = paste(cell_type.2.use,"_iNOA", sep=''),
            ident.2 = paste(cell_type.2.use,"_CTL", sep=''), 
            verbose = FALSE,
            feature=features,
            logfc.threshold=0,
            min.pct = 0)
layout <- '
#AAA#
#AAA#
BBBBB
'
mt = signif(mt, digits = 3)

wrap_plots(A = CVP , B = gridExtra::tableGrob(t(mt[, c('avg_logFC', 'p_val_adj')]), 
                                              theme = ttheme_minimal(padding = unit(c(8, 2), "mm"))), 
           design = layout)

pdf(paste(outdir,"Figure5_COL4_",cell_type.2.use,".pdf",sep=''),  width=12, height=6)
CVP
dev.off()
im = image_read_pdf(paste(outdir,"Figure5_COL4_",cell_type.2.use,".pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure5/Figure5_COL4_",cell_type.2.use,".tiff", sep=''), format = "tiff")
im

?ttheme_minimal()

library(grid)
library(gridExtra)

gene = c('AR','IGF1','IGF1R','IGF2','IGF2R',
          'COL1A1','COL1A2','COL4A1','COL4A2','COL4A3','COL4A4','COL4A5','COL4A6',
          'DLK1','NOTCH2','NOTCH3','HSD17B3')
#gene = c('ACE2')
#gene = c('MAGEB4', 'NANOS1', 'NR0B1','NR5A1','SOHLH1','SYCE1','TAF4B','TEX11','TEX15','WT1','ZMYND15')

cellTypes = c("STRO",
              'MCR',
              'END',
              "LEY",
              "MYD",
              "SRT")

#gene = 'INSL3'
#cellTypes = 'LEY'

Seurat.object = somatic.integrated.new
azoospermia.modulated.gene = data.frame()
for (cell_type in cellTypes) {
    
    azoospermia.response.gene <- FindMarkers(Seurat.object,
                                        ident.1 = paste(cell_type,"_iNOA", sep=''),
                                        ident.2 = paste(cell_type,"_CTL", sep=''), 
                                        verbose = FALSE,
                                        feature=gene,
                                        logfc.threshold=0,
                                        min.pct = 0)
    azoospermia.response.gene <- azoospermia.response.gene[gene,] 
    azoospermia.response.gene$geneID <- rownames(azoospermia.response.gene)
    azoospermia.response.gene$cellType <- cell_type 
    azoospermia.modulated.gene <- rbind(azoospermia.modulated.gene, azoospermia.response.gene)  
}   

myorder=c('cellType','geneID','avg_logFC','p_val_adj','p_val','pct.1','pct.2')
azoospermia.modulated.gene <- azoospermia.modulated.gene[,myorder]
write.xlsx(azoospermia.modulated.gene,
               file= paste(outdir,'azoospermia.modulated.gene.azoospermia.xlsx',sep=''), 
               row.names = F,
               asTable = T)

############ Paternal Genes #############
PD='/Users/tascini.annasofia/Dropbox (HSR Global)/Alfano_904_paperDraft/PaternalGenes/'

Pgene.table = read.xlsx(paste(PD,'Paternal_Genes.xlsx',sep=''),
                       sheet='NOLINK')
row.names(Pgene.table) <- Pgene.table$Gene

feature = row.names(somatic.integrated.new)
NF=setdiff(Pgene.table$Gene, feature)
NewTRY = unlist(strsplit(as.character(Pgene.table[NF,]$Aliases), ','))
a = setdiff(unique(NewTRY[is.na(NewTRY)==FALSE]), feature)
RESCUED = setdiff(unique(NewTRY[is.na(NewTRY)==FALSE]),a)

PG_azoospermia = c(intersect(Pgene.table$Gene, feature), RESCUED)

cellTypes = c("LEY",
              "MYD",
              "SRT",
              'MCR',
              'END',
              "STRO")

gene = PG_azoospermia
Seurat.object = somatic.integrated.new
azoospermia.modulated.gene = data.frame()
for (cell_type in cellTypes) {
    azoospermia.response.gene <- FindMarkers(Seurat.object,
                                        ident.1 = paste(cell_type,"_iNOA", sep=''),
                                        ident.2 = paste(cell_type,"_CTL", sep=''), 
                                        verbose = FALSE,
                                        feature=gene,
                                        logfc.threshold=0,
                                        min.pct = 0)
    azoospermia.response.gene <- azoospermia.response.gene[gene,] 
    azoospermia.response.gene$geneID <- rownames(azoospermia.response.gene)
    azoospermia.response.gene$cellType <- cell_type 
    azoospermia.modulated.gene <- rbind(azoospermia.modulated.gene, azoospermia.response.gene)  
}   

myorder=c('cellType','geneID','avg_logFC','p_val_adj','p_val','pct.1','pct.2')
azoospermia.modulated.gene <- azoospermia.modulated.gene[,myorder]
write.xlsx(azoospermia.modulated.gene,
               file= paste(PD,'azoospermia.modulated.gene.Paternal.xlsx',sep=''), 
               row.names = F,
               asTable = T)


Seurat_object <- somatic.integrated.new2

DefaultAssay(Seurat_object) <- "RNA"
markers.to.plot <- PG_azoospermia
dp = DotPlot(Seurat_object, 
        col.min = -1.5,
        col.max = 1.5,
        features = rev(markers.to.plot), 
        #cols = c("lightgrey", "red"), 
        cols = c('orange', 'dodgerblue2'),    
        dot.scale = 8,
            split.by ='condition') +
RotatedAxis() +
scale_size_continuous(range = c(1,10)) +
theme(axis.text.x = element_text(angle = 90,  color = "black", size=20, vjust =.5, hjust =1, family = 'mono'), 
      axis.title.x = element_text(face = "bold", color = "black", size = 20),
      axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=18),
      axis.title.y = element_text(face = "bold", color = "black", size = 20),
      legend.text = element_text(face = "bold", color = "black", size = 16),
      legend.position="right",
      panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
labs(x = "Genes", y = "Cell type")

dp

jpeg(paste(PD,'Paternal_genes_expression_inClusters_DotPlot_iNOAvsCTL.jpg',sep=''),width=1300, height=600, unit='px')
plot(dp)
dev.off()

pdf(paste(outdir,'Figure3_panel_F.pdf',sep=''),width=18, height=6)
plot(dp)
dev.off()

im = image_read_pdf(paste(outdir,"Figure3_panel_F.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure3_panel_F.tiff", sep=''), format = "tiff")
im

PG_count = as.data.frame(somatic.integrated.new@assays$RNA@counts[PG_azoospermia,])
df = data.frame()
for (gene in PG_azoospermia) {
    df.tmp = data.frame(count = as.numeric(PG_count[gene,]),
                       Gene = gene)
    df = rbind(df,df.tmp)
}

p <- ggplot(df, aes(x=Gene, y=count+1, fill = Gene)) + 
geom_boxplot() +
scale_y_continuous(trans='log10') +
theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
      axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=22, vjust =.5, hjust =1), 
      axis.title.x = element_text(face = "bold", color = "black", size = 24),
      axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
      axis.title.y = element_text(face = "bold", color = "black", size = 24),
      legend.text = element_text(face = "bold", color = "black", size = 12),
      legend.position="top",
      panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
labs(x = "Genes", y = "Log Expression Level") + NoLegend()
options(repr.plot.width=16, repr.plot.height=5)
p

jpeg(paste(PD,'Paternal_genes_expression_SOMATIC.jpg',sep=''),width=1600, height=500, unit='px')
plot(p)
dev.off()



######### Downregulated pathways #######
CND = 'down'
cell_types= c("LEY", "MYD", "SRT", "MCR","END", "STRO")
database = c('Reactome_2016','GO_Biological_Process_2018')
gene.l = list()
path.l = list()
gene.list = character()
indir = paste(outdir,'enrichR/', sep='')

for (cell_type in cell_types) {
  print(cell_type)
  gene.list = character()
  enrichR.file = paste(indir,'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
  
  
  enrichR.table = data.frame()
  for (dat in database) {
    Table <- read.xlsx(xlsxFile = enrichR.file, 
                       sheet = dat, 
                       startRow = 1, 
                       colNames = TRUE,
                       rowNames = TRUE, 
                       detectDates = FALSE, 
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA", 
                       fillMergedCells = FALSE)
    
    enrichR.table = rbind(enrichR.table , Table)
  }
  p = row.names(enrichR.table[enrichR.table$Adjusted.P.value < 0.05,])
  path.l[[cell_type]] = p
  
  for (pathway in p) {
    gene.list.tmp <- unlist(strsplit(enrichR.table[pathway,]$Genes, ';'))
    gene.list = c(gene.list, gene.list.tmp)
    
  }
  gene.l[[cell_type]] = unique(gene.list)
  write.table(gene.list, paste(outdir,'DW_',cell_type,'.txt',sep=''),quote=T, row.names=F, col.names=F)
}


outersect <- function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}

single.l=list()
only2.l = list()
for (i in cell_types){
    only1 = setdiff(gene.l[[i]],unique(as.character(unlist(gene.l[outersect(cell_types[1:5],i)]))))
    single.l[[i]] = only1
    for (j in cell_types) {
        if (j!=i) {
            only2 = setdiff(unique(as.character(unlist(gene.l[c(i,j)]))),
                            unique(as.character(unlist(gene.l[outersect(cell_types[1:5],c(i,j))]))))
            only2.l[[paste(i,j)]] = only2
        }
        }
}
prova = setdiff(unique(as.character(unlist(gene.l))),
                as.character(unlist(single.l)))
prova2 = setdiff(unique(as.character(unlist(gene.l))),
                 unique(as.character(unlist(only2.l))))

length(prova2)


Seurat_object <- somatic.integrated.new
order.2.plot = c('LEY_iNOA', 'LEY_CTL', 'MYD_iNOA','MYD_CTL','SRT_iNOA','SRT_CTL',
              'END_iNOA','END_CTL','MCR_iNOA','MCR_CTL','TCL_iNOA',
             'STRO_iNOA','STRO_CTL','UND_iNOA','UND_CTL')

levels(Seurat_object) <- order.2.plot

DefaultAssay(Seurat_object) <- "RNA"
#markers.to.plot <- AgeGenesClustersTable_dw$gene.name[AgeGenesClustersTable_dw$cluster==3]
dp_down = DotPlot(Seurat_object, 
        col.min = -1,
        col.max = 1,
        features = prova2, 
        cols = c("lightblue", "red"), 
        dot.scale = 8) +
scale_size_continuous(range = c(0.2,5)) +
RotatedAxis() +
coord_flip() +
theme(axis.text.x = element_text(angle = 45, face = "bold", color = "black", size=15,vjust =1.1), 
      axis.title.x = element_text(face = "bold", color = "black", size = 20),
      axis.text.y = element_text(angle = 0, color = "black", size=11, family = 'mono'), 
      axis.title.y = element_text(face = "bold", color = "black", size = 20),
      legend.text = element_text(face = "bold", color = "black", size = 16),
      legend.position="right",
      panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
labs(x = "Genes", y = "Cell type")

dp_down
pdf(paste(outdir,'DW_genes_exp.pdf',sep=''),width=8, height=19)
plot(dp_down)
dev.off()


p = rev(as.character(sort(Reduce(intersect, path.l))))
pathways.dataframe = data.frame(cellType = character(),
                                Pathway = character(),
                                gene.ratio = numeric(),
                                p.value = numeric(),
                                p.value.adj = numeric(),
                                stringsAsFactors=FALSE)
dir = paste(outdir,'enrichR/',sep='')
all.cell.types= c("LEY",
                  "MYD",
                  "SRT",
                  "MCR",
                  "END",
                  "STRO",
                  "UND")
database = c('Reactome_2016','GO_Biological_Process_2018')
N=length(p)
CND = 'DOWN'
fx <- function(x) eval(parse(text=enrichR.table$Overlap[x]))
fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))


for (cell_type in all.cell.types) {
  print(cell_type)
  enrichR.file = enrichR.file = paste(dir,'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
  enrichR.table = data.frame()
  for (dat in database) {
    Table <- read.xlsx(xlsxFile = enrichR.file, 
                       sheet = dat, 
                       startRow = 1, 
                       colNames = TRUE,
                       rowNames = TRUE, 
                       detectDates = FALSE, 
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA", 
                       fillMergedCells = FALSE)
    
    enrichR.table = rbind(enrichR.table , Table)
  }
  
  pathways.dataframe.entry = data.frame(cellType = rep(cell_type, length(p)),
                                        Pathway = p,
                                        gene.ratio = sapply(p, fx),
                                        p.value = enrichR.table[p,]$P.value,
                                        p.value.adj = enrichR.table[p,]$Adjusted.P.value)
  pathways.dataframe <- rbind(pathways.dataframe, pathways.dataframe.entry)       
}

removeHomo <- function(str, stopwords='Homo') {
  x <- unlist(strsplit(as.character(str), " "))
  paste(x[!x %in% stopwords], collapse = " ")
}
removeSapiens <- function(str, stopwords='sapiens') {
  x <- unlist(strsplit(as.character(str), " "))
  paste(x[!x %in% stopwords], collapse = " ")
}
pathways.dataframe$Pathway <- sapply(pathways.dataframe$Pathway, removeHomo)
pathways.dataframe$Pathway <- sapply(pathways.dataframe$Pathway, removeSapiens)


pathways.dataframe$Pathway <- factor(x = pathways.dataframe$Pathway, 
                                     levels = sort(unique(pathways.dataframe$Pathway),decreasing=T))


coloriamo = colors_Alf[all.cell.types,]$cols

patplot = ggplot(data=pathways.dataframe, aes(x=Pathway, y=-log10(p.value), fill=cellType)) +
geom_bar(stat="identity", position=position_dodge()) +
scale_fill_manual(values=as.character(coloriamo)) +
coord_flip() +
scale_x_discrete(position = 'top') + 
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = 'black', size=13, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=9, family = 'mono'),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 14),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))

pdf(paste(outdir,'DW_pathPLOT.pdf',sep=''),width=12, height=18)
plot(patplot)
dev.off()

patplot

layout <- '
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
'

p1 = wrap_plots(B = patplot, A = dp_down, design = layout) + plot_layout(guides = 'collect')

pdf(paste(outdir,'Figure7_main.pdf',sep=''),width=20, height=18)
p1
dev.off()

im = image_read_pdf(paste(outdir,"Figure7_main.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure7_main.pdf.tiff", sep=''), format = "tiff")
im

col = as.character(colors_Alf[names(gene.l),]$cols)

venn(gene.l, simplify=TRUE, opacity = 0.3, box = FALSE, ilab=TRUE, zcolor = col, ilcs = 1.5, sncs = 2)

pdf(paste(outdir,'Figure7_venn.pdf',sep=''),width=7, height=6.5)
venn(path.l, simplify=TRUE, opacity = 0.3, box = FALSE, elipse = T, ilab=TRUE, zcolor = col, ilcs = 1.5, sncs = 2.5)
dev.off()
im = image_read_pdf(paste(outdir,"Figure7_venn.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure7_venn.pdf.tiff", sep=''), format = "tiff")
im

a = venn(path.l, intersections = T)
list_intersection = attr(x = a, "intersections")
names(list_intersection) <- str_replace_all(names(list_intersection), ':', '_')
write.xlsx(list_intersection, paste(outdir,'intersection_Pathway.xlsx'))

df = data.frame()
for (element in names(list_intersection)) {
  entry = data.frame(Intersection = element, Pathways = list_intersection[[element]])
  df = rbind(df, entry)
}
write.xlsx(df, paste(outdir,'intersection_Pathway_2.xlsx'), asTable = T)


cell_type = 'MCR'
enrichR.file = paste(indir,'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
enrichR.table = data.frame()
for (dat in database) {
  Table <- read.xlsx(xlsxFile = enrichR.file, 
                     sheet = dat, 
                     startRow = 1, 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     detectDates = FALSE, 
                     skipEmptyRows = TRUE,
                     skipEmptyCols = TRUE,
                     na.strings = "NA", 
                     fillMergedCells = FALSE)
  
  enrichR.table = rbind(enrichR.table , Table)
}
enrichR.table <- enrichR.table[order(enrichR.table$Adjusted.P.value, decreasing = F),]

enrichR.table[grepl('amyloid', rownames(enrichR.table)),]



#### Pathway UP ####
CND = 'up'
database = c('Reactome_2016','GO_Biological_Process_2018')
#indir = paste(outdir, 'enrichR/',sep='')
indir =  enrichR.path
cell.type_2_plot = c("LEY", "MYD", "MCR")
pathway_to_plot <- list()
pathway_to_plot[["LEY"]] <- c("Cholesterol biosynthesis Homo sapiens R-HSA-191273",
                              "Extracellular matrix organization Homo sapiens R-HSA-1474244",
                              "Collagen formation Homo sapiens R-HSA-1474290",
                              "Collagen biosynthesis and modifying enzymes Homo sapiens R-HSA-1650814",
                              "collagen fibril organization (GO:0030199)",
                              "IRE1alpha activates chaperones Homo sapiens R-HSA-381070",
                              "regulation of cell migration (GO:0030334)",
                              "regulation of keratinocyte apoptotic process (GO:1902172)")

pathway_to_plot[["MYD"]] <- c("Cholesterol biosynthesis Homo sapiens R-HSA-191273",
                              "type I interferon signaling pathway (GO:0060337)",
                              "regulation of interferon-gamma-mediated signaling pathway (GO:0060334)",
                              "Extracellular matrix organization Homo sapiens R-HSA-1474244",
                              "Collagen biosynthesis and modifying enzymes Homo sapiens R-HSA-1650814",
                              "Scavenging by Class A Receptors Homo sapiens R-HSA-3000480",
                              "regulation of cell differentiation (GO:0045595)",
                              "glutathione derivative metabolic process (GO:1901685)",
                              "negative regulation of cellular macromolecule biosynthetic process (GO:2000113)",
                              "regulation of cell migration (GO:0030334)")


pathway_to_plot[["MCR"]] <- c("MHC class II antigen presentation Homo sapiens R-HSA-2132295",
                              "PD-1 signaling Homo sapiens R-HSA-389948",
                              "Interferon gamma signaling Homo sapiens R-HSA-877300",
                              "regulated exocytosis (GO:0045055)",
                              "regulation of amyloid-beta formation (GO:1902003)",
                              "neutrophil mediated immunity (GO:0002446)",
                              "homotypic cell-cell adhesion (GO:0034109)",
                             "cytokine-mediated signaling pathway (GO:0019221",
                             "production of miRNAs involved in gene silencing by miRNA (GO:0035196)",
                             "regulation of transcription from RNA polymerase II promoter (GO:0006357)",
                             "regulation of protein phosphorylation (GO:0001932)")

outdir = '/Users/tascini.annasofia/Dropbox (HSR Global)/Alfano_904_paperDraft/Final_updateiNOAvsCTL/'
for (i in cell.type_2_plot) {
  print(i)
  ds = somatic.integrated.new
    #ds = subset(somatic.integrated.new, downsample = 100)
  gl = 6
  if (i == 'MCR') {gl = 3}
    HP.col <- DGE_Heatmap(Seurat.object = ds, 
                             cell_type = i, 
                          enrichR.path = indir,
                             CND = 'up', 
                             database = c('Reactome_2016','GO_Biological_Process_2018'), 
                             pathway_to_plot[[i]],
                         gene.label = gl)
    
    pdf(paste(outdir,'HM_DS_',i,'.pdf',sep=''),  width=9.5, height=9)
    print(HP.col)
    dev.off()
    }
    
pathways2plot = pathway_to_plot
pathways.dataframe = data.frame(cellType = character(),
                                Pathway = character(),
                                gene.ratio = numeric(),
                                p.value = numeric(),
                                p.value.adj = numeric(),
                                stringsAsFactors=FALSE)
#dir = '/Users/tascini.annasofia/Dropbox (HSR Global)/Alfano_904_paperDraft/xlsx.table/somatics_DGE_enrichR/'
dir = indir
all.cell.types= c("LEY",
                 "MCR",                 "MYD")
database = c('Reactome_2016', 'GO_Biological_Process_2018')

CND = 'UP'
fx <- function(x) eval(parse(text=enrichR.table$Overlap[x]))
fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))

for (cell_type in all.cell.types) {
    print(cell_type)
    enrichR.file = paste(dir,"enrichR_DGE_",cell_type,'_',CND,'.xlsx',sep='')
    enrichR.table = data.frame()
    for (dat in database) {
        Table <- read.xlsx(xlsxFile = enrichR.file, 
                            sheet = dat, 
                            startRow = 1, 
                            colNames = TRUE,
                            rowNames = TRUE, 
                            detectDates = FALSE, 
                            skipEmptyRows = TRUE,
                            skipEmptyCols = TRUE,
                            na.strings = "NA", 
                            fillMergedCells = FALSE)
        
        enrichR.table = rbind(enrichR.table , Table)
    }
    enrichR.table <- enrichR.table[order(enrichR.table$Adjusted.P.value, decreasing = F),]
    N = length(pathways2plot[[cell_type]])
    pathways.dataframe.entry = data.frame(cellType = rep(cell_type, N),
                                          Pathway = pathways2plot[[cell_type]],
                                          gene.ratio = sapply(pathways2plot[[cell_type]], fx),
                                          p.value = enrichR.table[pathways2plot[[cell_type]],]$P.value,
                                          p.value.adj = enrichR.table[pathways2plot[[cell_type]],]$Adjusted.P.value)
    
    #pathways.dataframe <- rbind(pathways.dataframe, pathways.dataframe.entry)
    pathways.dataframe <- pathways.dataframe.entry[order(pathways.dataframe.entry$p.value, decreasing = F),]
    pathways.dataframe$Pathway.num = dim(pathways.dataframe)[1]:1
    pathways.dataframe$Pathway.num = as.factor(pathways.dataframe$Pathway.num )
    removeHomo <- function(str, stopwords='Homo') {
      x <- unlist(strsplit(as.character(str), " "))
      paste(x[!x %in% stopwords], collapse = " ")
    }
    removeSapiens <- function(str, stopwords='sapiens') {
      x <- unlist(strsplit(as.character(str), " "))
      paste(x[!x %in% stopwords], collapse = " ")
    }
    pathways.dataframe$Pathway <- sapply(pathways.dataframe$Pathway, removeHomo)
    pathways.dataframe$Pathway <- sapply(pathways.dataframe$Pathway, removeSapiens)
    pathways.dataframe$Pathway <- str_pad(pathways.dataframe$Pathway, width = 70, side = "right", pad = " ")
    
    changeName = TRUE
    if (changeName) {
        if (cell_type == 'MYD') {
            pathways.dataframe$Pathway[5] = "regulation of interferon-gamma-mediated signaling pathway \n (GO:0060334)"
            pathways.dataframe$Pathway[4] = "negative regulation of cellular macromolecule biosynthetic\n process (GO:2000113)"
        }
        }
    
    PP = ggplot(pathways.dataframe, aes(Pathway.num,-log10(p.value))) + 
      geom_point(aes(size = gene.ratio), color = colors_Alf[cell_type,]$cols) +
      scale_size_continuous(range = c(5,15), name = "Gene ratio") +
      ylim(c(2,10)) +
      coord_flip() +
      scale_x_discrete(breaks=pathways.dataframe$Pathway.num, 
                 labels=pathways.dataframe$Pathway,
                 position = "top") +
      theme(plot.title = element_text(color="black", size=18, face="bold.italic"),
      plot.subtitle = element_text(color="black", size=16, face="italic"),
      axis.text.x = element_text(angle = 90, face = "bold", color = 'black', size=13, hjust =1, family = "mono"), 
      axis.title.x = element_text(face = "bold", color = "black", size = 14),
      axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=13, family = "mono"),
      axis.title.y = element_text(face = "bold", color = "black", size = 14),
      legend.text = element_text(color = "black", size = 12),
      legend.title = element_text(face = "bold", color = "black", size = 14),
      panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
      labs(title = 'Upregulated Pathways in Reactome and GO Biological Process', 
         x = "Pathways", 
         y = "-Log10(pvalue)") 
    
    assign(paste('PP',cell_type,sep=''),PP)
    
    pdf(paste(outdir,'PP_',cell_type,'.pdf',sep=''),  width=11, height=9)
    print(PP)
    dev.off()    
}

for (cell_type in all.cell.types) {
  HM <- image_read_pdf(paste(outdir,"HM_DS_",cell_type,".pdf", sep=''),density = 160)
  PP <- image_read_pdf(paste(outdir,"PP_",cell_type,".pdf", sep=''), density = 160)
  im = image_append(c(HM, PP))
  image_write(im, path = paste(outdir,"HM_PP_",cell_type,".pdf", sep=''), format = "pdf")
  image_write(im, path = paste(outdir,"HM_PP_",cell_type,".tiff", sep=''), format = "tiff")
}


#### Pathway Down ####
CND = 'down'
cell_type = 'END'
enrichR.file = paste(indir,'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
enrichR.table = data.frame()
for (dat in database) {
  Table <- read.xlsx(xlsxFile = enrichR.file, 
                     sheet = dat, 
                     startRow = 1, 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     detectDates = FALSE, 
                     skipEmptyRows = TRUE,
                     skipEmptyCols = TRUE,
                     na.strings = "NA", 
                     fillMergedCells = FALSE)
  
  enrichR.table = rbind(enrichR.table , Table)
}
enrichR.table <- enrichR.table[order(enrichR.table$Adjusted.P.value, decreasing = F),]

enrichR.table[grepl("Transcriptional Regulation by TP53",rownames(enrichR.table)),]




CND = 'down'
database = c('Reactome_2016','GO_Biological_Process_2018')
indir = paste(outdir, 'enrichR/',sep='')

cell.type_2_plot = c("LEY", "MYD", "MCR", "END", "STRO")
pathway_to_plot <- list()
pathway_to_plot[["LEY"]] <- c("APC/C:Cdc20 mediated degradation of mitotic proteins Homo sapiens R-HSA-176409",
                              "mitochondrial respiratory chain complex I assembly (GO:0032981)", "mRNA Splicing Homo sapiens R-HSA-72172")

pathway_to_plot[["MYD"]] <- c("Respiratory electron transport Homo sapiens R-HSA-611105",
                              "mitochondrial respiratory chain complex I biogenesis (GO:0097031)",
                              "Smooth Muscle Contraction Homo sapiens R-HSA-445355",
                              "TP53 Regulates Metabolic Genes Homo sapiens R-HSA-5628897",
                              "Regulation of Hypoxia-inducible Factor (HIF) by oxygen Homo sapiens R-HSA-1234174")


pathway_to_plot[["MCR"]] <- c("Regulation of TP53 Activity through Methylation Homo sapiens R-HSA-6804760",
                              "APC/C:Cdc20 mediated degradation of mitotic proteins Homo sapiens R-HSA-176409",
                              "mitochondrial respiratory chain complex I assembly (GO:0032981)",
                              "Signaling by Type 1 Insulin-like Growth Factor 1 Receptor (IGF1R) Homo sapiens R-HSA-2404192",
                              "p53-Dependent G1 DNA Damage Response Homo sapiens R-HSA-69563")


pathway_to_plot[["STRO"]] <- c("Regulation of Hypoxia-inducible Factor (HIF) by oxygen Homo sapiens R-HSA-1234174")

pathway_to_plot[["END"]] <- c("mitochondrial respiratory chain complex I assembly (GO:0032981)",
                              "Regulation of Hypoxia-inducible Factor (HIF) by oxygen Homo sapiens R-HSA-1234174",
                              "p53-Dependent G1 DNA Damage Response Homo sapiens R-HSA-69563")

for (i in cell.type_2_plot) {
  print(i)
  ds = subset(somatic.integrated.new, downsample = 100)
  gl = 8
  #if (i == 'MCR') {gl = 3}
  HP.col <- DGE_Heatmap(Seurat.object = ds, 
                        cell_type = i, 
                        enrichR.path = indir,
                        CND = CND, 
                        database = c('Reactome_2016','GO_Biological_Process_2018'), 
                        pathway_to_plot[[i]],
                        gene.label = gl)
  
  pdf(paste(outdir,'HM_DW_',i,'.pdf',sep=''),  width=6, height=3.5)
  print(HP.col)
  dev.off()
}

pathways2plot = pathway_to_plot
pathways.dataframe = data.frame(cellType = character(),
                                Pathway = character(),
                                gene.ratio = numeric(),
                                p.value = numeric(),
                                p.value.adj = numeric(),
                                stringsAsFactors=FALSE)

CND = 'DOWN'
fx <- function(x) eval(parse(text=enrichR.table$Overlap[x]))
fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))
all.cell.types = c("LEY", "MYD", "MCR", "END")
for (cell_type in all.cell.types) {
  print(cell_type)
  enrichR.file = paste(dir,"enrichR_DGE_",cell_type,'_',CND,'.xlsx',sep='')
  enrichR.table = data.frame()
  for (dat in database) {
    Table <- read.xlsx(xlsxFile = enrichR.file, 
                       sheet = dat, 
                       startRow = 1, 
                       colNames = TRUE,
                       rowNames = TRUE, 
                       detectDates = FALSE, 
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA", 
                       fillMergedCells = FALSE)
    
    enrichR.table = rbind(enrichR.table , Table)
  }
  enrichR.table <- enrichR.table[order(enrichR.table$Adjusted.P.value, decreasing = F),]
  N = length(pathways2plot[[cell_type]])
  pathways.dataframe.entry = data.frame(cellType = rep(cell_type, N),
                                        Pathway = pathways2plot[[cell_type]],
                                        gene.ratio = sapply(pathways2plot[[cell_type]], fx),
                                        p.value = enrichR.table[pathways2plot[[cell_type]],]$P.value,
                                        p.value.adj = enrichR.table[pathways2plot[[cell_type]],]$Adjusted.P.value)
  
  #pathways.dataframe <- rbind(pathways.dataframe, pathways.dataframe.entry)
  pathways.dataframe <- pathways.dataframe.entry[order(pathways.dataframe.entry$p.value, decreasing = F),]
  pathways.dataframe$Pathway.num = dim(pathways.dataframe)[1]:1
  pathways.dataframe$Pathway.num = as.factor(pathways.dataframe$Pathway.num )
  removeHomo <- function(str, stopwords='Homo') {
    x <- unlist(strsplit(as.character(str), " "))
    paste(x[!x %in% stopwords], collapse = " ")
  }
  removeSapiens <- function(str, stopwords='sapiens') {
    x <- unlist(strsplit(as.character(str), " "))
    paste(x[!x %in% stopwords], collapse = " ")
  }
  pathways.dataframe$Pathway <- sapply(pathways.dataframe$Pathway, removeHomo)
  pathways.dataframe$Pathway <- sapply(pathways.dataframe$Pathway, removeSapiens)
  pathways.dataframe$Pathway <- str_pad(pathways.dataframe$Pathway, width = 70, side = "right", pad = " ")
  
  changeName = T
  if (changeName) {
    if (cell_type == 'MCR') {
      #pathways.dataframe$Pathway[5] = "regulation of interferon-gamma-mediated signaling pathway \n (GO:0060334)"
      pathways.dataframe$Pathway[4] = "Signaling by Type 1 Insulin-like Growth Factor 1 Receptor (IGF1R) \n R-HSA-2404192"
    }
  }
  
  PP = ggplot(pathways.dataframe, aes(Pathway.num,-log10(p.value))) + 
    geom_point(aes(size = gene.ratio), color = colors_Alf[cell_type,]$cols) +
    scale_size_continuous(range = c(5,15), name = "Gene ratio") +
    ylim(c(2,10)) +
    coord_flip() +
    scale_x_discrete(breaks=pathways.dataframe$Pathway.num, 
                     labels=pathways.dataframe$Pathway,
                     position = "top") +
    theme(plot.title = element_text(color="black", size=18, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = 'black', size=13, hjust =1, family = "mono"), 
          axis.title.x = element_text(face = "bold", color = "black", size = 14),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=13, family = "mono"),
          axis.title.y = element_text(face = "bold", color = "black", size = 14),
          legend.text = element_text(color = "black", size = 12),
          legend.title = element_text(face = "bold", color = "black", size = 14),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(title = 'Upregulated Pathways in Reactome and GO Biological Process', 
         x = "Pathways", 
         y = "-Log10(pvalue)") 
  
  assign(paste('PP',cell_type,sep=''),PP)
  
  pdf(paste(outdir,'PP_dw',cell_type,'.pdf',sep=''),  width=11, height=3.5)
  print(PP)
  dev.off()    
}

for (cell_type in all.cell.types) {
  HM <- image_read_pdf(paste(outdir,"HM_DW_",cell_type,".pdf", sep=''),density = 160)
  PP <- image_read_pdf(paste(outdir,"PP_dw",cell_type,".pdf", sep=''), density = 160)
  im = image_append(c(HM, PP))
  image_write(im, path = paste(outdir,"HM_PP_dw_",cell_type,".pdf", sep=''), format = "pdf")
  image_write(im, path = paste(outdir,"HM_PP_dw_",cell_type,".tiff", sep=''), format = "tiff")
}
######### otherplot######à

cell_type = 'SRT'
HM <- image_read_pdf(paste(outdir,"HM_DS_",cell_type,".pdf", sep=''),density = 140)
PP <- image_read_pdf(paste(outdir,'PP_',cell_type,".pdf", sep=''), density = 140)
VP <- image_read_pdf(paste(outdir,'VP_SRT_pathway.pdf', sep=''), density = 140)
im = image_append(c(PP, VP),stack = T)
im2 = image_append(c(HM, im),stack = F)
#im = image_annotate(im, caption, size = 75, location = "+50+650", color = "black")
im2
image_write(im2, path = paste(outdir,"HM_PP_",cell_type,".pdf", sep=''), format = "pdf")
image_write(im2, path = paste(outdir,"HM_PP_",cell_type,".tiff", sep=''), format = "tiff")

cell_type = 'STRO'
HM <- image_read_pdf(paste(outdir,"HM_DS_",cell_type,".pdf", sep=''),density = 140)
PP <- image_read_pdf(paste(outdir,'PP_',cell_type,".pdf", sep=''), density = 140)
VP <- image_read_pdf(paste(outdir,'VP_',cell_type,'_pathway.pdf', sep=''), density = 140)
im = image_append(c(PP, VP),stack = T)
im2 = image_append(c(HM, im),stack = F)
#im = image_annotate(im, caption, size = 75, location = "+50+650", color = "black")
im2
image_write(im, path = paste(outdir,"HM_PP_",cell_type,".pdf", sep=''), format = "pdf")
image_write(im, path = paste(outdir,"HM_PP_",cell_type,".tiff", sep=''), format = "tiff")

Seurat.object = subset(somatic.integrated.new2, idents = 'END')
DefaultAssay(Seurat.object) <- 'RNA'

genes = c('VAV3',"CAV1",'CALM1')
for (gene in genes) {
p = VlnPlot(Seurat.object,
        feature = gene,
        #slot = "counts", 
        log = TRUE,
        pt.size = 0,
       split.by = 'condition',
          split.plot = T ) +
ylab("Log Expression Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5),
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = 'top',
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c(legendiNOA,legendCTRL)) + NoLegend()
    assign(paste('p_',gene,sep=''),p) 
    }
options(repr.plot.width=12, repr.plot.height=5.5)

pdf(paste(outdir,'VP_END_pathway.pdf',sep=''),  width=12, height=5)
    print(p_VAV3 | p_CAV1 | p_CALM1)
dev.off()

?VlnPlot




path = as.character(pathway_to_plot[["END"]])
path <- sapply(path, removeHomo)
path <- sapply(path, removeSapiens)
caption = paste(as.character(path), 'p_adj_value = 1', sep =', ')

# END
cell_type = 'END'
HM <- image_read_pdf(paste(outdir,"HM_DS_",cell_type,".pdf", sep=''),density = 140)
VP <- image_read_pdf(paste(outdir,'VP_END_pathway.pdf', sep=''), density = 140)
im = image_append(c(HM, VP))
im = image_annotate(im, caption, size = 75, location = "+50+650", color = "black")
im
image_write(im, path = paste(outdir,"HM_PP_",cell_type,".pdf", sep=''), format = "pdf")
image_write(im, path = paste(outdir,"HM_PP_",cell_type,".tiff", sep=''), format = "tiff")

?image_append

pathways.dataframe$Pathway.num = dim(pathways.dataframe)[1]:1
pathways.dataframe$Pathway.num = as.factor(pathways.dataframe$Pathway.num )

removeHomo <- function(str, stopwords='Homo') {
  x <- unlist(strsplit(as.character(str), " "))
  paste(x[!x %in% stopwords], collapse = " ")
}
removeSapiens <- function(str, stopwords='sapiens') {
  x <- unlist(strsplit(as.character(str), " "))
  paste(x[!x %in% stopwords], collapse = " ")
}
pathways.dataframe$Pathway <- sapply(pathways.dataframe$Pathway, removeHomo)
pathways.dataframe$Pathway <- sapply(pathways.dataframe$Pathway, removeSapiens)

PP <- ggplot(pathways.dataframe, aes(Pathway.num,-log10(p.value))) + 
geom_point(aes(size = gene.ratio, col=cellType)) +
#scale_colour_manual(values=color) +
#geom_jitter(aes(col=cellType, size=gene.ratio)) +
coord_flip() +
scale_x_discrete(breaks=pathways.dataframe$Pathway.num, 
                 labels=pathways.dataframe$Pathway,
                 position = "top") +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = 'black', size=13, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 14),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=13),
          axis.title.y = element_text(face = "bold", color = "black", size = 13),
          legend.text = element_text(face = "bold", color = "black", size = 14),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) #+ NoLegend()
PP

cell.type_2_plot = c("END")
pathway_to_plot <- list()
pathway_to_plot[["END"]] <- c("Eukaryotic Translation Elongation Homo sapiens R-HSA-156842",
                                      "Selenocysteine synthesis Homo sapiens R-HSA-2408557", 
                                      "rRNA processing Homo sapiens R-HSA-72312")
for (i in cell.type_2_plot) {
    print(i)
    HP.col <- DGE_Heatmap(Seurat.object = somatic.integrated.new, 
                             cell_type = i, 
                             CND = 'down', 
                             database = 'Reactome_2016', 
                             pathway_to_plot[[i]])
    #pdf(paste(outdir,'HM_allcells_DW',i,'_','.pdf',sep=''),  width=10, height=10)
    plot(HP.col)
    #dev.off()
    }

### this could be removed ###
# import excel file
cell_type='LEY'
CND = 'up'
#pathway='Peptide chain elongation Homo sapiens R-HSA-156902'
#pathway="Selenocysteine synthesis Homo sapiens R-HSA-2408557"
pathway='Collagen biosynthesis and modifying enzymes Homo sapiens R-HSA-1650814'
#enrichR.file = paste('/Users/tascini.annasofia/Dropbox (HSR Global)/WORKSPACE/',
#                     'Alfano/',
#                     '904_infertilita_epigenetics/7_bioinfo/Paper_plot/xlsx.table/somatics_DGE_enrichR/',
#                     'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
enrichR.file = paste('/Users/tascini.annasofia/Dropbox (HSR Global)/',
                         'Alfano_904_paperDraft/',
                         'xlsx.table/somatics_DGE_enrichR/',
                         'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
Reactome.table <- read.xlsx(xlsxFile = enrichR.file, 
                            sheet = 'Reactome_2016', 
                            startRow = 1, 
                            colNames = TRUE,
                            rowNames = TRUE, 
                            detectDates = FALSE, 
                            skipEmptyRows = TRUE,
                            skipEmptyCols = TRUE,
                            na.strings = "NA", 
                            fillMergedCells = FALSE)
genes <- unlist(strsplit(Reactome.table[pathway,]$Genes, ';'))

#DoHeatmap
cell2select <- colnames(somatic.integrated.new)[somatic.integrated.new@meta.data$cell_type==cell_type]
cell.type.cells <- subset(somatic.integrated.new, cells = cell2select)
cell.type.cells.scaled <- ScaleData(cell.type.cells)
options(repr.plot.width=15, repr.plot.height=3)
DoHeatmap(cell.type.cells.scaled, 
          features = genes, 
          cells= cell2select,
          size = 3, 
          disp.min = -2, disp.max = 1.5,
          angle = 0,
          draw.lines =T) +
          scale_fill_gradientn(colours = coolwarm(100))
######

cell_type = c("MYD")
Seraut.object <- somatic.integrated.new
DefaultAssay(Seraut.object) <- "RNA"
SelectedCell = character()
    for (i in 1:length(cell_type)) {
        SelectedCell.tmp <- colnames(Seraut.object)[Seraut.object@meta.data$cell_type == cell_type[i]]
        SelectedCell <- c(SelectedCell, SelectedCell.tmp)
    }  
Seraut.object <- subset(Seraut.object, cells = SelectedCell)
Seraut.object <- ScaleData(Seraut.object)
options(repr.plot.width=18, repr.plot.height=7)
DoHeatmap(Seraut.object, 
          features = c("DLK1","NOTCH2"), 
          #cells= cell2select,
          size = 3, 
          disp.min = -2, disp.max = 1.5,
          angle = 0,
          draw.lines = T) +
          scale_fill_gradientn(colours = coolwarm(100))

# create cluster averages
Seurat.object <- somatic.integrated.new
DefaultAssay(Seurat.object) <- 'RNA'
cluster.averages <- AverageExpression(Seurat.object, return.seurat = TRUE)
somatic.cluster.averages <- cluster.averages

cell_type.2.use = 'SRT'
genes.Pregnenolone = c('GAMT', 'GATM','ODC1','CKB','SAT1')
options(repr.plot.width=10, repr.plot.height=6)
cell2select <- colnames(cluster.averages)[Idents(cluster.averages)==paste(cell_type.2.use,'_azoospermia',sep='') |  Idents(cluster.averages)==paste(cell_type.2.use,'_healthy',sep='')]
cluster.subset <- subset(somatic.cluster.averages, 
                         idents = c(paste(cell_type.2.use,'_iNOA',sep=''),
                                    paste(cell_type.2.use,'_CTL',sep='')))
HM <- DoHeatmap(cluster.subset, 
                features = genes.Pregnenolone, 
                    #cells = cell2select,
                    size = 2,
                    disp.min = -3,
                    disp.max = 3,
                    draw.lines = F) + 
scale_fill_gradientn(colours = coolwarm(100))
HM
#pdf(paste(outdir,'HM_ave_',cell_type.2.use,'2.pdf',sep=''),  width=5, height=5)
       # plot(HM)
# dev.off()

cell_type.2.use = c('LEY')
options(repr.plot.width=6, repr.plot.height=5)
cell2select <- colnames(cluster.averages)[Idents(cluster.averages)==paste(cell_type.2.use,'_iNOA',sep='') |  Idents(cluster.averages)==paste(cell_type.2.use,'_healthy',sep='')]
HM <- DoHeatmap(cluster.averages, 
                    features = c('DLK1', 'NOTCH2'), #genes.Pregnenolone,
                    cells = cell2select,
                    size = 2,
                    disp.min = -3,
                    disp.max = 3,
                    draw.lines = F) +
    scale_fill_gradientn(colours = coolwarm(100))

 #pdf(paste(outdir,'HM_ave_',cell_type.2.use,'2.pdf',sep=''),  width=5, height=5)
        plot(HM)
       # dev.off()

cell.type_2_plot = c("LEY", "MYD", "MCR","TCL","END")
#"Sertoli"
pathway_to_plot <- list()
pathway_to_plot[["LEY"]] <- c("Collagen biosynthesis and modifying enzymes Homo sapiens R-HSA-1650814",
                                 "Collagen formation Homo sapiens R-HSA-1474290",
                                 "IRE1alpha activates chaperones Homo sapiens R-HSA-381070")
pathway_to_plot[["MYD"]] <- c("Scavenging by Class A Receptors Homo sapiens R-HSA-3000480",
                                "Collagen biosynthesis and modifying enzymes Homo sapiens R-HSA-1650814",
                                "Cholesterol biosynthesis Homo sapiens R-HSA-191273")
pathway_to_plot[["SRT"]] <- c("Pregnenolone biosynthesis Homo sapiens R-HSA-196108",
                                  "Metabolism of steroid hormones Homo sapiens R-HSA-196071",
                                  "Steroid hormones Homo sapiens R-HSA-209943")
pathway_to_plot[["MCR"]] <- c("MHC class II antigen presentation Homo sapiens R-HSA-2132295",
                                     "Interferon gamma signaling Homo sapiens R-HSA-877300",
                                     "PD-1 signaling Homo sapiens R-HSA-389948")
pathway_to_plot[["TCL"]] <- c("Antigen Presentation: Folding, assembly and peptide loading of class I MHC Homo sapiens R-HSA-983170",
                                 "Interferon alpha/beta signaling Homo sapiens R-HSA-909733",
                                 "PD-1 signaling Homo sapiens R-HSA-389948")
pathway_to_plot[["END"]] <- c("VEGFR2 mediated vascular permeability Homo sapiens R-HSA-5218920")
pathway_to_plot2 <- list()


for (i in cell.type_2_plot) {
    print(i)
    HP.col <- DGE_HeatmapAVE(somatic.cluster.averages, 
                             cell_type = i, 
                             CND = 'up', 
                             database = 'Reactome_2016', 
                             pathway_to_plot[[i]])
    #pdf(paste(outdir,'HM_ave_',i,'_','.pdf',sep=''),  width=5, height=5)
    plot(HP.col)
    #dev.off()
    }
    

cell.type_2_plot = c("Endothelial")
pathway_to_plot <- list()
pathway_to_plot[["Endothelial"]] <- c("Eukaryotic Translation Elongation Homo sapiens R-HSA-156842",
                                      "Selenocysteine synthesis Homo sapiens R-HSA-2408557", 
                                      "rRNA processing Homo sapiens R-HSA-72312")
for (i in cell.type_2_plot) {
    print(i)
    for (j in pathway_to_plot[[i]]) {
        HP.col <- DGE_HeatmapAVE(cluster.averages, 
                           cell_type = i, 
                           CND = 'down', 
                           database = 'Reactome_2016', 
                           pathway = j)
        pdf(paste(outdir,'HM_ave_',i,'_',gsub("[[:punct:]]", "-", j),'.pdf',sep=''),  width=4.5, height=9)
        plot(HP.col)
        dev.off()
        }
    }

levels(somatic.integrated.new[[]]$condition)

DefaultAssay(somatic.integrated.new) <- 'RNA'
cell_type='SRT'
#signature= c('STAR','FDXR','CYP11A1')
signature = c('GAMT', 'GATM','ODC1','CKB','SAT1')
somatic.integrated.new[["Sign_exp"]] <- apply(FetchData(object = somatic.integrated.new, 
                                       vars = signature),
                                         1,
                                         mean)
cell2select <- colnames(somatic.integrated.new)[somatic.integrated.new@meta.data$cell_type==cell_type]
cell.type.cells <- subset(somatic.integrated.new, cells = cell2select)
cell.type.cells[["Sign_exp"]] <- apply(FetchData(object = cell.type.cells, 
                                       vars = signature),
                                         1,
                                         mean)
#VlnPlot(cell.type.cells, features = "Sign_exp", 
#                 split.by = "condition", 
#        group.by = "cell_type", 
#    pt.size = 0, combine = FALSE)

options(repr.plot.width=18, repr.plot.height=5)
plots <- VlnPlot(cell.type.cells, 
                 features = signature, 
                 split.by = "condition", 
                 group.by = "cell_type", 
    pt.size = 0, combine = FALSE) 
#+ scale_fill_manual(values = c('orange', 'dodgerblue2'))



#pdf(paste(outdir,'Sertoli_signature2.pdf',sep=''),width=10, height=6)
plot(CombinePlots(plots = plots, ncol = 5))
#dev.off()


Seurat.object = subset(somatic.integrated.new2, idents = 'SRT')
DefaultAssay(Seurat.object) <- 'RNA'


genes = c('STAR','FDXR','CYP11A1','AKR1B1','STARD3','GAMT', 'GATM','ODC1','CKB','SAT1')
for (gene in genes) {
p = VlnPlot(Seurat.object,
        feature = gene,
        #slot = "counts", 
        log = TRUE,
        pt.size = 0,
         split.plot = TRUE,
       split.by = 'condition') +
ylab("Log Expr. Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 14, hjust =.5),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=14),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = 'top',
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c(legendiNOA,legendCTRL)) + NoLegend()
    assign(paste('p_',gene,sep=''),p) 
    }


options(repr.plot.width=12, repr.plot.height=5.5)
print((p_STAR | p_CYP11A1 | p_AKR1B1 | p_FDXR) / (p_GAMT | p_GATM | p_ODC1 | p_CKB ))

pdf(paste(outdir,'VP_SRT_pathway.pdf',sep=''),  width=12, height=5)
    print((p_STAR | p_CYP11A1 | p_AKR1B1 | p_FDXR ) / (p_GAMT | p_GATM | p_ODC1 | p_CKB ))
dev.off()
p = (p_STAR | p_CYP11A1 | p_AKR1B1 | p_FDXR ) / (p_GAMT | p_GATM | p_ODC1 | p_CKB )

layout <- '
AACCCC
AADDDD
AADDDD
'
options(repr.plot.width=14, repr.plot.height=6) 
wrap_plots(A = HM, C = PPSTRO, D = p, 
           design = layout)



cell_type = 'SRT'
ds = subset(somatic.integrated.new, downsample = 100)
Seurat.object = ds
cell2select <- colnames(Seurat.object)[Seurat.object@meta.data$cell_type==cell_type]
cell.type.cells <- subset(Seurat.object, cells = cell2select)
cell.type.cells.scaled <- ScaleData(cell.type.cells)
gene.list = c('STAR','FDXR','CYP11A1','AKR1B1','FDX1','GAMT', 'GATM','ODC1','CKB','SAT1')
HM <- DoHeatmap(cell.type.cells.scaled,
                    group.colors = c('orange', 'dodgerblue2'),
                    features = gene.list, 
                    cells= cell2select,
                    size = 3, 
                    disp.min = -2, disp.max = 1.5,
                    angle = 45,
                    draw.lines =T,
                    label = F,
                    raster = F) +
scale_fill_gradientn(colours = coolwarm(200))
HM
pdf(paste(outdir,'HM_DS_',cell_type,'.pdf',sep=''),  width=7, height=5.5)
print(HM)
dev.off()

cell_type = 'SRT'
HM <- image_read_pdf(paste(outdir,"HM_DS_",cell_type,".pdf", sep=''),density = 140)
VP <- image_read_pdf(paste(outdir,'VP_SRT_pathway.pdf', sep=''), density = 140)
im = image_append(c(HM, VP))
#im = image_annotate(im, caption, size = 75, location = "+50+650", color = "black")
im
image_write(im, path = paste(outdir,"HM_PP_",cell_type,".pdf", sep=''), format = "pdf")
image_write(im, path = paste(outdir,"HM_PP_",cell_type,".tiff", sep=''), format = "tiff")

Seurat.object = subset(somatic.integrated.new2, idents = 'STRO')
DefaultAssay(Seurat.object) <- 'RNA'


genes = c('VIM','RACK1','RBM3','CLU','HSPA1A')
for (gene in genes) {
p = VlnPlot(Seurat.object,
        feature = gene,
        #slot = "counts", 
        log = TRUE,
        pt.size = 0,
         split.plot = T,
       split.by = 'condition') +
ylab("Log Expr Level") +
xlab("") +
theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 14, hjust =.5),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=14),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = 'top',
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
#geom_boxplot(width = 0.1, outlier.size=1) +
scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c(legendiNOA,legendCTRL)) + NoLegend()
    assign(paste('p_',gene,sep=''),p) 
    }

print((p_RACK1 | p_CLU | p_RBM3) / (p_VIM | p_HSPA1A))

pdf(paste(outdir,'STRO_genes.pdf',sep=''),  width=12, height=3)
    print(p_RACK1 | p_CLU | p_RBM3 | p_VIM | p_HSPA1A)
dev.off()

cell_type = 'STRO'
HM <- image_read_pdf(paste(outdir,"HM_DS_",cell_type,".pdf", sep=''),density = 140)
PP <- image_read_pdf(paste(outdir,"PP_",cell_type,".pdf", sep=''),density = 140)
VP <- image_read_pdf(paste(outdir,'STRO_genes.pdf', sep=''), density = 140)
im1 = image_append(c(PP, VP), stack=T)
im = image_append(c(HM, im1))
#im = image_annotate(im, caption, size = 75, location = "+50+650", color = "black")
im
image_write(im, path = paste(outdir,"HM_PP_",cell_type,".pdf", sep=''), format = "pdf")
image_write(im, path = paste(outdir,"HM_PP_",cell_type,".tiff", sep=''), format = "tiff")

options(repr.plot.width=8, repr.plot.height=3)
DefaultAssay(somatic.integrated.new2) <- "RNA"
genes.Pregnenolone = c("STAR", "CYP11A1","STARD4", "AKR1B1", "TSPO","FDX1","STARD3NL","STARD3","FDXR")#,"TSPOAP1")
markers.to.plot <- genes.Pregnenolone
#subset(pbmc, idents = c("NK", "B"))
#subset(somatic.integrated.3, idents = c("Leydig", "Myoid"))
#DotPlot(somatic.integrated.3, features = rev(markers.to.plot), 
#        cols = c("blue", "red"), dot.scale = 8, 
#        split.by = "condition") 
DotPlot(subset(somatic.integrated.new2, idents = c("SRT")), 
        features = rev(markers.to.plot), 
        cols = c("red", "red"), dot.scale = 8, 
        split.by = "condition") + RotatedAxis()

cell.type_2_plot = 'LEY'
pathway_to_plot = "Collagen biosynthesis and modifying enzymes Homo sapiens R-HSA-1650814"
pathway_to_plot2 = "Collagen biosynthesis and modifying enzymes Homo sapiens R-HSA-1650814"
HP.col <- DGE_HeatmapAVE(somatic.cluster.averages, 
            cell_type = cell.type_2_plot,
            CND = 'up',
            database = 'Reactome_2016',
            pathway = pathway_to_plot)
#pdf(paste(outdir,'HM_AVE_',cell.type_2_plot,'_',pathway_to_plot2,'.pdf',sep=''),  width=5, height=7)
HP.col
#dev.off()

genes <- c('VIM','DLK1','NOTCH2')
DotPlot(somatic.integrated.new2, 
        split.by = 'condition',
        assay = NULL, 
        features = genes, cols = coolwarm(13),
        col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
        group.by = NULL, 
        #split.by = 'cell_type', 
        scale.by = "radius",
        scale.min = NA, scale.max = NA)

options(repr.plot.width=6, repr.plot.height=6)
CellScatter(cluster.averages, 
            cell1 = "LEY_iNOA", 
            cell2 = "LEY_CTL", 
            highlight = genes,
            smooth = T, pt.size = 2) 
# + geom_pointdensity() + scale_color_viridis()

#Plot_sign(somatic.integrated.3, signature = genes)
options(repr.plot.width=12, repr.plot.height=6)
signature = c('STAR','CYP11A1','FDXR')
Seraut.object = somatic.integrated.new2
operator = sum
 x <- Seraut.object
    DefaultAssay(x) <- "RNA"
    x[["Sign_exp"]] <- apply(FetchData(object = x, 
                                       vars = signature),
                             1,
                             operator)
    FP <- FeaturePlot(x,
                      reduction = "umap", 
                      features = 'Sign_exp', 
                      pt.size = 1,
                      label = T, 
                      order = T,
                      cols = c("lightgrey", "red"),#)#,
                      split.by = 'condition')

FP

library(circlize)
library(ComplexHeatmap)
library(GetoptLong)

head(tcell.So@assays$RNA)

tcell.So <- subset(somatic.integrated.new2, idents = 'TCL')
tcell.matrix <- tcell.So@assays$RNA
expr <- tcell.matrix

Endothelial.So <- subset(somatic.integrated.3, idents = 'Endothelial')
Endothelial.So <- ScaleData(Endothelial.So)
Endothelial.matrix <- Endothelial.So@assays$RNA@data
expr <- Endothelial.matrix
expr <- as.matrix(expr)

#expr = expr[apply(expr, 1, function(x) sum(x > 0)/length(x) > 0.5), , drop = FALSE]
base_mean = rowMeans(expr)

cell_type = 'T-cell'
DGE.file = paste('/Users/tascini.annasofia/Dropbox (HSR Global)/WORKSPACE/',
                     'Alfano/',
                     '904_infertilita_epigenetics/7_bioinfo/Paper_plot/xlsx.table/somatics_DGE/',
                     'DGE_',cell_type,'.iNOA.vs.CTRL.xlsx', sep='')
DGE.table <- read.xlsx(xlsxFile = DGE.file, 
                            sheet = 1, 
                            startRow = 1, 
                            colNames = TRUE,
                            rowNames = TRUE, 
                            detectDates = FALSE, 
                            skipEmptyRows = TRUE,
                            skipEmptyCols = TRUE,
                            na.strings = "NA", 
                            fillMergedCells = FALSE)
DGE.gene = row.names(DGE.table[order(DGE.table$avg_logFC, decreasing = T),])

So <- subset(somatic.integrated.new, idents=c('TCL_iNOA','TCL_CTL')) 
So <- ScaleData(So)

Hp <- DoHeatmap(So, 
          features = DGE.gene, 
          size = 3, 
          disp.min = -2.5, disp.max = 2.5,
          angle = 0,
          draw.lines =T) 
Hp + scale_fill_gradientn(colours = coolwarm(100))

expr.DGE <- expr[rownames(expr) %in% DGE.gene,]
base_mean.DGE = rowMeans(expr.DGE)

library(GetoptLong)
options(repr.plot.width=16, repr.plot.height=16)
ht_list = Heatmap(expr.DGE, col = colorRamp2(c(-3.5, 0, 3), c("blue", "white", "red")), 
    name = "scaled_expr", column_title = qq("relative expression for @{nrow(expr.DGE)} genes"),
    show_column_names = T, width = unit(8, "cm"),
    heatmap_legend_param = list(title = "Scaled expr")) +
Heatmap(base_mean.DGE, name = "base_expr", width = unit(5, "mm"),
        heatmap_legend_param = list(title = "Base expr"))

ht_list

get_correlated_variable_genes = function(mat, 
                                         n = nrow(mat), 
                                         cor_cutoff = 0, 
                                         n_cutoff = 0) {
    ind = order(apply(mat, 1, function(x) {
            q = quantile(x, c(0.1, 0.9))
            x = x[x < q[1] & x > q[2]]
            var(x)/mean(x)
        }), decreasing = TRUE)[1:n]
    mat2 = mat[ind, , drop = FALSE]
    dt = cor(t(mat2), method = "spearman")
    diag(dt) = 0
    dt[abs(dt) < cor_cutoff] = 0
    dt[dt < 0] = -1
    dt[dt > 0] = 1

    i = colSums(abs(dt)) > n_cutoff

    mat3 = mat2[i, ,drop = FALSE]
    return(mat3)
}

mat = get_correlated_variable_genes(expr, cor_cutoff = 0.5, n_cutoff = 20)
mat2 = t(apply(mat, 1, function(x) {
    q10 = quantile(x, 0.1)
    q90 = quantile(x, 0.9)
    x[x < q10] = q10
    x[x > q90] = q90
    scale(x)
}))
colnames(mat2) = colnames(mat)

ribosomial.file <- read.table(file = 'ribosomial.genes.txt', sep = '\t')
ribosomial.genes <- ribosomial.file$V2

rpl = rownames(mat) %in% ribosomial.genes 

base_mean = rowMeans(mat)

base_mean_iNOA = rowMeans(mat[,startsWith(colnames(mat), "pz")])
base_mean_CTL = rowMeans(mat[,startsWith(colnames(mat), "Lit")])

library(GetoptLong)
options(repr.plot.width=16, repr.plot.height=16)
ht_list = Heatmap(mat2, col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")), 
    name = "scaled_expr", column_title = qq("relative expression for @{nrow(mat)} genes"),
    show_column_names = FALSE, width = unit(8, "cm"),
    heatmap_legend_param = list(title = "Scaled expr")) +
    Heatmap(base_mean, name = "base_expr", width = unit(5, "mm"),
        heatmap_legend_param = list(title = "Base expr")) +
    Heatmap(base_mean_iNOA, name = "base_expr_INOA", width = unit(5, "mm"),
        heatmap_legend_param = list(title = "Base expr")) +
    Heatmap(base_mean_CTL, name = "base_expr_CTL", width = unit(5, "mm"),
        heatmap_legend_param = list(title = "Base expr")) +
    Heatmap(rpl + 0, name = "ribonucleoprotein", col = c("0" = "white", "1" = "purple"), 
        show_heatmap_legend = FALSE, width = unit(5, "mm")) +
    Heatmap(cor(t(mat2)), name = "cor", 
        col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")), 
        show_row_names = FALSE, show_column_names = FALSE, row_dend_side = "right", 
        show_column_dend = FALSE, column_title = "pairwise correlation between genes",
        heatmap_legend_param = list(title = "Correlation"))
ht_list = draw(ht_list, main_heatmap = "cor")
decorate_column_dend("scaled_expr", {
    tree = column_dend(ht_list)$scaled_expr
    ind = cutree(as.hclust(tree), k = 2)[order.dendrogram(tree)]

    first_index = function(l) which(l)[1]
    last_index = function(l) { x = which(l); x[length(x)] }
    x1 = c(first_index(ind == 1), first_index(ind == 2)) - 1
    x2 = c(last_index(ind == 1), last_index(ind == 2))
    grid.rect(x = x1/length(ind), width = (x2 - x1)/length(ind), just = "left",
        default.units = "npc", gp = gpar(fill = c("#FF000040", "#00FF0040"), col = NA))
})

save(somatic.integrated.new, file = 'somatic.integrated')
save(somatic.integrated.new2, file = 'somatic.integrated.2')

save(somatic.integrated.new, file = paste(outdir, 'somatic.integrated.new', sep='') )
save(somatic.integrated.new2, file = paste(outdir, 'somatic.integrated.new2', sep='') )



