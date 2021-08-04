######## INTEGRATION ##########
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

########## load data ###########
OBJ_dir = '/Users/tascini.annasofia/Dropbox (HSR Global)/WORKSPACE/Alfano/AlfanoM_904_infertilita_epigenetics/7_bioinfo/SeuratObjects/'
SO_dir = "/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/Seurat_objects/"
dir.create(SO_dir)

start_time <- Sys.time()

load(paste(OBJ_dir,"Somatic_integration_beforeplots.RData", sep=''))

end_time <- Sys.time()
end_time - start_time

iNOA_donor1 <- RenameCells(iNOA_donor1, add.cell.id = 'iNOA1_')
iNOA_donor2 <- RenameCells(iNOA_donor2, add.cell.id = 'iNOA2_')
iNOA_donor3 <- RenameCells(iNOA_donor3, add.cell.id = 'iNOA3_')

#### functions ######
```{r functions, include = FALSE}
Gene_conversion = function(x, gtf_dictionary, gtf_OSR_dictionary){
  # convert gene name from two different versions
  gene.ensembl =  gtf_dictionary[which(gtf_dictionary$gene_name == x),]$gene_id
  genec = gtf_OSR_dictionary[gtf_OSR_dictionary$gene_id %in% gene.ensembl,]$gene_name
  if (x %in% genec | length(genec) == 0) {genec = x}
  newgene  = paste(genec, collapse = ';')
  return(newgene)
}


myDGE = function(Seurat.object, cell_type, outdir = './') {
  azoospermia.response <- FindMarkers(Seurat.object,
                                      ident.1 = paste(cell_type,"_iNOA", sep=''),
                                      ident.2 = paste(cell_type,"_CTL", sep=''), 
                                      verbose = FALSE)
  write.xlsx(azoospermia.response,
             file= paste(outdir,'DGE_', cell_type,'.iNOA.vs.CTRL.xlsx',sep=''), 
             row.names = T,
             asTable = T)
  return(azoospermia.response)
}

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


######### color #########
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
                                      "CTLonly",
                                      "LEY2",
                                      "MorrisOnly"), 
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
                                  'grey',
                                  "#FF0066",
                                  "#FF6600"))
rownames(colors_Alf) <- colors_Alf$celltype

############### prepare for integration ###############
testis_Morris$condition <- "CTL"
testis_Morris$source <- "Morris donor"
head(testis_Morris[[]])
testis_Morris <- RenameCells(testis_Morris, add.cell.id = 'Morris_')

dat = 'Morris'
dir = paste('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/',dat,'_iNOA/', sep ='')

somatic.list = list()
somatic.list = list(iNOA_donor1, iNOA_donor2, iNOA_donor3,
                    testis_Morris)
names(somatic.list) <- c('iNOA1', 'iNOA2','iNOA3', "Morris")

###### integration workflow #######
source('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/integration_function.R')
options(future.globals.maxSize = 3145728000) 

print(somatic.list)
int_obj = integration_workflow(somatic.list = somatic.list, 
                               dataset_iNOA_dir = dir, 
                               dataset = dat)



DimPlot(int_obj, label = T)
levels(int_obj)


############# Integration DGE enrichR and plots ########
cell_type <- c( "MYD",
                "LEY", 
                "LEY2",
                "UND",
                "MCR",
                "MorrisOnly",
                "SRT",
                "STRO",
                "END",
                "TCL")

dat = 'Morris2'
dir = paste('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/prj_COSR/AlfanoM_904_infertilita_epigenetics/',dat,'_iNOA/', sep ='')
dir.create(dir)
cell_type2 <- c( "MYD",
                "LEY", 
                "LEY",
                "UND",
                "MCR",
                "MorrisOnly",
                "SRT",
                "STRO",
                "END",
                "TCL")

library(enrichR)
DefaultAssay(int_obj) <- "integrated" 
final_int_obj = integration_postprocessing(SO = int_obj,
                                           new_cluster_id = cell_type2, 
                                           dataset = dat, 
                                           dataset_iNOA_dir = dir,
                                           res = 0.2, 
                                           TCL_DGE = TRUE,
                                           colors_Alf = colors_Alf)
save(final_int_obj, file = paste(SO_dir,'iNOA_',dat,'.integrated', sep =''))

####### N cells ##########

ND = 4

DimPlot(final_int_obj, group.by = 'condition', pt.size = 0.2) + scale_color_manual(values = c('orange', 'dodgerblue2'))
ggsave(filename = paste(dir, 'UMAP_bycondition.png', sep=''),
       width = 10, height = 9,
       units = 'cm')
DimPlot(final_int_obj, group.by = 'source', order = TRUE) + scale_color_manual(values = rainbow(ND))
ggsave(filename = paste(dir, 'UMAP_byDONOR.png', sep=''),
       width = 15, height = 9,
       units = 'cm')
DimPlot(final_int_obj, group.by = "celltype", split.by = 'condition')
table(Idents(final_int_obj))


A <- as.data.frame(prop.table(table(Idents(final_int_obj), 
                                    final_int_obj$condition),margin = 2))
A2 <- as.data.frame(table(Idents(final_int_obj), 
                          final_int_obj$condition),margin = 2)
A$N <- A2$Freq

col = as.character(colors_Alf[levels(final_int_obj),]$cols)
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

A <- as.data.frame(prop.table(table(Idents(final_int_obj), 
                                    final_int_obj$source),margin = 2))
A2 <- as.data.frame(table(Idents(final_int_obj), 
                          final_int_obj$source),margin = 2)
A$N <- A2$Freq

col = as.character(colors_Alf[levels(final_int_obj),]$cols)
N2 <- ggplot(data=A, aes(x=Var2, y=N, fill=Var1)) +
  geom_bar(stat="identity",  color="black", position="fill") + 
  scale_fill_manual(values = col) +
  #scale_x_discrete(labels=c('GSE124263_adult_Donor1'='GSE124263_#1', 
  #                         'GSE124263_adult_Donor2' = 'GSE124263_#2', 
  #                          'iNOA donor 1' = 'iNOA#1', 
  #                         'iNOA donor 2' ='iNOA#2', 
  #                        'iNOA donor 3' ='iNOA#3')) +
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

pdf(paste(dir,'N_iNOAvsCTL.pdf',sep=''),width=15, height=6)
print(N + N2 + plot_layout(guides = "collect"))
dev.off()



#COL = row.names(Seurat.object)[grep("^COL",row.names(Seurat.object))]
gene = c("DLK1", "NOTCH2", "IGF1", "IGF2", "INSL3", "HSD17B3",
         "CDKN2A",
         "COL1A1", "COL1A2", "COL4A1", "COL4A2", 
         'H19', 'MEG3', 'MEG8', 'HYMAI', 'PEG3', 'ERAP2', 'RTL1', "KCNQ1OT1")

Seurat.object = final_int_obj
DefaultAssay(Seurat.object) <-"RNA"
Idents(Seurat.object) <- "celltype.cond"

azoospermia.modulated.gene = data.frame()
for (ct in cell_type) {
  #print(ct)
  azoospermia.response.gene <- data.frame()
  azoospermia.response.gene <- FindMarkers(Seurat.object,
                                           ident.1 = paste(ct,"_iNOA", sep=''),
                                           ident.2 = paste(ct,"_CTL", sep=''), 
                                           verbose = FALSE,
                                           feature=gene,
                                           logfc.threshold=0,
                                           min.pct = 0)
  
  azoospermia.response.gene          <- azoospermia.response.gene[gene,] 
  azoospermia.response.gene$geneID   <- rownames(azoospermia.response.gene)
  azoospermia.response.gene$cellType <- ct 
  azoospermia.modulated.gene         <- rbind(azoospermia.modulated.gene, azoospermia.response.gene)  
}

#azoospermia.modulated.gene 
myorder=c('cellType','geneID','avg_logFC','p_val_adj','p_val','pct.1','pct.2')
azoospermia.modulated.gene <- azoospermia.modulated.gene[,myorder]
write.xlsx(azoospermia.modulated.gene,
           file= paste(dir,'azoospermia.modulated.gene.azoospermia.xlsx',sep=''), 
           row.names = FALSE,
           asTable = TRUE)

####### CDKN2A ########
Seurat.object = final_int_obj
DefaultAssay(Seurat.object) <- "RNA"
Seurat.object_s <- Seurat.object # subset(Seurat.object,  idents = c("LEY", "MYD","SRT"))
#Seurat.object_s = Seurat.object
vp_D <- VlnPlot(Seurat.object_s,
                feature = "CDKN2A",
                #slot = "counts", 
                #log = TRUE,
                split.plot = T,
                pt.size = 2,
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
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL')) +
  scale_color_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))

fp_D <- FeaturePlot(Seurat.object, 
                    reduction = "umap", 
                    features = c("CDKN2A"), 
                    label = T, 
                    order= T,
                    pt.size = 2,
                    max.cutoff = 'q99',
                    cols = c("lightgrey", "red"),
                    split.by='condition') + theme(legend.position = 'right')

vp_D
fp_D
ggsave(filename = paste(dir,"Featureplot_CDKN2A.png",sep=''),
       plot = fp_D,
       width = 20, height = 10,
       units = 'cm')
ggsave(filename = paste(dir,"Violinplot_CDKN2A.png",sep=''),
       plot = vp_D,
       width = 32, height = 16,
       units = 'cm')

##### DLK1 NOTCH2 #####

#DLK1
Seurat.object = final_int_obj
DefaultAssay(Seurat.object) <- "RNA"
#Seurat.object_s <- subset(Seurat.object,  idents = c("LEY", "LEY2", "MYD","SRT"))
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

Seurat.object$condition = plyr::revalue(Seurat.object$condition, c("iNOA"="iGCA", "CTL"="CTL"))

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


pdf(paste(dir,'Figure3_DLK1_NOTCH2_blend.pdf',sep=''),width=19, height=7)
wrap_plots(A = vp_D, B = vp_N, C = fp_blend, design = layout)
dev.off()

im = image_read_pdf(paste(dir,"Figure3_DLK1_NOTCH2_blend.pdf",sep=''),density = 140)
image_write(im, path = paste(dir,"Figure3_DLK1_NOTCH2_blend.tiff", sep=''), format = "tiff")

######## IGF1 IGF2##########
Seurat.object = final_int_obj
DefaultAssay(Seurat.object) <- "RNA"
#Seurat.object_s <- subset(Seurat.object, idents = c("LEY","LEY2", "MYD","SRT","END","STRO"))
#Seurat.object_s <- subset(Seurat.object, idents = c("LEY", "MYD"))
Seurat.object_s <- subset(Seurat.object, idents = c("LEY", "MYD","SRT","END","STRO"))
#

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

Seurat.object$condition = plyr::revalue(Seurat.object$condition, c("iNOA"="iGCA", "CTL"="CTL"))
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
vp_IG1

layout <- '
AACCCC
BBCCCC
' 
pdf(paste(dir,'Figure3_IGFs_blend.pdf',sep=''),width=19, height=7)
wrap_plots(A = vp_IG1, B = vp_IG2, C = fp_bl, design = layout)
dev.off()

im = image_read_pdf(paste(dir,"Figure3_IGFs_blend.pdf",sep=''),density = 140)
image_write(im, path = paste(dir,"Figure3_IGFs_blend.tiff", sep=''), format = "tiff")

#### COL1A1 COL1A2 ######
Seurat.object = final_int_obj
DefaultAssay(Seurat.object) <- "RNA"
#Seurat.object_s <- subset(Seurat.object,  idents = c("LEY", "LEY2","MYD","STRO"))
Seurat.object_s <- subset(Seurat.object,  idents = c("LEY","MYD","STRO"))
#Seurat.object_s <- subset(Seurat.object,  idents = c("LEY", "MYD"))

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

Seurat.object$condition = plyr::revalue(Seurat.object$condition, c("iNOA"="iGCA", "CTL"="CTL"))
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
pdf(paste(dir,'Figure5_COL1s_blend.pdf',sep=''),width=19, height=7)
wrap_plots(A = vp_IG1, B = vp_IG2, C = fp_bl, design = layout)
dev.off()

im = image_read_pdf(paste(dir,"Figure5_COL1s_blend.pdf",sep=''),density = 140)
image_write(im, path = paste(dir,"Figure5_COL1s_blend.tiff", sep=''), format = "tiff")


######### imprinted genes ##########
relevant_imprinted = c('H19', 'MEG3', 'MEG8','IGF2', 'DLK1', 'HYMAI', 'PEG3', 'ERAP2', 'RTL1', "KCNQ1OT1")
dir.create(paste(dir,'ImprintedGenes/', sep = ''))
dir.create(paste(dir,'ImprintedGenes/VlnPlot_r', sep = ''))
Seurat.object = final_int_obj
DefaultAssay(Seurat.object) <- "RNA"
#so = subset(Seurat.object, idents = c('LEY', 'LEY2','MYD'))
so = subset(Seurat.object, idents = c('LEY','MYD'))

for (gene in relevant_imprinted) {
  so[[gene]] <- FetchData(object = so, vars = gene)
  df= so[[]]
  df$condition = plyr::revalue(df$condition, c("iNOA"="iGCA", "CTL"="CTL"))
  v <- ggviolin(df, x = "condition", y = gene, 
                fill = "condition", color = "condition", width = 1,
                add = c('jitter'), add.params = list(size =.5, shape = 1),
                #              add.params = list(size = .5, shape = 20),
                alpha = 0.6,
                facet.by = 'celltype', short.panel.labs = TRUE,
                palette = c('orange', 'dodgerblue2'), 
                error.plot = "crossbar", draw_quantiles = c(0.5), trim = T,
                title = gene, 
                panel.labs.background = list(color = "white", fill = "white", size = 0.5),
                panel.labs.font = list(color = "black", face = "bold", size = 16))
  df$condition
  v = ggpar(v, ylab = 'Expression level', xlab = '')
  v = v + stat_compare_means(aes(group = condition),label = "p.signif", size =  8, 
                             label.x.npc = "center", label.y.npc = 0.85) +
    
    theme(plot.title = element_text(color="black", size=22, face="bold.italic", hjust = 0.5),
          plot.title.position =  'panel',
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 14, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 14),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=14),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = 'right',
          panel.spacing = unit('-0.1', "lines"),
          #legend.text = c('',''),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))
  assign(paste('v','c',gene,sep = '_'),v)
  pdf(paste(dir, 'ImprintedGenes/VlnPlot_r/','vlnPlot_',gene,'.pdf', sep=''), width = 5, height = 4.5)
  print(v)
  dev.off()
  
}

pdf(paste(dir, 'ImprintedGenes/VlnPlot_r/','vlnPlot_maternal.pdf', sep=''), width = 12, height = 4)
print(v_c_H19 + v_c_MEG3 + v_c_MEG8 + plot_layout(guides = "collect", ncol = 3) & 
        theme(legend.position = "right", legend.direction = 'vertical'))
dev.off()

pdf(paste(dir, 'ImprintedGenes/VlnPlot_r/','vlnPlot_paternal.pdf', sep=''), width = 28, height = 4)
print(v_c_IGF2 + v_c_DLK1 + v_c_HYMAI + v_c_PEG3 + v_c_ERAP2 + v_c_RTL1 + v_c_KCNQ1OT1 + plot_layout(guides = "collect", ncol = 7) & 
        theme(legend.position = "right", legend.direction = 'vertical'))
dev.off()


################ INSL3 #################
Seraut.object= final_int_obj

CELL_INSL3 <- WhichCells(Seraut.object, slot = 'counts', expression = INSL3 > 1)
somatic.integrated.new.INSL3 <- subset(Seraut.object, cells = CELL_INSL3)


legendiNOA='iGCA'
legendCTRL='CTL'

Seraut.object = somatic.integrated.new.INSL3
Idents(Seraut.object) <- "celltype"

cell_type.2.use = 'LEY'
#cell_type.2.use = 'LEY2'
cell_type = c(paste(cell_type.2.use,'_iNOA',sep=''),
              paste(cell_type.2.use,'_CTL',sep=''))
#cell_type = c(paste(cell_type.2.use,'_iNOA',sep=''),
#              paste(cell_type.2.use,'_OA',sep=''))
features = c('INSL3')


Cell.subset <- subset(Seraut.object, idents =  cell_type.2.use, invert=F)
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
                                          levels = c(paste(cell_type.2.use,"_iNOA",sep=''),paste(cell_type.2.use,"_CTL",sep='')))
#paste(cell_type.2.use,"_CTL",sep='')))

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
assign('INSL3_INSL3cell',CVP)


Seraut.object= final_int_obj
Idents(Seraut.object) <- "celltype"

cell_type = c(paste(cell_type.2.use,'_iNOA',sep=''),
              paste(cell_type.2.use,'_CTL',sep=''))
features = c('INSL3')


Cell.subset <- subset(Seraut.object, idents =  cell_type.2.use, invert=F)
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
                                          levels = c(paste(cell_type.2.use,"_iNOA",sep=''),paste(cell_type.2.use,"_CTL",sep='')))
#paste(cell_type.2.use,"_CTL",sep='')))

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

INSL3_allcell + plot_spacer() + INSL3_INSL3cell
p1 =INSL3_allcell; p2 =INSL3_INSL3cell;
(p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
  (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt")))
pdf(paste(dir,'Figure_INSL3.pdf',sep=''),width=8.5, height=3.5)
plot((p1 + theme(plot.margin = unit(c(0,40,0,15), "pt"))) +
       (p2 + theme(plot.margin = unit(c(0,15,0,40), "pt"))))
dev.off()

im = image_read_pdf(paste(dir,'Figure_INSL3.pdf',sep=''),density = 140)
image_write(im, path = paste(dir,'Figure_INSL3.tiff', sep=''), format = "tiff")


#### Pathway UP heatmaps ####
Seurat.object = final_int_obj
Seurat.object$cell_type = Seurat.object$celltype
Idents(Seurat.object) <- "celltype.cond"
DefaultAssay(Seurat.object) <- "RNA"

dir

CND = 'up'
database = c('Reactome_2016','GO_Biological_Process_2018')
indir =  paste(dir,'enrichR/', sep = '')
indir
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

pathway_to_plot[["LEY2"]] <- c("Cholesterol biosynthesis Homo sapiens R-HSA-191273",
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

indir = paste(dir, "enrichR/", sep='')
outdir = paste(dir, 'Pathways/', sep='') 
dir.create(outdir)

for (i in cell.type_2_plot) {
  print(i)
  ds = Seurat.object
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

######## Dotplots ##########à

pathways2plot = pathway_to_plot
pathways.dataframe = data.frame(cellType = character(),
                                Pathway = character(),
                                gene.ratio = numeric(),
                                p.value = numeric(),
                                p.value.adj = numeric(),
                                stringsAsFactors=FALSE)

all.cell.types= cell.type_2_plot

fx <- function(x) eval(parse(text=enrichR.table$Overlap[x]))
fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))
database = c('Reactome_2016','GO_Biological_Process_2018')
for (cell_type in all.cell.types) {
  print(cell_type)
  enrichR.file = paste(indir,"enrichR_DGE_",cell_type,'_',CND,'.xlsx',sep='')
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
    ylim(c(0,11)) +
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


######### Downregulated pathways #######
Seurat.object = final_int_obj

Seurat.object$cell_type = Seurat.object$celltype
Idents(Seurat.object) <- "celltype.cond"
DefaultAssay(Seurat.object) <- "RNA"

CND = 'down'
#cell_types= c("LEY", "LEY2", "MYD", "SRT","END", "STRO") # "MCR"
cell_types= c("LEY", "MYD", "SRT","END", "STRO")

database = c('Reactome_2016','GO_Biological_Process_2018')
gene.l = list()
path.l = list()
gene.list = character()
indir = paste(dir, "enrichR/", sep='')
outdir = paste(dir, "Pathways/", sep='')
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



order.2.plot = c('LEY_iNOA', 'LEY_CTL', 'LEY2_iNOA', 'LEY2_CTL',
                 'MorrisOnly_iNOA', 'MorrisOnly_CTL',
                 'MYD_iNOA','MYD_CTL','SRT_iNOA','SRT_CTL',
                 'END_iNOA','END_CTL','MCR_iNOA','MCR_CTL','TCL_iNOA', 'TCL_CTL',
                 'STRO_iNOA','STRO_CTL','UND_iNOA','UND_CTL')


order.2.plot = c('LEY_iNOA', 'LEY_CTL', 
                 'MorrisOnly_iNOA', 'MorrisOnly_CTL',
                 'MYD_iNOA','MYD_CTL','SRT_iNOA','SRT_CTL',
                 'END_iNOA','END_CTL','MCR_iNOA','MCR_CTL','TCL_iNOA', 'TCL_CTL',
                 'STRO_iNOA','STRO_CTL','UND_iNOA','UND_CTL')

#order.2.plot = c('LEY_iNOA', 'LEY_CTL', 'MYD_iNOA','MYD_CTL','SRT_iNOA','SRT_CTL',
#              'END_iNOA','END_CTL','MCR_iNOA','MCR_CTL','TCL_iNOA','TCL_CTL',
#             'STRO_iNOA','STRO_CTL','UND_iNOA','UND_CTL')

#order.2.plot = c('LEY_iNOA', 'LEY_OA', 'MYD_iNOA','MYD_OA','SRT_iNOA','SRT_OA')

levels(Seurat.object) <- order.2.plot
DefaultAssay(Seurat.object) <- "RNA"

#markers.to.plot <- AgeGenesClustersTable_dw$gene.name[AgeGenesClustersTable_dw$cluster==3]
dp_down = DotPlot(Seurat.object, 
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
p

pathways.dataframe = data.frame(cellType = character(),
                                Pathway = character(),
                                gene.ratio = numeric(),
                                p.value = numeric(),
                                p.value.adj = numeric(),
                                stringsAsFactors=FALSE)

all.cell.types= cell_types
database = c('Reactome_2016','GO_Biological_Process_2018')
N=length(p)
CND = 'DOWN'
fx <- function(x) eval(parse(text=enrichR.table$Overlap[x]))
fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))

for (cell_type in all.cell.types) {
  print(cell_type)
  enrichR.file = enrichR.file = paste(indir,'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
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

col = as.character(colors_Alf[names(gene.l),]$cols)

venn::venn(path.l, simplify=TRUE, opacity = 0.3, box = FALSE, ilab=TRUE, zcolor = col, ilcs = 1.5, sncs = 2)

pdf(paste(outdir,'Figure7_venn.pdf',sep=''),width=7, height=6.5)
print(venn::venn(path.l, simplify=TRUE, opacity = 0.3, box = FALSE, elipse = T, ilab=TRUE, zcolor = col, ilcs = 1.5, sncs = 2.5))
dev.off()
im = image_read_pdf(paste(outdir,"Figure7_venn.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure7_venn.pdf.tiff", sep=''), format = "tiff")




