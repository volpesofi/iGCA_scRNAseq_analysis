#' plots for paper 


######### libraries ##########
suppressMessages(library(plyr))
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
suppressMessages(library(enrichR))
suppressMessages(library(ggpubr))
suppressMessages(library(ggpubr))
suppressMessages(library("svglite"))

###### load data ########
datasets_dir = ('/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/')
OBJ_dir = '/Users/tascini.annasofia/Dropbox (HSR Global)/WORKSPACE/Alfano/AlfanoM_904_infertilita_epigenetics/7_bioinfo/SeuratObjects/'
SO_dir = '/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/Seurat_objects/'
outdir = '/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/pannelli_paper/'
outdir = '/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/Adult_iNOA_pannelli_perMassimo/'
dir.create(outdir, recursive = TRUE)

######## my functions #######
source("mySeuratfunctions.R")
source('./utility_functions.R')
source('./integration_function.R')

###### color #######
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



######## Load adult testis #########
load(paste(SO_dir, "cleanAdultIntegration", sep = ''))
cleanAdultIntegration
#Final_GSE124263_GSE112013_iNOA.integrated = Seurat.object
#Final_GSE124263_GSE112013_iNOA.integrated

Seurat.object = cleanAdultIntegration
table(Seurat.object$condition)
Seurat.object$condition = revalue(Seurat.object$condition, c("iNOA"="iGCA", "CTL"="CTL"))
Seurat.object$condition = factor(Seurat.object$condition, levels = c("iGCA", "CTL"))
Seurat.object$celltype.cond <- paste(Idents(Seurat.object), Seurat.object$condition, sep = "_")
Idents(Seurat.object) <- "celltype.cond"
order.2.plot = c("LEY_iGCA", "LEY_CTL", 
                 "MYD_iGCA", "MYD_CTL",
                 "SRT_iGCA", "SRT_CTL",
                 "END_iGCA","END_CTL",
                 "MCR_iGCA","MCR_CTL",
                 "TCL_iGCA",
                 "STRO_iGCA", "STRO_CTL",
                 "UND_iGCA", "UND_CTL")

levels(Seurat.object) <- order.2.plot
Seurat.object$celltype.cond = factor(Seurat.object$celltype.cond, levels = order.2.plot)


#---------------- UMAP split ----------------
Idents(Seurat.object) <- "celltype"
DefaultAssay(Seurat.object) <- "RNA"
col = as.character(colors_Alf[levels(Seurat.object),]$cols)
p.inf <- DimPlot(Seurat.object, 
                 reduction = "umap", 
                 label = T, 
                 order = F,
                 pt.size = 2, 
                 split.by = 'condition') +
  scale_color_manual(values = (col)) +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2") + NoLegend()

p.inf
df_dir = "/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/Dataframes_Figures/"
dir.create(df_dir)
xy_fig1G = Seurat.object[["umap"]]@cell.embeddings
write.xlsx(xy_fig1G, paste(df_dir, "Figure1_G.xlsx",sep =''), 
           row.names = T)

pdf(paste(outdir,"Figure1G.umap.pdf",sep=''),  width=12, height=6)
p.inf
dev.off()

ggsave(filename = paste(outdir,"Figure1G.umap.png",sep=''),p.inf, width=12, height=6)


#---------------- HM and MG ----------------
DefaultAssay(Seurat.object) <- "integrated"
Idents(Seurat.object) <- "celltype"
#cluster.markers.c = FindConservedMarkers(Seurat.object, grouping.var ="source", ident.1 = "TCL")
cluster.markers = FindAllMarkers(Seurat.object)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.xlsx(cluster.markers,
           file= paste(dataset_iNOA_dir, filename, sep =''), 
           row.names = T,
           asTable = T)
HM <- DoHeatmap(Seurat.object, 
                features = top10$gene, 
                group.colors = col,
                disp.min = -1.5,
                disp.max = 1.5,
                angle = 90) +
  scale_fill_gradientn(colours = coolwarm(200)) 

ggsave(paste(outdir, "SI/", "HM_integrated_IGCA_CTL.png", sep=''), HM, height = 15, width = 15)

# ----------- FeaturePlot of MF signature SI ----------------
DefaultAssay(Seurat.object) <- "RNA"
Idents(Seurat.object) <- "celltype"

DefaultAssay(Seurat.object) <- 'RNA'
pL <- Plot_sign(Seurat.object,
                signature= c('CFD','DLK1','LUM'), 
                operator = mean, title = 'LEY')
pM <- Plot_sign(Seurat.object,
                signature= c('ACTA2','MYH11','DES'), 
                operator = mean, title = 'MYD')
pS <- Plot_sign(Seurat.object,
                signature= c('FATE1','CITED1','SOX9', 'AMH'), 
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
EEHHGGDD
EEHHGG##
EEHHGG##
'
png(paste(outdir, "SI/","MG_integration_iGCA_CTL.png", sep=''),  width=1400, height=600)
print(wrap_plots(A = pL, B = pM, C = pS, D = pSTRO, E = pE, H = pMa, G =pT, design = layout))
dev.off()

pdf(paste(outdir, "SI/","MG_integration_iGCA_CTL.pdf", sep=''),  width=14, height=6)
print(wrap_plots(A = pL, B = pM, C = pS, D = pSTRO, E = pE, H = pMa, G =pT, design = layout))
dev.off()

##### EGR3 FOSB ##########

vp_IG1 <- VlnPlot(Seurat.object,
                  feature = c('FOSB'),
                  split.plot = TRUE,
                  pt.size = 0.2,
                  split.by = 'condition') +
theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
      plot.subtitle = element_text(color="black", size=16, face="italic"),
      axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
      axis.title.x = element_text(face = "bold", color = "black", size = 18),
      axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
      axis.title.y = element_text(face = "bold", color = "black", size = 18),
      legend.text = element_text(face = "bold", color = "black", size = 12),
      legend.position = "top",
      panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))
vp_IG1

fp_bl <- FeaturePlot(Seurat.object, 
                     reduction = "umap", 
                     pt.size = 1,
                     features = c('JUN'), 
                     label = T, 
                     order= T,
                     blend = F, 
                     repel = T,
                     cols = c('lightgrey','red'),
                     max.cutoff = 'q99',
                     split.by='condition') + theme(legend.position = 'right')

fp_bl
########## DLK1 feature ###########
DefaultAssay(Seurat.object) <- "RNA"
Idents(Seurat.object) <- "celltype"
Seurat.object_s <- subset(Seurat.object, idents = c("LEY", "MYD","SRT"))
vp_IG1 <- VlnPlot(Seurat.object_s,
                  feature = c('DLK1'),
                  split.plot = TRUE,
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
        legend.position = "top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))


vp_IG2 <- VlnPlot(Seurat.object_s,
                  feature = c('NOTCH2'),
                  pt.size = 0,
                  split.plot = TRUE,
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
        legend.position = "top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))

fig4Av1 = vp_IG1$data
fig4Av1 = fig4Av1[order(fig4Av1$ident),]
fig4Av1$cellID = rownames(fig4Av1)
fig4Av2 = vp_IG2$data
fig4Av2 = fig4Av2[order(fig4Av2$ident),]
fig4Av2$cellID = rownames(fig4Av2)
shcol = intersect(colnames(fig4Av1), colnames(fig4Av2))
fig4A = merge(fig4Av1, fig4Av2, by = shcol)
write.xlsx(fig4A,paste(df_dir,"Figure4_A.xlsx",sep=''),row.names = F)


fp_bl <- FeaturePlot(Seurat.object, 
                     reduction = "umap", 
                     pt.size = 1,
                     features = c('DLK1','NOTCH2'), 
                     label = T, 
                     order= T,
                     blend = T,       
                     max.cutoff = 'q99',
                     split.by='condition') + theme(legend.position = 'right')


layout <- '
AACCCC
BBCCCC
' 
pdf(paste(outdir,'Figure4_DLK1_NOTCH2_blend.pdf',sep=''),width=19, height=7)
wrap_plots(A = vp_IG1, B = vp_IG2, C = fp_bl, design = layout)
dev.off()

im = image_read_pdf(paste(outdir,"Figure4_DLK1_NOTCH2_blend.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure4_DLK1_NOTCH2_blend.tiff", sep=''), format = "tiff")


########## IGF2 feature ###########
DefaultAssay(Seurat.object) <- "RNA"
Seurat.object_s <- subset(Seurat.object, idents = c("LEY", "MYD","SRT"))
vp_IG1 <- VlnPlot(Seurat.object_s,
                  feature = c('IGF1'),
                  split.plot = TRUE,
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
        legend.position = "top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))

vp_IG2 <- VlnPlot(Seurat.object_s,
                  feature = c('IGF2'),
                  pt.size = 0,
                  split.plot = TRUE,
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
        legend.position = "top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))

fig4Av1 = vp_IG1$data
fig4Av1 = fig4Av1[order(fig4Av1$ident),]
fig4Av1$cellID = rownames(fig4Av1)
fig4Av2 = vp_IG2$data
fig4Av2 = fig4Av2[order(fig4Av2$ident),]
fig4Av2$cellID = rownames(fig4Av2)
shcol = intersect(colnames(fig4Av1), colnames(fig4Av2))
fig4A = merge(fig4Av1, fig4Av2, by = shcol)
write.xlsx(fig4A,paste(df_dir,"Figure4_B.xlsx",sep=''),row.names = F)



fp_bl <- FeaturePlot(Seurat.object, 
                     reduction = "umap", 
                     pt.size = 1,
                     features = c('IGF1','IGF2'), 
                     label = T, 
                     order= T,
                     blend = T,       
                     max.cutoff = 'q99',
                     split.by='condition') + theme(legend.position = 'right')

fp_bl
vp_IG1

layout <- '
AACCCC
BBCCCC
' 
pdf(paste(outdir,'Figure4_IGFs_blend.pdf',sep=''),width=19, height=7)
wrap_plots(A = vp_IG1, B = vp_IG2, C = fp_bl, design = layout)
dev.off()

im = image_read_pdf(paste(outdir,"Figure4_IGFs_blend.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure4_IGFs_blend.tiff", sep=''), format = "tiff")

# ---------- HSD17B3 -------------
DefaultAssay(Seurat.object) <- "RNA"
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

igca = fp_bl[[1]]$data
igca$cellID = row.names(igca)
write.xlsx(igca[,c("cellID", "UMAP_1", "UMAP_2","ident", "HSD17B3")],
           paste(df_dir,"Figure4_Dl.xlsx",sep=''),row.names = F)
igca = fp_bl[[2]]$data
igca$cellID = row.names(igca)
write.xlsx(igca[,c("cellID", "UMAP_1", "UMAP_2","ident", "HSD17B3")],
           paste(df_dir,"Figure4_Dr.xlsx",sep=''),row.names = F)


pdf(paste(outdir,'Figure4_HSD17B3.pdf',sep=''),width=8, height=4)
plot(fp_bl)
dev.off()

im = image_read_pdf(paste(outdir,'Figure4_HSD17B3.pdf',sep=''),density = 140)
image_write(im, path = paste(outdir,'Figure4_HSD17B3.tiff', sep=''), format = "tiff")

# ------------ figure4 statistics --------------

ct = "LEY"; c1 = "iGCA"; c2 = "CTL"
features = c("DLK1", "NOTCH2", "IGF2", "IGF1", "HSD17B3")
Idents(Seurat.object) <- "celltype.cond"
dge_value_for_figure = data.frame()
for (ct in levels(Seurat.object$celltype)[levels(Seurat.object$celltype)!="TCL"]) {
  print(ct)
  value <- FindMarkers(Seurat.object,
                       ident.1 = paste(ct,"_",c1, sep=''),
                       ident.2 = paste(ct,"_",c2, sep=''), 
                       verbose = FALSE,
                       feature = features,
                       logfc.threshold=0,
                       min.pct = 0)
  value$celltype = ct
  value$gene = row.names(value)
  dge_value_for_figure = rbind(dge_value_for_figure, value)
}
dge_value_for_figure = dge_value_for_figure[order(dge_value_for_figure$gene),
                     c("gene", "celltype", "avg_logFC", "p_val", "p_val_adj", "pct.1","pct.2")]

write.xlsx(dge_value_for_figure, file = paste(outdir, "Figure4_Genes_statistics.xlsx", sep=''), 
           row.names = FALSE, asTable = T)


################ INSL3 feature #################
Idents(Seurat.object) <- "celltype.cond"
CELL_INSL3 <- WhichCells(Seurat.object, slot = 'counts', expression = INSL3 > 1)
somatic.integrated.new.INSL3 <- subset(Seurat.object, cells = CELL_INSL3)

legendiNOA='iGCA'
legendCTRL='CTL'
Seurat.object.insl3 = somatic.integrated.new.INSL3

cell_type.2.use = 'LEY'
cell_type = c(paste(cell_type.2.use,'_iGCA',sep=''),
              paste(cell_type.2.use,'_CTL',sep=''))
features = c('INSL3')


Idents(Seurat.object.insl3) <- "celltype"
Cell.subset <- subset(Seurat.object.insl3, idents =  cell_type.2.use, invert=F)
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
                                           condition = Seurat.object.insl3[[]][names(gene.counts.tmp),
                                                                         'celltype.cond'],
                                           counts = (gene.counts.tmp))
  gene.counts.dataframe <- rbind(gene.counts.dataframe, gene.counts.dataframe.entry)
}
gene.counts.dataframe$condition <- factor(x = gene.counts.dataframe$condition, 
                                          levels = c(paste(cell_type.2.use,"_iGCA",sep=''),paste(cell_type.2.use,"_CTL",sep='')))
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
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c("iGCA", "CTL")) +
  labs(title = cell_type.2.use , 
       x = "Genes", 
       y = "Expression Level") 
CVP = CVP +labs(title= 'Cells expr. INSL3',
                x = "LEY", 
                y = "Expression Level") 
assign('INSL3_INSL3cell',CVP)
CVP
CVP$data
CVP$data$ident = "LEY" 
CVP$data$condition = revalue(CVP$data$condition, 
                             c("LEY_iGCA"="iGCA", "LEY_CTL"="CTL"))
write.xlsx(CVP$data[,c("ident", "condition", "cell_barcode", "counts")],
           paste(df_dir,"Figure4_C.xlsx",sep=''),row.names = F)

# ------ INSL3 all cell --------

Idents(Seurat.object) <- "celltype.cond"

cell_type.2.use = 'LEY'
cell_type = c(paste(cell_type.2.use,'_iGCA',sep=''),
              paste(cell_type.2.use,'_CTL',sep=''))
features = c('INSL3'); ct = 'LEY'; c1='iGCA'; c2= 'CTL'
INS3_all <- FindMarkers(Seurat.object,
                        ident.1 = paste(ct,"_",c1, sep=''),
                        ident.2 = paste(ct,"_",c2, sep=''), 
                        verbose = FALSE,
                        feature = features,
                        logfc.threshold=0,
                        min.pct = 0)
Idents(Seurat.object) <- "celltype"
Cell.subset <- subset(Seurat.object, idents =  cell_type.2.use, invert=F)
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
                                           condition = Seurat.object[[]][names(gene.counts.tmp),
                                                                         'celltype.cond'],
                                           counts = (gene.counts.tmp))
  gene.counts.dataframe <- rbind(gene.counts.dataframe, gene.counts.dataframe.entry)
}
gene.counts.dataframe$condition <- factor(x = gene.counts.dataframe$condition, 
                                          levels = c(paste(cell_type.2.use,"_iGCA",sep=''),paste(cell_type.2.use,"_CTL",sep='')))
#paste(cell_type.2.use,"_CTL",sep='')))

CVP = ggplot(gene.counts.dataframe, aes(x=gene, y = counts, fill=condition)) +
  geom_split_violin() + 
  scale_x_discrete('LEY') +
  scale_fill_manual(values = c('purple', 'orange')) +
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
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c("iGCA", "CTL")) +
  labs(title = cell_type.2.use , 
       x = "Genes", 
       y = "Expression Level") 
CVP = CVP +labs(title= 'All cells',
                x = "LEY", 
                y = "Expression Level") 
CVP

assign('INSL3_allcell',CVP) 

INSL3_allcell + plot_spacer() + INSL3_INSL3cell
p1 =INSL3_allcell; p2 =INSL3_INSL3cell;
(p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
  (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt")))
pdf(paste(outdir,'Figure4_INSL3.pdf',sep=''),width=8.5, height=3.5)
plot((p1 + theme(plot.margin = unit(c(0,40,0,15), "pt"))) +
       (p2 + theme(plot.margin = unit(c(0,15,0,40), "pt"))))
dev.off()

im = image_read_pdf(paste(outdir,'Figure4_INSL3.pdf',sep=''),density = 140)
image_write(im, path = paste(outdir,'Figure4_INSL3.tiff', sep=''), format = "tiff")

ct = "LEY"; c1 = "iGCA"; c2 = "CTL"
Idents(Seurat.object) <- "celltype.cond"
INS3_exp <- FindMarkers(somatic.integrated.new.INSL3,
                        ident.1 = paste(ct,"_",c1, sep=''),
                        ident.2 = paste(ct,"_",c2, sep=''), 
                        verbose = FALSE,
                        feature = features,
                        logfc.threshold=0,
                        min.pct = 0)
INS3_all <- FindMarkers(Seurat.object,
                        ident.1 = paste(ct,"_",c1, sep=''),
                        ident.2 = paste(ct,"_",c2, sep=''), 
                        verbose = FALSE,
                        feature = features,
                        logfc.threshold=0,
                        min.pct = 0)


INSL3 = rbind(INS3_all, INS3_exp)
INSL3$cells = c("all", "expr. INS3")
write.xlsx(INSL3, file = paste(outdir, "ISNL3_statistics.xlsx", sep=''), row.names = TRUE)

####### % cells #########
A <- as.data.frame(prop.table(table(Seurat.object$celltype, 
                                    Seurat.object$condition),margin = 2))
A2 <- as.data.frame(table(Seurat.object$celltype, 
                          Seurat.object$condition),margin = 2)
A$N <- A2$Freq

col = as.character(colors_Alf[levels(Seurat.object$celltype),]$cols)
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
ggsave(filename = paste(outdir, "Figure1_H.png", sep = ''), N, width = 6, height = 4)

########## Aging genes ##########
PD="/Users/tascini.annasofia/Dropbox (HSR Global)/Alfano_904_paperDraft/"
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

markers.to.plot = AgeGenesClustersTable_dw$gene.name[AgeGenesClustersTable_dw$cluster==3]

Idents(Seurat.object) <- "celltype.cond"
DefaultAssay(Seurat.object) <- "RNA"

dp_RPL = DotPlot(Seurat.object, 
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

write.xlsx(dp_RPL$data[,c("features.plot","id", "avg.exp", "avg.exp.scaled", "pct.exp")], 
           paste(df_dir,"Fig6_dotplot.xlsx",sep=''),row.names = F)


pdf(paste(outdir,'Figure6_DotPlot_RPLgenes.pdf',sep=''),width=14, height=7)
plot(dp_RPL)
dev.off()

im = image_read_pdf(paste(outdir,"Figure6_DotPlot_RPLgenes.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure6_DotPlot_RPLgene.tiff", sep=''), format = "tiff")


######### Pathways ##########
###### pathway #######
PGenes = list()
PPath = list()
PGenes2 = list()
PPath2 = list()

gene_pathway = function(pathway, enrichRfolder, CND) {
  PGenes = list()
  cellTypes = names(table(Seurat.object$celltype))
  cellTypes = cellTypes[cellTypes !="TCL"]
  cellTypes = cellTypes[cellTypes !="UND"]
  for (celltype in cellTypes) {
    print(paste("searching", celltype))
    enrich_file = paste(enrichRfolder,'enrichR_DGE_',celltype,'_',CND,'.xlsx', sep='')
    #  
    enrichR.table = data.frame()
    for (dat in getSheetNames(enrich_file)) {
      Table <- read.xlsx(xlsxFile = enrich_file, 
                         sheet = dat, 
                         startRow = 1, 
                         colNames = TRUE,
                         rowNames = F, 
                         detectDates = FALSE, 
                         skipEmptyRows = TRUE,
                         skipEmptyCols = TRUE,
                         na.strings = "NA", 
                         fillMergedCells = FALSE)
      Table$database = dat
      enrichR.table = rbind(enrichR.table , Table)
    }
    
    Seleno_pathways = enrichR.table[enrichR.table$Term %in% pathway,]
    Genes = sort(unique(unlist(strsplit(Seleno_pathways$Genes,';'))))
    PGenes[[celltype]] = Genes
  }
  return(PGenes)
}

###### mit respiratory chain #######
path_file = paste("../Final_GSE124263_GSE112013_iNOA/enrichR/",sep='')
#path_file = paste("../cleanAdultIntegration/enrichR/",sep='')
pathway = 'mitochondrial respiratory chain complex assembly (GO:0033108)'
CND = 'down'
MITrespiratory_list = gene_pathway(pathway = pathway, enrichRfolder = path_file, CND = CND)
MITrespiratory_genes = sort(unique(unlist(MITrespiratory_list)))
markers.to.plot <- MITrespiratory_genes
markers.to.plot

dp1 = DotPlot(Seurat.object, 
              features = rev(markers.to.plot), 
              cols = c("lightblue", "red"), 
              dot.scale = 8) +
  RotatedAxis() +
  scale_size_continuous(range = c(1,8)) +
  theme(axis.text.x = element_text(angle = 90,  color = "black", size=20, vjust =.5, hjust =1, family = 'mono'), 
        axis.title.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=18),
        axis.title.y = element_text(face = "bold", color = "black", size = 20),
        legend.text = element_text(face = "bold", color = "black", size = 16),
        legend.position="right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "Genes", y = "Cell type") + NoLegend()
dp1

pdf(paste(outdir, 'Dot_mitResp.pdf', sep =''), width = 7, height = 6)
dp1
dev.off()

write.xlsx(dp1$data[,c("features.plot","id", "avg.exp", "avg.exp.scaled", "pct.exp")], 
           paste(df_dir,"Mith_dotplot.xlsx",sep=''),row.names = F)

###### HIF #######
pathway = 'Oxygen-dependent proline hydroxylation of Hypoxia-inducible Factor Alpha Homo sapiens R-HSA-1234176'
pathway = 'Regulation of Hypoxia-inducible Factor (HIF) by oxygen Homo sapiens R-HSA-1234174'
CND = 'down'
HIF_list = gene_pathway(pathway = pathway, enrichRfolder = path_file, CND = CND)
HIF_genes = sort(unique(unlist(HIF_list)))
markers.to.plot <- HIF_genes
markers.to.plot

dp2 = DotPlot(Seurat.object, 
              features = rev(markers.to.plot), 
              cols =c("lightblue", "red"),
              dot.scale = 8) +
  RotatedAxis() +
  scale_size_continuous(range = c(1,8)) +
  theme(axis.text.x = element_text(angle = 90,  color = "black", size=20, vjust =.5, hjust =1, family = 'mono'), 
        axis.title.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.text = element_text(face = "bold", color = "black", size = 16),
        legend.position="right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "Genes", y = "Cell type") 

pdf(paste(outdir, 'Dot_HIF.pdf', sep =''), width = 4, height = 6)
dp2
dev.off()

write.xlsx(dp2$data[,c("features.plot","id", "avg.exp", "avg.exp.scaled", "pct.exp")], 
           paste(df_dir,"HIF_dotplot.xlsx",sep=''),row.names = F)

# ----------- down pathway statistics ------------
features = c(MITrespiratory_genes, HIF_genes)
c1 = "iGCA"; c2 = "CTL"
Idents(Seurat.object) <- "celltype.cond"
dge_value_for_figure = data.frame()
for (ct in levels(Seurat.object$celltype)[levels(Seurat.object$celltype)!="TCL"]) {
  print(ct)
  value <- FindMarkers(Seurat.object,
                       ident.1 = paste(ct,"_",c1, sep=''),
                       ident.2 = paste(ct,"_",c2, sep=''), 
                       verbose = FALSE,
                       feature = features,
                       logfc.threshold=0,
                       min.pct = 0)
  value$celltype = ct
  value$gene = row.names(value)
  dge_value_for_figure = rbind(dge_value_for_figure, value)
}

dge_value_for_figure = dge_value_for_figure[order(dge_value_for_figure$gene),
                                            c("gene", "celltype", "avg_logFC", "p_val", "p_val_adj", "pct.1","pct.2")]

write.xlsx(dge_value_for_figure, file = paste(outdir, "Figure3_HIF_MIT_Genes_statistics.xlsx", sep=''), 
           row.names = FALSE, asTable = T)



###### Selenoprotein #######
SEL = row.names(Seurat.object)[grep('^sel', row.names(Seurat.object), ignore.case = T)]
SELENO = row.names(Seurat.object)[grep('^seleno', row.names(Seurat.object), ignore.case = T)]
SELENO
cell_type <- c( "LEY",
                "MYD", 
                "END",
                "MCR",
                "STRO",
                "SRT")

#cell_type <- c( "LEY",
#               "MYD", 
#                "SRT",
#                "UND")

#COL = row.names(Seurat.object)[grep("^COL",row.names(Seurat.object))]
gene = c(SELENO, c("INMT","GPX1","GPX3","GPX4"))
DefaultAssay(Seurat.object) <-"RNA"
Idents(Seurat.object) <- "celltype.cond"
order.2.plot = c("LEY_iGCA", "LEY_CTL", 
                 "MYD_iGCA", "MYD_CTL",
                 "SRT_iGCA", "SRT_CTL",
                 "END_iGCA","END_CTL",
                 "MCR_iGCA","MCR_CTL",
                 "TCL_iGCA",
                 "STRO_iGCA", "STRO_CTL",
                 "UND_iGCA", "UND_CTL")

levels(Seurat.object) <- order.2.plot
DefaultAssay(Seurat.object) <- "RNA"
azoospermia.modulated.gene = data.frame()
for (ct in cell_type) {
  #print(ct)
  azoospermia.response.gene <- data.frame()
  azoospermia.response.gene <- FindMarkers(Seurat.object,
                                           ident.1 = paste(ct,"_iGCA", sep=''),
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
           file= paste(outdir,'seleno.xlsx',sep=''), 
           row.names = FALSE,
           asTable = TRUE)

genes = c("INMT","GPX1","GPX3","GPX4",'SELENOP', 'SELENBP1', 'SELENOH')
genes = c(SELENO, c("INMT","GPX1","GPX3","GPX4"))
genes = genes[genes != "SELENOOLP"]
#genes
markers.to.plot <- genes 

Idents(Seurat.object) <- "celltype.cond"
order.2.plot = c("LEY_iGCA", "LEY_CTL", 
                 "MYD_iGCA", "MYD_CTL",
                 "SRT_iGCA", "SRT_CTL",
                 "END_iGCA","END_CTL",
                 "MCR_iGCA","MCR_CTL",
                 "TCL_iGCA",
                 "STRO_iGCA", "STRO_CTL",
                 "UND_iGCA", "UND_CTL")

levels(Seurat.object) <- order.2.plot

dp3 = DotPlot(Seurat.object, 
              features = rev(markers.to.plot), 
              cols =c("lightblue", "red"),
              dot.scale = 8) +
  #split.by ='condition') +
  RotatedAxis() +
  #coord_flip() +
  scale_size_continuous(range = c(1,8)) +
  theme(axis.text.x = element_text(angle = 90,  color = "black", size=20, vjust =.5, hjust =1, family = 'mono'), 
        axis.title.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=18),
        axis.title.y = element_text(face = "bold", color = "black", size = 20),
        legend.text = element_text(face = "bold", color = "black", size = 16),
        legend.position="right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "Genes", y = "Cell type") 

pdf(paste(outdir, 'Figure3_Dot_SP.pdf', sep =''), width = 8, height = 6)
dp3
dev.off()

write.xlsx(dp3$data[,c("features.plot","id", "avg.exp", "avg.exp.scaled", "pct.exp")], 
           paste(df_dir,"Redoxin_dotplot.xlsx",sep=''),row.names = F)


features = genes
Idents(Seurat.object) <- "celltype.cond"
dge_value_for_figure = data.frame()
c1 = 'iGCA'; c2='CTL'
for (ct in levels(Seurat.object$celltype)[levels(Seurat.object$celltype)!="TCL"]) {
  print(ct)
  value <- FindMarkers(Seurat.object,
                       ident.1 = paste(ct,"_",c1, sep=''),
                       ident.2 = paste(ct,"_",c2, sep=''), 
                       verbose = FALSE,
                       feature = features,
                       logfc.threshold=0,
                       min.pct = 0)
  value$celltype = ct
  value$gene = row.names(value)
  dge_value_for_figure = rbind(dge_value_for_figure, value)
}
dge_value_for_figure = dge_value_for_figure[order(dge_value_for_figure$gene),
                                            c("gene", "celltype", "avg_logFC", "p_val", "p_val_adj", "pct.1","pct.2")]

write.xlsx(dge_value_for_figure, file = paste(outdir, "SI/", "SELENOproteins_statistics.xlsx", sep=''), 
           row.names = FALSE, asTable = T)

# GPX1
Idents(Seurat.object) <- "celltype"
Seurat.object_s <- subset(Seurat.object, idents = c("MCR", "SRT"))
g = 'GPX1'
vp <- VlnPlot(Seurat.object_s,
              feature = c(g),
              split.plot = TRUE,
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
        legend.position = "top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))
ggsave(paste(outdir, g, "_VP.png", sep=''), vp, height = 4, width = 5)
  
Idents(Seurat.object) <- "celltype"
Seurat.object_s <- subset(Seurat.object, idents = c("LEY"))
g = 'INMT'
vp <- VlnPlot(Seurat.object_s,
              feature = c(g),
              split.plot = TRUE,
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
        legend.position = "top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))
ggsave(paste(outdir, g, "_VP.png", sep=''), vp, height = 4, width = 4)


Idents(Seurat.object) <- "celltype"
Seurat.object_s <- subset(Seurat.object, idents = c("LEY",
                                                    "MYD",
                                                    "END",
                                                    "MCR",

                                                                                                        "SRT"))
for (g  in SELENO) {
  vp <- VlnPlot(Seurat.object_s,
                feature = c(g),
                split.plot = TRUE,
                pt.size = 0.1,
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
          legend.position = "top",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))
  vp
  ggsave(paste(outdir, g, "_VP.png", sep=''), vp, height = 4, width = 10)
  }


########## Pathways ###########

#### Pathway UP heatmaps ####
Seurat.object$cell_type = Seurat.object$celltype
Idents(Seurat.object) <- "celltype.cond"
DefaultAssay(Seurat.object) <- "RNA"

dir = '/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/cleanAdultIntegration/'

CND = 'up'
database = c('Reactome_2016','GO_Biological_Process_2018')
indir =  paste(dir,'enrichR/', sep = '')
indir
cell.type_2_plot = c("LEY", "MYD", "STRO", "MCR", "SRT")

pathway_to_plot <- list()
pathway_to_plot[["LEY"]] <- c("nuclear-transcribed mRNA catabolic process, deadenylation-independent decay (GO:0031086)",
                              "Collagen biosynthesis and modifying enzymes Homo sapiens R-HSA-1650814",
                              "Scavenging by Class A Receptors Homo sapiens R-HSA-3000480",
                              "Cholesterol biosynthesis Homo sapiens R-HSA-191273",
                              "cellular response to corticosteroid stimulus (GO:0071384)")


pathway_to_plot[["MYD"]] <- c("Collagen biosynthesis and modifying enzymes Homo sapiens R-HSA-1650814",
                              "Antigen Presentation: Folding, assembly and peptide loading of class I MHC Homo sapiens R-HSA-983170",
                              "Cholesterol biosynthesis Homo sapiens R-HSA-191273",
                              "Pre-NOTCH Transcription and Translation Homo sapiens R-HSA-1912408",
                              "Scavenging by Class A Receptors Homo sapiens R-HSA-3000480")

pathway_to_plot[["STRO"]] <- c("Extracellular matrix organization Homo sapiens R-HSA-1474244",
                               "Degradation of the extracellular matrix Homo sapiens R-HSA-1474228",
                              "Collagen biosynthesis and modifying enzymes Homo sapiens R-HSA-1650814",
                              "Response to elevated platelet cytosolic Ca2+ Homo sapiens R-HSA-76005",
                              "cellular response to transforming growth factor beta stimulus (GO:0071560)")

pathway_to_plot[["MCR"]] <- c("Immune System Homo sapiens R-HSA-168256",
                              "MHC class II antigen presentation Homo sapiens R-HSA-2132295",
                              "Interferon gamma signaling Homo sapiens R-HSA-877300",
                              "neutrophil activation involved in immune response (GO:0002283)")

pathway_to_plot[["SRT"]] <- c("substantia nigra development (GO:0021762)",
                              "glyceraldehyde-3-phosphate metabolic process (GO:0019682)",
                              "Metabolism Homo sapiens R-HSA-1430728",
                              "energy coupled proton transport, down electrochemical gradient (GO:0015985)")


DGE_Heatmap2 = function(Seurat.object, 
                        cell_type, 
                        enrichR.path = './',
                        CND = 'up', 
                        database = c('Reactome_2016','GO_Biological_Process_2018'), 
                        pathways.2.plot,
                        gene.label = 6){

  DefaultAssay(Seurat.object) <- 'RNA'
  enrichR.file = paste(enrichR.path,
                       'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
  
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
  
  gene.list = character()
  for (pathway in pathways.2.plot) {
    gene.list.tmp <- unlist(strsplit(enrichR.table[pathway,]$Genes, ';'))
    gene.list = c(gene.list, gene.list.tmp)
  }
  #gene.list <- unlist(strsplit(enrichR.table[pathway,]$Genes, ';'))
  cell2select <- colnames(Seurat.object)[Seurat.object@meta.data$cell_type==cell_type]
  cell.type.cells <- subset(Seurat.object, cells = cell2select)
  cell.type.cells.scaled <- ScaleData(cell.type.cells)
  #options(repr.plot.width=15, repr.plot.height=3)
  HM <- DoHeatmap(cell.type.cells.scaled,
                  group.colors = c('orange', 'dodgerblue2'),
                  features = sort(gene.list), 
                  cells= cell2select,
                  size = 3, 
                  disp.min = -2, disp.max = 1.5,
                  angle = 45,
                  draw.lines =T,
                  label = F,
                  raster = TRUE) +
    scale_fill_gradientn(colours = coolwarm(200)) +
    theme(axis.text.y = element_text(color = 'black', size=gene.label))
  return(HM)
}

for (i in cell.type_2_plot) {
  print(i)
  ds = Seurat.object
  gl = 8; h = 6
  if (i == 'MCR') {gl = 2; h = 6.5}
  HP.col <- DGE_Heatmap2(Seurat.object = ds, 
                        cell_type = i, 
                        enrichR.path = indir,
                        CND = 'up', 
                        database = c('Reactome_2016','GO_Biological_Process_2018'), 
                        pathway_to_plot[[i]],
                        gene.label = gl)
  pdf(paste(outdir,'HM_DS_',i,'_2.pdf',sep=''),  width=9, height=h)
  print(HP.col)
  dev.off()
  Sdata = HP.col$data 
  Sdata = Sdata[!is.na(Sdata$Expression),] 
  Sdata[str_sub(Sdata$Cell, 1, 5) == "Donor" & str_sub(Sdata$Identity, -5, -1) == "_iGCA","Identity"] = 
    paste(i, "CTL", sep="_")
  write.xlsx(Sdata, paste(df_dir, 'SData_HM_DS_',i,'.xlsx',sep=''), 
             row.names = F)
  
  }


######## Dotplots ##########

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
    if (cell_type == 'LEY') {
      pathways.dataframe$Pathway[3] = "nuclear-transcribed mRNA catabolic process, deadenylation-independent \n decay (GO:0031086)"
    }
    if (cell_type == 'MYD') {
      pathways.dataframe$Pathway[4] = "Antigen Presentation: Folding, assembly and peptide loading of \n class I MHC R-HSA-983170"
    }
    if (cell_type == 'STRO') {
      pathways.dataframe$Pathway[5] = "cellular response to transforming growth factor beta stimulus \n (GO:0071560)"
    }
    if (cell_type == 'SRT') {
      pathways.dataframe$Pathway[1] = "energy coupled proton transport, down electrochemical gradient \n (GO:0015985)"
    }
    
  }
  
  PP = ggplot(pathways.dataframe, aes(Pathway.num,-log10(p.value))) + 
    geom_point(aes(size = gene.ratio), color = colors_Alf[cell_type,]$cols) +
    scale_size_continuous(range = c(5,15), name = "Gene ratio") +
    ylim(c(0,16)) +
    coord_flip() +
    scale_x_discrete(breaks=pathways.dataframe$Pathway.num, 
                     labels=pathways.dataframe$Pathway,
                     position = "top") +
    theme(plot.title = element_text(color="black", size=18, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = 'black', size=20,family = "mono"), 
          axis.title.x = element_text(face = "bold", color = "black", size = 14),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=25, family = "mono"),
          axis.title.y = element_text(face = "bold", color = "black", size = 14),
          legend.text = element_text(color = "black", size = 12),
          legend.title = element_text(face = "bold", color = "black", size = 14),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(title = 'Upregulated Pathways in Reactome and GO Biological Process', 
         x = "Pathways", 
         y = "-Log10(pvalue)") 
  
  assign(paste('PP',cell_type,sep=''),PP)
  
  pdf(paste(outdir,'PP_',cell_type,'.pdf',sep=''),  width=18, height=6)
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
#Seurat.object = Final_GSE124263_GSE112013_iNOA.integrated
Idents(Seurat.object) <- "celltype"
Seurat.object$cell_type = Seurat.object$celltype
Seurat.object$celltype.cond <- paste(Idents(Seurat.object), Seurat.object$condition, sep = "_")
Idents(Seurat.object) <- "celltype.cond"
DefaultAssay(Seurat.object) <- "RNA"

CND = 'down'
cell_types= c("LEY",  "MYD",  "MCR", "SRT","END", "STRO")


database = c('Reactome_2016','GO_Biological_Process_2018')
gene.l = list()
path.l = list()
gene.list = character()
indir = paste(dir, "enrichR/", sep='')
#outdir = paste(dir, "Pathways/", sep='')
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

prova2 = intersect(gene.l$LEY,gene.l$MYD)
length(prova2)


DefaultAssay(Seurat.object) <- "RNA"

#markers.to.plot <- AgeGenesClustersTable_dw$gene.name[AgeGenesClustersTable_dw$cluster==3]
dp_down = DotPlot(Seurat.object, 
                  col.min = -1,
                  col.max = 1,
                  features = prova2, 
                  cols = c("lightblue", "red"), 
                  dot.scale = 2) +
  #scale_size_continuous(range = c(0.2,5)) +
  RotatedAxis() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, face = "bold", color = "black", size=15,vjust =1.1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(angle = 0, color = "black", size=9, family = 'mono'), 
        axis.title.y = element_text(face = "bold", color = "black", size = 20),
        legend.text = element_text(face = "bold", color = "black", size = 16),
        legend.position="right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "Genes", y = "Cell type")

dp_down

write.xlsx(dp_down$data[,c("features.plot","id", "avg.exp", "avg.exp.scaled", "pct.exp")], 
           paste(df_dir,"DownregGene_dotplot.xlsx",sep=''),row.names = F)


pdf(paste(outdir,'DW_genes_exp.pdf',sep=''),width=8, height=18)
plot(dp_down)
dev.off()


p = rev(as.character(sort(Reduce(intersect, path.l))))


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
        axis.text.x = element_text(angle = 0, face = "bold", color = 'black', size=13, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16, family = 'mono'),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 18),
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))

pdf(paste(outdir,'DW_pathPLOT.pdf',sep=''),width=20, height=16)
plot(patplot)
dev.off()

outdir
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

pdf(paste(outdir,'Figure3_main.pdf',sep=''),width=25, height=18)
p1
dev.off()

im = image_read_pdf(paste(outdir,"Figure3_main.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure3_main.pdf.tiff", sep=''), format = "tiff")

col = as.character(colors_Alf[names(path.l),]$cols)

venn::venn(path.l, simplify=TRUE, opacity = 0.3, box = FALSE, ilab=TRUE, zcolor = col, ilcs = 1.5, sncs = 2)

pdf(paste(outdir,'Figure3_venn.pdf',sep=''),width=7, height=6.5)
print(venn::venn(path.l, simplify=TRUE, opacity = 0.3, box = FALSE, elipse = T, ilab=TRUE, zcolor = col, ilcs = 1.5, sncs = 2.5))
dev.off()
im = image_read_pdf(paste(outdir,"Figure3_venn.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure3_venn.pdf.tiff", sep=''), format = "tiff")

# ------ Table down pathway intersection --------
a = venn::venn(path.l, intersections = T)
list_intersection = attr(x = a, "intersections")
names(list_intersection) <- str_replace_all(names(list_intersection), ':', '_')
write.xlsx(list_intersection, paste(outdir,'intersection_Pathway.xlsx', sep=''))

df = data.frame()
for (element in names(list_intersection)) {
  entry = data.frame(Intersection = element, Pathways = list_intersection[[element]])
  df = rbind(df, entry)
}
write.xlsx(df, paste(outdir,'intersection_Pathway_2.xlsx'), asTable = T)


#### Tcell ####
# Blend plot
DimPlot(Seurat.object)
DefaultAssay(Seurat.object) <- 'RNA'
Idents(Seurat.object) <- "celltype"
fp_blend <- FeaturePlot(subset(Seurat.object, idents = 'TCL'), 
                        reduction = "umap", 
                        pt.size = 3,
                        features = c('CD8A','CD69'),
                        label = F, 
                        repel = T,
                        order= T,
                        blend.threshold = 0.5,
                        #slot = 'counts',
                        blend = T,       
                        max.cutoff = 'q99',
                        cols = c("lightgrey", "red","green")) + theme(legend.position = 'right')

pdf(paste(outdir, 'Figure2_Tcell_CD8_CD69_blend_tcellonly.pdf',sep =''), width = 12, height = 4)
fp_blend
dev.off()

fp_blend <- FeaturePlot(subset(Seurat.object, idents = 'TCL'), 
                        reduction = "umap", 
                        pt.size = 3,
                        features = c('CCL5','CCL4'),
                        label = F, 
                        repel = T,
                        order= T,
                        blend.threshold = 0.5,
                        #slot = 'counts',
                        blend = T,       
                        max.cutoff = 'q99',
                        cols = c("lightgrey", "red","green")) + theme(legend.position = 'right')
fp_blend

fp_blend <- FeaturePlot(subset(Seurat.object, idents = 'TCL'), 
                        reduction = "umap", 
                        pt.size = 3,
                        features = c('GZMM','GZMK'),
                        label = F, 
                        repel = T,
                        order= T,
                        blend.threshold = 0.5,
                        #slot = 'counts',
                        blend = T,       
                        max.cutoff = 'q99',
                        cols = c("lightgrey", "red","green")) + theme(legend.position = 'right')

pdf(paste(outdir, 'Figure2_Tcell_GZMM_GZMK_blend_tcellonly.pdf',sep =''),width = 12, height = 4)
fp_blend
dev.off()

vccl4 = VlnPlot(subset(Seurat.object, idents = 'TCL'), c("CCL4"), pt.size = 0) + 
  xlab("") + scale_fill_manual(values = c('orange')) + 
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) 
ggsave(paste(outdir, "Figure2_TCL_CCL4.png", sep=''), vccl4, height = 5, width = 4)

vccl5 = VlnPlot(subset(Seurat.object, idents = 'TCL'), c("CCL5"), pt.size = 0) + 
  xlab("") + scale_fill_manual(values = c('orange')) + 
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 18),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
        axis.title.y = element_text(face = "bold", color = "black", size = 18),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) 
ggsave(paste(outdir, "Figure2_TCL_CCL5.png", sep=''), vccl5, height = 5, width = 4)


path_file =outdir
fl = list.files(path_file, pattern=glob2rx("*.pdf"))
fl = paste("Figure1/","Figure1G.umap.pdf", sep='')
for (file in fl){
  im = image_read_pdf(paste(path_file,file,sep=''),density = 200)
  image_write(im, path = paste(path_file,file,'.tiff', sep=''), format = "tiff")
}

Idents(Seurat.object) <- "celltype"
TCo = subset(Seurat.object, idents = 'TCL')
head(TCo@assays$RNA@counts)[,1:5]

myTCL = TCo
FeaturePlot(myTCL, features = c("TUBA1A", "TUBA1B", "TUBB", "TUBB4B"))
Plot_sign(myTCL, signature = c("TUBA1A", "TUBA1B", "TUBB", "TUBB4B"), operator = sum)

save(TCo, file = paste(SO_dir, "TCL_OSR", sep=''))

matrixTcells = GetAssayData(TCo)
dat="cleanAdultIntegration"
write.table(matrixTcells, file = paste(dat,'_ScaleData_matrix_TCL_allintegrated.tsv',sep =''), quote = F)
TCL.pred = read.table(file = "SingleR_testis_FinalIntegration_TCL.tsv", quote = "", sep = '\t', header = T)
head(TCL.pred)
TCo$SingleR_BlueprintEncodeData = TCL.pred$labels
DefaultAssay(TCo) <- "RNA"
TCo_2 = RunUMAP(TCo, dims = 1:7)
DimPlot(TCo_2, group.by = "SingleR_BlueprintEncodeData", pt.size = 4, cols = viridis(4)) + 
  ggtitle("SingleR BlueprintEncodeData")
p = DimPlot(TCo, group.by = "SingleR_BlueprintEncodeData", pt.size = 4, cols = viridis(4)) + 
  ggtitle("SingleR BlueprintEncodeData") + 
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=12, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
        axis.title.y = element_text(face = "bold", color = "black", size = 14),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) 

p
A = as.data.frame(table(TCo$SingleR_BlueprintEncodeData, 
      TCo$celltype)[,"TCL"])
A2 = as.data.frame(prop.table(table(TCo$SingleR_BlueprintEncodeData, 
                    TCo$celltype))[,"TCL"])
A$N = A2$`prop.table(table(TCo$SingleR_BlueprintEncodeData, TCo$celltype))[, "TCL"]`
A$prediction = row.names(A)
A$celltype = "TCL"
col = viridis(4)
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
layout <- "
AAAABB
AAAABB
"
wrap_plots(A = p, B = N,  design = layout)
ggsave(paste(outdir, "Tcell_SingleR_mapping.png", sep=''), wrap_plots(A = p, B = N,  design = layout), 
       width = 10, height = 4)

#---------- Macrophages ---------------
GM = c("HLA-DRA", "CD36", "CD209", "CCL2", "CCL8")

c1 = "iGCA"; c2 = "CTL"
features = GM
Idents(Seurat.object) <- "celltype.cond"
dge_value_for_figure = data.frame()
for (ct in "MCR") {
  print(ct)
  value <- FindMarkers(Seurat.object,
                       ident.1 = paste(ct,"_",c1, sep=''),
                       ident.2 = paste(ct,"_",c2, sep=''), 
                       verbose = FALSE,
                       feature = features,
                       logfc.threshold=0,
                       min.pct = 0)
  value$celltype = ct
  value$gene = row.names(value)
  dge_value_for_figure = rbind(dge_value_for_figure, value)
}
dge_value_for_figure = dge_value_for_figure[order(dge_value_for_figure$gene),
                                            c("gene", "celltype", "avg_logFC", "p_val", "p_val_adj", "pct.1","pct.2")]

dge_value_for_figure
write.xlsx(dge_value_for_figure, file = paste(outdir, "Figure2_MCR_Genes_statistics.xlsx", sep=''), 
           row.names = FALSE, asTable = T)

GM = c( "CD86", "TNF", "ARG1")
GM = c("SOCS3" , "CCL5","NOS2", "TNF")
GM = c("CCL8", "CD84", "CD36",
       "CD59",
       "CD209",
       "CD81",
       "CD46",
       "CD74",
       "CD93",
       "CD22")
GM = c("CD163")
GM = c("TNF", "IL1B", "CCL4", "STAB1", "F13A1")
DefaultAssay(Seurat.object) <- "RNA"
Idents(Seurat.object) <- "celltype"
for (g in GM) {
  vp <- VlnPlot(subset(Seurat.object, idents = "MCR"),
                feature = g,
                split.plot = TRUE,
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
          legend.position = "top",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))
  print(vp)
  ggsave(paste(outdir, "Figure2_MCR_",g,".png", sep=''), vp, height = 5, width = 5)
}
Plot_sign(Seurat.object, signature = GM,operator = mean, title = "Macrophage")

x <- Seurat.object
DefaultAssay(x) <- "RNA"
signature = GM; operator = sum
x[["Sign_exp"]] <- apply(FetchData(object = x, 
                                   vars = signature),
                         1,
                         operator)
FP <- FeaturePlot(x, reduction = "umap", 
                  features = 'Sign_exp', 
                  label = T, 
                  order=T,
                  split.by = "condition",
                  cols = c("lightgrey", "red")) 
FP

ggsave(paste(outdir, "Figure2_TCL_CCL4.png", sep=''), vccl4, height = 5, width = 4)


vp <- VlnPlot(subset(Seurat.object, idents = "MCR"),
              feature = "CD74",
              split.plot = TRUE,
              pt.size = 0,
              split.by = 'condition')
fig2K = vp$data
fig2K = fig2K[order(fig2K$ident),]
fig2K$cellID = rownames(fig2K)
head(fig2K)
write.xlsx(fig2K[,colnames(fig2K)[c(4,2,3,1)]],paste(df_dir,"Figure2K.xlsx",sep=''),row.names = F)

DefaultAssay(Seurat.object) <- "RNA"
Idents(Seurat.object) <- "celltype"
Seurat.object_s <- subset(Seurat.object,  idents = c("LEY", "MYD","STRO"))
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


fig4Av1 = vp_IG1$data
fig4Av1 = fig4Av1[order(fig4Av1$ident),]
fig4Av1$cellID = rownames(fig4Av1)
fig4Av2 = vp_IG2$data
fig4Av2 = fig4Av2[order(fig4Av2$ident),]
fig4Av2$cellID = rownames(fig4Av2)
shcol = intersect(colnames(fig4Av1), colnames(fig4Av2))
fig4A = merge(fig4Av1, fig4Av2, by = shcol)
write.xlsx(fig4A,paste(df_dir,"Figure5_Aviolin.xlsx",sep=''),row.names = F)



# table(Seurat.object$condition)
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

fp_bl
fp_bl_df = list()
for (i in 1:8) {
  fp_bl_df[[i]] = fp_bl[[i]]$data
  fp_bl_df[[i]]$cellID = row.names(fp_bl_df[[i]])
}
iGCA = Reduce(merge, fp_bl_df[1:3])
CTL = Reduce(merge, fp_bl_df[5:7])
legend = fp_bl_df[[4]]
write.xlsx(iGCA[,colnames(iGCA)[c(4,1:3,5:7)]],paste(df_dir,"Figure5_Ablend_igca.xlsx",sep=''),row.names = F)
write.xlsx(CTL[,colnames(CTL)[c(4,1:3,5:7)]],paste(df_dir,"Figure5_Ablend_ctl.xlsx",sep=''),row.names = F)
write.xlsx(legend,paste(df_dir,"Figure5_Ablend_leg.xlsx",sep=''),row.names = F)
colnames(iGCA)


layout <- '
AACCCC
BBCCCC
' 
pdf(paste(outdir,'Figure5_COL1s_blend.pdf',sep=''),width=19, height=7)
wrap_plots(A = vp_IG1, B = vp_IG2, C = fp_bl, design = layout)
dev.off()

im = image_read_pdf(paste(outdir,"Figure5_COL1s_blend.pdf",sep=''),density = 140)
image_write(im, path = paste(outdir,"Figure5_COL1s_blend.tiff", sep=''), format = "tiff")

# ------------ figure5 statistics -------------- #

row.names(cleanAdultIntegration)[grep("COL",row.names(cleanAdultIntegration))]
c1 = "iGCA"; c2 = "CTL"
features = c("COL1A1", "COL1A2", "COL4A1", "COL4A2","COL4A3","COL4A4", "OGN")
Idents(Seurat.object) <- "celltype.cond"
dge_value_for_figure = data.frame()
for (ct in levels(Seurat.object$celltype)[levels(Seurat.object$celltype)!="TCL"]) {
  print(ct)
  value <- FindMarkers(Seurat.object,
                       ident.1 = paste(ct,"_",c1, sep=''),
                       ident.2 = paste(ct,"_",c2, sep=''), 
                       verbose = FALSE,
                       feature = features,
                       logfc.threshold=0,
                       min.pct = 0)
  value$celltype = ct
  value$gene = row.names(value)
  dge_value_for_figure = rbind(dge_value_for_figure, value)
}
dge_value_for_figure = dge_value_for_figure[order(dge_value_for_figure$gene),
                                            c("gene", "celltype", "avg_logFC", "p_val", "p_val_adj", "pct.1","pct.2")]

dge_value_for_figure
write.xlsx(dge_value_for_figure, file = paste(outdir, "Figure5_Genes_statistics.xlsx", sep=''), 
           row.names = FALSE, asTable = T)


# ------------- Supplementary ----------------
# --------------Jacard Heatmap ---------------
library('philentropy')
library('IntClust')
library('pheatmap')

dir = '/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/cleanAdultIntegration/'

indir = paste(dir, 'enrichR/', sep ='')
outdirJ = paste(outdir, "SI/", sep ='')
dir.create(outdirJ, recursive = TRUE, showWarnings = FALSE)

# --------------Jacard Heatmap UP ---------------

CND = 'up'
database = c('Reactome_2016','GO_Biological_Process_2018')
p.value.thr = 0.05
d = 'Jaccard_up/'
dir.create(paste(outdirJ,d,sep=''))

cell_types <- c("LEY",
                "MYD",
                "STRO",
                "MCR", "SRT")
#"STRO")

breaksList = seq(0, 1, by = 0.001)
cutree_rows_values = c(6,6,7,10)

CND = 'up'
database = c('Reactome_2016','GO_Biological_Process_2018')
p.value.thr = 0.05
for (N in 1:length(cell_types)) { 
  cell_type = cell_types[N]
  cutree_rows_N = cutree_rows_values[N]
  enrichR.file = paste(indir,'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
  
  Pathways.Table = data.frame()
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
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  gene.list = list()
  gene.all = character()
  pathways = unlist(row.names(Pathways.Table[Pathways.Table$Adjusted.P.value < p.value.thr,]))
  
  for (p in pathways) {
    gene.list[[p]] <- unlist(strsplit(Pathways.Table[p,]$Genes, ';')) 
    gene.all = c(gene.all, unlist(strsplit(Pathways.Table[p,]$Genes, ';')))
  }
  gene.all = unique(gene.all)
  
  # MAtrix
  M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
  row.names(M) = pathways 
  colnames(M) = gene.all 
  
  for (pat in pathways) {
    for (gene in gene.all) {
      if (gene %in% gene.list[[pat]]) {
        M[pat,gene] <- 1 
      }            
    }    
  }
  
  # Jaccard dist
  Jacard.Matrix <- distance(M, method = "jaccard")
  if (length(pathways)==2) {
    Jacard.Matrix_new = as.matrix(rbind(c(0,Jacard.Matrix),c(Jacard.Matrix,0)))
    Jacard.Matrix = Jacard.Matrix_new
  }
  row.names(Jacard.Matrix) <- pathways
  colnames(Jacard.Matrix) <- pathways
  
  pheatmap::pheatmap(Jacard.Matrix,
                     border_color = 'darkgrey',
                     color = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(breaksList)),   
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     cellwidth = 6, cellheight = 6,
                     cutree_rows = cutree_rows_N,
                     show_colnames = FALSE,
                     #main = paste(cell_type,'cells pathways - Jaccard distance heatmap'),
                     fontsize = 12,
                     fontsize_row = 6,
                     filename = paste(outdirJ,d,'JaccardDist_RP_GOBP_',cell_type,'_',CND,".pdf", sep='')) 
}

path_file = paste(outdirJ,d,sep='')
fl = list.files(path_file, pattern=glob2rx("*.pdf"))
fl
#fl = paste("DevelopmentGene/statistcs/","Figure1G.umap.pdf", sep='')
for (file in fl){
  im = image_read_pdf(paste(path_file,file,sep=''),density = 200)
  image_write(im, path = paste(path_file,file,'.tiff', sep=''), format = "tiff")
}

# --------------Jacard Heatmap DW ---------------
d = 'JACCARD_dw/'
dir.create(paste(outdirJ,d,sep=''))

cell_types <- c("LEY",
                "MYD",
                "SRT",
                "MCR",
                "STRO", "END")

breaksList = seq(0, 1, by = 0.001)

cutree_rows_values = c(7,8,5,5,4,5)

CND = 'down'
database = c('Reactome_2016','GO_Biological_Process_2018')
p.value.thr = 0.05
for (N in 1:length(cell_types)) { 
  cell_type = cell_types[N]
  cutree_rows_N = cutree_rows_values[N]
  enrichR.file = paste(indir,'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
  
  Pathways.Table = data.frame()
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
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  gene.list = list()
  gene.all = character()
  pathways = unlist(row.names(Pathways.Table[Pathways.Table$Adjusted.P.value < p.value.thr,]))
  
  for (p in pathways) {
    gene.list[[p]] <- unlist(strsplit(Pathways.Table[p,]$Genes, ';')) 
    gene.all = c(gene.all, unlist(strsplit(Pathways.Table[p,]$Genes, ';')))
  }
  gene.all = unique(gene.all)
  
  # MAtrix
  M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
  row.names(M) = pathways 
  colnames(M) = gene.all 
  
  for (pat in pathways) {
    for (gene in gene.all) {
      if (gene %in% gene.list[[pat]]) {
        M[pat,gene] <- 1 
      }            
    }    
  }
  
  # Jaccard dist
  Jacard.Matrix <- distance(M, method = "jaccard")
  if (length(pathways)==2) {
    Jacard.Matrix_new = as.matrix(rbind(c(0,Jacard.Matrix),c(Jacard.Matrix,0)))
    Jacard.Matrix = Jacard.Matrix_new
  }
  row.names(Jacard.Matrix) <- pathways
  colnames(Jacard.Matrix) <- pathways
  
  pheatmap::pheatmap(Jacard.Matrix,
                     border_color = 'darkgrey',
                     color = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(breaksList)),   
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     cellwidth = 1, cellheight = 1,
                     cutree_rows = cutree_rows_N,
                     show_colnames = FALSE,
                     #main = paste(cell_type,'cells pathways - Jaccard distance heatmap'),
                     fontsize = 12,
                     fontsize_row = 1,
                     filename = paste(outdirJ,d,'JaccardDist_RP_GOBP_',cell_type,'_',CND,".pdf", sep=''))
  #dev.off() 
}


path_file = paste(outdirJ,d,sep='')
fl = list.files(path_file, pattern=glob2rx("*.pdf"))
fl
#fl = paste("DevelopmentGene/statistcs/","Figure1G.umap.pdf", sep='')
for (file in fl){
  im = image_read_pdf(paste(path_file,file,sep=''),density = 200)
  image_write(im, path = paste(path_file,file,'.tiff', sep=''), format = "tiff")
}

# ------------- Imprinted --------------------
library(biomaRt)
ensembl97 = useMart("ensembl",
                    dataset="hsapiens_gene_ensembl", 
                    host = 'http://jul2019.archive.ensembl.org')

IG_directory='/Users/tascini.annasofia/Dropbox (HSR Global)/Alfano_904_paperDraft/PaternalGenes/'

Pgene.table = read.xlsx(paste(IG_directory,'Paternal_Genes.xlsx',sep=''),
                        sheet='Paternal')
Mgene.table = read.xlsx(paste(IG_directory,'Paternal_Genes.xlsx',sep=''),
                        sheet='Maternal')
row.names(Mgene.table) <- Mgene.table$Gene

DefaultAssay(Seurat.object) <- "RNA"
features = row.names(Seurat.object)
PG_expressed = intersect(Pgene.table$Gene, features)
fraction = (length(PG_expressed)/length(Pgene.table$Gene))*100
print(paste('% parental genes expressed in azoospermia', as.character(fraction),'%',sep = ' '))

MG_expressed = intersect(Mgene.table$Gene, features)
MG_expressed = MG_expressed[MG_expressed!='KCNK9']
fraction = (length(MG_expressed)/length(Mgene.table$Gene))*100
print(paste('% maternal genes expressed in azoospermia', as.character(fraction),'%',sep = ' '))


# -------- Imprinted DotPlots ---------
markers.to.plot <- PG_expressed
dp = DotPlot(Seurat.object, 
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

pdf(paste(outdir,'SI/','paternal_Dotplot.pdf',sep=''),width=18, height=6.5)
plot(dp)
dev.off()

markers.to.plot <- MG_expressed
dp_m = DotPlot(Seurat.object, 
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

dp_m

pdf(paste(outdir,'SI/','maternal_Dotplot.pdf',sep=''),width=17, height=6.5)
plot(dp_m)
dev.off()

# ------- Imprinted violin  ------
relevant_imprinted = c('H19', 'MEG3', 'MEG8','IGF2', 'DLK1', 'HYMAI', 'PEG3', 'ERAP2', 'RTL1', "KCNQ1OT1")
dir.create(paste('ImprintedGenes/', sep = ''))
dir.create(paste('ImprintedGenes/VlnPlot_r', sep = ''))

DefaultAssay(Seurat.object) <- "RNA"
Idents(Seurat.object) <- "celltype"
so = subset(Seurat.object, idents = c('LEY', 'MYD'))
Sdata = list()
for (gene in relevant_imprinted) {
  so[[gene]] <- FetchData(object = so, vars = gene)
  df= so[[]]
  #df$condition = revalue(df$condition, c("iNOA"="iGCA", "CTL"="CTL"))
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
  
  Sdata[[gene]] = v$data[,c("celltype", "condition", gene)]
  Sdata[[gene]]$cellID = row.names(Sdata[[gene]])
  
  pdf(paste('../ImprintedGenes/VlnPlot_r/','vlnPlot_',gene,'.pdf', sep=''), width = 5, height = 4.5)
  print(v)
  dev.off()
  
}

fig4G = Reduce(merge, Sdata)
write.xlsx(fig4G, paste(df_dir,"Figure4_G.xlsx",sep=''),row.names = F)

pdf(paste(dir, '../ImprintedGenes/VlnPlot_r/','vlnPlot_maternal.pdf', sep=''), width = 12, height = 4)
print(v_c_H19 + v_c_MEG3 + v_c_MEG8 + plot_layout(guides = "collect", ncol = 3) & 
        theme(legend.position = "right", legend.direction = 'vertical'))
dev.off()

pdf(paste(dir, '../ImprintedGenes/VlnPlot_r/','vlnPlot_paternal.pdf', sep=''), width = 28, height = 4)
print(v_c_IGF2 + v_c_DLK1 + v_c_HYMAI + v_c_PEG3 + v_c_ERAP2 + v_c_RTL1 + v_c_KCNQ1OT1 + plot_layout(guides = "collect", ncol = 7) & 
        theme(legend.position = "right", legend.direction = 'vertical'))
dev.off()


# --------- statistics of Imprinted by cluster --------------
cellTypes = c("LEY",
              "MYD",
              "SRT",
              'MCR',
              'END',
              "STRO",
              "UND")
signature = PG_expressed
azoospermia.modulated.Pimprinted= list()
biomartCacheClear()
ensembl97 = useMart("ensembl",
                    dataset="hsapiens_gene_ensembl", 
                    host = 'http://jul2019.archive.ensembl.org')

details_PG = biomaRt::getBM(attributes= c('external_gene_name','chromosome_name',
                                 'start_position', 'end_position','description'),
                            filters = "external_gene_name",
                            values = Pgene.table$Gene,
                            mart = ensembl97,
                            useCache = FALSE)

setdiff(Pgene.table$Gene, details_PG$external_gene_name)
details_PG = details_PG[order(details_PG$external_gene_name),]
details_MG = getBM(attributes= c('external_gene_name','chromosome_name',
                                 'start_position', 'end_position','description'),
                   filter = 'external_gene_name',
                   values = Mgene.table$Gene,
                   mart = ensembl97,
                   useCache = FALSE)
setdiff(Mgene.table$Gene, details_MG$external_gene_name)

Idents(Seurat.object) <- "celltype.cond"
for (cell_type in cellTypes) {
  print(cell_type)
  azoospermia.response.gene <- FindMarkers(Seurat.object,
                                           ident.1 = paste(cell_type,"_iGCA", sep=''),
                                           ident.2 = paste(cell_type,"_CTL", sep=''), 
                                           verbose = FALSE,
                                           feature = signature,
                                           logfc.threshold=0,
                                           min.pct = 0)
  azoospermia.response.gene$geneID <- row.names(azoospermia.response.gene)
  myorder=c('geneID','avg_logFC','p_val_adj','p_val','pct.1','pct.2')
  azoospermia.response.gene <- azoospermia.response.gene[,myorder]
  colnames(azoospermia.response.gene) <- c('geneID','avg_logFC','p_val_adj','p_val','pct.iNOA','pct.CTL')
  azoospermia.modulated.Pimprinted[[cell_type]] =  azoospermia.response.gene
  #azoospermia.modulated.Pimprinted[[cell_type]] = 
  #  azoospermia.modulated.Pimprinted[[cell_type]][!is.na(azoospermia.modulated.Pimprinted[[cell_type]]$avg_logFC),]
  azoospermia.modulated.Pimprinted[[cell_type]] = 
    azoospermia.modulated.Pimprinted[[cell_type]][order(azoospermia.modulated.Pimprinted[[cell_type]]$geneID),]
  descriptions = character()
  for (g in azoospermia.modulated.Pimprinted[[cell_type]]$geneID) {
    d = unique(details_PG[details_PG$external_gene_name == g,]$description)
    descriptions = c(descriptions, d)
  }
  azoospermia.modulated.Pimprinted[[cell_type]]$description = descriptions
  
}

write.xlsx(x = azoospermia.modulated.Pimprinted,
           file = paste(outdir, 'DGE_paternal.xlsx',sep=''),
           quote = F, row.names = F, colnames= T, asTable = T)


signature = MG_expressed
azoospermia.modulated.imprinted= list()
for (cell_type in cellTypes) {
  print(cell_type)
  azoospermia.response.gene <- FindMarkers(Seurat.object,
                                           ident.1 = paste(cell_type,"_iGCA", sep=''),
                                           ident.2 = paste(cell_type,"_CTL", sep=''), 
                                           verbose = FALSE,
                                           feature = signature,
                                           logfc.threshold=0,
                                           min.pct = 0)
  azoospermia.response.gene$geneID <- row.names(azoospermia.response.gene)
  myorder=c('geneID','avg_logFC','p_val_adj','p_val','pct.1','pct.2')
  azoospermia.response.gene <- azoospermia.response.gene[,myorder]
  colnames(azoospermia.response.gene) <- c('geneID','avg_logFC','p_val_adj','p_val','pct.iNOA','pct.CTL')
  azoospermia.modulated.imprinted[[cell_type]] =  azoospermia.response.gene
  azoospermia.modulated.imprinted[[cell_type]] = 
    azoospermia.modulated.imprinted[[cell_type]][order(azoospermia.modulated.imprinted[[cell_type]]$geneID),]
  descriptions = character()
  for (g in azoospermia.modulated.imprinted[[cell_type]]$geneID) {
    d = unique(details_MG[details_MG$external_gene_name == g,]$description)
    descriptions = c(descriptions, d)
  }
  azoospermia.modulated.imprinted[[cell_type]]$description = descriptions
  
}

write.xlsx(x = azoospermia.modulated.imprinted,
           file = paste(outdir, 'DGE_maternal.xlsx',sep=''),
           quote = F, row.names = F, colnames= T, asTable = T)



path_file = paste(outdir,"SI/",sep='')
fl = list.files(path_file, pattern=glob2rx("*.pdf"))
fl
#fl = paste("DevelopmentGene/statistcs/","Figure1G.umap.pdf", sep='')
for (file in fl){
  im = image_read_pdf(paste(path_file,file,sep=''),density = 200)
  image_write(im, path = paste(path_file,file,'.tiff', sep=''), format = "tiff")
}

######### SRT iNOA ####
dpaper = paste(outdir, 'SRT_iNOA/', sep = '')
dir.create(dpaper)
file = '/Users/tascini.annasofia/Downloads/41467_2020_19414_MOESM8_ESM.txt'
x = file
raw <- read.delim(x, check.names=FALSE, stringsAsFactors=FALSE, sep ='\t')
#raw <- raw[-1, ]
res <- as.list(raw)
res <- lapply(res, function(x) x[!is.na(x) & nchar(x) > 0])
names(res) <- colnames(raw)

r2 = list()
for (ls in names(res)) {
  r2[[ls]] = intersect(res[[ls]], row.names(Seurat.object))
  print(setdiff(res[[ls]],r2[[ls]]))
}
lengths(res)
lengths(r2)


for (ls in names(r2)) {
  print(ls)
  x <- Seurat.object
  DefaultAssay(x) <- "RNA"
  x[[ls]] <- apply(FetchData(object = x, 
                             vars = r2[[ls]]),
                   1,
                   median)
  VlnPlot(x,
          features = ls,
          split.by = 'condition', 
          pt.size = 0,
          split.plot = T) +
    ylab("Median of Expression") +
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
    #scale_fill_manual(values = viridis(8))
    scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))
  ggsave(paste(dpaper, "VlnPlot_",ls,".png", sep = ''), width = 12, height = 4)
}

gene = c('FOSB', 'EGR3', 'HOPX', 'TSC22D1', 'JUN', 'NR4A1', 'ENO1','BEX1','DEFB119')
for (g in gene) {
  VlnPlot(x,
         features = g,
          split.by = 'condition', 
          pt.size = 0.2,
          split.plot = T) +
    ylab("Sum of Expression") +
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
    scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c('iGCA','CTL'))
  ggsave(paste(dpaper, "VlnPlot_",g,".png", sep = ''), width = 12, height = 4)
}

table(cleanAdultIntegration$celltype)

# ----------- ribosomal -------------
dir.create(paste(outdir, "SI_Ribosomal/", sep =''))

load("/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/AlfanoM_904_infertilita_epigenetics/testis_GSE124263")
Idents(testis) = "dev"
testis_ad = subset(testis, idents = "adult")
SO = testis_ad
SO <- NormalizeData(SO, verbose = FALSE)
SO <- FindVariableFeatures(SO, 
                           selection.method = 'vst',
                           nfeatures = 2000,
                           verbose = FALSE)
SO <- ScaleData(SO, features = row.names(SO), verbose = FALSE)
SO <- RunPCA(SO, features = VariableFeatures(object = SO), verbose = FALSE)
SO <- FindNeighbors(SO, dims = 1:20, verbose = FALSE)
SO <- FindClusters(SO, resolution = 0.1, verbose = FALSE)
SO <- RunUMAP(SO, dims = 1:20, verbose = FALSE)
DimPlot(SO, label = T)

new.cluster.ids <- c( "LEY",
                      "S'gonia", 
                      "MYD",
                      "SSCs",
                      "END",
                      "S'tidis",
                      "S'cytes",
                      "MCR",
                      "STRO",
                      "Sperm",
                      "SRT")

names(new.cluster.ids) <- levels(SO)
SO2 <- RenameIdents(SO, new.cluster.ids)
DimPlot(SO2)
Idents(SO) <- "RNA_snn_res.0.1"
SO2[["perc.RP"]] <- PercentageFeatureSet(SO, pattern = "^RP")
RP_plot <- FeaturePlot(SO2,
                       repel = F,
                       reduction = "umap",
                       features = "perc.RP", 
                       pt.size = 1,
                       order = T,
                       label = T, 
                       max.cutoff = 'q98') + 
  ggtitle("") + scale_color_continuous(name = "%RP", 
                                      type = "gradient",
                                      low = "lightgrey",
                                      high = "red")
print(RP_plot)
pdf(file = paste(outdir, "SI_Ribosomal/RP_GSE124263",
                 ".pdf", sep =''),
    width = 7, height = 6)
print(RP_plot)
dev.off()
listSO = SplitObject(cleanAdultIntegration, split.by = "source")
for (i in names(listSO)) {
  SO = listSO[[i]]
  DefaultAssay(SO) <- "RNA"
  SO <- NormalizeData(SO, verbose = FALSE)
  SO <- FindVariableFeatures(SO, 
                                  selection.method = 'vst',
                                  nfeatures = 2000,
                                  verbose = FALSE)
  SO <- ScaleData(SO, features = row.names(SO), verbose = FALSE)
  SO <- RunPCA(SO, features = VariableFeatures(object = SO), verbose = FALSE)
  SO <- FindNeighbors(SO, dims = 1:20, verbose = FALSE)
  SO <- FindClusters(SO, resolution = 0.5, verbose = FALSE)
  SO <- RunUMAP(SO, dims = 1:10, verbose = FALSE)
  DimPlot(SO, reduction = "umap", group.by = "celltype")
  SO[["perc.RP"]] <- PercentageFeatureSet(SO, pattern = "^RP")
  Idents(SO) <- "celltype"
  RP_plot <- FeaturePlot(SO,
                        repel = T,
                        reduction = "umap",
                        features = "perc.RP", 
                        pt.size = 1,
                        order = T,
                        label = T, 
                        max.cutoff = 'q98') + 
    ggtitle(i) + scale_color_continuous(name = "%RP", 
                                        type = "gradient",
                                        low = "lightgrey",
                                        high = "red")
  print(RP_plot)
  pdf(file = paste(outdir, "SI_Ribosomal/RP_",str_remove_all(string = i, pattern = " "),
                   ".pdf", sep =''),
      width = 5, height = 5)
  print(RP_plot)
  dev.off()
}

# ----------- save celltypes -------------
listA = SplitObject(cleanAdultIntegration, split.by = "source")
celltype_excel = list()
for (c in names(listA)[1:3]) {
  metadata = listA[[c]]@meta.data
  celltypeDF = metadata[,c("source", "cell_type")]
  celltypeDF$Barcode = str_sub(row.names(celltypeDF), 7,-1)
  celltype_excel[[unique(celltypeDF$source)]] = celltypeDF[, c("Barcode","cell_type")]
}
write.xlsx(x = celltype_excel, file = paste(outdir, "Table_celltypes_iNOA.xlsx", sep=''), row.names= F)
