#cellphoneDB prepare and plot
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
library(scales)
library(viridis)
library(SeuratWrappers)
###library(slingshot)
#require(BiocStyle)
library(SingleCellExperiment)

setwd("/beegfs/scratch/ric.cosr/ric.alfano/AlfanoM_904_infertilita_epigenetics/CellPhoneDB/")

indir = "/beegfs/scratch/ric.cosr/ric.alfano/AlfanoM_904_infertilita_epigenetics/SeuratObjects/"
outdir = "SeuratOutputs/"
dir.create(outdir, recursive = T)

load(paste(indir, "cleanAdultIntegration", sep=''))
DefaultAssay(cleanAdultIntegration) <- "RNA"

listAlfano = SplitObject(cleanAdultIntegration, split.by = "condition")

# ---------- prepare iGCA for cellphoneDB analysis -------------- 
count_raw_iNOA =as.matrix(GetAssayData(object = listAlfano$iNOA, slot = "counts"))
head(rownames(count_raw_iNOA))
#normalize 
count_norm_iNOA <- apply(count_raw_iNOA, 2, function(x) (x/sum(x))*10000)
Gene <- rownames(count_norm_iNOA)
count_norm_iNOA <- cbind(Gene,count_norm_iNOA)
head(count_norm_iNOA[,1:5])

meta_data_iNOA <- cbind(rownames(listAlfano$iNOA@meta.data),as.data.frame(Idents(listAlfano$iNOA)))
head(meta_data_iNOA)
colnames(meta_data_iNOA)<- c("Cell","cell_type")
keeps <- colnames(count_norm_iNOA)
mySubset_iNOA <- meta_data_iNOA[meta_data_iNOA$Cell %in% keeps, ]

dim(count_norm_iNOA)
dim(mySubset_iNOA)

# writeOutputs to use with cellphoneDB
write.table(count_norm_iNOA, paste(outdir, 'cellphonedb_count_iGCA.txt', sep=''), 
            sep='\t', quote=F, row.names = F)
write.table(mySubset_iNOA, paste(outdir,'cellphonedb_meta_iGCA.txt', sep=''), 
            sep='\t', quote=F, row.names=F)
setdiff(colnames(count_norm_iNOA), mySubset_iNOA$Cell)

# ------- iGCA cellphoneDB ---------
# cellphonedb method statistical_analysis ./SeuratOutputs/cellphonedb_meta_iGCA.txt SeuratOutputs/cellphonedb_count_iGCA.txt --counts-data=gene_name --threads=4 --iterations=100
# cellphonedb plot dot_plot
# cellphonedb plot heatmap_plot SeuratOutputs/cellphonedb_meta_iGCA.txt
integration_heatmap_list = list()
toplot_matrix = read.table("./out_iGCA/count_network.txt", header=T)
#toplot_matrix = read.table("./out/count_network.txt", header=T)
integration_heatmap_list[["iGCA"]] = toplot_matrix 
toplot_matrix1 <- toplot_matrix[order(toplot_matrix$SOURCE, toplot_matrix$TARGET),]
toplot_matrix1

library(pheatmap)
library(reshape2)
matrix_to_plot= acast(toplot_matrix1, SOURCE~TARGET, value.var="count")
head(matrix_to_plot)
oredr = c("LEY", "MYD", "STRO","END", "MCR", "TCL", "SRT","UND") 
matrix_to_plot_n = matrix_to_plot[oredr,oredr]
head(matrix_to_plot_n)
# create heatmap using pheatmap
library(rcartocolor)
my_colors = rev(carto_pal(7, "Earth"))
col1 = "dodgerblue4"; col2 = 'peachpuff'; col3 = 'deeppink4'
col.heatmap <- colorRampPalette(c(col1,col2,col3))( 1000 )
max(matrix_to_plot_n)
Breaks <- seq(min(0), max(120), length = 1000)
family='Arial'
pheatmap(matrix_to_plot_n, 
         cluster_rows = F,
         cluster_cols = F, color = col.heatmap,
         show_rownames = T,
         breaks = Breaks,
         family = family,
         fontsize = 18,
         filename = "iGCA_heatmap.pdf",
         main = 'iGCA',
         legend = TRUE)

dev.off()

# ---------- prepare CTL for cellphoneDB analysis -------------- 

#row counts
count_raw_CTL = as.matrix(GetAssayData(object = listAlfano$CTL, slot = "counts"))
#normalize 
count_norm_CTL <- apply(count_raw_CTL, 2, function(x) (x/sum(x))*10000)
Gene <- rownames(count_norm_CTL)
count_norm_CTL <- cbind(Gene,count_norm_CTL)
head(count_norm_CTL[,1:5])

meta_data_CTL <- cbind(rownames(listAlfano$CTL@meta.data),as.data.frame(Idents(listAlfano$CTL)))
head(meta_data_CTL)
colnames(meta_data_CTL)<- c("Cell","cell_type")
keeps <- colnames(count_norm_CTL)
mySubset_CTL <- meta_data_CTL[meta_data_CTL$Cell %in% keeps, ]

dim(count_norm_CTL)
dim(mySubset_CTL)

# writeOutputs to use with cellphoneDB
write.table(count_norm_CTL, paste(outdir, 'cellphonedb_count_CTL.txt', sep=''), 
            sep='\t', quote=F, row.names = F)
write.table(mySubset_CTL, paste(outdir,'cellphonedb_meta_CTL.txt', sep=''), 
            sep='\t', quote=F, row.names=F)
setdiff(colnames(count_norm_CTL), mySubset_CTL$Cell)

# ------- CTL cellphoneDB ---------
# cellphonedb method statistical_analysis ./SeuratOutputs/cellphonedb_meta_iGCA.txt SeuratOutputs/cellphonedb_count_iGCA.txt --counts-data=gene_name --threads=4 --iterations=100
# cellphonedb plot dot_plot
# cellphonedb plot heatmap_plot SeuratOutputs/cellphonedb_meta_iGCA.txt
toplot_matrix = read.table("./out_CTL/count_network.txt", header=T)
integration_heatmap_list[["CTL"]] = toplot_matrix 
toplot_matrix1 <- toplot_matrix[order(toplot_matrix$SOURCE, toplot_matrix$TARGET),]
toplot_matrix1
library(pheatmap)
library(reshape2)
matrix_to_plot= acast(toplot_matrix1, SOURCE~TARGET, value.var="count")
oredr = c("LEY", "MYD", "STRO","END", "MCR","SRT","UND") 
matrix_to_plot_n_CTL = matrix_to_plot[oredr,oredr]
matrix_to_plot_n_CTL
#matrix_to_plot_n_CTL$TCL = 0
A = rbind(matrix_to_plot_n_CTL, c(0,0,0,0,0,0,0))
B = cbind(A, c(0,0,0,0,0,0,0,0))
B
row.names(B)[8] = "TCL"
colnames(B)[8] = "TCL"
oredr = c("LEY", "MYD", "STRO","END",  "MCR", "TCL", "SRT","UND") 
B2 <- B[oredr,oredr]
B2
oredr
# create heatmap using pheatmap
col1 = "dodgerblue4"; col2 = 'peachpuff'; col3 = 'deeppink4'
col.heatmap <- colorRampPalette(c(col1,col2,col3))( 1000 )
max(matrix_to_plot_n_CTL)
Breaks <- seq(min(0), max(120), length = 1000)
family='Arial'
pheatmap(matrix_to_plot_n_CTL, 
         cluster_rows = F,
         cluster_cols = F, color = col.heatmap,
         show_rownames = T,
         breaks = Breaks,
         family = family,
         fontsize = 18,
         main = 'CTL', 
         filename = "CTL_heatmap.pdf",
         legend = TRUE)
B2
pheatmap(B2, 
         cluster_rows = F,
         cluster_cols = F, color = col.heatmap,
         show_rownames = T,
         breaks = Breaks,
         family = family,
         fontsize = 18,
         main = 'CTL', 
         filename = "CTL_heatmap_TCL.pdf",
         legend = TRUE)
dev.off()

write.xlsx(integration_heatmap_list, "SourceData_integration_heatmap.xlsx", row.names = FALSE)

# -------- Dot plots ---------

toplot_matrix = read.table("../Cellphone/out_iGCA/significant_means.txt", 
                           header=T, sep = "\t", check.names = FALSE,
                           as.is = TRUE)
head(toplot_matrix)
A = colnames(toplot_matrix[, grep("LEY", names(toplot_matrix))])

library(ggplot2)
dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = 'plot.pdf',
                    width = 8,
                    height = 10,
                    minscale = -3,
                    maxscale = 3,
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){
  
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }
  
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  
  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  #,
  pl = ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette, 
                          limits = c(minscale,maxscale)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
  write.xlsx(pl$data, 
             paste(df_dir,filename,"_Cellphone_dotplot.xlsx",sep=''),row.names = F)
  
  if (output_extension == '.pdf') {
    ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
    ggsave(filename, width = width, height = height, limitsize=F)
  }
}

sr_all = c()
sr= toplot_matrix[grep("CD74",toplot_matrix$interacting_pair),]$interacting_pair
sr_all = c(sr_all, sr)
sc = colnames(toplot_matrix[, grep("MCR", names(toplot_matrix))])
dot_plot(selected_rows = sr,
         selected_columns = sc,
         width = 8, height = 3,
         minscale = -2,
         maxscale = 7,
         filename = "CD74_dotplot_iGCA.pdf",
         means_path = "../Cellphone/out_iGCA/means.txt",
         pvalues_path = "../Cellphone/out_iGCA/pvalues.txt")

dot_plot(selected_rows = sr,
         selected_columns = sc[!grepl("TCL",sc)],
         width = 8, height = 3,
         minscale = -2,
         maxscale = 7,
         filename = "CD74_dotplot_CTL.pdf",
         means_path = "out_CTL/means.txt",
         pvalues_path = "out_CTL/pvalues.txt")

sr= toplot_matrix[grep("DLK1",toplot_matrix$interacting_pair),]$interacting_pair
sr_all = c(sr_all, sr)
sc = colnames(toplot_matrix[, grep("LEY", names(toplot_matrix))])
dot_plot(selected_rows = sr,
         selected_columns = sc,
         width = 10, height = 4,
         minscale = -4,
         maxscale = 4,
         filename = "DLK1_dotplot_iGCA.pdf",
         means_path = "out_iGCA/means.txt",
         pvalues_path = "out_iGCA/pvalues.txt")


dot_plot(selected_rows = sr,
         selected_columns = sc[!grepl("TCL",sc)],
         width = 10, height = 4,
         minscale = -4,
         maxscale = 4,
         filename = "DLK1_dotplot_CTL.pdf",
         means_path = "out_CTL/means.txt",
         pvalues_path = "out_CTL/pvalues.txt")
#C5AR1_RPS19


sr= toplot_matrix[grep("COL1A",toplot_matrix$interacting_pair),]$interacting_pair
sr_all = c(sr_all, sr)
sr
sc = colnames(toplot_matrix[, grep("LEY", names(toplot_matrix))])
dot_plot(selected_rows = sr,
         selected_columns = sc,
         width = 10, height = 5,
         minscale = -4,
         maxscale = 4,
         filename = "COL1_dotplot_iGCA.pdf",
         means_path = "out_iGCA/means.txt",
         pvalues_path = "out_iGCA/pvalues.txt")

sc = colnames(toplot_matrix[, grep("LEY", names(toplot_matrix))])
dot_plot(selected_rows = sr,
         selected_columns = sc[!grepl("TCL",sc)],
         width = 10, height = 5,
         minscale = -4,
         maxscale = 4,
         filename = "COL1_dotplot_CTL.pdf",
         means_path = "out_CTL/means.txt",
         pvalues_path = "out_CTL/pvalues.txt")


sr= toplot_matrix[grep("COL4A",toplot_matrix$interacting_pair),]$interacting_pair
sr
sc = colnames(toplot_matrix[, grep("LEY", names(toplot_matrix))])
dot_plot(selected_rows = sr,
         selected_columns = sc,
         width = 10, height = 5,
         minscale = -2,
         maxscale = 2,
         filename = "COL4_dotplot_iGCA.pdf",
         means_path = "out_iGCA/means.txt",
         pvalues_path = "out_iGCA/pvalues.txt")

sc = colnames(toplot_matrix[, grep("LEY", names(toplot_matrix))])
dot_plot(selected_rows = sr,
         selected_columns = sc[!grepl("TCL",sc)],
         width = 10, height = 5,
         minscale = -2,
         maxscale = 2,
         filename = "COL4_dotplot_CTL.pdf",
         means_path = "out_CTL/means.txt",
         pvalues_path = "out_CTL/pvalues.txt")

sr= toplot_matrix[grep("IGF2",toplot_matrix$interacting_pair),]$interacting_pair
sr_all = c(sr_all, sr)
sr
sc = colnames(toplot_matrix[, grep("LEY", names(toplot_matrix))])
dot_plot(selected_rows = sr,
         selected_columns = sc,
         width = 10, height = 4,
         minscale = -3,
         maxscale = 3,
         filename = "IGF2_dotplot_iGCA.pdf",
         means_path = "out_iGCA/means.txt",
         pvalues_path = "out_iGCA/pvalues.txt")
dot_plot(selected_rows = sr,
         selected_columns = sc[!grepl("TCL",sc)],
         width = 10, height = 4,
         minscale = -3,
         maxscale = 3,
         filename = "IGF2_dotplot_CTL.pdf",
         means_path = "out_CTL/means.txt",
         pvalues_path = "out_CTL/pvalues.txt")


sr_all = c(sr_all, "CD40_INSL3", "HLA-C_FAM3C", "HLA-DRB1_OGN", "HLA-DPA1_TNFSF9", "C5AR1_RPS19")
sr= toplot_matrix[grep("C5AR1",toplot_matrix$interacting_pair),]$interacting_pair
sr
#sr_all = c(sr_all, sr)
sr_all = unique(sr_all)
sr_all = sr_all[!sr_all=="C5AR1 RPS19"]
dot_plot(selected_rows = sr_all,
         selected_columns = NULL,
         width = 18, height = 10,
         minscale = -8,
         maxscale = 8,
         filename = "Selected_dotplot_iGCA.pdf",
         means_path = "../Cellphone/out_iGCA/means.txt",
         pvalues_path = "../Cellphone/out_iGCA/pvalues.txt")
dot_plot(selected_rows = sr_all,
         selected_columns = NULL,
         width = 18, height = 10,
         minscale = -8,
         maxscale = 8,
         filename = "Selected_dotplot_CTL.pdf",
         means_path = "../Cellphone/out_CTL/means.txt",
         pvalues_path = "../Cellphone/out_CTL/pvalues.txt")

#CD40 INSL3
#HLA-DRB1 OGN
#C5AR1 RPS19
# HLA-C_FAM3C
Seurat.object = cleanAdultIntegration
cleanAdultIntegration
DimPlot(cleanAdultIntegration)
Seurat.object$condition = plyr::revalue(Seurat.object$condition, c("iNOA"="iGCA", "CTL"="CTL"))
Seurat.object$condition = factor(Seurat.object$condition, levels = c("iGCA", "CTL"))
VlnPlot(Seurat.object, idents = c("LEY", "MYD", "STRO"), 
        features = c("COPA", "APP", "MIF"), split.by = "condition", 
        pt.size = 0, cols = c('orange','dodgerblue2'))


VlnPlot(Seurat.object, idents = "MCR", 
        features = c("CD74"), group.by = "condition", pt.size = 0, cols = c('orange','dodgerblue2'))

c1 = "iNOA"; c2 = "CTL"
features = c("COPA", "APP", "MIF", "CD74")
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
