#mySeuratfunctions.R

##### CUSTOM VIOLIN PLOT #######

# function to plot a violin plot of a list of genes in a cell type subset
VlnPlt_subset <- function(Seraut.object, cell_type, features){   
    SelectedCell = character()
    for (i in 1:length(cell_type)) {
        SelectedCell.tmp <- colnames(Seraut.object)[Seraut.object@meta.data$cell_type == cell_type[i]]
        SelectedCell <- c(SelectedCell, SelectedCell.tmp)
    }   
    
    
    Cell.subset <- subset(Seraut.object, cells = SelectedCell, invert=F)
    subset.counts <- Cell.subset@assays$RNA@counts
    gene.counts.dataframe = data.frame(gene = character(),
                                       cell_barcode = character(),
                                       counts = integer(),
                                       stringsAsFactors=FALSE)
    for (gene in features) {
        gene.counts.tmp = subset.counts[gene,]
        gene.counts.dataframe.entry = data.frame(gene = rep(gene, length(gene.counts.tmp)),
                                                 cell_barcode =  names(gene.counts.tmp), 
                                                 counts = (gene.counts.tmp))
        gene.counts.dataframe <- rbind(gene.counts.dataframe, gene.counts.dataframe.entry)
        }
    vp <- ggplot(gene.counts.dataframe, aes(x=gene, y = counts+1, fill=gene)) +
            geom_violin(trim=FALSE) + 
            scale_y_continuous(trans='log10') +
            geom_jitter(shape=16, position=position_jitter(0.4)) +  
      theme(plot.title = element_text(color="blue", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 60, face = "bold", color = 'dodgerblue4', size=22, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "dodgerblue2", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = 'dodgerblue4', size=22),
        axis.title.y = element_text(face = "bold", color = "dodgerblue2", size = 24),
        legend.text = element_text(face = "bold", color = "dodgerblue2", size = 22),
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
        labs(title = cat('Expression subset ', cell_type, sep = ''), 
               x = "Genes", 
               y = "Log Expression Level") 
    return(vp)
}

#vp <- VlnPlt_subset(Seraut.object = azoospermia.integrated.subset2,
#                   cell_type = 'Leydig',
#                   features = collagene_genes )
#vp
##################################

# split by condition 

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
                           
###################### ############
                           
##### Vary resolution #######                           
VaryResolution <- function(Seurat.object, 
                           nPC = 20) {
    
    nPC=20
    for (resolu in 1:10) {
        print(resolu/10)
        SO_tmp <- Seurat.object
        SO_tmp <- FindNeighbors(SO_tmp, dims = 1:nPC)
        SO_tmp <- FindClusters(SO_tmp, resolution = resolu/10) 
        RunUMAP(SO_tmp, dims = 1:nPC)
        assign(paste('plot_res', toString(resolu/10), sep=''), 
               DimPlot(SO_tmp, reduction = "umap", label = TRUE)) 
    }
    combined.plot <- CombinePlots(plots = list(plot_res0.1, plot_res0.2, plot_res0.3, 
                                               plot_res0.4, plot_res0.5, plot_res0.6, 
                                               plot_res0.7, plot_res0.8, plot_res0.9))
    return(combined.plot)
}


#res.combined.plot <- VaryResolution(azoospermia.integrated.subset2)
#options(repr.plot.width=17, repr.plot.height=17)
#res.combined.plot  
##################################

##### Save markers in a file and plot Heatmap #######
SaveMarkers <- function(Seurat.object, 
                        filename_xlsx,
                        thresh.use = 0.25,
                        min.pct = 0.25,
                        min.diff.pct = -Inf,
                        test.use = "wilcox",
                        LogFC.onlypos = FALSE) {
    # calculate and save markers and plot Heatmap
    cluster.markers = FindAllMarkers(Seurat.object, 
                                     thresh.use = thresh.use, 
                                     test.use=test.use, 
                                     min.pct=min.pct, 
                                     min.diff.pct=min.diff.pct, 
                                     only.pos=LogFC.onlypos)
    write.xlsx(cluster.markers,
               file= filename_xlsx, 
               row.names = T,
               asTable = T)
    top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    top3  <- cluster.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
    HM <- DoHeatmap(Seurat.object,                     
                    features = top10$gene,
                    disp.min = -2,
                    disp.max = 2,
                    angle = 90,
                   raster = F) +
    scale_fill_gradientn(colours = coolwarm(200)) 
    return(top10)
        
}

#HM_integrated.subset <- SaveMarkers(somatic.integrated,
#                                  'somatic_markers.xlsx')
                           
##################################

##### Plot signature of genes in Feature plot #######
Plot_sign <- function(Seraut.object, signature, operator = sum) {
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
                      cols = c("lightgrey", "red")) +
                      #cols = as.vector(coolwarm(5))) +
    theme(plot.title = element_text(color="black", size=18, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=10, face="italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = 'black', size=16, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 18),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 18),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(title = "Signature plot", subtitle = paste('MarkerGenes: ',toString(signature), sep=''), 
         x = "UMAP 1", y = "UMAP 2") 
    return(FP)
    }
##################################

##### Plot Violin Plot and Feature plot of a gene/list of genes #######
PlotFeature <- function(Seurat.object, features) {    
    DefaultAssay(Seurat.object) <- "RNA"
    vp <- VlnPlot(Seurat.object, 
                  feature = features, 
                  slot = "counts", 
                  log = TRUE)
    fp <- FeaturePlot(Seurat.object, 
                      reduction = "umap", 
                      pt.size = 2,
                      features = features, 
                      label = T, 
                      cols = c("lightgrey", "red"))
    plot <- vp | fp
    return(plot)
}
##################################

##### DGE #####
##### Perform DGE #####
myDGE <- function(Seurat.object, cell_type, outdir = './') {
    azoospermia.response <- FindMarkers(Seurat.object,
                                        ident.1 = paste(cell_type,"_iNOA", sep=''),
                                        ident.2 = paste(cell_type,"_CTL", sep=''), 
                                        verbose = FALSE)
    
    write.xlsx(azoospermia.response,
           file= paste(outdir,'DGE_', cell_type,'.iNOA.vs.CTRL.xlsx',sep=''), 
           row.names = T,
           asTable = T)
    return()    
}


myDGE_Deseq2 <- function(Seurat.object, cell_type) {
    azoospermia.response <- FindMarkers(Seurat.object,
                                        ident.1 = paste(cell_type,"_iNOA", sep=''),
                                        ident.2 = paste(cell_type,"_CTL", sep=''), 
                                        test.use = "DESeq2",
                                        verbose = FALSE)
    
    write.xlsx(azoospermia.response,
           file= paste('DGE_DESeq2_', cell_type,'.azoospermia.xlsx',sep=''), 
           row.names = T,
           asTable = T)
    return()    
}

##### DGE Plot #####
myDGEplot <- function(Seurat.object, cell_type) {
    DefaultAssay(Seurat.object) <- 'RNA'
    Cell.type <- subset(Seurat.object, idents = cell_type)
    Idents(Cell.type) <- "condition"
    avg.CellType <- log1p(AverageExpression(Cell.type, verbose = FALSE)$RNA)
    avg.CellType$gene <- rownames(avg.CellType)
    plot <- ggplot(avg.CellType, aes(azoospermia,healthy)) + geom_point() + ggtitle(cell_type)
    result <- list()
    result[['avg.CellType']] <- avg.CellType
    result[['plot']] <- plot
    return(result)
}

############## end ####################  

##### HEATmap functions #####
##### all cells HEATmap #####
DGE_Heatmap <- function(Seurat.object, 
                        cell_type, 
                        enrichR.path = './',
                        CND = 'up', 
                        database = c('Reactome_2016','GO_Biological_Process_2018'), 
                        pathways.2.plot,
                       gene.label = 6){
#    enrichR.file = paste('/Users/tascini.annasofia/Dropbox (HSR Global)/WORKSPACE/',
#                     'Alfano/',
#                     '904_infertilita_epigenetics/7_bioinfo/Paper_plot/xlsx.table/somatics_DGE_enrichR/',
#                     'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
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
                    raster = F) +
    scale_fill_gradientn(colours = coolwarm(200)) +
    theme(axis.text.y = element_text(color = 'black', size=gene.label))
    return(HM)
}

##### average HEATmap #####
DGE_HeatmapAVE <- function(cluster.averages, 
                        cell_type, 
                        CND = 'up', 
                        database = 'Reactome_2016', 
                        pathways.2.plot){
    
    #enrichR.file = paste('/Users/tascini.annasofia/Dropbox (HSR Global)/WORKSPACE/',
    #                 'Alfano/',
    #                 '904_infertilita_epigenetics/7_bioinfo/Paper_plot/xlsx.table/somatics_DGE_enrichR/',
    #                 'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
    
    enrichR.file = paste('/Users/tascini.annasofia/Dropbox (HSR Global)/',
                         'Alfano_904_paperDraft/',
                         'xlsx.table/somatics_DGE_enrichR/',
                         'enrichR_DGE_',cell_type,'_',CND,'.xlsx', sep='')
    
    enrichR.table <- read.xlsx(xlsxFile = enrichR.file, 
                            sheet = database, 
                            startRow = 1, 
                            colNames = TRUE,
                            rowNames = TRUE, 
                            detectDates = FALSE, 
                            skipEmptyRows = TRUE,
                            skipEmptyCols = TRUE,
                            na.strings = "NA", 
                            fillMergedCells = FALSE)
    gene.list = character()
        for (pathway in pathways.2.plot) {
            gene.list.tmp <- unlist(strsplit(enrichR.table[pathway,]$Genes, ';'))
            gene.list = c(gene.list, gene.list.tmp)
        }
    cell_type.2.use = cell_type
    cell2select <- colnames(cluster.averages)[Idents(cluster.averages) == 
                                              paste(cell_type.2.use,'_azoospermia',sep='') |
                                              Idents(cluster.averages) ==
                                              paste(cell_type.2.use,'_healthy',sep='')]
    cluster.averages.subset <- subset(cluster.averages, 
                                      idents = c(paste(cell_type.2.use,'_azoospermia',sep=''),
                                                 paste(cell_type.2.use,'_healthy',sep='')))
    
    HM <- DoHeatmap(cluster.averages.subset, 
                    features = gene.list,
                    #cells = cell2select,
                    size = 2,
                    disp.min = -3,
                    disp.max = 3,
                    draw.lines = F, 
                    raster = F) +
    scale_fill_gradientn(colours = coolwarm(200))
    return(HM)
}
                           
# subset cell and normalise
SubsetFunction <- function(Seraut.object, 
                           cell_type, 
                           min_nFeature_RNA = 500, 
                           max_nFeature_RNA = 6500, 
                           max_percent_MT = 20,
                          nDIM=10){
    
    SelectedCell = character()
    for (i in 1:length(cell_type)) {
        SelectedCell.tmp <- colnames(Seraut.object)[Seraut.object@meta.data$cell_type == cell_type[i]]
        SelectedCell <- c(SelectedCell, SelectedCell.tmp)
    }        
            
    Cell.subset <- subset(Seraut.object, cells = SelectedCell, invert=F)
    DimPlot(Seraut.object, 
            cells.highlight = SelectedCell, 
            order = T, pt.size = 2, 
            split.by = 'sample')
    
    DefaultAssay(Cell.subset) <- "RNA"
    Cell.subset <- SplitObject(Cell.subset,split.by = "sample")
    print(Cell.subset)
    for (i in 1:length(Cell.subset)) {
        Cell.subset[[i]][["percent.mt"]] <- PercentageFeatureSet(Cell.subset[[i]], pattern = "^MT-")
        Cell.subset[[i]] <-  subset(Cell.subset[[i]],
                                    subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)
        Cell.subset[[i]] <- NormalizeData(Cell.subset[[i]], verbose = FALSE)
        Cell.subset[[i]] <- FindVariableFeatures(Cell.subset[[i]], 
                                                 selection.method = 'vst',
                                                 nfeatures = 2000,
                                                 verbose = FALSE)
        }
    anchors.subset <- FindIntegrationAnchors(object.list = Cell.subset, dims = 1:nDIM, verbose=FALSE)
    integrated.subset <- IntegrateData(anchorset = anchors.subset, dims = 1:nDIM, verbose=FALSE)                                          
 
    DefaultAssay(integrated.subset) <- "RNA"
    all.genes <- rownames(integrated.subset)
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    integrated.subset <- CellCycleScoring(integrated.subset, 
                                           s.features = s.genes, 
                                           g2m.features = g2m.genes, 
                                           set.ident = TRUE,
                                           verbose = FALSE)
    
    DefaultAssay(integrated.subset) <- "integrated"
    integrated.subset <- ScaleData(integrated.subset, 
                                   verbose = FALSE, 
                                   vars.to.regress = c("percent.mt", "nFeature_RNA"),
                                   features = all.genes)

    integrated.subset <- ScaleData(integrated.subset,
                                   verbose = FALSE,
                                   vars.to.regress = c("S.Score", "G2M.Score"), 
                                   features = all.genes) 
    return(integrated.subset)
    }
                           
