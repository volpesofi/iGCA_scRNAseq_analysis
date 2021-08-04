# integration function
integration_workflow <- function(somatic.list,   
                                 dataset_iNOA_dir = '.', 
                                 dataset = '',
                                 min_nFeature_RNA = 500, 
                                 max_nFeature_RNA = 6500, 
                                 max_percent_MT = 20,   
                                 res = 0.2, nPC = 20){
  dir.create(dataset_iNOA_dir, recursive = TRUE)
  
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
  somatic.anchors <- FindIntegrationAnchors(object.list = somatic.list, dims = 1:30, verbose = TRUE)
  somatic.integrated <- IntegrateData(anchorset = somatic.anchors, dims = 1:30, verbose = TRUE)
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
                                         verbose = FALSE)
  
  
  DefaultAssay(somatic.integrated) <- 'integrated'
  
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
  
  E = ElbowPlot(somatic.integrated)
  ggsave(filename = paste(dataset_iNOA_dir,"ElbowPlot_",dataset,".png",sep=''),
         plot = E,
         width = 10, height = 7,
         units = 'cm')
  
  #UMAP
  somatic.integrated <- RunUMAP(somatic.integrated, reduction = "pca", dims = 1:20, verbose = FALSE)
  # clusterinig
  somatic.integrated <- FindNeighbors(somatic.integrated, dims = 1:nPC, verbose = FALSE)
  somatic.integrated <- FindClusters(somatic.integrated, resolution = res, verbose = FALSE)
  somatic.integrated <- RunUMAP(somatic.integrated, dims = 1:nPC, verbose = FALSE)
  
  
  bys = DimPlot(somatic.integrated, reduction = "umap", label = T, pt.size = 1, split.by = 'source', order = TRUE) + 
    theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=22, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 24),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
          axis.title.y = element_text(face = "bold", color = "black", size = 24),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position="top",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(subtitle = paste('Res = ', res,', nPC = ',nPC, sep = ''), x = "UMAP 1", y = "UMAP 2") + NoLegend()
  
  ggsave(filename = paste(dataset_iNOA_dir,"UMAP_bysources_allcell_integration_",dataset,"_iNOA_res",res,".png",sep=''),
         plot = bys,
         width = 45, height = 15, 
         units = 'cm')
  
  
  filename = paste("integration_",dataset,"_iNOA_res",res,".xlsx", sep = '')
  existing_file = filename %in% list.files(path = dataset_iNOA_dir, pattern = ".xlsx")
  if (!existing_file) {
    
    cluster.markers = FindAllMarkers(somatic.integrated, verbose=FALSE)
    top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    write.xlsx(cluster.markers,
               file = paste(dataset_iNOA_dir, filename, sep =''), 
               row.names = TRUE,
               asTable = TRUE)
    
    write.csv(cluster.markers, 
              paste(dataset_iNOA_dir, 'seurat_MG_integration_',dataset,'_iNOA_res',res,'.csv', sep =''))
    
    HM =  DoHeatmap(somatic.integrated, features = top10$gene, angle = 90) + NoLegend() + 
      scale_fill_gradientn(colours = coolwarm(200)) 
    ggsave(filename = paste(dataset_iNOA_dir,"HM_integration_",dataset,"_iNOA_res",res,".png",sep=''),
           plot = HM,
           width = 14, height = 18, 
           units = 'cm')
  }
  
  
  leydig_signature = c('CFD','DLK1','LUM','CALB2')
  myoid_signature = c('ACTA2','MYH11','DES', 'MYL9')
  sertoli_signature = c('FATE1','CITED1','SOX9','AMH','CLDN11')
  macrophage_signature = c('CD14','CD74','HLA-DRA','HLA-DRB1')
  endothelial_signature = c('VWF', "EGFL7",'CD34', 'PRSS23', "RBP7")
  t_cell_signature =  c('GZMA','GZMK','CD2','CCL5','NKG7')
  
  DefaultAssay(somatic.integrated) <- 'RNA'
  pL <- Plot_sign(somatic.integrated,
                  signature= c('CFD','DLK1','LUM'), 
                  operator = mean, title = 'LEY')
  pM <- Plot_sign(somatic.integrated,
                  signature= c('ACTA2','MYH11','DES'), 
                  operator = mean, title = 'MYD')
  pS <- Plot_sign(somatic.integrated,
                  signature= c('FATE1','CITED1','SOX9', 'AMH'), 
                  operator = mean, title = 'SRT')
  pMa <- Plot_sign(somatic.integrated,
                   signature= c('CD14','CD74','HLA-DRA'), 
                   operator = mean, title = 'MCR')
  pE <- Plot_sign(somatic.integrated,
                  signature= c('VWF',"EGFL7", 'PRSS23'), 
                  operator = mean, title = 'END')
  pT <- Plot_sign(somatic.integrated,
                  signature= c('GZMA','CD2','CCL5'), 
                  operator = mean, title = 'TCL')
  pSTRO <- Plot_sign(somatic.integrated,
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
  png(paste(dataset_iNOA_dir, "MG_integration_",dataset,"_iNOA_res",res,".png", sep=''),  width=1400, height=600)
  print(wrap_plots(A = pL, B = pM, C = pS, D = pSTRO, E = pE, H = pMa, G =pT, design = layout))
  dev.off()
  
  DefaultAssay(somatic.integrated) <- 'integrated'
  return(somatic.integrated)
}



integration_postprocessing <- function(SO, new_cluster_id, 
                                       dataset = '', dataset_iNOA_dir = '.',
                                       TCL_DGE = FALSE,
                                       plotTCL = FALSE,
                                       c1 = "iNOA", c2 = "CTL",
                                       res = 0.2, colors_Alf) {
  
  DefaultAssay(SO) <- "integrated"
  names(new_cluster_id) <- levels(SO)
  SO <- RenameIdents(SO, new_cluster_id)
  
  col = as.character(colors_Alf[levels(SO),]$cols)
  
  png(paste(dataset_iNOA_dir,"UMAP_cellAssigned_",dataset,"_res",res,".png",sep=''),  
      width = 500, height = 400)
  print(DimPlot(SO, reduction = "umap", label = TRUE, pt.size = 1, order = TRUE) +
    scale_color_manual(values = (col)) +
    theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=22, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 24),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
          axis.title.y = element_text(face = "bold", color = "black", size = 24),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "UMAP 1", y = "UMAP 2"))
  dev.off()
  
  filename = paste("integration_",dataset,"_res",res,"_celltype.xlsx", sep = '')
  existing_file = filename %in% list.files(path = dataset_iNOA_dir, pattern = ".xlsx")
  
  if (!existing_file) {
    cluster.markers = FindAllMarkers(SO)
    top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    write.xlsx(cluster.markers,
               file= paste(dataset_iNOA_dir, filename, sep =''), 
               row.names = T,
               asTable = T)
    HM =  DoHeatmap(SO, features = top10$gene, angle = 90) + NoLegend() + 
      scale_fill_gradientn(colours = coolwarm(200)) 
    
    png(paste(dataset_iNOA_dir,"HM_celltype_integration_",dataset,"_res",res,".png",sep=''),  width = 1400, height = 1600)
    plot(HM)
    dev.off()
  }
  
  ######## DGE ###########
  dir_DGE = paste(dataset_iNOA_dir,'DGE/', sep='')
  dir.create(dir_DGE)
  
  DefaultAssay(SO) <- "RNA"
  # relabel cells to perform DGE
  SO$celltype.cond <- paste(Idents(SO), SO$condition, sep = "_")
  SO$celltype <- Idents(SO)
  Idents(SO) <- "celltype.cond"
  if (!TCL_DGE) {new_cluster_id = new_cluster_id[new_cluster_id != 'TCL']}
  DGEresults = lapply(new_cluster_id, 
                      myDGE, 
                      Seurat.object = SO, 
                      outdir = dir_DGE,
                      c1 = c1,
                      c2 = c2)
  filename = paste(dataset_iNOA_dir,"DGEall_integration_",dataset,"_res",res,".xlsx", sep = '')
  write.xlsx(x = DGEresults, file = filename, asTable = T, row.names = T)
  
  ######## enrichment ###########
  dir_enrichR = paste(dataset_iNOA_dir,'enrichR/', sep='')
  dir.create(dir_enrichR)
  
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
  
  
  cell_types = new_cluster_id
  DGE_gene_up_list <- list()
  DGE_gene_dw_list <- list()
  
  for (celltype in cell_types) {
    DGE_file = paste(dir_DGE,"DGE_", celltype, ".",c1,".vs.",c2,".xlsx", sep= '')
    A <- read.xlsx(DGE_file)
    
    up.genes <- A$row.names[A$avg_logFC > 0 & A$p_val_adj <= 0.05]
    down.genes <- A$row.names[A$avg_logFC < 0 & A$p_val_adj <= 0.05]
    both.genes <- A$row.names[A$p_val_adj <= 0.05]
    
    DGE_gene_up_list[[celltype]] <- up.genes
    DGE_gene_dw_list[[celltype]] <- down.genes
    
    enrichr.list <- list()
    enrichr.list <- lapply(list(up.genes,down.genes,both.genes), function(x) {
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
      if (!is.null(enrichr.list[[i]])) {write.xlsx(x = enrichr.list[[i]], file = filename, asTable = T)}
    }
  }
  
  if (plotTCL) {
  Idents(SO) <- "celltype"
  SO$condition <- factor(x = SO$condition, levels = c(c1, c2))
  #SO  = somatic.integrated.new2
  for (g in  c('CD3E','TRAC', 'GZMK')) {
    v = VlnPlot(subset(SO, idents = c('MCR','TCL')), features =g, 
                split.by = 'condition', assay = 'RNA', split.plot = F, pt.size = 0.5) + 
      theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
            plot.subtitle = element_text(color="black", size=16, face="italic"),
            axis.text.x = element_text(angle = 0, face = "bold", color = "black", size = 16, hjust =.5), 
            axis.title.x = element_text(face = "bold", color = "black", size = 18),
            axis.text.y = element_text(angle = 0, face = "bold", color = 'black', size=16),
            axis.title.y = element_text(face = "bold", color = "black", size = 18),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.position = c(0.1, 0.92),
            panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) #+
      #geom_boxplot(width = 0.1, outlier.size=1) +
      #scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c(c1,c2))
    assign(paste('v_', g, sep =''),v)
  }
  v_CD3E | v_TRAC | v_GZMK
  
  png(paste(dataset_iNOA_dir,"MG_Tcell_",dataset,"_res",res,".png",sep=''),  width = 1400, height = 600)
  plot(v_CD3E | v_TRAC | v_GZMK)
  dev.off()
  
  
  for (g in  c('HLA-DRA','CD74', 'HLA-DPA1')) {
    v = VlnPlot(subset(SO , idents = c('MCR','TCL')), features =g, 
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
      scale_fill_manual(values = c('orange', 'dodgerblue2'), labels = c(c1,c2))
    assign(paste('v_', g, sep =''),v)
  }
  
  `v_HLA-DRA` | v_CD74 | `v_HLA-DPA1`
  png(paste(dataset_iNOA_dir,"MG_MCR_",dataset,"_res",res,".png",sep=''),  width = 1400, height = 600)
  plot(`v_HLA-DRA` | v_CD74 | `v_HLA-DPA1`)
  dev.off()
  }
  
  return(SO)
}