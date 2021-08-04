# functions

Gene_conversion = function(x, gtf_dictionary, gtf_OSR_dictionary){
  # convert gene name from two different versions
  gene.ensembl =  gtf_dictionary[which(gtf_dictionary$gene_name == x),]$gene_id
  genec = gtf_OSR_dictionary[gtf_OSR_dictionary$gene_id %in% gene.ensembl,]$gene_name
  if (x %in% genec | length(genec) == 0) {genec = x}
  newgene  = paste(genec, collapse = ';')
  return(newgene)
}


myDGE = function(Seurat.object, cell_type, outdir = './', 
                 c1 = "iNOA", c2 = "CTL") {
  azoospermia.response <- FindMarkers(Seurat.object,
                                      ident.1 = paste(cell_type,"_", c1, sep=''),
                                      ident.2 = paste(cell_type,"_", c2, sep=''), 
                                      verbose = FALSE)
  write.xlsx(azoospermia.response,
             file= paste(outdir,'DGE_', cell_type,'.',c1,'.vs.',c2,'.xlsx',sep=''), 
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
    labs(title = title, subtitle = paste('MarkerGenes: ',toString(signature), sep=''), 
         x = "UMAP 1", y = "UMAP 2") 

  return(FP)
}

