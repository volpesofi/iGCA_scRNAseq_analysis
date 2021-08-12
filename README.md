# iGCA_scRNAseq_analysis
Analysis script of "Aging, inflammation and DNA damage in the somatic testicular niche with idiopathic germ cell aplasia"

The single testis samples were analysed following the standard Seurat workflow [https://satijalab.org/seurat/articles/pbmc3k_tutorial.html]
and integrated with the integration workflow [https://satijalab.org/seurat/articles/integration_introduction.html].

We used Seurat program (v.3.1.5), within the R environment (v3.6.1).

### Seurat analysis and integration ###

In the ``Seurat`` folder you can find the main analysis script:
  * Single samples iNOA/OA analysis: ``seurat.analysis.singleSamples.ipynb``
  * integration of iNOA samples: ``azoospermia.integration.seurat.ipynb``
  * integration of adults samples (3 iGCA + CTL from GSE): ``iGCA_CTL_integration.R``

### Analysis of testis samples from the literature ###
In the ``GEOdatasets`` you can find a script for analysis of testis samples from the literature (several GSE accession codes) and identification of somatic cell: ``testis_dataset_all.R``

### Leydig cell development ###
In the ``Development`` folder you can find:
  * the script for the integration of neonatal, prepubertal, adults and iGCA samples: ``development.R``
  * the script with the definition of the Stage signatures A, B, C and their expression levels in iGCA and CTL samples: ``development_StagesLEY.R``

### T cell characterisation ###
In ``Tcell_characterisation`` folder you can find the script for Tcell characterisation according to GSE126030: ``Tcell.R``

### CellPhoneDB analysis ###
In the ``Receptor_Ligand_interactions`` folder you can find the script to prepare files for CellPhoneDB analysis ans plot the CellPhoneDB results: ``CellPhoneDB_prepare_plots.R``

## Plotting with ggplot2 ##
In the ``PaperPlot`` folder you can find the script for reproducing the paper figures:  ``plotsNatComm.R``

## Utility functions ##
The scripts ``integration_function.R``, ``mySeuratfunctions.R``, ``utility_functions.R`` are a set of functions (mainly wrapper of Seurat functions) I use for the analysis workflow.
