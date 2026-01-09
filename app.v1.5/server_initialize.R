
# the annotation table 
setwd(wd)

## aws variables
region = "us-east-1"
bucket_name= "forjazul.galaxy"
ora_tab_list <- NULL
norm_and_de_tab_list <- NULL
custom_upset_terms_list <- c()
volcano_plot_data_list <- list()
heatmap_plot_data_list <- list()
pca_data_list  <- list()
fdr_filter <- 0
fc_filter_low <- 0
fc_filter_high <- 0
correct_metadata_input = T
precomp_metadata <- data.frame()
fastqc_files <- data.frame()
galaxy_outputs_choices_list <- list()
for(run_id in list.files("metadata", pattern = ".csv$", full.names = T)){
  precomp_metadata <- rbindlist(list(precomp_metadata, 
                                     as.data.frame(read.table(run_id, header = T , sep = "," )) %>%
                                       select(BioProject, Run, SampleName) 
  ))
}
precomp_metadata <- as.data.frame(precomp_metadata)


MSI_data_folder = "/mnt/storage/alexandros/bcro/data/central_efs/data/users/Public/tmp/kallisto/v_union2_transcriptome.fasta/MSI/tx_level_DE/"


degs <- c()
rv <- reactiveValues(gsea_plots = NULL,
                     IP = NULL, 
                     show_menu_fastqc = F,
                     show_norm_and_de_mout = NULL,
                     show_clustering_mout = NULL,
                     show_gso_mout = NULL,
                     show_KeggScape_mout = NULL,
                     show_gsea_mout = NULL,
                     show_gsea_plots_mout = NULL,
                     show_ora_mout = NULL,
                     
                     runs_galaxy_df = data.frame(),
                     galaxy_outputs = data.frame(),
                     
                     condition = "",
                     conditions = c(),
                     
                     update_annotations = 0,
                     update_msi_df = 0,
                     update_msi_clustering = 0,
                     
                     norm_and_de_signal =0,
                     data = NULL, 
                     text = NULL, 
                     okplot = FALSE, 
                     #okplot_KEGGmap_ora = FALSE,
                     okplot_KEGGmap_gso = FALSE,
                     okplot_KEGGmap_gsea = FALSE,
                     okplot_gseaplot = FALSE,
                     map = "",
                     map_gsea = "",
                     go_term = "",
                     run_counter = 0, 
                     
                     gse = NULL,
                     gse_desc = NULL,
                     
                     # barplot1 = NULL,
                     # barplot2 = NULL,
                     # dotplot = NULL,
                     # emapplot = NULL,
                     # cnetplot = NULL,
                     # heatplot = NULL,
                     # treeplot = NULL,
                     # upsetplot = NULL,
                     
                     # dotplot_gsea = NULL,
                     # ridgeplot_gsea = NULL,
                     
                     filtered = data.frame(), 
                     #pheatmap = NULL,
                     #volcano0 = NULL,
                     pca_plot = NULL,
                     heatmap_overview_plot = NULL,
                     
                     all_cluster_heatmap = NULL,
                     all_together_tcpm = c(),
                     
                     gsea_list = list(),
                     #gsea = NULL,
                     current_gsea = NULL,
                     #gseaplot = NULL,
                     
                     enriched_kegg_pathways = NULL,
                     #enriched_kegg_pathways_desc = NULL,
                     #enriched_kegg_pathways_id = NULL,
                     submission_by = "",
                     input_genes = c(),
                     input_data = NULL,
                     non_redundant_input_data_size = NULL,
                     input_genes_with_annotations = c(),
                     
                     kegg_ora_vector = NULL,
                     kegg_gsea_vector = NULL,
                     
                     
                     ora = data.frame(), 
                     msi_df = data.frame(), 
                     gsea_input = NULL,
                     gsea = data.frame(),
                     gsea_kegg = data.frame(),
                     go = NULL,
                     enricher_plots = NULL,
                     
                     fc = c(),
                     
                     dim_msi_df = NULL ,
                     
                     kegg_gso_col_low = "red",
                     kegg_gso_col_high = "green",
                     kegg_gsea_col_low =  "red",
                     kegg_gsea_col_high = "green",
                     
                     aws_files = c(),
                     metadata_df = data.frame(),
                     
                     fastqc_run_acc = 0
                     
)