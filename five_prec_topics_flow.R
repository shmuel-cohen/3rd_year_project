
################################################################################
## Globals stuff and libreries   
################################################################################

K <- 10
NC <- 6

library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(hexbin)
library(viridis)
library(Seurat)
library(SeuratDisk)
library(SeuratData)

library(gridExtra)
install.packages("Seurat")

#install.packages("hexbin")
################################################################################
## open objects   
################################################################################

obj <- LoadH5Seurat("/Users/itamar_shahar/Downloads/objects/SuperAgersRandomSubset0.05.h5seurat")
## fsad
five_prec_clustered <- LoadH5Seurat("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/five_prec.h5seurat")
## libreraries 
set.seed(1)

################################################################################
## transpose count 
################################################################################
trans <- five_prec_clustered@assays$RNA@counts
trans <- t(trans)
five_prec_clustered@assays$RNA@counts <- trans

################################################################################
## feel the data: 
################################################################################

counts <- five_prec_clustered@assays$RNA@counts
#head(five_prec_clustered@assays$RNA@counts, 20)
#dim(five_prec_clustered@assays$RNA@counts)

#dim(counts)

#mean(counts>0)


################################################################################
## fit the topic model
################################################################################
fit <- fit_topic_model(counts, 
                       k = K,
                       control.init = list(nc = NC),
                       control.main = list(nc = NC),
                       control.refine = list(nc = NC)
                       )
#fit_topic_model
#fitted_k_10 <- obj@fit

saveRDS(fit, 
        paste0("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/fitted_k_", K, ".rds"))
        #head(fit$F)  



################################################################################
## upload the fiitted model 
################################################################################

# Load the Seurat object from the RDS file
fit <- readRDS("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/fitted_k_10.rds")


################################################################################
## Adding the topics as column to the obj@meta.data  
################################################################################



new_meta_data <- cbind(five_prec_clustered@meta.data, fit$L)

five_prec_clustered@meta.data <- new_meta_data

five_prec_clustered@meta.data$celltype <- five_prec_clustered@meta.data$new_col
five_prec_clustered@meta.data$new_col <- NULL
SaveH5Seurat(five_prec_clustered,
             filename="five_prec_clustered_and_topics", 
             overwrite = TRUE)
################################################################################  
## visualize the plans (topics) for each cell 
################################################################################
df <- data.frame(row.names = seq(from = 1, to = length(five_prec_clustered$seurat_clusters)),
                 CellID = rownames(five_prec_clustered@meta.data),
                 Cluster = as.factor(five_prec_clustered$new_col))
additional_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                       "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                       "#cab2d6", "#6a3d9a")
structure_plot(fit,
               topics = 1:K,
               colors = additional_colors,
               #colors =scale_fill_viridis_d,
               gap = 25,
               #ticks = 
               #grouping = five_prec_clustered$seurat_clusters,
              #+  scale_colour_viridis_c()
               grouping = df$Cluster
               )

ggsave(paste0("topics_k_10", ".pdf"),
       #plot = g,
       #width = 35, 
       #height = 15,
       limitsize = FALSE,
       path="/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/plots")
      
################################################################################  
## FeaturePlot - umap 
################################################################################
custom_feature_names <- paste0("k", 1:K)
list(custom_feature_names)

g<-FeaturePlot(object = five_prec_clustered,
            #features = list(custom_feature_names)
            features =custom_feature_names,
            label = TRUE, 
            combine = F,
            #blend = TRUE,
            )# + labs(title ="Cells Type and the Average Expression of Over Different Topics") 
remove_axes <- function(plot) {
  plot + theme_void()
}

# Apply the function to each plot in the list
g_no_axes <- lapply(g, remove_axes)

# Combine the modified plots into a single plot
g <- do.call(grid.arrange, c(g_no_axes, ncol = 3))  # Change ncol as needed



?FeaturePlot

ggsave(paste0("topics_k_", K, "_distribution_over_clusters.pdf"),
       plot = g,
       width = 35, 
       height = 15,
       limitsize = FALSE,
       path="/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/plots")  

################################################################################  
## DotPlot
################################################################################


DotPlot(object = five_prec_clustered,
        assay	= "RNA",
        features = custom_feature_names,
        scale.min = 0,
        scale.max = 100,
        scale = F,
        #group.by =  "seurat_clusters",       #cols =  c("#f7fbff", "#deebf7", "#c6dbef", "#9ecae1","#6baed6", "#4292c6", "#2171b5", "#08519c","#08306b", "#08305b"))
        group.by =  "celltype",       #cols =  c("#f7fbff", "#deebf7", "#c6dbef", "#9ecae1","#6baed6", "#4292c6", "#2171b5", "#08519c","#08306b", "#08305b"))
) + labs(x = "Topics", y= "Cell Type", title="Cells Type and the Average Expression of Over Different Topics ") + theme(text = element_text(family = "Arial"))



?DotPlot
ggsave(paste0("topics_k_", K, "_DotPlot.pdf"),
       #plot = g,
       #width = 35, 
       #height = 15,
       limitsize = FALSE,
       path="/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/plots")


head(five_prec_clustered@assays$RNA@counts[,"MIR1302-2HG"],10)

################################################################################  
## ViolinPlot
################################################################################
