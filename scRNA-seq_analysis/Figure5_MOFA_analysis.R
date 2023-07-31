library(data.table)
library(purrr)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(MOFA2)
library(ggrepel)
library(ggpubr)
library(GGally)

colours=c("#FFFF00","#1CE6FF", "#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43",
          "#00C2A0","#FFAA92","#FF90C9", "#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8","#013349","#00846F",
          "#372101","#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75",
          "#B77B68","#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD",
          "#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757")

## Load old RNA matrix
old<- read.csv("E85_old_curated.csv", header=TRUE, row.names=1)

## Create Seurat object
seurat_old <- CreateSeuratObject(counts = old, min.cells = 3, assay = "RNA",names.field = 1,
                            names.delim = "-",
                            project = "old")

## Load metadata
metadata<- read.csv("nascent_old_new_metadata.csv", header =1)

seurat_old@meta.data$cell_type <-metadata$cell_type

## Load nascent RNA matrix
new<- read.csv("E85_new_curated.csv", header=TRUE, row.names=1)

seurat_new <- CreateSeuratObject(counts = new, min.cells = 3, assay = "RNA",names.field = 1,
                            names.delim = "-",
                            project = "new")

## Normalize and pre-process independently
seurat_old <- NormalizeData(seurat_old, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_new <- NormalizeData(seurat_new, normalization.method = "LogNormalize", scale.factor = 10000)

barcodes <- intersect(colnames(seurat_new), colnames(seurat_old))
seurat_new <- seurat_new[,barcodes]
seurat_old <- seurat_old[,barcodes]

seurat_old <- FindVariableFeatures(seurat_old, 
                                   selection.method = "vst", 
                                   nfeatures = 2500
)
old_rna.features <- seurat_old@assays$RNA@var.features

seurat_new <- FindVariableFeatures(seurat_new, 
                                   selection.method = "vst", 
                                   nfeatures = 2500
)
new_rna.features <- seurat_new@assays$RNA@var.features

plot_variance_explained(model, x="view", y="factor")


## Compute factors
mofa <- create_mofa(list(
  "oldRNA" = as.matrix( seurat_old@assays$RNA@data[old_rna.features,] ),
  "newRNA" = as.matrix( seurat_new@assays$RNA@data[new_rna.features,] )
))

mofa

## Cap the number of factors at 6
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 6

# Training options: let's use default options
train_opts <- get_default_training_options(mofa)
train_opts$seed <- 42

mofa <- prepare_mofa(
  object = mofa,
  model_options = model_opts,
  training_options = train_opts
)


mofa <- run_mofa(mofa,"old_new_MOFA.hdf5",use_basilisk = TRUE)

# transfer cell type labels
mofa@samples_metadata$cell_type <- metadata$cell_type

plot_variance_explained(mofa)
plot_variance_explained(mofa, plot_total = TRUE)[[2]]


factors <- 1:get_dimensions(mofa)[["K"]]
factors <- factors[!factors%in%c(4,7)]

mofa <- run_umap(mofa)
plot_dimred(mofa, method="UMAP",color_by="cell_type")+ scale_fill_manual(values=colours)

## plot all weights for all factors
plot_weights(mofa, 
             view = "newRNA", 
             factors = 2, 
             nfeatures = 5, 
             text_size = 2
)

plot_top_weights(mofa,
                 view = "newRNA",
                 factor = 2,
                 nfeatures = 10,     
                 scale = T           
)

plot_data_heatmap(mofa, 
                  view = "newRNA",
                  factor = 2,  
                  features = 25,
                  denoise=TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

plot_factor(mofa, factors=1, group_by = "cell_type", color_by="cell_type") +
  theme(
    axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1)
  )

plot_factors(mofa, factors=1:4, color_by = "cell_type")

## Plot the variance explained by each factor
plot_variance_explained(model, x="view", y="factor")