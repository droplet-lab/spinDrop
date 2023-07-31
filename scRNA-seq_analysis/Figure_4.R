library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(reshape2)
library(biomaRt)
library(ggrepel)
library(corrplot)
library(colorspace)

## Read seurat objects with 10x, sciRNA-seq, spinDrop sorted_not_sorted, technology stored in metadata$data
obj <- readRDS("040722_integrated_data_brain_annotated.rds")

## Integrate the datasets using Seurat v3
obj.list <- SplitObject(obj, split.by = "data")

for (i in names(obj.list)) {
  obj.list[[i]] <- SCTransform(obj.list[[i]], verbose = FALSE)
}

obj.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = obj.features)

# chose 10x as the reference dataset containing the annotations
reference_dataset <- which(names(obj.list) == "10X")

obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", 
                                       anchor.features = obj.features, reference = reference_dataset)
obj.integrated2 <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT")

# dimensional reduction
obj.integrated2 <- RunPCA(object = obj.integrated2, verbose = FALSE)
obj.integrated2 <- RunUMAP(object = obj.integrated2, dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% FindClusters(resolution =0.5) %>% identity()

# cell type annotation
Idents(obj.integrated) <- obj.integrated@meta.data$data
obj.subset <- subset(obj.integrated, idents = "10X")

#Plot correlation heatmap
x <- AverageExpression(obj.subset, group.by="Class")
y <- AverageExpression(obj.subset, group.by="seurat_clusters")
coeff_cor <-cor(x$SCT, y$SCT)

Idents(obj.subset)<- obj.subset@meta.data$Class
DimPlot(obj.subset, cols = colours, group.by=c("Class","seurat_clusters"))

coeff_cor2 <- as.data.frame(coeff_cor)

obj.integrated@meta.data$celltype <- obj.integrated@meta.data$seurat_clusters
Idents(obj.integrated) <- "celltype"

#transfer cell type labels to the combined object
celltypes <- rownames(coeff_cor2)[apply(coeff_cor2, 2, which.max)]
names(celltypes) <- levels(obj.integrated)
obj.integrated <- RenameIdents(obj.integrated, celltypes)

obj.integrated@meta.data$celltype_n <-obj.integrated@active.ident
DimPlot(obj.integrated, split.by="data",  cols=colours, pt.size = 0.2)

## Plot populations of cells
d2 <- obj@meta.data %>% 
  group_by(data, celltype_n) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))


colours_rand <-sample(colours)

ggplot(d2, aes(x = factor(data), y = perc*100, fill = factor(celltype_n))) +
  geom_bar(stat="identity", width = 0.7, colour="black") +
  labs(x = "type", y = "percent (%)", fill = "Cell type") +
  theme_classic() + scale_fill_manual(values=colours)

##DEG computation
obj.list2 <- SplitObject(obj.integrated, split.by = "data")
for (i in names(obj.list2)) {
  obj.list2[[i]] <- SCTransform(obj.list2[[i]], verbose = FALSE)
}


obj.comp <- merge(x= obj.list2[["10X"]], y= obj.list2[["Brain_E10_5+"]])

obj.integrated2 <- SCTransform(obj.comp, verbose = FALSE)


ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")


annotations <- getBM(attributes = c('external_gene_name','gene_biotype'),
                mart = ensembl)

cells <-unique(obj.integrated2@meta.data$celltype_n)
list_diff_genes <-  rep(NA, length(unique(obj.integrated2@meta.data$celltype_n)))

for (i in 1:length(unique(obj.integrated2@meta.data$celltype_n))) {
  cm <- subset(obj.integrated2, idents= unique(obj.integrated2@meta.data$celltype_n)[[i]])
  print(unique(obj.integrated2@meta.data$celltype_n)[[i]])
  Idents(cm) <- "data"
  ld <- FindMarkers(cm, ident.1 = "Brain_E10_5+", ident.2 = "10X",logfc.threshold =0.05,min.pct=0.1, min.cells.group=10)
  ld$biotype <- annotations$gene_biotype[match(rownames(ld), annotations$external_gene_name)]
  write.csv(ld,paste(unique(obj.integrated2@meta.data$celltype_n)[[i]],".csv",sep=""))
  list_diff_genes[i+1] <- length(ld$p_val)
}


## check differences between all technologies
downsampled_dataset <- merge(x=spindrop_down,y=tenx_down)
ls <- FindMarkers(obj, ident.1 = "10X", indent.2 ="Brain_E10_5+")
ls$biotype <- annotations$gene_biotype[match(rownames(ls), annotations$external_gene_name)]

## Plot Neuroblast data
rs <- read.csv("Neuroblast.csv")

rs$diffexpressed <- "no"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
rs$diffexpressed[(rs$avg_log2FC > 1 | rs$avg_log2FC < (-1))] <- "yes"

rs$delabel <- NA
rs$delabel[rs$diffexpressed != "no"] <- rs$X[rs$diffexpressed != "no"]

ggplot(data=rs, aes(x=avg_log2FC, y=-log10(p_val_adj), col=biotype, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 10)+
  scale_color_manual(values=colours) +
  geom_vline(xintercept=c(-1, 1), col="red") 

ggplot(data=rs, aes(x=avg_log2FC, y=-log10(p_val_adj))) + 
  geom_point() + 
  theme_minimal() 

## biotype plot
Idents(obj)<- obj@meta.data$data

tenx <- subset(obj, ident="10X")
spindrop<- subset(obj, ident="Brain_E10_5+")
sciRNA <- subset(obj, ident="sciRNA")

tenx_down <- tenx[, WhichCells(tenx,  downsample = 2000)]
spindrop_down <- spindrop[, WhichCells(spindrop,  downsample = 2000)]
sciRNA_down <- sciRNA[, WhichCells(sciRNA,  downsample = 2000)]

xz <- data.frame(spindrop = AverageExpression(spindrop_down)$SCT,
                 tenx = AverageExpression(tenx_down)$SCT,
                 sciRNA = AverageExpression(sciRNA_down)$SCT)
colnames(xz) <- c("spindrop", "tenx","sciRNA")

xz$genes <- rownames(xz)
xz$biotype <- annotations$gene_biotype[match(rownames(xz), annotations$external_gene_name)]

tenx_bio <- xz %>%
  group_by(biotype) %>%
  summarise(n = n_distinct(tenx))

spindrop_bio <- xz %>%
  group_by(biotype) %>%
  summarise(n = n_distinct(spindrop))

sciRNA_bio <- xz %>%
  group_by(biotype) %>%
  summarise(n = sum(n_distinct(sciRNA)))

bioc_count <- xz %>%
  group_by(biotype) %>%
  summarise(tenx= sum(tenx), sciRNA=sum(sciRNA), spindrop=sum(spindrop))
