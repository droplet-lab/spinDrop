library(harmony)
library(Seurat)
library(biomaRt)
library(ggrepel)

setwd("~/Desktop/scRNA-seq_analysis/matrices")
# load seurat objects and annotate genes for FADS-sorted samples
obj <- readRDS('B_cell.dgecounts.rds')
matrix <- as.matrix(obj$umicount$inex$all)

annot <- getBM(attributes = c('external_gene_name','gene_biotype'),
                     mart = ensembl)

rownames(matrix) <- annot$external_gene_name[match(rownames(matrix), annot$ensembl_gene_id)]
B_cells <- CreateSeuratObject(matrix,min.cells = 3, assay = "RNA",
                              project = "B_cell", min.features = 300)

B_cells@meta.data$type <- "indrop"
B_cells@meta.data$subtype <- "B cells"

# load mouse PBMC 10x object
tenx <- Read10X("mouse_PBMC_10X")

tenx_s <- CreateSeuratObject(tenx,min.cells = 3, assay = "RNA",names.field = 1,
                             names.delim = "_",
                             project = "tenx", min.features = 500)

tenx_s@meta.data$type <- "10x"
tenx_s@meta.data$subtype <- "10x"

pbmc <- merge(x= tenx_s ,y= c(B_cells))


# integrate datasets
obj.list <- SplitObject(pbmc, split.by = "type")

for (i in names(obj.list)) {
  obj.list[[i]] <- SCTransform(obj.list[[i]], verbose = FALSE)
}

obj.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = obj.features)

# Use 10x as the reference dataset
reference_dataset <- which(names(obj.list) == "10X")

obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", 
                                      anchor.features = obj.features)
obj.integrated2 <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT")

#dimernsional reduction
obj.integrated2 <- RunPCA(object = obj.integrated2, verbose = FALSE)
obj.integrated2 <- RunHarmony(obj.integrated2, c("subtype"),assay.use = "integrated")
obj.integrated2 <- RunUMAP(obj.integrated2, reduction = "harmony", dims = 1:10) %>% FindNeighbors(reduction = "harmony", dims = 1:50) %>% FindClusters(resolution=2) %>% identity() 

#Cluster ID	Markers	Cell Type
#0	IL7R, CCR7	Naive CD4+ T
#1	CD14, LYZ	CD14+ Mono
#2	IL7R, S100A4	Memory CD4+
#3	MS4A1	B
#4	CD8A	CD8+ T
#5	FCGR3A, MS4A7	FCGR3A+ Mono
#6	GNLY, NKG7	NK
#7	FCER1A, CST3	DC
#8	PPBP	Platelet
DotPlot(obj.integrated2, features=c("Il7r","Ccr7","Cd14","S100a4","Ms4a1","Cd8a","Fgr3a","Ms4a7","Nkg7","Fcer1a","Cst3","Ppbp", "S100a9", "Lst1",
                                    "Hla-dpb1"))

DotPlot(obj.integrated2, features=c("Cst3", "Il7r","Lst1","Nkg7","S100a9","Cd8a","Ms4a1","Ccr7"))

DimPlot(obj.integrated2, cols=colours)

new.cluster.ids <- c("Naive CD4+ T cells", "B cells", "B cells", "B cells", "CD8+ T cells",
                     "B cells", "CD14+ Mono", "NK cells", "FCGR3A+ Mono", "Memory CD4+ T cells", "DC cells","B cells", "CD8+ T cells", "B cells")

## Annotate the clusters
names(new.cluster.ids) <- levels(obj.integrated2)
obj.integrated2 <- RenameIdents(obj.integrated2, new.cluster.ids)

