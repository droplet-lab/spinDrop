##Scripts to run DropletQC on the FADS sorted and unsorted mouse brain samples
## address queries to marcin.tabaka@gmail.com


require(ggplot2)
require(data.table)
require(DropletQC)
require(Seurat)
require(harmony)
require(ggExtra)

tab10 <- c("26 98 165", "253 105 16", "40 147 34", "203 17 30", "129 78 175", 
           "121 67 58", "218 93 182", "108 108 108", "175 179 28", "30 178 197")
tab10<-sapply(strsplit(tab10, " "), function(x)
  rgb(x[1], x[2], x[3], maxColorValue=255))

FADS_sorted<-readRDS("FADS_sorted.dgecounts.rds")
intron<-FADS_sorted$umicount$intron$all
all<-FADS_sorted$umicount$inex$all
commonBCs<-intersect(colnames(intron),colnames(all))
intron<-intron[,commonBCs]
all<-all[,commonBCs]

total_intron<-Matrix::colSums(intron)
total_all<-Matrix::colSums(all)

NF<-total_intron/total_all

df_sorted<-data.frame(UMIs_detected=total_all,NF=NF,Sample="FADS_obj_plus")

FADS_unsorted<-readRDS("FADS_unsorted.dgecounts.rds")
intron<-FADS_unsorted$umicount$intron$all
all<-FADS_unsorted$umicount$inex$all
commonBCs<-intersect(colnames(intron),colnames(all))
intron<-intron[,commonBCs]
all<-all[,commonBCs]

total_intron<-Matrix::colSums(intron)
total_all<-Matrix::colSums(all)
NF<-total_intron/total_all

df_unsorted<-data.frame(UMIs_detected=total_all,NF=NF,Sample="FADS_obj_min")

df<-rbind(df_sorted,df_unsorted)

df<-subset(df,UMIs_detected>=100)
ggplot(df, aes(x=NF,y=log10(UMIs_detected)))+geom_point(size=0.5)

gbm.nf.umi <- data.frame(row.names = rownames(df),
                         nf=df$NF,
                         umi=df$UMIs_detected)
gbm.ed <- identify_empty_drops(nf_umi=gbm.nf.umi,include_plot = T,nf_rescue = 0.07,umi_rescue = 2000)

ggplot(gbm.ed , aes(x=nf,y=log10(umi),color=cell_status))+geom_point(size=0.5)

FADS_sorted<-FADS_sorted$umicount$inex$all
FADS_unsorted<-FADS_unsorted$umicount$inex$all


gene.meta<-data.frame(fread("~/gene.conversion.txt"))
rownames(gene.meta)<-gene.meta$ID

gene.meta$Gene_Symbol<-make.unique(gene.meta$Gene_Symbol)

id.1<-intersect(rownames(gene.meta),rownames(FADS_sorted))
FADS_sorted<-FADS_sorted[id.1,]
rownames(FADS_sorted)<-gene.meta[id.1,"Gene_Symbol"]

id.2<-intersect(rownames(gene.meta),rownames(FADS_unsorted))
FADS_unsorted<-FADS_unsorted[id.2,]
rownames(FADS_unsorted)<-gene.meta[id.2,"Gene_Symbol"]

FADS_sorted <- CreateSeuratObject(counts = FADS_sorted[,intersect(rownames(df_sorted),rownames(gbm.ed))], 
                                  project = "FADS_sorted",min.cells = 10)
FADS_unsorted <- CreateSeuratObject(counts = FADS_unsorted[,intersect(rownames(df_unsorted),rownames(gbm.ed))], 
                                 project = "FADS_unsorted",min.cells = 10)

FADS_obj <- merge(FADS_sorted, y = FADS_unsorted)

FADS_obj <- NormalizeData(FADS_obj)
FADS_obj <- FindVariableFeatures(FADS_obj, selection.method = "vst", nfeatures = 1000)

top10 <- head(VariableFeatures(FADS_obj), 10)
plot1 <- VariableFeaturePlot(FADS_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)


require(future)
FADS_obj <- ScaleData(FADS_obj, 
                    verbose = T)

print("Running PCA")
FADS_obj <- RunPCA(object      = FADS_obj,  npcs = 50,
                 do.print    = FALSE)

ElbowPlot(FADS_obj,ndims = 50)
nDIM=14

FADS_obj <- RunUMAP(object        = FADS_obj, n.neighbors = 50L,min.dist = 0.0001,
                  reduction = "pca", 
                  dims      = 1:nDIM)

pdf(file = "Fig_100UMIs.pdf",width = 5.25,height=3.5)
DimPlot(FADS_obj, reduction = "umap",cols = tab10,label = F,label.size = 5,repel = T,pt.size = 0.01,shuffle = T)+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("UMAP 1")+ylab("UMAP 2")
dev.off()


p1<-ggplot(FADS_obj@meta.data[sample(1:nrow(FADS_obj@meta.data)),],aes(log10(nCount_RNA),log10(nFeature_RNA),color=orig.ident))+
  geom_point(size=0.1)+scale_color_manual(values = tab10)
ggMarginal(p1, type = "density", groupFill = TRUE)


FADS_obj <- RunHarmony(FADS_obj, group.by.vars = "orig.ident")

FADS_obj <- RunUMAP(object        = FADS_obj, n.neighbors = 50L,min.dist = 0.01, 
                  reduction = "harmony", 
                  dims      = 1:nDIM)

DimPlot(FADS_obj, reduction = "umap",group.by = "orig.ident",pt.size = 0.05,cols = tab10,shuffle = T)


FADS_obj <- FindNeighbors(FADS_obj, dims = 1:nDIM,
                        k.param = 15,
                        reduction = "harmony")
FADS_obj <- FindClusters(FADS_obj, resolution = 0.05)

DimPlot(FADS_obj, reduction = "umap",cols = tab10,label = T,label.size = 10)

meta<-FADS_obj@meta.data

ggplot(meta, aes(log10(nCount_RNA),color=RNA_snn_res.0.05)) + 
  stat_ecdf(geom = "line",show.legend = T,size=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank())+
  scale_color_manual(values=tab10)+
  labs(y="ECDF",x="log10(#Transcripts)")


gbm.ed.flt<-cbind(gbm.ed,cell_type=FADS_obj@meta.data[rownames(gbm.ed),]$RNA_snn_res.0.05)

p1<-ggplot(gbm.ed.flt , aes(x=nf,y=log10(umi),color=cell_type))+
  geom_point(size=1.5)+scale_color_manual(values=tab10)+theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

gbm.ed.dc <- identify_damaged_cells(gbm.ed.flt,nf_sep = 0.07,umi_sep_perc = 20)

m<-gbm.ed.dc$df 
m$cell_status<-factor(m$cell_status,levels=c("damaged_cell","empty_droplet","cell"))

p1<-ggplot(m , aes(x=nf,y=log10(umi),color=cell_status))+theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_manual(labels = c("Damaged cells", "Empty droplets","Cells"),values=tab10[c(3,2,1)])+
  geom_point(size=1.0)+
  coord_trans(y="log10")+
 # scale_y_continuous(breaks = c(1000,10000),labels= c(1000,10000))+
  xlab("Nuclear RNA fraction")+ylab("UMIs detected, log10")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  #theme(legend.position = c(0.8, 0.9))+
  guides(color=guide_legend(title="",override.aes = list(size=2)))
pdf(file = "Fig_DropletQC.pdf",width = 5,height=3.5)
ggMarginal(p1, type = "density", groupFill = F)
dev.off()

