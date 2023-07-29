# spinDrop
This is the repository for the publication: De Jonghe, Kaminski et al. spinDrop: a modular droplet microfluidic platform to maximise single-cell RNA sequencing information content, Nature Communications (2023)

Droplet microfluidic methods have massively increased the throughput of single-cell RNA sequencing campaigns. The benefit of scale-up is, however, accompanied by increased background noise when processing challenging samples as well as lower overall RNA capture efficiency. These drawbacks stem from the lack of strategies to enrich for high-quality material at the moment of cell encapsulation and the lack of implementable multi-step enzymatic processes that increase RNA capture. Here we alleviate both bottlenecks using fluorescence-activated droplet sorting to enrich for droplets that contain single viable cells, intact nuclei or fixed cells and use reagent addition to droplets by picoinjection to perform multi-step lysis and reverse transcription. Our methodology increases gene detection rates fivefold, while reducing background noise by up to half, depending on sample quality. We harness these unique properties to deliver a high-quality molecular atlas of mouse brain development using highly damaged input material. Our method is broadly applicable to other droplet-based workflows to deliver sensitive and accurate single-cell profiling at a reduced cost.

## scRNA-seq pre-processing
In this repository, scripts to pre-process the R1, R2, R3 and R4 are shown. Briefly, the files are pre-processed by zUMIs (Parekh et al., Gigascience (2018)). Next, the concatenated intron and exon matrices are exctracted (inex) and annotations are added using bioMart (Durinck et al. Nature Protocols (2009)).

## LabVIEW FADS
In this repository, one can find the LabVIEW FADS files used for sorting, for more information regarding the scripts, please contact fh111@cam.ac.uk

## scRNA-seq analysis
This repository contains the scripts for data analysis. Data filtering, integration, differential expression analysis and global benchmarking was performed using Seurat (Hao, et al., bioRxiv (2022) [Seurat v5], Hao*, Hao*, et al., Cell (2021) [Seurat v4], Stuart*, Butler*, et al., Cell (2019) [Seurat v3], Butler, et al., Nat Biotechnology (2018) [Seurat v2], Satija*, Farrell*, et al., Nat Biotechnology (2015) [Seurat v1]. RNA velocity analysis, trajectory inference and nascent RNA analyses were carried out using scanpy (Wolf et al. Genome Biology (2018)) and scvelo (Bergen et al. Nature Biotechnology (2020)).
