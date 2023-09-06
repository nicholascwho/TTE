library(Seurat)
library(dplyr)
library(patchwork)
library(harmony)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(GenomicFeatures)
library(ComplexHeatmap)

create_seurat_for_integration <- function(directory, fn_contain, output=c("Seurat","Matrix_DF")){
  # directory: path to directory containing cell ranger output
  # fn_contain: process only samples containing this pattern
  # output: desired output class
  folders <- list.dirs(directory, recursive = F, full.names = F)
  tmp <- 0
  for(i in 1:length(folders)){
    if(any(sapply(fn_contain, function(x){grepl(x,folders[i])}))){
      tmp <- tmp + 1
      final_folder <- paste(directory, "\\", folders[i], "\\filtered_feature_bc_matrix\\", sep = "")
      print(folders[i])
      if(tmp == 1){
        counts_mat <- Matrix::readMM(gzfile(paste(final_folder, "matrix.mtx.gz", sep = "")))
        features <- read.csv(gzfile(paste(final_folder, "features.tsv.gz", sep = "")), header = F, sep = "\t")
        barcodes <- read.csv(gzfile(paste(final_folder, "barcodes.tsv.gz", sep = "")), header = F, sep = "\t")
        barcodes <- apply(barcodes,2,function(x){paste("S",i,"_",x,sep="")})
        meta <- data.frame(origin=rep(folders[i], nrow(barcodes)))
      } else {
        barcodes1 <- read.csv(gzfile(paste(final_folder, "barcodes.tsv.gz", sep = "")), header = F, sep = "\t")
        barcodes1 <- apply(barcodes1,2,function(x){paste("S",i,"_",x,sep="")})
        barcodes <- rbind(barcodes, barcodes1)
        features1 <- read.csv(gzfile(paste(final_folder, "features.tsv.gz", sep = "")), header = F, sep = "\t")
        if(!all(features1[,1]==features[,1])){
          features1$number <- 1:nrow(features1)
          row.names(features1) <- featuress1[,1]
          features1 <- features1[features[,1],]
          counts_mat1 <- Matrix::readMM(gzfile(paste(final_folder, "matrix.mtx.gz", sep = "")))
          counts_mat1 <- counts_mat1[features1$number,]
          counts_mat <- cbind(counts_mat, counts_mat1)
        } else {
          counts_mat <- cbind(counts_mat, Matrix::readMM(gzfile(paste(final_folder, "matrix.mtx.gz", sep = ""))))
        }
        meta <- rbind(meta, data.frame(origin=rep(folders[i], nrow(barcodes1))))
      }
      counts_mat@Dimnames[[2]] <- as.character(barcodes[,1])
      counts_mat@Dimnames[[1]] <- as.character(features[,2])
    }
  }
  print("Read in all datasets")
  row.names(meta) <- counts_mat@Dimnames[[2]]
  dup_gene_names <- duplicated(counts_mat@Dimnames[[1]])
  dup_gene_id <- which(dup_gene_names==T)
  for(i in dup_gene_id){
    dup_gene <- counts_mat@Dimnames[[1]][i]
    first_gene_id <- which(counts_mat@Dimnames[[1]]==dup_gene)[1]
    counts_mat[first_gene_id,] <- counts_mat[first_gene_id,] + counts_mat[i,]
  }
  counts_mat <- counts_mat[-dup_gene_id,]
  print("Removed gene variants")
  if(output=="Seurat"){
    result <- CreateSeuratObject(counts = counts_mat, meta.data = meta, min.cells = 3, min.features = 100)
    result <- subset(result, nCount_RNA > 0)
  } else if(output=="Matrix_DF"){
    result <- list(Counts=counts_mat, Meta=meta)
  }
  return(result)
}
pdac <- create_seurat_for_integration("C:\\Users\\nicho\\Desktop\\NUS_PhD\\230128 PDAC Analysis\\cellranger_out", "Stellate", "Seurat")
pdac[["percent.mt"]] <- PercentageFeatureSet(pdac, pattern = "^MT-")

pdac$PatientID <- strsplit(pdac$origin, split=" ") %>% lapply(function(x){x[[2]]}) %>% unlist()
pdac$Sample_type <- strsplit(pdac$origin, split=" ") %>% lapply(function(x){x[[3]]}) %>% unlist()
pdac$Sample_name <- strsplit(pdac$origin, split=" ") %>% lapply(function(x){x[[1]]}) %>% unlist()

svg("C:/Users/nicho/Desktop/NUS_PhD/230128 PDAC Analysis/Stellate samples QC.svg", width = 12)
VlnPlot(pdac, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Sample_name")
dev.off()
pdac_filt <- subset(pdac, subset = nFeature_RNA > 100 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 10000)

pdac_filt <- SCTransform(pdac_filt, verbose = T, assay = "RNA")
pdac_filt <- RunPCA(pdac_filt, verbose = T, assay = "SCT")
pdac_filt <- RunHarmony(pdac_filt, group.by.vars = c("PatientID"), verbose = T, assay.use = "SCT", max.iter.harmony = 10)
ElbowPlot(pdac_filt, reduction = "harmony")
pdac_filt <- RunUMAP(pdac_filt, reduction = "harmony", verbose = T, assay = "SCT", dims = 1:10)
pdac_filt <- FindNeighbors(pdac_filt, dims = 1:2, reduction = "umap")
pdac_filt <- FindClusters(pdac_filt, res = 0.003)
pdac_filt_umap <- cbind(pdac_filt@meta.data,Embeddings(pdac_filt,reduction="umap"),t(pdac_filt@assays$SCT@data[c("COL1A1","ACTA2","CTHRC1","TAGLN",
                                                                                                                 "IL6","IL11","CXCL1","CXCL8"),]))
pdac_filt_umap$origin <- strsplit(pdac_filt_umap$origin,split=" ") %>% lapply(function(x){x[1]}) %>% unlist()

# to visualize patient sample distribution
pdf(filename, width = 15)  #replace 'filename' with desired filename
ggarrange(ggplot(pdac_filt_umap,aes(UMAP_1,UMAP_2,color=SCT_snn_res.0.003)) + geom_point(size=0.75) + 
            scale_color_manual("Cluster",values=c("#1dbde6","#f1515e")) +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_text(size=10), panel.grid.major = element_blank(),title = element_text(size=14,face="bold"),
                  axis.ticks=element_blank(),axis.text.y = element_text(size=10),
                  axis.title = element_text(size=12),legend.key=element_rect(fill=NA)) +
            xlab("") + ylab("") + ggtitle("Cluster"),
          ggplot(pdac_filt_umap[sample(1:nrow(pdac_filt_umap)),],aes(UMAP_1,UMAP_2,color=origin)) + geom_point(size=0.75) + 
            scale_color_manual("Sample",values=c("#373f51","#355070","#840032","#e59500","#ffbf69","#e5dada")) +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_text(size=10), panel.grid.major = element_blank(),title = element_text(size=14,face="bold"),
                  axis.ticks=element_blank(),axis.text.y = element_text(size=10),
                  axis.title = element_text(size=12),legend.key=element_rect(fill=NA)) +
            xlab("") + ylab("") + ggtitle("Sample"),
          nrow=1,ncol=2)
dev.off()

pdac_filt_genes_2plot <- cbind(pdac_filt@meta.data, t(pdac_filt@assays$SCT@data[c("COL1A1","COL1A2","COL3A1","FN1","THY1",
                                                                                  "RGS5","PDGFRB","MYL9","NDUFA4L2","SOD3"),]))
pdac_filt_genes_2plot$origin <- strsplit(pdac_filt_genes_2plot$origin,split=" ") %>% lapply(function(x){x[1]}) %>% unlist()

#to visualize expression of fibroblast and stellate cell genes by sample
pdf(filename, width = 20, height=10) #replace 'filename' with desired filename
ggarrange(ggplot(pdac_filt_genes_2plot,aes(origin,COL1A1,fill=origin)) + geom_violin() + 
            scale_fill_manual("Sample",values=c("#373f51","#355070","#840032","#e59500","#ffbf69","#e5dada")) +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_blank(), panel.grid.major.x = element_blank(),title = element_text(size=14,face="bold"),
                  panel.grid.major.y = element_line(color="black",linetype="dashed", linewidth=0.5),
                  axis.ticks.x=element_blank(),axis.text.y = element_text(size=11),legend.key=element_rect(fill=NA),
                  axis.title = element_text(size=12)) + geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
            xlab("") + ylab("") + ggtitle("COL1A1"),
          ggplot(pdac_filt_genes_2plot,aes(origin,COL1A2,fill=origin)) + geom_violin() + 
            scale_fill_manual("Sample",values=c("#373f51","#355070","#840032","#e59500","#ffbf69","#e5dada")) +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_blank(), panel.grid.major.x = element_blank(),title = element_text(size=14,face="bold"),
                  panel.grid.major.y = element_line(color="black",linetype="dashed", linewidth=0.5),
                  axis.ticks.x=element_blank(),axis.text.y = element_text(size=11),legend.key=element_rect(fill=NA),
                  axis.title = element_text(size=12)) + geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
            xlab("") + ylab("") + ggtitle("COL1A2"),
          ggplot(pdac_filt_genes_2plot,aes(origin,COL3A1,fill=origin)) + geom_violin() + 
            scale_fill_manual("Sample",values=c("#373f51","#355070","#840032","#e59500","#ffbf69","#e5dada")) +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_blank(), panel.grid.major.x = element_blank(),title = element_text(size=14,face="bold"),
                  panel.grid.major.y = element_line(color="black",linetype="dashed", linewidth=0.5),
                  axis.ticks.x=element_blank(),axis.text.y = element_text(size=11),legend.key=element_rect(fill=NA),
                  axis.title = element_text(size=12)) + geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
            xlab("") + ylab("") + ggtitle("COL3A1"),
          ggplot(pdac_filt_genes_2plot,aes(origin,FN1,fill=origin)) + geom_violin() + 
            scale_fill_manual("Sample",values=c("#373f51","#355070","#840032","#e59500","#ffbf69","#e5dada")) +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_blank(), panel.grid.major.x = element_blank(),title = element_text(size=14,face="bold"),
                  panel.grid.major.y = element_line(color="black",linetype="dashed", linewidth=0.5),
                  axis.ticks.x=element_blank(),axis.text.y = element_text(size=11),legend.key=element_rect(fill=NA),
                  axis.title = element_text(size=12)) + geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
            xlab("") + ylab("") + ggtitle("FN1"),
          ggplot(pdac_filt_genes_2plot,aes(origin,THY1,fill=origin)) + geom_violin() + 
            scale_fill_manual("Sample",values=c("#373f51","#355070","#840032","#e59500","#ffbf69","#e5dada")) +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_blank(), panel.grid.major.x = element_blank(),title = element_text(size=14,face="bold"),
                  panel.grid.major.y = element_line(color="black",linetype="dashed", linewidth=0.5),
                  axis.ticks.x=element_blank(),axis.text.y = element_text(size=11),legend.key=element_rect(fill=NA),
                  axis.title = element_text(size=12)) + geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
            xlab("") + ylab("") + ggtitle("THY1"),
          ggplot(pdac_filt_genes_2plot,aes(origin,RGS5,fill=origin)) + geom_violin() + 
            scale_fill_manual("Sample",values=c("#373f51","#355070","#840032","#e59500","#ffbf69","#e5dada")) +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_blank(), panel.grid.major.x = element_blank(),title = element_text(size=14,face="bold"),
                  panel.grid.major.y = element_line(color="black",linetype="dashed", linewidth=0.5),
                  axis.ticks.x=element_blank(),axis.text.y = element_text(size=11),legend.key=element_rect(fill=NA),
                  axis.title = element_text(size=12)) + geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
            xlab("") + ylab("") + ggtitle("RGS5"),
          ggplot(pdac_filt_genes_2plot,aes(origin,PDGFRB,fill=origin)) + geom_violin() + 
            scale_fill_manual("Sample",values=c("#373f51","#355070","#840032","#e59500","#ffbf69","#e5dada")) +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_blank(), panel.grid.major.x = element_blank(),title = element_text(size=14,face="bold"),
                  panel.grid.major.y = element_line(color="black",linetype="dashed", linewidth=0.5),
                  axis.ticks.x=element_blank(),axis.text.y = element_text(size=11),legend.key=element_rect(fill=NA),
                  axis.title = element_text(size=12)) + geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
            xlab("") + ylab("") + ggtitle("PDGFRB"),
          ggplot(pdac_filt_genes_2plot,aes(origin,MYL9,fill=origin)) + geom_violin() + 
            scale_fill_manual("Sample",values=c("#373f51","#355070","#840032","#e59500","#ffbf69","#e5dada")) +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_blank(), panel.grid.major.x = element_blank(),title = element_text(size=14,face="bold"),
                  panel.grid.major.y = element_line(color="black",linetype="dashed", linewidth=0.5),
                  axis.ticks.x=element_blank(),axis.text.y = element_text(size=11),legend.key=element_rect(fill=NA),
                  axis.title = element_text(size=12)) + geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
            xlab("") + ylab("") + ggtitle("MYL9"),
          ggplot(pdac_filt_genes_2plot,aes(origin,NDUFA4L2,fill=origin)) + geom_violin() + 
            scale_fill_manual("Sample",values=c("#373f51","#355070","#840032","#e59500","#ffbf69","#e5dada")) +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_blank(), panel.grid.major.x = element_blank(),title = element_text(size=14,face="bold"),
                  panel.grid.major.y = element_line(color="black",linetype="dashed", linewidth=0.5),
                  axis.ticks.x=element_blank(),axis.text.y = element_text(size=11),legend.key=element_rect(fill=NA),
                  axis.title = element_text(size=12)) + geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
            xlab("") + ylab("") + ggtitle("NDUFA4L2"),
          ggplot(pdac_filt_genes_2plot,aes(origin,SOD3,fill=origin)) + geom_violin() + 
            scale_fill_manual("Sample",values=c("#373f51","#355070","#840032","#e59500","#ffbf69","#e5dada")) +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_blank(), panel.grid.major.x = element_blank(),title = element_text(size=14,face="bold"),
                  panel.grid.major.y = element_line(color="black",linetype="dashed", linewidth=0.5),
                  axis.ticks.x=element_blank(),axis.text.y = element_text(size=11),legend.key=element_rect(fill=NA),
                  axis.title = element_text(size=12)) + geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
            xlab("") + ylab("") + ggtitle("SOD3"), nrow=2, ncol=5, common.legend = T, legend="right")

dev.off()

#to visualize expression of fibroblast and stellate cell genes on umap
pdf(filename, width = 16, height = 10) #replace 'filename' with desired filename
ggarrange(ggplot(pdac_filt_umap,aes(UMAP_1,UMAP_2,color=COL1A1)) + geom_point(size=0.75) + 
            scale_color_gradient("Cluster",low="#c1c1c1",high="#db2b39") +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_text(size=10), panel.grid.major = element_blank(),title = element_text(size=14,face="bold"),
                  axis.ticks=element_blank(),axis.text.y = element_text(size=10),
                  axis.title = element_text(size=12),legend.key=element_rect(fill=NA)) +
            xlab("") + ylab("") + ggtitle("COL1A1"),
          ggplot(pdac_filt_umap,aes(UMAP_1,UMAP_2,color=CTHRC1)) + geom_point(size=0.75) + 
            scale_color_gradient("Cluster",low="#c1c1c1",high="#db2b39") +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_text(size=10), panel.grid.major = element_blank(),title = element_text(size=14,face="bold"),
                  axis.ticks=element_blank(),axis.text.y = element_text(size=10),
                  axis.title = element_text(size=12),legend.key=element_rect(fill=NA)) +
            xlab("") + ylab("") + ggtitle("CTHRC1"),
          ggplot(pdac_filt_umap,aes(UMAP_1,UMAP_2,color=ACTA2)) + geom_point(size=0.75) + 
            scale_color_gradient("Cluster",low="#c1c1c1",high="#db2b39") +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_text(size=10), panel.grid.major = element_blank(),title = element_text(size=14,face="bold"),
                  axis.ticks=element_blank(),axis.text.y = element_text(size=10),
                  axis.title = element_text(size=12),legend.key=element_rect(fill=NA)) +
            xlab("") + ylab("") + ggtitle("ACTA2"),
          ggplot(pdac_filt_umap,aes(UMAP_1,UMAP_2,color=TAGLN)) + geom_point(size=0.75) + 
            scale_color_gradient("Cluster",low="#c1c1c1",high="#db2b39") +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_text(size=10), panel.grid.major = element_blank(),title = element_text(size=14,face="bold"),
                  axis.ticks=element_blank(),axis.text.y = element_text(size=10),
                  axis.title = element_text(size=12),legend.key=element_rect(fill=NA)) +
            xlab("") + ylab("") + ggtitle("TAGLN"),
          ggplot(pdac_filt_umap,aes(UMAP_1,UMAP_2,color=IL6)) + geom_point(size=0.75) + 
            scale_color_gradient("Cluster",low="#c1c1c1",high="#db2b39") +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_text(size=10), panel.grid.major = element_blank(),title = element_text(size=14,face="bold"),
                  axis.ticks=element_blank(),axis.text.y = element_text(size=10),
                  axis.title = element_text(size=12),legend.key=element_rect(fill=NA)) +
            xlab("") + ylab("") + ggtitle("IL6"),
          ggplot(pdac_filt_umap,aes(UMAP_1,UMAP_2,color=IL11)) + geom_point(size=0.75) + 
            scale_color_gradient("Cluster",low="#c1c1c1",high="#db2b39") +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_text(size=10), panel.grid.major = element_blank(),title = element_text(size=14,face="bold"),
                  axis.ticks=element_blank(),axis.text.y = element_text(size=10),
                  axis.title = element_text(size=12),legend.key=element_rect(fill=NA)) +
            xlab("") + ylab("") + ggtitle("IL11"),
          ggplot(pdac_filt_umap,aes(UMAP_1,UMAP_2,color=CXCL1)) + geom_point(size=0.75) + 
            scale_color_gradient("Cluster",low="#c1c1c1",high="#db2b39") +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_text(size=10), panel.grid.major = element_blank(),title = element_text(size=14,face="bold"),
                  axis.ticks=element_blank(),axis.text.y = element_text(size=10),
                  axis.title = element_text(size=12),legend.key=element_rect(fill=NA)) +
            xlab("") + ylab("") + ggtitle("CXCL1"),
          ggplot(pdac_filt_umap,aes(UMAP_1,UMAP_2,color=CXCL8)) + geom_point(size=0.75) + 
            scale_color_gradient("Cluster",low="#c1c1c1",high="#db2b39") +
            theme(panel.background = element_rect(color="black",linewidth=1,fill=NA),
                  panel.grid.minor = element_blank(), plot.background = element_blank(),
                  axis.text.x=element_text(size=10), panel.grid.major = element_blank(),title = element_text(size=14,face="bold"),
                  axis.ticks=element_blank(),axis.text.y = element_text(size=10),
                  axis.title = element_text(size=12),legend.key=element_rect(fill=NA)) +
            xlab("") + ylab("") + ggtitle("CXCL8"), nrow=2, ncol=4, common.legend = T, legend = "right")
dev.off()

pdac_filt_pseudobulk_column_anno <- HeatmapAnnotation(Cluster = c("0","1"), col = list(Cluster = c("0"="#1dbde6","1"="#f1515e")))

# to visualize expression of CAF cluster DEGs
pdf(filename, width = 7, height = 10) #replace 'filename' with desired filename
Heatmap(pdac_filt_pseudobulk[c(row.names(pdac_filt_deg)[1:30],row.names(pdac_filt_deg)[5585:5614]),], cluster_rows = F,
        col = circlize::colorRamp2(c(4.673033, 7.399862, 10.08504), c("#0d3b66", "#faf0ca", "#db2b39")), name = "Normalized\nCPM", cluster_columns = F,
        rect_gp = gpar(col = "black", lwd = 0.2), row_split = rep(c("Cluster 0 Markers","Cluster 1 Markers"),each=30),
        heatmap_legend_param = list(legend_height = unit(3,"cm")), show_row_names = T, show_column_names = F,
        top_annotation = pdac_filt_pseudobulk_column_anno, row_gap = unit(3,"mm"))
dev.off()

# to visualize fibroblast and stellate genes by sample
pdf(filename, width = 25, height=10) #replace 'filename' with desired filename
VlnPlot(pdac_filt, c("COL1A1","COL1A2","COL3A1","FN1","THY1",
                     "RGS5","PDGFRB","ARIRF","NDUFA4L2","SOD3"), ncol=5, group.by = "origin")
dev.off()

fib_stellate_gs <- read.table(filename, header = T) #replace filename with path to fibroblast and stellate gene signatures

pseudo_bulk <- function(seu, cluster_col){
  # seu: seurat object
  # cluster_col: column name containing cluster information to group cells by
  cell_types <- unique(as.character(seu@meta.data[,cluster_col]))
  cluster_list <- as.list(cell_types)
  cluster_list <- lapply(cluster_list, function(x){rowSums(seu@assays$SCT@counts[,which(seu@meta.data[,cluster_col]==x)])})
  names(cluster_list) <- cell_types
  cluster_list <- data.frame(cluster_list)
  to_rm <- which(rowSums(cluster_list)==0)
  if(length(to_rm)>0){
    cluster_list <- cluster_list[-to_rm,]
    row.names(cluster_list) <- row.names(seu@assays$SCT@counts)[-to_rm]
  } else {
    row.names(cluster_list) <- row.names(seu@assays$SCT@counts)
  }
  cluster_list <- t(asinh(t(cluster_list*1e6)/colSums(cluster_list)))
  # cluster_list <- cluster_list/apply(cluster_list,1,max)
  return(cluster_list)
}
pdac_filt_pseudobulk <- pseudo_bulk(pdac_filt, "SCT_snn_res.0.003")

#get transcript lengths
refseq_ids <- as.list(org.Hs.egREFSEQ[mappedkeys(org.Hs.egREFSEQ)])
mart<- biomaRt::useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
refseq_mapping <- biomaRt::getBM(attributes = c("refseq_mrna","hgnc_symbol","transcript_length"), 
                                 filters="refseq_mrna", # you swap out of this filter for whatever your input is
                                 values=unlist(refseq_ids), # vector of your NMf
                                 mart=mart)
refseq_mapping <- refseq_mapping[which(sapply(refseq_mapping$hgnc_symbol,nchar)>0),]
refseq_mapping <- refseq_mapping[order(refseq_mapping$transcript_length, decreasing=T),] %>% .[-which(duplicated(.$hgnc_symbol)),]

counts_to_tpm <- function(counts, transcript_length_df){
  # counts: counts data
  # transcript_length_df: transcript lengths
  lengths <- sapply(row.names(counts),function(x){if(x %in% transcript_length_df$hgnc_symbol){transcript_length_df$transcript_length[which(transcript_length_df$hgnc_symbol==x)]}else{-1}})
  to_remove <- which(lengths==-1)
  if(length(to_remove)>0){
    lengths <- lengths[-to_remove]
    counts <- counts[-to_remove,]
  }
  lengths <- as.numeric(lengths)/1000
  rpk <- counts/lengths
  lib_size <- colSums(rpk)/1e6
  tpm <- t(t(rpk)/lib_size)
  return(tpm)
}
pdac_stellate_tpm <- counts_to_tpm(pdac_filt_pseudobulk,refseq_mapping)
pdac_stellate_tpm_ssgsea <- cbind(data.frame(Genes=row.names(pdac_stellate_tpm)),pdac_stellate_tpm)
pdac_stellate_tpm_ssgsea <- rbind(c(-666,rep("PDAC",ncol(pdac_stellate_tpm))),pdac_stellate_tpm_ssgsea)
pdac_stellate_tpm_ssgsea <- cbind(data.frame(id=c("Cancer_type",row.names(pdac_stellate_tpm))),pdac_stellate_tpm_ssgsea)
pdac_stellate_tpm1 <- counts_to_tpm(pdac_stellatebyclus_bulk,refseq_mapping)
pdac_stellate_tpm_ssgsea1 <- cbind(data.frame(Genes=row.names(pdac_stellate_tpm1)),pdac_stellate_tpm1)
pdac_stellate_tpm_ssgsea1 <- rbind(c(-666,rep("PDAC",ncol(pdac_stellate_tpm1))),pdac_stellate_tpm_ssgsea1)
pdac_stellate_tpm_ssgsea1 <- cbind(data.frame(id=c("Cancer_type",row.names(pdac_stellate_tpm1))),pdac_stellate_tpm_ssgsea1)
pdac_tumor_tpm <- counts_to_tpm(pdac_tumor_bulk,refseq_mapping)
pdac_tumor_tpm_ssgsea <- cbind(data.frame(Genes=row.names(pdac_tumor_tpm)),pdac_tumor_tpm)
pdac_tumor_tpm_ssgsea <- rbind(c(-666,rep("PDAC",ncol(pdac_tumor_tpm))),pdac_tumor_tpm_ssgsea)
pdac_tumor_tpm_ssgsea <- cbind(data.frame(id=c("Cancer_type",row.names(pdac_tumor_tpm))),pdac_tumor_tpm_ssgsea)
write.table("#1.3", 
            file=filename, # replace filename with desired path
            col.names = F, sep = "\t", quote=F, row.names=F, append=T, eol="\n")
write.table("15809\t6\t1\t1", 
            file=filename, # replace filename with desired path
            col.names = F, sep = "\t", quote=F, row.names=F, append=T, eol="\n")
write.table(pdac_stellate_tpm_ssgsea, 
            file=filename, # replace filename with desired path
            col.names = T, sep = "\t", quote=F, row.names=F, append=T)

ssgsea <- read.table("path to ssGSEA NES results", sep="\t", header=T, row.names=1, skip=2)
get_scores_df <- function(ssgsea_nes){
  res <- data.frame(set=rep(row.names(ssgsea_nes),ncol(ssgsea_nes)),
                    sample=rep(colnames(ssgsea_nes),each=nrow(ssgsea_nes)),
                    NES=unlist(ssgsea_nes))
  return(res)
}
ssgsea_nes_df <- get_scores_df(ssgsea[2:13,16:21])
ssgsea_nes_df$NES <- as.numeric(ssgsea_nes_df$NES)

# to visualize ssGSEA results
svg(filename, # replace filename with desired path
    width=9, height=6)
ggplot(ssgsea_clus_nes_df, aes(x=cell_type,y=as.numeric(NES),fill=cell_type)) +
  geom_violin(color="black", position=position_dodge(0.5), linewidth=0.5) +
  scale_fill_manual(values=c("#006e90","#f18f01")) +
  theme(panel.background = element_rect(color="black",linewidth=2,fill=NA),
        panel.grid.minor = element_blank(), plot.background = element_blank(),
        axis.text.x=element_text(size=12,color="black",angle=90,hjust=1), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="black",linetype="dashed", linewidth=0.5),
        axis.ticks.x=element_blank(), legend.position="none",
        axis.title = element_text(size=16), strip.background = element_blank(), strip.text = element_text(size=14)) +
  ylab("NES") + xlab("Cell Type") + geom_boxplot(width=0.1, fill="#e3e3e3", linewidth=0.5) + ylim(c(0,20)) +
  facet_grid(cols = vars(sample))
dev.off()
