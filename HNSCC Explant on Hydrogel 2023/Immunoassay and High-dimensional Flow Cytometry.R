# for immunoassay
suppressPackageStartupMessages({library(dendextend)})
legendplex <- read.xlsx("# filename of data containing immunoassay results", colNames = T, rowNames = T)
row_dend <- as.dendrogram(hclust(dist(scale(asinh(legendplex)))))
row_dend <- color_branches(row_dend, k = 2, col = c("#BA181B", "#023E8A"))
Heatmap(t(scale(asinh(legendplex))), name = "z-score", cluster_columns = row_dend, column_split = 2,
        col = colorRamp2(c(-2, 0, 2), c("#023E8A", "#fdf0d5", "#BA181B")), row_names_side = "left", row_dend_side = "right",
        rect_gp = gpar(col = "black", lwd = 0.05), column_dend_gp = gpar(lwd=1), row_dend_gp = gpar(lwd=0.5, col="black"))

#for high-dimensional flow cytometry
suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyverse)
  library(data.table)
  library(circlize)
  library(SingleCellExperiment)
  library(Rphenograph)
  library(flowAI)
  library(ComplexHeatmap)
})
set.seed(123)

daf_pheno <- Rphenograph(t(assay(daf)),k=30)
clusters <- rep(0,length(daf$cluster_PhenoGraph_30))
for(i in 1:21){clusters[as.numeric(daf_pheno[[2]][[i]])] <- i}
lymphoid <- as.data.frame(cbind(reducedDim(lymphoid_sce, "UMAP", withDimnames=TRUE),data.frame(cluster=as.character(lymphoid_sce$cluster_PhenoGraph_30),patient=lymphoid_sce$patient_id,condition=lymphoid_sce$condition)))
# lymphoid <- as.data.frame(cbind(reducedDim(daf, "UMAP", withDimnames=TRUE),clusters))
colnames(lymphoid)[1:2] <- c("UMAP1","UMAP2")
lymphoid_ct <- c("NK cell","Treg","NKT-like cell","CD4 T cell","NK cell","Vd2+ T cell","Myeloid","B cell","CD8 T cell","CD69+ CD103+ CD8 T cell","Myeloid","LAG3+ CD4 T cell","Myeloid","CCR5+ CD69+ B cell","Plasma cell","Myeloid","CD39+ CD103+ CD8 T cell","CD4 T cell","CD4 T cell","Unknown HN396 cell","Myeloid")
imm_col <- c("CD8 T cell"="#DCD316","CD69+ CD103+ CD8 T cell"="#ff8d00","CD39+ CD103+ CD8 T cell"="#ffee00","CD4 T cell"="#55cdfc","LAG3+ CD4 T cell"="#008121","Treg"="#f7a8b8","NKT-like cell"="#000000","Vd2+ T cell"="#794e10",
             "NK cell"="#e50000","B cell"="#760188","CCR5+ CD69+ B cell"="#004cff","Plasma cell"="#9661ff","Myeloid"="#9ef01a","Unknown HN396 cell"="#2f6690")
lymphoid$cluster <- factor(sapply(lymphoid$cluster, function(x){lymphoid_ct[x]}), levels=c(names(imm_col)))
ggplot(lymphoid, aes(UMAP1,UMAP2,color=cluster)) + geom_point(size=0.4, show.legend = F) +
  theme(panel.background=element_rect(color="black",fill=NA),axis.title=element_blank(),axis.text=element_blank(),
        axis.line=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank())

#function to determine proportion of each cell type by patient
flow_proportion <- function(meta, group_by, split_by, combine_by, group_levels, combine_levels){
  g <- unique(as.character(meta[,group_by]))
  cb <- unique(as.character(meta[,combine_by]))
  s <- unique(as.character(meta[,split_by]))
  result <- list()
  for(ss in s){
    result[[ss]] <- list(group=c(),split=c(),proportion=c())
    names(result[[ss]]) <- c(group_by, combine_by, "Proportion")
    meta_s <- meta[which(meta[,split_by]==ss),]
    for(cc in cb){
      meta_sc <- meta_s[which(meta_s[,combine_by]==cc),]
      ncell_sc <- nrow(meta_sc)
      for(gg in g){
        result[[ss]][[group_by]] <- c(result[[ss]][[group_by]], gg)
        result[[ss]][[combine_by]] <- c(result[[ss]][[combine_by]], cc)
        result[[ss]]$Proportion <- c(result[[ss]]$Proportion, length(which(meta_sc[,group_by]==gg))*100/ncell_sc)
      }
    }
    result[[ss]] <- data.frame(result[[ss]])
    result[[ss]][,group_by] <- factor(result[[ss]][,group_by], levels = group_levels)
    result[[ss]][,combine_by] <- factor(result[[ss]][,combine_by], levels = combine_levels)
  }
  result1 <- result[[1]]
  if(length(result)>1){
    for(i in 2:length(result)){
      result1 <- rbind(result1, result[[i]])
    }
  }
  result1[,split_by] <- rep(s, each=length(g)*length(cb))
  return(result1)
}
lymphoid_prop <- flow_proportion(lymphoid, "cluster", "patient", "condition", names(imm_col), c("D0", "D7 HA"))
ggplot(lymphoid_prop, aes(x=condition, y=Proportion, fill=cluster)) +
  geom_bar(colour = "black", stat = "identity",show.legend = F) + ylab("Proportion (%)") + xlab("") + #coord_polar("y",start=0) +
  scale_fill_manual(values = imm_col, name = "") + theme_void() +
  theme(axis.text = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(),
        axis.title = element_blank(), strip.text.y.left = element_blank(),
        strip.background = element_blank(), legend.text = element_blank(), legend.title = element_blank(),
        legend.key.height = unit(0,"cm")) + 
  facet_grid(.~patient, switch="y")

# function to determine differentially expressed markers/proteins. Outputs used to generate barplots representing cell type proportions in Figures 2 and 4
diff_exp_flow_markers <- function(sce, clus){
  clus_label <- unique(as.numeric(clus)) %>% .[order(.)]
  dat <- assay(sce)
  result <- list(Ref_mean=c(),Rest_mean=c(),Mean_diff=c(),Ref_median=c(),Rest_median=c(),Median_diff=c(),Cluster=rep(clus_label,each=nrow(dat)),Marker=rep(row.names(dat),length(clus_label)))
  for(i in clus_label){
    ref_id <- which(clus==i)
    result$Ref_mean <- c(result$Ref_mean, rowMeans(dat[,ref_id]))
    result$Rest_mean <- c(result$Rest_mean, rowMeans(dat[,-ref_id]))
    result$Ref_median <- c(result$Ref_median, rowMedians(dat[,ref_id]))
    result$Rest_median <- c(result$Rest_median, rowMedians(dat[,-ref_id]))
  }
  result$Mean_diff <- result$Ref_mean - result$Rest_mean
  result$Median_diff <- result$Ref_median - result$Rest_median
  result <- data.frame(result) %>% .[order(.$Median_diff, decreasing = T),] %>% .[order(.$Cluster),]
  return(result)
}
fibroblast_dem <- diff_exp_flow_markers(fibroblast_sce, as.numeric(membership(fibroblast_pheno[[2]])))
lymphoid_dem <- diff_exp_flow_markers(lymphoid_sce, lymphoid_sce@colData$cluster_PhenoGraph_30)
# above used to generate Supplementary Item 4
