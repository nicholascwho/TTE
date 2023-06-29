suppressPackageStartupMessages({
  library(openxlsx)
  library(DESeq2)
  library(ggplot2)
  library(ComplexHeatmap)
  library(dplyr)
  library(FactoMineR)
  library(factoextra)
  library(apeglm)
  library(patchwork)
  library(sva)
  library(Seurat)
  library(KEGGREST)
  library(msigdbr)
  library(clusterProfiler)
})

# function used to combine raw counts files from different samples
create_counts_df <- function(dir, start_row, counts_col, counts_file_id, sep = "_"){
  filenames <- list.files(dir)
  x <- 0
  for(i in 1:length(filenames)){
    if(grepl(counts_file_id, filenames[i])){
      x <- x + 1
      file <- read.table(paste(dir, filenames[i], sep = "\\")) %>% .[start_row:nrow(.),]
      genes <- file[,1]
      file <- as.data.frame(matrix(as.numeric(file[,counts_col]), ncol = 1))
      row.names(file) <- genes
      if(x == 1){
        result <- file
      } else {
        res_not_file <- setdiff(rownames(result), rownames(file))
        file_not_res <- setdiff(rownames(file), rownames(result))
        if(length(res_not_file) > 0){
          tmp <- matrix(0, nrow = length(res_not_file), ncol = 1)
          row.names(tmp) <- res_not_file
          file <- rbind(file, tmp) %>% .[order(row.names(.)),]
        }
        if(length(file_not_res) > 0){
          tmp <- matrix(0, nrow = length(file_not_res), ncol = result)
          row.names(tmp) <- file_not_res
          result <- rbind(result, tmp) %>% .[order(row.names(.)),]
        }
        result <- cbind(result, file)
      }
    }
  }
  colnames(result) <- strsplit(filenames, split = sep) %>% lapply(function(x){x[1]}) %>% unlist()
  no_exp <- which(as.numeric(rowSums(result))==0)
  if(length(no_exp)>0){
    result <- result[-which(rowSums(result)==0),]
  }
  return(result)
}

p12_14_16 <- create_counts_df("# insert folder containing STAR alignment output", 5, 4, "ReadsPerGene")
p12_14_16_meta <- data.frame()# create data frame containing 3 columns named Patient, Condition, and Batch
p12_14_16_cpm <- t(t(p12_14_16*10e6)/colSums(p12_14_16))
hist(log(rowMeans(p12_14_16_cpm)+1,2),breaks = 100)
p12_14_16_cpm_filt <- p12_14_16_cpm[which(rowMeans(p12_14_16_cpm)>31),]
p12_14_16_filt <- p12_14_16[which(rowMeans(p12_14_16_cpm)>31),]

p12_14_16_filt_adj <- ComBat_seq(as.matrix(p12_14_16_filt), batch = p12_14_16_meta$Batch)
p12_14_16_filt_adj2 <- ComBat_seq(as.matrix(p12_14_16_filt_adj), batch = p12_14_16_meta$Patient)
p12_14_16_filt_adj2_cpm <- t(t(p12_14_16_filt_adj2*10e6)/colSums(p12_14_16_filt_adj2))

# function to calculate average expression within same condition
average_expression <- function(dat, met){
  lst_res <- list()
  for(i in levels(met)){
    if(length(which(met==i))==1){
      lst_res[[i]] <- dat[,which(met==i)]
    } else {
      # lst_res[[i]] <- apply(dat[,which(met==i)],1,median)
      lst_res[[i]] <- rowMeans(dat[,which(met==i)])
    }
  }
  res <- data.frame(lst_res)
  row.names(res) <- row.names(dat)
  return(res)
}
p12_14_16_filt_adj_cpm_mean <- average_expression(p12_14_16_filt_adj2_cpm, factor(p12_14_16_meta$Condition,levels=c("D0 Pre-Slicing","D0 Post-Slicing","RGD","RDG")))

paths_ha_row <- rowAnnotation(df = data.frame("Pathway"=factor(c(rep("Notch",7),rep("JAK/STAT",20),rep("TGFβ",20)),levels=c("JAK/STAT","Notch","TGFβ"))),
                              col = list(Pathway=c("JAK/STAT"="#ff218c","Notch"="#ffd800",
                                                   "TGFβ"="#21b1ff")))
# used to plot heatmap in Figure 1F
Heatmap(t(scale(t(asinh(p12_14_16_filt_adj_cpm_mean[c(rgd_kegg_gsea2@result["KEGG_NOTCH_SIGNALING_PATHWAY","core_enrichment"] %>% strsplit("/") %>% unlist(),
                                                      rgd_kegg_gsea2@result["KEGG_JAK_STAT_SIGNALING_PATHWAY","core_enrichment"] %>% strsplit("/") %>% unlist(),
                                                      rgd_kegg_gsea2@result["KEGG_TGF_BETA_SIGNALING_PATHWAY","core_enrichment"] %>% strsplit("/") %>% unlist()),c(1:4)])))),
        cluster_rows = T, row_gap = unit(3,"mm"), right_annotation = paths_ha_row, row_split = c(rep("Notch",7),rep("JAK/STAT",20),rep("TGFβ",20)), 
        col = circlize::colorRamp2(c(-2, 0, 2), c("#023E8A", "#fdf0d5", "#BA181B")), name = "Relative\nExpression",
        rect_gp = gpar(col = "black", lwd = 0.2), heatmap_legend_param = list(legend_height = unit(3,"cm")), column_names_rot = 90,
        cluster_columns = T, column_names_gp = gpar(fontsize = 15), row_names_gp = gpar(fontface = "italic",cex=1),
        column_labels = c("D0 Pre-slicing", "D0 Post-slicing", "RGD", "RDG"))
