suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
  library(ComplexHeatmap)
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(clusterProfiler)
  library(harmony)
  library(circlize)
  library(msigdbr)
  library(org.Hs.eg.db)
})

# function to read in and combine data from multiple samples
create_seurat_for_integration <- function(directory, fn_contain, output=c("Seurat","Matrix_DF")){
  folders <- list.dirs(directory, recursive = F, full.names = F)
  tmp <- 0
  for(i in 1:length(folders)){
    if(any(sapply(fn_contain, function(x){grepl(x,folders[i])}))){
      tmp <- tmp + 1
      final_folder <- paste(directory, "\\", folders[i], "\\outs\\filtered_feature_bc_matrix\\", sep = "")
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
all_slices <- create_seurat_for_integration("#insert folder name containing 10x scRNAseq output", " ", output="Seurat")
all_slices[["percent_mt"]] <- PercentageFeatureSet(all_slices, pattern="^MT-")
all_slices@meta.data[,"PatientID"] <- strsplit(all_slices$origin, split = " ") %>% lapply(function(x){x[1]}) %>% unlist()
all_slices@meta.data[,"Day"] <- strsplit(all_slices$origin, split = " ") %>% lapply(function(x){x[2]}) %>% unlist() %>% factor(levels = c("D0","D7"))
all_slices@meta.data[,"Condition"] <- strsplit(all_slices$origin, split = " ") %>% lapply(function(x){x[3]}) %>% unlist() %>% factor(levels = c("Original","PhaseX"))
all_slices <- subset(all_slices, nCount_RNA<50000 & nFeature_RNA<7000 & percent_mt<15)

#Extended Figure 2A
VlnPlot(all_slices, features=c("nCount_RNA","nFeature_RNA","percent_mt"), group.by = "origin")

all_slices <- SCTransform(all_slices)
all_slices <- RunPCA(all_slices)
all_slices <- RunHarmony(all_slices, group.by.vars = "origin", verbose = T, assay.use = "SCT", max.iter.harmony = 15, tau = 300)
ElbowPlot(all_slices, reduction = "harmony")
all_slices <- RunUMAP(all_slices, reduction = "harmony", verbose = F, assay = "SCT", dims = 1:15)
all_slices <- FindNeighbors(all_slices, dims = 1:2, reduction = "umap")

cell_type_gs <- msigdbr(category = "C8") %>% .[, c("gs_name", "human_gene_symbol")]
go_gs <- msigdbr(category = "C5", subcategory = "BP") %>% .[, c("gs_name", "human_gene_symbol")]
kegg_gs <- msigdbr(category = "C2", subcategory = "KEGG") %>% .[, c("gs_name", "human_gene_symbol")]
hallmark_gs <- msigdbr(category = "H") %>% .[, c("gs_name", "human_gene_symbol")]
reactome_gs <- msigdbr(category = "C2", subcategory = "REACTOME") %>% .[, c("gs_name", "human_gene_symbol")]

all_slices_filt <- subset(all_slices, SCT_snn_res.0.01 %in% c(0:8))
all_slices_global_deg <- FindAllMarkers(all_slices_filt, logfc.threshold = 0, return.thresh = 1) %>% .[order(.$avg_log2FC, decreasing = T),] %>% .[order(.$cluster),]
# above used to produce Supplementary Item 2 (sheet name: Global), which contains a filtered list where average log2 fold change is >0.25 and p value <0.01
all_slices_cell_types <- lapply(as.list(c(0:8)),function(x){GSEA(setNames(all_slices_global_deg$avg_log2FC[which(all_slices_global_deg$cluster==x)],nm=all_slices_global_deg$gene[which(all_slices_global_deg$cluster==x)]),TERM2GENE = cell_type_gs, pvalueCutoff = 1, eps = 0, seed = 123, pAdjustMethod = "fdr")})
all_slices_filt$cell_type <- sapply(as.numeric(Idents(all_slices_filt)), function(x){c("Lymphoid","Epithelial","Mesenchymal","Myeloid","Lymphoid","Lymphoid","Myocyte","Endothelial","Myeloid")[x]}) %>% factor(levels = unique(.)[order(.)])

all_slices_sub <- list("Lymphoid"=subset(all_slices_filt, cell_type == "Lymphoid"),
                       "Myeloid"=subset(all_slices_filt, cell_type == "Myeloid"),
                       "Mesenchymal"=subset(all_slices_filt, cell_type == "Mesenchymal"),
                       "Epithelial"=subset(all_slices_filt, cell_type == "Epithelial"))
all_slices_sub <- lapply(all_slices_sub, function(x){FindNeighbors(x, dims = 1:2, reduction = "umap")})
all_slices_sub$Lymphoid <- FindClusters(all_slices_sub$Lymphoid, resolution = 0.15)
all_slices_sub$Myeloid <- FindClusters(all_slices_sub$Myeloid, resolution = 0.1)
all_slices_sub$Mesenchymal <- FindClusters(all_slices_sub$Mesenchymal, resolution = 0.02)
all_slices_sub$Epithelial <- FindClusters(all_slices_sub$Epithelial, resolution = 0.07)

all_slices_sub_deg_nonaming <- lapply(all_slices_sub, function(x){FindAllMarkers(x, return.thresh = 0.01, logfc.threshold = 0.1, only.pos = T) %>%  .[order(.$avg_log2FC, decreasing = T),] %>% .[order(.$cluster),]})
all_slices_sub$Lymphoid$sub_cell_type <- sapply(as.numeric(Idents(all_slices_sub$Lymphoid)), function(x){c("Treg","B cell","Naive/Memory T cell","CD8 T cell","CD8 T cell","Plasma cell","B cell","B cell","CD8 T cell","NKT cell","CD4 T cell","Plasma cell","Myocyte","Plasma cell","B cell","Myocyte","Plasma cell")[x]}) %>% factor(levels = unique(.)[order(.)])
all_slices_sub$Myeloid$sub_cell_type <- sapply(as.numeric(Idents(all_slices_sub$Myeloid)), function(x){c("IL1B high Macrophage","Neutrophil","C1Q high Macrophage","Monocyte","Mast cell","cDC","Neutrophil","Mast cell","Neutrophil")[x]}) %>% factor(levels = unique(.)[order(.)])
all_slices_sub$Fibroblast$sub_cell_type <- sapply(as.numeric(Idents(all_slices_sub$Fibroblast)), function(x){c("myCAF","Pericyte","iCAF","apCAF")[x]}) %>% factor(levels = unique(.)[order(.)])
all_slices_sub$Epithelial$sub_cell_type <- sapply(as.numeric(Idents(all_slices_sub$Epithelial)), function(x){c("KRT high Epithelial","pEMT Epithelial","KRT high Epithelial","Cycling Epithelial","Stress Response Epithelial","Inteferon Response Epithelial","KRT high Epithelial","cEMT Epithelial","KRT high Epithelial","Stress Response Epithelial","Stem Cell-like Epithelial")[x]}) %>% factor(levels = unique(.)[order(.)])
# for epi: https://www.nature.com/articles/s41588-022-01141-9
all_slices_sub$Myeloid$sub_cell_type1 <- sapply(as.numeric(Idents(all_slices_sub$Myeloid)), function(x){c("Macrophage","Neutrophil","Macrophage","Monocyte","Mast cell","cDC","Neutrophil","Mast cell","Neutrophil")[x]}) %>% factor(levels = unique(.)[order(.)])
all_slices_sub$Fibroblast$sub_cell_type1 <- sapply(as.numeric(Idents(all_slices_sub$Fibroblast)), function(x){c("CAF","Pericyte","CAF","CAF")[x]}) %>% factor(levels = unique(.)[order(.)])
all_slices_sub$Lymphoid$sub_cell_type1 <- sapply(as.numeric(Idents(all_slices_sub$Lymphoid)), function(x){c("T cell","B cell","T cell","T cell","T cell","Plasma cell","B cell","B cell","T cell","T cell","T cell","Plasma cell","Myocyte","Plasma cell","B cell","Myocyte","Plasma cell")[x]}) %>% factor(levels = unique(.)[order(.)])
for(i in 1:4){Idents(all_slices_sub[[i]]) <- all_slices_sub[[i]]$sub_cell_type}
all_slices_sub_deg <- lapply(all_slices_sub, function(x){FindAllMarkers(x, return.thresh = 0.01, logfc.threshold = 0.1) %>%  .[order(.$avg_log2FC, decreasing = T),] %>% .[order(.$cluster),]})

# function to update original Seurat object with new cell subtype information
rename_ref_seu <- function(ref_seu, sub_seu_lst, new_col_name){
  final_ct <- as.character(Idents(ref_seu))
  names(final_ct) <- names(Idents(ref_seu))
  lvls <- levels(Idents(ref_seu))
  for(i in 1:length(sub_seu_lst)){
    final_ct[row.names(sub_seu_lst[[i]]$seu@meta.data)] <- as.character(sub_seu_lst[[i]]$seu@meta.data[,sub_seu_lst[[i]]$col])
  }
  final_ct <- factor(final_ct, levels=unique(final_ct))
  ref_seu[[new_col_name]] <- final_ct
  Idents(ref_seu) <- ref_seu[[new_col_name]]
  return(ref_seu)
}

Idents(all_slices_filt) <- all_slices_filt$cell_type
all_slices_sub$`Myeloid (non-Mast cell)`$sub_cell_type1 <- as.character(all_slices_sub$`Myeloid (non-Mast cell)`$sub_cell_type) %>% sapply(function(x){if(grepl("Macrophage",x)){"Macrophage"}else{x}}) %>% factor(., levels=c("Macrophage","Monocyte","Dendritic cell","Neutrophil","Mast cell"))
all_slices_filt <- rename_ref_seu(all_slices_filt,
                                  list(list(seu=all_slices_sub$Fibroblast, col="sub_cell_type", type="Fibroblast"),
                                       list(seu=all_slices_sub$`T cell`, col="sub_cell_type", type="T cell"),
                                       list(seu=all_slices_sub$`Myeloid (non-Mast cell)`, col="sub_cell_type", type="Myeloid (non-Mast cell)"),
                                       list(seu=all_slices_sub$`Cancer cell`, col="sub_cell_type", type="Cancer cell")), "final_cell_type")
all_slices_filt$final_cell_type1 <- as.character(all_slices_filt$final_cell_type) %>% sapply(function(x){
  if(grepl("Macrophage",x)){"Macrophage"}
  else if(grepl("Cancer",x)){"Cancer cell"}
  else if(grepl("CAF",x)){"Fibroblast"}
  else{x}
} %>% factor(levels=c("Cancer cell","Fibroblast","Pericyte","Myocyte","Endothelial","CD8 T cell","CD4 T cell","Treg","Naive/Memory T cell","B cell","Plasma cell","Macrophage","Dendritic cell","Monocyte","Neutrophil","Mast cell")))
all_slices_filt$final_cell_type2 <- as.character(all_slices_filt$final_cell_type1) %>% sapply(function(x){
  if(grepl("T",x)){"T cell"}
  else{x}
} %>% factor(levels=c("Cancer cell","Fibroblast","Pericyte","Myocyte","Endothelial","T cell","B cell","Plasma cell","Macrophage","Dendritic cell","Monocyte","Neutrophil","Mast cell")))
all_slices$final_cell_type2 <- rep("Other",nrow(all_slices@meta.data))
all_slices$final_cell_type2[row.names(all_slices_filt@meta.data)] <- as.character(all_slices_filt$final_cell_type2)
all_slices$final_cell_type2 <- factor(all_slices$final_cell_type2, levels=c("Cancer cell","Fibroblast","Pericyte","Myocyte","Endothelial","T cell","B cell","Plasma cell","Macrophage","Dendritic cell","Monocyte","Neutrophil","Mast cell","Other"))

Idents(all_slices) <- all_slices$final_cell_type2
# to generate UMAP plots across Figures 2 to 4 and Extended Figures 2 to 3. Vary same code for respective plots
ggplot(cbind(Embeddings(subset(all_slices_filt,final_cell_type %in% levels(all_slices_filt$final_cell_type)[c(3,4,6,7,10:13,15,17,18,21,22)]&Condition=="Original"),"umap"),
             subset(all_slices_filt,final_cell_type %in% levels(all_slices_filt$final_cell_type)[c(3,4,6,7,10:13,15,17,18,21,22)]&Condition=="Original")[["final_cell_type"]]), aes(UMAP_1,UMAP_2,color=final_cell_type)) +
  geom_point(size=0.4) + scale_color_manual(values=c("CD8 T cell"="#acff35","CD4 T cell"="#55cdfc",
                                                     "Treg"="#f7a8b8","Naive/Memory T cell"="#000000",
                                                     "SPP1 high Macrophage"="#794e10","C1Q high Macrophage"="#e50000",
                                                     "APOE high Macrophage"="#ff8d00","B cell"="#ffee00",
                                                     "Plasma cell"="#008121","Neutrophil"="#004cff",
                                                     "cDC"="#760188","Monocyte"="#9661ff",
                                                     "Mast cell"="#e0e1dd")) +
  theme(panel.background = element_rect(fill=NA,color="black"), axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        legend.position = "none") + xlim(-15,7) + ylim(-15,7)

all_slices_filt$imm_ct <- sapply(as.character(all_slices_filt$final_cell_type2), function(x){if(x %in% c("Macrophage","T cell","Mast cell","Neutrophil","B cell","Plasma cell","Dendritic cell","Monocyte")){"Immune"}else{"Non_immune"}})
# function to create data frame suitable for ggplot
prop_barplot <- function(df, cond, ct, patient){
  conds <- unique(as.character(df[,cond]))
  no_cond <- length(conds)
  cts <- unique(as.character(df[,ct]))
  no_ct <- length(cts)
  patients <- unique(as.character(df[,patient]))
  no_patient <- length(patients)
  result <- list(Condition=rep(conds, each = no_patient*no_ct),
                 Patient=rep(rep(patients, each = no_ct),no_cond),
                 Cell_type=rep(rep(cts,no_patient),no_cond),
                 Proportion=c())
  for(i in conds){
    for(j in patients){
      total_cells <- length(which(df[,patient]==j & df[,cond]==i))
      for(k in cts){
        curr_cell_type <- length(which(df[,patient]==j & df[,cond]==i & df[,ct]==k))
        result$Proportion <- c(result$Proportion, curr_cell_type*100/total_cells)
      }
    }
  }
  return(data.frame(result) %>% .[order(.$Cell_type),] %>% .[order(.$Patient),] %>% .[order(.$Cell_type),])
}
all_slices_prop <- prop_barplot(all_slices_filt@meta.data,"Condition","cell_type","PatientID")

all_slices_prop_nonimmune <- prop_barplot(subset(all_slices_filt,imm_ct=="Non_immune")@meta.data,"Condition","final_cell_type2","PatientID")
all_slices_prop_immune <- prop_barplot(subset(all_slices_filt,imm_ct=="Immune")@meta.data,"Condition","final_cell_type2","PatientID")
all_slices_prop_immune$Condition <- factor(all_slices_prop_immune$Condition, levels=c("Original","PhaseX"))
all_slices_prop_immune$Patient <- factor(all_slices_prop_immune$Patient, levels=c("HN377","HN385","HN390","HN396"))
all_slices_prop_immune$Cell_type <- factor(all_slices_prop_immune$Cell_type, levels(all_slices_filt$final_cell_type2)[6:13])
all_slices_prop_nonimmune$Condition <- factor(all_slices_prop_nonimmune$Condition, levels=c("Original","PhaseX"))
all_slices_prop_nonimmune$Patient <- factor(all_slices_prop_nonimmune$Patient, levels=c("HN377","HN385","HN390","HN396"))
all_slices_prop_nonimmune$Cell_type <- factor(all_slices_prop_nonimmune$Cell_type, levels=levels(all_slices_filt$final_cell_type2)[1:5])
all_slices_prop_epi <- prop_barplot(all_slices_sub$`Cancer cell`@meta.data,"Condition","sub_cell_type","PatientID")
all_slices_prop_caf <- prop_barplot(all_slices_sub$Fibroblast@meta.data,"Condition","sub_cell_type","PatientID")
all_slices_prop_immune_sub <- prop_barplot(subset(all_slices_filt,imm_ct == "Immune")@meta.data,"Condition","final_cell_type","PatientID")

# used to produce proportion barplots in Figures 2 to 4 and Extended Figure 3. Vary same code to generate respective plots.
ggplot(all_slices_prop_nonimmune, aes(x=Condition, y=Proportion, fill=Cell_type)) +
  geom_bar(colour = "black", stat = "identity", show.legend = F) + ylab("Proportion (%)") + xlab("") + #coord_polar("y",start=0) +
  scale_fill_manual(values = c("Cancer cell"="#FF0000","Fibroblast"="#FF8000","Pericyte"="#16B251",
                               "Myocyte"="#0066FF","Endothelial"="#F74ED6")) +
  theme(axis.text = element_text(size=10), axis.ticks = element_line(linewidth=1),
        panel.background = element_rect(fill=NA,color="black",linewidth=1),
        panel.grid = element_blank(),
        axis.title = element_blank(), strip.text.y.left = element_blank(), strip.text = element_blank(),
        strip.background = element_blank(), legend.text = element_blank(), legend.title = element_blank()) +
  facet_wrap(.~Patient, nrow=1)

#run infercnv
tmp <- all_slices_filt@assays$SCT@counts
colnames(tmp) <- sapply(colnames(all_slices_filt@assays$SCT@counts),function(x){gsub("-",".",x)})
all_slices_filt$final_cell_type2a <- as.character(all_slices_filt$final_cell_type2)
all_slices_filt$final_cell_type2a[which(all_slices_filt$final_cell_type2a=="Cancer cell")] <- paste("malignant",all_slices_filt$PatientID[which(all_slices_filt$final_cell_type2a=="Cancer cell")],sep="_")
write.table(data.frame(cell=sapply(colnames(all_slices_filt@assays$SCT@counts),function(x){gsub("-",".",x)}),
                       type=all_slices_filt$final_cell_type2a),
            file="C:/Users/nicho/Desktop/NUS_PhD/220204 scRNAseq_slices/all_slices_meta.txt", sep = "\t", quote = F, row.names = F, col.names = F)
all_slices_infercnv <- infercnv::CreateInfercnvObject(raw_counts_matrix=tmp,
                                                      annotations_file="C:/Users/nicho/Desktop/NUS_PhD/220204 scRNAseq_slices/all_slices_meta.txt",
                                                      delim="\t",
                                                      gene_order_file="C:/Users/nicho/Desktop/NUS_PhD/220204 scRNAseq_slices/Infercnv genomic positions.txt",
                                                      ref_group_names=NULL) #this line removes assumption of malignancy
options(scipen = 100)
all_slices_infercnv <- infercnv::run(all_slices_infercnv,
                                     cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                     out_dir="C:/Users/nicho/Desktop/NUS_PhD/220204 scRNAseq_slices/infercnv_230410",
                                     cluster_by_groups=TRUE,
                                     denoise=TRUE,
                                     HMM=TRUE, num_threads=15,no_plot=TRUE)
tmp <- all_slices_filt
row.names(tmp@meta.data) <- sapply(row.names(tmp@meta.data), function(x){gsub("-","\\.",x)})
tmp <- infercnv::add_to_seurat(tmp, "C:/Users/nicho/Desktop/NUS_PhD/220204 scRNAseq_slices/infercnv_230410")
cnv_results <- read.table("C:/Users/nicho/Desktop/NUS_PhD/220204 scRNAseq_slices/infercnv_230410/map_metadata_from_infercnv.txt")
row.names(cnv_results) <- sapply(row.names(cnv_results), function(x){gsub("\\.","-",x)})
rm(tmp)
all_slices_filt[["total_mutation"]] <- rowSums(cnv_results[,which(grepl("proportion_scaled_cnv",colnames(cnv_results)))])
length(which(all_slices_filt$total_mutation>0.5561163 & all_slices_filt$final_cell_type2=="Cancer cell"))/length(which(all_slices_filt$final_cell_type2=="Cancer cell"))
length(which(all_slices_filt$total_mutation>0.5561163 & (!all_slices_filt$final_cell_type2=="Cancer cell")))/length(which(!all_slices_filt$final_cell_type2=="Cancer cell"))
all_orig_col <- cell_pal(as.character(Idents(subset(all_slices,Condition=="Original"))), c("Epithelial"="#FF0000","CAF"="#FF8000","Endothelial"="#F74ED6","Pericyte"="#16B251",
                                                                                           "Myocyte"="#0066FF","T cell"="#00FE7F","B cell"="#8C56F8","Plasma cell"="#A6752C","cDC"="#E19FD3",
                                                                                           "Macrophage"="#F2F808","Neutrophil"="#50D2FA","Mast cell"="#023047","Monocyte"="#9A4452","Not annotated"="#A9A9A9"))

# used to produce infercnv plot in Extended Figure 2
ggplot(cbind(Embeddings(all_slices_filt,"umap"),all_slices_filt[["total_mutation"]]), aes(UMAP_1,UMAP_2,color=total_mutation)) +
  geom_point(size=0.4) + scale_color_gradient(low="grey",high="red",name="") +
  theme(panel.background = element_rect(fill=NA,color="black"), axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
  xlim(-16.5,15.5) + ylim(-17,15)

# function used to process NATMI output
create_natmi_diffedge_op <- function(named_folders, top_pct, weight_type){
  all_edge_name <- paste("All_edges_",weight_type,".csv",sep="")
  op <- c("Appeared_", "Disappeared_", "UP-regulated_", "DOWN-regulated_", "Stable_")
  op <- paste(op,weight_type,".csv",sep="")
  names(op) <- c("Appeared","Disappeared","Upregulated","Downregulated","Stable")
  weights <- c()
  result <- list()
  for(i in 1:length(named_folders)){
    all_edges <- read.csv(paste(named_folders[i],all_edge_name,sep="\\"))
    if(i == 1){weights <- c(weights, all_edges$Edge.expression.weight.in.condition.1)}
    weights <- c(weights, all_edges$Edge.expression.weight.in.condition.2)
  }
  if(top_pct>=1){top_pct <- top_pct/100; print("top_pct detected to be in % units, converting to decimal")}
  thresh <- quantile(weights, probs = 1-top_pct)
  for(j in names(named_folders)){
    result[[j]] <- list()
    op_fn <- sapply(op, function(x){paste(named_folders[j],x,sep="\\")})
    for(k in names(op)){
      df <- read.csv(op_fn[k])
      print(paste("Initial edge count for ",j," ",k,": ",as.character(nrow(df)),sep=""))
      if(nrow(df)>0){
        if(k %in% c("Appeared", "Disappeared")){df <- subset(df, Delta.edge.expression.weight > thresh)}
        else if(k %in% c("Upregulated","Downregulated","Stable")){df <- df[which(df$Edge.expression.weight.in.condition.2 > thresh | df$Edge.expression.weight.in.condition.1 > thresh),]}
      }
      print(paste("Final edge count for ",j," ",k,": ",as.character(nrow(df)),sep=""))
      result[[j]][[k]] <- df
    }
  }
  return(result)
}
all_diff_natmi_mean <- create_natmi_diffedge_op("PhaseX"="#insert filename of NATMI output", 20, "mean")

# function used to find ligand-receptor interactions important in cancer
find_cancer_ligands <- function(cancer_fn, edges_fn, senders, targets, top_pct, cancer.type=NA){
  cancer_file <- read.xlsx(cancer_fn)
  if(!is.na(cancer.type[1])){cancer_file <- subset(cancer_file, cancer_type %in% cancer.type)}
  cancer_genes <- unique(cancer_file$symbol)
  edges <- read.csv(edges_fn) #%>% subset(Edge.average.expression.weight > 0)
  result <- list(cell_pair=c(),gene=c(),"function"=c())
  for(j in senders){
    edges1 <- subset(edges, Sending.cluster == j)
    for(i in targets){
      edges2 <- subset(edges1, Target.cluster == i)
      weights <- (edges2$Edge.average.expression.weight * edges2$Edge.average.expression.derived.specificity)^0.5
      thresh <- quantile(weights, probs = 1-top_pct)
      edges2 <- edges2[which(weights>thresh),]
      ligands <- unique(edges2$Ligand.symbol) %>% .[which(. %in% cancer_genes)]
      receptors <- unique(edges2$Receptor.symbol) %>% .[which(. %in% cancer_genes)]
      print(paste(j,"-",i," has ",as.character(length(ligands))," ligands and ",as.character(length(receptors))," receptors", sep=""))
      result$cell_pair <- c(result$cell_pair, rep(paste(j,i,sep="-"),length(ligands)+length(receptors)))
      result$gene <- c(result$gene, ligands, receptors)
      result$`function` <- c(result$`function`, rep("Ligand",length(ligands)), rep("Receptor",length(receptors)))
    }
  }
  return(data.frame(result))
}
cancer_interactions_original <- find_cancer_ligands("# insert filename of excel downloaded from  Network of Cancer Genes containing cancer drivers with annotation and supporting evidence",
                                                    "# insert filename of NATMI output generated for Original tumor",
                                                    c("Epithelial","CAF","Macrophage","CD8 T cell"),c("Epithelial","CAF","Macrophage","CD8 T cell"),0.05,
                                                    cancer.type=c("pan-cancer_adult","squamous_head_and_neck_cancer","esophageal_squamous_carcinoma","nasopharyngeal_carcinoma","oral_squamous_cell_carcinoma"))
# function used to find *preserved* ligand-receptor interactions important in cancer
find_shared_interactions <- function(ligands_df, shared_inter, top_pct){
  result <- list(cell_type_pair=c(), lig_rec_pair=c())
  cond1_weights <- shared_inter$Edge.expression.weight.in.condition.1 * shared_inter$Edge.specificity.weight.in.condition.1
  cond1_thresh <- quantile(cond1_weights, probs=0.6)
  cond2_weights <- shared_inter$Edge.expression.weight.in.condition.2 * shared_inter$Edge.specificity.weight.in.condition.2
  cond2_thresh <- quantile(cond2_weights, probs=0.6)
  shared_inter <- shared_inter[intersect(which(cond1_weights>cond1_thresh),which(cond2_weights>cond2_thresh)),]
  for(i in unique(ligands_df[,1])){
    ligands_df_st <- subset(ligands_df, cell_pair == i)
    sen_tar <- strsplit(i,split="-") %>% .[[1]]
    shared_inter_st <- subset(shared_inter, Sending.cluster==sen_tar[1] & Target.cluster==sen_tar[2])
    weights <- shared_inter_st$Edge.expression.weight.in.condition.1 * shared_inter_st$Edge.specificity.weight.in.condition.1 * shared_inter_st$Edge.expression.weight.in.condition.2 * shared_inter_st$Edge.specificity.weight.in.condition.2
    thresh <- quantile(weights, probs = 1-top_pct)
    shared_inter_st <- shared_inter_st[which(weights>thresh),]
    shared_inter_st <- subset(shared_inter_st,Ligand.symbol %in% ligands_df_st[,2] | Receptor.symbol %in% ligands_df_st[,2])
    result$lig_rec_pair <- c(result$lig_rec_pair, unlist(paste(shared_inter_st$Ligand.symbol, shared_inter_st$Receptor.symbol, sep = "-")))
    result$cell_type_pair <- c(result$cell_type_pair, rep(i,nrow(shared_inter_st)))
  }
  return(data.frame(result))
}
shared_phasex_orig_interactions <- read.csv("# insert filename of NATMI output from comparing PhaseX to Original") %>% subset(., Edge.expression.weight.in.condition.1>0 & Edge.expression.weight.in.condition.2>0)
shared_phasex_orig_cancer_interactions <- find_shared_interactions(cancer_interactions_original,shared_phasex_orig_interactions,0.05)
# above used to generate Supplementary Item 5

# function to determine degree of preservation of ligand-receptor interactions
natmi_diffedge_adjmat_with_weight <- function(processed_diffedge, cell_types, pct, cell_type_col1, cell_type_col2,scale_by=1){
  result <- matrix(NA, nrow=length(cell_types), ncol=length(cell_types))
  change <- list() #columns are senders, rows are targets
  change1 <- list()
  row.names(result) <- cell_types #sender
  colnames(result) <- cell_types #receiver
  for(sen in cell_types){
    list_sen <- lapply(processed_diffedge, function(x){subset(x, Sending.cluster == sen)})
    change[[sen]] <- c()
    change1[[sen]] <- c()
    sen_id <- which(cell_types==sen)
    for(tar in cell_types){
      list_sen_tar <- lapply(list_sen, function(x){subset(x, Target.cluster == tar)})
      cond1 <- 0
      cond2 <- 0
      shared <- 0
      if("Appeared" %in% names(processed_diffedge) && nrow(list_sen_tar$Appeared) > 0){cond2 <- cond2 + nrow(list_sen_tar$Appeared)}
      if("Disappeared" %in% names(processed_diffedge) && nrow(list_sen_tar$Disappeared) > 0){cond1 <- cond1 + nrow(list_sen_tar$Disappeared)}
      if("Upregulated" %in% names(processed_diffedge) && nrow(list_sen_tar$Upregulated) > 0){shared <- shared + nrow(list_sen_tar$Upregulated)}
      if("Downregulated" %in% names(processed_diffedge) && nrow(list_sen_tar$Downregulated) > 0){shared <- shared + nrow(list_sen_tar$Downregulated)}
      if("Stable" %in% names(processed_diffedge) && nrow(list_sen_tar$Stable) > 0){shared <- shared + nrow(list_sen_tar$Stable)}
      cond1 <- cond1 + shared
      # cond2 <- cond2 + shared
      # log2fc <- log2(cond2/cond1)
      pct_preserved <- round((shared/(cond1+cond2))*100)
      pct_preserved1 <- pct_preserved
      pct_preserved <- pct_preserved-(2*(100-pct_preserved))
      if(pct_preserved<10){paste("0",as.character(pct_preserved),sep="")}else if(pct_preserved==100){pct_preserved<-"99"}else{as.character(pct_preserved)}
      # if(log2fc < -0.5849625){change[[sen]] <- c(change[[sen]],"#d3d3d3")} else if(log2fc > 0.5849625){change[[sen]] <- c(change[[sen]],"#e9d8a6")} else {change[[sen]] <- c(change[[sen]],"#660708")}
      # result[sen, tar] <- abs(log2fc)
      # if(pct_preserved < pct){change[[sen]] <- c(change[[sen]],paste(cell_type_col2[sen_id],"99",sep=""))} else {change[[sen]] <- c(change[[sen]],paste(cell_type_col1[sen_id],"80",sep=""))}
      change[[sen]] <- c(change[[sen]],paste(cell_type_col1[sen_id],pct_preserved,sep=""))
      change1[[sen]] <- c(change1[[sen]],pct_preserved1)
      result[sen, tar] <- cond1^scale_by
    }
  }
  change <- data.frame(change) %>% t() #rows are senders, columns are targets; same as result
  change1 <- data.frame(change1) %>% t()
  colnames(change) <- cell_types
  colnames(change1) <- cell_types
  return(list(fc=result, direction=change, actual_direction=change1))
}
all_phasex_d0_adj_weight_mean <- natmi_diffedge_adjmat_with_weight(all_diff_natmi_mean$PhaseX[-5],
                                                                   levels(all_slices_filt$final_cell_type1),0.7,
                                                                   pals::kelly(19)[c(3,4,10,6:8,19,5,11:15,17:19)],#c("#FFE98B","#C9ADCF","#00FAA1","#CDE3F7","#FF799C","#E0D8BE","#CDCDCD","#FFC681","#F2C6D5","#5BC1FF","#FCCABC","#A397C9","#FFD889","#FFF869","#E26F54"),
                                                                   rep("#000000",16))
# circos plot in Figure 4
circlize::chordDiagram(all_phasex_d0_adj_weight_mean$fc, directional=1, col=all_phasex_d0_adj_weight_mean$direction,
                       annotationTrack = c("name","grid"), link.arr.type = "big.arrow", link.arr.length = 0.07,
                       grid.col=pals::kelly(19)[c(3,4,10,6:8,19,5,11:15,17:19)], direction.type = c("arrows"))

# find DEGs of cell subtypes using only cells from original tumor
# Used to produce remaining sheets in Supplementary Item 2
Idents(all_slices_sub$Fibroblast) <- all_slices_sub$Fibroblast$sub_cell_type
all_slices_orig_fib_deg <- FindAllMarkers(subset(all_slices_sub$Fibroblast, Condition=="Original"), logfc.threshold = 0, return.thresh = 1) %>% .[order(.$avg_log2FC, decreasing=T),] %>% .[order(.$cluster),]
Idents(all_slices_sub$`Myeloid (non-Mast cell)`) <- all_slices_sub$`Myeloid (non-Mast cell)`$sub_cell_type
all_slices_orig_macro_deg <- FindAllMarkers(subset(all_slices_sub$`Myeloid (non-Mast cell)`, Condition=="Original" & sub_cell_type %in% c("C1Q high Macrophage","SPP1 high Macrophage","IL1B high Macrophage")), logfc.threshold = 0, return.thresh = 1) %>% .[order(.$avg_log2FC, decreasing=T),] %>% .[order(.$cluster),]
Idents(all_slices_sub$`T cell`) <- all_slices_sub$`T cell`$sub_cell_type
all_slices_orig_tcell_deg <- FindAllMarkers(subset(all_slices_sub$`T cell`, Condition=="Original"), logfc.threshold = 0, return.thresh = 1) %>% .[order(.$avg_log2FC, decreasing=T),] %>% .[order(.$cluster),]
Idents(all_slices_sub$Epithelial) <- all_slices_sub$Epithelial$sub_cell_type
all_slices_orig_epi_deg <- FindAllMarkers(all_slices_sub$Epithelial, logfc.threshold = 0, return.thresh = 1) %>% .[order(.$avg_log2FC, decreasing=T),] %>% .[order(.$cluster),]

orig_fib_gsea_go <- lapply(list(myCAF="myCAF",iCAF="iCAF",apCAF="apCAF"), function(x){GSEA(setNames(all_slices_orig_fib_deg$avg_log2FC[which(all_slices_orig_fib_deg$cluster==x)],nm=all_slices_orig_fib_deg$gene[which(all_slices_orig_fib_deg$cluster==x)]),TERM2GENE = go_gs, pvalueCutoff = 1, eps = 0, seed = 123, pAdjustMethod = "fdr")})
orig_macro_gsea_go <- lapply(list(`IL1B high Macrophage`="IL1B high Macrophage",`C1Q high Macrophage`="C1Q high Macrophage",`SPP1 high Macrophage`="SPP1 high Macrophage"), function(x){GSEA(setNames(all_slices_orig_macro_deg$avg_log2FC[which(all_slices_orig_macro_deg$cluster==x)],nm=all_slices_orig_macro_deg$gene[which(all_slices_orig_macro_deg$cluster==x)]),TERM2GENE = go_gs, pvalueCutoff = 1, eps = 0, seed = 123, pAdjustMethod = "fdr")})
orig_tcell_gsea_go <- lapply(list(`CD8 T cell`="CD8 T cell",`CD4 T cell`="CD4 T cell",Treg="Treg"), function(x){GSEA(setNames(all_slices_orig_tcell_deg$avg_log2FC[which(all_slices_orig_tcell_deg$cluster==x)],nm=all_slices_orig_tcell_deg$gene[which(all_slices_orig_tcell_deg$cluster==x)]),TERM2GENE = go_gs, pvalueCutoff = 1, eps = 0, seed = 123, pAdjustMethod = "fdr")})

orig_fib_gsea_kegg <- lapply(list(myCAF="myCAF",iCAF="iCAF",apCAF="apCAF"), function(x){GSEA(setNames(all_slices_orig_fib_deg$avg_log2FC[which(all_slices_orig_fib_deg$cluster==x)],nm=all_slices_orig_fib_deg$gene[which(all_slices_orig_fib_deg$cluster==x)]),TERM2GENE = kegg_gs, pvalueCutoff = 1, eps = 0, seed = 123, pAdjustMethod = "fdr")})
orig_macro_gsea_kegg <- lapply(list(`IL1B high Macrophage`="IL1B high Macrophage",`C1Q high Macrophage`="C1Q high Macrophage",`SPP1 high Macrophage`="SPP1 high Macrophage"), function(x){GSEA(setNames(all_slices_orig_macro_deg$avg_log2FC[which(all_slices_orig_macro_deg$cluster==x)],nm=all_slices_orig_macro_deg$gene[which(all_slices_orig_macro_deg$cluster==x)]),TERM2GENE = kegg_gs, pvalueCutoff = 1, eps = 0, seed = 123, pAdjustMethod = "fdr")})
orig_tcell_gsea_kegg <- lapply(list(`CD8 T cell`="CD8 T cell",`CD4 T cell`="CD4 T cell",Treg="Treg"), function(x){GSEA(setNames(all_slices_orig_tcell_deg$avg_log2FC[which(all_slices_orig_tcell_deg$cluster==x)],nm=all_slices_orig_tcell_deg$gene[which(all_slices_orig_tcell_deg$cluster==x)]),TERM2GENE = kegg_gs, pvalueCutoff = 1, eps = 0, seed = 123, pAdjustMethod = "fdr")})

gsea_oi_xlsx <- read.xlsx("# filename of excel sheet containing important pathways from each cell subtype", sheet = "Sheet1")
# data is composed of 3 columns, named cell_type, path, and new_path: cell_type - cell subtype name; path - original gene set term name from msigdb; new_path - name used for plotting

# function used to calculate the correlation of gene expression using genes from each cell subtype's enriched biological pathways
gene_expression_corr_subtype <- function(patients, cell_type, gsea, gsea_oi, condition_nm="Condition", cell_type_nm="final_cell_type", baseline){
  counts <- lapply(patients, function(x){as.matrix(x@assays$SCT@counts[,which(x[[cell_type_nm]]==cell_type)])})
  condition <- lapply(patients, function(x){as.character(x[[condition_nm]][which(x[[cell_type_nm]]==cell_type),1])})
  gsea <- subset(gsea, ID %in% gsea_oi$path & qvalue<0.05 & NES>0)
  conditions <- as.list(as.character(unique(unlist(lapply(patients, function(x){x[[condition_nm]]})))))
  names(conditions) <- unlist(conditions)
  count <- as.list(names(patients))
  names(count) <- unlist(count)
  cpm_by_condition <- lapply(count, function(x){lapply(conditions, function(y){rowSums(counts[[x]][,which(condition[[x]]==y)])})})
  cpm_by_condition <- lapply(cpm_by_condition, function(y){lapply(y, function(x){asinh((x*1e6)/sum(x))})})
  thresh <- lapply(cpm_by_condition, function(y){quantile(unlist(y) %>% .[which(.>0)], probs = 0.05)})
  keep_lst <- lapply(count, function(x){unique(unlist(lapply(conditions, function(y){which(cpm_by_condition[[x]][[y]]>thresh[[x]])})))})
  keep <- row.names(counts[[1]])[keep_lst[[1]]]
  for(i in 2:length(keep_lst)){keep <- intersect(keep, row.names(counts[[i]])[keep_lst[[i]]])}
  keep <- keep[order(keep)]
  cpm_by_condition <- lapply(cpm_by_condition, function(y){lapply(y, function(x){x[keep]})})
  corr <- c()
  pathnames <- c()
  rm_path <- 0
  cond_oi <- names(conditions)[-which(names(conditions)==baseline)]
  for(path in gsea_oi$path){
    genes_oi <- unlist(strsplit(gsea[which(gsea$ID==path),"core_enrichment"], split = "/")) %>% intersect(., names(cpm_by_condition[[1]][[1]]))
    if(length(genes_oi)>0){
      no_genes <- paste(" - ",as.character(length(genes_oi)), sep = "")
      pathnames <- c(pathnames, rep(paste(gsea_oi$new_path[which(gsea_oi$path==path)], no_genes, sep = ""),(length(cond_oi)*length(counts))))
      for(cpm_by_condition1 in cpm_by_condition){
        for(cond in cond_oi){
          corr <- c(corr, cor(cpm_by_condition1[[cond]][genes_oi], cpm_by_condition1[[baseline]][genes_oi], method="pearson"))
        }
      }
    } else {rm_path <- rm_path + 1}
  }
  cond_res <- rep(cond_oi, nrow(gsea_oi)-rm_path)
  result <- data.frame(Path=pathnames, Condition=cond_res, Patient=rep(rep(names(counts),each=length(cond_oi)), length(corr)/(length(names(counts))*length(cond_oi))), Correlation=corr)
  return(result)
}

impt_subtype_gsea <- list(myCAF=rbind(orig_fib_gsea_go$myCAF@result,orig_fib_gsea_kegg$myCAF@result),iCAF=rbind(orig_fib_gsea_go$iCAF@result,orig_fib_gsea_kegg$iCAF@result),apCAF=rbind(orig_fib_gsea_go$apCAF@result,orig_fib_gsea_kegg$apCAF@result),
                          `IL1B high Macrophage`=rbind(orig_macro_gsea_go$`IL1B high Macrophage`@result,orig_macro_gsea_kegg$`IL1B high Macrophage`@result),
                          `C1Q high Macrophage`=rbind(orig_macro_gsea_go$`C1Q high Macrophage`@result,orig_macro_gsea_kegg$`C1Q high Macrophage`@result),
                          `SPP1 high Macrophage`=rbind(orig_macro_gsea_go$`SPP1 high Macrophage`@result,orig_macro_gsea_kegg$`SPP1 high Macrophage`@result),
                          `CD8 T cell`=rbind(orig_tcell_gsea_go$`CD8 T cell`@result,orig_tcell_gsea_kegg$`CD8 T cell`@result),
                          `CD4 T cell`=rbind(orig_tcell_gsea_go$`CD4 T cell`@result,orig_tcell_gsea_kegg$`CD4 T cell`@result),
                          Treg=rbind(orig_tcell_gsea_go$Treg@result,orig_tcell_gsea_kegg$Treg@result)) %>% lapply(function(x){subset(x,NES>0 & p.adjust<0.05)}) %>% lapply(function(x){x[order(x$NES,decreasing = T),]})
gene_exp_corr_subtype <- lapply(list(`myCAF`="myCAF", `iCAF`="iCAF", `apCAF`="apCAF", `C1Q high Macrophage`="C1Q high Macrophage",
                                     `IL1B high Macrophage`="IL1B high Macrophage",`SPP1 high Macrophage`="SPP1 high Macrophage",
                                     `CD8 T cell`="CD8 T cell",`CD4 T cell`="CD4 T cell", `Treg`="Treg"),
                                function(x){gene_expression_corr_subtype(list("HN377"=subset(all_slices_filt,PatientID=="HN377"),"HN385"=subset(all_slices_filt,PatientID=="HN385"), "HN390"=subset(all_slices_filt,PatientID=="HN390"), "HN396"=subset(all_slices_filt,PatientID=="HN396")), x,
                                                                         impt_subtype_gsea[[x]], subset(gsea_oi_xlsx, cell_type == x),
                                                                         "Condition", "final_cell_type", "Original")})

# same as above but focusing on cancer cells, and pathways known to be involved in carcinogenesis
gene_expression_corr_epi <- function(patients, cell_type, gs, gs_oi, condition_nm="Condition", cell_type_nm="final_cell_type", baseline){
  counts <- lapply(patients, function(x){as.matrix(x@assays$SCT@counts[,which(x[[cell_type_nm]]==cell_type)])})
  condition <- lapply(patients, function(x){as.character(x[[condition_nm]][which(x[[cell_type_nm]]==cell_type),1])})
  conditions <- as.list(as.character(unique(unlist(lapply(patients, function(x){x[[condition_nm]]})))))
  names(conditions) <- unlist(conditions)
  count <- as.list(names(patients))
  names(count) <- unlist(count)
  cpm_by_condition <- lapply(count, function(x){lapply(conditions, function(y){rowSums(counts[[x]][,which(condition[[x]]==y)])})})
  cpm_by_condition <- lapply(cpm_by_condition, function(y){lapply(y, function(x){asinh((x*1e6)/sum(x))})})
  thresh <- lapply(cpm_by_condition, function(y){quantile(unlist(y) %>% .[which(.>0)], probs = 0.05)})
  keep_lst <- lapply(count, function(x){unique(unlist(lapply(conditions, function(y){which(cpm_by_condition[[x]][[y]]>thresh[[x]])})))})
  keep <- row.names(counts[[1]])[keep_lst[[1]]]
  for(i in 2:length(keep_lst)){keep <- intersect(keep, row.names(counts[[i]])[keep_lst[[i]]])}
  keep <- keep[order(keep)]
  cpm_by_condition <- lapply(cpm_by_condition, function(y){lapply(y, function(x){x[keep]})})
  corr <- c()
  pathnames <- c()
  rm_path <- 0
  cond_oi <- names(conditions)[-which(names(conditions)==baseline)]
  for(path in names(gs_oi)){
    genes_oi <- gs$human_gene_symbol[which(gs$gs_name==path)] %>% intersect(., names(cpm_by_condition[[1]][[1]]))
    if(length(genes_oi)>0){
      no_genes <- paste(" - ",as.character(length(genes_oi)), sep = "")
      pathnames <- c(pathnames, rep(paste(gs_oi[path], no_genes, sep = ""),(length(cond_oi)*length(counts))))
      for(cpm_by_condition1 in cpm_by_condition){
        for(cond in cond_oi){
          corr <- c(corr, cor(cpm_by_condition1[[cond]][genes_oi], cpm_by_condition1[[baseline]][genes_oi], method="pearson"))
        }
      }
    } else {rm_path <- rm_path + 1}
  }
  cond_res <- rep(cond_oi, length(gs_oi)-rm_path)
  result <- data.frame(Path=pathnames, Condition=cond_res, Patient=rep(rep(names(counts),each=length(cond_oi)), length(corr)/(length(names(counts))*length(cond_oi))), Correlation=corr)
  return(result)
}
epi_paths <- c("KEGG_CELL_CYCLE"="Cell Cycle",
               "KEGG_HEDGEHOG_SIGNALING_PATHWAY"="Hedgehog",
               "KEGG_MAPK_SIGNALING_PATHWAY"="MAPK",
               "KEGG_NOTCH_SIGNALING_PATHWAY"="NOTCH",
               "KEGG_P53_SIGNALING_PATHWAY"="p53",
               "KEGG_TGF_BETA_SIGNALING_PATHWAY"="TGFβ",
               "KEGG_VEGF_SIGNALING_PATHWAY"="VEGF",
               "KEGG_WNT_SIGNALING_PATHWAY"="WNT",
               "REACTOME_SIGNALING_BY_HIPPO"="Hippo",
               "REACTOME_SIGNALING_BY_EGFR_IN_CANCER"="EGFR",
               "HALLMARK_IL6_JAK_STAT3_SIGNALING"="IL6-JAK/STAT3",
               "HALLMARK_KRAS_SIGNALING_UP"="KRAS",
               "HALLMARK_TNFA_SIGNALING_VIA_NFKB"="TNFα-NF-κB",
               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"="EMT",
               "HALLMARK_HYPOXIA"="Hypoxia")
Idents(all_slices_filt) <- all_slices_filt$final_cell_type
gene_exp_corr_subtype$Epithelial <- gene_expression_corr_epi(list("HN377"=subset(all_slices_filt,PatientID=="HN377"),"HN385"=subset(all_slices_filt,PatientID=="HN385"), "HN390"=subset(all_slices_filt,PatientID=="HN390"), "HN396"=subset(all_slices_filt,PatientID=="HN396")), "Cancer cell",
                                                             rbind(kegg_gs[which(kegg_gs$gs_name %in% names(epi_paths)),],reactome_gs[which(reactome_gs$gs_name %in% names(epi_paths)),],hallmark_gs[which(hallmark_gs$gs_name %in% names(epi_paths)),]),
                                                             epi_paths,"Condition", "final_cell_type1", "Original")
for(i in 1:length(gene_exp_corr_subtype)){gene_exp_corr_subtype[[i]]$Patient <- factor(gene_exp_corr_subtype[[i]]$Patient, levels=c("HN377","HN385","HN390","HN396"))}

#used to generate radar plots from Figures 3-4 and Extended Figure 3
fmsb::radarchart(gene_exp_corr_subtype_mat$`Epithelial`$Matrix, axistype = 1,
                 # Customize the polygon
                 pcol = c("HN377"="#8ecae6","HN385"="#ff5a5f","HN390"="#023047","HN396"="#e98a15"),  plwd = 2, plty = 1, #pfcol = scales::alpha(c("#457b9d","#ff006e","#ffc857"), c(0.7,0.5,0.3)),
                 # Customize the grid
                 cglcol = "grey", cglty = 5, cglwd = 1,
                 # Customize the axis
                 axislabcol = "black", 
                 # Variable labels
                 vlcex = 1.2, vlabels = NA,
                 caxislabels = rep(NA,5), pty = c(15:17,4) #seq(0.2,1,by=0.2), title = "Epithelial"
)

# function to determine preservation of specific ligand-receptor interactions
natmi_specific_inter_change <- function(int_df, diff_edge_file){
  diff_edge <- lapply(as.list(diff_edge_file), read.csv)
  result <- list(lig_rec=c(), sen_tar=c(), weight=c(), "function"=c())
  sen <- list(sen=as.character(int_df[[2]][which(int_df[[2]][,2]=="Sender"),1]), no_sen=length(which(int_df[[2]][,2]=="Sender")))
  tar <- list(tar=as.character(int_df[[2]][which(int_df[[2]][,2]=="Target"),1]), no_tar=length(which(int_df[[2]][,2]=="Target")))
  sen_tar <- c()
  
  for(lr in 1:nrow(int_df[[1]])){
    lig <- int_df[[1]]$Ligand[lr]
    rec <- int_df[[1]]$Receptor[lr]
    result$lig_rec <- c(result$lig_rec, rep(paste(lig, "-", rec, sep = ""), sen$no_sen * tar$no_tar))
    result$`function` <- c(result$`function`, rep(int_df[[1]]$`Function`[lr], sen$no_sen * tar$no_tar))
    diff_edge_lr <- lapply(diff_edge, function(x){subset(x, Ligand.symbol==lig & Receptor.symbol==rec)})
    for(ss in sen$sen){
      diff_edge_lr_s <- lapply(diff_edge_lr, function(x){subset(x, Sending.cluster==ss)})
      for(tt in tar$tar){
        diff_edge_lr_st <- lapply(diff_edge_lr_s, function(x){subset(x, Target.cluster==tt)})
        result$weight <- c(result$weight, unlist(lapply(diff_edge_lr_st, function(x){
          if(nrow(x)>0){
            sapply(1:nrow(x),function(y){if(x$Edge.expression.weight.in.condition.1[y]>0){
              if(x$Edge.expression.weight.in.condition.2[y]>0){"Preserved"}else{"Lost"}
            } else if(x$Edge.expression.weight.in.condition.2[y]>0){"Gained"}else{"Not detected"}
            })
          } else {"Not detected"}})))
        if(lr == 1){
          sen_tar <- c(sen_tar, paste(ss, "-", tt, sep = ""))
        }
      }
    }
  }
  result$sen_tar <- rep(sen_tar, nrow(int_df[[1]]))
  result$weight <- factor(result$weight, levels=c("Preserved","Gained","Lost","Not detected"))
  return(data.frame(result))
}
combine_natmi_specific_inter_change <- function(int_df, diff_edge_file_lst){
  result <- natmi_specific_inter_change(int_df = int_df, diff_edge_file = diff_edge_file_lst[[1]])
  if(length(diff_edge_file_lst)>1){
    for(i in 2:length(diff_edge_file_lst)){
      result <- rbind(result, natmi_specific_inter_change(int_df = int_df, diff_edge_file = diff_edge_file_lst[[i]]))
    }
  }
  result$Condition <- factor(rep(names(diff_edge_file_lst), each = nrow(result)/length(diff_edge_file_lst)), levels = names(diff_edge_file_lst))
  return(result)
}
all_immune_inter_plot_mean <- combine_natmi_specific_inter_change(immune_interactions, #data frame of immune checkpoint interactions (has 2 columns with names Ligand and Receptor. Each containing respective gene names of interacting proteins involved in immune checkpoints)
                                                                  list(PhaseX="C:\\Users\\nicho\\Desktop\\NUS_PhD\\220204 scRNAseq_slices\\NATMI\\NATMI\\data_230410\\All_PhaseX-Original_mean\\Delta_edges_lrc2p\\All_edges_mean.csv"))

# used to generate Figure 5B
subset(all_immune_inter_plot_mean,Condition=="PhaseX") %>% ggplot(aes(x=lig_rec, y=sen_tar, shape=weight, fill=weight)) +
  geom_point(stroke = 0.1, size = 2) + xlab("") + ylab("Sender-Target") +
  scale_shape_manual(values = c(23,25,1), name = "Change") + scale_color_manual(values = c("white","white","black")) +
  scale_fill_manual(values = c("#fe4a49","#009fb7","#ffffff"), name = "Change") +
  theme(axis.text.y = element_text(angle=0, family="sans", hjust=1, size=8, color="black"), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.1, linetype = "solid"),
        panel.grid = element_blank(), axis.line = element_blank(), axis.text.x = element_text(angle=90, family="sans", hjust=1, size=8, color="black",face="italic"),
        axis.title = element_text(family = "sans", size=9), strip.text.y.left = element_blank(),
        strip.background = element_blank(), legend.text = element_text(size=8), legend.title = element_text(size=8),
        legend.key.height = unit(0.3,"cm"), legend.position = "none")

# function used to generate pseudobulk RNA data from scRNA-seq data
pseudo_bulk <- function(seu, cluster_col){
  cell_types <- unique(as.character(seu@meta.data[,cluster_col]))
  cluster_list <- as.list(cell_types)
  cluster_list <- lapply(cluster_list, function(x){rowSums(seu@assays$SCT@counts[,which(seu@meta.data[,cluster_col]==x)])})
  names(cluster_list) <- cell_types
  cluster_list <- data.frame(cluster_list)
  to_rm <- which(rowSums(cluster_list)==0)
  cluster_list <- cluster_list[-to_rm,]
  row.names(cluster_list) <- row.names(seu@assays$SCT@counts)[-to_rm]
  cluster_list <- t(asinh(t(cluster_list*1e6)/colSums(cluster_list)))
  cluster_list <- cluster_list/apply(cluster_list,1,max)
  return(cluster_list)
}
all_slices_filt_psuedobulk_ha <- pseudo_bulk(subset(all_slices_filt,Condition=="PhaseX"), "final_cell_type1")
all_slices_filt_psuedobulk_orig <- pseudo_bulk(subset(all_slices_filt,Condition=="Original"), "final_cell_type1")
all_slices_filt_psuedobulk_ha <- all_slices_filt_psuedobulk_ha[,c("Cancer.cell","Fibroblast","Pericyte","Myocyte","Endothelial","CD8.T.cell","CD4.T.cell","Treg","Naive.Memory.T.cell","B.cell","Plasma.cell","Macrophage","Dendritic.cell","Monocyte","Neutrophil","Mast.cell")]
all_slices_filt_psuedobulk_orig <- all_slices_filt_psuedobulk_orig[,c("Cancer.cell","Fibroblast","Pericyte","Myocyte","Endothelial","CD8.T.cell","CD4.T.cell","Treg","Naive.Memory.T.cell","B.cell","Plasma.cell","Macrophage","Dendritic.cell","Monocyte","Neutrophil","Mast.cell")]
# used to generate Figure 5A
Heatmap(all_slices_filt_psuedobulk_orig[unique(c(immune_interactions$LR$Ligand,immune_interactions$LR$Receptor)),], cluster_rows = F,
        col = circlize::colorRamp2(c(0, 0.5, 1), c("#000022", "#4b3f72", "#ffd800")), name = "Relative Expression", cluster_columns = F,
        #rect_gp = gpar(col = "black", lwd = 0.2),
        heatmap_legend_param = list(legend_height = unit(3,"cm")), column_names_rot = 90, show_row_names = T, show_column_names = T)
Heatmap(all_slices_filt_psuedobulk_ha[unique(c(immune_interactions$LR$Ligand,immune_interactions$LR$Receptor)),], cluster_rows = F,
        col = circlize::colorRamp2(c(0, 0.5, 1), c("#000022", "#4b3f72", "#ffd800")), name = "Relative Expression", cluster_columns = F,
        #rect_gp = gpar(col = "black", lwd = 0.2),
        heatmap_legend_param = list(legend_height = unit(3,"cm")), column_names_rot = 90, show_row_names = T, show_column_names = T)

# function used to create matrix suitable for heatmap plotting
my_sc_heatmap <- function(seu, cluster_col, condition_col, color){
  clusters <- unique(as.character(seu@meta.data[,cluster_col]))
  meta <- seu@meta.data
  for(i in length(color):1){
    meta <- meta[order(meta[,names(color)[i]]),]
  }
  # arranged_ids <- unlist(lapply(as.list(clusters),function(x){which(seu@meta.data[,cluster_col]==x)}))
  mat <- seu@assays$SCT@counts[,row.names(meta)]
  mat <- t(asinh(t(mat*1e6)/colSums(mat)))
  if(min(mat)<0){mat <- mat+min(mat)}
  mat <- mat/apply(mat,1,max)
  annot_df <- data.frame(name=meta[,cluster_col],name1=meta[,condition_col])
  colnames(annot_df) <- names(color)
  annot <- rowAnnotation(df = annot_df, col = color)
  return(list("Matrix"=as.matrix(mat),"Annotation"=annot))
}
all_slices_fib_deg_filt <- all_slices_orig_fib_deg[order(all_slices_orig_fib_deg$avg_log2FC, decreasing=T),] %>% .[-which(duplicated(.$gene)),] %>% .[order(.$cluster),] %>% subset(avg_log2FC>0.58 & p_val_adj<0.01 & pct.1 > 0.15)
all_slices_filt_caf_hm <- my_sc_heatmap(subset(all_slices_sub$Fibroblast,sub_cell_type1=="CAF"),
                                        "sub_cell_type","Condition",
                                        color = list("sub_cell_type"=c("myCAF"="#fed700","iCAF"="#21b0fe","apCAF"="#ff218c"),
                                                     "Condition"=c("Original"="#003049","PhaseX"="#d62828")))
# used to plot Figure 3G
Heatmap(t(all_slices_filt_caf_hm$Matrix[c(subset(all_slices_fib_deg_filt,cluster=="myCAF")$gene[1:15],subset(all_slices_fib_deg_filt,cluster=="iCAF")$gene[1:15],subset(all_slices_fib_deg_filt,cluster=="apCAF")$gene[1:15]),]), cluster_rows = F,
        col = circlize::colorRamp2(c(0, 0.5, 1), c("#000022", "#4b3f72", "#ffd800")), name = "Relative\nExpression", cluster_columns = F,
        column_split = factor(rep(c("myCAF","iCAF","apCAF"),each=8),levels=c("myCAF","iCAF","apCAF")), column_gap = unit(2,"mm"),
        heatmap_legend_param = list(legend_height = unit(3,"cm")), column_names_rot = 90, left_annotation = all_slices_filt_caf_hm$Annotation,
        show_column_names = T, show_row_names = F, row_names_side = "left", use_raster=T, column_names_gp = gpar(fontface = "italic"),
        row_split = factor(c(rep("myCAF",3914),rep("iCAF",1240),rep("apCAF",494)),levels=c("myCAF","iCAF","apCAF")), row_gap = unit(2,"mm"))

