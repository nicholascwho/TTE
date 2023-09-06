suppressPackageStartupMessages({
  library(openxlsx)
  library(DESeq2)
  library(ggplot2)
  library(ComplexHeatmap)
  library(dplyr)
  library(FactoMineR)
  library(factoextra)
  library(apeglm)
  library(Seurat)
  library(patchwork)
  library(sva)
  library(tidyverse)
  library(clusterProfiler)
  library(circlize)
  library(enrichplot)
  library(slingshot)
  library(tidyverse)
  library(hrbrthemes)
  library(viridis)
  library(networkD3)
  library(nichenetr)
  library(singleCellHaystack)
  library(scriabin)
  library(loomR)
  library(sf)
  library(msigdbr)
  library(scales)
  library(survival)
  library(survminer)})

#scRNA-seq
all <- CreateSeuratObject(counts, meta, min.cells = 3) # counts and meta represent counts and metadata of downloaded data
all[["percent_mt"]] <- PercentageFeatureSet(all, pattern="^MT-")
all <- subset(all, nCount_RNA > 200 & nCount_RNA < 50000 & nFeature_RNA > 100 & nFeature_RNA < 6000 & percen_mt < 0.15)
all <- SCTransform(all)
all <- RunPCA(all)
all <- RunHarmony(all, group.by.vars = "origin", verbose = T, assay.use = "SCT", max.iter.harmony = 15, tau = 300)
ElbowPlot(all, reduction = "harmony")
all <- RunUMAP(all, reduction = "harmony", verbose = F, assay = "SCT", dims = 1:14)
all <- FindNeighbors(all, dims = 1:2, reduction = "umap")
all <- FindClusters(all, res = 0.2)
pri <- subset(all, site == "P") # subset cells for primary site cells
pri_deg <- FindAllMarkers(pri, logfc.threshold = 0, return.thresh = 1) %>% .[order(.$avg_log2FC, decreasing=T),] %>% .[order(.$cluster),]
## rename clusters under column name cell_type
pri_fib <- subset(pri, cell_type == "CAF")
pri_fib1 <- RunUMAP(pri_fib, reduction = "harmony", dims = 1:14)
set.seed(123)
pri_fib_haystack <- haystack(x = Embeddings(pri_fib1, reduction = "harmony")[,1:14], expression = pri_fib1@assays$SCT@data[which(rowSums(pri_fib1@assays$SCT@data)>0),])
pri_fib_haystack_sorted_filt <- show_result_haystack(res.haystack = pri_fib_haystack, n = length(which(pri_fib_haystack$results$log.p.adj < -3))) #p.adj<0.01
pri_fib_haystack_hclust <- hclust_haystack(Embeddings(pri_fib1, reduction = "harmony")[,1:14], pri_fib1@assays$SCT@data[row.names(pri_fib_haystack_sorted_filt), ], grid.coordinates=pri_fib_haystack$info$grid.coordinates)
pri_fib_haystack_hclusters <- cutree(pri_fib_haystack_hclust, k=4) # check till 10
d <- cbind(Embeddings(pri_fib1, reduction = "umap"), t(as.matrix(pri_fib1@assays$SCT@data[row.names(pri_fib_haystack_sorted_filt),]))) %>% as.data.frame()
for (cluster in unique(pri_fib_haystack_hclusters)) {
  d[[paste0("cluster_", cluster)]] <- colMeans(pri_fib1@assays$SCT@data[names(which(pri_fib_haystack_hclusters == cluster)), ])
}

# to visualize expression of haystack gene sets
svg(filename, height = 7, width = 7) # replace filename with desired path
lapply(tail(colnames(d),4), function(cluster) {
  ggplot(d, aes(UMAP_1, UMAP_2, color=.data[[cluster]])) +
    geom_point(size=0.4) +
    scale_color_gradient("Cluster",low="#c1c1c1",high="#db2b39") +
    theme(panel.background = element_rect(fill=NA,color="black"), axis.line = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          legend.position = "none")
}) |> patchwork::wrap_plots()
dev.off()
pri_fib1 <- FindNeighbors(pri_fib1, dims = 1:2, reduction = "umap")
pri_fib1 <- FindClusters(pri_fib1, res = 0.04)
DimPlot(pri_fib1, label=T)
pri_fib_deg_final <- FindAllMarkers(pri_fib1, logfc.threshold = 0, return.thresh = 1, min.pct = 0.001) %>% .[order(.$avg_log2FC,decreasing=T), ] %>% .[order(.$cluster), ]
get_haystack_deg <- function(deg_df, haystack_clusters, cluster_annot, pval_thresh=0.05, log2fc_thresh=1){
  result <- list()
  for(deg_clus in levels(deg_df$cluster)){
    result[[deg_clus]] <- list()
    result[[deg_clus]]$up <- intersect(subset(deg_df, cluster==deg_clus & avg_log2FC > log2fc_thresh & p_val_adj < pval_thresh)$gene, 
                                       names(haystack_clusters)[which(haystack_clusters==cluster_annot[deg_clus])])
    result[[deg_clus]]$down <- intersect(subset(deg_df, cluster==deg_clus & avg_log2FC < 0-log2fc_thresh & p_val_adj < pval_thresh)$gene, 
                                         names(haystack_clusters)[which(!haystack_clusters==cluster_annot[deg_clus])])
  }
  return(result)
}

fib_top_deg <- get_haystack_deg(pri_fib_deg_final, pri_fib_haystack_hclusters, c("0"=1,"1"=3,"2"=4,"3"=2), log2fc_thresh = 0.25) 

prop_barplot <- function(seu, ct, patient){
  # seu: seurat object
  # ct: column name containing info on cluster of interest
  # patient: column name containing info on patient of origin
  result <- list(Type=c(),Patient=c(),Proportion=c())
  for(i in unique(as.character(seu@meta.data[,patient]))){
    tmp <- seu@meta.data[which(seu@meta.data[,patient]==i),]
    for(j in unique(as.character(seu@meta.data[,ct]))){
      tmp1 <- tmp[which(tmp[,ct]==j),]
      result$Type <- c(result$Type, j)
      result$Patient <- c(result$Patient, i)
      result$Proportion <- c(result$Proportion, nrow(tmp1)*100/nrow(tmp))
    }
  }
  return(as.data.frame(result))
}
pri_fib_prop <- prop_barplot(pri_fib1,"caf_cell_type","patientID")
pri_fib_prop$Type <- factor(pri_fib_prop$Type,levels=names(caf_colours))
pri_fib_prop$Patient <- factor(pri_fib_prop$Patient,levels=unique(pri_fib_prop$Patient) %>% .[order(.)])

# to visualize CAF state distribution
svg(filename, # replace filename with desired path
    width=9, height=6)
ggplot(pri_fib_prop, aes(x=Patient, y=Proportion, fill=Type)) +
  geom_bar(colour = "black", stat = "identity", show.legend = T) + ylab("Proportion (%)") + xlab("") + #coord_polar("y",start=0) +
  scale_fill_manual(values = caf_colours) +
  theme(axis.text = element_text(size=14), axis.ticks = element_line(linewidth=1),
        panel.background = element_rect(fill=NA,color="black",linewidth=1),
        panel.grid = element_blank(),
        axis.title = element_blank(), legend.text = element_text(size=10))
dev.off()

create_ranked_named_gene_for_gsea <- function(deg_df, clus_col, logfc_col, gene_col){
  # deg_df: DEG output
  # clus_col: name of column containing cluster information
  # logfc_col: name of column containing log2FC information
  # gene_col: name of column containing gene information
  result <- list()
  deg_df <- deg_df[order(deg_df[,logfc_col], decreasing = T),] %>% .[order(.[,clus_col], decreasing = F),]
  for(i in unique(deg_df[[clus_col]])){
    result[[i]] <- as.numeric(deg_df[which(deg_df$cluster==i),logfc_col])
    names(result[[i]]) <- as.character(deg_df[which(deg_df$cluster==i),gene_col])
  }
  return(result)
}

go_gs <- msigdbr(category = "C5", subcategory = "BP") %>% select(., gs_name, human_gene_symbol)
kegg_gs <- msigdbr(category = "C2", subcategory = "KEGG") %>% select(., gs_name, human_gene_symbol)
reactome_gs <- msigdbr::msigdbr(category = "C2", subcategory = "REACTOME") %>% .[, c("gs_name", "human_gene_symbol")]
other_caf_markers1 <= read.xlsx("path to collection of gene signatures")

ranked_gopal_fib_genes <- create_ranked_named_gene_for_gsea(pri_fib_deg_final,"cluster","avg_log2FC","gene")
gsea_other_caf0 <- GSEA(ranked_gopal_fib_genes$`0`, pAdjustMethod = "fdr", TERM2GENE = other_caf_markers1[,2:3], seed = 123, pvalueCutoff = 1, eps = 0, minGSSize = 5)
gsea_other_caf1 <- GSEA(ranked_gopal_fib_genes$`1`, pAdjustMethod = "fdr", TERM2GENE = other_caf_markers1[,2:3], seed = 123, pvalueCutoff = 1, eps = 0, minGSSize = 5)
gsea_other_caf2 <- GSEA(ranked_gopal_fib_genes$`2`, pAdjustMethod = "fdr", TERM2GENE = other_caf_markers1[,2:3], seed = 123, pvalueCutoff = 1, eps = 0, minGSSize = 5)
gsea_other_caf3 <- GSEA(ranked_gopal_fib_genes$`3`, pAdjustMethod = "fdr", TERM2GENE = other_caf_markers1[,2:3], seed = 123, pvalueCutoff = 1, eps = 0, minGSSize = 5)
gsea_path_caf0 <- GSEA(ranked_gopal_fib_genes$`0`, pAdjustMethod = "fdr", TERM2GENE = rbind(go_gs,rbind(kegg_gs,reactome_gs)), seed = 123, pvalueCutoff = 1, eps = 0)
gsea_path_caf1 <- GSEA(ranked_gopal_fib_genes$`1`, pAdjustMethod = "fdr", TERM2GENE = rbind(go_gs,rbind(kegg_gs,reactome_gs)), seed = 123, pvalueCutoff = 1, eps = 0)
gsea_path_caf2 <- GSEA(ranked_gopal_fib_genes$`2`, pAdjustMethod = "fdr", TERM2GENE = rbind(go_gs,rbind(kegg_gs,reactome_gs)), seed = 123, pvalueCutoff = 1, eps = 0)
gsea_path_caf3 <- GSEA(ranked_gopal_fib_genes$`3`, pAdjustMethod = "fdr", TERM2GENE = rbind(go_gs,rbind(kegg_gs,reactome_gs)), seed = 123, pvalueCutoff = 1, eps = 0)

dotplot_for_gsea <- function(gsea_lst, p_thresh=0.001, term_start=6){
  # gsea_lst: named list of gsea output
  # p_thresh: adjusted p value threshold
  # term_start: when the actual name of the gene set starts
  result <- list("Term"=c(),"NES"=c(),"padj"=c(),"Cluster"=c())
  for(i in 1:length(gsea_lst)){
    filt <- subset(gsea_lst[[i]], p.adjust<p_thresh)
    result$Term <- c(result$Term, sapply(filt$Description, function(x){gsub("_", " ", x) %>% substr(., term_start, nchar(.))}) %>% as.character())
    result$NES <- c(result$NES, filt$NES)
    result$padj <- c(result$padj, filt$p.adjust)
    result$Cluster <- c(result$Cluster, rep(names(gsea_lst)[i], nrow(filt)))
  }
  result <- data.frame(result)
  unique_terms <- unique(result$Term)
  for(i in unique_terms){
    filt <- which(result$Term==i)
    if(all(result[filt,"NES"] < rep(0,length(filt)))){
      result <- result[-filt,]
    }
  }
  unique_terms <- unique(result$Term)
  result <- result[order(result$NES, decreasing = T),]
  term_order <- c()
  for(i in unique_terms){
    term_order <- c(term_order, which(result$Term==i)[1])
  }
  tmp <- result[term_order,]
  result <- rbind(tmp, result[-term_order,])
  return(result)
}
vivo_other_gsea_dotplot <- dotplot_for_gsea(list("CAF-0"=gsea_other_caf0@result,"CAF-1"=gsea_other_caf1@result,"CAF-2"=gsea_other_caf2@result,"CAF-3"=gsea_other_caf3@result),p_thresh = 0.05,term_start = 1)
vivo_other_gsea_dotplot$Cluster <- factor(vivo_other_gsea_dotplot$Cluster, levels = caf_names)
vivo_other_gsea_dotplot$color <- map_color(vivo_other_gsea_dotplot$Cluster, list("CAF-0"="#ee6055","CAF-1"="#ffd97d","CAF-2"="#90c2e7","CAF-3"="#4e8098"))
vivo_other_gsea_dotplot$Term1 <- factor(vivo_other_gsea_dotplot$Term1, levels = c("Puram CAF1 - HNSCC", "Elyada myCAF - PDAC", "Dominguez LRRC15 TGFB CAF - PDAC", "Kieffer ecm myCAF - BC", "Kieffer TGFb myCAF - BC", "Kieffer wound myCAF - BC", "Wang myCAF - PDAC", "Wu myCAF - BC","Chung myCAF - NSCLC","Hanley Myofibroblast - NSCLC",
                                                                                  "Puram CAF2 - HNSCC", "Elyada iCAF - PDAC", "Dominguez IL1 CAF - PDAC", "Wang iCAF - PDAC", "Wu iCAF - BC", "Kieffer detox iCAF - BC", "Kieffer IL iCAF - BC","Chung iCAF - NSCLC",
                                                                                  "Wang meCAF - PDAC","Dominguez Nonmalignant Fibroblast - PDAC","Kieffer Normal fibroblast - BC","Wu ImmaturePVL - BC","Givel CAF S1 - Ovarian","Givel CAF S4 - Ovarian","Hanley Adventitial - NSCLC","Hanley Alveolar - NSCLC","Elyada KPC apCAF - PDAC","Zhu Adipose Stromal Cell - Pan"))
vivo_other_gsea_dotplot$CAF <- factor(vivo_other_gsea_dotplot$CAF, levels = c("myCAF", "iCAF", "Other"))
vivo_other_gsea_dotplot_final <- subset(vivo_other_gsea_dotplot, NES > 0)

circos_mat_for_gsea <- function(df, sector_column){
  # df: output of dotplot_for_gsea
  # sector_column: name of column containing info gene sets that belong to the same sector
  result <- list()
  for(i in unique(df[[sector_column]])){
    tmp <- df[which(df[[sector_column]]==i),]
    my_terms <- unique(tmp$Term) %>% .[order(.)]
    my_clusters <- unique(tmp$Cluster) %>% .[order(.)]
    no_row <- length(my_terms)
    no_col <- length(my_clusters)
    curr_mat <- matrix(0, nrow = no_row, ncol = no_col)
    row.names(curr_mat) <- my_terms
    colnames(curr_mat) <- my_clusters
    for(j in 1:no_row){
      for(k in 1:no_col){
        tmp1 <- subset(tmp, Cluster == my_clusters[k] & Term == my_terms[j])
        if(nrow(tmp1) > 0){
          curr_mat[j,k] <- tmp1$NES[1]
        }
      }
    }
    result[[i]] <- curr_mat
  }
  return(result)
}
vivo_other_gsea_circos <- circos_mat_for_gsea(vivo_other_gsea_dotplot_final, "CAF")
# to visualize circos plot comparing in vivo caf states to CAFs in other cancer types
svg(filename, height = 7, width = 7) # replace filename with desired path
circos.clear()
circos.par("track.height" = 0.5, gap.after=6)
circos.initialize(sectors = vivo_other_gsea_dotplot_final$CAF, x = vivo_other_gsea_dotplot_final$Term1)
circos.track(sectors = vivo_other_gsea_dotplot_final$CAF, ylim = c(0,4), bg.col = c("#ffcfd2","#e3f2fd","#fefae0"), bg.lwd = 5, bg.border = "white", panel.fun = function(x,y) {
  sector.index = get.cell.meta.data("sector.index")
  xcenter = get.cell.meta.data("xcenter")
  ycenter = get.cell.meta.data("ycenter")
  # circos.text(xcenter, ycenter, sector.index,
  #             cex = 1.5,
  #             niceFacing =T,
  #             facing="bending.inside", adj = c(0.9,-3))
})
circos.barplot(vivo_other_gsea_circos$myCAF, seq(1.47,9.54,1.151429), col = c("#ef233c"), sector.index = "myCAF", track.index = 1, lwd = 1, border = "black")
circos.yaxis("right", labels.niceFacing=T, labels.cex=0.7, sector.index = "myCAF")
circos.barplot(vivo_other_gsea_circos$iCAF, seq(11.48,17.5,0.86), col = c("#ffba08","#3f88c5"), sector.index = "iCAF", track.index = 1, lwd = 1, border = "black")
circos.yaxis("right", labels.niceFacing=F, labels.cex=0.7, sector.index = "iCAF")
circos.barplot(vivo_other_gsea_circos$Other, seq(19.4,27.6,0.9111111), col = caf_colours, sector.index = "Other", track.index = 1, lwd = 1, border = "black")
circos.yaxis("right", labels.niceFacing=F, labels.cex=0.7, sector.index = "Other")
circos.axis(h = "bottom", labels = 1:8, major.at = seq(1.47,9.54,1.151429), sector.index = "myCAF", direction = "inside", labels.facing = "inside", labels.niceFacing = T, major.tick = F, col = NA, labels.col = "#000000", labels.cex = 0.7)
circos.axis(h = "bottom", labels = 9:16, major.at = seq(11.48,17.5,0.86), sector.index = "iCAF", direction = "inside", labels.facing = "inside", labels.niceFacing = T, major.tick = F, col = NA, labels.col = "#000000", labels.cex = 0.7)
circos.axis(h = "bottom", labels = 17:26, major.at = seq(19.4,27.6,0.9111111), sector.index = "Other", direction = "inside", labels.facing = "inside", labels.niceFacing = F, major.tick = F, col = NA, labels.col = "#000000", labels.cex = 0.7)
circos.clear()
dev.off()

create_gsea_sankey_df <- function(filename, sheetnames, cluster_names){
  # filename: path to excel containing GSEA output, including column that groups gene sets into larger groups
  # sheetnames: which sheets to include, with each sheet containing GSEA output for each state
  # cluster_names: what to name the cell state/type; should correspond to sheetnames
  for(i in 1:length(sheetnames)){
    if(i==1){
      df <- read.xlsx(filename, sheetnames[i]) %>% .[which(!is.na(.$Type)),]
      curr_nrow <- nrow(df)
      clus <- rep(cluster_names[i], curr_nrow)
    } else {
      tmp <- read.xlsx(filename, sheetnames[i]) %>% .[which(!is.na(.$Type)),]
      df <- rbind(df, tmp)
      clus <- c(clus, rep(cluster_names[i], nrow(df)-curr_nrow))
      curr_nrow <- nrow(df)
    }
  }
  df$Cluster <- clus
  result <- list(Cluster=c(), Process=c(), No_Pathways=c())
  for(i in cluster_names){
    for(j in unique(df$Type)){
      tmp <- subset(df, Cluster == i & Type == j)
      if(!nrow(tmp) == 0){
        result$Cluster <- c(result$Cluster, i)
        result$Process <- c(result$Process, j)
        result$No_Pathways <- c(result$No_Pathways, nrow(tmp))
      }
    }
  }
  result <- data.frame(result)
  node <- data.frame(Names = c(unique(result$Cluster), unique(result$Process)))
  result$ClusterID <- match(result$Cluster, node$Names) - 1
  result$ProcessID <- match(result$Process, node$Names) - 1
  result$link_type <- node[result$ClusterID + 1, "Names"]
  return(list(df=result, nodes=node))
}
gsea_df_for_sankey <- create_gsea_sankey_df("C:\\Users\\nicho\\Desktop\\NUS_PhD\\220117 HN294P CAF RNAseq\\230606\\GSEA 4 clusters\\CAF Results.xlsx",
                                            c("CAF0", "CAF1", "CAF2", "CAF3"),
                                            c("CAF-0", "CAF-1", "CAF-2", "CAF-3"))
color_scale <- 'd3.scaleOrdinal() .range(["#ef233c","#ffba08","#3f88c5","#032b43","#ffe66d","#ffe66d","#ffe66d","#ffe66d","#ffe66d","#ffe66d","#ffe66d","#ffe66d","#ffe66d","#ffe66d","#ffe66d","#ffe66d","#ffe66d","#ffe66d","#ffe66d"])'
sn <- sankeyNetwork(Links = gsea_df_for_sankey$df, Nodes = gsea_df_for_sankey$nodes,
                    Source = "ClusterID", Target = "ProcessID",
                    Value = "No_Pathways", NodeID = "Names", 
                    sinksRight=T, nodeWidth=20, fontSize=15, nodePadding=10, height = 850, width = 700, margin= list("left"=200),
                    fontFamily = "Arial", LinkGroup = "link_type", colourScale = color_scale)
sn <- htmlwidgets::onRender(sn,
                            '
  function(el, x) {
    d3.selectAll(".node text").attr("text-anchor", "begin").attr("x", 30);
  }
  '
)
# to visualize sankey diagram for biological pathway GSEA
saveNetwork(sn, "desired filename")
webshot::webshot("desired filename","desired filename", vwidth = 950, vheight = 700)


## for pySCENIC; replace filename with desired path
tmp <- pri_fib1@assays$SCT@counts[which(rowSums(pri_fib1@assays$SCT@counts > 0) > 3),]
create(filename = filename, data = tmp, transpose = T, feature.attrs = list(Genes=row.names(tmp)), cell.attrs = list(ID=colnames(tmp),nCount=colSums(tmp),nFeature=apply(tmp,2,function(x){length(which(x>0))})), overwrite = T)
loom_file <- connect(filename = filename, mode = "r+")
loom_file$close_all()
write.table(log1p(tmp), file = filename, quote = F, sep = "\t")
pri_fib_pyscenic_auc <- read.csv("path to pySCENIC AUC output", row.names = 1)
colnames(pri_fib_pyscenic_auc) <- unlist(sapply(colnames(pri_fib_pyscenic_auc),function(x){substr(x,1,nchar(x)-3)}))
pri_fib_pyscenic_seu <- CreateSeuratObject(counts = Matrix::Matrix(t(pri_fib_pyscenic_auc), sparse = TRUE), meta.data = pri_fib1@meta.data[,-c(1:3)], min.cells = 3)
VlnPlot(pri_fib_pyscenic_seu, "nFeature_RNA")
pri_fib_pyscenic_seu <- subset(pri_fib_pyscenic_seu, nFeature_RNA>190 & nFeature_RNA<225)
pri_fib_pyscenic_seu <- NormalizeData(pri_fib_pyscenic_seu, normalization.method = "RC", scale.factor = 100, verbose = F, assay = "RNA")
pri_fib_pyscenic_seu <- ScaleData(pri_fib_pyscenic_seu, assay = "RNA", features = row.names(pri_fib_pyscenic_seu@assays$RNA@data))
pri_fib_pyscenic_seu <- RunPCA(pri_fib_pyscenic_seu, assay = "RNA", features = row.names(pri_fib_pyscenic_seu@assays$RNA@data))
pri_fib_pyscenic_seu <- harmony::RunHarmony(pri_fib_pyscenic_seu, group.by.vars = "origin", verbose = T,max.iter.harmony = 10, assay.use = "RNA", tau = 500)
ElbowPlot(pri_fib_pyscenic_seu, reduction = "harmony")
pri_fib_pyscenic_seu <- RunUMAP(pri_fib_pyscenic_seu, reduction = "harmony", verbose = T, assay = "RNA", dims = 1:7)
pri_fib_pyscenic_seu <- FindNeighbors(pri_fib_pyscenic_seu, dims = 1:2, reduction = "umap")
pri_fib_pyscenic_seu <- FindClusters(pri_fib_pyscenic_seu, res = 0.04)
DimPlot(pri_fib_pyscenic_seu, label = "T", group.by = "SCT_snn_res.0.2") + DimPlot(pri_fib_pyscenic_seu, label = "T", group.by = "caf_cell_type") + DimPlot(pri_fib_pyscenic_seu, group.by = "origin")
Idents(pri_fib_pyscenic_seu) <- pri_fib_pyscenic_seu$caf_cell_type
pri_fib_pyscenic_datf <- FindAllMarkers(pri_fib_pyscenic_seu, assay = "RNA", logfc.threshold = 0, return.thresh = 1, min.pct = 0.01, slot = "counts") %>% .[order(.$avg_log2FC, decreasing = T),] %>% .[order(.$cluster),]
pyscenic_caf_cluster_colors <- cell_pal(Idents(pri_fib_pyscenic_seu), c("myCAF"="#ef233c","IFN-iCAF"="#ffba08","IL-iCAF"="#3f88c5","Adipose-like CAF"="#032b43"))
pyscenic_cluster_colors <- cell_pal(as.character(pri_fib_pyscenic_seu$SCT_snn_res.0.2), c("0"="#33a8c7","1"="#52e3e1","2"="#a0e426","3"="#fdf148","4"="#ffab00","5"="#f77976","6"="#f050ae","7"="#00043a","8"="#d883ff","9"="#9336fd"))
svg("C:/Users/nicho/Desktop/NUS_PhD/220117 HN294P CAF RNAseq/230606/Pri Fib pyscenic UMAP.svg", height = 7, width = 14)
par(mfrow=c(1,2), mai = c(0.1,0.1,0.1,0.1))
plot(Embeddings(pri_fib_pyscenic_seu, "umap"), col = pyscenic_cluster_colors, pch = 16, cex = 0.5, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
plot(Embeddings(pri_fib_pyscenic_seu, "umap"), col = pyscenic_caf_cluster_colors, pch = 16, cex = 0.5, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
dev.off()
fib_sling_tf <- slingshot(Embeddings(pri_fib_pyscenic_seu, "umap"), clusterLabels = as.character(Idents(pri_fib_pyscenic_seu)), start.clus = 5, stretch = 0)
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun, categories)
    return(pal[cell_vars])
  }
}
cell_colors1 <- cell_pal(as.character(pri_fib_pyscenic_seu$SCT_snn_res.0.04), c("0"="#ef233c","1"="#ffba08","2"="#3f88c5","3"="#032b43"))
svg("path to store plot of UMAP based on TF, including pseudotime trajectory", height = 7, width = 7)
plot(Embeddings(pri_fib_pyscenic_seu, "umap"), col = cell_colors1, pch = 16, cex = 0.5, ylab = "", xlab = "", xaxt = "n", yaxt = "n")
lines(fib_sling_tf@metadata$curves$Lineage1, lwd = 2, col = 'black')
lines(fib_sling_tf@metadata$curves$Lineage2, lwd = 2, col = 'black')
lines(fib_sling_tf@metadata$curves$Lineage3, lwd = 2, col = 'black')
dev.off()

## ARACNe-AP
aracne_op <- read.table("path to aracne output", header = T, sep = "\t") %>% subset(pvalue < 1E-8)
gsea_aracne_caf0 <- GSEA(ranked_gopal_fib_genes$`0`, pAdjustMethod = "fdr", TERM2GENE = aracne_op[,1:2], pvalueCutoff = 1, minGSSize = 2, seed = 123, eps = 0)
gsea_aracne_caf1 <- GSEA(ranked_gopal_fib_genes$`1`, pAdjustMethod = "fdr", TERM2GENE = aracne_op[,1:2], pvalueCutoff = 1, minGSSize = 2, seed = 123, eps = 0)
gsea_aracne_caf2 <- GSEA(ranked_gopal_fib_genes$`2`, pAdjustMethod = "fdr", TERM2GENE = aracne_op[,1:2], pvalueCutoff = 1, minGSSize = 2, seed = 123, eps = 0)
gsea_aracne_caf3 <- GSEA(ranked_gopal_fib_genes$`3`, pAdjustMethod = "fdr", TERM2GENE = aracne_op[,1:2], pvalueCutoff = 1, minGSSize = 2, seed = 123, eps = 0)
nes_compare_df <- function(gsea_lst){
  # gsea_lst: named list of gsea outputs
  # gene sets must all be same
  result <- list(Gene=gsea_lst[[1]]@result$ID)
  for(i in 1:length(gsea_lst)){
    result[[names(gsea_lst)[i]]] <- gsea_lst[[i]]@result[result$Gene,"NES"]
  }
  for(i in 1:(length(gsea_lst)-1)){
    for(j in (i+1):length(gsea_lst)){
      curr_name <- paste(names(gsea_lst)[i],names(gsea_lst)[j],sep="_")
      result[[curr_name]] <- c()
      for(k in 1:length(result$Gene)){
        if(result[[names(gsea_lst)[i]]][k] > 0 & result[[names(gsea_lst)[j]]][k] > 0){
          result[[curr_name]] <- c(result[[curr_name]], "Both")
        } else if(result[[names(gsea_lst)[i]]][k] < 0 & result[[names(gsea_lst)[j]]][k] < 0){
          result[[curr_name]] <- c(result[[curr_name]], "Neither")
        } else if(result[[names(gsea_lst)[i]]][k] < 0 & result[[names(gsea_lst)[j]]][k] > 0){
          result[[curr_name]] <- c(result[[curr_name]], names(gsea_lst)[j])
        }else if(result[[names(gsea_lst)[i]]][k] > 0 & result[[names(gsea_lst)[j]]][k] < 0){
          result[[curr_name]] <- c(result[[curr_name]], names(gsea_lst)[i])
        }
      }
      result[[curr_name]] <- factor(result[[curr_name]], levels = c("Both", names(gsea_lst)[i],
                                                                    names(gsea_lst)[j], "Neither"))
    }
  }
  return(data.frame(result))
}
aracne_nes_df <- nes_compare_df(list("myCAF"=gsea_aracne_caf0, "IL-iCAF"=gsea_aracne_caf2,
                                     "IFN-iCAF"=gsea_aracne_caf1, "Adipose-like CAF"=gsea_aracne_caf3))
aracne_nes_df <- cbind(aracne_nes_df, data.frame("myCAF_TF"=rep("",nrow(aracne_nes_df)),
                                                 "IL.iCAF_TF"=rep(NA,nrow(aracne_nes_df)),
                                                 "IFN.iCAF_TF"=rep(NA,nrow(aracne_nes_df)),
                                                 "Adipose.like.CAF_TF"=rep(NA,nrow(aracne_nes_df))))
aracne_nes_df$myCAF_TF[which(aracne_nes_df$Gene %in% c("FOSL2","CREB3L1","TEAD1"))] <- aracne_nes_df$Gene[which(aracne_nes_df$Gene %in% c("FOSL2","CREB3L1","TEAD1"))]
aracne_nes_df$IL.iCAF_TF[which(aracne_nes_df$Gene %in% c("POU2F2","RELB"))] <- aracne_nes_df$Gene[which(aracne_nes_df$Gene %in% c("POU2F2","RELB"))]
aracne_nes_df$IFN.iCAF_TF[which(aracne_nes_df$Gene %in% c("NR2F2","IRF7"))] <- aracne_nes_df$Gene[which(aracne_nes_df$Gene %in% c("NR2F2","IRF7"))]
aracne_nes_df$Adipose.like.CAF_TF[which(aracne_nes_df$Gene %in% c("NR2F1"))] <- aracne_nes_df$Gene[which(aracne_nes_df$Gene %in% c("NR2F1"))]

# to visualize comparison of ARACNe-AP-predicted transcription factor activities
svg("desired filename", height = 7, width = 7)
ggplot(aracne_nes_df, aes(x=Adipose.like.CAF, y=IFN.iCAF, color = IFN.iCAF_Adipose.like.CAF)) +
  geom_point(size=2) + scale_color_manual(values=c("Adipose-like CAF"="#032b43","IFN-iCAF"="#ffba08","Both"="#817326","Neither"="black")) +
  theme(text=element_text(size=10,family="sans"), axis.text=element_text(size=10),
        axis.title=element_text(size=15),axis.title.y=element_text(size=15),
        panel.background = element_rect(fill = NA, color = "black", linewidth = 1),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("Adipose-like CAF") + ylab("IFN-iCAF") +
  ggrepel::geom_text_repel(label = aracne_nes_df$Adipose.like.CAF_TF, color = "black", size = 6) +
  ggrepel::geom_text_repel(label = aracne_nes_df$IFN.iCAF_TF, color = "black", size = 6) +
  geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') +
  annotate("text",label="Both\n4",x=0.1,y=0.3, color = "#817326", hjust = 0, size = 7) +
  annotate("text",label="48\nNeither",x=-0.1,y=-0.3, color = "black", hjust = 1, size = 7) +
  annotate("text",label="15\nAdipose-like CAF",x=0.1,y=-0.3, color = "#032b43", hjust = 0, size = 7) +
  annotate("text",label="IFN-iCAF\n45",x=-0.1,y=0.3, color = "#ffba08", hjust = 1, size = 7)
dev.off()
svg("desired filename", height = 7, width = 7)
ggplot(aracne_nes_df, aes(x=IFN.iCAF, y=myCAF, color = myCAF_IFN.iCAF)) +
  geom_point(size=2) + scale_color_manual(values=c("myCAF"="#ef233c","IFN-iCAF"="#ffba08","Both"="#f76f22","Neither"="black")) +
  theme(text=element_text(size=10,family="sans"), axis.text=element_text(size=10),
        axis.title=element_text(size=15),axis.title.y=element_text(size=15),
        panel.background = element_rect(fill = NA, color = "black", linewidth = 1),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("IFN-iCAF") + ylab("myCAF") +
  ggrepel::geom_text_repel(label = aracne_nes_df$myCAF_TF, color = "black", size = 6, max.overlaps = 50) +
  ggrepel::geom_text_repel(label = aracne_nes_df$IFN.iCAF_TF, color = "black", size = 6, max.overlaps = 50) +
  geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') +
  annotate("text",label="Both\n1",x=0.1,y=0.3, color = "#f76f22", hjust = 0, size = 7) +
  annotate("text",label="25\nNeither",x=-0.1,y=-0.3, color = "black", hjust = 1, size = 7) +
  annotate("text",label="myCAF\n38",x=-0.1,y=0.3, color = "#ef233c", hjust = 1, size = 7) +
  annotate("text",label="48\nIFN-iCAF",x=0.1,y=-0.3, color = "#ffba08", hjust = 0, size = 7)
dev.off()

## scriabin, refer to scriabin github for tutorial
scriabin_compare_patients <- function(scriabin_dai_lst, min_patients){
  all_int <- c()
  all_clus <- c()
  scriabin_dai_lst <- lapply(scriabin_dai_lst, function(x){subset(x, avg_log2FC > log2(1.2) & p_val_adj < 0.01)})
  result <- list(Interaction=c(),Cluster=c(),avg_avg_log2FC=c(),avg_padj=c())
  for(i in scriabin_dai_lst){
    all_int <- c(all_int, i$gene)
    all_clus <- c(all_clus, levels(i$cluster))
  }
  all_int <- unique(all_int)
  all_clus <- unique(all_clus)
  for(int in all_int){
    for(clus in all_clus){
      fc <- c()
      pval <- c()
      for(pt in scriabin_dai_lst){
        if(int %in% pt$gene){
          fc <- c(fc, pt$avg_log2FC[which(pt$gene == int & pt$cluster == clus)])
          pval <- c(pval, pt$p_val_adj[which(pt$gene == int & pt$cluster == clus)])
        }
      }
      if(length(fc)>=min_patients){
        result$Interaction <- c(result$Interaction, int)
        result$Cluster <- c(result$Cluster, clus)
        result$avg_avg_log2FC <- c(result$avg_avg_log2FC, mean(fc))
        result$avg_padj <- c(result$avg_padj, mean(pval))
      }
    }
  }
  return(data.frame(result) %>% .[order(.$avg_avg_log2FC, decreasing = T),] %>% .[order(.$Cluster),])
}
scriabin_conserved_caf_interactions <- scriabin_compare_patients(scriabin_all_diff_interactions, 3)
scriabin_conserved_caf_interactions[,c("Ligand","Receptor")] <- t(data.frame(strsplit(scriabin_conserved_caf_interactions$Interaction, "=")))

## for NATMI and CellPhoneDB
diff_edge <- function(op, lig_col, rec_col, sen_col, tar_col, val_col, tar_oi){
  # op: NATMI output
  # lig_col, rec_col, sen_col, tar_col, val_col: names of columns representing ligands, receptors, senders, receivers, and desired edge parameter (eg, weight, specificity, geometric mean of both), respectively
  op1 <- op[which(op[,tar_col] %in% tar_oi),]
  result <- list("Sender"=c(), "Ligand"=c(), "Receptor"=c(), "Target"=c(),
                 "Edge weight"=c(), "Rest edge weight"=c(), "Difference"=c())
  sen_lig_rec <- unique(paste(op1[,sen_col],op1[,lig_col],op1[,rec_col], sep = "__"))
  n_tar <- length(tar_oi)
  for(slr in sen_lig_rec){
    slr1 <- unlist(strsplit(slr, split = "__"))
    op2 <- op1[which(op1[,sen_col]==slr1[1] & op1[,lig_col]==slr1[2] & op1[,rec_col]==slr1[3]),]
    nrow_op2 <- nrow(op2)
    if(nrow_op2 > 0){
      result$Sender <- c(result$Sender, rep(slr1[1],nrow_op2))
      result$Ligand <- c(result$Ligand, rep(slr1[2],nrow_op2))
      result$Receptor <- c(result$Receptor, rep(slr1[3],nrow_op2))
      result$Target <- c(result$Target, as.character(op2[,tar_col]))
      for(tar in 1:nrow_op2){
        curr_weight <- as.numeric(op2[tar,val_col])
        if(nrow_op2>1){rest_weight <- mean(as.numeric(op2[-tar,val_col]))} else {rest_weight <- 0}
        result$`Edge weight` <- c(result$`Edge weight`, curr_weight)
        result$`Rest edge weight` <- c(result$`Rest edge weight`, rest_weight)
        result$Difference <- c(result$Difference, curr_weight-rest_weight)
      }
    }
  }
  result <- data.frame(result) %>% .[order(.$Difference, decreasing = T),] %>% .[order(.$Target),]
  return(result)
}
cpdb_for_diff_edge <- function(cpdb_op, cell_types, ignore_first_ncol){
  # to convert CellPhoneDB output suitable for diff_edge
  # cpdb_op: CellPhoneDB output
  # cell_types: cell types present in dataset, in order of appearance
  # ignore_first_ncol: number of columns to ignore such that only the edge weight data are included
  process_cpdb <- function(cpdb, no_cell_type, ignore_first_col){
    secreted_swap <- which(cpdb$receptor_a == "True" & cpdb$receptor_b == "False")
    integrin_swap <- which(cpdb$is_integrin == "True")
    for(swap in unique(c(secreted_swap, integrin_swap))){
      if(cpdb$is_integrin[swap] == "True"){
        if(grepl("complex", cpdb$partner_a[swap])){
          cpdb[swap,c("partner_a", "partner_b", "gene_a", "gene_b")] <- cpdb[swap,c("partner_b", "partner_a", "gene_b", "gene_a")]
          cpdb[swap,c("receptor_a", "receptor_b")] <- c("False", "True")
          to_swap <- T
        } else {
          cpdb[swap,c("receptor_a", "receptor_b")] <- c("False", "True")
          to_swap <- F
        }
      } else {
        cpdb[swap,c("partner_a", "partner_b", "gene_a", "gene_b")] <- cpdb[swap,c("partner_b", "partner_a", "gene_b", "gene_a")]
        cpdb[swap,c("receptor_a", "receptor_b")] <- c("False", "True")
        to_swap <- T
      }
      if(to_swap){
        for(i in 1:(no_cell_type-1)){
          for(j in (i+1):no_cell_type){
            cpdb[swap, c(((i-1)*no_cell_type)+j+ignore_first_col, ((j-1)*no_cell_type)+i+ignore_first_col)] <- cpdb[swap, c(((j-1)*no_cell_type)+i+ignore_first_col, ((i-1)*no_cell_type)+j+ignore_first_col)]
          }
        }
      }
    }
    return(cpdb)
  }
  cpdb_op <- process_cpdb(cpdb_op, length(cell_types), ignore_first_ncol)
  result <- list("Sender"=c(), "Ligand"=c(), "Receptor"=c(), "Target"=c(),
                 "Edge mean"=c())
  sen_tar <- c()
  for(sen in cell_types){
    for(tar in cell_types){
      sen_tar <- c(sen_tar, paste(sen, tar, sep = "__"))
    }
  }
  no_sen_tar <- length(sen_tar)
  sen_tar <- strsplit(sen_tar, split = "__")
  sens <- unlist(lapply(sen_tar, function(x){x[1]}))
  tars <- unlist(lapply(sen_tar, function(x){x[2]}))
  cpdb_op1 <- cpdb_op[,(ignore_first_ncol+1):ncol(cpdb_op)]
  for(lr in 1:nrow(cpdb_op)){
    if(cpdb_op$gene_a[lr]==""){lig <- strsplit(cpdb_op$partner_a[lr], split = ":") %>% unlist() %>% .[2]} else {lig <- cpdb_op$gene_a[lr]}
    if(cpdb_op$gene_b[lr]==""){rec <- strsplit(cpdb_op$partner_b[lr], split = ":") %>% unlist() %>% .[2]} else {rec <- cpdb_op$gene_b[lr]}
    non_zero <- which(!is.na(cpdb_op1[lr,]))
    no_non_zero <- length(non_zero)
    result$Ligand <- c(result$Ligand, rep(lig, no_non_zero))
    result$Receptor <- c(result$Receptor, rep(rec, no_non_zero))
    result$Sender <- c(result$Sender, sens[non_zero])
    result$Target <- c(result$Target, tars[non_zero])
    result$`Edge mean` <- c(result$`Edge mean`, as.numeric(cpdb_op1[lr,non_zero]))
  }
  return(data.frame(result))
}

pri_natmi <- read.csv("path to NATMI results; used excel to calculate geometric mean of edge weight and specificity")
pri_natmi <- pri_natmi[-which(pri_natmi$Receptor.symbol=="COL13A1"),]
pri_natmi$Average.edge.weight.and.specificity <- sqrt(pri_natmi$Edge.average.expression.derived.specificity * pri_natmi$Edge.average.expression.weight)
caf_natmi <- subset(pri_natmi, grepl("CAF", Target.cluster)) %>% .[.$Average.edge.weight.and.specificity > quantile(.$Average.edge.weight.and.specificity, probs = 0.8),]
caf_natmi <- diff_edge(pri_natmi, "Ligand.symbol", "Receptor.symbol", "Sending.cluster", "Target.cluster",
                       "Average.edge.weight.and.specificity", c("Adipose-like CAF", "myCAF", "IFN-iCAF","IL-iCAF"))
pri_cpdb <- read.table("path to CellPhoneDB output", sep = "\t", header = T)
pri_cpdb1 <- cpdb_for_diff_edge(pri_cpdb, c("Adipose-like CAF","B cell","CD4 T cell","CD8 T cell","Endothelial",
                                            "IFN-iCAF","IL-iCAF","Malignant","Mast cell","Monocyte","NK cell","NKT cell",
                                            "Naive/Memory T cell","Neutrophil","Pericyte","Plasma cell",
                                            "TAM","Treg","myCAF","Myeloid DC",
                                            "pDC"), 12)
# pri_cpdb1$Sender[which(pri_cpdb1$Sender == "myeloid DC")] <- rep("Myeloid DC", length(which(pri_cpdb1$Sender == "myeloid DC")))
# pri_cpdb1$Target[which(pri_cpdb1$Target == "myeloid DC")] <- rep("Myeloid DC", length(which(pri_cpdb1$Target == "myeloid DC")))
caf_cpdb <- diff_edge(pri_cpdb1, "Ligand", "Receptor", "Sender", "Target", "Edge.mean",
                      c("myCAF", "Adipose-like CAF", "IL-iCAF", "IFN-iCAF"))

circos_for_lr <- function(diff_edge_op, target_clus, top_n_lr, sender_colours, target_colour){
  # diff_edge_op: output of diff_edge
  # target_clus: which receiving cell type to focus on
  # top_n_lr: top percentage of interactions to keep
  # sender_colours: color to represent each sender cell type
  # target_colour: color to represent each receiver cell type
  diff_edge_op1 <- diff_edge_op[order(diff_edge_op$Difference, decreasing = T),] %>% subset(Target==target_clus)
  lig_rec <- paste(diff_edge_op1$Ligand, diff_edge_op1$Receptor, sep = "__")
  lig_rec <- which(!duplicated(lig_rec)) %>% .[1:min(length(.), top_n_lr)]
  diff_edge_op1 <- diff_edge_op1[lig_rec,]
  sen_col <- unlist(sapply(diff_edge_op1$Sender, function(x){sender_colours[x]}))
  names(sen_col) <- diff_edge_op1$Ligand
  tar_col <- rep(target_colour, nrow(diff_edge_op1))
  names(tar_col) <- diff_edge_op1$Receptor
  diff_edge_op1$Sender <- factor(diff_edge_op1$Sender, levels = names(sender_colours))
  diff_edge_op1 <- diff_edge_op1[order(diff_edge_op1$Sender),]
  grouping <- as.character(diff_edge_op1$Ligand)
  names(grouping) <- as.character(diff_edge_op1$Sender)
  diff_edge_op1 <- diff_edge_op1[order(diff_edge_op1$Difference, decreasing = T),] %>% .[order(.$Target),]
  result <- list("df"=diff_edge_op1,"col"=c(sen_col,tar_col))
  for(i in unique(names(grouping))){
    result[[i]] <- grouping[which(names(grouping)==i)]
  }
  return(result)
}
caf_colours <- c("myCAF"="#ef233c","IL-iCAF"="#3f88c5","IFN-iCAF"="#ffba08","Adipose-like CAF"="#032b43")
cell_type_colours <- c("Malignant"="#FF0000","myCAF"="#ef233c","IL-iCAF"="#3f88c5","IFN-iCAF"="#ffba08","Adipose-like CAF"="#032b43",
                       "Pericyte"="#16B251","Endothelial"="#F74ED6",
                       "CD8 T cell"="#00FE7F","B cell"="#8C56F8","Plasma cell"="#A6752C",
                       "NK cell"="#0066FF","TAM"="#ffd500","Myeloid DC"="#E19FD3","NKT cell"="#7209b7",
                       "pDC"="#a6a2a2","Monocyte"="#9A4452",
                       "Neutrophil"="#50D2FA","Mast cell"="#023047", "CD4 T cell"="#ee6c4d",
                       "Naive/Memory T cell"="black")
ilicaf_cpdb_circos <- circos_for_lr(caf_cpdb, "IL-iCAF", 20, cell_type_colours, caf_colours["IL-iCAF"])
ifnicaf_cpdb_circos <- circos_for_lr(caf_cpdb, "IFN-iCAF", 20, cell_type_colours, caf_colours["IFN-iCAF"])
adipocaf_cpdb_circos <- circos_for_lr(caf_cpdb, "Adipose-like CAF", 20, cell_type_colours, caf_colours["Adipose-like CAF"])
mycaf_cpdb_circos <- circos_for_lr(caf_cpdb, "myCAF", 20, cell_type_colours, caf_colours["myCAF"])
ilicaf_natmi_circos <- circos_for_lr(caf_natmi, "IL-iCAF", 20, cell_type_colours, caf_colours["IL-iCAF"])
ifnicaf_natmi_circos <- circos_for_lr(caf_natmi, "IFN-iCAF", 20, cell_type_colours, caf_colours["IFN-iCAF"])
adipocaf_natmi_circos <- circos_for_lr(caf_natmi, "Adipose-like CAF", 20, cell_type_colours, caf_colours["Adipose-like CAF"])
mycaf_natmi_circos <- circos_for_lr(caf_natmi, "myCAF", 20, cell_type_colours, caf_colours["myCAF"])

cpdb_split_complex <- function(circos_for_lr_op, recep_annot_lst, ligand_annot_lst=NULL){
  # circos_for_lr_op: output of circos_for_lr
  # recep_annot_lst, ligand_annot_lst: list of complex receptors or ligands, respectively, that need to be split into individual subunits
  split_complex <- function(df, column2change, annot){
    to_change <- which(df[,column2change] == names(annot))
    to_bind <- df[to_change,]
    no_row <- nrow(to_bind)
    if(length(annot[[1]])>1){
      for(i in 1:(length(annot[[1]])-1)){
        to_bind <- rbind(to_bind, to_bind)
      }
    }
    to_bind[,column2change] <- rep(annot[[1]], each=no_row)
    df <- rbind(df[-to_change,],to_bind) %>% .[order(.$Difference, decreasing=T),] %>% .[order(.$Target),]
    return(df)
  }
  for(com in 1:length(recep_annot_lst)){
    circos_for_lr_op$df <- split_complex(circos_for_lr_op$df, "Receptor", recep_annot_lst[com])
  }
  if(length(ligand_annot_lst)>0){
    for(com in 1:length(ligand_annot_lst)){
      circos_for_lr_op$df <- split_complex(circos_for_lr_op$df, "Ligand", ligand_annot_lst[com])
    }
  }
  return(circos_for_lr_op)
}
ilicaf_cpdb_circos <- cpdb_split_complex(ilicaf_cpdb_circos, recep_annot_lst = list("integrin_aVb3_complex"=c("ITGAV","ITGB3")))
ifnicaf_cpdb_circos <- cpdb_split_complex(ifnicaf_cpdb_circos, recep_annot_lst = list("integrin_a1b1_complex"=c("ITGA1","ITGB1")), ligand_annot_lst = list("Adenosine_byNT5E_and_SLC29A1"=c("NT5E","SLC29A1")))
mycaf_cpdb_circos <- cpdb_split_complex(mycaf_cpdb_circos, recep_annot_lst = list("integrin_a11b1_complex"=c("ITGA11","ITGB1"),
                                                                                  "integrin_aVb1_complex"=c("ITGAV","ITGB1"),
                                                                                  "integrin_a4b1_complex"=c("ITGA4","ITGB1"),
                                                                                  "integrin_aVb5_complex"=c("ITGAV","ITGB5"),
                                                                                  "integrin_a5b1_complex"=c("ITGA5","ITGB1")),
                                        ligand_annot_lst = list("ProstaglandinE2_byPTGES3"="PTGES3"))
adipocaf_cpdb_circos <- cpdb_split_complex(adipocaf_cpdb_circos, recep_annot_lst = list("IGF1R_enhancerComGPC3"=c("IGFR1","GPC3")),
                                           ligand_annot_lst = list("ProstaglandinF2a_byAKR1B1"="AKR1B1",
                                                                   "ProstaglandinF2a_byCBR1"="CBR1"))
caf_cpdb_df <- rbind(ilicaf_cpdb_circos$df,ifnicaf_cpdb_circos$df,adipocaf_cpdb_circos$df,mycaf_cpdb_circos$df)
caf_natmi_df <- rbind(ilicaf_natmi_circos$df,ifnicaf_natmi_circos$df,adipocaf_natmi_circos$df,mycaf_natmi_circos$df)

combine_cpdb_natmi_scriabin <- function(cpdb_circos_df, natmi_circos_df, scriabin_dai, targets){
  # cpdb_circos_df, natmi_circos_df: output of circos_for_lr from CellPhoneDB and NATMI analyses, respectively
  # scriabin_dai: differential ligand-receptor interaction analysis for Scriabin analysis
  # targets: receiving cell types of interest
  lig <- unique(c(cpdb_circos_df$Ligand, natmi_circos_df$Ligand, scriabin_dai$Ligand)) %>% .[order(.)]
  rec <- unique(c(cpdb_circos_df$Receptor, natmi_circos_df$Receptor, scriabin_dai$Receptor)) %>% .[order(.)]
  result <- matrix("No",nrow=length(lig)*length(targets)/2,ncol=length(rec)*length(targets)/2)
  row.names(result) <- paste(rep(lig,each=length(targets)/2),targets[1:(length(targets)/2)],sep="_")
  colnames(result) <- paste(rep(rec,each=length(targets)/2),targets[((length(targets)/2)+1):length(targets)],sep="_")
  for(l in lig){
    for(r in rec){
      count <- rep(0,length(targets))
      names(count) <- targets
      tmp <- subset(cpdb_circos_df, Ligand == l & Receptor == r)
      if(nrow(tmp)>0){
        count[tmp$Target] <- count[tmp$Target] + 1
      }
      tmp <- subset(natmi_circos_df, Ligand == l & Receptor == r)
      if(nrow(tmp)>0){
        count[tmp$Target] <- count[tmp$Target] + 1
      }
      tmp <- subset(scriabin_dai, Ligand == l & Receptor == r)
      if(nrow(tmp)>0){
        count[tmp$Cluster] <- count[tmp$Cluster] + 1
      }
      for(i in 1:length(count)){
        if(count[i]>1 & i<=(length(count)/2)){
          result[paste(l,names(count)[i],sep="_"),paste(r,names(count)[3],sep="_")] <- paste(names(count)[i],as.character(count[i]),sep="_")
        } else if(count[i]>1 & i>(length(count)/2)){
          result[paste(l,names(count)[i-2],sep="_"),paste(r,names(count)[4],sep="_")] <- paste(names(count)[i],as.character(count[i]),sep="_")
        }
      }
    }
  }
  rem_row <- c()
  rem_col <- c()
  for(i in seq(1,nrow(result),2)){
    if(length(unique(as.character(result[i:(i+1),])))==1){rem_row <- c(rem_row, i, i+1)}
  }
  result <- result[-rem_row,]
  for(i in seq(1,ncol(result),2)){
    if(length(unique(as.character(result[,i:(i+1)])))==1){rem_col <- c(rem_col, i, i+1)}
  }
  result <- result[,-rem_col]
  # rem_row <- apply(result,1,unique) %>% lapply(length) %>% unlist()==1
  # rem_col <- apply(result,2,unique) %>% lapply(length) %>% unlist()==1
  # return(result[-which(rem_row),-which(rem_col)])
  return(result)
}
lr_mat <- combine_cpdb_natmi_scriabin(caf_cpdb_df, caf_natmi_df, scriabin_conserved_caf_interactions,
                                      c("myCAF","IL-iCAF","IFN-iCAF","Adipose-like CAF"))

## get CAF state-specific markers for subsequent TCGA analyses
mycaf_v_all_deg <- FindMarkers(pri1, ident.1 = "myCAF")
ilcaf_v_all_deg <- FindMarkers(pri1, ident.1 = "IL-iCAF")
ifncaf_v_all_deg <- FindMarkers(pri1, ident.1 = "IFN-iCAF")
adipocaf_v_all_deg <- FindMarkers(pri1, ident.1 = "Adipose-like CAF")
caf_v_all_deg <- cbind(rbind(mycaf_v_all_deg,ilcaf_v_all_deg,ifncaf_v_all_deg,adipocaf_v_all_deg),
                       data.frame(cluster=c(rep("myCAF",nrow(mycaf_v_all_deg)),
                                            rep("IL-iCAF",nrow(ilcaf_v_all_deg)),
                                            rep("IFN-iCAF",nrow(ifncaf_v_all_deg)),
                                            rep("Adipose-like CAF",nrow(adipocaf_v_all_deg))),
                                  gene=c(row.names(mycaf_v_all_deg),
                                         row.names(ilcaf_v_all_deg),
                                         row.names(ifncaf_v_all_deg),
                                         row.names(adipocaf_v_all_deg)))) %>%
  .[order(.$avg_log2FC, decreasing=T),] %>% .[order(.$cluster),]
fib_v_all_top_deg <- get_top_markers(caf_v_all_deg, log2fc_thresh = 0.25)

# bulk RNA seq
create_counts_df <- function(dir, start_row, counts_col, counts_file_id){
  # dir: path to directory containing counts output from STAR
  # start_row: row where data starts
  # counts_col: column to use (depends on library prep method)
  # counts_file_id: include only files whose names contain this pattern
  filenames <- list.files(dir)
  path <- paste(dir, "\\", sep = "")
  sample_id <- c()
  x <- 0
  for(i in 1:length(filenames)){
    if(grepl(counts_file_id, filenames[i])){
      x <- x + 1
      file <- read.table(paste(path, filenames[i], sep = "")) %>% .[start_row:nrow(.),]
      genes <- file[,1]
      file <- as.data.frame(matrix(as.numeric(file[,counts_col]), ncol = 1))
      row.names(file) <- genes
      i_split <- strsplit(filenames[i], split = "_")
      sample_id <- c(sample_id, paste(i_split[[1]][1], i_split[[1]][2], sep = "_"))
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
  colnames(result) <- sample_id
  no_exp <- which(as.numeric(rowSums(result))==0)
  if(length(no_exp)>0){
    result <- result[-which(rowSums(result)==0),]
  }
  return(result)
}

combined_cafs <- create_counts_df(filename, 5, 4, "ReadsPerGene") # replace filename with path to directory containing counts
combined_cafs_meta <- data.frame(condition = factor(c("RGD", rep(c("2D", "RDG", "isoPQ", "RGD"), 3), rep(c("2D", "isoPQ", "RDG", "RGD"), 6)), levels = c("RGD", "RDG", "isoPQ", "2D")),
                                 patient = factor(c(rep("HN217P", 13), rep("HN294P", 12), rep("HN338P", 12)), levels = c("HN217P", "HN294P","HN338P")),
                                 passage = c(rep(7,5), rep(8,4), rep(9,4), rep(10,4), rep(8,4), rep(9,4), rep(7,4), rep(8,4), rep(9,4)))
combined_cafs_cpm <- t(t(combined_cafs*1e6)/colSums(combined_cafs))
hist(log2(unlist(combined_cafs_cpm)), breaks = 100)
combined_cafs_cpm_filt <- combined_cafs_cpm[which(rowMeans(combined_cafs_cpm)>4),]

# correct counts for batch effects for use in visualization and deconvolution
combined_cafs_counts_adj <- sva::ComBat_seq(as.matrix(combined_cafs[row.names(combined_cafs_cpm_filt),]), batch = c(rep(1,13), rep(2,12), rep(3,12)))
combined_cafs_cpm_adj <- t(t(combined_cafs_counts_adj*1e6)/colSums(combined_cafs_counts_adj))

# CibersortX parameters
## Signature matrix: disable quantile normalization = false, min expression = 0.1, sampling = 0.9
## Impute cell fractions: batch correction= S-mode, disable quantile normalization = true, permuatations = 100
ciber <- read.csv("path to CibersortX output")
ciber$Condition <- combined_cafs_meta$condition %>% factor(levels=c("2D","RGD","RDG","isoPQ"))
ciber$Patient <- combined_cafs_meta$patient
ciber <- ciber[order(ciber$Condition),]
ciber4gg <- function(df, prop_columns, condition_column){
  no_prop <- length(prop_columns)
  result <- data.frame(Mixture=rep(df$Mixture,each=no_prop) %>% factor(levels=df$Mixture),
                       Condition=rep(df[,condition_column],each=no_prop),
                       Type=rep(prop_columns,nrow(df)),
                       Proportion=unlist(as.data.frame(t(df[,prop_columns]))))
  return(result)
}
ciber_gg <- ciber4gg(ciber,colnames(ciber)[2:5],"Condition")
ciber_gg$x <- c(rep(1:9,each=4),rep(11:20,each=4),rep(22:30,each=4),rep(32:40,each=4))
caf_colours1 <- caf_colours
names(caf_colours1) <- colnames(ciber)[2:5]
# to visualize proportion of CAF states in in vitro samples
svg("desired filename",
    height = 5, width=10)
ggplot(ciber_gg, aes(x = x, y = Proportion, fill = Type)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=caf_colours1,name="State") + xlab("Sample") +
  theme(text=element_text(size=10,family="sans"), axis.text.x=element_blank(), axis.text.y=element_text(size=10),
        legend.text=element_text(size=10), legend.title=element_text(size=10), axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15), panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.ticks.x = element_blank())
dev.off()

combined_cafs_meta_deseq <- data.frame(patient = factor(c(rep("HN217P", 13), rep("HN294P", 12), rep("HN338P", 12)), levels = c("HN217P", "HN294P", "HN338P")),
                                       twoD = factor(c(0,rep(c(1,0,0,0),9)), levels = c(0,1)),
                                       RGD.PQ = factor(c(1,rep(c(0,0,0,1),9)), levels = c(0,1)),
                                       RDG.PQ = factor(c(0,rep(c(0,1,0,0),3),rep(c(0,0,1,0),6)), levels = c(0,1)),
                                       RGD.isoPQ = factor(c(0,rep(c(0,0,1,0),3),rep(c(0,1,0,0),6)), levels = c(0,1)))
row.names(combined_cafs_meta_deseq) <- colnames(combined_cafs_cpm_filt)

#exclude 2D for DEG analysis
dds_rgd <- DESeq2::DESeqDataSetFromMatrix(countData = combined_cafs[row.names(combined_cafs_cpm_filt),which(combined_cafs_meta_deseq$twoD==0)],
                                          colData = combined_cafs_meta_deseq[which(combined_cafs_meta_deseq$twoD==0),],
                                          design= ~ patient + RGD.PQ)
dds_rgd <- DESeq2::DESeq(dds_rgd)
res_rgd <- as.data.frame(DESeq2::lfcShrink(dds_rgd, coef="RGD.PQ_1_vs_0", type="apeglm")) %>% .[order(.$log2FoldChange, decreasing = T),]
res_rgd$gene <- row.names(res_rgd)

dds_rdg <- DESeq2::DESeqDataSetFromMatrix(countData = combined_cafs[row.names(combined_cafs_cpm_filt),which(combined_cafs_meta_deseq$twoD==0)],
                                          colData = combined_cafs_meta_deseq[which(combined_cafs_meta_deseq$twoD==0),],
                                          design= ~ patient + RDG.PQ)
dds_rdg <- DESeq2::DESeq(dds_rdg)
res_rdg <- as.data.frame(DESeq2::lfcShrink(dds_rdg, coef="RDG.PQ_1_vs_0", type="apeglm")) %>% .[order(.$log2FoldChange, decreasing = T),]
res_rdg$gene <- row.names(res_rdg)

dds_iso <- DESeq2::DESeqDataSetFromMatrix(countData = combined_cafs[row.names(combined_cafs_cpm_filt),which(combined_cafs_meta_deseq$twoD==0)],
                                          colData = combined_cafs_meta_deseq[which(combined_cafs_meta_deseq$twoD==0),],
                                          design= ~ patient + RGD.isoPQ)
dds_iso <- DESeq2::DESeq(dds_iso)
res_iso <- as.data.frame(DESeq2::lfcShrink(dds_iso, coef="RGD.isoPQ_1_vs_0", type="apeglm")) %>% .[order(.$log2FoldChange, decreasing = T),]
res_iso$gene <- row.names(res_iso)

ranked_rgd_deg <- res_rgd$log2FoldChange
names(ranked_rgd_deg) <- res_rgd$gene
ranked_rdg_deg <- res_rdg$log2FoldChange
names(ranked_rdg_deg) <- res_rdg$gene
ranked_iso_deg <- res_iso$log2FoldChange
names(ranked_iso_deg) <- res_iso$gene

gsea_gopalcafs_rgd <- GSEA(ranked_rgd_deg, pAdjustMethod = "fdr", TERM2GENE = gopal_fib_markers[,2:3], seed = 123, pvalueCutoff = 1, eps = 0)
gsea_gopalcafs_rdg <- GSEA(ranked_rdg_deg, pAdjustMethod = "fdr", TERM2GENE = gopal_fib_markers[,2:3], seed = 123, pvalueCutoff = 1, eps = 0)
gsea_gopalcafs_iso <- GSEA(ranked_iso_deg, pAdjustMethod = "fdr", TERM2GENE = gopal_fib_markers[,2:3], seed = 123, pvalueCutoff = 1, eps = 0)

pool_samples <- function(mtx, cluster, pool_fn){
  # mtx: matrix of genes in rows and samples in columns
  # cluster: vector indicating how samples should be grouped
  # pool_fn: function such as mean or median to find average of expression per cluster
  result <- list()
  for(i in unique(cluster)){
    clus_id <- which(cluster==i)
    result[[i]] <- as.vector(apply(mtx[,clus_id], 1, pool_fn))
  }
  result <- as.data.frame(result)
  row.names(result) <- row.names(mtx)
  return(result)
}
pooled_combined_caf_cpm <- pool_samples(combined_cafs_cpm_filt_adj[,-seq(2,37,4)], c("RGD", rep(c("RDG", "isoPQ", "RGD"), 3), rep(c("isoPQ", "RDG", "RGD"), 6)), mean)

gsea_dotplot <- function(lst, group_vec, rs_fn, split){
  clus <- strsplit(lst[[1]][,"ID"], split = split) %>% lapply(function(x){if(length(x)==2){x[1]}else{paste(x[1],x[2])}}) %>% unlist() %>% unique()
  result <- list(Cluster=rep(clus,each=length(lst)),Condition=rep(names(lst),length(clus)),NES=c(),"-log10(FDR)"=c())
  for(i in clus){
    for(j in names(lst)){
      tmp <- lst[[j]][which(grepl(i,lst[[j]][,"ID"])),]
      id <- sapply(group_vec, function(x){which(grepl(x,tmp$ID))})
      result$NES <- c(result$NES, as.numeric(tmp[id[1],"NES"]))
      for(k in 2:length(id)){
        result$NES[length(result$NES)] <- rs_fn(result$NES[length(result$NES)], tmp[id[k],"NES"])
      }
      result$`-log10(FDR)` <- c(result$`-log10(FDR)`, mean(sum(-log10(tmp$qvalue))))
    }
  }
  result <- data.frame(result)
  result$Cluster <- paste(result$Cluster, rep("Score",nrow(result)), sep = " ")
  result <- result[order(result$Cluster),]
  return(result)
}
vitro_vivo_gsea <- gsea_dotplot(lst=list("Integrin-sensitive and Degradable"=gsea_gopalcafs_rgd@result,
                                         "- Integrin-sensitive"=gsea_gopalcafs_rdg@result,
                                         "- Degradable"=gsea_gopalcafs_iso@result),
                                group_vec=c("up", "down"), rs_fn=function(x,y){return(x-y)},
                                split=" ")
vitro_vivo_gsea$Condition <- factor(vitro_vivo_gsea$Condition, levels = c("Integrin-sensitive and Degradable","- Degradable","- Integrin-sensitive"))

# to visualize how in vitro samples compared to in vivo
svg(filename, height = 15, width = 15, bg = NA) # replace filename with desired path
circos.clear()
circos.par("track.height" = 0.7, gap.after = 10)
circos.initialize(sectors = vitro_vivo_gsea$Cluster, x = vitro_vivo_gsea$Condition)
circos.track(sectors = vitro_vivo_gsea$Cluster, ylim = c(-5,5), bg.col = caf_colours, bg.lwd = 0.1, bg.border = "white", panel.fun = function(x,y) {
  sector.index = get.cell.meta.data("sector.index")
  xcenter = get.cell.meta.data("xcenter")
  ycenter = get.cell.meta.data("ycenter")
  # circos.text(xcenter, ycenter, sector.index,
  #             cex = 2.5,
  #             niceFacing =T,
  #             facing="bending.inside", adj=c(0,5))
  circos.yaxis(side = "left", at = seq(-4,4,2), labels = seq(-4,4,2), col = "black", labels.cex = 2.5)
  circos.barplot(vitro_vivo_gsea[which(vitro_vivo_gsea$Cluster==sector.index),"NES"], seq(1.4,2.6,0.6), col = c("#d00000","#eca400","#0d3b66"), sector.index = sector.index, track.index = 1, lwd = 2, border = "white")
  # circos.xaxis(h="bottom",labels=1:3,major.tick=F,minor.tick=F,major.at=F,col="black",direction = "inside", labels.facing = "inside")
})
circos.axis(h = "bottom", labels = 1:3, major.at = c(1.35,2,2.65), sector.index = "myCAF Score", direction = "inside", labels.facing = "inside", labels.niceFacing = T, major.tick = F, col = NA, labels.col = "#000000", labels.cex = 2.5)
circos.axis(h = "bottom", labels = 1:3, major.at = c(1.35,2,2.65), sector.index = "IL-iCAF Score", direction = "inside", labels.facing = "inside", labels.niceFacing = T, major.tick = F, col = NA, labels.col = "#000000", labels.cex = 2.5)
circos.axis(h = "bottom", labels = 1:3, major.at = c(1.35,2,2.65), sector.index = "IFN-iCAF Score", direction = "inside", labels.facing = "inside", labels.niceFacing = T, major.tick = F, col = NA, labels.col = "#000000", labels.cex = 2.5)
circos.axis(h = "bottom", labels = 1:3, major.at = c(1.35,2,2.65), sector.index = "Adipose-like CAF Score", direction = "inside", labels.facing = "inside", labels.niceFacing = T, major.tick = F, col = NA, labels.col = "#000000", labels.cex = 2.5)
circos.clear()
dev.off()

gsea_aracne_rgd <- GSEA(ranked_rgd_deg, pAdjustMethod = "fdr", TERM2GENE = aracne_bulk_op[,1:2], pvalueCutoff = 1, minGSSize = 2, seed = 123, eps = 0)
gsea_aracne_rdg <- GSEA(ranked_rdg_deg, pAdjustMethod = "fdr", TERM2GENE = aracne_bulk_op[,1:2], pvalueCutoff = 1, minGSSize = 2, seed = 123, eps = 0)
gsea_aracne_iso <- GSEA(ranked_iso_deg, pAdjustMethod = "fdr", TERM2GENE = aracne_bulk_op[,1:2], pvalueCutoff = 1, minGSSize = 2, seed = 123, eps = 0)
gsea_aracne_rgdv0 <- data.frame(myCAF=gsea_aracne_caf0@result[intersect(gsea_aracne_caf0@result$ID,gsea_aracne_rgd@result$ID),"NES"],
                                RGD=gsea_aracne_rgd@result[intersect(gsea_aracne_caf0@result$ID,gsea_aracne_rgd@result$ID),"NES"],
                                common=intersect(gsea_aracne_caf0@result$ID,gsea_aracne_rgd@result$ID) %>% sapply(function(x){if(x %in% c("CREB3L1","TEAD1","FOSL2")){x}else{""}}))
gsea_aracne_isov2 <- data.frame(ILiCAF=gsea_aracne_caf2@result[intersect(gsea_aracne_caf2@result$ID,gsea_aracne_iso@result$ID),"NES"],
                                isoPQ=gsea_aracne_iso@result[intersect(gsea_aracne_caf2@result$ID,gsea_aracne_iso@result$ID),"NES"],
                                common=intersect(gsea_aracne_caf2@result$ID,gsea_aracne_iso@result$ID) %>% sapply(function(x){if(x %in% c("RELB","POU2F2")){x}else{""}}))
svg(filename, # replace filename with desired path
    height=5, width=5)
ggplot(gsea_aracne_rgdv0, aes(myCAF,RGD)) + geom_point() +
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid = element_blank(),
        axis.text = element_text(size=12), axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.5), axis.title = element_text(size=15)) +
  geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") +
  annotate("text",label="26.3%",x=0.1,y=0.3, color = "#ff0022", hjust = 0, size = 7) +
  annotate("text",label="43.4%",x=-0.1,y=-0.3, color = "#ff0022", hjust = 1, size = 7) +
  annotate("text",label="20.2%",x=-0.1,y=0.3, color = "#011627", hjust = 1, size = 7) +
  annotate("text",label="10.1%",x=0.1,y=-0.3, color = "#011627", hjust = 0, size = 7) +
  ggrepel::geom_text_repel(label = gsea_aracne_rgdv0$common, color = "black", size = 6) +
  xlab("In vivo myCAF") + ylab("In vitro myCAF")
dev.off()
svg(filename, # replace filename with desired path
    height=5, width=5)
ggplot(gsea_aracne_isov2, aes(ILiCAF,isoPQ)) + geom_point() +
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid = element_blank(),
        axis.text = element_text(size=12), axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.5), axis.title = element_text(size=15)) +
  geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") +
  annotate("text",label="28.3%",x=0.1,y=0.3, color = "#ff0022", hjust = 0, size = 7) +
  annotate("text",label="26.3%",x=-0.1,y=-0.3, color = "#ff0022", hjust = 1, size = 7) +
  annotate("text",label="17.2%",x=-0.1,y=0.3, color = "#011627", hjust = 1, size = 7) +
  annotate("text",label="28.3%",x=0.1,y=-0.3, color = "#011627", hjust = 0, size = 7) +
  ggrepel::geom_text_repel(label = gsea_aracne_isov2$common, color = "black", size = 6) +
  xlab("In vivo IL-iCAF") + ylab("In vitro IL-iCAF")
dev.off()

#for bulk aracne
write.table(log1p(combined_cafs_counts_adj[row.names(combined_cafs_cpm_filt_adj),-seq(2,37,4)]), quote = F, sep = "\t",
            file = "C:/Users/nicho/ARACNe-AP-master/230710_caf_bulk/combined_cafs_counts_adj_filt.txt")
aracne_bulk_op <- read.table("C:/Users/nicho/ARACNe-AP-master/230710_caf_bulk/output/network.txt", header = T)

## to visualize comparison of in vivo and in vitro CAF states based on ARACNe-AP-predicted transcription factor activities
svg("desired filename",
    height=5, width=5)
ggplot(gsea_aracne_rgdv0, aes(myCAF,RGD)) + geom_point() +
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid = element_blank(),
        axis.text = element_text(size=12), axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.5), axis.title = element_text(size=15)) +
  geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") +
  annotate("text",label="26.3%",x=0.1,y=0.3, color = "#ff0022", hjust = 0, size = 7) +
  annotate("text",label="43.4%",x=-0.1,y=-0.3, color = "#ff0022", hjust = 1, size = 7) +
  annotate("text",label="20.2%",x=-0.1,y=0.3, color = "#011627", hjust = 1, size = 7) +
  annotate("text",label="10.1%",x=0.1,y=-0.3, color = "#011627", hjust = 0, size = 7) +
  ggrepel::geom_text_repel(label = gsea_aracne_rgdv0$common, color = "black", size = 6) +
  xlab("In vivo myCAF") + ylab("In vitro myCAF")
dev.off()
svg("desired filename",
    height=5, width=5)
ggplot(gsea_aracne_isov2, aes(ILiCAF,isoPQ)) + geom_point() +
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid = element_blank(),
        axis.text = element_text(size=12), axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.5), axis.title = element_text(size=15)) +
  geom_hline(yintercept=0,linetype="dashed") + geom_vline(xintercept=0,linetype="dashed") +
  annotate("text",label="28.3%",x=0.1,y=0.3, color = "#ff0022", hjust = 0, size = 7) +
  annotate("text",label="26.3%",x=-0.1,y=-0.3, color = "#ff0022", hjust = 1, size = 7) +
  annotate("text",label="17.2%",x=-0.1,y=0.3, color = "#011627", hjust = 1, size = 7) +
  annotate("text",label="28.3%",x=0.1,y=-0.3, color = "#011627", hjust = 0, size = 7) +
  ggrepel::geom_text_repel(label = gsea_aracne_isov2$common, color = "black", size = 6) +
  xlab("In vivo IL-iCAF") + ylab("In vitro IL-iCAF")
dev.off()

# TCGA paclitaxel/docetaxel analysis
ensg2eg <- as.list(org.Hs.egENSEMBL2EG) %>% lapply(function(x){x[1]})
eg2symbol <- as.list(org.Hs.egSYMBOL[mappedkeys(org.Hs.egSYMBOL)]) %>% lapply(function(x){x[1]})
sum_isoforms <- function(mat, gene_vec){
  dup_genes <- which(duplicated(gene_vec))
  finished <- c()
  for(i in dup_genes){
    goi <- gene_vec[i]
    if(!goi %in% finished){
      finished <- c(goi, finished)
      goi_id <- which(gene_vec==goi)
      mat[goi_id[1],] <- colSums(mat[goi_id,])
    }
  }
  mat <- mat[-dup_genes,]
  row.names(mat) <- gene_vec[-dup_genes]
  return(mat)
}

hnsc_tpm <- read.table("path to TPM data of TCGA-HNSC cohort",
                       header = T, sep = "\t")
hnsc_tpm <- hnsc_tpm[which(hnsc_tpm$Entrez_Gene_Id %in% names(eg2symbol)),] %>% .[-which(.$Hugo_Symbol %in% c(""," ")),]
tmp <- duplicated(hnsc_tpm$Hugo_Symbol)
row.names(hnsc_tpm) <- sapply(1:nrow(hnsc_tpm),function(x){if(tmp[x]){paste(hnsc_tpm$Hugo_Symbol[x],".1",sep="")}else{hnsc_tpm$Hugo_Symbol[x]}})
hnsc_tpm1 <- sum_isoforms(hnsc_tpm[,3:517],hnsc_tpm$Hugo_Symbol)
hist(log2(apply(hnsc_tpm1,1,median)+1), breaks=1000)
hnsc_tpm_filt <- hnsc_tpm1[which(apply(hnsc_tpm1,1,median)>2),]
hnsc_tpm_ssgsea <- cbind(data.frame(Genes=row.names(hnsc_tpm_filt)),hnsc_tpm_filt)
hnsc_tpm_ssgsea <- rbind(c(-666,rep("HNSC",ncol(hnsc_tpm_filt))),hnsc_tpm_ssgsea)
hnsc_tpm_ssgsea <- cbind(data.frame(id=c("Cancer_type",row.names(hnsc_tpm_filt))),hnsc_tpm_ssgsea)
write.table(hnsc_tpm_ssgsea, 
            file="desired filename",
            col.names = T, sep = "\t", quote=F, row.names=F)
hnsc_ssgsea <- read.table("ssGSEA output path", sep="\t", header=T, row.names=1, skip=2)
hnsc_ssgsea_nes <- hnsc_ssgsea[2:9,1034:1548]
for(i in 1:4){
  hnsc_ssgsea_nes[(i*2)-1,] <- as.numeric(hnsc_ssgsea_nes[(i*2)-1,])-as.numeric(hnsc_ssgsea_nes[(i*2),])
}
hnsc_ssgsea_nes <- hnsc_ssgsea_nes[c(7,3,5,1),]
for(i in 1:3){
  for(j in (i+1):4){
    hnsc_ssgsea_nes <- rbind(hnsc_ssgsea_nes,as.numeric(hnsc_ssgsea_nes[i,])-as.numeric(hnsc_ssgsea_nes[j,]))
  }
}
hnsc_ssgsea_nes <- rbind(hnsc_ssgsea_nes, as.numeric(hnsc_ssgsea_nes[2,])+as.numeric(hnsc_ssgsea_nes[3,]))
hnsc_ssgsea_nes <- rbind(hnsc_ssgsea_nes, as.numeric(hnsc_ssgsea_nes[1,])-as.numeric(hnsc_ssgsea_nes[11,]))
hnsc_ssgsea_nes <- rbind(hnsc_ssgsea_nes, as.numeric(hnsc_ssgsea_nes[11,])-as.numeric(hnsc_ssgsea_nes[4,]))
row.names(hnsc_ssgsea_nes) <- c("myCAF","ILiCAF","IFNiCAF","AdipoCAF","myCAFvILiCAF","myCAFvIFNiCAF","myCAFvAdipoCAF","ILiCAFvIFNiCAF","ILiCAFvAdipoCAF","IFNiCAFvAdipoCAF","iCAF","myCAFviCAF","iCAFvAdipoCAF")
hnsc_meta_all <- read.table("C:/Users/nicho/Desktop/NUS_PhD/230118 Microtubule TCGA/HNSC TCGA/hnsc_tcga_pan_can_atlas_2018_clinical_data_mrna.tsv",header=T,sep="\t")
hnsc_meta_all$Sample.ID <- sapply(hnsc_meta_all$Sample.ID,function(x){gsub("-",".",x)})
hnsc_meta_hpvneg <- subset(hnsc_meta_all, Subtype=="HNSC_HPV-")
hnsc_meta_hpvneg_clean <- hnsc_meta_hpvneg[-unique(c(which(is.na(hnsc_meta_hpvneg$Months.of.disease.specific.survival)),
                                                     which(is.na(hnsc_meta_hpvneg$Disease.specific.Survival.status)),
                                                     which(is.na(hnsc_meta_hpvneg$Overall.Survival..Months.)),
                                                     which(is.na(hnsc_meta_hpvneg$Overall.Survival.Status)),
                                                     which(is.na(hnsc_meta_hpvneg$Progress.Free.Survival..Months.)),
                                                     which(is.na(hnsc_meta_hpvneg$Progression.Free.Status)))),]
hnsc_meta_hpvneg_clean$Months.of.disease.specific.survival <- as.numeric(hnsc_meta_hpvneg_clean$Months.of.disease.specific.survival)
hnsc_meta_hpvneg_clean$Overall.Survival..Months. <- as.numeric(hnsc_meta_hpvneg_clean$Overall.Survival..Months.)
hnsc_meta_hpvneg_clean$Progress.Free.Survival..Months. <- as.numeric(hnsc_meta_hpvneg_clean$Progress.Free.Survival..Months.)
hnsc_meta_hpvneg_clean$Disease.specific.Survival.status <- sapply(hnsc_meta_hpvneg_clean$Disease.specific.Survival.status,function(x){if(grepl("0",x)){0}else{1}})
hnsc_meta_hpvneg_clean$Overall.Survival.Status <- sapply(hnsc_meta_hpvneg_clean$Overall.Survival.Status,function(x){if(grepl("0",x)){0}else{1}})
hnsc_meta_hpvneg_clean$Progression.Free.Status <- sapply(hnsc_meta_hpvneg_clean$Progression.Free.Status,function(x){if(grepl("0",x)){0}else{1}})
hnsc_ssgsea_nes_hpvneg <- hnsc_ssgsea_nes[,which(colnames(hnsc_ssgsea_nes) %in% hnsc_meta_hpvneg_clean$Sample.ID)]
for(i in 1:nrow(hnsc_ssgsea_nes_hpvneg)){
  tmp <- quantile(as.numeric(hnsc_ssgsea_nes_hpvneg[i,]),probs = c(1/3,2/3))
  hnsc_meta_hpvneg_clean[[row.names(hnsc_ssgsea_nes_hpvneg)[i]]] <- sapply(as.numeric(hnsc_ssgsea_nes_hpvneg[i,]),function(x){if(x<=tmp[1]){0}else if(x<=tmp[2]){1}else{2}})
  rm(tmp)
}
hnsc_pd_meta <- read.table("path to file containing metadata of patients treated with paclitaxel/docetaxel",
                           header = T, sep = "\t")

#DSS
hnsc_hpvneg_pd_dss_km <- Surv(time = hnsc_meta_hpvneg_clean$Months.of.disease.specific.survival[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID)],
                              event = hnsc_meta_hpvneg_clean$Disease.specific.Survival.status[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID)])
hnsc_hpvneg_pd_dss_km_treatment <- list(myCAF=survfit(hnsc_hpvneg_pd_dss_km~myCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        ILiCAF=survfit(hnsc_hpvneg_pd_dss_km~ILiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        IFNiCAF=survfit(hnsc_hpvneg_pd_dss_km~IFNiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        AdipoCAF=survfit(hnsc_hpvneg_pd_dss_km~AdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        myCAFvILiCAF=survfit(hnsc_hpvneg_pd_dss_km~myCAFvILiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        myCAFvIFNiCAF=survfit(hnsc_hpvneg_pd_dss_km~myCAFvIFNiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        myCAFvAdipoCAF=survfit(hnsc_hpvneg_pd_dss_km~myCAFvAdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        ILiCAFvIFNiCAF=survfit(hnsc_hpvneg_pd_dss_km~ILiCAFvIFNiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        ILiCAFvAdipoCAF=survfit(hnsc_hpvneg_pd_dss_km~ILiCAFvAdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        IFNiCAFvAdipoCAF=survfit(hnsc_hpvneg_pd_dss_km~IFNiCAFvAdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        iCAF=survfit(hnsc_hpvneg_pd_dss_km~iCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        myCAFviCAF=survfit(hnsc_hpvneg_pd_dss_km~myCAFviCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        iCAFvAdipoCAF=survfit(hnsc_hpvneg_pd_dss_km~iCAFvAdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'))

#OS
hnsc_hpvneg_pd_os_km <- Surv(time = hnsc_meta_hpvneg_clean$Overall.Survival..Months.[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID)],
                             event = hnsc_meta_hpvneg_clean$Overall.Survival.Status[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID)])
hnsc_hpvneg_pd_os_km_treatment <- list(myCAF=survfit(hnsc_hpvneg_pd_os_km~myCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                       ILiCAF=survfit(hnsc_hpvneg_pd_os_km~ILiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                       IFNiCAF=survfit(hnsc_hpvneg_pd_os_km~IFNiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                       AdipoCAF=survfit(hnsc_hpvneg_pd_os_km~AdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                       myCAFvILiCAF=survfit(hnsc_hpvneg_pd_os_km~myCAFvILiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                       myCAFvIFNiCAF=survfit(hnsc_hpvneg_pd_os_km~myCAFvIFNiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                       myCAFvAdipoCAF=survfit(hnsc_hpvneg_pd_os_km~myCAFvAdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                       ILiCAFvIFNiCAF=survfit(hnsc_hpvneg_pd_os_km~ILiCAFvIFNiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                       ILiCAFvAdipoCAF=survfit(hnsc_hpvneg_pd_os_km~ILiCAFvAdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                       IFNiCAFvAdipoCAF=survfit(hnsc_hpvneg_pd_os_km~IFNiCAFvAdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                       iCAF=survfit(hnsc_hpvneg_pd_os_km~iCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                       myCAFviCAF=survfit(hnsc_hpvneg_pd_os_km~myCAFviCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                       iCAFvAdipoCAF=survfit(hnsc_hpvneg_pd_os_km~iCAFvAdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'))

#PFS
hnsc_hpvneg_pd_pfs_km <- Surv(time = hnsc_meta_hpvneg_clean$Progress.Free.Survival..Months.[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID)],
                              event = hnsc_meta_hpvneg_clean$Progression.Free.Status[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID)])
hnsc_hpvneg_pd_pfs_km_treatment <- list(myCAF=survfit(hnsc_hpvneg_pd_pfs_km~myCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        ILiCAF=survfit(hnsc_hpvneg_pd_pfs_km~ILiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        IFNiCAF=survfit(hnsc_hpvneg_pd_pfs_km~IFNiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        AdipoCAF=survfit(hnsc_hpvneg_pd_pfs_km~AdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        myCAFvILiCAF=survfit(hnsc_hpvneg_pd_pfs_km~myCAFvILiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        myCAFvIFNiCAF=survfit(hnsc_hpvneg_pd_pfs_km~myCAFvIFNiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        myCAFvAdipoCAF=survfit(hnsc_hpvneg_pd_pfs_km~myCAFvAdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        ILiCAFvIFNiCAF=survfit(hnsc_hpvneg_pd_pfs_km~ILiCAFvIFNiCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        ILiCAFvAdipoCAF=survfit(hnsc_hpvneg_pd_pfs_km~ILiCAFvAdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        IFNiCAFvAdipoCAF=survfit(hnsc_hpvneg_pd_pfs_km~IFNiCAFvAdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        iCAF=survfit(hnsc_hpvneg_pd_pfs_km~iCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        myCAFviCAF=survfit(hnsc_hpvneg_pd_pfs_km~myCAFviCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'),
                                        iCAFvAdipoCAF=survfit(hnsc_hpvneg_pd_pfs_km~iCAFvAdipoCAF,data=hnsc_meta_hpvneg_clean[which(hnsc_meta_hpvneg_clean$Patient.ID %in% hnsc_pd_meta$Patient.ID),],type='kaplan-meier',conf.type='log'))
