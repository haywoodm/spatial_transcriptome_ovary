## Load Libraries
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(stringi)
library(RColorBrewer)
library(pheatmap)

## Read in file from 3_new_cluster_ids.R script
slide.integrated <- readRDS(file = "slide.integrated.newIDs.rds")

## subset data into subcluster of interest (CL-R, Stroma, Follicle, or CL-P)
d <- subset(slide.integrated, idents = c("Follicle"))

## Proceed with downstream analysis on the integrated dataset
d <- RunPCA(object = d, verbose = FALSE, npcs = 30, assay = "integrated")
d <- FindNeighbors(d, reduction = "pca", dims = 1:30, k.param = 20)
d <- FindClusters(d, verbose = FALSE, resolution = 0.7)
d <- RunUMAP(d, reduction = "pca", dims = 1:30)

#####################################################
## Figure 7                                        ##
#####################################################
## Plot UMAPs
plot1 <- DimPlot(d, reduction = "umap", label = TRUE, label.size = 6) +
  plot_annotation(title = '') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
plot2 <- DimPlot(d, reduction = "umap", group.by = "group") +
  plot_annotation(title = '') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
plotUMAP <- plot1 + plot2 &
  labs(x = "UMAP1", y = "UMAP2", title = "") &
  theme(text = element_text(size=20), axis.text = element_text(size = 20))

## group proportions per cluster
clusters <- as.data.frame((d@meta.data) %>%
                            group_by(group, seurat_clusters, .drop = FALSE) %>%
                            summarize(count = n()) %>%
                            group_by(group) %>%
                            mutate(Proportion = count/sum(count)))

aclusters <- filter(clusters, group == "aged")
yclusters <- filter(clusters, group == "young")
aclusters$remaining <- sum(aclusters$count)-aclusters$count
yclusters$remaining <- sum(yclusters$count)-yclusters$count
clusters <- left_join(yclusters, aclusters, by = "seurat_clusters")

## fisher.exact on all rows
fisher <- apply(clusters, 1, 
                function(x) {
                  tbl <- matrix(as.numeric(x[c(3,5,7,9)]), ncol=2, byrow=T)
                  fisher.test(tbl, alternative="two.sided")
                })
clusters$p <- unlist(lapply(fisher, function(i) i$p.value))
clusters$OR <- unlist(lapply(fisher, function(i) i$estimate))
clusters$FDR <- p.adjust(clusters$p, method = "BH")
clusters$logOR <- log10(clusters$OR) 
clusters$logFDR <- -log10(clusters$FDR)

## replace infinite ORs with +-2 for plotting purposes
clusters$logOR <- ifelse(clusters$logOR == "Inf", 2,
                         ifelse(clusters$logOR == "-Inf", -2, clusters$logOR))

clusters$total <- clusters$count.x + clusters$count.y
clusters$total_prop <- clusters$total/(sum(clusters$count.x + clusters$count.y))

## plot FDR/OR
plot3 <- ggplot(clusters, aes(x = logOR, y = logFDR, color = seurat_clusters, label = seurat_clusters)) +
  geom_point(size = 5) +
  geom_text_repel(min.segment.length = 0.5,
                  box.padding = 0.5, point.padding = 0.5, max.overlaps = 20,
                  size = 4, fontface = "bold", segment.size = 1, nudge_y = .4) +
  geom_hline(yintercept = 1.3) +
  geom_vline(xintercept = c((-log10(0.67)), (-log10(1.5))), linetype = 2) +
  labs(y = "-log10 FDR", x = "log10 Odds Ratio") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20))

## Plot Figure 7a/b/c/d
plotUMAP + plot3

#####################################################
## Table 3 information                             ##
#####################################################
names <- ifelse(clusters$count.x != 0 & clusters$count.y != 0, levels(d@meta.data$seurat_clusters), "")
names <- stri_remove_empty(names)
gene_list <- list()

for (i in 1:length(names)) {
  gene_list[[i]] <- FindMarkers(d, assay = "Spatial", ident.1 = "aged", group.by = "group", subset.ident = names[i],
                                logfc.threshold = 0, latent.vars = "slide", test.use = "LR")
}
names(gene_list) <- names

#####################################################
## Supplementary Figure 3 (Follicles only)         ##
#####################################################
counts <- GetAssayData(d, assay = "Spatial", slot = "data")

## Genes of interest
genes <- c("Amh", "Fshr", "Cyp11a1", "Cyp19a1", "Hsd17b1", "Esr1", "Esr2", "Ar", "Pgr",
           "Sohlh1", "Nobox", "Kitlg", "Bmp15", "Gdf9", "Tbp2", "Zar1", "Gja1", "Gja4",
           "Lhcgr", "Foxo3", "Inha", "Inhbb")

counts <- as.data.frame((counts[rownames(counts) %in% genes, ]))

cluster.averages <- as.data.frame(AverageExpression(d, assay = "Spatial", slot = "counts"))
cluster.averages <- as.data.frame((cluster.averages[rownames(cluster.averages) %in% genes, ]))

names(cluster.averages) <- names

## heatmap
hr <- hclust(as.dist(1-cor(t(cluster.averages), method="spearman")), method="average")
hc <- hclust(as.dist(1-cor(cluster.averages, method="spearman")), method="average")
breaksList = seq(-1, 1, by = 0.1)

## Plot Supplementary Figure 3
pheatmap(cluster.averages, scale = "row",
         cluster_cols = hc, cluster_rows = hr, show_colnames = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         annotation_names_col = FALSE, show_rownames = TRUE, breaks = breaksList,
         legend_breaks = c(-0.9,0,0.9,1), legend_labels = c("-1", "0", "1","Z-score\n\n"))

##############################################################
## Supplementary Figures 4 and 5 (Follicle and Oocyte only) ##
##############################################################

images <- c("D655", "D656", "D657", "D658")

## Plot Figures 4 and 5
for (i in 1:length(images)) {
  plot <- SpatialDimPlot(d, images = images[i], crop = FALSE, pt.size.factor = 1.5) +
    labs(fill = "")
  print(plot)
}