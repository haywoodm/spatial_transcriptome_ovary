## Load Libraries
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(stringi)

## Read in file from 3_new_cluster_ids.R script
slide.integrated <- readRDS(file = "slide.integrated.global.rds")

## Subset out spots expressing Gdf9 and Zp3
DefaultAssay(slide.integrated) <- "Spatial"
oocyte <- subset(x = slide.integrated, subset = (Gdf9 > 0 & Zp3 > 0))

## Proceed with downstream analysis on the integrated dataset
DefaultAssay(oocyte) <- "integrated"
oocyte <- RunPCA(object = oocyte, verbose = FALSE, npcs = 10, assay = "integrated")
oocyte <- FindNeighbors(oocyte, reduction = "pca", dims = 1:10, k.param = 5)
oocyte <- FindClusters(oocyte, verbose = FALSE, resolution = 0.3)
oocyte <- RunUMAP(oocyte, reduction = "pca", dims = 1:10)

#####################################################
## Figure 8a and 8b                                ##
#####################################################
## Plot Figure 8a and 8b
plot1 <- DimPlot(oocyte, reduction = "umap", label = TRUE, label.size = 6) +
  plot_annotation(title = '') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
plot2 <- DimPlot(oocyte, reduction = "umap", group.by = "group") +
  plot_annotation(title = '') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
plot1 + plot2 &
  labs(x = "UMAP1", y = "UMAP2", title = "") &
  theme(text = element_text(size=20), axis.text = element_text(size = 20))

#####################################################
## Figure 8c                                       ##
#####################################################
## find cluster markers
allMarkersSpatial <- FindAllMarkers(oocyte, assay = "Spatial", only.pos = TRUE, latent.vars = "slide", test.use = "LR")
top10 <- allMarkersSpatial %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

## Plot Figure 8c
DoHeatmap(oocyte, features = top10$gene, size = 4.5)

###################################################
## Figure 8d                                     ##
###################################################
## Assign new cluster IDs
new.cluster.ids <- c("oocyte", "oocyte", "oocyte", "oocyte", "not oocyte", "oocyte", "not oocyte")
names(new.cluster.ids) <- levels(oocyte)
oocyte <- RenameIdents(oocyte, new.cluster.ids)
oocyte[["new.cluster.ids"]] <- Idents(object = oocyte)

## Plot Figure 8d
plot1 <- DimPlot(oocyte, reduction = "umap", label = TRUE, label.size = 6) +
  plot_annotation(title = 'Oocyte') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))

plot1 + labs(x = "UMAP1", y = "UMAP2", title = "") +  theme(text = element_text(size=20), axis.text = element_text(size = 20))

###################################################
## Figure 8e and 8f                              ##
###################################################
## Subset oocyte clusters and proceed with downstream analysis
d <- subset(oocyte, idents = c("oocyte"))
d <- RunPCA(object = d, verbose = FALSE, npcs = 10, assay = "integrated")
d <- FindNeighbors(d, reduction = "pca", dims = 1:10, k.param = 10)
d <- FindClusters(d, verbose = FALSE, resolution = 0.3)
d <- RunUMAP(d, reduction = "pca", dims = 1:10)

## Plot Figure 8e and 8f
plot1 <- DimPlot(d, reduction = "umap", label = TRUE, label.size = 6) +
  plot_annotation(title = 'Oocyte') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
plot2 <- DimPlot(d, reduction = "umap", group.by = "group") +
  plot_annotation(title = '') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))

plot1 + labs(x = "UMAP1", y = "UMAP2", title = "") +  theme(text = element_text(size=20), axis.text = element_text(size = 20))
plot2 + labs(x = "UMAP1", y = "UMAP2", title = "") +  theme(text = element_text(size=20), axis.text = element_text(size = 20))

###################################################
## Figure 8g                                     ##
###################################################
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

## Plot Figure 8g
ggplot(clusters, aes(x = logOR, y = logFDR, color = seurat_clusters, label = seurat_clusters)) +
  geom_point(size = 5) +
  geom_text_repel(min.segment.length = 0.5,
                  box.padding = 0.5, point.padding = 0.5, max.overlaps = 20,
                  size = 4, fontface = "bold", segment.size = 1) +
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