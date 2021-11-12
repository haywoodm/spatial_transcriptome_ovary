## Load Libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)

## Read in file from 1_preprocessing.R script
slide.integrated <- readRDS(file = "slide.integrated.preprocessed.rds")

## Proceed with downstream analysis on the integrated dataset
slide.integrated <- RunPCA(object = slide.integrated, verbose = FALSE, npcs = 50)
slide.integrated <- FindNeighbors(slide.integrated, reduction = "pca", dims = 1:50, k.param = 10)
slide.integrated <- FindClusters(slide.integrated, verbose = FALSE, resolution = 0.9)
slide.integrated <- RunUMAP(slide.integrated, reduction = "pca", dims = 1:50)

## Save integrated analysis
saveRDS(slide.integrated, file = "slide.integrated.global.rds")

#####################################################
## Figure 1b and 1c                                ##
#####################################################
## Plot Figure 1b and 1c
plot1 <- DimPlot(slide.integrated, reduction = "umap", label = TRUE, label.size = 6) +
  plot_annotation(title = '')
plot2 <- DimPlot(slide.integrated, reduction = "umap", group.by = "group", label.size = 6) +
  plot_annotation(title = '')
plot1 + plot2 &
  labs(x = "UMAP1", y = "UMAP2", title = "") &
  theme(text = element_text(size=20), axis.text = element_text(size = 20))

#####################################################
## Figure 1d                                       ##
#####################################################
## summarizing mouse proportions per cluster
clusters <- as.data.frame(slide.integrated@meta.data) %>%
  group_by(mouse, seurat_clusters, .drop = FALSE) %>%
  summarize(count = n()) %>%
  group_by(mouse) %>%
  mutate(Proportion = count/sum(count))

## Plot Figure 1d
ggplot(clusters, aes(x = mouse, y = Proportion)) +
  geom_col(aes(fill = seurat_clusters)) +
  labs(fill = "Seurat Clusters", x = "Sample") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(text = element_text(size=20))

#####################################################
## Figure 1e                                       ##
#####################################################
## group proportions per cluster
clusters <- as.data.frame((slide.integrated@meta.data) %>%
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

## Plot Figure 1e
ggplot(clusters, aes(x = logOR, y = logFDR, color = seurat_clusters, label = seurat_clusters)) +
  geom_point(size = 5) +
  geom_text_repel(min.segment.length = 0.5,
                  box.padding = 0.4, point.padding = 0.5, max.overlaps = 20,
                  size = 6, fontface = "bold", segment.size = 1, nudge_y = 2) +
  geom_hline(yintercept = 1.3) +
  xlim(-2, 1.4) +
  geom_vline(xintercept = c((-log10(0.5)), (-log10(2))), linetype = 2) +
  labs(y = "-log10 FDR", x = "log10 Odds Ratio") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size=20))

#####################################################
## Figure 2a                                       ##
#####################################################

## find cluster markers
allMarkersSpatial <- FindAllMarkers(slide.integrated, assay = "Spatial", only.pos = TRUE, latent.vars = "slide", test.use = "LR")
top10 <- allMarkersSpatial %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

## Plot Figure 2a
DoHeatmap(slide.integrated, features = top10$gene, size = 4.5) +
  theme(axis.text.y = element_blank())