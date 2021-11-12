## Load Libraries
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(ggrepel)

## Read in file from 2_global_analysis.R script
slide.integrated <- readRDS(file = "slide.integrated.global.rds")

## Assign new cluster IDs
new.cluster.ids <- c("CL-R", "Stroma", "Stroma", "CL-R", "Stroma", "Stroma",
                     "CL-R", "Follicle", "CL-R", "Follicle", "CL-P",
                     "Epithelium 1", "CL-P", "Stroma", "Stroma", "Follicle",
                     "Follicle", "CL-P", "CL-R", "Follicle", "Follicle",
                     "CL-P", "Epithelium 2", "Epithelium 3")
names(new.cluster.ids) <- levels(slide.integrated)
slide.integrated <- RenameIdents(slide.integrated, new.cluster.ids)
slide.integrated[["new.cluster.ids"]] <- Idents(object = slide.integrated)

## Save new ID analysis
saveRDS(slide.integrated, file = "slide.integrated.newIDs.rds")

#####################################################
## Figure 3a                                       ##
#####################################################
## Plot Figure 3a
DimPlot(slide.integrated, reduction = "umap", label = TRUE, label.size = 4.75) +
  plot_annotation(title = '') +
  labs(x = "UMAP1", y = "UMAP2", title = "") +
  NoLegend() +
  expand_limits(x = 11) +
  theme(text = element_text(size=20),
        axis.text = element_text(size = 20))

#####################################################
## Figure 3b  and 3d                               ##
#####################################################
## group proportions per cluster
clusters <- as.data.frame((slide.integrated@meta.data) %>%
                            group_by(group, new.cluster.ids, .drop = FALSE) %>%
                            summarize(count = n()) %>%
                            group_by(group) %>%
                            mutate(Proportion = count/sum(count)))

aclusters <- filter(clusters, group == "aged")
yclusters <- filter(clusters, group == "young")
aclusters$remaining <- sum(aclusters$count)-aclusters$count
yclusters$remaining <- sum(yclusters$count)-yclusters$count
clusters <- left_join(yclusters, aclusters, by = "new.cluster.ids")

clusters$total <- clusters$count.x + clusters$count.y
clusters$total_prop <- clusters$total/(sum(clusters$count.x + clusters$count.y))

## Plot Figure 3b
ggplot(clusters, aes(x = new.cluster.ids, y = total_prop, fill = new.cluster.ids)) +
  geom_col() +
  labs(x = "", y = "Proportion of Total Spots") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=20)) +
  theme(plot.margin = margin(10, 10, 10, 50))

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

## Plot Figure 3d
ggplot(clusters, aes(x = logOR, y = logFDR, color = new.cluster.ids, label = new.cluster.ids)) +
  geom_point(size = 5) +
  geom_text_repel(min.segment.length = 0.5,
                  box.padding = 0.5, point.padding = 0.5, max.overlaps = 20,
                  size = 4, fontface = "bold", segment.size = 1, nudge_y = .4) +
  geom_hline(yintercept = 1.3) +
  geom_vline(xintercept = c((-log10(0.67)), (-log10(1.5))), linetype = 2) +
  labs(y = "-log10 FDR", x = "log10 Odds Ratio") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size=20))

#####################################################
## Figure 3c                                       ##
#####################################################
## summarizing mouse proportions per cluster
clusters <- as.data.frame(slide.integrated@meta.data) %>%
  group_by(mouse, new.cluster.ids, .drop = FALSE) %>%
  summarize(count = n()) %>%
  group_by(mouse) %>%
  mutate(Proportion = count/sum(count))

## Plot Figure 3c
ggplot(clusters, aes(x = mouse, y = Proportion)) +
  geom_col(aes(fill = new.cluster.ids)) +
  labs(fill = "Seurat Clusters", x = "Sample") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(text = element_text(size=20))

#####################################################
## Figure 2b                                       ##
#####################################################
## find cluster markers
allMarkersSpatial <- FindAllMarkers(slide.integrated, assay = "Spatial", only.pos = TRUE, latent.vars = "slide", test.use = "LR")
top10 <- allMarkersSpatial %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

## Plot Figure 2b
DoHeatmap(slide.integrated, features = top10$gene, size = 4.5) +
  theme(plot.margin = margin(50, 10, 10, 10))

#####################################################
## Figure 4                                        ##
#####################################################
images <- c("D655", "D656", "D657", "D658")

## Plot Figure 4
for (i in 1:length(images)) {
  plot <- SpatialDimPlot(slide.integrated, images = images[i], crop = FALSE, pt.size.factor = 1.5) +
    labs(fill = "")
  print(plot)
}

#####################################################
## Table 1 information                             ##
#####################################################
names <- unique(new.cluster.ids)
gene_list <- list()

for (i in 1:length(names)) {
  gene_list[[i]] <- FindMarkers(slide.integrated, assay = "Spatial", ident.1 = "aged", group.by = "group", subset.ident = names[i],
                           logfc.threshold = 0, latent.vars = "slide", test.use = "LR")
}
names(gene_list) <- names

#####################################################
## Supplementary Table 1 information               ##
#####################################################
names <- unique(new.cluster.ids)
gene_list <- list()

for (i in 1:length(names)) {
  gene_list[[i]] <- FindMarkers(slide.integrated, assay = "Spatial", ident.1 = names[i], ident.2 = NULL,
                                logfc.threshold = 0, latent.vars = "slide", test.use = "LR")
}
names(gene_list) <- names