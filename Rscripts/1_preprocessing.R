## Load Libraries
library(dplyr)
library(Seurat)

## Load the data as outlined in the 10x Genomics pipeline
## https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/rkit
bcs_merge <- readRDS(file = "bcs_merge.rds")

## divide each image by young and aged sections
## for slide D21-655: young is on the top and aged is on the bottom
## for slides D21-656, 657, 658: young is on the left and aged is on the right
bcs_merge$group <- as.factor(
  ifelse((bcs_merge$slide == "D21-655" & bcs_merge$row < 37), "aged",
         ifelse((bcs_merge$slide == "D21-656" & bcs_merge$col > 64), "aged",
                ifelse((bcs_merge$slide == "D21-657" & bcs_merge$col > 74), "aged",
                       ## even row number 0-30, column > 56 = aged
                       ifelse((bcs_merge$slide == "D21-658" &
                                 (bcs_merge$row %in% c(seq(0,30, by = 2))) & bcs_merge$col > 56), "aged",
                              ## odd row number 0-30, column > 55 = aged
                              ifelse((bcs_merge$slide == "D21-658" &
                                        (bcs_merge$row %in% c(seq(1,29, by = 2))) &
                                        bcs_merge$col > 55), "aged",
                                     ## even row number 31-77, column > 54 = aged
                                     ifelse((bcs_merge$slide == "D21-658" &
                                               (bcs_merge$row %in% c(seq(32,76, by = 2))) & bcs_merge$col > 54), "aged",
                                            ## odd row number 31-77, column > 55 = aged
                                            ifelse((bcs_merge$slide == "D21-658" &
                                                      (bcs_merge$row %in% c(seq(31,77, by = 2))) &
                                                      bcs_merge$col > 55), "aged", "young"))))))))

## Add in metadata information about group, slide, and mouse
bcs_merge$mouse <- ifelse(bcs_merge$group == "aged" & bcs_merge$slide == "D21-655", "A1",
                          ifelse(bcs_merge$group == "aged" & bcs_merge$slide == "D21-656", "A2",
                                 ifelse(bcs_merge$group == "aged" & bcs_merge$slide == "D21-657", "A3",
                                        ifelse(bcs_merge$group == "aged" & bcs_merge$slide == "D21-658", "A4",
                                               ifelse(bcs_merge$group == "young" & bcs_merge$slide == "D21-655", "Y1",
                                                      ifelse(bcs_merge$group == "young" & bcs_merge$slide == "D21-656", "Y2",
                                                             ifelse(bcs_merge$group == "young" & bcs_merge$slide == "D21-657", "Y3", "Y4")))))))

## delineate each metadata set for later use
D21_655_meta <- filter(bcs_merge, slide == "D21-655")
D21_656_meta <- filter(bcs_merge, slide == "D21-656")
D21_657_meta <- filter(bcs_merge, slide == "D21-657")
D21_658_meta <- filter(bcs_merge, slide == "D21-658")
meta <- list(D21_655_meta, D21_656_meta, D21_657_meta, D21_658_meta)

## Load each Visium slide
slides <- c("D21-655", "D21-656", "D21-657", "D21-658")
slices <- c("D655", "D656", "D657", "D658")
slide.list <- list()

for (i in 1:length(slides)) {
  slide.list[[i]] <- Load10X_Spatial(paste("C:/Path/to/slides/",slides[i],"/outs/", sep =""),
                                     filename = "filtered_feature_bc_matrix.h5",
                                     slice = paste(slices[i]))
  
  ## Calculate percent mitochondria content of each spot for QC
  slide.list[[i]][["percent.mt"]] <- PercentageFeatureSet(slide.list[[i]], pattern = "^mt-")
  
}

names(slide.list) <- slices
rm(slices)

## Add in metadata to Seurat object
for (i in 1:length(slide.list)){
  seurat <- slide.list[[i]]@meta.data
  seurat <- tibble::rownames_to_column(seurat, var = "barcode")
  metadata <- left_join(meta[[i]], seurat)
  metadata <- tibble::column_to_rownames(metadata, var = "barcode")
  slide.list[[i]] <- AddMetaData(slide.list[[i]], metadata)
  rm(seurat)
  rm(metadata)
}

## Subset out specific non-tissue spots (< 25% tissue coverage)
slide.list$D655 <- subset(
  x = slide.list$D655,
  subset = ((row != 46 | col != 52) &
              (row != 43 | col != 15) &
              (row != 47 | col != 53) &
              (row != 49 | col != 55) &
              (row != 48 | col != 58) &
              (row != 48 | col != 60) &
              (row != 48 | col != 54) &
              (row != 47 | col != 59) &
              (row != 46 | col != 54) &
              (row != 48 | col != 52)))

slide.list$D656 <- subset(
  x = slide.list$D656,
  subset = ((row != 22 | col != 22) &
              (row != 60 | col != 18) &
              (row != 60 | col != 20) &
              (row != 72 | col != 20) &
              (row != 64 | col != 30) &
              (row != 40 | col != 30) &
              (row != 41 | col != 31) &
              (row != 49 | col != 41) &
              (row != 52 | col != 42) &
              (row != 52 | col != 16) &
              (row != 52 | col != 18) &
              (row != 53 | col != 23) &
              (row != 53 | col != 25) &
              (row != 57 | col != 35) &
              (row != 40 | col != 96) &
              (row != 69 | col != 99) &
              (row != 44 | col != 102) &
              (row != 40 | col != 94) &
              (row != 38 | col != 98) &
              (row != 39 | col != 99) &
              (row != 70 | col != 100) &
              (row != 39 | col != 93) &
              (row != 9 | col != 101) &
              (row != 39 | col != 95) &
              (row != 38 | col != 92) &
              (row != 38 | col != 96) &
              (row != 40 | col != 98) &
              (row != 40 | col != 108) &
              (row != 38 | col != 94) &
              (row != 39 | col != 97) &
              (row != 69 | col != 101) &
              (row != 0 | col != 4)))

slide.list$D657 <- subset(
  x = slide.list$D657,
  subset = ((row != 41 | col != 107) &
              (row != 10 | col != 126) &
              (row != 9 | col != 125) &
              (row != 7 | col != 127) &
              (row != 9 | col != 127) &
              (row != 8 | col != 126) &
              (row != 42 | col != 106) &
              (row != 41 | col != 105) &
              (row != 52 | col != 48) &
              (row != 47 | col != 63) &
              (row != 77 | col != 63) &
              (row != 51 | col != 47) &
              (row != 32 | col != 40) &
              (row != 77 | col != 67) &
              (row != 9 | col != 21)))

slide.list$D658 <- subset(
  x = slide.list$D658,
  subset = ((row != 62 | col != 98) &
              (row != 5 | col != 103) &
              (row != 45 | col != 73) &
              (row != 63 | col != 99) &
              (row != 42 | col != 68) &
              (row != 44 | col != 70) &
              (row != 16 | col != 28) &
              (row != 16 | col != 32) &
              (row != 51 | col != 47) &
              (row != 16 | col != 30) &
              (row != 65 | col != 13) &
              (row != 29 | col != 23) &
              (row != 15 | col != 31) &
              (row != 10 | col != 14) &
              (row != 20 | col != 0)))


## Setup the Seurat object list, and run SCTransform on each object individually
options(future.globals.maxSize = 4000 * 1024^2)
for (i in 1:length(slide.list)) {
  slide.list[[i]] <- SCTransform(slide.list[[i]], assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
}

## Select features for downstream integration, and prepare to integrate datasets
slide.features <- SelectIntegrationFeatures(object.list = slide.list, nfeatures = 3000)
slide.list <- PrepSCTIntegration(object.list = slide.list, anchor.features = slide.features, 
                                 verbose = FALSE)

## Identify anchors and integrate the datasets
slide.anchors <- FindIntegrationAnchors(object.list = slide.list, normalization.method = "SCT", 
                                        anchor.features = slide.features, verbose = FALSE)
slide.integrated <- IntegrateData(anchorset = slide.anchors, normalization.method = "SCT", 
                                  verbose = FALSE)

## Save integrated object
saveRDS(slide.integrated, file = "slide.integrated.preprocessed.rds")