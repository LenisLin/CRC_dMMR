library(Seurat)
library(Matrix)

dataPath <- "/mnt/raid5/ProjectData/HorizontalProject/zsz/Data/"

sce_pre_tumor <- readRDS(paste0(dataPath, "zsz_updated_sce.rds"))
sce_post_tumor <- readRDS(paste0(dataPath, "zsz_Postreat-Tumor-NONPCR(includeZZM16).rds"))
sce_dyh_normal <- readRDS(paste0(dataPath, "DYH_N.rds"))
sce_zzm_normal <- readRDS(paste0(dataPath, "ZZM_N.rds"))

## view and process each element

### sce_zzm_normal
if (1) {
  rowSums(sce_zzm_normal) ## sce_zzm_normal is a matrix with counts
  sce_zzm_normal[1:5, 1:5]
}

### sce_dyh_normal
if (1) {
  head(sce_dyh_normal@meta.data)
  sce_dyh_normal_meta <- sce_dyh_normal@meta.data

  layers <- names(sce_dyh_normal@assays$RNA@layers)
  count_matrices <- lapply(layers, function(layer) {
    GetAssayData(sce_dyh_normal, slot = layer)
  })

  ## filter features and combine
  features <- rownames(count_matrices[[1]])
  for (i in 2:length(count_matrices)) {
    features <- intersect(features, rownames(count_matrices[[i]]))
  }

  subset_matrices <- lapply(count_matrices, function(mat) {
    mat[features, , drop = FALSE]
  })

  sce_dyh_normal_exp <- do.call(cbind, subset_matrices)
  dim(sce_dyh_normal_exp)

  rm(subset_matrices, count_matrices, layers)
}

### sce_post_tumor
head(sce_post_tumor@meta.data)
sce_post_tumor_meta <- sce_post_tumor@meta.data
sce_post_tumor_exp <- sce_post_tumor@assays$RNA@counts

### sce_pre_tumor
tail(sce_pre_tumor@meta.data)
sce_pre_tumor_meta <- sce_pre_tumor@meta.data
sce_pre_tumor_exp <- sce_pre_tumor@assays$RNA@counts

### Combine all features
if (1) {
rownames <- c(colnames(sce_zzm_normal),rownames(sce_dyh_normal_meta),rownames(sce_post_tumor_meta),rownames(sce_pre_tumor_meta))

tissues <- c(rep("Normal",length(c(colnames(sce_zzm_normal),rownames(sce_dyh_normal_meta)))),rep("Tumor",length(c(rownames(sce_post_tumor_meta),rownames(sce_pre_tumor_meta)))))
timepoints <- c(rep(NA,length(tissues[tissues=="Normal"])),rep("Post",length(rownames(sce_post_tumor_meta))),rep("Pre",length(rownames(sce_pre_tumor_meta))))

duplicateCell <- rownames[duplicated(rownames)]
removeidx <- rownames %in% duplicateCell

all_meta <- as.data.frame(matrix(data = NA,nrow = length(rownames),ncol = 3))
colnames(all_meta) <- c("cellid","tissue","timepoint")
all_meta$cellid <- rownames
all_meta$tissue <- tissues
all_meta$timepoint <- timepoints

all_meta$orig.ident <- sapply(all_meta$cellid,function(x){
  return(unname(strsplit(x,"_")[[1]][1]))
})

all_meta$group <- sce_pre_tumor_meta[match(all_meta$orig.ident,sce_pre_tumor_meta$orig.ident),"group"]
head(all_meta)
dim(all_meta)

all_features <- intersect(rownames(sce_zzm_normal),rownames(sce_dyh_normal_exp))
all_features <- intersect(all_features,rownames(sce_post_tumor_exp))
all_features <- intersect(all_features,rownames(sce_pre_tumor_exp))
all_exp <- cbind(
  sce_zzm_normal[match(all_features,rownames(sce_zzm_normal)),],
  sce_dyh_normal_exp[match(all_features,rownames(sce_dyh_normal_exp)),],
  sce_post_tumor_exp[match(all_features,rownames(sce_post_tumor_exp)),],
  sce_pre_tumor_exp[match(all_features,rownames(sce_pre_tumor_exp)),]
)

dim(all_exp)

identical(colnames(all_exp),all_meta$cellid)

all_exp <- all_exp[,!removeidx]
all_meta <- all_meta[!removeidx,]

rownames(all_meta) <- all_meta$cellid
all_meta <- all_meta[,-1]

}

seu <- CreateSeuratObject(counts = all_exp,meta.data = all_meta)

# Save
saveRDS(seu,"/mnt/raid5/ProjectData/HorizontalProject/zsz/Data/all_sce.rds")

# Save the row and column names to separate files
savePath <- "/home/lenislin/Experiment/projects/HorizontalProject/zsz_CRC/Data/Step1_preProcessing/"
writeLines(rownames(all_exp), paste0(savePath,"exp_rownames.tx"))
writeLines(colnames(all_exp), paste0(savePath,"exp_colnames.txt"))

writeMM(all_exp, paste0(savePath,"exp.mtx"))
write.csv(all_meta, paste0(savePath,"meta.csv"), row.names = T)
