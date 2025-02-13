# from https://github.com/gbouland/binary-differential-analysis

library(magrittr)
counts <- rio::import("data/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz")

samplesheet <- rio::import("data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz")

rownames(counts) <- counts$Index
counts$Index <- NULL

samplesheet <- samplesheet[,c("Index","Cell_type","Class","Patient")]
colnames(samplesheet) <- c("IDs","cellype","status","patient")
selection <- table(samplesheet$cellype,samplesheet$status) %>% as.matrix()
selection <- selection[selection[,1]>5  & selection[,2]>5,] %>% rownames()
samplesheet <- samplesheet[samplesheet$cellype %in% selection,]
counts <- counts[,samplesheet$IDs]
sums <- rowSums(counts)
exclude <- names(sums[sums == 0])
counts <- counts[!rownames(counts) %in% exclude,]


colorectalCancer_set <- list("counts" = counts, "samplesheet" = samplesheet)
saveRDS(object = colorectalCancer_set, file = "data/colorectalCancer.rds")

set <- readRDS("./data/colorectalCancer.rds")#Load processed CR dataset
samplesheet <- set$samplesheet
counts <- set$counts
samplesheet$celltype <- "all"

source("src/DDs_function_parallel.R")
source("src/DEGs_function.R")

## Execute in parallel
DDs <- DDsPar(data = counts,
              samplesheet = samplesheet,
              cells = "all",
              contrast = c("Normal","Tumor"),
              cores = 4,
              chunks = 48,
              out = "./results/CR_res_DDs.rds")

## no funciona
#DEGs <- DEGsPar(data = counts,
#                samplesheet = samplesheet,
#                cells = "all",
#                contrast = c("Normal","Tumor"),
#                test="wilcox",
#                cores=6,
#                out = "./results/CR_res_wilcox.rds",
#                norm=TRUE)

## igual, pero mas lento
DEGs_result <- DEGs(counts, samplesheet, cells = "all", contrast = c("Normal","Tumor"), test="wilcox")
saveRDS(object = DEGs_result, file = "results/CR_res_wilcox.rds")
