if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("apeglm"))

rm(list = ls())
library(tidyverse) 
library(DESeq2)

#Carga de datos
dds_GSEA <- readRDS("Apartado2/input/dds2.rds")

# DEGs
res <- results(dds_GSEA, alpha = 0.05, contrast = c("group", "DPN24h", "Control24h"))
summary(res)

# Shrunken LFC
res.ape <- lfcShrink(dds = dds_GSEA, coef = "group_DPN24h_vs_Control24h", type = "apeglm",
                     res = res)
summary(res.ape)

#Creación del .rnk
rnk <- data.frame(Feature = rownames(res.ape), LFC = res.ape$log2FoldChange)

#Guardar el .rnk
write.table(rnk, file = "Apartado2/input/DPN-Control_24h.rnk", sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)
