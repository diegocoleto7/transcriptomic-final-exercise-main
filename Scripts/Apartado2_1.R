if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "tidyverse", "pheatmap", "RColorBrewer"))

rm(list = ls())

library(DESeq2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)


# Cargar datos
counts_data <- read.csv(file = "Apartado2/input/rawcounts.tsv", sep = "\t", row.names = 1)
experiment_data <- read.csv(file = "Apartado2/input/metadata.tsv", sep = "\t", row.names = 1)

# Preprocesamiento de datos (factorización)
experiment_data$patient<-factor(experiment_data$patient)
experiment_data$agent<-factor(experiment_data$agent)
experiment_data$time<-factor(experiment_data$time)

# Crear variable 'group'
experiment_data$group <- as.factor(paste0(experiment_data$agent, experiment_data$time))

#Comprobación
all(rownames(experiment_data) == colnames(counts_data))
all(rownames(experiment_data) %in% colnames(counts_data))

# Crear DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = experiment_data,
                              design = ~ patient + group)
summary(dds)
# Filtrar genes con bajo recuento
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
summary(dds)

# Transformación VST
vsd <- vst(dds, blind = TRUE)

# PCA
plotPCA(vsd, intgroup = c("patient", "group", "agent", "time"))

#Eliminar outlayer
experiment_data = subset(experiment_data, patient != 4)
counts_data <- counts_data[, rownames(experiment_data)]

#Comprobación
all(rownames(experiment_data) == colnames(counts_data))
all(rownames(experiment_data) %in% colnames(counts_data))

# Crear DESeqDataSet de nuevo
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = experiment_data,
                              design = ~ patient + group)
summary(dds)

# Filtrar genes con bajo recuento
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
summary(dds)

# Transformación VST
vsd <- vst(dds, blind = TRUE)

# PCA
plotPCA(vsd, intgroup = c("patient", "group", "agent", "time"))

# Matriz de distancias
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$patient, vsd$group, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# Análisis de expresión diferencial
dds2 <- DESeq(dds, test = "Wald")
plotDispEsts(dds2)

# Resultados DPN
res_DPN <- results(dds2,
                   contrast = c("group", "Control24h", "DPN24h"),
                   alpha = 0.05,
                   pAdjustMethod = "BH")
summary(res_DPN)
#Probamos la seleccion de un corte de logfold para el primer análisis
my_results_DPN <- results(object = dds2,
                          contrast = c("group", "Control24h", "DPN24h"),
                          alpha = 0.05,
                          lfcThreshold = 1,
                          pAdjustMethod = "BH")
summary(my_results_DPN)
# Resultados OHT
res_OHT <- results(dds2,
                   contrast = c("group", "Control24h", "OHT24h"),
                   alpha = 0.05,
                   pAdjustMethod = "BH")

summary(res_OHT)

#Probamos la seleccion de un corte de logfold para el segundo análisis
my_results_OHT <- results(object = dds2,
                          contrast = c("group", "Control24h", "OHT24h"),
                          alpha = 0.05,
                          lfcThreshold = 1,
                          pAdjustMethod = "BH")
summary(my_results_OHT)

#Visualización de los resultados
mat_DPN <- assay(vsd)[head(order(my_results_DPN$padj), 30), ]
mat_DPN <- mat_DPN - rowMeans(mat_DPN)
anno <- as.data.frame(colData(vsd)[ , c("patient", "agent", "time")])
ann_colors <- list(
    patient = c("1" ="red", "2"= "yellow", "3" = "blue"),
    agent = c(Control ="gold1", DPN ="turquoise2", OHT = "purple"),
    time = c("24h" = "red2", "48h" ="aquamarine")
)

mat_OHT <- assay(vsd)[head(order(my_results_OHT$padj), 30), ]
mat_OHT <- mat_OHT - rowMeans(mat_OHT)
anno <- as.data.frame(colData(vsd)[ , c("patient", "agent", "time")])
ann_colors <- list(
    patient = c("1" ="red", "2"= "yellow", "3" = "blue"),
    agent = c(Control ="gold1", DPN ="turquoise2", OHT = "purple"),
    time = c("24h" = "red2", "48h" ="aquamarine")
)

par(mfrow = c(1, 2))
heatmap_DPN <- pheatmap(mat = mat_DPN, annotation_col = anno, show_colnames = F, annotation_colors = ann_colors, main = "Heatmap Control vs DPN -24 hours")
heatmap_OHT <- pheatmap(mat = mat_OHT, annotation_col = anno, show_colnames = F, annotation_colors = ann_colors, main = "Heatmap Control vs OHT -24 hours")

#Guardar los resutados
resOrdered_1 <- subset(my_results_DPN[order(my_results_DPN$padj),], padj < 0.05)
resOrdered_2 <- subset(my_results_OHT[order(my_results_OHT$padj),], padj < 0.05)

dir.create(path = "Apartado2/out")

write.csv(resOrdered_1, file = "Apartado2/out/res_Ctrl-DPN_24h.csv")
write.csv(resOrdered_2, file = "Apartado2/out/res_Ctrl-OHT_24h.csv")

saveRDS(object = dds2, file = "Apartado2/input/dds2.rds")
