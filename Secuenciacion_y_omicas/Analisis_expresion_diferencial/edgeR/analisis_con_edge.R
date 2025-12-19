# ------------------------------------------------------------------------------
# 1. Configuración inicial
# ------------------------------------------------------------------------------
library(tximport)

BiocManager::install("edgeR")
library(edgeR)

library(ggplot2)

if(!require(ggrepel)) install.packages("ggrepel")
library(ggrepel)

library(pheatmap)
library(EnhancedVolcano)

# ------------------------------------------------------------------------------
# 2. Establecemos el directorio de trabajo
# ------------------------------------------------------------------------------

## Fijamos el directorio de trabajo en la carpeta donde tengamos los ficheros tsv y csv
setwd("C:/Users/usuario/Desktop/master_marta/SECUENCIACION_Y_OMICAS/ACTIVIDADES/ACTIVIDAD 2/mubio03_act2/Salmon")

# ------------------------------------------------------------------------------
# 3. Cargamos el diseño experimental y la traducción de los transcritos a genes
# ------------------------------------------------------------------------------

sample_info = read_csv("Design.csv", locale = locale(encoding = "UTF-8"))
tx2gene = read_tsv("Transcrito_a_Gen.tsv", col_names = FALSE)
colnames(tx2gene) = c("TXNAME", "GENEID")

sample_info
tx2gene

# ------------------------------------------------------------------------------
# 4. Definimos rutas a los ficheros de cuantificación de Salmon
# ------------------------------------------------------------------------------

## Cambiamos la ruta a la carpeta donde se encuentre alojada la carpeta denominada Salmon
setwd("C:/Users/usuario/Desktop/master_marta/SECUENCIACION_Y_OMICAS/ACTIVIDADES/ACTIVIDAD 2/mubio03_act2")

files = file.path("Salmon", sample_info$Sample, "quant.sf")
names(files) = sample_info$Sample

files

# ------------------------------------------------------------------------------
# 5. Solo conservar los archivos que existen
# ------------------------------------------------------------------------------

file.exists(files)
files_existentes = files[file.exists(files)]

# ------------------------------------------------------------------------------
#  6. Leemos los datos de expresión con tximport
# ------------------------------------------------------------------------------
txi <- tximport(files_existentes, type="salmon", tx2gene=tx2gene, countsFromAbundance="lengthScaledTPM")
txi$counts

# ------------------------------------------------------------------------------
#  7. Preparar los datos
# ------------------------------------------------------------------------------

# txi$counts ya tiene los conteos agregados a nivel de genes
counts <- round(txi$counts)  # edgeR requiere enteros

## Eliminamos las muestras que no son de interés en nuestro análisis 
sample_info
sample_info <- sample_info[-c(1,2),]
sample_info

## Nombramos las filas con el nombre de las muestras
sample_info <- as.data.frame(sample_info)
rownames(sample_info) = sample_info$Sample
head(sample_info)

## Simplificamos las columnas que van a ser de interés para el análisis
sample_info <- sample_info[, c("Condition"), drop=FALSE]
sample_info

## Simplificamos el nombre de las condiciones
sample_info$Condition[sample_info$Condition == 'Normopeso'] = 'normopeso'
sample_info$Condition[sample_info$Condition == 'Sobrepeso/Obeso2'] = 'obeso2'
sample_info

## Formateamos los valores de las condiciones para convertirlos en factor
sample_info$Condition = factor(sample_info$Condition, levels=c("normopeso", "obeso2"))
sample_info$Condition
sample_info

# ------------------------------------------------------------------------------
#  8.  Crear objeto DGEList
# ------------------------------------------------------------------------------

dge <- DGEList(counts = counts, group = sample_info$Condition)

## Normalización (TMM)
dge <- calcNormFactors(dge)

## Diseño y estimación de dispersión
design <- model.matrix(~Condition, data=sample_info)
dge <- estimateDisp(dge, design)

## Ajuste del modelo y test exacto 
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef=2)  
## coef=2 compara el segundo nivel de Condition frente al primero

## Resultados
res <- topTags(qlf, n=Inf)  # todos los genes
res_df <- as.data.frame(res)

# ------------------------------------------------------------------------------
#  9. Representaciones gráficas
# ------------------------------------------------------------------------------

## MA plot con edgeR
et <- exactTest(dge, pair = c("normopeso", "obeso2"))
res_edger <- topTags(et, n = Inf)$table
sig <- res_edger$FDR < 0.05

et$table$FDR <- p.adjust(et$table$PValue, method = "BH")
top <- order(et$table$FDR)[1:10]

plotMD(et,
       main = "MA plot: Sobrepeso/Obeso2 vs Normopeso")


text(et$table$logCPM[top],
     et$table$logFC[top],
     labels = rownames(et$table)[top],
     pos = 3,
     cex = 0.9)

## Volcano plot
EnhancedVolcano(res_df,
                lab = rownames(res_df),
                x = 'logFC',
                y = 'PValue',
                pCutoff = 0.05,
                FCcutoff = 1,
                labSize = 3,
                axisLabSize = 10)


## Heatmap de los genes diferencialmente expresados
## Seleccionamos genes con FDR < 0.05
sig_genes <- rownames(res_df)[res_df$FDR < 0.05]
mat <- counts[sig_genes, ]
mat_scaled <- t(scale(t(log2(mat + 1))))  # normalización simple para visualización


pheatmap(mat_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(sample_info),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 8,
         color = colorRampPalette(c("blue", "white", "red"))(50))

## Extraer top genes
top10_genes <- head(res_df[order(res_df$FDR), ], 10)
top10_genes

