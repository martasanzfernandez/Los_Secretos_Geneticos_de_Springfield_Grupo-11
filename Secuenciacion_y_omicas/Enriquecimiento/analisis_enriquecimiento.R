# ------------------------------------------------------------------------------
# 1. Configuración inicial
# ------------------------------------------------------------------------------
library(readr)
BiocManager::install("tximport")
library(tximport)
library(DESeq2)


# ------------------------------------------------------------------------------
# 2. Establecemos el directorio de trabajo
# ------------------------------------------------------------------------------
setwd("~/MASTER BIOINFORMÁTICA UNIR/1ER CUATRIMESTRE/Secuenciación y Ómicas de Próxima Generación/ACTIVIDADES/Actividad_2_Grupal/TallerGrupal_Ficheros")

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
files = file.path("Salmon", sample_info$Sample, "quant.sf")
names(files) = sample_info$Sample
files


# ------------------------------------------------------------------------------
# 5. Solo conservar los archivos que existen
# ------------------------------------------------------------------------------
file.exists(files)
files_existentes = files[file.exists(files)]
files_existentes


# ------------------------------------------------------------------------------
#  6. Leemos los datos de expresión con tximport
# ------------------------------------------------------------------------------
txi <- tximport(files_existentes, type="salmon", tx2gene=tx2gene, countsFromAbundance="lengthScaledTPM")
txi$counts


# ------------------------------------------------------------------------------
#  7. Ajustamos la tabla de metadatos
# ------------------------------------------------------------------------------

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
# 8. Verificación de la expresión génica
# ------------------------------------------------------------------------------

# Seleccionamos el gen de interés, por ejemplo "LEP"
gene_counts <- txi$counts["LEP", ]

# Convertimos en dataframe combinando con metadatos
gene_data <- data.frame(
  Sample = colnames(txi$counts),
  Condition = sample_info$Condition,
  Counts = as.numeric(gene_counts)
)

gene_data

# Visualizamos con ggplot las diferencias de expresión entre genotipos

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
ggplot(gene_data, aes(x = Condition, y = Counts, fill = Condition)) + geom_boxplot() 


# ------------------------------------------------------------------------------
# 9. Análisis diferencial con DESeq2
# ------------------------------------------------------------------------------

# Redondeamos los conteos a enteros. DESeq2 solo acepta enteros puros
count_data_int <- round(txi$counts)

## Creamos un objeto DESeqDataSet con el conteo y los metadatos, especificando el genotipo como diseño experimental
dds <- DESeqDataSetFromMatrix(countData=count_data_int, colData=sample_info, design=~Condition)

## En este caso, como tenemos muy pocos genes, no vamos a filtrar los genes con bajo conteo

# Ejecutamos DESeq2
dds <- DESeq(dds)

# Extraemos la información de los resultados específicos comparando las condiciones del diseño
res = results(dds, contrast=c("Condition", "obeso2", "normopeso"), alpha=0.05)

head(res)


# ------------------------------------------------------------------------------
# 10. Visualización de los resultados
# ------------------------------------------------------------------------------
## MA plot
# Definimos un umbral de significancia
alpha_threshold <- 0.05

# MA plot básico con DESeq2, resaltando en azul los genes significativos
plotMA(res, alpha = alpha_threshold, main="MA plot Sobrepeso/Obeso2 vs Normopeso")

# Opcional: añadir etiquetas a los genes con padj < 0.05 y log2FC grande
library(ggrepel) # para etiquetas más legibles

# Creamos un dataframe de genes significativos
sig_genes <- as.data.frame(res)
sig_genes$Gene <- rownames(sig_genes)
sig_genes <- sig_genes[!is.na(sig_genes$padj) & sig_genes$padj < alpha_threshold, ]

# Dibujar MA plot con ggplot2 para controlar etiquetas
if(!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)

ggplot(as.data.frame(res), aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = padj < alpha_threshold), alpha=0.6) +
  scale_x_log10() +
  scale_color_manual(values = c("black", "blue")) +
  geom_text_repel(data = sig_genes, aes(label=Gene),
                  size = 3, max.overlaps = 20) +
  theme_bw() +
  labs(title="MA plot", x="Mean of normalized counts", y="Log2 Fold Change") +
  theme(legend.position = "none") 

# Volcano plot
if(!require(EnhancedVolcano)){
  BiocManager::install("EnhancedVolcano")
  library(EnhancedVolcano)
}

# Usamos los nombres de los genes desde rownames
res$Gene.Name <- rownames(res)

EnhancedVolcano(res,
                lab=res$Gene.Name,
                x='log2FoldChange',
                y='pvalue',
                pCutoff=0.05,   # nivel de significancia
                FCcutoff=1,     # log2 fold change mínimo
                labSize = 3,
                axisLabSize = 10)


# HeatMap de todos los genes
if(!require(pheatmap)){
  install.packages("pheatmap")
  library(pheatmap)
}

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
mat <- assay(vsd)[rownames(res), ] 
rownames(mat) <- rownames(res)  # nombres de los genes
mat_scaled <- t(scale(t(mat)))  # escala por fila (genes)

pheatmap(mat_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(dds)),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 3,
         fontsize_col = 8,
         color = colorRampPalette(c("blue", "white", "red"))(50))

# Top genes
res_ordered <- res[order(res$padj), ]   # ordenar por padj
top10_genes <- head(res_ordered, 10)    # 10 genes más significativos
top10_genes

# Heatmap de los top 10 genes
mat_subset <- mat_scaled[head(order(res$padj), 10), ]
pheatmap(mat_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(dds)),
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50))


# ------------------------------------------------------------------------------
# 11. Análisis de enriquecimiento GO enrichment (ORA)
# ------------------------------------------------------------------------------

# Preparación de paquetes necesarios
pkgs_bioc <- c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot", "dplyr")
for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p)
  library(p, character.only = TRUE)
}


# Preparamos los resultados del DESqe2 para el análisis de enriquecimiento
res_df <- as.data.frame(res)
res_df$SYMBOL <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj), ] # Eliminamos NA
res_df


# Selección de genes significativos entre condiciones para GO ORA
genes_sig <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

genes_sig_symbols <- genes_sig$SYMBOL
length(genes_sig_symbols)


# Conversión de SYMBOL a ENTREZID
gene_conv <- bitr(
  genes_sig_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

genes_entrez <- gene_conv$ENTREZID
genes_entrez


# GO enrichment (ORA): Biological Process (BP)
ego_bp <- enrichGO(
  gene          = genes_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_bp)


#  Visualización de GO
dotplot(ego_bp, showCategory = 20) + 
  ggtitle("GO Biological Process enrichment") + 
  theme(axis.text.y=element_text(size=8))


# Reactome ORA
reactome_ora <- enrichPathway(
  gene          = genes_entrez,
  organism      = "human",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(reactome_ora)


#  Visualización de Reactome ORA
dotplot(reactome_ora, showCategory = 20) +
  ggtitle("Reactome pathway enrichment (ORA)") +
  theme(axis.text.y=element_text(size=6))

