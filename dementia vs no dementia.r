library(DESeq2)
library(reshape2)
library(ggplot2)
library('RColorBrewer')
library('pheatmap')
setwd("/media/bioinfo/Geral/Downloads/TBI/gene_expression_matrix_2016-03-03")
tbi = read.csv("fpkm_table_unnormalized.csv", header=TRUE, row.names =1)
#View(tbi)
tbigenes = read.csv("rows-genes.csv")
tbisamples = read.csv("columns-samples.csv")
TBI = tbi
rownames(TBI) = tbigenes$gene_symbol
colnames(TBI) = tbisamples$rnaseq_profile_id

#Incluindo coluna de condições neurodegenerativas no dataset

#noDementiaGroup = c('326765665', '326765657', '309335443', '309335475', '326765671', '326765689', '326765667', '467056405')
alzGroup = c('309335447', '467056409', '309335479', '309335491', '326765669', '326765685', '309335481', '326765675',
             '309335439', '309335459', '326765673', '309335471', '309335455', '309335493', '326765655', '309335477')
otherDementiaGroup = c('326765654', '326765676', '326765680', '326765670', '309335496', '309335463', '326765650',
                       '309335468')
dementiaGroup = c(alzGroup, otherDementiaGroup)

samples = tbisamples
samples$condition = 'NoDementia'

#samples$condition[(which(samples$donor_id %in% noDementiaGroup))] = 'No Dementia'
#samples$condition[(which(samples$donor_id %in% alzGroup))] = "Alzheimers"
#samples$condition[(which(samples$donor_id %in% otherDementiaGroup))] = 'Other Dementia'
samples$condition[(which(samples$donor_id %in% dementiaGroup))] = 'Dementia'

#which(samples$condition == '')
length(samples$condition)
unique(samples$condition)
length(which(samples$condition == 'NoDementia'))

#Selecionando apenas as amostras de interesse no estudo

#newSamples = samples[which(samples$condition == c('No Dementia', 'Alzheimers', 'Other Dementia')),]
#newSamples = newSamples[order(newSamples$condition, newSamples$structure_id),]


#View(samples)

#Genes de interesse
#goi = c('SLC20A2', 'PDGFB', 'PDGFRB', 'XPR1', 'KIAA1161', 'JAM2')

#Focando o estudo nos genes de interesse
#newDataset = TBI[which(rownames(TBI) %in% goi),which(colnames(TBI) %in% newSamples$rnaseq_profile_id)]
#newDataset = TBI[which(rownames(TBI) %in% goi),]

#View(tbisamples)
#ncol(newDataset)
#nrow(newDataset)
nrow(TBI)

df = data.frame(matrix(unlist(TBI), nrow=50281, byrow=T),stringsAsFactors=False)

rownames(df) = rownames(TBI)
colnames(df) = colnames(TBI)

pseudoCount = log2(df +1)
#boxplot(pseudoCount)
#View(pseudoCount)

multiplyby100 = function(x){
  return(x*100)
}

df1 = lapply(df, multiplyby100)
#View(df1)
df1 = lapply(df1, as.integer)

df1 = data.frame(matrix(unlist(df1), nrow=50281, byrow=T),stringsAsFactors=True)
#View(df1)
rownames(df1) = rownames(newDataset)
colnames(df1) = colnames(newDataset)
#TESTE
#View(df1)

#rownames(df1) = tbigenes$gene_symbol
#colnames(df1) = tbisamples$donor_id
#colData = as.data.frame(samples$condition)
#View(colData)

condition = as.factor(samples$condition)
mycols = data.frame(row.names = samples$rnaseq_profile_id, condition)

dds = DESeqDataSetFromMatrix(countData = df1, colData = mycols, design = ~ condition)
dds
#Funcionou... kinda.


#Removendo valores 0 ?
table(colSums(counts(dds) == 0))
idx = which.max(colSums(counts(dds) == 0))
dds2 = dds[ , -idx]

table(colSums(counts(dds2) == 0))
dds2 = dds2[, -(which.max(colSums(counts(dds) == 0)))]
dds2
#Não funcionou.

#Normalização.
dds = estimateSizeFactors(dds, type = 'poscounts')
sizeFactors(dds)

norm_counts = counts(dds, normalized = TRUE)
pseudoCounts = log2(norm_counts + 1)
df = melt(pseudoCount)
df = data.frame(df, Condition = substr(df$variable, 1, 4))
#ggplot(df, aes(x = df$variable, y = value, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count +1)))

dds = DESeq(dds)

res = results(dds)

mcols(res, use.names = TRUE)
summary(res)

'''
#MA plot
plotMA(res, ylim=c(-7,7))

resShrink = lfcShrink(dds, coef = 2)
plotMA(resShrink, ylim=c(-5,5))

#Extraindo genes diferenciados
resSig = subset(res, padj < 0.1)
resSig

#10 most upregulated genes in aboral vs oral
head(resShrink[order(resShrink$log2FoldChange),],10)

plotCounts(dds,"SLC20A2", "condition")
plotCounts(dds,"XPR1", "condition")
plotCounts(dds,"PDGFB", "condition")
plotCounts(dds,"PDGFRB", "condition")
plotCounts(dds,"KIAA1161", "condition")
plotCounts(dds,"JAM2", "condition")

#10 most down-regulated
head(resShrink[order(resShrink$log2FoldChange, decreasing=TRUE),],10)
#plotCounts(dds, "SLC20A2", "condition")

#Presenting highly significant expressed genes
head(res[order(res$padj),],5)
plotCounts(dds, "SLC20A2", "condition")
plotCounts(dds,"XPR1", "condition")
plotCounts(dds,"PDGFB", "condition")
plotCounts(dds,"PDGFRB", "condition")
plotCounts(dds,"KIAA1161", "condition")
plotCounts(dds,"JAM2", "condition")

#Filtering DE genes to fold change > (-2, 2)
resLFC = results(dds, lfcThreshold = 1)
summary(resLFC)

plotMA(resLFC, ylim=c(-7,7)) +
  abline(h = 1, col = "blue") +
  abline(h = -1, col = "blue")

resLFCShrunk = lfcShrink(dds, coef=2, res=resLFC)
plotMA(resLFCShrunk, ylim=c(-5,5))+
  abline(h = 1, col = "blue") +
  abline(h = -1, col = "blue")

'''
#Adjusting cutoff by alpha (FDR)
res.05 = results(dds, alpha=.05)
table(res.05$padj < .05)
res = res[order(res$padj),]
summary(res.05)

res.05Shrunk = lfcShrink(dds, coef=2, res=res.05)
plotMA(res.05Shrunk, alpha = 0.05, ylim=c(-5,5)) +
  abline(h = 1, col = "blue") +
  abline(h = -1, col = "blue")

#5 most significantly up-regulated genes in the aboral organ
resUpReg = res[which(res$log2FoldChange <0), ]
head(resUpReg[order(resUpReg$padj),],5)

plotCounts(dds, "SLC20A2", "condition")


#Heatmap of the significant DE genes
rld = rlogTransformation(dds)
mat = assay(rld)[head(order(res$padj),30),]
mat = mat - rowMeans(mat) #substracting row means from each value

df = as.data.frame(colData(rld)[,c("condition")])
colnames(df) = "condition"
rownames(df) = colnames(mat)
pheatmap(mat,annotation_col=df, show_rownames = T, show_colnames = T, fontsize_col = 6, angle_col = c('45'), cluster_cols = T, labels_col = samples$structure_acronym[(which(rownames(df) %in% samples$rnaseq_profile_id))], filename = "Heatmap in wholebrain.png")
#Leva muito tempo com essa quantidade de dados.

#Accessing data for a specific gene
assay(dds)["PRNP",]

#Indo além
#DE results
table(res$padj<0.05)
res <- res[order(res$padj), ]

resdata = merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)), by="row.names", sorte = FALSE)
names(resdata)[1] = 'Gene'
write.csv(resdata, "Nat dementia vs no dementia RNASeq analysis.csv")
head(resdata)
#write.csv(resdata, 'PFBC genes expression in nat dementia vs no dementia.csv')

######################################################################
############## Analisando por estrutura ##############################
######################################################################

TCxData = which(samples$structure_acronym == 'TCx')
TCxData = samples$rnaseq_profile_id[TCxData]
FWMData = which(samples$structure_acronym == 'FWM')
FWMData = samples$rnaseq_profile_id[FWMData]
HIPData = which(samples$structure_acronym == 'HIP')
HIPData = samples$rnaseq_profile_id[HIPData]
PCxData = which(samples$structure_acronym == 'PCx')
PCxData = samples$rnaseq_profile_id[PCxData]

TCxDataTest = TBI[,which(colnames(TBI) %in% TCxData)]
FWMDataTest = TBI[,which(colnames(TBI) %in% FWMData)]
HipDataTest = TBI[,which(colnames(TBI) %in% HIPData)]
PCxDataTest = TBI[,which(colnames(TBI) %in% PCxData)]

#newDataTest = newDataTest[order(newDataTest$structure),]
#View(newDataTest)

TCxSamples = samples[which(samples$structure_acronym == 'TCx'),]
FWMSamples = samples[which(samples$structure_acronym == 'FWM'),]
HipSamples = samples[which(samples$structure_acronym == 'HIP'),]
PCxSamples = samples[which(samples$structure_acronym == 'PCx'),]

ncol(TCxDataTest) == nrow(TCxSamples)
ncol(FWMDataTest) == nrow(FWMSamples)
ncol(HipDataTest) == nrow(HipSamples)
ncol(PCxDataTest) == nrow(PCxSamples)

######################################################################
############## Analisando os GOI no NeoCortex ##############################
######################################################################

df = data.frame(matrix(unlist(TCxDataTest), nrow=50281, byrow=T),stringsAsFactors=False)

rownames(df) = rownames(TCxDataTest)
colnames(df) = colnames(TCxDataTest)

pseudoCount = log2(df +1)
#boxplot(pseudoCount)
#View(pseudoCount)

df1 = lapply(df, multiplyby100)
#View(df1)
df1 = lapply(df1, as.integer)

df1 = data.frame(matrix(unlist(df1), nrow=50281, byrow=T),stringsAsFactors=True)
#View(df1)
rownames(df1) = rownames(TCxDataTest)
colnames(df1) = colnames(TCxDataTest)
#TESTE
#View(df1)

condition = as.factor(TCxSamples$condition)
mycols = data.frame(row.names = TCxSamples$rnaseq_profile_id, condition)


dds = DESeqDataSetFromMatrix(countData = df1, colData = mycols, design = ~ condition)
dds

#Removendo valores 0 ?
table(colSums(counts(dds) == 0))
idx = which.max(colSums(counts(dds) == 0))
dds2 = dds[ , -idx]

table(colSums(counts(dds2) == 0))
dds2 = dds2[, -(which.max(colSums(counts(dds) == 0)))]
dds2

#Normalização.
dds = estimateSizeFactors(dds, type = 'poscounts')
sizeFactors(dds)

norm_counts = counts(dds, normalized = TRUE)
pseudoCounts = log2(norm_counts + 1)
df = melt(pseudoCount)
df = data.frame(df, Condition = substr(df$variable, 1, 4))
#ggplot(df, aes(x = df$variable, y = value, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count +1)))

dds = DESeq(dds)

res = results(dds)

mcols(res, use.names = TRUE)
summary(res)

#write.csv(res, 'expressao genica no TCx - individuos selecionados.csv')

#Heatmap of the significant DE genes
rld = rlogTransformation(dds)
mat = assay(rld)[head(order(res$padj),30),]
mat = mat - rowMeans(mat) #substracting row means from each value

df = as.data.frame(colData(rld)[,c("condition")])
colnames(df) = "condition"
rownames(df) = colnames(mat)
pheatmap(mat,annotation_col=df, show_rownames = T, show_colnames = T, fontsize_col = 6, angle_col = c('45'), cluster_cols = T, filename = "Heatmap in TCx.png")


######################################################################
############## Analisando os GOI no Massa branca frontal (FWM) ##############################
######################################################################

df = data.frame(matrix(unlist(FWMDataTest), nrow=50281, byrow=T),stringsAsFactors=False)

rownames(df) = rownames(FWMDataTest)
colnames(df) = colnames(FWMDataTest)

pseudoCount = log2(df +1)
#boxplot(pseudoCount)
#View(pseudoCount)

df1 = lapply(df, multiplyby100)
#View(df1)
df1 = lapply(df1, as.integer)

df1 = data.frame(matrix(unlist(df1), nrow=50281, byrow=T),stringsAsFactors=True)
#View(df1)
rownames(df1) = rownames(FWMDataTest)
colnames(df1) = colnames(FWMDataTest)
#TESTE
#View(df1)

condition = as.factor(FWMSamples$condition)
mycols = data.frame(row.names = FWMSamples$rnaseq_profile_id, condition)


dds = DESeqDataSetFromMatrix(countData = df1, colData = mycols, design = ~ condition)
dds

#Removendo valores 0 ?
table(colSums(counts(dds) == 0))
idx = which.max(colSums(counts(dds) == 0))
dds2 = dds[ , -idx]

table(colSums(counts(dds2) == 0))
dds2 = dds2[, -(which.max(colSums(counts(dds) == 0)))]
dds2

#Normalização.
dds = estimateSizeFactors(dds, type = 'poscounts')
sizeFactors(dds)

norm_counts = counts(dds, normalized = TRUE)
pseudoCounts = log2(norm_counts + 1)
df = melt(pseudoCount)
df = data.frame(df, Condition = substr(df$variable, 1, 4))
#ggplot(df, aes(x = df$variable, y = value, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count +1)))

dds = DESeq(dds)

res = results(dds)

mcols(res, use.names = TRUE)
summary(res)

#write.csv(res, 'expressao genica no TCx - individuos selecionados.csv')

#Heatmap of the significant DE genes
rld = rlogTransformation(dds)
mat = assay(rld)[head(order(res$padj),30),]
mat = mat - rowMeans(mat) #substracting row means from each value

df = as.data.frame(colData(rld)[,c("condition")])
colnames(df) = "condition"
rownames(df) = colnames(mat)
pheatmap(mat,annotation_col=df, show_rownames = T, show_colnames = T, fontsize_col = 6, angle_col = c('45'), cluster_cols = T, filename = "Heatmap in FWM.png")

######################################################################
############## Analisando os GOI no Hipotálamo ##############################
######################################################################

df = data.frame(matrix(unlist(HipDataTest), nrow=50281, byrow=T),stringsAsFactors=False)

rownames(df) = rownames(HipDataTest)
colnames(df) = colnames(HipDataTest)

pseudoCount = log2(df +1)
#boxplot(pseudoCount)
#View(pseudoCount)

df1 = lapply(df, multiplyby100)
#View(df1)
df1 = lapply(df1, as.integer)

df1 = data.frame(matrix(unlist(df1), nrow=50281, byrow=T),stringsAsFactors=True)
#View(df1)
rownames(df1) = rownames(HipDataTest)
colnames(df1) = colnames(HipDataTest)
#TESTE
#View(df1)

condition = as.factor(HipSamples$condition)
mycols = data.frame(row.names = HipSamples$rnaseq_profile_id, condition)


dds = DESeqDataSetFromMatrix(countData = df1, colData = mycols, design = ~ condition)
dds

#Removendo valores 0 ?
table(colSums(counts(dds) == 0))
idx = which.max(colSums(counts(dds) == 0))
dds2 = dds[ , -idx]

table(colSums(counts(dds2) == 0))
dds2 = dds2[, -(which.max(colSums(counts(dds) == 0)))]
dds2

#Normalização.
dds = estimateSizeFactors(dds, type = 'poscounts')
sizeFactors(dds)

norm_counts = counts(dds, normalized = TRUE)
pseudoCounts = log2(norm_counts + 1)
df = melt(pseudoCount)
df = data.frame(df, Condition = substr(df$variable, 1, 4))
#ggplot(df, aes(x = df$variable, y = value, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count +1)))

dds = DESeq(dds)

res = results(dds)

mcols(res, use.names = TRUE)
summary(res)

#write.csv(res, 'expressao genica no TCx - individuos selecionados.csv')

#Heatmap of the significant DE genes
rld = rlogTransformation(dds)
mat = assay(rld)[head(order(res$padj),30),]
mat = mat - rowMeans(mat) #substracting row means from each value

df = as.data.frame(colData(rld)[,c("condition")])
colnames(df) = "condition"
rownames(df) = colnames(mat)
pheatmap(mat,annotation_col=df, show_rownames = T, show_colnames = T, fontsize_col = 6, angle_col = c('45'), cluster_cols = T, filename = "Heatmap in Hip.png")

######################################################################
############## Analisando os GOI no Córtex Parietal ##############################
######################################################################
df = data.frame(matrix(unlist(PCxDataTest), nrow=50281, byrow=T),stringsAsFactors=False)

rownames(df) = rownames(PCxDataTest)
colnames(df) = colnames(PCxDataTest)

pseudoCount = log2(df +1)
#boxplot(pseudoCount)
#View(pseudoCount)

df1 = lapply(df, multiplyby100)
#View(df1)
df1 = lapply(df1, as.integer)

df1 = data.frame(matrix(unlist(df1), nrow=50281, byrow=T),stringsAsFactors=True)
#View(df1)
rownames(df1) = rownames(PCxDataTest)
colnames(df1) = colnames(PCxDataTest)
#TESTE
#View(df1)

condition = as.factor(PCxSamples$condition)
mycols = data.frame(row.names = PCxSamples$rnaseq_profile_id, condition)


dds = DESeqDataSetFromMatrix(countData = df1, colData = mycols, design = ~ condition)
dds

#Removendo valores 0 ?
table(colSums(counts(dds) == 0))
idx = which.max(colSums(counts(dds) == 0))
dds2 = dds[ , -idx]

table(colSums(counts(dds2) == 0))
dds2 = dds2[, -(which.max(colSums(counts(dds) == 0)))]
dds2

#Normalização.
dds = estimateSizeFactors(dds, type = 'poscounts')
sizeFactors(dds)

norm_counts = counts(dds, normalized = TRUE)
pseudoCounts = log2(norm_counts + 1)
df = melt(pseudoCount)
df = data.frame(df, Condition = substr(df$variable, 1, 4))
#ggplot(df, aes(x = df$variable, y = value, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count +1)))

dds = DESeq(dds)

res = results(dds)

mcols(res, use.names = TRUE)
summary(res)

#write.csv(res, 'expressao genica no TCx - individuos selecionados.csv')

#Heatmap of the significant DE genes
rld = rlogTransformation(dds)
mat = assay(rld)[head(order(res$padj),30),]
mat = mat - rowMeans(mat) #substracting row means from each value

df = as.data.frame(colData(rld)[,c("condition")])
colnames(df) = "condition"
rownames(df) = colnames(mat)
pheatmap(mat,annotation_col=df, show_rownames = T, show_colnames = T, fontsize_col = 6, angle_col = c('45'), cluster_cols = T, filename = "Heatmap in PCx.png")