### Dementia vs No Dementia
#### A Script By Denis Moura
##### drmouradap@gmail.com
###### July, 2020
####################

### Description
  """
    This Script treats the Aging Dementia and Traumatic Brain Injury Study
    transcriptome dataset from Allen Institute for Brain Science
    (http://aging.brain-map.org/). It's goal is to use this data to
    create machine learning models to identify gene expression patterns
    linked to brains with dementia.
  """
#####################

### Loading libraries
#####################

library(DESeq2)
library(reshape2)
library(ggplot2)
library('RColorBrewer')
library('pheatmap')
library(e1071)
library(glue)
library(rpart)
library(data.table)
library(randomForest)
library('EnhancedVolcano')


### Loading data
#####################

setwd("/home/denis/Lab/Allen/TBI")

tbi = read.csv("fpkm_table_unnormalized.csv", header=TRUE, row.names =1)
length(tbi)
tbigenes = read.csv("rows-genes.csv")
tbisamples = read.csv("columns-samples.csv")

TCxDEGenes = read.csv('TCxDEGenes.csv')
PCxDEGenes = read.csv('PCxDEGenes.csv')
HipDEGenes = read.csv('HipDEGenes.csv')
FWMDEGenes = read.csv('FWMDEGenes.csv')
whole_brain = merge(FWMDEGenes, TCxDEGenes, all.x = T, all.y = T)
whole_brain = merge(whole_brain, HipDEGenes, all.x = T, all.y = T)
whole_brain = merge(whole_brain, PCxDEGenes, all.x = T, all.y = T)


### Treating the data
#####################

TBI = tbi
rownames(TBI) = tbigenes$gene_symbol
colnames(TBI) = tbisamples$rnaseq_profile_id

alzGroup = c('309335447', '467056409', '309335479', '309335491', '326765669', '326765685', '309335481', '326765675',
             '309335439', '309335459', '326765673', '309335471', '309335455', '309335493', '326765655', '309335477')
otherDementiaGroup = c('326765654', '326765676', '326765680', '326765670', '309335496', '309335463', '326765650',
                       '309335468')
dementiaGroup = c(alzGroup, otherDementiaGroup)

samples = tbisamples
samples$condition = 'NoDementia'

samples$condition[(which(samples$donor_id %in% dementiaGroup))] = 'Dementia'

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

TCxSamples = samples[which(samples$structure_acronym == 'TCx'),]
FWMSamples = samples[which(samples$structure_acronym == 'FWM'),]
HipSamples = samples[which(samples$structure_acronym == 'HIP'),]
PCxSamples = samples[which(samples$structure_acronym == 'PCx'),]

whole_brain = merge(FWMDEGenes, TCxDEGenes, all.x = T, all.y = T)
whole_brain = merge(whole_brain, HipDEGenes, all.x = T, all.y = T)
whole_brain = merge(whole_brain, PCxDEGenes, all.x = T, all.y = T)

### Creating functions
#####################

multiplyby100 = function(x){
  return(x*100)
}

analyze_structure = function(dataframe, sample_data) {
  df = data.frame(matrix(unlist(dataframe), nrow=50281, byrow=T),stringsAsFactors=False)
  rownames(df) = rownames(dataframe)
  colnames(df) = colnames(dataframe)
  pseudoCount = log2(df +1)
  df1 = lapply(df, multiplyby100)
  df1 = lapply(df1, as.integer)
  df1 = data.frame(matrix(unlist(df1), nrow=50281, byrow=T),stringsAsFactors=True)
  rownames(df1) = rownames(dataframe)
  colnames(df1) = colnames(dataframe)
  condition = as.factor(sample_data$condition)
  mycols = data.frame(row.names = sample_data$rnaseq_profile_id, condition)
  dds = DESeqDataSetFromMatrix(countData = df1, colData = mycols, design = ~ condition)
  
  #Normalização.
  dds = estimateSizeFactors(dds, type = 'poscounts')
  sizeFactors(dds)
  norm_counts = counts(dds, normalized = TRUE)
  pseudoCounts = log2(norm_counts + 1)
  df = melt(pseudoCount)
  df = data.frame(df, Condition = substr(df$variable, 1, 4))
  dds = DESeq(dds)
  res = results(dds)
  mcols(res, use.names = TRUE)
  print(summary(res))
  return (dds)
}

create_heatmap = function(area, dds) {
  
  res = results(dds)
  vsn = vst(dds)
  mat = assay(vsn)[head(order(res$padj),1000),]
  mat = mat - rowMeans(mat) #substracting row means from each value
  
  df = as.data.frame(colData(vsn)[,c("condition")])
  colnames(df) = "condition"
  rownames(df) = colnames(mat)
  pheatmap(mat,annotation_col=df, show_rownames = F,
           show_colnames = F, fontsize_col = 6,
           angle_col = c('45'), cluster_cols = T,
           cluster_rows = F,
           filename = glue('{area} Differentially Expressed genes in Dementia.png'))
}

create_classification_tree = function(dataframe, area) {
  set.seed(length(dataframe))
  amostra = sample(2, length(dataframe), replace = T, prob = c(0.7,0.3))
  
  training_set = dataframe[amostra==1,]
  testing_set = dataframe[amostra==2,]
  
  model = rpart(condition ~., data = dataframe, method = 'class')
  printcp(model)
  plotcp(model)
  summary(model)
  
  pfit = prune(model, cp = model$cptable[which.min(model$cptable[,'xerror']),'CP'])
  plot(pfit, uniform = T, main = 'Classification Tree for Dementia')
  text(pfit, use.n=T, all=T, cex=.8)
  
  prediction = predict(model, dataframe)

  comparison = cbind(prediction, dataframe$condition)
  
  print(comparison)
  
  test = predict(model, testing_set)
  dementia = cbind(testing_set, test)
  dementia['Result'] = ifelse(dementia$Dementia>=0.5, "Dementia", "No Dementia")
  
  confusao = table(dementia$condition, dementia$Result)
  
  write.table(confusao, file = glue('{area}_CT_Confusion.csv'), col.names = NA)
  
  success_rate = (confusao[1] + confusao[4]) / sum(confusao)
  print(glue('Model success rate: {success_rate}'))
  
  return (model)
}

plot_decision_tree = function(tree) {
  rpart.plot(tree, type = 5, roundint = F)
}

save_decision_tree = function(tree, name) {
  png(glue("{name} Decision Tree.png"), width = 800,
           height = 600)
  rpart.plot(tree, type = 5, roundint = F)
  dev.off()
}

multidimensional_scaling = function(model, dataframe) {
  distance.matrix = dist(1-model$proximity)
  
  mds.stuff = cmdscale(distance.matrix, eig=TRUE, x.ret = TRUE)
  mds.var.per = round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
  
  mds.values = mds.stuff$points
  mds.data = data.frame(Sample = rownames(mds.values),
                        X=mds.values[,1],
                        Y=mds.values[,2],
                        Status=dataframe$condition)
  
  ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
    geom_text(aes(color=Status)) +
    theme_bw() +
    xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
    ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
    ggtitle("MDS plot using (1 - Random Forest Proximities)")
}

create_random_forest = function(dataframe, area) {
  set.seed(length(dataframe))
  # Filtering features by importance
  importance = random.forest.importance(condition ~., dataframe)
  importance = setDT(importance, keep.rownames = TRUE)[]
  rf_features = importance[tail(order(importance$attr_importance),length(which(importance$attr_importance > 3)))]
  
  # Finding the optimal number of splits
  
  opt_mtry = tuneRF(x = dataframe[,-ncol(dataframe)], y = as.factor(dataframe[,ncol(dataframe)]), data = dataframe, mtryStart = 3, ntreeTry=1000, improve = 0.05, doBest = T)
  
  # Creating the model with all features
  All_Features_Model = randomForest(x = dataframe[,-ncol(dataframe)], y = as.factor(dataframe[,ncol(dataframe)]), data = dataframe, ntree= 1000, mtry = opt_mtry$mtry, proximity = TRUE)
  print(All_Features_Model$confusion)
  write.table(All_Features_Model$confusion, glue("{area}_RF_All_Confusion.csv"), sep = '\t', col.names = NA)
  
  # Creating the model with selected features
  Selected_Model = randomForest(x = dataframe[,which(importance$attr_importance > 3)], y = as.factor(dataframe[,ncol(dataframe)]), data = dataframe, ntree = 1000, mtry = opt_mtry$mtry, proximity = T)
  print(Selected_Model$confusion)
  write.table(Selected_Model$confusion, glue("{area}_RF_Selected_Confusion.csv"), sep = '\t', col.names = NA)
  
  #Testing the distance matrices for class correlation
  plot(multidimensional_scaling(All_Features_Model, dataframe))
  plot(multidimensional_scaling(Selected_Model, dataframe))
  
  return (Selected_Model)
}

### Running the scripts
#####################

#Analyzing unnormalized read-count data from ADTBI

TC_data = analyze_structure(TCxDataTest, TCxSamples)
PC_data = analyze_structure(PCxDataTest, PCxSamples)
Hyp_data = analyze_structure(HipDataTest, HipSamples)
FWM_data = analyze_structure(FWMDataTest, FWMSamples)

# Creating heatmaps for Differentially Expressed Genes

create_heatmap('Temporal Cortex', TC_data)
create_heatmap('Parietal Cortex', PC_data)
create_heatmap('Hypothalamus', Hyp_data)
create_heatmap('Frontotemporal White Matter', FWM_data)

# Creating Classification Trees Models

TC_CT_Model = create_classification_tree(TCxDEGenes, 'TC')
PC_CT_Model = create_classification_tree(PCxDEGenes, 'PC')
FWM_CT_Model = create_classification_tree(FWMDEGenes, 'FWM')
Hyp_CT_Model = create_classification_tree(HipDEGenes, 'Hip')
WB_CT_Model = create_classification_tree(whole_brain, 'WB')

# Saving Classification Trees Plots

save_decision_tree(TC_CT_Model, "Temporal Cortex")
save_decision_tree(PC_CT_Model, "Parietal Cortex")
save_decision_tree(Hyp_CT_Model, "Hippocampus")
save_decision_tree(FWM_CT_Model, "Frontotemporal White Matter")
save_decision_tree(WB_CT_Model, "Whole Brain")

# Creating Random Forest Models

TC_RF_Model = create_random_forest(TCxDEGenes, 'TC')
PC_RF_Model = create_random_forest(PCxDEGenes, 'PC')
FWM_RF_Model = create_random_forest(FWMDEGenes, "FWM")
Hip_RF_Model = create_random_forest(HipDEGenes, "Hip")

# Saving Models
save(TC_CT_Model, file = 'TC_CT_Model.RData')
save(TC_RF_Model, file = 'TC_RF_Model.RData')

save(PC_CT_Model, file = 'PC_CT_Model.RData')
save(PC_RF_Model, file = 'PC_RF_Model.RData')

save(Hyp_CT_Model, file = 'Hip_CT_Model.RData')
save(Hip_RF_Model, file = 'Hip_RF_Model.RData')

save(FWM_CT_Model, file = 'FWM_CT_Model.RData')
save(FWM_RF_Model, file = 'FWM_RF_Model.RData')

# Extracting selected features

write.csv(TC_CT_Model$variable.importance[tail(order(TC_CT_Model$variable.importance), 10)], file = 'TC_CT_Model_Features.csv')
write.csv(PC_CT_Model$variable.importance[tail(order(PC_CT_Model$variable.importance), 10)], file = 'PC_CT_Model_Features.csv')
write.csv(Hyp_CT_Model$variable.importance[tail(order(Hyp_CT_Model$variable.importance), 10)], file = 'Hyp_CT_Model_Features.csv')
write.csv(FWM_CT_Model$variable.importance[tail(order(FWM_CT_Model$variable.importance), 10)], file = 'FWM_CT_Model_Features.csv')

write.csv(TC_RF_Model$importance[tail(order(TC_RF_Model$importance), 10),], file = 'TC_RF_Model_Features.csv')
write.csv(PC_RF_Model$importance[tail(order(PC_RF_Model$importance), 10),], file = 'PC_RF_Model_Features.csv')
write.csv(Hip_RF_Model$importance[tail(order(Hip_RF_Model$importance), 10),], file = 'Hip_RF_Model_Features.csv')
write.csv(FWM_RF_Model$importance[tail(order(FWM_RF_Model$importance), 10),], file = 'FWM_RF_Model_Features.csv')


### Testing
#####################


res = results(TC_data)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')
