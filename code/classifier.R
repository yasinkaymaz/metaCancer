source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")


load("~/data/TCGA/LungCancerGDC.sub.seurat.Robj")
head(Lc.sub@meta.data$definition,3)
classSet <- Lc.sub@meta.data %>%
  select(shortLetterCode,definition,name) %>%
  mutate(ClassType = if_else( shortLetterCode == 'TP', as.character(name), as.character(definition) )) %>%
  mutate(ClassType = if_else( ClassType == "Recurrent Solid Tumor", "Lung Adenocarcinoma" , as.character(ClassType) )) %>%
  select(ClassType)


rownames(classSet) <- rownames(Lc.sub@meta.data)

Lc.sub@meta.data <- cbind(Lc.sub@meta.data, classSet)
head(Lc.sub@meta.data,3)
#Train a model:
Lc.rf <- CellTyperTrainer2(ExpressionData = Lc.sub@data,
                             CellLabels = Lc.sub@meta.data$ClassType,
                             run.name = "Lc.rf",
                             do.splitTest = F,
                             PCs = 5,
                             improve = F)


pdf("Feature_importance.pdf",width = 4, height = 7)
plot(varImp(Lc.rf, scale=FALSE,type=2),top=50)
dev.off()


#Tumor stage Predictors:
classSet <- NULL
  Lc.sub@meta.data %>%
  select(shortLetterCode,definition,name,tumor_stage) %>% head()


  mutate(ClassType = if_else( shortLetterCode == 'TP', as.character(name), as.character(definition) )) %>%
  mutate(ClassType = if_else( ClassType == "Recurrent Solid Tumor", "Lung Adenocarcinoma" , as.character(ClassType) )) %>%
  select(ClassType)



impSet <- varImp(Lc.rf, scale=FALSE)

impSetTop <- impSet$importance %>%
  add_column(Gene=rownames(impSet$importance)) %>%
  arrange(-`Lung Adenocarcinoma`) %>% head(10)

FeaturePlot(Lc.sub, features.plot = impSetTop$Gene,cols.use = c("lightblue","red"), pt.size = 2, reduction.use = "umap")


impSetTop <- impSet$importance %>%
  add_column(Gene=rownames(impSet$importance)) %>%
  arrange(-`Lung Squamous Cell Carcinoma`) %>% head(10)

FeaturePlot(Lc.sub, features.plot = impSetTop$Gene,cols.use = c("lightblue","red"), pt.size = 2, reduction.use = "umap")

impSetTop <- impSet$importance %>%
  add_column(Gene=rownames(impSet$importance)) %>%
  arrange(-`Solid Tissue Normal`) %>% head(10)

FeaturePlot(Lc.sub, features.plot = impSetTop$Gene,cols.use = c("lightblue","red"), pt.size = 2, reduction.use = "umap")
