source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")

#Get the data and process:
library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(plotly)

query.exp <- GDCquery(project = c("TCGA-LUSC","TCGA-LUAD"),
                          legacy = FALSE,
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "HTSeq - Counts",
                          access = "open")
setwd("~/data/TCGA/")
GDCdownload(query.exp)

exp.lc <- GDCprepare(query = query.exp, save = TRUE, save.filename = "lcExp.rda")

# get expression matrix
data <- assay(exp.lc)
data[1:5,1:5]
dim(data)
# get genes information
genes.info <- rowRanges(exp.lc)
genes.info
rownames(data) <- make.names(genes.info$external_gene_name,unique = T)

# get sample information
sample.info <- colData(exp.lc)
head(sample.info)
dim(sample.info)

#Create a Seurat object
Lc <- SeuratWrapper(ExpData = data, NewMeta = as.data.frame(sample.info), ProjectLabel = "LungCancerGDC", Normalize = T, scale.only.var = T, PCs = 20, dump.files = F)
Lc <- RunUMAP(object = Lc, dims = 1:10)
save(Lc, file="~/data/TCGA/LungCancerGDC.seurat.Robj")

dim(Lc@meta.data)
colnames(Lc@meta.data)
#Subset
classes <- c("Squamous cell carcinoma, NOS","Adenocarcinoma, NOS")
samps <- rownames(Lc@meta.data[which(Lc@meta.data$primary_diagnosis %in% classes),])

Lc.sub <- SubsetData(object = Lc, cells.use = samps, do.clean=T)
Lc.sub <- QuickSeurat(Lc.sub)#To add variable genes to var.genes
Lc.sub <- RunUMAP(object = Lc.sub, dims = 1:10)
save(Lc.sub, file="~/data/TCGA/LungCancerGDC.sub.seurat.Robj")
# TODO: Create the plots/data for general study statistics

DimPlot(object = Lc.sub, reduction = "umap", group.by="gender",pt.size=2)
DimPlot(object = Lc.sub, reduction = "umap", group.by="name")
DimPlot(object = Lc.sub, reduction = "umap", group.by="definition")
DimPlot(object = Lc.sub, reduction = "umap", group.by="ClassType")#--> The most informative
DimPlot(object = Lc.sub, reduction = "umap", group.by="race")
DimPlot(object = Lc.sub, reduction = "umap", group.by="ethnicity")
DimPlot(object = Lc.sub, reduction = "umap", group.by="primary_diagnosis")
DimPlot(object = Lc.sub, reduction = "umap", group.by="tumor_stage")
DimPlot(object = Lc.sub, reduction = "umap", group.by="vital_status")
DimPlot(object = Lc.sub, reduction = "umap", group.by="tissue_or_organ_of_origin")

ctp <- DimPlot(object = Lc.sub, reduction = "umap", group.by="ClassType")
ggplotly(ctp)
Lc@meta.data$cigarettes_per_day[is.na(Lc@meta.data$cigarettes_per_day)] <- 0
FeaturePlot(Lc,features.plot = "cigarettes_per_day",cols.use = c("lightblue","red"),pt.size = 2)

Lc@meta.data$years_smoked[is.na(Lc@meta.data$years_smoked)] <- 0
FeaturePlot(Lc,features.plot = "years_smoked",cols.use = c("lightblue","red"),pt.size = 2)


FeaturePlot(Lc,features.plot = "",cols.use = c("lightblue","red"),pt.size = 2)


FeaturePlot(Lc,features.plot = c("SFTPA1","TMPRSS4","MLST8","TP63","TTF1","NAPSA","PROM1", "CD44", "THY1", "ABCG2", "ALDH1A1"),cols.use = c("lightblue","red"),pt.size = 2)
#PROM1 = CD133
#THY1 = CD90
FeaturePlot(Lc,features.plot = c("CD5L", "BLACAT1", "WNT7A", "FAM107A", "SPP1", "GPM6A", "AGER", "CD36"),cols.use = c("lightblue","red"),pt.size = 2)
FeaturePlot(Lc,features.plot = c("CAV1", "HS6ST2", "HMGB3", "FABP4", "TEK", "RTKN2", "GRK5", "GCNT3", "SH3GL3", "SCUBE1"),cols.use = c("lightblue","red"),pt.size = 2)

#Novel Markers:
FeaturePlot(Lc,features.plot = "MIR205HG",cols.use = c("lightblue","red"),pt.size = 2,reduction.use = "umap")
FeaturePlot(Lc,features.plot = "DSG3",cols.use = c("lightblue","red"),pt.size = 2)
FeaturePlot(Lc,features.plot = "NECTIN1",cols.use = c("lightblue","red"),pt.size = 2)
FeaturePlot(Lc,features.plot = "CLCA2",cols.use = c("lightblue","red"),pt.size = 2)
FeaturePlot(Lc,features.plot = "DPYSL2",cols.use = c("lightblue","red"),pt.size = 2)
FeaturePlot(Lc,features.plot = "CELF2",cols.use = c("lightblue","red"),pt.size = 2)


#DE
Lc.sub <- SetAllIdent(Lc.sub, id="ClassType")
#SCC vs Normal tissues
scc.vs.norm <- FindMarkers(object = Lc.sub, ident.1 = "Lung Squamous Cell Carcinoma",
                           ident.2 = "Solid Tissue Normal",
                           logfc.threshold = 0.25,
                           min.pct = 0.25)
scc.vs.norm <- scc.vs.norm %>% add_column(gene=rownames(scc.vs.norm)) %>% filter(p_val_adj < 0.05)
dim(scc.vs.norm)

#AC vs Normal tissues
ac.vs.norm <- FindMarkers(object = Lc.sub, ident.1 = "Solid Tissue Normal",
                          ident.2 = "Lung Adenocarcinoma",
                          logfc.threshold = 0.25,
                          min.pct = 0.25)
ac.vs.norm <- ac.vs.norm %>% add_column(gene=rownames(ac.vs.norm)) %>% filter(p_val_adj < 0.05)
dim(ac.vs.norm)

#AC vs SCC
ac.vs.scc <- FindMarkers(object = Lc.sub, ident.1 = "Lung Adenocarcinoma",
                         ident.2 = "Lung Squamous Cell Carcinoma",
                         logfc.threshold = 0.25,
                         min.pct = 0.25)
ac.vs.scc <- ac.vs.scc %>% add_column(gene=rownames(ac.vs.scc)) %>% filter(p_val_adj < 0.05)
dim(ac.vs.scc)

save(scc.vs.norm, ac.vs.norm, ac.vs.scc, file="DEGs.paired.Rdata")



ac.vs.norm.sub <- ac.vs.norm %>% filter(abs(avg_logFC) > 1.5) %>% top_n(wt = avg_logFC, n = 10)
VlnPlot(object = Lc.sub, nCol = 1,single.legend = T,
        features = ac.vs.norm.sub$gene,
        ident.include = c("Lung Adenocarcinoma", "Solid Tissue Normal"))


Lc.markers <- FindAllMarkers(object = Lc.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(Lc.markers, file="~/data/TCGA/LungCancerGDC.sub.DE_testResults.Rdata")

all.markers <- Lc.markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05)
markers <- Lc.markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% top_n(n = 30, wt = avg_logFC)
VlnPlot(object = Lc.sub, features = markers$gene)

DoHeatmap(object = Lc.sub, genes.use =  markers$gene,cex.col = 0)

VlnPlot(object = Lc.sub, features =c("SFTPA1","TMPRSS4","MLST8","TP63","TTF1","NAPSA","PROM1", "CD44", "THY1", "ABCG2", "ALDH1A1"))

scc.markers <- Lc.markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05, cluster == "Lung Squamous Cell Carcinoma")
ac.markers <- Lc.markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05, cluster == "Lung Adenocarcinoma")
norm.markers <- Lc.markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05, cluster == "Solid Tissue Normal")

ggplot(scc.markers , aes(avg_logFC, -log10( p_val_adj)))+geom_point(aes(color=avg_logFC ))+scale_color_gradient(low="blue", high="red")
ggplot(ac.markers , aes(avg_logFC, -log10( p_val_adj) ))+geom_point(aes(color=avg_logFC ))+scale_color_gradient(low="blue", high="red")
ggplot(norm.markers , aes(avg_logFC, -log10( p_val_adj) ))+geom_point(aes(color=avg_logFC ))+scale_color_gradient(low="blue", high="red")

p <- ggplot(norm.markers , aes(avg_logFC, -log10( p_val_adj) ))+
  geom_point(aes(color=avg_logFC,text=sprintf("Gene: %s<br>Sample Type: %s", gene, cluster) ))+
  scale_color_gradient(low="blue", high="red")
ggplotly(p)




#Gene ontology and Pathway Enrichment Analysis:
ansEA.scc <- TCGAanalyze_EAcomplete(TFname="DEA genes of SCC", RegulonList = scc.markers$gene)
ansEA.ac <- TCGAanalyze_EAcomplete(TFname="DEA genes of Adenocarcinoma", RegulonList = ac.markers$gene)
ansEA.norm <- TCGAanalyze_EAcomplete(TFname="DEA genes of Normal Tissue", RegulonList = norm.markers$gene)

TCGAvisualize_EAbarplot(tf = rownames(ansEA.ac$ResBP),
                        filename = NULL,
                        GOBPTab = ansEA.ac$ResBP,
                        GOCCTab = ansEA.ac$ResCC,
                        GOMFTab = ansEA.ac$ResMF,
                        PathTab = ansEA.ac$ResPat,
                        nRGTab = ac.markers$gene,
                        nBar = 20)

TCGAvisualize_EAbarplot(tf = rownames(ansEA.scc$ResBP),
                        filename = NULL,
                        GOBPTab = ansEA.scc$ResBP,
                        GOCCTab = ansEA.scc$ResCC,
                        GOMFTab = ansEA.scc$ResMF,
                        PathTab = ansEA.scc$ResPat,
                        nRGTab = scc.markers$gene,
                        nBar = 20)

TCGAvisualize_EAbarplot(tf = rownames(ansEA.norm$ResBP),
                        filename = NULL,
                        GOBPTab = ansEA.norm$ResBP,
                        GOCCTab = ansEA.norm$ResCC,
                        GOMFTab = ansEA.norm$ResMF,
                        PathTab = ansEA.norm$ResPat,
                        nRGTab = ansEA.norm$gene,
                        nBar = 20)
