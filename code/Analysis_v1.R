library(tidyverse)
#wget https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.4.0_OSX_x64_10.12.6.zip
#unzip https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.4.0_OSX_x64_10.12.6.zip
#./gdc-client download -m gdc_manifest_20190408_194451.txt
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")


sample <- read.delim("data/sample.tsv",header=T)
exposure <- read.delim("data/exposure.tsv", header = T)
portion <- read.delim("data/portion.tsv", header = T)
slide <- read.delim("data/slide.tsv", header=T)
analyte <- read.delim("data/analyte.tsv",header = T)
clinical <- read.delim("data/clinical.tsv", header=T)
samplesheet <- read.delim("data/gdc_sample_sheet.2019-04-08.tsv", header = T)
manifest <- read.delim("data/gdc_manifest_20190408_194451.txt", header = T)

comb <- samplesheet %>% left_join(clinical,by=c("Case.ID" = "submitter_id") ) %>% left_join(exposure, by=c("Case.ID" = "submitter_id"))


rownames(comb) <- make.names(comb$Sample.ID,unique = T)
comb$Sample.ID <- make.names(comb$Sample.ID,unique = T)

exp <- read.delim(paste("~/data/TCGA/data/",comb[1,]$File.ID,"/",comb[1,]$File.Name,sep = ""), header=F,row.names = 1)
colnames(exp) <- c(as.character(comb[1,]$Sample.ID))
head(exp)
for(i in 2:length(comb[,1])){print(i)
s2 <- read.delim(paste("~/data/TCGA/data/",comb[i,]$File.ID,"/",comb[i,]$File.Name,sep = ""), header=F,row.names = 1)
colnames(s2) <- c(as.character(comb[i,]$Sample.ID))
exp <- cbind(exp,s2)
}


save(exp, file="data/expression.FPKM-UQ.txt")


Lc <- SeuratWrapper(ExpData = exp,NewMeta = comb, ProjectLabel = "LungCancer", Normalize = F, scale.only.var = T, PCs = 20, dump.files = F)
head(Lc@meta.data)
TSNEPlot(Lc, group.by="Sample.Type")
TSNEPlot(Lc, group.by="gender")
TSNEPlot(Lc, group.by="race")
TSNEPlot(Lc, group.by="ethnicity")
TSNEPlot(Lc, group.by="primary_diagnosis")
TSNEPlot(Lc, group.by="tumor_stage")
TSNEPlot(Lc, group.by="vital_status")
TSNEPlot(Lc, group.by="tissue_or_organ_of_origin")
TSNEPlot(Lc, group.by="project_id.y")
TSNEPlot(Lc, group.by="res.1")

FeaturePlot(Lc,features.plot = "years_smoked")

#Study cohort infographics
Lc@meta.data %>% ggplot(aes(x=primary_diagnosis,color=tumor_stage))+geom_bar()







#DE
load("~/data/TCGA/LungCancer.seurat.Robj")
Lc <- SeuratObj
rm(SeuratObj)

Lc = SetAllIdent(Lc, id = 'res.1')
ers <- FindAllMarkers(Lc, only.pos = T,logfc.threshold=.1)
head(Lc@meta.data)
Lc <- RunUMAP(object = Lc, dims = 1:10)
DimPlot(object = Lc, reduction = "umap")
