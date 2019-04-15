source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")

#Get the data and process:
library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(plotly)
library(gridExtra)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

load("~/data/TCGA/LungCancerGDC.sub.seurat.Robj")

p.d <- DimPlot(object = Lc.sub, reduction = "umap", group.by="ClassType",pt.size = 1.5, cols.use = c("darkgreen","red","blue"))#--> The most informative
ggsave(plot = p.d, file="UMAP-ClassType.pdf",width = 8,height = 6)

Lc.sub@meta.data %>% as_tibble()

#Random forest classifier for Major classes:
load("~/data/TCGA/Lc.rf.RF_model.Robj")

pdf("Feature_importance.pdf",width = 4, height = 7)
plot(varImp(Lc.rf, scale=FALSE,type=2),top=50)
dev.off()

#Paired DEGs
load("DEGs.paired.Rdata")
#Plots
p1 <- ggplot(scc.vs.norm, aes(avg_logFC, -log10( p_val_adj) ))+
  geom_point(aes(color=avg_logFC,text=sprintf("Gene: %s", gene) ))+
  scale_color_gradient(low="lightblue", high="red")+
  labs(title="Normal vs SCC")
ggplotly(p1)

p2 <- ggplot(ac.vs.norm, aes(avg_logFC, -log10( p_val_adj) ))+
  geom_point(aes(color=avg_logFC,text=sprintf("Gene: %s", gene) ))+
  scale_color_gradient(low="darkgreen", high="lightblue")+
  labs(title="AC vs Normal")
ggplotly(p2)

p3 <- ggplot(ac.vs.scc, aes(avg_logFC, -log10( p_val_adj) ))+
  geom_point(aes(color=avg_logFC,text=sprintf("Gene: %s", gene) ))+
  scale_color_gradient(low="red", high="darkgreen")+
  labs(title="SCC vs AC")
ggplotly(p3)

plots <- list()

plots[[1]] <- p1
plots[[2]] <- p2
plots[[3]] <- p3

ggsave(multiplot(plotlist = plots), file="DEGs_paired.plots.pdf",width = 6,height = 12)

#DE genes Violin plot
ac.vs.norm.sub <- ac.vs.norm %>% filter(abs(avg_logFC) > 1.5) %>% top_n(wt = avg_logFC, n = 10)
i=0; plots <- list()
for (gene in ac.vs.norm.sub$gene){print(gene)
  i=i+1
  p1 <- VlnPlot(object = Lc.sub, nCol = 1,single.legend = T,
        features = gene, remove.legend = T,size.y.use = 0,
        ident.include = c("Lung Adenocarcinoma", "Solid Tissue Normal"),same.y.lims = T)+
  coord_flip()
  plots[[i]] <- p1
}
ggsave(multiplot(plotlist = plots), file="DEGs_paired.AC-norm.plots.pdf",width = 6,height = 20)
