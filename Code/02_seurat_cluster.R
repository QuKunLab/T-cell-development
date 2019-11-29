install.packages('devtools')

# Replace '2.3.0' with your desired version
devtools::install_version(package = 'Seurat', version = package_version('2.3.0'))

library(Seurat)
library(dplyr)
library(Matrix)
rawdata=read.table('~/Desktop/single/T/10x/Figure_layout/gene/10x_2003_log2_q_norm.txt',header=TRUE,sep='\t',row.name='gene_id')
t=CreateSeuratObject(raw.data=rawdata,min.cells=3,min.genes=500,project='T cell development')
t@data = Matrix(as.matrix(rawdata), sparse = TRUE)

mito.genes=grep(pattern='mt-',x=rownames(x=t@data),value=TRUE)
percent.mito=Matrix::colSums(t@raw.data[mito.genes,])/Matrix::colSums(t@raw.data)
t=AddMetaData(object=t,metadata=percent.mito,col.name='percent.mito')
VlnPlot(object=t,features.plot=c('nGene',"nUMI",'percent.mito'),nCol= 3 )
par(mfrow=c(1,2))
GenePlot(object=t,gene1='nUMI',gene2='percent.mito')
GenePlot(object=t,gene1='nUMI',gene2='nGene')
t=FilterCells(object=t,subset.names=c('nGene','percent.mito'),low.thresholds=c(500,-Inf),high.thresholds=c(4500,0.4))

t <- FindVariableGenes(object = t, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.025, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = t@var.genes)


t <- ScaleData(object = t, vars.to.regress = c("nUMI", "percent.mito"))
t <- RunPCA(object = t, pc.genes = t@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print=5)
PCAPlot(object = t, dim.1 = 1, dim.2 = 2)
t <- ProjectPCA(object = t, do.print = FALSE)
PCHeatmap(object = t, pc.use = 1, cells.use = 50, do.balanced = TRUE, label.columns = FALSE)

t <- FindClusters(object = t, reduction.type = "pca", dims.use = 1:20,resolution=2,print.output=0,save.SNN=TRUE)
t <- RunTSNE(object = t, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = t)

