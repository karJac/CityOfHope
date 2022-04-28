library(DoubletFinder)
library(Matrix)
library(Seurat)
library(readxl)
library(tidyverse)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(parallel)
library(R.devices)
library(pbmcapply)
library(data.table)
library(glmGamPoi)

wait(10000)

source("~/all/customFunctionsKJ.R")
setwd("/home/rstudio/all/COH/pilot")

filesList <- list.files(".",pattern="h5")
projectNames <- c("Naive1","PBS1","PBS2","CpG1","CpG2","Naive2","PBS3","PBS4","CpG3","CpG4")
seuObjList <- mclapply(filesList,Read10X_h5,mc.cores = 10)
seuObjList <- mcmapply(CreateSeuratObject,seuObjList,projectNames,mc.cores = 10)


samplesListQC <- pbmclapply(seuObjList, function(sample) {
     sample <- PercentageFeatureSet(sample, "^mt-", col.name = "percent_mito")
     sample <- PercentageFeatureSet(sample, "^Rp[sl]", col.name = "percent_ribo")
     sample <- PercentageFeatureSet(sample, "^Hb[^(p)]", col.name = "percent_hb")
     sample <- PercentageFeatureSet(sample, "Pecam1|Pf4", col.name = "percent_plat")
     
     feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb", "percent_plat")
     sample$orig.ident <- Project(sample)
     
     nCountFeaturePlot <- FeatureScatter(sample, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
     
     # filter the extreme outliers to get properl violin plot
     selected_mito <- WhichCells(sample, expression = percent_mito < 15)
     sample <- subset(sample, cells = selected_mito)
     
     QCplot <- VlnPlot(sample, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
          NoLegend()
     
     
     
     exp <- sample@assays$RNA@counts
     
     # mouse_mito <- read_excel("Mouse.MitoCarta3.0.xls") # lepsze przyrownanie niż "^mt-" ale nie skonczylem go pisac
     # mouseMito <- mouse_mito[,3]
     
     
     selected_c <- WhichCells(sample, expression = nFeature_RNA > 200)
     selected_f <- rownames(sample)[Matrix::rowSums(sample) > 20]
     sample.filt <- subset(sample, features = selected_f, cells = selected_c)
     select_UMI <- WhichCells(sample.filt, expression = nCount_RNA  > 3000)
     sample.filt <- subset(sample, cells = select_UMI)
     dim(sample.filt)
     remove(sample)
     
     selected_mito <- WhichCells(sample.filt, expression = percent_mito < 4) # normaly ~10-15%
     selected_ribo <- WhichCells(sample.filt, expression = percent_ribo > 5) # normaly ~5 (cuz it correlates with high levels of mito genes)
     #tryCatch(select_over <- WhichCells(sample.filt, expression = nCount_RNA < 80000), error = function(e) NULL)
     tryCatch(select_under <- WhichCells(sample.filt, expression = nCount_RNA > 200), error = function(e) NULL)
     tryCatch(selected_hb <- WhichCells(sample.filt, expression = percent_hb < 2), error = function(e) NULL)
     tryCatch(selected_plat <- WhichCells(sample.filt, expression = percent_plat < 1), error = function(e) NULL)
     
     sample.filt <- subset(sample.filt, cells = selected_mito)
     sample.filt <- subset(sample.filt, cells = selected_ribo)
     #sample.filt <- subset(sample.filt, cells = select_over)
     sample.filt <- subset(sample.filt, cells = select_under)
     if (exists("selected_hb") == TRUE) {
          sample.filt <- subset(sample.filt, cells = selected_hb)
     }
     if (exists("selected_plat") == TRUE) {
          sample.filt <- subset(sample.filt, cells = selected_plat)
     }
     
     dim(sample.filt)
     
     table(sample.filt$orig.ident)
     
     QCfilteredPlot <- VlnPlot(sample.filt, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
          NoLegend()
     
     # Compute the relative expression of each gene per cell Use sparse matrix
     # operations, if your dataset is large, doing matrix devisions the regular way
     # will take a very long time.
     par(mar = c(4, 8, 2, 1))
     C <- sample.filt@assays$RNA@counts
     C <- Matrix::t(Matrix::t(C) / Matrix::colSums(C)) * 100
     most_expressed <- order(apply(C, 1, median), decreasing = T)[30:1]
     tmp <- data.frame(t(C[most_expressed, ]))
     tmp2 <- pivot_longer(tmp, cols = colnames(tmp))
     tmp2$name <- factor(tmp2$name, levels = colnames(tmp), ordered = TRUE)
     highestExprsPlot <- ggplot(tmp2, aes(x = value, y = name, fill = name)) +
          geom_boxplot()
     
     
     png(paste("QCplot_", Project(sample.filt), ".png", sep = ""), height = 1200, width = 1200)
     print(plot_grid(nCountFeaturePlot, highestExprsPlot,
                     nrow = 2,
                     labels = c(paste("QC for sample: ", Project(sample.filt)), "most expressed"),
                     vjust = 0.85, greedy = FALSE
     )) # improve annotations of graphs because they are overlapping with figures
     dev.off()
     return(sample.filt)
}, mc.cores = length(seuObjList))



samplesListFilt <- pbmclapply(samplesListQC, function(sample.filt) {
     
     tmp <- sample.filt
     # Filter MALAT1
     sample.filt <- sample.filt[!grepl("Malat1", rownames(sample.filt)), ]
     # Filter Mitocondrial
     # sample.filt <- sample.filt[!grepl("^MT-", rownames(sample.filt)), ]
     # Filter Ribosomal
     # sample.filt <- sample.filt[!grepl("^RP[SL][[:digit:]]", rownames(sample.filt)), ]
     # Filter Ribosomal rRNA
     sample.filt <- sample.filt[!grepl("rRNA", ignore.case = TRUE, rownames(sample.filt)), ]
     
     return(sample.filt)
}, mc.cores = length(samplesListQC))



httr::set_config(httr::config(ssl_verifypeer = FALSE))
m.s.genes <- NULL
while(is.null(m.s.genes)){  #infinite loop that will run until converhumangenelist wont crash (it crashes because of connction issues with ensmeble)
     try(m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes))
}

m.g2m.genes <- NULL
while(is.null(m.g2m.genes)){
     try(m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes))
}


samplesListNorm <- mclapply(samplesListFilt, function(sample.filt) {
     
     # Before running CellCycleScoring the data need to be normalized and
     # logtransformed.
     #sample.filt@active.assay <- "RNA"
     sample.filt <- NormalizeData(sample.filt)
     
     
     sample.filt <- CellCycleScoring(object = sample.filt, g2m.features = m.g2m.genes,
                                     s.features = m.s.genes)
     cellCycleGraph <- VlnPlot(sample.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
                               ncol = 4, pt.size = 0.1)
     
     
     png(paste("CellCycleplot_",Project(sample.filt),".png",sep=''))
     print(cellCycleGraph)
     dev.off()
     
     
     sample.filt <- FindVariableFeatures(sample.filt, verbose = F)
     sample.filt <- ScaleData(sample.filt,
                              vars.to.regress = c("nFeature_RNA", "percent_mito"),
                              verbose = F
     )
     sample.filt <- RunPCA(sample.filt, verbose = F, npcs = 40)
     sample.filt <- RunUMAP(sample.filt, dims = 1:40, verbose = F)
     
     return(sample.filt)
}, mc.cores = length(samplesListFilt))



samplesAfterDF2 <- pbmclapply(samplesListNorm, function(sample.filt) {
     sweep.res <- paramSweep_v3(sample.filt)
     sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
     bcmvn <- find.pK(sweep.stats)
     
     pK <- as.numeric(as.character(bcmvn$pK))
     BCmetric <- bcmvn$BCmetric
     pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
     
     par(mar = c(5, 4, 4, 8) + 1, cex.main = 1.2, font.main = 2)
     # plot(x = pK, y = BCmetric, pch = 16,type="b",  #visualtion of BCmetrics (but function finds the maximum automaticly)
     #     col = "blue",lty=1)
     # abline(v=pK_choose,lwd=2,col='red',lty=2)
     # title("The BCmvn distributions")
     # text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
     
     # define the expected number of doublet cellscells.
     nExp <-  round(ncol(sample.filt) * 800 / 10000) 
     sample.filt <- doubletFinder_v3(sample.filt, pN = 0.25, pK = pK_choose, nExp = nExp, PCs = 1:10)
     n <- dim(sample.filt@meta.data)[[2]] - 1
     m <- dim(sample.filt@meta.data)[[2]]
     colnames(sample.filt@meta.data)[[n]] <- "pANN"
     colnames(sample.filt@meta.data)[[m]] <- "DoubletFinder"
     
     return(sample.filt)
}, mc.cores = length(samplesListNorm))


alldata <- merge(samplesAfterDF2[[1]], c(samplesAfterDF2[2:length(samplesAfterDF2)]))
DefaultAssay(alldata) <- 'RNA'
sctransform <- SCTransform(alldata,ncells=10000,method = "glmGamPoi",vars.to.regress = "percent_mito") %>% RunPCA(npc=50, assay="SCT")
#alldata <- NormalizeData(alldata,assay = "RNA", normalization.method = "CLR", margin = 2) %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA(reduction.name = 'pca')
#JackStrawPlot(sctransform, dims = 1:45,reduction = "pca")


#saveRDS(samplesAfterDF2,"samplesAfterDF2.rds")
#saveRDS(sctransform,"sctransform.rds")
#saveRDS(alldata,"alldata.rds")



library(pheatmap)
library(enrichR)
library(rafalib)
library(clustree)
library(leiden)


sctransform@active.assay <- "SCT"

alldata <- sctransform
rm(sctransform)


alldata <- subset(alldata, cells=WhichCells(alldata,idents = "Singlet"))

alldata <- FindNeighbors(cdc, dims = 1:45, k.param = 60, prune.SNN = 1/15)
names(alldata@graphs)

# Clustering with Leiden (algorithm 4)
for (res in c(0.5,1,1.5,2,2.5)) {
     alldata <- FindClusters(alldata, graph.name = names(alldata@graphs)[[2]], resolution = res, algorithm = 1) # 4 = Leiden ALGORITHM
}

alldata <- RunUMAP(alldata, reduction = 'pca', dims = 1:47, assay = 'SCT', 
                   reduction.name = 'umap', reduction.key = 'UMAP_', min.dist = 0.3, n.neighbors = 40, seed.use=69)




######## Labels prediction with ACTINN made in python ######
# https://github.com/mafeiyang/ACTINN?fbclid=IwAR0yC7GrS11L-2whWwt76sjifX87SKpvNKgOdAvDNXG-uenEkt2PcwDOZL8

###########################################################



at_num <- matrix(as.numeric(prob),    # Convert to numeric matrix
                 ncol = ncol(prob))
colnames(at_num) <- colnames(prob)
rownames(at_num) <- rownames(prob)



maxList <- mclapply(1:dim(prob)[[2]],function(i){
     max(prob[,i])
},mc.cores=100)
maxList <- unlist(maxList)

probLowerThan07 <- prob[,sct$probMax < 0.7]
m <- sample(dim(probLowerThan07)[[2]], 5000)
lastTwo <- mclapply(m, function(i) {
     n <- 36
     return(sort(prob[, i])[c(n, n - 1)])
}, mc.cores = 100)

LTdf <- do.call(rbind,lastTwo)


mynames <- c("q","w","e","r")

mat <- matrix(c(c("q","r"),c("w","q"),c("w","q"),c("r","q")),ncol=2,byrow=TRUE)


zmat <-matrix(rep(rep(0,4),4),ncol=4) 
rownames(zmat) <- mynames
colnames(zmat) <- mynames

for ( i in 1:dim(mat)[[1]]){
     zmat[mat[i,1],mat[i,2]] <- zmat[mat[i,1],mat[i,2]]+1
}



zmat <-matrix(rep(rep(0,36),36),ncol=36) 
rownames(zmat) <- rownames(prob)
colnames(zmat) <- rownames(prob)

for ( i in 1:dim(LTdf)[[1]]){
     zmat[LTdf[i,1],LTdf[i,2]] <- zmat[LTdf[i,1],LTdf[i,2]]+1
}

zmat<-log(zmat+1)

pheatmap(zmat,cluster_rows = TRUE,cluster_cols = TRUE)



redsct <- RunPCA(redsct, verbose = F, npcs = 50)
redsct <- RunUMAP(redsct, reduction = 'pca', dims = 1:47, assay = 'SCT', reduction.name = 'umap2', reduction.key = 'UMAP_', min.dist = 0.3, n.neighbors = 40, seed.use=69)





JackStrawPlot(redsct, dims = 1:50,reduction = "pca") 


redsct <- FindNeighbors(redsct, dims = 1:40, k.param = 60, prune.SNN = 1/15)
names(redsct@graphs)


for (res in c(0.5,1,1.5,2)) {
     redsct <- FindClusters(redsct, graph.name = names(redsct@graphs)[[2]], resolution = res, algorithm = 4) # 4 = Leiden ALGORITHM
}

tmp <- LTdf[LTdf[,1] < 0.6 & LTdf[,1] > 0.4,]


lol <- mclapply(names(table(alldata$labels)),function(i){
     return(greyClusters(alldata,i,red="umap2"))
},mc.cores=36)
do.call(grid.arrange,c(lol,ncol=8))






################## scnym prediction ####################


SaveH5Seurat(alldata, filename = "train_set.h5Seurat")
Convert("train_set.h5Seurat", dest = "h5ad")
SaveH5Seurat(redsct, filename = "test_set.h5Seurat")
Convert("test_set.h5Seurat", dest = "h5ad")


########################################################
     
     
     


tmp1 <- mak$test2[mak$cond=="PBS"] %>% table
tmp2 <- mak$test2[mak$cond=="CpG"] %>% table


##### Pie charts ####

myClusters <- names(tmp1)
names(tmp1) <- rep("PBS",length(tmp1))
names(tmp2) <- rep("CpG",length(tmp2))
pie(tmp1[1],tmp2[1],main=names[1])
#### 
     
     
     

















     
     
