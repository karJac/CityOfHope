labels <- read_tsv("predicted_label.txt")
alldata$predlab <- labels$celltype
lymph <- subset(alldata, cells=WhichCells(alldata,idents=c(27,30,18,22,25,39)))
lymph$predlab <- as.factor(lymph$predlab)
lymph@active.ident <- lymph$predlab
redlymph <- subset(lymph, cells=WhichCells(lymph,idents=levels(lymph$predlab)[-c(1,2,3,5,6,11,13,14,16,18,19,20,22,23,25,27,29,30,34)]))
redlymph@reductions$pca <- NULL
redlymph@reductions$umap <- NULL
redlymph$SCT_snn_res.2.5 <- NULL
redlymph$SCT_snn_res.2 <- NULL
redlymph$SCT_snn_res.1.5 <- NULL
redlymph$SCT_snn_res.1 <- NULL
redlymph$SCT_snn_res.0.5 <- NULL
redlymph <- RunPCA(redlymph, npc=50, assay="SCT")
JackStrawPlot(redlymph, dims = 49,reduction = "pca") 
redlymph <- FindNeighbors(redlymph, dims = 1:40, k.param = 60, prune.SNN = 1/15)
redlymph <- RunUMAP(redlymph, reduction = 'pca', dims = 1:47, assay = 'SCT', 
                   reduction.name = 'umap', reduction.key = 'UMAP_', min.dist = 0.3, n.neighbors = 40, seed.use=69)
for (res in c(3)) {
     redlymph <- FindClusters(redlymph, graph.name = names(redlymph@graphs)[[2]], resolution = res, algorithm = 4) # 4 = Leiden ALGORITHM
}
redlymph$cond[which(redlymph$cond=="Naive1")] <- "Naive"
redlymph$cond[which(redlymph$cond=="Naive2")] <- "Naive"
redlymph$cond[which(redlymph$cond=="PBS1")] <- "PBS"
redlymph$cond[which(redlymph$cond=="PBS2")] <- "PBS"
redlymph$cond[which(redlymph$cond=="PBS3")] <- "PBS"
redlymph$cond[which(redlymph$cond=="PBS4")] <- "PBS"
redlymph$cond[which(redlymph$cond=="CpG1")] <- "CpG"
redlymph$cond[which(redlymph$cond=="CpG2")] <- "CpG"
redlymph$cond[which(redlymph$cond=="CpG3")] <- "CpG"
redlymph$cond[which(redlymph$cond=="CpG4")] <- "CpG"


redlymph <- subset(rl,cells=WhichCells(rl,idents=c(2,5,6,4)))
redlymph@reductions$pca <- NULL
redlymph@reductions$umap <- NULL
redlymph$SCT_snn_res.2.5 <- NULL
redlymph$SCT_snn_res.2 <- NULL
redlymph$SCT_snn_res.1.5 <- NULL
redlymph$SCT_snn_res.1 <- NULL
redlymph$SCT_snn_res.0.5 <- NULL
redlymph <- RunPCA(redlymph, npc=50, assay="SCT")
JackStrawPlot(redlymph, dims = 49,reduction = "pca") 
redlymph <- FindNeighbors(redlymph, dims = 1:40, k.param = 60, prune.SNN = 1/15)
redlymph <- RunUMAP(redlymph, reduction = 'pca', dims = 1:47, assay = 'SCT', 
                    reduction.name = 'umap', reduction.key = 'UMAP_', min.dist = 0.3, n.neighbors = 40, seed.use=69)
for (res in c(2)) {
     redlymph <- FindClusters(redlymph, graph.name = names(redlymph@graphs)[[2]], resolution = res, algorithm = 4) # 4 = Leiden ALGORITHM
}




Idents(rl) <- rl$lab
rl$lab[rl$lab=="3"] <- "NKT"
rl$lab[rl$lab=="1"] <- "Tregs"
rl$lab[rl$lab=="7"] <- "Th1 CpG"
rl$lab[rl$lab=="6"] <- "Th1 PBS"
rl$lab[rl$lab=="4"] <- "Effector/Prexhausted Cd8+"
rl$lab[rl$lab=="5"] <- "5"
rl$lab[rl$lab=="2"] <- "2"

saveRDS(rl,"rlRDS")

cpalldata$test[which(cpalldata$test=="32")] <- "cDC1 + preDC"
cpalldata$test[which(cpalldata$test=="21")] <- "cDC2"
cpalldata$test[which(cpalldata$test=="21")] <- "Oligodendrocytes"
cpalldata$test[which(cpalldata$test=="37")] <- "Neutrophiles"
cpalldata$test[which(cpalldata$test=="25")] <- "Bcells"
cpalldata$test[which(cpalldata$test=="39")] <- "Bcells"
cpalldata$test[which(cpalldata$test=="21")] <- "cDC2"
cpalldata$test[which(cpalldata$test=="21")] <- "cDC2"
cpalldata$test[which(cpalldata$test=="27")] <- "NKcells"
cpalldata$test[which(cpalldata$test=="30")] <- "Cd8+"
cpalldata$test[which(cpalldata$test=="18")] <- "Cd4+"
cpalldata$test[which(cpalldata$test=="22")] <- "Proliferating Mg & Cd8"





cpalldata$cond[which(cpalldata$cond=="Naive1")] <- "Naive"
cpalldata$cond[which(cpalldata$cond=="Naive2")] <- "Naive"
cpalldata$cond[which(cpalldata$cond=="PBS1")] <- "PBS"
cpalldata$cond[which(cpalldata$cond=="PBS2")] <- "PBS"
cpalldata$cond[which(cpalldata$cond=="PBS3")] <- "PBS"
cpalldata$cond[which(cpalldata$cond=="PBS4")] <- "PBS"
cpalldata$cond[which(cpalldata$cond=="CpG1")] <- "CpG"
cpalldata$cond[which(cpalldata$cond=="CpG2")] <- "CpG"
cpalldata$cond[which(cpalldata$cond=="CpG3")] <- "CpG"
cpalldata$cond[which(cpalldata$cond=="CpG4")] <- "CpG"


mak <- subset(cpalldata, cells= WhichCells(cpalldata,idents=c("10","4","3","11","9","14","13")))
mak@reductions$pca <- NULL
mak@reductions$umap <- NULL
mak$SCT_snn_res.2.5 <- NULL
mak$SCT_snn_res.2 <- NULL
mak$SCT_snn_res.1.5 <- NULL
mak$SCT_snn_res.1 <- NULL
mak$SCT_snn_res.0.5 <- NULL
mak <- RunPCA(mak, npc=50, assay="SCT")
mak <- FindNeighbors(mak, dims = 1:40, k.param = 60, prune.SNN = 1/15)
mak <- RunUMAP(mak, reduction = 'pca', dims = 1:47, assay = 'SCT', 
                    reduction.name = 'umap', reduction.key = 'UMAP_', min.dist = 0.3, n.neighbors = 40, seed.use=69)
for (res in c(0.5,1,1.5,2,2.5)) {
     mak <- FindClusters(mak, graph.name = names(redlymph@graphs)[[2]], resolution = res, algorithm = 4) # 4 = Leiden ALGORITHM
}





lapply(1:12,function(i){
     pie(table(mak$meh[mak$SCT_snn_res.0.5==i]), main=i)
})

#1 <- Mf
#2 <- not induced Mf
#3 <- Mf
#4 <- wszystko i nic
#5 <- Mo/Mf
#6 <- not ind Mf
#7 <- phag Mf
#8 <- cDC2
#9 <- Mo
#10 <- cDC1 + preDC
#11 <- wszystko i nic
#12 <- Bam



mak$test <- mak$SCT_snn_res.0.5
mak$test[mak$test==1] <- "Mf"
mak$test[mak$test==2] <- "niMf"
mak$test[mak$test==3] <- "Mf"
mak$test[mak$test==4] <- "?"
mak$test[mak$test==5] <- "Mo/Mf"
mak$test[mak$test==6] <- "niMF"
mak$test[mak$test==7] <- "phMf"
mak$test[mak$test==8] <- "cDC2"
mak$test[mak$test==9] <- "Mo"
mak$test[mak$test==10] <- "cDC1+preDC"
mak$test[mak$test==11] <- "??"
mak$test[mak$test==12] <- "BAMs"
mak$test[mak$SCT_snn_res.0.5==1] <- "Mf1CpG"
mak$test[mak$SCT_snn_res.0.5==3] <- "Mf2CpG"
mak$test[mak$SCT_snn_res.0.5==2] <- "niMfPBS"
mak$test[mak$SCT_snn_res.0.5==6] <- "niMf"

Idents(mak) <- mak$cond
p <- DimPlot(mak, label=T, pt.size=0.6)
Idents(mak) <- as.factor(mak$test)
DimPlot(mak, label=T, pt.size=0.6) + p


Idents(mak) <- mak$cond
DimPlot(mak, cells.highlight = WhichCells(mak,idents="PBS")) + DimPlot(mak, cells.highlight = WhichCells(mak,idents="CpG"))





comp <- readRDS("/home/rstudio/all/Opus/redo_pilot/sctalldata3.rds")


clustNum <- DimPlot(mak, label=T, pt.size=0.6)

clustNames <- DimPlot(mak, label=T, pt.size=0.6)

clustComp <- DimPlot(mak, label=T, pt.size=0.6)





############## MARKERS Microglia ################

Idents(mak) <- mak$test


### only.pos=TRUE do sciezek, do CSV only.pos=FALSE ###
mark1 <- FindMarkers(mak,"ActMg1",c(c("ActMg2","Proliferating Mg","BAMs")),only.pos = TRUE)
mark1 <- mark1[which(mark1$p_val_adj < 0.05),]
mark2 <- FindMarkers(mak,"ActMg2",c("ActMg1","Proliferating Mg","BAMs"),only.pos = TRUE)
mark2 <- mark2[which(mark2$p_val_adj < 0.05),]
mark3 <- FindMarkers(mak,"Proliferating Mg",c("ActMg1","ActMg2","BAMs"),only.pos = TRUE)
mark3 <- mark3[which(mark3$p_val_adj < 0.05),]
mark4 <- FindMarkers(mak,"BAMs",c("Proliferating Mg","ActMg1","ActMg2"),only.pos = TRUE)
mark4 <- mark4[which(mark4$p_val_adj < 0.05),]

markers <- list(mark1,mark2,mark3,mark4)

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][match(sort(markers[[i]]$avg_log2FC, decreasing = TRUE),markers[[i]]$avg_log2FC),]
}

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][c(1:25,(length(markers[[i]][[1]])-25):(length(markers[[i]][[1]]))),]
}

clusters <- c("?","niMfPBS","niMf Mg","BAMs")
for (i in 1:(length(markers))){
     markers[[i]]$cluster <- rep(clusters[[i]],dim(markers[[i]])[1])
}

for (i in 1:(length(markers))){
     markers[[i]]$gene <- rownames(markers[[i]])
}

markers <- do.call(rbind,markers)

markers <- markers[which(markers$p_val_adj < 0.05),]

write_csv(markers,"markersActMicroglia.csv")

library("rio")

myfiles <- list.files(pattern = "top50upregGenesMysteryCluster.csv")

for (i in myfiles){
     convert(i,paste(gsub('.{4}$','',i),".xlsx",sep=""))
}


#############################################################


lol <- names(table(mak$test))
for (i in lol){
     print(i)
     print(mean(mak@assays$SCT@data["Tgfb1",WhichCells(mak,idents=i)]))
     print("\n")
}


saveRDS(mak,"mak.rds")




############## MARKERS Lymphocytes ################

Idents(mak) <- mak$test


### only.pos=TRUE do sciezek, do CSV only.pos=FALSE ###
mark1 <- FindMarkers(mak,"ActMg1",c(c("ActMg2","Proliferating Mg","BAMs")),only.pos = TRUE)
mark1 <- mark1[which(mark1$p_val_adj < 0.05),]
mark2 <- FindMarkers(mak,"ActMg2",c("ActMg1","Proliferating Mg","BAMs"),only.pos = TRUE)
mark2 <- mark2[which(mark2$p_val_adj < 0.05),]
mark3 <- FindMarkers(mak,"Proliferating Mg",c("ActMg1","ActMg2","BAMs"),only.pos = TRUE)
mark3 <- mark3[which(mark3$p_val_adj < 0.05),]
mark4 <- FindMarkers(mak,"BAMs",c("Proliferating Mg","ActMg1","ActMg2"),only.pos = TRUE)
mark4 <- mark4[which(mark4$p_val_adj < 0.05),]

markers <- list(mark1,mark2,mark3,mark4)

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][match(sort(markers[[i]]$avg_log2FC, decreasing = TRUE),markers[[i]]$avg_log2FC),]
}

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][c(1:25,(length(markers[[i]][[1]])-25):(length(markers[[i]][[1]]))),]
}

clusters <- c("?","niMfPBS","niMf Mg","BAMs")
for (i in 1:(length(markers))){
     markers[[i]]$cluster <- rep(clusters[[i]],dim(markers[[i]])[1])
}

for (i in 1:(length(markers))){
     markers[[i]]$gene <- rownames(markers[[i]])
}

markers <- do.call(rbind,markers)

markers <- markers[which(markers$p_val_adj < 0.05),]

write_csv(markers,"markersActMicroglia.csv")

library("rio")

myfiles <- list.files(pattern = "markersActMicroglia.csv")

for (i in myfiles){
     convert(i,paste(gsub('.{4}$','',i),".xlsx",sep=""))
}


#############################################################


######################## MARKERS MACROPHAGES ###############################

Idents(mak) <- mak$test


### only.pos=TRUE do sciezek, do CSV only.pos=FALSE ###
mark1 <- FindMarkers(mak,"?",c(c("niMfPBS","niMf","Mo/Mf","protumor_CpG","M0_CpG","exhMf_CpG")),only.pos = FALSE)
mark1 <- mark1[which(mark1$p_val_adj < 0.05),]
mark2 <- FindMarkers(mak,"niMfPBS",c("?","niMf","Mo/Mf","protumor_CpG","M0_CpG","exhMf_CpG"),only.pos = FALSE)
mark2 <- mark2[which(mark2$p_val_adj < 0.05),]
mark3 <- FindMarkers(mak,"niMf",c("?","niMfPBS","Mo/Mf","protumor_CpG","M0_CpG","exhMf_CpG"),only.pos = FALSE)
mark3 <- mark3[which(mark3$p_val_adj < 0.05),]


markers <- list(mark1,mark2,mark3)

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][match(sort(markers[[i]]$avg_log2FC, decreasing = TRUE),markers[[i]]$avg_log2FC),]
}

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][c(1:25,(length(markers[[i]][[1]])-25):(length(markers[[i]][[1]]))),]
}

clusters <- c("?","niMfPBS","niMf")
for (i in 1:(length(markers))){
     markers[[i]]$cluster <- rep(clusters[[i]],dim(markers[[i]])[1])
}

for (i in 1:(length(markers))){
     markers[[i]]$gene <- rownames(markers[[i]])
}

markers <- do.call(rbind,markers)

markers <- markers[which(markers$p_val_adj < 0.05),]

desc <- read_csv("../../Opus/redo_pilot/makMarkDesc.csv")
desc <- as.data.frame(desc)
colnames(desc)[[2]] <- "desc"
rownames(desc) <- make.unique(desc$gene)

markers$desc <- rep(1,length(markers$gene))
for (i in 1:length(markers$gene)){
     if (markers$gene[i] %in% desc$gene){
          markers$desc[i] <- desc[i,][[2]]
     }
}


write_csv(markers,"markersMfPBS.csv")

library("rio")

myfiles <- list.files(pattern = "markersMfPBS.csv")

for (i in myfiles){
     convert(i,paste(gsub('.{4}$','',i),".xlsx",sep=""))
}





sciezki <- mclapply(markers,function(i){
     kk2 <- enrichKEGG(gene         = gen$entrezgene_id[match(i$gene,gen$external_gene_name)],
                       organism     = 'mmu',
                       pvalueCutoff = 0.05)
     return(kk2)
},mc.cores=length(markers))

for (i in 1:length(sciezki)){
     sciezki[[i]]@result <- sciezki[[i]]@result[sort(sciezki[[i]]@result$Count, decreasing = TRUE, index.return=TRUE)[[2]],]
}

myplots <- lapply(1:length(sciezki),function(i){
     dotplot(sciezki[[i]], showCategory=25) + ggtitle(clusters[i])
})


png("macrophagesPBSSciezkiKEGGSorted.png",width=1200,height=2000)
do.call(ggarrange,myplots)
dev.off()






##### PBS vs CpG #####

Idents(mak) <- mak$cond

mark1 <- FindMarkers(mak,"CpG",c(c("PBS")),only.pos = FALSE)
mark1 <- mark1[which(mark1$p_val_adj < 0.05),]


markers <- list(mark1)

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][match(sort(markers[[i]]$avg_log2FC, decreasing = TRUE),markers[[i]]$avg_log2FC),]
}

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][c(1:25,(length(markers[[i]][[1]])-25):(length(markers[[i]][[1]]))),]
}

clusters <- c("CpG")
for (i in 1:(length(markers))){
     markers[[i]]$cluster <- rep(clusters[[i]],dim(markers[[i]])[1])
}

for (i in 1:(length(markers))){
     markers[[i]]$gene <- rownames(markers[[i]])
}

markers <- do.call(rbind,markers)

markers <- markers[which(markers$p_val_adj < 0.05),]

desc <- read_csv("../../Opus/redo_pilot/makMarkDesc.csv")
desc <- as.data.frame(desc)
colnames(desc)[[2]] <- "desc"
rownames(desc) <- make.unique(desc$gene)

markers$desc <- rep(1,length(markers$gene))
for (i in 1:length(markers$gene)){
     if (markers$gene[i] %in% desc$gene){
          markers$desc[i] <- desc[i,][[2]]
     }
}


write_csv(markers,"markersMfCpgvsPBS.csv")

library("rio")

myfiles <- list.files(pattern = "markersMfCpgvsPBS.csv")

for (i in myfiles){
     convert(i,paste(gsub('.{4}$','',i),".xlsx",sep=""))
}





sciezki <- mclapply(markers,function(i){
     kk1 <- enrichKEGG(gene         = gen$entrezgene_id[match(i$gene,gen$external_gene_name)],
                       organism     = 'mmu',
                       pvalueCutoff = 0.05)
     # kk1 <- enrichGO(
     #      gene = na.omit(gen$entrezgene_id[match(i$gene, gen$external_gene_name)]), OrgDb = "org.Mm.eg.db",
     #      pvalueCutoff = 0.05, pAdjustMethod = "fdr", ont = "BP", readable = TRUE
     # )
     return(kk1)
},mc.cores=length(markers))

for (i in 1:length(sciezki)){
     sciezki[[i]]@result <- sciezki[[i]]@result[sort(sciezki[[i]]@result$Count, decreasing = TRUE, index.return=TRUE)[[2]],]
}

myplots <- lapply(1:length(sciezki),function(i){
     dotplot(sciezki[[i]], showCategory=25) + ggtitle(clusters[i])
})


png("macrophagesPBSvsCpGSciezkiKEGGSorted.png",width=600,height=1000)
do.call(ggarrange,myplots)
dev.off()


#################################




mak <- readRDS("mak.rds")

mak$cond[mak$cond=="CpG"] <- "Treated"
mak$cond[mak$cond=="PBS"] <- "Shams"
Idents(mak) <- mak$cond
mak <- subset(mak,idents=c("Shams","Treated"))

mak$test2[mak$test2=="Mystery PBS"] <- "Inert MΦ_1"
mak$test2[mak$test2=="Inert Mph"] <- "Inert MΦ_1"
mak$test2[mak$test2=="Mystery CpG"] <- "Inert MΦ_2"
mak$test2[mak$test2=="phag CpG"] <- "Phag MΦ"
mak$test2[mak$test2=="protumor CpG"] <- "Protumor MΦ"
mak$test2[mak$test2=="X PBS"] <- "UN MΦ_1"
mak$test2[mak$test2=="Y PBS"] <- "UN MΦ_2"
mak$test2[mak$test2=="tmp CpG"] <- "UN MΦ_3"
mak$test2[mak$test2=="?"] <- "UN MΦ_1"
mak$test2[mak$test2=="niMf"] <- "UN MΦ_2"
mak$test2[mak$test2=="protumor_CpG"] <- "Protumor MΦ"
mak$test2[mak$test2=="Mo/Mf"] <- "Mo/MΦ"


mak@assays$SCT@data["Arg1",mak$test2=="Shams"] %>% mean
#[1] 0.5937672
mak@assays$SCT@data["Arg1",mak$cond=="Treated"] %>% mean
#[1] 0.5015015 
mak@assays$SCT@data["Stat3",mak$test2=="Inert Mph_2"] %>% mean
#[1] 0.5937672
mak@assays$SCT@data["Stat3",mak$test2=="Inert Mph_1"] %>% mean
#[1] 0.5937672
mak@assays$SCT@data["Stat3",mak$test2=="Phag Mph"] %>% mean
#[1] 0.5937672




tmp1 <- mak$test2[mak$cond=="PBS"] %>% table
tmp2 <- mak$test2[mak$cond=="CpG"] %>% table






mak@assays$SCT@data["Cd86",mak$test2=="Inert Mph_2"] %>% mean




for (i in names(table(mak$test2))){
     print(i)
print(mak@assays$SCT@data["Stat3",mak$test2==i] %>% mean)
}



######################### DC ###########################

cdc <- subset(cpalldata, idents=c("cDC1","cDC2"))
cdc@reductions$umap@cell.embeddings[,2][cdc@reductions$umap@cell.embeddings[,2] < 0.3] <- cdc@reductions$umap@cell.embeddings[,2][cdc@reductions$umap@cell.embeddings[,2] < 0.3] + 3.2

cdc <- FindNeighbors(cdc, dims = 1:40, k.param = 60, prune.SNN = 1/15)
cdc <- FindClusters(cdc, graph.name = names(cdc@graphs)[[2]], resolution = 1.5, algorithm = 4)
#mature Cdc1 == cluster "6" -- rest- immature

##################################################



tmp <- subset(rl,idents="Treated")


#             Th1 Treated Proliferating Lymphocytes                     Tregs 
#151                        77                                           39 
#Effector2 Cd8+               Th1 Control             Effector Cd8+ 
#     73                        14                        62 
#NKT                         8 
#21                         4 

# Eff/Reg ratio:
(73+62)/39
#[1] 3.461538


tmp <- subset(rl,idents="Shams")

#              Th1 Control                       NKT            Effector2 Cd8+             Effector Cd8+ 
#                   161                           221                202                        158 
#Th1 Treated         Tregs         Proliferating Lymphocytes                8 
#   9                 286                        77                        44

# Eff/Reg ratio:
(202 + 158)/286
#1.258741






