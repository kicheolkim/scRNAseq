### differential expression test using Seurat
library(xlsx); library(multcomp)


#####  add time + celltype
immune.cells <- SetAllIdent(immune.cells, id = "res.2")    # select resolution of tree
table(immune.cells@meta.data$batchid, immune.cells@meta.data$res.2)
immune.cells <- RenameIdent(
  object = immune.cells,
  old.ident.name = c('0'),
  new.ident.name = 'Mp1'
)
immune.cells <- RenameIdent(
  object = immune.cells,
  old.ident.name = c('1'),
  new.ident.name = 'Mg1'
)
immune.cells <- RenameIdent(
  object = immune.cells,
  old.ident.name = c('2'),
  new.ident.name = 'Mg2'
)
immune.cells <- RenameIdent(
  object = immune.cells,
  old.ident.name = c('3'),
  new.ident.name = 'Mg3'
)
immune.cells <- RenameIdent(
  object = immune.cells,
  old.ident.name = c('4'),
  new.ident.name = 'Mp2'
)
immune.cells <- RenameIdent(
  object = immune.cells,
  old.ident.name = c('5'),
  new.ident.name = 'Mg4'
)
immune.cells <- RenameIdent(
  object = immune.cells,
  old.ident.name = c('6'),
  new.ident.name = 'Mg5'
)
immune.cells <- RenameIdent(
  object = immune.cells,
  old.ident.name = c('7'),
  new.ident.name = 'Nk'
)
immune.cells <- RenameIdent(
  object = immune.cells,
  old.ident.name = c('8'),
  new.ident.name = 'Mg6'
)
immune.cells <- RenameIdent(
  object = immune.cells,
  old.ident.name = c('9'),
  new.ident.name = 'Mc'
)
table(immune.cells@meta.data$batchid, immune.cells@ident)
samplecluster = paste(immune.cells@ident, immune.cells@meta.data$batchid, sep = '_')
names(samplecluster) = rownames(immune.cells@meta.data)
head(samplecluster)

immune.cells <- AddMetaData(
  object = immune.cells,
  metadata = samplecluster,
  col.name = "TimeCell2")
table(immune.cells@meta.data$TimeCell)

save(immune.cells,file="Seurat_object_immune_final clusters 2.RData")


################################################################################################################################
library(Seurat)
immune.cells <- SetAllIdent(immune.cells, id = "finalcluster")
TSNEPlot(object = immune.cells, pt.size=1.5, do.label = T)

#####  DEG between Macrophage1 and Macrophage2
setwd("C:/Users/kicheol/Desktop/Projects/Andres_scRNA/Results_Seurat_new_NoCorrection/ImmuneCells/DEG/Mp1-Mp2 - clusters")
comp_deseq2 <- FindMarkers(immune.cells, ident.1 = "Macrophage1", ident.2= "Macrophage2", test.use="DESeq2")
comp_mast <- FindMarkers(immune.cells, ident.1 = "Macrophage1", ident.2= "Macrophage2", test.use="MAST")
comp_wilcox <- FindMarkers(immune.cells, ident.1 = "Macrophage1", ident.2= "Macrophage2")
write.csv(comp_deseq2, "deg_deseq2.csv")
write.csv(comp_mast, "deg_mast.csv")
write.csv(comp_wilcox, "deg_wilcox.csv")

geneList <- intersect(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_mast, p_val_adj<0.05)))
write.csv(data.frame(geneList), "significant_gene_list_overlap.csv")

pdf("DEG_heatmap-Mp1_Mp2.pdf", paper='letter')
DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$finalcluster %in% c("Macrophage1", "Macrophage2")],
          genes.use = geneList,
          group.label.rot=TRUE, group.cex = 10, cex.row = 6, slim.col.label = TRUE)
dev.off()

DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$finalcluster %in% c("Macrophage1", "Macrophage2")],
          genes.use = row.names(subset(comp_wilcox, p_val_adj<0.05)),
          group.label.rot=FALSE, group.cex = 10, cex.row = 4, slim.col.label = TRUE)



### Venn diagram
library("VennDiagram")
library("RColorBrewer")
venn <- list("deseq2"=row.names(subset(comp_deseq2, p_val_adj<0.05)), 
             "mast"=row.names(subset(comp_mast, p_val_adj<0.05)),
             "wilcox"=row.names(subset(comp_wilcox, p_val_adj<0.05)))
grid.newpage()
venn.plot <- venn.diagram(venn, NULL, resolution=600, fill=brewer.pal(3, "Dark2"), col="transparent", alpha=c(0.5,0.5,0.5), 
                          main.cex = 1.5, sub.cex=2, cat.cex=1,
                          category=c("deseq2","mast","wilcox"), main="(p<0.05)")
grid.draw(venn.plot)




#############################################################################################################################
### set the identity to the new variable 
immune.cells <- SetAllIdent(immune.cells, id = "TimeCell")
TSNEPlot(object = immune.cells, pt.size=1.5, do.label = T)
table(immune.cells@meta.data$TimeCell)

#####  DEG among time-point
library(dplyr)
## using DE based on a model using the negative binomial distribution (DESeq2)

### Microglia
setwd("C:/Users/kicheol/Desktop/Projects/Andres_scRNA/Results_Seurat_new_NoCorrection/ImmuneCells/DEG/Microglia_all - Time")
#  Mg_1_Baseline, Mg_2_PreOnset, Mg_3_Chronic
comp_deseq1 <- FindMarkers(immune.cells, ident.1 = "Mg_1_Baseline", ident.2= "Mg_2_PreOnset", test.use="DESeq2")
comp_deseq2 <- FindMarkers(immune.cells, ident.1 = "Mg_1_Baseline", ident.2= "Mg_3_Chronic", test.use="DESeq2")
comp_deseq3 <- FindMarkers(immune.cells, ident.1 = "Mg_2_PreOnset", ident.2= "Mg_3_Chronic", test.use="DESeq2")
write.csv(comp_deseq1, "comp_Mg-1Base_2Pre-DESeq2.csv")
write.csv(comp_deseq2, "comp_Mg-1Base_3Chro-DESeq2.csv")
write.csv(comp_deseq3, "comp_Mg-2Pre_3Chro-DESeq2.csv")

comp_wilcox1 <- FindMarkers(immune.cells, ident.1 = "Mg_1_Baseline", ident.2= "Mg_2_PreOnset", logfc.threshold = 0.1)
comp_wilcox2 <- FindMarkers(immune.cells, ident.1 = "Mg_1_Baseline", ident.2= "Mg_3_Chronic", logfc.threshold = 0.1)
comp_wilcox3 <- FindMarkers(immune.cells, ident.1 = "Mg_2_PreOnset", ident.2= "Mg_3_Chronic", logfc.threshold = 0.1)
write.csv(comp_wilcox1, "comp_Mg-1Base_2Pre-wilcox.csv")
write.csv(comp_wilcox2, "comp_Mg-1Base_3Chro-wilcox.csv")
write.csv(comp_wilcox3, "comp_Mg-2Pre_3Chro-wilcox.csv")

comp1 <- intersect(row.names(subset(comp_deseq1, p_val_adj<0.05)),row.names(subset(comp_mast1, p_val_adj<0.05)))
comp2 <- intersect(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_mast2, p_val_adj<0.05)))
comp3 <- intersect(row.names(subset(comp_deseq3, p_val_adj<0.05)),row.names(subset(comp_mast3, p_val_adj<0.05)))
comp <- c(comp1, comp2, comp3)

comp1 <- row.names(subset(comp_wilcox1, p_val_adj<0.05))
comp2 <- row.names(subset(comp_wilcox2, p_val_adj<0.05))
comp3 <- row.names(subset(comp_wilcox3, p_val_adj<0.05))
comp <- c(comp1, comp2, comp3)

# heatmap
pdf("DEG_heatmap-Microglia.pdf", paper='letter')
DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$TimeCell %in% c("Mg_1_Baseline", "Mg_2_PreOnset", "Mg_3_Chronic")],
          genes.use = comp,
          group.label.rot=FALSE, group.cex = 10, cex.row = 6, slim.col.label = TRUE)
dev.off()
rm(comp_deseq1, comp_deseq2, comp_deseq3, comp_mast1, comp_mast2, comp_mast3)

# boxplot
imm_data <- as.data.frame(as.matrix(immune.cells@data))
imm_meta <- as.data.frame(immune.cells@meta.data)

tmp_meta <- subset(imm_meta, imm_meta$finalcluster == "Microglia")
tmp_meta <- tmp_meta[order(tmp_meta$batchid),]
tmp_data <- imm_data[,colnames(imm_data) %in% row.names(tmp_meta)]

for (i in 1:length(comp)){
  tmp_plot <- data.frame(grp=tmp_meta$batchid, expr=t(tmp_data[row.names(tmp_data)==comp[i],]))
  colnames(tmp_plot) <- c("grp","expr")
  png(filename=paste0(comp[i],"_gene.png"))
    gp <- ggplot(tmp_plot, aes(grp, expr)) + 
      geom_boxplot(outlier.colour="Black", outlier.shape=1, aes(colour=grp), show.legend = FALSE) +
      geom_jitter(width=0.2, shape = 16, size=1, aes(colour = grp), show.legend = FALSE) +
      ggtitle(paste0("Microglia: ", comp[i]," gene")) +                              ##### name of cell type
      xlab("") + ylab("Expression level")
    print(gp)
  dev.off()
}



### Macrophage 1
#  Mp1_1_Baseline, Mp1_2_PreOnset, Mp1_3_Chronic
comp_deseq1 <- FindMarkers(immune.cells, ident.1 = "Mp1_1_Baseline", ident.2= "Mp1_2_PreOnset", test.use="DESeq2")
comp_deseq2 <- FindMarkers(immune.cells, ident.1 = "Mp1_1_Baseline", ident.2= "Mp1_3_Chronic", test.use="DESeq2")
comp_deseq3 <- FindMarkers(immune.cells, ident.1 = "Mp1_2_PreOnset", ident.2= "Mp1_3_Chronic", test.use="DESeq2")
write.csv(comp_deseq1, "comp_Mp1-1Base_2Pre-DESeq2.csv")
write.csv(comp_deseq2, "comp_Mp1-1Base_3Chro-DESeq2.csv")
write.csv(comp_deseq3, "comp_Mp1-2Pre_3Chro-DESeq2.csv")

comp_wilcox1 <- FindMarkers(immune.cells, ident.1 = "Mp1_1_Baseline", ident.2= "Mp1_2_PreOnset", logfc.threshold = 0.1)
comp_wilcox2 <- FindMarkers(immune.cells, ident.1 = "Mp1_1_Baseline", ident.2= "Mp1_3_Chronic", logfc.threshold = 0.1)
comp_wilcox3 <- FindMarkers(immune.cells, ident.1 = "Mp1_2_PreOnset", ident.2= "Mp1_3_Chronic", logfc.threshold = 0.1)
write.csv(comp_wilcox1, "comp_Mp1-1Base_2Pre-wilcox.csv")
write.csv(comp_wilcox2, "comp_Mp1-1Base_3Chro-wilcox.csv")
write.csv(comp_wilcox3, "comp_Mp1-2Pre_3Chro-wilcox.csv")

comp1 <- intersect(row.names(subset(comp_deseq1, p_val_adj<0.05)),row.names(subset(comp_mast1, p_val_adj<0.05)))
comp2 <- intersect(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_mast2, p_val_adj<0.05)))
comp3 <- intersect(row.names(subset(comp_deseq3, p_val_adj<0.05)),row.names(subset(comp_mast3, p_val_adj<0.05)))
comp <- c("Hbb-bs","Pcp2","Ccl4")

comp1 <- row.names(subset(comp_wilcox1, p_val_adj<0.05))
comp2 <- row.names(subset(comp_wilcox2, p_val_adj<0.05))
comp3 <- row.names(subset(comp_wilcox3, p_val_adj<0.05))
comp <- c(comp1, comp2, comp3)


# heatmap
pdf("DEG_heatmap-Macrophage1.pdf", paper='letter')
DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$TimeCell %in% c("Mp1_1_Baseline", "Mp1_2_PreOnset", "Mp1_3_Chronic")],
          genes.use = comp,
          group.label.rot=FALSE, group.cex = 10, cex.row = 8, slim.col.label = TRUE)
dev.off()

# boxplot
imm_data <- as.data.frame(as.matrix(immune.cells@data))
imm_meta <- as.data.frame(immune.cells@meta.data)

tmp_meta <- subset(imm_meta, imm_meta$finalcluster == "Macrophage1")
tmp_meta <- tmp_meta[order(tmp_meta$batchid),]
tmp_data <- imm_data[,colnames(imm_data) %in% row.names(tmp_meta)]

for (i in 1:length(comp)){
  tmp_plot <- data.frame(grp=tmp_meta$batchid, expr=t(tmp_data[row.names(tmp_data)==comp[i],]))
  colnames(tmp_plot) <- c("grp","expr")
  png(filename=paste0(comp[i],"_gene.png"))
  gp <- ggplot(tmp_plot, aes(grp, expr)) + 
    geom_boxplot(outlier.colour="Black", outlier.shape=1, aes(colour=grp), show.legend = FALSE) +
    geom_jitter(width=0.2, shape = 16, size=1, aes(colour = grp), show.legend = FALSE) +
    ggtitle(paste0("Macrophage1: ", comp[i]," gene")) +                              ##### name of cell type
    xlab("") + ylab("Expression level")
  print(gp)
  dev.off()
}

rm(comp_deseq1, comp_deseq2, comp_deseq3, comp_mast1, comp_mast2, comp_mast3)


## Macrophage 2
#  Mp2_1_Baseline, Mp2_2_PreOnset, Mp2_3_Chronic
comp_deseq1 <- FindMarkers(immune.cells, ident.1 = "Mp2_1_Baseline", ident.2= "Mp2_2_PreOnset", test.use="DESeq2")
comp_deseq2 <- FindMarkers(immune.cells, ident.1 = "Mp2_1_Baseline", ident.2= "Mp2_3_Chronic", test.use="DESeq2")
comp_deseq3 <- FindMarkers(immune.cells, ident.1 = "Mp2_2_PreOnset", ident.2= "Mp2_3_Chronic", test.use="DESeq2")
write.csv(comp_deseq1, "comp_Mp2-1Base_2Pre-DESeq2.csv")
write.csv(comp_deseq2, "comp_Mp2-1Base_3Chro-DESeq2.csv")
write.csv(comp_deseq3, "comp_Mp2-2Pre_3Chro-DESeq2.csv")

comp_wilcox1 <- FindMarkers(immune.cells, ident.1 = "Mp2_1_Baseline", ident.2= "Mp2_2_PreOnset", logfc.threshold = 0.1)
comp_wilcox2 <- FindMarkers(immune.cells, ident.1 = "Mp2_1_Baseline", ident.2= "Mp2_3_Chronic", logfc.threshold = 0.1)
comp_wilcox3 <- FindMarkers(immune.cells, ident.1 = "Mp2_2_PreOnset", ident.2= "Mp2_3_Chronic", logfc.threshold = 0.1)
write.csv(comp_wilcox1, "comp_Mp2-1Base_2Pre-wilcox.csv")
write.csv(comp_wilcox2, "comp_Mp2-1Base_3Chro-wilcox.csv")
write.csv(comp_wilcox3, "comp_Mp2-2Pre_3Chro-wilcox.csv")

comp1 <- intersect(row.names(subset(comp_deseq1, p_val_adj<0.05)),row.names(subset(comp_mast1, p_val_adj<0.05)))
comp2 <- intersect(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_mast2, p_val_adj<0.05)))
comp3 <- intersect(row.names(subset(comp_deseq3, p_val_adj<0.05)),row.names(subset(comp_mast3, p_val_adj<0.05)))
#comp <- c(comp1, comp2, comp3)
comp <- c("Hbb-bs","Hba-a1","Pcp2","Pde6h","Rho","Hspa1a","Opn1sw","mt-Co3","mt-Co1")

comp1 <- row.names(subset(comp_wilcox1, p_val_adj<0.05))
comp2 <- row.names(subset(comp_wilcox2, p_val_adj<0.05))
comp3 <- row.names(subset(comp_wilcox3, p_val_adj<0.05))
comp <- c(comp1, comp2, comp3)

# heatmap
pdf("DEG_heatmap-Macrophage2.pdf", paper='letter')
DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$TimeCell %in% c("Mp2_1_Baseline", "Mp2_2_PreOnset", "Mp2_3_Chronic")],
          genes.use = comp,
          group.label.rot=TRUE, group.cex = 10, cex.row = 8, slim.col.label = TRUE)
dev.off()

# boxplot
imm_data <- as.data.frame(as.matrix(immune.cells@data))
imm_meta <- as.data.frame(immune.cells@meta.data)

tmp_meta <- subset(imm_meta, imm_meta$finalcluster == "Macrophage2")
tmp_meta <- tmp_meta[order(tmp_meta$batchid),]
tmp_data <- imm_data[,colnames(imm_data) %in% row.names(tmp_meta)]

for (i in 1:length(comp)){
  tmp_plot <- data.frame(grp=tmp_meta$batchid, expr=t(tmp_data[row.names(tmp_data)==comp[i],]))
  colnames(tmp_plot) <- c("grp","expr")
  png(filename=paste0(comp[i],"_gene.png"))
  gp <- ggplot(tmp_plot, aes(grp, expr)) + 
    geom_boxplot(outlier.colour="Black", outlier.shape=1, aes(colour=grp), show.legend = FALSE) +
    geom_jitter(width=0.2, shape = 16, size=1, aes(colour = grp), show.legend = FALSE) +
    ggtitle(paste0("Macrophage2: ", comp[i]," gene")) +                              ##### name of cell type
    xlab("") + ylab("Expression level")
  print(gp)
  dev.off()
}

rm(comp_deseq1, comp_deseq2, comp_deseq3, comp_mast1, comp_mast2, comp_mast3)




##################  Within Microglia  ###############################################################################################
library(Seurat)
setwd("C:/Users/kicheol/Desktop/Projects/Andres_scRNA/Results_Seurat_new_NoCorrection/ImmuneCells/DEG/Microglia - clusters")
immune.cells <- SetAllIdent(immune.cells, id = "res.2")
TSNEPlot(object = immune.cells, pt.size=1.5, do.label = T)


#####  DEG between Microglia 1 and 2
setwd("C:/Users/kicheol/Desktop/Projects/Andres_scRNA/Results_Seurat_new_NoCorrection/ImmuneCells/DEG/Microglia - clusters/Microglia - cluster 1 and 2")
comp_deseq2 <- FindMarkers(immune.cells, ident.1 = "1", ident.2= "2", test.use="DESeq2")
comp_wilcox <- FindMarkers(immune.cells, ident.1 = "1", ident.2= "2")
write.csv(comp_deseq2, "deg_deseq2.csv")
write.csv(comp_wilcox, "deg_wilcox.csv")
geneList <- c(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
geneList_all <- intersect(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
write.csv(data.frame(geneList_all), "significant_gene_list_overlap.csv")

pdf("DEG_heatmap.pdf", paper='letter')
DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$res.2 %in% c("1", "2")],
          genes.use = geneList,
          group.label.rot=TRUE, group.cex = 10, cex.row = 6, slim.col.label = TRUE)
dev.off()

# boxplot
imm_data <- as.data.frame(as.matrix(immune.cells@data))
imm_meta <- as.data.frame(immune.cells@meta.data)
tmp_meta <- subset(imm_meta, imm_meta$res.2 == "1" | imm_meta$res.2 == "2")
tmp_meta <- tmp_meta[order(tmp_meta$batchid),]
tmp_data <- imm_data[,colnames(imm_data) %in% row.names(tmp_meta)]

for (i in 1:length(geneList)){
  tmp_plot <- data.frame(grp=tmp_meta$batchid, expr=t(tmp_data[row.names(tmp_data)==geneList[i],]))
  colnames(tmp_plot) <- c("grp","expr")
  png(filename=paste0("boxplot/",geneList[i],"_gene.png"))
  gp <- ggplot(tmp_plot, aes(grp, expr)) + 
    geom_boxplot(outlier.colour="Black", outlier.shape=1, aes(colour=grp), show.legend = FALSE) +
    geom_jitter(width=0.2, shape = 16, size=1, aes(colour = grp), show.legend = FALSE) +
    ggtitle(paste0(geneList[i]," gene")) +                              ##### name of cell type
    xlab("") + ylab("Expression level")
  print(gp)
  dev.off()
}



#####  DEG between Microglia 1 and 3
setwd("C:/Users/kicheol/Desktop/Projects/Andres_scRNA/Results_Seurat_new_NoCorrection/ImmuneCells/DEG/Microglia - clusters/Microglia - cluster 1 and 3")
comp_deseq2 <- FindMarkers(immune.cells, ident.1 = "1", ident.2= "3", test.use="DESeq2")
comp_wilcox <- FindMarkers(immune.cells, ident.1 = "1", ident.2= "3")
write.csv(comp_deseq2, "deg_deseq2.csv")
write.csv(comp_wilcox, "deg_wilcox.csv")
geneList <- c(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
geneList_all <- intersect(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
write.csv(data.frame(geneList_all), "significant_gene_list_overlap.csv")

pdf("DEG_heatmap.pdf", paper='letter')
DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$res.2 %in% c("1", "3")],
          genes.use = geneList,
          group.label.rot=TRUE, group.cex = 10, cex.row = 6, slim.col.label = TRUE)
dev.off()

# boxplot
imm_data <- as.data.frame(as.matrix(immune.cells@data))
imm_meta <- as.data.frame(immune.cells@meta.data)
tmp_meta <- subset(imm_meta, imm_meta$res.2 == "1" | imm_meta$res.2 == "3")
tmp_meta <- tmp_meta[order(tmp_meta$batchid),]
tmp_data <- imm_data[,colnames(imm_data) %in% row.names(tmp_meta)]

for (i in 1:length(geneList)){
  tmp_plot <- data.frame(grp=tmp_meta$batchid, expr=t(tmp_data[row.names(tmp_data)==geneList[i],]))
  colnames(tmp_plot) <- c("grp","expr")
  png(filename=paste0("boxplot/",geneList[i],"_gene.png"))
  gp <- ggplot(tmp_plot, aes(grp, expr)) + 
    geom_boxplot(outlier.colour="Black", outlier.shape=1, aes(colour=grp), show.legend = FALSE) +
    geom_jitter(width=0.2, shape = 16, size=1, aes(colour = grp), show.legend = FALSE) +
    ggtitle(paste0(geneList[i]," gene")) +                              ##### name of cell type
    xlab("") + ylab("Expression level")
  print(gp)
  dev.off()
}



#####  DEG between Microglia 1 and 5
setwd("C:/Users/kicheol/Desktop/Projects/Andres_scRNA/Results_Seurat_new_NoCorrection/ImmuneCells/DEG/Microglia - clusters/Microglia - cluster 1 and 5")
comp_deseq2 <- FindMarkers(immune.cells, ident.1 = "1", ident.2= "5", test.use="DESeq2")
comp_wilcox <- FindMarkers(immune.cells, ident.1 = "1", ident.2= "5")
write.csv(comp_deseq2, "deg_deseq2.csv")
write.csv(comp_wilcox, "deg_wilcox.csv")
geneList <- c(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
geneList_all <- intersect(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
write.csv(data.frame(geneList_all), "significant_gene_list_overlap.csv")

pdf("DEG_heatmap.pdf", paper='letter')
DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$res.2 %in% c("1", "5")],
          genes.use = geneList,
          group.label.rot=TRUE, group.cex = 10, cex.row = 6, slim.col.label = TRUE)
dev.off()

# boxplot
imm_data <- as.data.frame(as.matrix(immune.cells@data))
imm_meta <- as.data.frame(immune.cells@meta.data)
tmp_meta <- subset(imm_meta, imm_meta$res.2 == "1" | imm_meta$res.2 == "5")
tmp_meta <- tmp_meta[order(tmp_meta$batchid),]
tmp_data <- imm_data[,colnames(imm_data) %in% row.names(tmp_meta)]

for (i in 1:length(geneList)){
  tmp_plot <- data.frame(grp=tmp_meta$batchid, expr=t(tmp_data[row.names(tmp_data)==geneList[i],]))
  colnames(tmp_plot) <- c("grp","expr")
  png(filename=paste0("boxplot/",geneList[i],"_gene.png"))
  gp <- ggplot(tmp_plot, aes(grp, expr)) + 
    geom_boxplot(outlier.colour="Black", outlier.shape=1, aes(colour=grp), show.legend = FALSE) +
    geom_jitter(width=0.2, shape = 16, size=1, aes(colour = grp), show.legend = FALSE) +
    ggtitle(paste0(geneList[i]," gene")) +                              ##### name of cell type
    xlab("") + ylab("Expression level")
  print(gp)
  dev.off()
}




#####  DEG between Microglia 1 and 6
setwd("C:/Users/kicheol/Desktop/Projects/Andres_scRNA/Results_Seurat_new_NoCorrection/ImmuneCells/DEG/Microglia - clusters/Microglia - cluster 1 and 6")
comp_deseq2 <- FindMarkers(immune.cells, ident.1 = "1", ident.2= "6", test.use="DESeq2")
comp_wilcox <- FindMarkers(immune.cells, ident.1 = "1", ident.2= "6")
write.csv(comp_deseq2, "deg_deseq2.csv")
write.csv(comp_wilcox, "deg_wilcox.csv")
geneList <- c(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
geneList_all <- intersect(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
write.csv(data.frame(geneList_all), "significant_gene_list_overlap.csv")

pdf("DEG_heatmap.pdf", paper='letter')
DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$res.2 %in% c("1", "6")],
          genes.use = geneList,
          group.label.rot=TRUE, group.cex = 10, cex.row = 6, slim.col.label = TRUE)
dev.off()

# boxplot
imm_data <- as.data.frame(as.matrix(immune.cells@data))
imm_meta <- as.data.frame(immune.cells@meta.data)
tmp_meta <- subset(imm_meta, imm_meta$res.2 == "1" | imm_meta$res.2 == "6")
tmp_meta <- tmp_meta[order(tmp_meta$batchid),]
tmp_data <- imm_data[,colnames(imm_data) %in% row.names(tmp_meta)]

for (i in 1:length(geneList)){
  tmp_plot <- data.frame(grp=tmp_meta$batchid, expr=t(tmp_data[row.names(tmp_data)==geneList[i],]))
  colnames(tmp_plot) <- c("grp","expr")
  png(filename=paste0("boxplot/",geneList[i],"_gene.png"))
  gp <- ggplot(tmp_plot, aes(grp, expr)) + 
    geom_boxplot(outlier.colour="Black", outlier.shape=1, aes(colour=grp), show.legend = FALSE) +
    geom_jitter(width=0.2, shape = 16, size=1, aes(colour = grp), show.legend = FALSE) +
    ggtitle(paste0(geneList[i]," gene")) +                              ##### name of cell type
    xlab("") + ylab("Expression level")
  print(gp)
  dev.off()
}



#####  DEG between Microglia 1 and 8
setwd("C:/Users/kicheol/Desktop/Projects/Andres_scRNA/Results_Seurat_new_NoCorrection/ImmuneCells/DEG/Microglia - clusters/Microglia - cluster 1 and 8")
comp_deseq2 <- FindMarkers(immune.cells, ident.1 = "1", ident.2= "8", test.use="DESeq2")
comp_wilcox <- FindMarkers(immune.cells, ident.1 = "1", ident.2= "8")
write.csv(comp_deseq2, "deg_deseq2.csv")
write.csv(comp_wilcox, "deg_wilcox.csv")
geneList <- c(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
geneList_all <- intersect(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
write.csv(data.frame(geneList_all), "significant_gene_list_overlap.csv")

pdf("DEG_heatmap.pdf", paper='letter')
DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$res.2 %in% c("1", "8")],
          genes.use = geneList,
          group.label.rot=TRUE, group.cex = 10, cex.row = 6, slim.col.label = TRUE)
dev.off()

# boxplot
imm_data <- as.data.frame(as.matrix(immune.cells@data))
imm_meta <- as.data.frame(immune.cells@meta.data)
tmp_meta <- subset(imm_meta, imm_meta$res.2 == "1" | imm_meta$res.2 == "8")
tmp_meta <- tmp_meta[order(tmp_meta$batchid),]
tmp_data <- imm_data[,colnames(imm_data) %in% row.names(tmp_meta)]

for (i in 1:length(geneList)){
  tmp_plot <- data.frame(grp=tmp_meta$batchid, expr=t(tmp_data[row.names(tmp_data)==geneList[i],]))
  colnames(tmp_plot) <- c("grp","expr")
  png(filename=paste0("boxplot/",geneList[i],"_gene.png"))
  gp <- ggplot(tmp_plot, aes(grp, expr)) + 
    geom_boxplot(outlier.colour="Black", outlier.shape=1, aes(colour=grp), show.legend = FALSE) +
    geom_jitter(width=0.2, shape = 16, size=1, aes(colour = grp), show.legend = FALSE) +
    ggtitle(paste0(geneList[i]," gene")) +                              ##### name of cell type
    xlab("") + ylab("Expression level")
  print(gp)
  dev.off()
}



#####  DEG between Microglia 2 and 5
setwd("C:/Users/kicheol/Desktop/Projects/Andres_scRNA/Results_Seurat_new_NoCorrection/ImmuneCells/DEG/Microglia - clusters/Microglia - cluster 2 and 5")
comp_deseq2 <- FindMarkers(immune.cells, ident.1 = "2", ident.2= "5", test.use="DESeq2")
comp_wilcox <- FindMarkers(immune.cells, ident.1 = "2", ident.2= "5")
write.csv(comp_deseq2, "deg_deseq2.csv")
write.csv(comp_wilcox, "deg_wilcox.csv")
geneList <- c(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
geneList_all <- intersect(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
write.csv(data.frame(geneList_all), "significant_gene_list_overlap.csv")

pdf("DEG_heatmap.pdf", paper='letter')
DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$res.2 %in% c("2", "5")],
          genes.use = geneList,
          group.label.rot=TRUE, group.cex = 10, cex.row = 6, slim.col.label = TRUE)
dev.off()

# boxplot
imm_data <- as.data.frame(as.matrix(immune.cells@data))
imm_meta <- as.data.frame(immune.cells@meta.data)
tmp_meta <- subset(imm_meta, imm_meta$res.2 == "2" | imm_meta$res.2 == "5")
tmp_meta <- tmp_meta[order(tmp_meta$batchid),]
tmp_data <- imm_data[,colnames(imm_data) %in% row.names(tmp_meta)]

for (i in 1:length(geneList)){
  tmp_plot <- data.frame(grp=tmp_meta$batchid, expr=t(tmp_data[row.names(tmp_data)==geneList[i],]))
  colnames(tmp_plot) <- c("grp","expr")
  png(filename=paste0("boxplot/",geneList[i],"_gene.png"))
  gp <- ggplot(tmp_plot, aes(grp, expr)) + 
    geom_boxplot(outlier.colour="Black", outlier.shape=1, aes(colour=grp), show.legend = FALSE) +
    geom_jitter(width=0.2, shape = 16, size=1, aes(colour = grp), show.legend = FALSE) +
    ggtitle(paste0(geneList[i]," gene")) +                              ##### name of cell type
    xlab("") + ylab("Expression level")
  print(gp)
  dev.off()
}



#####  DEG between Microglia 3 and 6
setwd("C:/Users/kicheol/Desktop/Projects/Andres_scRNA/Results_Seurat_new_NoCorrection/ImmuneCells/DEG/Microglia - clusters/Microglia - cluster 3 and 6")
comp_deseq2 <- FindMarkers(immune.cells, ident.1 = "3", ident.2= "6", test.use="DESeq2")
comp_wilcox <- FindMarkers(immune.cells, ident.1 = "3", ident.2= "6")
write.csv(comp_deseq2, "deg_deseq2.csv")
write.csv(comp_wilcox, "deg_wilcox.csv")
geneList <- c(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
geneList_all <- intersect(row.names(subset(comp_deseq2, p_val_adj<0.05)),row.names(subset(comp_wilcox, p_val_adj<0.05)))
write.csv(data.frame(geneList_all), "significant_gene_list_overlap.csv")

pdf("DEG_heatmap.pdf", paper='letter')
DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$res.2 %in% c("3", "6")],
          genes.use = geneList,
          group.label.rot=TRUE, group.cex = 10, cex.row = 6, slim.col.label = TRUE)
dev.off()

# boxplot
imm_data <- as.data.frame(as.matrix(immune.cells@data))
imm_meta <- as.data.frame(immune.cells@meta.data)
tmp_meta <- subset(imm_meta, imm_meta$res.2 == "3" | imm_meta$res.2 == "6")
tmp_meta <- tmp_meta[order(tmp_meta$batchid),]
tmp_data <- imm_data[,colnames(imm_data) %in% row.names(tmp_meta)]

for (i in 1:length(geneList)){
  tmp_plot <- data.frame(grp=tmp_meta$batchid, expr=t(tmp_data[row.names(tmp_data)==geneList[i],]))
  colnames(tmp_plot) <- c("grp","expr")
  png(filename=paste0("boxplot/",geneList[i],"_gene.png"))
  gp <- ggplot(tmp_plot, aes(grp, expr)) + 
    geom_boxplot(outlier.colour="Black", outlier.shape=1, aes(colour=grp), show.legend = FALSE) +
    geom_jitter(width=0.2, shape = 16, size=1, aes(colour = grp), show.legend = FALSE) +
    ggtitle(paste0(geneList[i]," gene")) +                              ##### name of cell type
    xlab("") + ylab("Expression level")
  print(gp)
  dev.off()
}


################################################################################################
# time-point comparison for all cells
setwd("C:/Users/kicheol/Desktop/Projects/Andres_scRNA/Results_Seurat_new_NoCorrection/ImmuneCells/DEG/Time - All immune cells")
immune.cells <- SetAllIdent(immune.cells, id = "batchid")

comp_wilcox1 <- FindMarkers(immune.cells, ident.1 = "1_Baseline", ident.2= "2_PreOnset", logfc.threshold = 0.1)
comp_wilcox2 <- FindMarkers(immune.cells, ident.1 = "1_Baseline", ident.2= "3_Chronic", logfc.threshold = 0.1)
comp_wilcox3 <- FindMarkers(immune.cells, ident.1 = "2_PreOnset", ident.2= "3_Chronic", logfc.threshold = 0.1)
write.csv(comp_wilcox1, "comp_ALL-1Base_2Pre-wilcox.csv")
write.csv(comp_wilcox2, "comp_ALL-1Base_3Chro-wilcox.csv")
write.csv(comp_wilcox3, "comp_ALL-2Pre_3Chro-wilcox.csv")


comp1 <- row.names(subset(comp_wilcox1, p_val_adj<0.05))
comp2 <- row.names(subset(comp_wilcox2, p_val_adj<0.05))
comp3 <- row.names(subset(comp_wilcox3, p_val_adj<0.05))
comp <- c(comp1, comp2, comp3)


# heatmap
pdf("DEG_heatmap-ImmuneCells.pdf", paper='letter')
DoHeatmap(immune.cells, use.scaled = TRUE,
          cells.use=rownames(immune.cells@meta.data)[immune.cells@meta.data$batchid %in% c("1_Baseline", "2_PreOnset", "3_Chronic")],
          genes.use = comp,
          group.label.rot=FALSE, group.cex = 10, cex.row = 2.5, slim.col.label = TRUE)
dev.off()
