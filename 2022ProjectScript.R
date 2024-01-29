library(TCGAbiolinks)
install.packages("digest")
install.packages("rlang")
install.packages("RSQLite")
install.packages("yaml")
library(TCGAbiolinks)

projectString <- "TCGA-STAD"
query <- GDCquery(project = projectString,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts",
                  experimental.strategy = "RNA-Seq")

table(query$results[[1]]$sample_type)

library("SummarizedExperiment")

dim(query$results[[1]]) #Check results

query2 <- query
dupeIDsLogic <- duplicated(query2$results[[1]]$cases.submitter_id)
dupeIDs <- query2$results[[1]]$cases.submitter_id[dupeIDsLogic]
normalLogic <- query2$results[[1]]$sample_type=="Solid Tissue Normal"|
  query2$results[[1]]$sample_type=="Recurrent Tumor"
normalIDs          <-query2$results[[1]]$cases.submitter_id[normalLogic]
dupeIDs.inNormal   <-dupeIDs[dupeIDs%in%normalIDs]
dupeIDs.RowsLogic  <-query2$results[[1]]$cases.submitter_id%in%dupeIDs.inNormal
query2$results[[1]]<-query2$results[[1]][dupeIDs.RowsLogic,]


GDCdownload(query2,files.per.chunk = 50)
cnt_og <- GDCprepare(query2)

save(cnt_og,file = paste0("cnt_og.",projectString,".Rdata"))
save.image(paste0("Image_01AfterPreparingCntFile_",projectString,".Rimage"))

#-------------------------------------------------------------------------------
cnt_woNA <- cnt_og[,!is.na(cnt_og$age_at_index)]
table (cnt_woNA$age_at_index,useNA="ifany")
View (as.data.frame(colData(cnt_woNA)))

cntPT<-cnt_woNA[,cnt_woNA$shortLetterCode=="TP"]
View(as.data.frame(colData(cntPT)))

cnt_sub_below68<-cntPT[,cntPT$age_at_index<68]
cnt_sub_68nabove<-cntPT[,cntPT$age_at_index>=68]

cnt_zero_below68<-assays(cnt_sub_below68)[[1]]==0
cnt_zero_68nabove<-assays(cnt_sub_68nabove)[[1]]==0
meanCnt0_below68<-rowMeans(cnt_zero_below68,na.rm=T)
meanCnt0_68nabove<-rowMeans(cnt_zero_68nabove,na.rm=1)
zeroFreqLogic<-meanCnt0_below68<0.50 & meanCnt0_68nabove<0.50
cntPT<-cntPT[zeroFreqLogic,]

save.image(paste0("Image_02PriorToDESeq2_",projectString,".Rimage"))
save(cntPT, file = paste0("cnt_sub_",projectString,".rData"))

cnt_de <- cntPT
cnt_de$comp <- cnt_de$age_at_index<68
cnt_de$comp <- as.factor(cnt_de$comp)
levels(cnt_de$comp) 
cnt_de$comp <- relevel(cnt_de$comp, "TRUE")
levels(cnt_de$comp)
levels(cnt_de$comp)<-gsub(" ","_",levels(cnt_de$comp))

ddsSE <- DESeqDataSet(cnt_de, design = ~ comp) 
dds   <- DESeq(ddsSE) 
res   <- results(dds)

resOutput <- cbind(as.data.frame(rowRanges(cnt_de)),
                   as.data.frame(res))
fdr.cut.off <- 0.01 
lfc.cut.off <- 1.50
resOutput$padjLteCutoff <- resOutput$padj <= fdr.cut.off
resOutput$absFoldGteCutoff <- abs(resOutput$log2FoldChange) >= lfc.cut.off
resOutput$bothCutoffs <- resOutput$absFoldGteCutoff & resOutput$padjLteCutoff

table(paste0("pAdj_log2Change:",resOutput$padjLteCutoff,"_",resOutput$absFoldGteCutoff))
table(resOutput$bothCutoffs)
table(resOutput$bothCutoffs)/sum(table(resOutput$bothCutoffs))

resPlotted <- resOutput[sample(1:nrow(resOutput),size = nrow(resOutput)),] 
volc1 <- ggplot(resPlotted, aes(x     = log2FoldChange,
                                y     = -log10(padj),
                                color = baseMean ,
                                shape = bothCutoffs)) +
  geom_point() +
  geom_vline(aes(xintercept = -lfc.cut.off)       , color = "blue" , linetype = "dashed") +
  geom_vline(aes(xintercept =  lfc.cut.off)       , color = "red"  , linetype = "dashed") +
  geom_hline(aes(yintercept = -log10(fdr.cut.off)), color = "black", linetype = "dashed") +
  labs(title = "TCGA-STAD, Differential Gene Expression\nBetween PT Patients <68 Y/O and >=68 Y/O",
       y     = "-log10(FDR-Adjusted p-values)",
       x     = "log2(Fold Change)",
       color = "Base mean count",
       shape = "Sig. DE") +
  theme_bw() +
  geom_text_repel(data = head(resPlotted[order(resPlotted$padj),],10),
                  aes(label = gene_name),
                  show.legend = FALSE,color = "black",alpha = 0.8) +
  scale_color_viridis_c(trans ="log10")
volc1

geneInv.index    <- which.min(resPlotted$padj)
geneInv.name     <- resPlotted$gene_name[geneInv.index]
geneInv.ensemble <- rownames(resPlotted)[geneInv.index]
normCounts <- as.data.frame(counts(dds, normalized = T))
geneInv<- data.frame(
  sample     = rownames(colData(dds)),
  group      = cnt_de$comp,
  normCounts = unlist(normCounts[match(geneInv.ensemble,rownames(normCounts)),])
)
geneInv$rank <- rank(geneInv$normCounts,ties.method = "random")
geneInv$rank <- (geneInv$rank-min(geneInv$rank)) / (max(geneInv$rank) - min(geneInv$rank))
geneInv$rank <- geneInv$rank*0.8 - 0.4
geneInv$rankWithFactor <- geneInv$rank + as.numeric(geneInv$group)

oneGenePlot <- ggplot(geneInv, aes(group, normCounts))+
  geom_boxplot(outlier.shape = NA, color = "blue", fill="yellow", alpha=0.2 )+
  geom_point(mapping = aes(x = rankWithFactor),
             size = 3,shape = 8)+
  labs(title=paste0(geneInv.name))+
  scale_color_viridis_c()+theme_bw()

pca.sum <-summary(pca)
pctPC1<-paste0("PC1 (",round(pca.sum$importance[2,1]*100,1),"% variance explained)")
pctPC2<-paste0("PC2 (",round(pca.sum$importance[2,2]*100,1),"% variance explained)")

pca <- prcomp(t(normCounts),scale. = T)
plot(pca$x[,1],pca$x[,2], xlab = pctPC1, ylab = pctPC2,asp = 1,col = dds$comp)



plottedDf<-data.frame(pca$x)
plottedDf <- cbind(plottedDf,colData(dds))
scatterD3::scatterD3(PC1,PC2,data = plottedDf,
                     lab = sample_submitter_id,
                     labels_size = 0,
                     col_var = comp,
                     fixed = T)

ggsave(filename = "Figure_Volcano1.png",volc1,
       dpi = 600, width = 84*2,height = 84,units = "mm")
ggsave(filename = paste0("Figure_",geneInv.name,".png"),oneGenePlot,
       dpi = 600, width = 84*2,height = 84,units = "mm")
write.csv(table(resOutput$padjLteCutoff,resOutput$absFoldGteCutoff),
          file = "Table_DEAnalysisSummary.csv")
write.csv(resOutput[order(resOutput$padj,resOutput$gene_name),],
          file = "SupplementalTable_DEAnalysisResults.csv")
write.csv(colData(cnt_de),file = "SupplementalTable_SampleMetadata.csv")

save.image(paste0("Image_03AfterDESeq2_",project,".Rimage"))







