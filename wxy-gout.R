
#1.数据清洗

library(devtools)
library(GEOmirror)
library(idmap3)
library(GEOquery)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(dplyr)
library(limma)
library(tidyr)

# 下载GEO数据和注释信息 ---------------------------------------------------------
data_use <- read.csv("mRNA.csv")
rownames(data_use) <- data_use[,1]
data_use <- data_use[,-1]
data_use <- as.data.frame(data_use)

group <- read.csv("分组.csv",header = T)
group <- group[,-1]

# 判断表达矩阵是否需要对数处理 --------------------------------------------------
boxplot(data_use,las=2)#las完整展示行名
qx <- as.numeric(quantile(data_use, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))  
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) 
LogC
if (LogC) { 
  data_use[which(data_use <= 0)] <- NaN
  data_use <- log2(data_use) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
}   #进行log2取对数

# 数据检验
boxplot(data_use,las=2)#las完整展示行名

PCA_result <- PCA(t(data_use), graph = F)
fviz_pca_ind(PCA_result,
             geom.ind = c("point","text"), 
             mean.point = F, 
             repel = T, 
             col.ind = group$group, 
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             legend.title = "Groups")

Sample_clust <- dist(t(data_use)) #计算变量间距离
hc <- hclust(Sample_clust) #hclust进行聚类
plot(hc,
     hang = -1,
     cex = 0.8)

# 获取注释信息 ------------------------------------------------------------------
id3 <- idmap3::get_pipe_IDs("GPL21827")
colnames(id3)
# 探针ID转换为gene_symbol -------------------------------------------------------
data_use$probe_id <- rownames(data_use)
data_with_name <- merge(data_use,
                        id3,
                        by.x = "probe_id",
                        by.y = "probe_id") 

dim(data_with_name)
data_with_name <- data_with_name[data_with_name$symbol != "",] 
dim(data_with_name)
table(duplicated(data_with_name[,ncol(data_with_name)]))
data_with_name <- avereps(data_with_name[,-c(1,ncol(data_with_name))],
                          ID = data_with_name$symbol) 

table(duplicated(rownames(data_with_name))) 
data_with_name['GAPDH',] 
data_with_name['ACTB',]
write.table(data.frame(gene_symbol=rownames(data_with_name),data_with_name),
            file = "Matrix_of_expression.txt",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
write.csv(data_with_name,"Matrix_of_expression.csv")



#2. 差异表达

# 加载所需要的R包 ----------------------------------------------------------------
library(limma)
library(ggplot2) 
library(pheatmap) 

# 输入表达矩阵和分组文件 -------------------------------------------------------------
expr_data<-read.table("Matrix_of_expression.txt",header = T,
                      row.names = 1,sep = "\t")
group<-read.csv("分组.csv",header = T,row.names = 1,sep = ",")

# #构建分组矩阵--design ---------------------------------------------------------
design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr_data)

# #构建比较矩阵——contrast -------------------------------------------------------
contrast.matrix <- makeContrasts(Gout-Control,levels = design)

# #线性拟合模型构建 ---------------------------------------------------------------
fit <- lmFit(expr_data,design) #非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)  
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = 1,n = Inf)
DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
                       ifelse(DEG$logFC > 2, "up-regulated",
                              ifelse(DEG$logFC < -2, "down-regulated", "unchanged")))
table(DEG$regulate)
write.table(data.frame(gene_symbol=rownames(DEG),DEG),file = "DEG_result.txt",
            sep = "\t",quote = F,row.names = F,col.names = T)

# 区分上下调基因 -----------------------------------------------------------------
DE_1_0.05 <- DEG[DEG$P.Value<0.05&abs(DEG$logFC)>2,]
upGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "up-regulated",]
downGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "down-regulated",]
write.csv(upGene_1_0.05,"upGene_1_0.05.csv")
write.csv(downGene_1_0.05,"down-regulated.csv")



# 火山图的绘制 ------------------------------------------------------------------
pdf("volcano.pdf")
ggplot(DEG,aes(x=logFC,y=-log10(P.Value)))+ 
  geom_point(alpha=0.6,size=3.5,aes(color=regulate))+ 
  ylab("-log10(P.Value)")+ 
  scale_color_manual(values = c("blue", "grey", "red"))+ 
  geom_vline(xintercept = c(-1,1),lty=4,col ="black",lwd=0.8)+ 
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ 
  theme_bw()  
dev.off()





#3. WGCNA

library("WGCNA")

clinical <- read.csv("clinical.csv",header=T)
rownames(clinical) <- clinical[,1]
clinical <- clinical[,-1]
class(clinical)

fpkm <- read.table("Matrix_of_expression.txt",header=T,
                   sep="\t",dec=".")
rownames(fpkm) <- fpkm[,1]
fpkm <- fpkm[,-1]
fpkm <- as.matrix(fpkm)

WGCNA_matrix = t(fpkm[order(apply(fpkm,1,mad), decreasing = T)[1:10000],])

datExpr <- WGCNA_matrix  
sampleNames = rownames(datExpr)
##对样本聚类
sampleTree = hclust(dist(datExpr), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", 
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

############################################################
##########################软阈值ֵ############################
############################################################
powers = c(seq(from = 1, to=10, by=2),c(11:20))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
par(mfrow = c(1,2));
cex1 = 0.90;###可选0.85/0.90
# Scale-free topology fit index as a function of the soft-thresholding power
png(file = "SoftThreshold.png",width = 800, height = 600)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.80,col="red") #0.85/0.90太少可选择0.8，不能低于0.8

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()
power = sft$powerEstimate
power
####################################################################
###########################一步式网络构建###########################
####################################################################
net = blockwiseModules(
  datExpr,
  power = 15,      
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 100,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "AS-green-FPKM-TOM",
  verbose = 3
)
table(net$colors)


png(file = "moduleCluster.png", width = 1200, height = 800)

mergedColors = labels2colors(net$colors)

plotDendroAndColors(dendro = net$dendrograms[[1]], 
                    colors = mergedColors[net$blockGenes[[1]]],
                    groupLabels = "Module colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

MEs = moduleEigengenes(datExpr, mergedColors)$eigengenes

MEs = orderMEs(MEs)
####################################################################
############################与表型数据相结合########################
####################################################################

datTraits <- clinical
moduleTraitCor=cor(MEs, datTraits, use="p")
#write.table(file="Step04-modPhysiological.cor.xls",moduleTraitCor,sep="\t",quote=F)
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, nSamples)
#write.table(file="Step04-modPhysiological.p.xls",moduleTraitPvalue,sep="\t",quote=F)



png(file="Module_trait_relationships.png",width=800,height=900)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)

labeledHeatmap(Matrix=moduleTraitCor,
               xLabels=colnames(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=TRUE,
               cex.text=1.2,
               cex.lab=0.9,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))


dev.off()



########################################################################
MF_stage = as.data.frame(datTraits$Gout)

names(MF_stage) = "MF_failure"

modNames = substring(names(MEs), 3) 


geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));


names(geneModuleMembership) = paste("MF", modNames, sep="");
names(MMPvalue) = paste("p.MF", modNames, sep="");


geneTraitSignificance = as.data.frame(cor(datExpr, MF_stage, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));


names(geneTraitSignificance) = paste("GS.", names(MF_stage), sep="")
names(GSPvalue) = paste("p.GS.", names(MF_stage), sep="")


module = "turquoise" #选定颜色

column = match(module, modNames);

moduleColors=mergedColors

moduleGenes = moduleColors==module;



png(file="Module_membership_vs_gene_significance.png")
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

turquoise_gene <- colnames(datExpr)[moduleGenes]

turquoise <- as.data.frame(turquoise_gene)
write.csv(turquoise,file = "turquoise_gene.csv")




#4.KEGG富集分析
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("pathview")
library("ggnewscale")
library("DOSE")
library(stringr)

pvalueFilter=0.05        
qvalueFilter=1        
showNum=20


rt=read.table("target.txt",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)

kk@result$Description= gsub(' - Homo sapiens \\(human\\)','',kk@result$Description)

KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

#KEGG所有结果
write.table(KEGG,file="KEGG.xls",sep="\t",quote=F,row.names = F)

#筛选前20
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

pdf(file="KEGG_barplot.pdf",width = 9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel) +scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

pdf(file="KEGG_bubble.pdf",width = 9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

pdf(file="KEGG_cnet.pdf",width = 9,height = 8)
af=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(af, showCategory = showNum, categorySize="pvalue",circular = TRUE,colorEdge = TRUE,cex_label_category=0.65,cex_label_gene=0.6)
dev.off()

pdf(file="KEGG_net.pdf",width = 9,height = 7)
x2 <- pairwise_termsim(kk)
emapplot(x2,showCategory = showNum,cex_label_category=0.65,color = "pvalue",layout ="nicely")
dev.off()  


#5.GO富集分析
library("org.Hs.eg.db")  
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("ggnewscale")
library("enrichplot")
library("DOSE")
library(stringr)

pvalueFilter=0.05         
qvalueFilter=1  
showNum=8

rt=read.table("target.txt",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}


kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

write.table(GO,file="GO.xls",sep="\t",quote=F,row.names = F)


if(nrow(GO)<30){
  showNum=nrow(GO)
}

pdf(file="GO_barplot.pdf",width = 9,height =7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
print(bar)
dev.off()


pdf(file="GO_bubble.pdf",width = 9,height =7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
print(bub)
dev.off()

pdf(file="GO_cnet.pdf",width = 10,height = 5)
af=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kk, showCategory = 10, categorySize="pvalue",circular = TRUE,colorEdge = TRUE,cex_label_category=0.65,cex_label_gene=0.6)
dev.off()

pdf(file="GO_net.pdf",width = 9,height = 7)
x2 <- pairwise_termsim(kk)
emapplot(x2,showCategory = 20,cex_label_category=0.65,color = "pvalue",layout ="nicely" )
dev.off()



#6.SSGSEA

library(GSVA)#算法
library(tidyverse)#数据处理
library(ggpubr)#绘图
library(ggplot2)#绘图
library(pheatmap)#绘制热图

library(tinyarray)
library(dplyr)
library(Hmisc)


DEG_expr <- read.table("Matrix_of_expression.txt",header=T,
                       sep="\t",dec=".",row.names = 1)
group <- read.csv("分组.csv")
group <- group[,-1]
colnames(group)[1] <- "Samples"
colnames(group)[2] <- "group"
markergenes <- read.csv("markergenes.csv")
table(markergenes$Cell.type)
#数据检验
boxplot(DEG_expr,outline=F, notch=F , las=2)
qx <- as.numeric(quantile(DEG_expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))  #数据的分布，样本分位数
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)#判断是否进行log的标准
LogC
if (LogC) { 
  DEG_expr[which(DEG_expr <= 0)] <- NaN
  DEG_expr <- log2(DEG_expr) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
}  



geneset <- split(markergenes,markergenes$Cell.type)

im_geneset <- lapply(geneset, function(x){
  gene = x$Metagene
  unique(gene)
})
lapply(im_geneset[1:3], head)

save(im_geneset,file = "im_geneset.Rdata")

DEG_expr <- as.matrix(DEG_expr) 



result <- gsva(DEG_expr,im_geneset,method = "ssgsea")
result1 <- as.data.frame(t(result))

write.csv(result1,"ssGSEA_result.csv")


data <- cbind(group,result1)
colnames(data)
data <- pivot_longer(data = data,
                     cols = 3:30,
                     names_to = "celltype",
                     values_to = "proportion")
#开始绘图
pdf(file="分组箱式图+显著性P值.pdf",width = 10,height = 8)
ggboxplot(data = data,
          x = "celltype",
          y = "proportion",
          combine = TRUE,
          merge = FALSE,
          color = "black",
          fill = "group",
          palette = c("#1C3EDF","#DF1C26"),
          title = NULL,
          xlab = "ssGSEA",
          ylab = "Expression",
          bxp.errorbar = FALSE,
          bxp.errorbar.width = 0.2,
          facet.by = NULL,
          panel.labs = NULL,
          short.panel.labs = TRUE,
          linetype = "solid",
          size = NULL,
          width = 0.8,
          notch = FALSE,
          outlier.shape = 20,
          select = NULL,
          remove = NULL,
          order = NULL,
          error.plot = "pointrange",
          label = NULL,
          font.label = list(size = 12, color = "black"),
          label.select = NULL,
          repel = TRUE,
          label.rectangle = TRUE, ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1)) +  
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.",hide.ns = F,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("***", "**", "*", "ns")))
dev.off()  


# 基因相关性分析 ---------------------------------------------------------------
choose_gene_1se <- c("CXCL8","CXCL1","CXCL2","IL1B","IL6")
nc = t(rbind(result,DEG_expr[choose_gene_1se,]))
m = rcorr(nc)$r[1:nrow(result),(ncol(nc)-length(choose_gene_1se)+1):ncol(nc)]

p = rcorr(nc)$P[1:nrow(result),(ncol(nc)-length(choose_gene_1se)+1):ncol(nc)]
p
p=as.data.frame(p)

tmp <- mutate_all(p, funs(case_when(. < 0.0001 ~ "****",
                                    . < 0.001 ~ "***",
                                    . < 0.01 ~ "**",
                                    . < 0.05 ~ "*",
                                    TRUE ~ "")))


p1 <- pheatmap(t(m),
               display_numbers =t(tmp),
               angle_col =90,
               color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
               border_color = "white",
               cellwidth = 15, 
               cellheight = 15,
               width = 7, 
               height=9.1,
               treeheight_col = 0,
               treeheight_row = 0)



#7. CERNA
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="http://mirrors.cloud.tencent.com/CRAN/"))
options(download.file.method = 'libcurl')
options(url.method='libcurl')


library(multiMiR)

### 1.从mRNA得到miRNA,9hubgenes.txt。
x = read.table("hubgenes.txt",stringsAsFactors = F)$V1;x


gene2mir <- get_multimir(org     = 'hsa',
                         target  = x,
                         table   = 'validated',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff      = 500000)
mit = gene2mir@data[gene2mir@data$database=="mirtarbase",];dim(mit)


table(mit$support_type)

miRNAs = unique(mit$mature_mirna_id)

### 2.从miRNA得到lncRNA


starbase = data.table::fread("starBaseV3_hg19_CLIP-seq_lncRNA_all.txt");dim(starbase)


load("anno.Rdata")
lnc_anno$gene_id = stringr::str_remove(lnc_anno$gene_id,"\\.\\d")
p1 = starbase$geneName %in% lnc_anno$gene_name;table(p1)
p2 = starbase$geneID %in% lnc_anno$gene_id;table(p2)

starbase = starbase[p2,];dim(starbase)
lnc_mi = starbase[starbase$miRNAname %in% miRNAs,]
colnames(lnc_mi)


p2 = lnc_mi$pancancerNum >10  &lnc_mi$clipExpNum>5;table(p2)
p3 = lnc_mi$pancancerNum >10  &lnc_mi$clipExpNum>5 & lnc_mi$degraExpNum >0;table(p3)
lnc_mi$geneName[p2]
lnc_mi$geneName[p3]

### 3.准备cytoscape的输入文件
ez1 = mit[,3:4]
ez2 = lnc_mi[p2,c(2,4)]
colnames(ez1) = colnames(ez2)

library(dplyr)
ez = rbind(ez1,ez2);dim(ez)
ez = distinct(ez,miRNAname,geneName);dim(ez)


tp = data.frame(nodes = c(ez1$miRNAname,
                          ez2$miRNAname,
                          ez1$geneName,
                          ez2$geneName),
                type = rep(c("mi","pc","lnc"),
                           times = c(nrow(ez1)+nrow(ez2),
                                     nrow(ez1),
                                     nrow(ez2))
                )
)
dim(tp);head(tp)

tp = distinct(tp,nodes,.keep_all = T)
table(tp$type)

write.csv(ez,file = "ez.csv")
write.csv(tp,file = "tp.csv")











