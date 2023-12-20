rm(list = ls())  
options(stringsAsFactors = F)
#加载相关的R包
library(FactoMineR)
library(factoextra)  
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
library(limma)

setwd('C:\\Users\\user\\Desktop\\乳腺癌\\Deseq2')
#读入count文件
gene <- read.table('../data_mrna_seq_v2_rsem.txt',sep = '\t',header = T,check.names = F)
head(gene$Hugo_Symbol)
gene<-gene[!gene$Hugo_Symbol%in%"",]
rowna<-gene$Hugo_Symbol
gene<-gene[,c(3:ncol(gene))]
gene<-as.matrix(gene)

rownames(gene)<-rowna
#去除重复的行名
gene<-avereps(gene)
gene<-as.data.frame(gene)
a<-colnames(gene)

#剔除全为0的

name1<-rownames(gene)
gene=as.data.frame(lapply(gene,as.numeric))
rownames(gene)<-name1
gene<-round(gene,digits = 0)
all <- apply(gene,1,function(x)all(x==0))
newdata <- gene[!all,]

colnames(newdata)<-a

group=sapply(strsplit(colnames(newdata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
table(group)

#分组
clinical<-read.table('../data_clinical_patient.txt',sep = '\t',header = T)

table(clinical$SUBTYPE)
table(clinical$CANCER_TYPE_ACRONYM)

co_name<-intersect(substr(colnames(newdata),1,12),clinical$PATIENT_ID)
table(substr(colnames(newdata),1,12)%in%co_name)

clinical<-clinical[match(substr(colnames(newdata),1,12),clinical$PATIENT_ID),]


doms<-c("PATIENT_ID","SUBTYPE","AGE","SEX",
      "PATH_M_STAGE","PATH_N_STAGE","PATH_T_STAGE","OS_STATUS","OS_MONTHS")

clinical<-clinical[,doms]
clinical$group<-ifelse(clinical$SUBTYPE=='BRCA_Her2','BRCA_Her2','BRCA_Other')
table(clinical$group)
identical(substr(colnames(newdata),1,12),clinical$PATIENT_ID)
write.csv(clinical,'clinical.csv')


coldata <- data.frame(group =factor( clinical$group))
coldata <- data.frame(rom.names=colnames(newdata),coldata)
coldata
str(newdata)


#构建Deseq2对象
dds<-DESeqDataSetFromMatrix(countData = round(newdata),colData=coldata,design=~group)
rld <- vst(dds, blind=FALSE)
write.table(assay(rld),  file="Deseq2_rld.txt", sep="\t", quote=F, col.names=NA)
dat <- as.data.frame(assay(rld))#获得Deseq2标准化的数据     


#去除在75%样本中不表达的基因
keep <- rowSums(dat>0) >= floor(0.75*ncol(dat))
table(keep)
dat<-dat[keep,]

#画箱图
color<-brewer.pal(7,'Paired')
p<-boxplot(dat, ylab="dat", col=color,main=" normalized data ",
           outline = F, notch = F)
dev.off()

#画树图
sampleDists <- dist(t(dat))   #dist默认计算矩阵行与行的距离， 因此需要转置
sampleDistMatrix <- as.matrix(sampleDists)  
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)  #选取热图的颜色
p0 <- pheatmap::pheatmap(sampleDistMatrix,
                         fontsize=7,
                         clustering_distance_rows=sampleDists,
                         clustering_distance_cols=sampleDists,
                         angle_col=45,
                         col=colors)
ggsave(p0,filename = 'check_dist.pdf',width = 7.5,height =6)
dev.off()
pdf("check_hclust.pdf")
plot(hclust(sampleDists))
dev.off()


#若用 rld 数据，还可使用DESeq2自带函数 
pca <- plotPCA(rld, ntop = 500, intgroup=c("group"))
ggsave(pca, filename = 'check_PCA.pdf',width = 7.5,height =6)
dev.off()




#预处理，过滤低表达基因
dds <- dds[rowSums(counts(dds) >= 1) >= 100, ] #至少有2个样品count>=1



#利用PCA主成分分析法评估数据质量
rld <- vst(dds, blind=FALSE) 
plotPCA(rld, intgroup="group",ntop=500)


#####################未进行sva预测的差异#########################
#差异基因分析
#1.标准方法
dds2<-DESeq(dds)
suppressMessages(dds2)
res<- results(dds2,contrast =c('group','BRCA_Her2','BRCA_Other'),pAdjustMethod='fdr',alpha=0.05)

#查看差异基因数目
res2 <- as.data.frame(res)
table(res2$padj<0.05)
#按P值排序
deseq_res2 <- as.data.frame(res2[order(res2$padj),])
deseq_res2$genename <-rownames(deseq_res2)
sig_deseq_res2 <-subset(deseq_res2,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < (-1)))
#先输出第7列，再输出前6列
write.csv(deseq_res2[c(7,1:6)],'DESeq2_res2.csv')
write.csv(sig_deseq_res2[c(7,1:6)],'sig DESeq22.csv')

diff_up2 = subset(deseq_res2,padj<0.05 & (log2FoldChange >1))
write.csv(diff_up2,file="diff_up2.csv",row.names = F)
diff_down2 = subset(deseq_res2,padj<0.05 & (log2FoldChange < -1))
write.csv(diff_down2,file="diff_down2.csv",row.names = F)



