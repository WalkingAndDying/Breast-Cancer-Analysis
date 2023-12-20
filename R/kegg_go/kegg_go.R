
rm(list = ls())
#相关R包载入：
library(dplyr)#数据清洗
library(org.Hs.eg.db)#ID转换
library(clusterProfiler)#富集分析
library(ggplot2)#绘图
library(RColorBrewer)#配色调整

setwd('C:\\Users\\user\\Desktop\\乳腺癌\\kegg_go')
#加载差异基因文件
DESeq2<-read.csv('../Deseq2//sig DESeq22.csv',header = T,row.names = 1)
head(DESeq2)

#添加上下调基因分组标签：
DESeq2$group <- case_when(
  DESeq2$log2FoldChange >=1 & DESeq2$padj <=0.05 ~ "up",
  DESeq2$log2FoldChange <=-1 & DESeq2$padj <= 0.05 ~ "down",
  abs(DESeq2$log2FoldChange) <= 2 ~ "none",
  DESeq2$padj >=0.05 ~ "none")
head(DESeq2)

table(DESeq2$group)

#分别筛选上调基因、下调基因或所有差异基因（上调+下调）：
up <- rownames(DESeq2)[DESeq2$group=="up"]#差异上调
down <- rownames(DESeq2)[DESeq2$group=="down"]#差异下调
diff <- c(up,down)#所有差异基因
head(up)
head(down)


#ID转换：
#查看可转换的ID类型：
columns(org.Hs.eg.db)


#########################先做上调基因的富集分析#############################################
#使用函数bitr(基于org.Hs.eg.db包)：
#up：
up_entrez <- bitr(up,
                  fromType = "SYMBOL",#现有的ID类型
                  toType = "ENTREZID",#需转换的ID类型
                  OrgDb = "org.Hs.eg.db")
#提示1.36%的基因无法映射是正常的，一般不同数据库ID转换会存在缺失
head(up_entrez)


#KEGG富集分析，我们以总差异基因(diff_entrez)为例：
KEGG_diff <- enrichKEGG(gene = up_entrez$ENTREZID,
                        organism = "hsa",#物种，Homo sapiens (human)
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

#将ENTREZID转化为可读的gene symbol
KEGG_diff<- setReadable(KEGG_diff, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

#保存富集结果
write.csv(KEGG_diff@result,file = 'KEGG_up.CSV')



library(enrichplot)
#富集条形图绘制：
#绘图


#这里快速绘图，会发现直接使用p值或q值绘图颜色不太能看出差别，这种时候可以考虑使用-log10p或者-log10q来绘图

#使用ggplot2进行更加自由的可视化呈现：
#先提取富集结果表前Top20：

KEGG_result<-KEGG_diff@result
KEGG_top20 <- KEGG_result[1:20,]

#指定绘图顺序（转换为因子）：
KEGG_top20$pathway <- factor(KEGG_top20$Description,levels = rev(KEGG_top20$Description))

#Top20富集数目条形图：
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))


#富集因子（Rich Factor）计算：
colnames(KEGG_top20) #方便复制所需列名

#Rich Factor = GenRatio/BgRatio
top20_rf <- apply(KEGG_top20,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  KEGG_top20_rf <- round(GeneRatio/BgRatio,2)
  KEGG_top20_rf
})
head(top20_rf)
KEGG_top20$Rich_Factor <- top20_rf


#富集因子版显著性气泡图绘制：
p3 <- ggplot(data = KEGG_top20,
             aes(x = Rich_Factor, # X轴用富集因子来映射
                 y = pathway))+
  geom_point(aes(size = Count,
                 color = -log10(pvalue)))+
  theme_bw()+
  scale_color_distiller(palette = "Spectral",direction = -1) +
  labs(x = "Rich Factor",
       y = "",
       title = "Dotplot of Enriched KEGG Pathways",
       size = "Gene Number") +
  mytheme
p3

pdf(file = 'up_kegg.pdf',width = 12,height = 15)
p3
dev.off()
############################################################


#GO(Gene Ontology)富集分析：
##MF(我们以总差异基因的GO富集为例):
GO_MF_diff <- enrichGO(gene = up_entrez$ENTREZID, #用来富集的差异基因
                       OrgDb = org.Hs.eg.db, #指定包含该物种注释信息的org包
                       ont = "MF", #可以三选一分别富集,或者"ALL"合并
                       pAdjustMethod = "BH", #多重假设检验矫正方法
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE) #是否将gene ID映射到gene name
#提取结果表格：
GO_MF_result <- GO_MF_diff@result



#CC:
GO_CC_diff <- enrichGO(gene = up_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
#提取结果表格：
GO_CC_result <- GO_CC_diff@result

#BP:
GO_BP_diff <- enrichGO(gene = up_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
#提取结果表格：
GO_BP_result <- GO_BP_diff@result



#MF、CC、BP三合一：
GO_all_diff <- enrichGO(gene = up_entrez$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = "ALL", #三合一选择“ALL”
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
#提取结果表格：
GO_all_result <- GO_all_diff@result
write.csv(GO_all_result,'GO_up_result.csv')

# 绘制GO富集分析条形图，结果默认按qvalue升序，分别选出前十的term进行绘图即可
goBP <- subset(GO_all_result,subset = (ONTOLOGY == "BP"))[1:10,]
goCC <- subset(GO_all_result,subset = (ONTOLOGY == "CC"))[1:7,]
goMF <- subset(GO_all_result,subset = (ONTOLOGY == "MF"))[1:10,]
go.df <- rbind(goBP,goCC,goMF)
# 使画出的GO term的顺序与输入一致
library(tidyverse)

go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
# 绘图
go_bar <- ggplot(data = go.df, # 绘图使用的数据
                 aes(x = Description, y = Count,fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 11), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
go_bar
ggsave(go_bar,filename = "GO_up_Barplot.pdf",width = 12,height = 15)
















###################下调基因的富集分析####################################

rm(list = ls())
#相关R包载入：
library(dplyr)#数据清洗
library(org.Hs.eg.db)#ID转换
library(clusterProfiler)#富集分析
library(ggplot2)#绘图
library(RColorBrewer)#配色调整



#加载差异基因文件
DESeq2<-read.csv('../Deseq2//sig DESeq22.csv',header = T,row.names = 1)
head(DESeq2)

#添加上下调基因分组标签：
DESeq2$group <- case_when(
  DESeq2$log2FoldChange >=1 & DESeq2$padj <=0.05 ~ "up",
  DESeq2$log2FoldChange <=-1 & DESeq2$padj <= 0.05 ~ "down",
  abs(DESeq2$log2FoldChange) <= 2 ~ "none",
  DESeq2$padj >=0.05 ~ "none")
head(DESeq2)

table(DESeq2$group)

#分别筛选上调基因、下调基因或所有差异基因（上调+下调）：
up <- rownames(DESeq2)[DESeq2$group=="up"]#差异上调
down <- rownames(DESeq2)[DESeq2$group=="down"]#差异下调
diff <- c(up,down)#所有差异基因
head(up)
head(down)


#ID转换：
#查看可转换的ID类型：
columns(org.Hs.eg.db)


down_entrez <- bitr(down,
                    fromType = "SYMBOL",#现有的ID类型
                    toType = "ENTREZID",#需转换的ID类型
                    OrgDb = "org.Hs.eg.db")
#提示1.36%的基因无法映射是正常的，一般不同数据库ID转换会存在缺失
head(down_entrez)


#KEGG富集分析，我们以总差异基因(diff_entrez)为例：
KEGG_diff <- enrichKEGG(gene = down_entrez$ENTREZID,
                        organism = "hsa",#物种，Homo sapiens (human)
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

#将ENTREZID转化为可读的gene symbol
KEGG_diff<- setReadable(KEGG_diff, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

#保存富集结果
write.csv(KEGG_diff@result,file = 'KEGG_down.CSV')

#使用ggplot2进行更加自由的可视化呈现：
#先提取富集结果表前Top20：
KEGG_result<-KEGG_diff@result
KEGG_top20 <- KEGG_result[1:20,]

#指定绘图顺序（转换为因子）：
KEGG_top20$pathway <- factor(KEGG_top20$Description,levels = rev(KEGG_top20$Description))

#Top20富集数目条形图：
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))





#富集因子（Rich Factor）计算：
colnames(KEGG_top20) #方便复制所需列名

#Rich Factor = GenRatio/BgRatio
top20_rf <- apply(KEGG_top20,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  KEGG_top20_rf <- round(GeneRatio/BgRatio,2)
  KEGG_top20_rf
})
head(top20_rf)
KEGG_top20$Rich_Factor <- top20_rf


#富集因子版显著性气泡图绘制：
p3 <- ggplot(data = KEGG_top20,
             aes(x = Rich_Factor, # X轴用富集因子来映射
                 y = pathway))+
  geom_point(aes(size = Count,
                 color = -log10(pvalue)))+
  theme_bw()+
  scale_color_distiller(palette = "Spectral",direction = -1) +
  labs(x = "Rich Factor",
       y = "",
       title = "Dotplot of Enriched KEGG Pathways",
       size = "Gene Number") +
  mytheme
p3



pdf(file = 'down_kegg.pdf',width = 12,height = 15)
p3
dev.off()

############################################################


#GO(Gene Ontology)富集分析：
##MF(我们以总差异基因的GO富集为例):
GO_MF_diff <- enrichGO(gene = down_entrez$ENTREZID, #用来富集的差异基因
                       OrgDb = org.Hs.eg.db, #指定包含该物种注释信息的org包
                       ont = "MF", #可以三选一分别富集,或者"ALL"合并
                       pAdjustMethod = "BH", #多重假设检验矫正方法
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE) #是否将gene ID映射到gene name
#提取结果表格：
GO_MF_result <- GO_MF_diff@result



#CC:
GO_CC_diff <- enrichGO(gene = down_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
#提取结果表格：
GO_CC_result <- GO_CC_diff@result

#BP:
GO_BP_diff <- enrichGO(gene = down_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
#提取结果表格：
GO_BP_result <- GO_BP_diff@result



#MF、CC、BP三合一：
GO_all_diff <- enrichGO(gene = down_entrez$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = "ALL", #三合一选择“ALL”
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
#提取结果表格：
GO_all_result <- GO_all_diff@result
write.csv(GO_all_result,'GO_down_result.csv')

# 绘制GO富集分析条形图+，结果默认按qvalue升序，分别选出前十的term进行绘图即可
goBP <- subset(GO_all_result,subset = (ONTOLOGY == "BP"))[1:10,]
goCC <- subset(GO_all_result,subset = (ONTOLOGY == "CC"))[1:10,]
goMF <- subset(GO_all_result,subset = (ONTOLOGY == "MF"))[1:10,]
go.df <- rbind(goBP,goCC,goMF)
# 使画出的GO term的顺序与输入一致
go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
# 绘图
go_bar <- ggplot(data = go.df, # 绘图使用的数据
                 aes(x = Description, y = Count,fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 11), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
go_bar
ggsave(go_bar,filename = "GO_down_Barplot.pdf",width = 12,height = 15)


