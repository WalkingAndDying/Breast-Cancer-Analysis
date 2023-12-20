rm(list = ls())


#相关R包和数据载入：
library(stringr)
library(caret)
library(survminer)
library(survival)
library(pheatmap)
#载入TCGA表达矩阵和临床信息表格(第1期整理)：
setwd('C:\\Users\\user\\Desktop\\乳腺癌\\kegg_go')
exp<-read.table('../Deseq2/Deseq2_rld.txt',sep = '\t',header = T,row.names = 1)

diff<-read.csv('../Deseq2/sig DESeq22.csv')
table(abs(diff$log2FoldChange)>2)
diff<-diff[abs(diff$log2FoldChange)>2,]

#矩阵自检，如果是character则需转换为numeric：

exp<-exp[diff$genename,]

clinical<-read.csv('../Deseq2/clinical.csv')
clinical<-clinical[,c(2,9:11)]




clinical$OS_STATUS<-ifelse(clinical$OS_STATUS=="0:LIVING",0,1)
colnames(clinical)[2:3]<-c('status','OS')

clinical$OS<-clinical$OS/12
a<-!clinical$OS==0
clinical<-clinical[a,]

exp<-exp[,a]

cox <- apply(
  exp,1,function(x){
    clinical$gene <- as.numeric(x)
    cox_genes <- coxph(Surv(OS, status) ~ gene, data =clinical)
    coef <- coef(cox_genes) #回归系数
    SE <- sqrt(diag(vcov(cox_genes))) #标准误
    HR <- exp(coef) #风险比
    cox_need <- cbind(HR = HR,
                      HR.95L = exp(coef - qnorm(.975, 0, 1) * SE),
                      HR.95H = exp(coef + qnorm(.975, 0, 1) * SE),
                      pvalue = 1 - pchisq((coef/SE)^2, 1))
    return(cox_need['gene',])
  }
)
unicox <- t(cox)
head(unicox)
#提取预后差异基因列表(p<0.05)：
diff_unicox <- unicox[unicox[,4]<0.01,]
dim(diff_unicox) #从66个差异基因中筛选出了25个候选预后DEHGs
table(diff_unicox[,1]<1) #18个风险因子，7个保护因子
head(diff_unicox)




#继续载入lasso回归所需R包：
library(glmnet)


#构建候选预后DEHGs表达矩阵：
exp <- exp[rownames(diff_unicox),]
#使用lasso回归进一步收缩单因素cox筛选的25个候选DEHGs：
x <- t(exp)
y <- data.matrix(Surv(time = clinical$OS,event = clinical$status))



x[1:6,1:6]
head(y)

#构建模型：
fit <- glmnet(x, y, family = 'cox', type.measure = "deviance", nfolds = 10)
plot(fit,xvar = 'lambda',label = T) #候选DEHGs的lasso系数



#十折交叉检验筛选最佳lambda：
set.seed(007)
lasso_fit <- cv.glmnet(x, y, family = 'cox', type.measure = 'deviance', nfolds = 10)
plot(lasso_fit)


#提取最佳λ值(这里选择1se对应lambda)：
lambda.1se <- lasso_fit$lambda.1se
lambda.1se


#使用1se的lambda重新建模：
model_lasso_1se <- glmnet(x, y, family = 'cox', type.measure = 'deviance', nfolds = 10,lambda = lambda.1se)
#拎出建模使用基因：
gene_1se <- rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]#as.numeric后"."会转化为0
gene_1se #筛选出7个


#######################多因素COX#########################

exp<-as.data.frame(t(exp))
exp <- as.data.frame(exp[,gene_1se])
dt <- cbind(clinical,exp)
#构建最优多因素cox比例风险回归模型：
colnames(dt)[6]<-'MAPT_AS1'
train_cox <- coxph(Surv(OS, status) ~ NTRK3 + MAPT_AS1 + SEMA3B + FABP7 + RPH3A+CEL+TFPI2+CLCA4+IYD+CDHR4+QRFPR+PAX7+OR52E6, data = dt)
train_cox

coef <- coef(train_cox)
coef

write.table(coef,file = 'coef.txt',sep = '\t')

#Step1：先将每个基因表达量*对应系数（注意顺序，基因和对应系数不能乱）：
x <- data.frame(exp$NTRK3*coef[1],
                dt$MAPT_AS1*coef[2],
                exp$SEMA3B**coef[3],
                exp$FABP7*coef[4],
                exp$RPH3A*coef[5],
                exp$CEL*coef[6],
                exp$TFPI2*coef[7],
                exp$CLCA4*coef[8],
                exp$IYD*coef[9],
                exp$CDHR4*coef[10],
                exp$QRFPR*coef[11],
                exp$PAX7*coef[12],
                exp$OR52E6*coef[13])
colnames(x) <- names(coef)
head(x)


#Step2:将每行相加即为每个样本的风险评分：
dt$score <- apply(x,1,sum) #相加，并将风险评分列添加到训练集矩阵备用



dt <- dt[order(dt$score,decreasing = F),] #按风险评分从低到高排序
dt$id <- c(1:length(dt$score)) #根据调整后顺序建立编号id
dt$status <- ifelse(dt$status==0,'alive','death') #0/1转换为字符生死
dt$status <- factor(dt$status,levels = c("death","alive")) #指定因子调整顺序
dt$Risk_Group <- ifelse(dt$score<median(dt$score),'Low Risk','High Risk') #将风险评分按中位数拆分为高/低风险两组
dt$Risk_Group <- factor(dt$Risk_Group,levels = c('Low Risk','High Risk')) #指定顺序
head(dt)


exp <- dt[,5:17] #提取表达矩阵,并调整一下顺序(按风险因子和保护因子排列)
head(exp)


p1 <- ggplot(dt,aes(x = id,y = score)) +
  geom_point(aes(col = Risk_Group)) +
  scale_colour_manual(values = c("blue","red")) +
  geom_hline(yintercept = median(dt$score), colour="grey", linetype="dashed", size=0.8) +
  geom_vline(xintercept = sum(dt$Risk_Group == "Low Risk"), colour="grey", linetype = "dashed", size = 0.8) +
  theme_bw()
p1



#1.2.患者生存时间散点图绘制：
p2 <- ggplot(dt,aes(x = id,y = OS)) +
  geom_point(aes(col = status)) +
  scale_colour_manual(values = c("red","blue")) +
  geom_vline(xintercept = sum(dt$Risk_Group == "Low Risk"), colour = "grey", linetype = "dashed", size = 0.8) +
  theme_bw()
p2



#1.3.表达量热图绘制：
mycol <- colorRampPalette(c("blue","yellow","red"))(100) #自定义颜色
exp2 <- t(scale(exp)) #矩阵标准化：
exp2[1:5,1:4]


#添加分组信息：
annotation <- data.frame(Type = as.vector(dt$Risk_Group))
rownames(annotation) <- rownames(exp)
annotation$Type <- factor(annotation$Type,levels = c('Low Risk','High Risk'))
head(annotation)


ann_colors <- list(Type = c('Low Risk' = "blue",
                            'High Risk' = "red")) #添加分组颜色信息
#绘图：
pheatmap(exp2,
         col= mycol,
         cluster_rows = F,
         cluster_cols = F,
         show_colnames = F,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         annotation_legend = F
)


#将热图转化为ggplot2对象：
library(ggplotify)
p3 <- as.ggplot(as.grob(pheatmap(exp2,
                                 col= mycol,
                                 cluster_rows = F,
                                 cluster_cols = F,
                                 show_colnames = F,
                                 annotation_col = annotation,
                                 annotation_colors = ann_colors,
                                 annotation_legend = F
)))

p4 <- p3+theme(plot.margin=unit(rep(0.5,4),'cm')) #自行调整尝试，尽量和p1p2中心对齐
p4
#拼图：

library(cowplot)
plot_grid(p1,p2,p4, nrow = 3, align = "v", axis = "tlbr") #可以看到直接拼在一起热图是对不齐的




#相关R包载入：
library(survminer)
library(survival)
#重新载入数据：



dt$risk <- ifelse(dt$score > median(dt$score),"High","Low") #将风险评分按中位数拆分为高/低风险两组
#2.生存分析绘制：
dt$status <- ifelse(dt$status=='alive',0,1) 

fit <- survfit(Surv(OS, status) ~ risk, data = dt)
fit
ggsurvplot(
  fit,
  data = dt,
  censor = T, #是否绘制删失点
  censor.shape = "|", censor.size = 4,
  conf.int = TRUE, #是否显示置信区间
  conf.int.style = "ribbon",
  conf.int.alpha = 0.3,
  pval = TRUE, #是否显示P值
  pval.size = 5,
  legend = "top", #图例位置
  legend.title = 'Risk Score',
  legend.labs = c("High Risk","Low Risk"),
  xlab = "Years",
  ylab = "Survival probablity",
  title = "Discovery TCGA Cohort",
  palette = c('#ed0000','#00468b'), #调用色板or自行创建
  ggtheme = theme_bw(), #主题修改
  risk.table = TRUE, #是否风险表添加
  risk.table.col = "strata", #颜色跟随曲线
  risk.table.title = 'Number at risk',
  fontsize = 4,
  risk.table.y.text = FALSE, #是否显示风险表y轴标签
  risk.table.height = 0.2,
)

