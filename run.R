library(randomForest)
library(ggplot2)
library(pheatmap)
# 训练模型
# 读取分组
design = read.table("3750group3-1.txt",header = T, row.names = 1)
#将$Type列转换成因子factor
design$Type=as.factor(design$Type)
# 读取domain矩阵
otu_table = read.table("train_matrix_fixed.txt",header = T, row.names = 1)
otu_table = t(otu_table)
# 取出"group1"(训练集)
design_sub = subset(design, Group %in% c("group1"))
summary(design_sub)


#判断"design_sub"的行名（细菌编号）是否在结构域矩阵"out_table"中
idx = rownames(design_sub) %in% colnames(otu_table)
# 取出Group1对应的数据集
design_sub = design_sub[idx,]
otu_sub = otu_table[, rownames(design_sub)]
summary(design_sub)


# 训练构建模型
set.seed(3)
rf = randomForest(t(otu_sub), design_sub$Type, importance=TRUE, proximity=T, ntree = 1000)
print(rf)
# 保存模型
save(rf,file = 'bestmodel_group3_1.RData')
# 直接加载模型
load('bestmodel_group3_1.RData')
# 交叉验证选择Features
set.seed(827) # 随机数据保证结果可重复，必须
# rfcv是随机森林交叉验证函数：Random Forest Cross Validation
result = rfcv(t(otu_sub), design_sub$Type, cv.fold=5)
save(result,file = 'best_rfcv_group3_1.RData')

# 查看错误率表，68时错误率最低，为最佳模型
result$error.cv
error_data <- as.data.frame(result$error.cv)
write.table(error_data,file = 'best_error_group3_1.txt',sep = '\t',row.names = T,
            quote = F,col.names = NA)
# 绘制验证结果 
with(result,plot(n.var, error.cv, log="x", type="o", lwd=1))# 交叉验证的结果建议多做5-6次，将结果统一在一张图上

# 导出训练集观测值与预测结果
train.p = predict(rf, type = "response")
df = data.frame(observed = design_sub$Type, predict = train.p)  #提取数据的子集作为另一部分数据再预测一下
# 保存预测结果与真实结果比较
write.table(df,file = "train_predict_group3_1.txt",quote = F,sep = '\t', row.names = T, col.names = T)



# 导出feature重要性
imp= as.data.frame(rf$importance)
imp = imp[order(imp[,1],decreasing = T),]
head(imp,n=10)
write.table(imp,file = "best_importance_class_group3_1.txt",quote = F,sep = '\t', row.names = T, col.names = T)
# 简单可视化
varImpPlot(rf, main = "Top 10 - Feature importance",n.var = 10, bg = par("bg"), color = par("fg"), gcolor = par("fg"), lcolor = "gray" )



# ggplot2美化feature贡献度
# 读取所有feature贡献度
imp = read.table("best_importance_class_group3_1.txt", header=T, row.names= 1, sep="\t") 
# 分析选择top23分组效果最好
imp = head(imp, n=23)
# 反向排序X轴，让柱状图从上往下画
imp = imp[order(1:23,decreasing = T),]
#将imp的按第三列从小到大排序
imp = imp[order(imp[,3]),]
# 取出列名
imp$Domain = gsub("","",rownames(imp),perl=TRUE) 

imp$Domain=factor(imp$Domain,levels = imp$Domain)

# 图1. feature重要性柱状图
library(ggplot2)
p=ggplot(data = imp, mapping = aes(x=Domain,y=MeanDecreaseAccuracy,fill=Domain)) + 
  geom_bar(stat="identity")+coord_flip()+theme_bw()
ggsave(p,filename = "imp_shape.pdf",width = 16,height = 9)


# group2验证
# design = read.table("group_Gr_all.txt",header = T, row.names = 1)
design_test = subset(design, Group %in% c("group2")) 
summary(design_test)
idx = rownames(design_test) %in% colnames(otu_table)
design_test = design_test[idx,]
otu_sub = otu_table[,rownames(design_test)]
summary(design_test)
# 转置，并添加分组信息
otutab_t = as.data.frame(t(otu_sub))
# 将Group2的分组信息添加到domain矩阵中
# 表示按照行名将矩阵"design"中"Group"列的内容添加到矩阵"otutab_t"中并将列名设置为"otutab_t$Group中的"Group",
otutab_t$Type = design[rownames(otutab_t),]$Type

set.seed(13)
otutab.pred = predict(rf, t(otu_sub) )  
pre_tab = table(observed=otutab_t[,"Type"],predicted=otutab.pred) 
pre_tab
# 整理样本原始分组和预测分类
predict = data.frame(Type = otutab_t[,"Type"], predicted=otutab.pred)
# 保存预测结果表
write.table("SampleID\t", file=paste("RF_prediction_binary.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
write.table(predict, file = "RF_prediction_binary.txt",append = T, quote = F, row.names = T, col.names = T, sep = "\t")
# 此处应有报错
# Warning message:
# In write.table(predict, file = "RF_prediction_binary.txt", append = T,  :
#                appending column names to file

