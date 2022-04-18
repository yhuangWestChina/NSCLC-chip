#univariate receiver operating characteristic curve (ROC) analysis
library(pROC)
dat<-read.table("variable_uROC.txt",head=T,row.names=1,sep="\\t")
dat<-data.frame(dat)
uroc<-auc(dat\$Group,dat\$$name)
write.table(file="variable.uROC.result.txt",x=uroc)

#multivariate ROC analysis and machine learning analysis
library(randomForest)
library(e1071)
library(ROCR)
library(AUCRF)
library(dplyr)
dat <- read.table("gwas.result.txt",head=T,row.names=1,sep="\t")
dat <- read.table("gwas.result.txt",head=T,sep="\t")
dat <- data.frame(dat)
control_data <- data.frame(subset(dat,Group==1))
case_data <- data.frame(subset(dat,Group==0))
control_data$Group=as.factor(control_data$Group)
case_data$Group=as.factor(case_data$Group)

## step 1 - spliting data into two parts
train_case_data=na.omit(sample_n(case_data,size=138))
test_case_data <- na.omit(case_data[!case_data$Sample %in% train_case_data$Sample,])
train_control_data=na.omit(sample_n(control_data,size=226))
test_control_data <- na.omit(control_data[!control_data$Sample %in% train_control_data$Sample,])
train_data <- na.omit(rbind(train_case_data,train_control_data))
rownames(train_control_data)<- train_control_data$Sample
train_data$Sample <- NULL
test_data <- na.omit(rbind(test_case_data,test_control_data))
rownames(test_data)<- test_data$Sample
test_data$Sample <- NULL

## AUCRF to evaluate the significance of single variable
aucrf_dat <- dat
aucrf_dat$Group <- as.factor(aucrf_dat$Group)
fit <- AUCRF(Group~.,data=aucrf_dat,na.action=na.omit,ntree=20000)
fit$ranking
write.table(file="all_variables.sig.snps.AUCRF.ranking.txt",fit$ranking,sep="\t")

#randomForest:

nsclc_rf <- randomForest(Group~.,data=train_data,ntree=20000,mtry=3,importance=T,proximity=T,na.action = na.omit)
pdf("all.variable_importance.pdf")
varImpPlot(nsclc_rf,main="variable importance")
dev.off()
pre_ran <- predict(nsclc_rf,newdata=test_data)
obs_p_ran <- data.frame(prob=pre_ran,obs=test_data$Group)
table(test_data$Group,pre_ran,dnn=c("True","Prediction"))
temp_conf <- confusionMatrix(pre_ran,test_data$Group)
write.table(file="all.rf.conf.txt",x=temp_conf$overall,sep="\t",quote=F)
pre_ran <- predict(nsclc_rf,newdata=test_data,type="prob")
pred <- prediction(pre_ran[,2],test_data$Group)
perf <- performance(pred,"tpr","fpr")
pdf("all.AUCRF.rf.perf.AUC.pdf")
plot(perf,colorize=T,lwd=10,colorkey.relwidth=0.1)
dev.off()
write.table(file="all.AUCRF.rf.prob.txt",x=pre_ran,sep="\t")

#SVM:
nsclc_svm <- svm(Group ~.,data=train_data,type = 'C',kernel = 'radial',probability=TRUE)
pre_svm <- predict(nsclc_svm,newdata=test_data,probability=T)
obs_p_svm = data.frame(prob=pre_svm,obs=test_data$Group)
svm_roc <- roc(test_data$Group,as.numeric(pre_svm))
temp_conf <- confusionMatrix(pre_svm,test_data$Group)
write.table(file="all.svm.conf.txt",x=temp_conf$overall,sep="\t",quote=F)
pred.prob <- attr(pre_svm,"probabilities")
pred <- prediction(pred.prob[,2],test_data$Group)
perf <- performance(pred,"tpr","fpr")
pdf("all.AUCRF.svm.perf.AUC.pdf")
plot(perf,colorize=T,lwd=10,colorkey.relwidth=0.1)
dev.off()
write.table(file="all.AUCRF.svm.prob.txt",x=pred.prob,sep="\t")
