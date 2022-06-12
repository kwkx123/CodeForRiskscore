#UTF-8 encoding
#首先先把所有数据集的所有模型基因提取出来，总共有68个gene
library(limma)
setwd("E:\\10分文章新训练集\\表达谱汇总")
gene_list=read.table("allmodel_genelist.txt", header=T, sep="\t", check.names=F)
fix(gene_list)
#首先是37642——96
rt=read.table("beat.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
fix(data)
ncol(data)
data=data[which(row.names(data)[]%in%gene_list[,1]),]
nrow(data)
row.names(data)
gene_list[,1]
fix(data)

# for (i in 1:ncol(data)) {
#   data[,i]=log2(data[,i]+1)
# }
write.table(data,
            file="beatdata68.txt",
            sep="\t",
            quote=F,
            row.names=T)









setwd("E:\\10分文章新训练集\\代码上传\\IPRP")####please change this path to the "IPRP" file!
load(".\\IPRPmodel.RData")
library(limma)
library(survival)
library(survivalROC)
library(survminer)
namedata=c("37642test","tcga","12417","106291","beat","fimm","all","allsva")# 7 test sets
namedata=c("beatnomo")
pvaluesurvival=c() ###the p values in K-M curves
for (ID in namedata) {
  cox_gene_list=read.table("cox_gene_list.txt", header=T, sep="\t", check.names=F)
  cox_pair_list=read.table("cox_pair_list.txt", header=T, sep="\t", check.names=F)
  expdata=read.table(paste0(ID,"data114.txt"), header=T, sep="\t", check.names=F,row.names = 1)
  out=data.frame()
  for(i in 1:nrow(cox_gene_list))
  {
    num=which(row.names(expdata)[]==cox_gene_list[i,1])
    out=rbind(out,expdata[num,])
    
  }
  #expression data is in "out", then we are going to pair these genes
  tcgaPair=data.frame()
  rt=out
  sampleNum=ncol(rt)
  for(i in 1:(nrow(rt)-1)){
    for(j in (i+1):nrow(rt)){
      pair=ifelse(rt[i,]>rt[j,], 1, 0)
      rownames(pair)=paste0(rownames(rt)[i],"|",rownames(rt)[j])
      tcgaPair=rbind(tcgaPair, pair)
    }
  }
  tcgaOut=rbind(ID=colnames(tcgaPair), tcgaPair)
  tcgaOut=tcgaOut[which(rownames(tcgaOut)[]%in%cox_pair_list[,1]),]
  tcgaOut=t(tcgaOut)
  cli=read.table(paste0(ID,"cli.txt"), header=T, sep="\t", check.names=F) ####input survival data of these sets
  cli[1,2]
  intername=intersect(cli[,1],row.names(tcgaOut))
  cli=cli[which(cli[,1]%in%intername),]
  tcgaOut=tcgaOut[intername,]
  fustat=vector(length = nrow(tcgaOut))
  futime=vector(length = nrow(tcgaOut))
  for(i in 1:nrow(tcgaOut))
  {
    
    
    num=which(cli[,1]==rownames(tcgaOut)[i])
    fustat[i]=cli[num,3]
    futime[i]=cli[num,2]
    
  }
  tcgaOut=cbind(tcgaOut,fustat,futime)
  tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365 ####we used "year" in survival time
  tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
  tcgaOut=as.data.frame(tcgaOut)
  for(i in 1:ncol(tcgaOut))
  {
    tcgaOut[,i]=as.numeric(as.character(tcgaOut[,i]))
  }
  riskScoreTest=predict(multiCox,type="risk",newdata=tcgaOut)      #we obtained the RGP risk sore of this test set
  medianTrainRisk=0.6844
  riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
  riskfile=cbind(id=rownames(cbind(tcgaOut,riskScoreTest,riskTest)),cbind(tcgaOut,riskScore=riskScoreTest,risk=riskTest))
  write.table(riskfile,
              file=paste0(".\\IPRP score\\",ID,"riskTest.txt"),
              sep="\t",
              quote=F,
              row.names=F)#save the RGP risk score of this test set in "RGP score" file.
  
  

  #to draw K-M curves using below codes
  rt=riskfile
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  
  
  
  pdf(file=paste0(".\\survival\\",ID,"survivalTest.pdf"),width=5.5,height=5)
  plot(fit, 
       lwd=2,
       col=c("red","blue"),
       xlab="Time (year)",
       ylab="Survival rate",
       main=paste("Survival curve (p=", pValue ,")",sep=""),
       mark.time=T)
  legend("topright", 
         c("high risk", "low risk"),
         lwd=2,
         col=c("red","blue"))
  dev.off()
  
  
  pvaluesurvival=c(pvaluesurvival,pValue)
  
  
  #to draw ROC curves using below codes
  rocCol=rainbow(3)
  aucText=c()
  
  pdf(file=paste0(".\\ROC\\",ID,"multiROC.pdf"),width=6,height=6)
  
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
       xlab="False positive rate", ylab="True positive rate",
       main=paste0("IPRP score in", ID),
       lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  aucText=c(aucText,paste0("1st year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  
  
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =2, method="KM")
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],
        lwd = 2)
  aucText=c(aucText,paste0("2nd year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  
  
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =3, method="KM")
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],
        lwd = 2)
  aucText=c(aucText,paste0("3rd year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  
  
  
  abline(0,1)
  
  
  legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
  dev.off()
  print(ID)

}
outdata=cbind(namedata,pvaluesurvival)
fix(outdata)






#######下面开始处理7个其他模型#############################################
namedata=c("37642test","tcga","12417","106291","beat","fimm","all","allsva")#
namedata=c("all","allsva")#
library(limma)
library(survival)
library(survivalROC)
pvaluesurvival=c()
for (ID in namedata) {
  setwd("E:\\10分文章新训练集\\代码上传\\PS29MRC")
  cox_gene_list=read.table("genelist.txt", header=T, sep="\t", check.names=F)
  #cox_gene_list=read.table("106291securitygenelist.txt", header=T, sep="\t", check.names=F)
  rawexpression=read.table(paste0(ID,"data114.txt"), header=T, sep="\t", check.names=F,row.names = 1)
  #rawexpression=read.table("GSEbeatsecurity.txt", header=T, sep="\t", check.names=F,row.names = 1)
  out=data.frame()
  for(i in 1:nrow(cox_gene_list))
  {
    num=which(row.names(rawexpression)[]==cox_gene_list[i,1])
    out=rbind(out,rawexpression[num,])
    
  }
  
  tcgaOut=out

  tcgaOut=t(tcgaOut)
  tcgaOut=as.data.frame(tcgaOut)
  cli=read.table(paste0(ID,"cli.txt"), header=T, sep="\t", check.names=F) ####input survival data of these sets
  intername=intersect(cli[,1],row.names(tcgaOut))
  cli=cli[which(cli[,1]%in%intername),]
  tcgaOut=tcgaOut[intername,]
  
  Targetcli=cli
  fix(Targetcli)
  fustat=vector(length = nrow(tcgaOut))
  futime=vector(length = nrow(tcgaOut))
  for(i in 1:nrow(tcgaOut))
  {
    
    num=which(Targetcli[,1]==rownames(tcgaOut)[i])
    fustat[i]=Targetcli[num,3]
    futime[i]=Targetcli[num,2]
    
    #print(i)
  }
  tcgaOut=cbind(tcgaOut,fustat,futime)
  tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
  tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
  tcgaOut=as.data.frame(tcgaOut)
  for(i in 1:ncol(tcgaOut))
  {
    tcgaOut[,i]=as.numeric(tcgaOut[,i])
  }
  
  #####下面开始添加riskScore
  riskScore=vector(length = nrow(tcgaOut))
  for(i in 1:nrow(tcgaOut))
  {
    sumvalue=0
    for(j in 1:nrow(cox_gene_list))
    {
      sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
    }
    riskScore[i]=sumvalue
    
  }
  
  ############分组的中位数为42.3010
  tcgaOut=cbind(tcgaOut,riskScore)
  risk=vector(length = nrow(tcgaOut))
  for(i in 1:nrow(tcgaOut))
  {
    if(tcgaOut[i,(nrow(cox_gene_list)+3)]<= median(riskScore))
    {
      risk[i]="low"
    }
    else
    {
      risk[i]="high"
    }
  }
  tcgaOut=cbind(tcgaOut,risk)
  write.table(tcgaOut,
              file=paste0(ID,"riskTest.txt"),
              sep="\t",
              quote = F,
              row.names=T)
  
  
  
  
  #绘制test组生存曲线
  rt=tcgaOut
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     title=ID,
                     pval=pValue,
                     pval.size=6,
                     conf.int=F,
                     legend.title="risk group",
                     legend.labs=c("high","low"),
                     font.legend=12,
                     fontsize=4,
                     xlab="Time(years)",
                     ylab="Overall survival",
                     break.time.by = 2,
                     palette=c("red","blue"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.height=.25)
  pdf(file=paste0(ID,"survivalTest.pdf"),width=5.5,height=5)#we can get K-M curve of this test set in "survival" file.
  print(surPlot)
  dev.off()
  
  
  pvaluesurvival=c(pvaluesurvival,pValue)
  
  
  
  #rt$futime=rt$futime/365
  rocCol=rainbow(3)
  aucText=c()
  
  #????risk score??ROC????
  pdf(file=paste0(ID,"multiROC.pdf"),width=6,height=6)
  
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
       xlab="False positive rate", ylab="True positive rate",
       main="Signature in TCGA",
       lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  aucText=c(aucText,paste0("1st year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  
  
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =2, method="KM")
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],
        lwd = 2)
  aucText=c(aucText,paste0("2nd year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  
  
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =3, method="KM")
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],
        lwd = 2)
  aucText=c(aucText,paste0("3rd year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  
  
  
  abline(0,1)
  
  
  legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
  dev.off()
  print(ID)
  
}
outdata=cbind(namedata,pvaluesurvival)
fix(outdata)



####下面使用IPRP 对treatment分组进行测试
setwd("E:\\10分文章新训练集\\代码上传\\IPRPtreatment作图专用")####please change this path to the "IPRP" file!
load(".\\IPRPmodel.RData")
library(limma)
library(survival)
library(survivalROC)
library(survminer)
library(timeROC)
namedata=c("37642test","tcga","12417","106291","beat","fimm","all","allsva")# 7 test sets
namedata=c("tcga")
classche="induct"
pvaluesurvival=c() ###the p values in K-M curves
for (ID in namedata) {
  cox_gene_list=read.table("cox_gene_list.txt", header=T, sep="\t", check.names=F)
  cox_pair_list=read.table("cox_pair_list.txt", header=T, sep="\t", check.names=F)
  expdata=read.table(paste0(ID,"data114.txt"), header=T, sep="\t", check.names=F,row.names = 1)
  out=data.frame()
  for(i in 1:nrow(cox_gene_list))
  {
    num=which(row.names(expdata)[]==cox_gene_list[i,1])
    out=rbind(out,expdata[num,])
    
  }
  #expression data is in "out", then we are going to pair these genes
  tcgaPair=data.frame()
  rt=out
  sampleNum=ncol(rt)
  for(i in 1:(nrow(rt)-1)){
    for(j in (i+1):nrow(rt)){
      pair=ifelse(rt[i,]>rt[j,], 1, 0)
      rownames(pair)=paste0(rownames(rt)[i],"|",rownames(rt)[j])
      tcgaPair=rbind(tcgaPair, pair)
    }
  }
  tcgaOut=rbind(ID=colnames(tcgaPair), tcgaPair)
  tcgaOut=tcgaOut[which(rownames(tcgaOut)[]%in%cox_pair_list[,1]),]
  tcgaOut=t(tcgaOut)
  cli=read.table(paste0(ID,"cli.txt"), header=T, sep="\t", check.names=F) ####input survival data of these sets
  cli[1,2]
  chemoid=read.table(paste0(ID,classche,".csv"), header=T, sep=",",row.names = 1, check.names=F) 
  intername0=intersect(cli[,1],row.names(tcgaOut))
  intername=intersect(intername0,row.names(chemoid))
  cli=cli[which(cli[,1]%in%intername),]
  tcgaOut=tcgaOut[intername,]
  fustat=vector(length = nrow(tcgaOut))
  futime=vector(length = nrow(tcgaOut))
  for(i in 1:nrow(tcgaOut))
  {
    
    
    num=which(cli[,1]==rownames(tcgaOut)[i])
    fustat[i]=cli[num,3]
    futime[i]=cli[num,2]
    
  }
  tcgaOut=cbind(tcgaOut,fustat,futime)
  tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365 ####we used "year" in survival time
  tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
  tcgaOut=as.data.frame(tcgaOut)
  for(i in 1:ncol(tcgaOut))
  {
    tcgaOut[,i]=as.numeric(as.character(tcgaOut[,i]))
  }
  riskScoreTest=predict(multiCox,type="risk",newdata=tcgaOut)      #we obtained the RGP risk sore of this test set
  medianTrainRisk=0.6844
  riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
  riskfile=cbind(id=rownames(cbind(tcgaOut,riskScoreTest,riskTest)),cbind(tcgaOut,riskScore=riskScoreTest,risk=riskTest))
  write.table(riskfile,
              file=paste0(".\\IPRP score\\",ID,classche,"riskTest.txt"),
              sep="\t",
              quote=F,
              row.names=F)#save the RGP risk score of this test set in "RGP score" file.
  
  
  
  #to draw K-M curves using below codes
  rt=riskfile
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  
  
  
  pdf(file=paste0(".\\survival\\",ID,classche,"survivalTest.pdf"),width=5.5,height=5)
  plot(fit, 
       lwd=2,
       col=c("red","blue"),
       xlab="Time (year)",
       ylab="Survival rate",
       main=paste("Survival curve (p=", pValue ,")",sep=""),
       mark.time=T)
  legend("topright", 
         c("high risk", "low risk"),
         lwd=2,
         col=c("red","blue"))
  dev.off()
  
  
  pvaluesurvival=c(pvaluesurvival,pValue)
  
  
  #to draw ROC curves using below codes

  ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
                 marker=rt$riskScore, cause=1,
                 weighting='aalen',
                 times=c(1,2,3), ROC=TRUE)
  pdf(file=paste0(".\\ROC\\",ID,classche,"multiROC.pdf"),width=6,height=6)
  
  plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
  plot(ROC_rt,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
          c(paste0('1st year AUC=',sprintf("%.03f",ROC_rt$AUC[1]),"*"),
            paste0('2nd year AUC=',sprintf("%.03f",ROC_rt$AUC[2]),"*"),
            paste0('3rd year AUC=',sprintf("%.03f",ROC_rt$AUC[3]),"*")
          ),
         col=c("green","blue",'red'),lwd=2,bty = 'n')
  
  dev.off()
  print(ID)
  # c(paste0('1st year AUC=',sprintf("%.03f",ROC_rt$AUC[1]),"*"),
  #   paste0('2nd year AUC=',sprintf("%.03f",ROC_rt$AUC[2]),"*"),
  #   paste0('3rd year AUC=',sprintf("%.03f",ROC_rt$AUC[3]),"*")
  # )
}
outdata=cbind(namedata,pvaluesurvival)
fix(outdata)



