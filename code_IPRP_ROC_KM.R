#UTF-8 encoding
setwd("E:\\white\\LSC17")
library(survival)
library(survivalROC)
library(survminer)
library(limma)
library(survival)
library(survivalROC)
cox_gene_list=read.table("genelist.txt", header=T, sep="\t", check.names=F)

rt=read.table("BeatAML65patients.txt", header=T, sep="\t", check.names=F)#read the expression data
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rawexpression=avereps(data)
fix(rawexpression)




out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(row.names(rawexpression)[]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
}
  
tcgaOut=out
tcgaOut=t(tcgaOut)
row.names(tcgaOut)=colnames(rawexpression)
tcgaOut=as.data.frame(tcgaOut)
cli=read.table("clinical_65patients.txt", header=T, sep="\t", check.names=F) ####input survival data of these sets
intername=intersect(cli[,1],row.names(tcgaOut))
cli=cli[which(cli[,1]%in%intername),]
tcgaOut=tcgaOut[intername,]
Targetcli=cli
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  num=which(Targetcli[,1]==rownames(tcgaOut)[i])
  fustat[i]=Targetcli[num,3]
  futime[i]=Targetcli[num,2]
}
tcgaOut=cbind(tcgaOut,fustat,futime)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
for(i in 1:ncol(tcgaOut))
{
  tcgaOut[,i]=as.numeric(tcgaOut[,i])
}
  
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
              file="riskTest.txt",
              sep="\t",
              quote = F,
              row.names=T)
  
  
  
  
  #survival plot
rt=tcgaOut
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

surPlot=ggsurvplot(fit, 
                   data=rt,
                   title="BeatAML",
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
pdf(file="survivalTest.pdf",width=5.5,height=5)#we can get K-M curve of this test set in "survival" file.
print(surPlot)
dev.off()
  

###ROC curve
rocCol=rainbow(3)
aucText=c()

pdf(file="multiROC.pdf",width=6,height=6)

par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)

roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     main="Signature in BeatAML",
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


