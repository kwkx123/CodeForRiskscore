####set path to the "XX" file, please change the path manually
setwd("E:\\prognosis\\code_for_riskscore")
#######Firstly using limma to average duplicate genes
library(limma)
rt=read.table("expression.csv", header=T, sep=",", check.names=F)#read the expression data
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)###56633 genes left
fix(data)####have a look at the data


cox_gene_list=read.table("genelist.txt", header=T, sep="\t", check.names=F)
###read the gene and coeffficient of LSC17.Note! check whether all genes of ths signature were in the expression data



out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(row.names(data)[]==cox_gene_list[i,1])
  out=rbind(out,data[num,])
}
row.names(out)=cox_gene_list[,1]
colnames(out)=colnames(data)
  
out=t(out)
out=as.data.frame(out)
  
  
riskScore=vector(length = nrow(out))  ###TO save riskscore in "riskScore" 
for(i in 1:nrow(out))
{
  sumvalue=0
  for(j in 1:nrow(cox_gene_list))
  {
    sumvalue=sumvalue+out[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
}
  
out=cbind(out,riskScore)
write.table(out,
              file="riskScore.csv",
              sep=",",
              quote = F,
              row.names=T)#save the riskScore
  
