library(Mime1)
library(limma)
library(tidyverse)
library(data.table)
library(survival)
library(survminer)



TCGA = read.table("TCGA_TPM.txt", header = T, sep = "\t", check.names = F,row.names = 1)
TCGA.time = read.table("time.TCGA.txt", header = T, sep = "\t", check.names = F,row.names = 1)

GSE76427 = read.table("GSE76427.txt", header = T, sep = "\t", check.names = F,row.names = 1)
GSE76427.time = read.table("time.GSE76427.txt", header = T, sep = "\t", check.names = F,row.names = 1)
#标准化
GSE76427=normalizeBetweenArrays(GSE76427)

genelist = read.table("Genelist.txt", header = F, sep = "\t", check.names = F)


TCGA.gene = rownames(TCGA)[rowMeans(TCGA)>1]
GSE76427.gene = rownames(GSE76427)[rowMeans(GSE76427)>1]
samegene <- Reduce(intersect, list(TCGA.gene,GSE76427.gene,genelist[,1]))
length(samegene)

group=sapply(strsplit(colnames(TCGA),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
TCGA = TCGA[,group == 0]

TCGA=t(TCGA)
rownames(TCGA)=substr(rownames(TCGA),1,12)
sameSample=intersect(row.names(TCGA), row.names(TCGA.time))
TCGA=TCGA[sameSample,]
TCGA.time=TCGA.time[sameSample,]
TCGA = log2(TCGA +1)
TCGA=cbind(rownames(TCGA),TCGA.time,TCGA)
colnames(TCGA)[1] = "ID"
TCGA[1:7,1:7]

GSE76427 = t(GSE76427)
sameSample=intersect(row.names(GSE76427), row.names(GSE76427.time))
GSE76427=GSE76427[sameSample,]
GSE76427.time=GSE76427.time[sameSample,]
GSE76427=cbind(rownames(GSE76427),GSE76427.time,GSE76427)
colnames(GSE76427)[1] = "ID"
GSE76427[1:7,1:7]

list_train_vali_Data = list(TCGA,GSE76427)
names(list_train_vali_Data) = c("TCGA","GSE76427")

res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data[[1]],
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = samegene,
                       mode = 'all',nodesize =5,seed = 5201314)


saveRDS(res,"res.rds")

res = readRDS("res.rds")


pdf(file="1.cindex_dis_all.pdf", width=8, height=15)
cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],order =names(list_train_vali_Data),width = 0.35)
dev.off()

final_model = "RSF + SuperPC"

write.table(res[["Cindex.res"]], file="Cindex.res.txt", sep="\t", quote=F, col.names=T,row.names = F)
write.table(res[["Sig.genes"]], file="unicox.txt", sep="\t", quote=F, col.names=F,row.names = F)

RS = data.frame()
for (i in seq_along(res[["riskscore"]])) {
  for (l in c(1:length(list_train_vali_Data))) {
    rt = cbind(res[["riskscore"]][[i]][[l]][,c(1,4)],rep(names(res[["riskscore"]][i]),length(rownames(res[["riskscore"]][[i]][[l]][,c(1,4)]))))
    colnames(rt)[3] = "model"
    RS = rbind(RS,rt)
  }
}
write.table(RS, file="riskscore.txt", sep="\t", quote=F, col.names=T,row.names = F)

pdf(file="2.cindex_dis_select.pdf", width=8, height=5)
cindex_dis_select(res,
                  model=final_model,
                  order= names(list_train_vali_Data))
dev.off()


pdf(file=("6auc_dis_select.pdf"), width=8, height=4)
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name=final_model,
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,3,5))
dev.off()

