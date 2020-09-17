setwd("~/Dropbox/Kim/cuffnormformycellsJuly2018/")
require(RColorBrewer)
library(RColorBrewer)
require(pcaMethods)
require(qvalue)
require(limma)
require(marray)
require(sva)
require(Gviz)
require(genefu)

genes = read.delim(file = "/Users/staufferkm/Dropbox/Kim/cuffnormformycellsJuly2018/mycells.July2018.genes.entrez.attr_table", header = T, as.is = T, sep = "\t", quote = "") ###fix this to be entrez
gene.exp= read.table(file ="/Users/staufferkm/Dropbox/Kim/cuffnormformycellsJuly2018/mycells.July2018.genes.fpkm_table.fixed", header = T, as.is = T, sep = "\t", quote = "")
tar = read.delim(file = "/Users/staufferkm/Dropbox/Kim/cuffnormformycellsJuly2018/mycells.July2018.samples.table", header = T, as.is = T, sep = "\t", quote = "")

####fix this from here on
gene.exp1=gene.exp[,-1]
rownames(gene.exp1) = gene.exp[,1]
m=apply(gene.exp1, 1, mean)
hist(log2(m+.5), breaks =250)
high=which(m>1)
genes1 = genes[high,]
gene.exp2 = gene.exp1[high,]
gene.exp2 = data.matrix(gene.exp2)
gene.exp3 = log2(gene.exp2+.5)

# names = c("MCF7shMLL2P4", "MCF7shMLL2P5", "MCF7shLucifP10", "MCF7shKDM6AP624", "ZR751shMLL2P3", "ZR751shMLL3P7", "MCF7shMLL3P8", "T47DshKDM6AP8", "MCF7shMLL3P7", "ZR751shMLL2P6", "MCF7shLucifP11", "T47DshLucifP13", "ZR751shKDM6AP4", "T47DshMLL2P6", "ZR751shKDM6AP2", "ZR751shMLL3P5", "T47DshLucifP8", "T47DshMLL2P5", "ZR751shLucifP3", "T47DshKDM6AP5", "ZR751shLucifP4", "MCF7shKDM6AP6", "T47DshMLL3P13", "T47DshMLL3P12")
# colnames(gene.exp3) <- names
# head(gene.exp3)
# 
# MCF7 = c(1,2,3,4,7,9,11,22)
# ZR751 = c(5,6,10,13,15,16,19,21)
# T47D = c(8,12,14,17,18,20,23,24)
# Z.MLL3.Lucif = c(6,16,19,21)
gene.exp3.Z.MLL3.Lucif = gene.exp3[,Z.MLL3.Lucif]
names = c("ZR751shMLL3P7", "ZR751shMLL3P5", "ZR751shLucifP3", "ZR751shLucifP3")
colnames(gene.exp3.Z.MLL3.Lucif) <- names
tar.Z.MLL3.Lucif = tar[Z.MLL3.Lucif,]
dim(gene.exp3.Z.MLL3.Lucif)
dim (tar.Z.MLL3.Lucif)
######################PICK UP HERE ON AUG 30 2018

# gene.exp3.MCF7 = gene.exp3[,MCF7]
# gene.exp3.ZR751 = gene.exp3[,ZR751]
# gene.exp3.T47D = gene.exp3[,T47D]
# gene.exp3.Z.MLL3.Lucif = gene.exp3[,Z.MLL3.Lucif]
# tar.MCF7 = tar[MCF7,]
# tar.ZR751 = tar[ZR751,]
# tar.T47D = tar[T47D,]
# tar.Z.MLL3.Lucif = tar[Z.MLL3.Lucif,]
tar.Z.MLL3.Lucif$KD = NA
tar.Z.MLL3.Lucif$KD = c("KD", "KD", "Con", "Con")



cd=cor(gene.exp3.Z.MLL3.Lucif)
cna=apply(cd,2,function(x) replace(x,which(x==1),NA))
dt=as.dist(1-cd)
h=hclust(dt,method = "ward.D")
lab.palette=colorRampPalette(rev(brewer.pal(9,"Spectral")), space = "Lab")
lp<-lab.palette(100)
heatmap(cna,Colv=as.dendrogram(h),Rowv=as.dendrogram(h),scale="row", col = lp, labRow=colnames(gene.exp3.Z.MLL3.Lucif), labCol=colnames(gene.exp1.Z.MLL3.Lucif))


# ###let's see if I can get rid of more low values so that SVA runs without error
# m1=apply(gene.exp3.Z.MLL3.Lucif, 1, mean)
# hist(m, breaks =1000,xlim=c(-2,15))
# highs=which(m1>(-1))
# gene.exp4.Z.MLL3.Lucif=gene.exp3.Z.MLL3.Lucif[highs,]
# genes2=genes1[highs,]
# ###SUCCESS!

mod = model.matrix(~ as.factor(tar.Z.MLL3.Lucif$KD)) ###do I need to use the internal scale as a factor
mod0 = model.matrix(~1, data=tar.Z.MLL3.Lucif)
svaobj = sva(gene.exp3.Z.MLL3.Lucif,mod,mod0)
save(svaobj, file = "ZR751.shLsh3.SVA")


cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

fcleanZR751shLsh3 = cleanY(gene.exp3.Z.MLL3.Lucif,mod,svaobj$sv)
save(fcleanZR751shLsh3, file = "stauffer.ZR751shLsh3.fclean")





#load("nonpara.batch.correct.arteaga.BRE.final.tss.expression")
mat=c()
bet = c()
pv  =c()
for (i in 1:dim(fcleanZR751shLsh3)[1]) {
  ref=as.character(rownames(fcleanZR751shLsh3)[i])
  print(ref)
  y=fcleanZR751shLsh3[i,]
  yv=as.vector(unlist(y))
  s=as.factor(tar.Z.MLL3.Lucif$KD)
  #s = tar1[,31]
  g=lm(yv~s)
  an=anova(g)
  summ=summary(g)
  v=c(ref,an[,5]) ###p-value from anova overall between the two categories
  beta = c(i,summ$coefficients[,1]) ###Estimate
  pval=c(i,summ$coefficients[,4]) ###p value for KD status specifically
  tot.pval = c(ref, pval)
  mat=rbind(mat,v)
  bet=rbind(bet, beta)
  pv=rbind(pv, tot.pval)
}

qpv_otherKD=qvalue(as.numeric(mat[,2])) ###this is the p-value for the ANOVA
summary(qpv_otherKD)
sigresisKD = which(qpv_otherKD$qvalue<0.05)  
length(sigresisKD) 
matKDSignificant = mat[sigresisKD,]
betKDSignificant = bet[sigresisKD,]
pvKDSignificant = pv[sigresisKD,]
genes1Sig = genes1[sigresisKD,]

##pvfclean[,4] is the column for p-value for MutStatus in DEGs
###betfclean[,3] is the column for estimate
###matfclean[,2] is the column for p-value for that category overall in assigning variance
refgene = genes1Sig[,10]
matKDSignificant.1 = matKDSignificant[,2] ###### THE Q-VALUES, AND THE COLUMN NUMBERS FOR THIS ANALYSIS, this should be correct
betKDSignificant.1 = betKDSignificant[,3]
pvKDSignificant.1 = pvKDSignificant[,4] 
MLL3KDANOVAtablefclean = c()
MLL3KDANOVAtablefclean= cbind(refgene, matKDSignificant.1, betKDSignificant.1, pvKDSignificant.1) 
MLL3KDANOVAtablefclean.1 = as.table(MLL3KDANOVAtablefclean)
length(which(MLL3KDANOVAtablefclean.1[,4]<0.05))
### going to keep only the genes that were differentially expressed based on MLL3 Mut Status
KeepforSigfclean = which(MLL3KDANOVAtablefclean.1[,4]<0.05)
MLL3KDANOVAtablefclean.2 = MLL3KDANOVAtablefclean.1[KeepforSigfclean,]
dim(MLL3KDANOVAtablefclean.1)
dim(MLL3KDANOVAtablefclean.2)
### now going to divide the genes into those that are upregulated vs downregulated
length(which(MLL3KDANOVAtablefclean.2[,3]<0)) ###693
length(which(MLL3KDANOVAtablefclean.2[,3]>0)) ###688 
MLL3KDUpfclean = which(MLL3KDANOVAtablefclean.2[,3]>0)
MLL3KDDownfclean = which(MLL3KDANOVAtablefclean.2[,3]<0)
length(MLL3KDDownfclean)
length(MLL3KDUpfclean)
MLL3KDANOVAtableUpfclean.2 = MLL3KDANOVAtablefclean.2[MLL3KDUpfclean,]
MLL3KDANOVAtableDownfclean.2 = MLL3KDANOVAtablefclean.2[MLL3KDDownfclean,]
write.table(MLL3KDANOVAtableUpfclean.2, file="sigUPMLL3KDfcleanJuly2018.txt", col.name = F, row.names = F, quote = F, sep = "\t")
write.table(MLL3KDANOVAtableDownfclean.2, file="sigDOWNMLL3KDfcleanJuly2018.txt", col.name = F, row.names = F, quote = F, sep = "\t")
write.table(MLL3KDANOVAtablefclean.1, file= "allgenesMLL3KDfcleanJuly2018.txt", col.name = F, row.names = F, quote = F, sep = "\t")




