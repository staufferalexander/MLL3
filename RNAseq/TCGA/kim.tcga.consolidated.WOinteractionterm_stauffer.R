setwd("~/Dropbox/Kim/cuffnorm/")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("genefilter")
require(RColorBrewer)
require(pcaMethods)
require(qvalue)
require(limma)
require(marray)
require(sva)
require(Gviz)
require(genefu)
require(genefilter)
tcga.clin = read.delim(file = "/Users/staufferkm/Desktop/Data/Sequencing/RNA-seq/TCGA Analysis/TCGA files/tcga.clinical.data.052414.txt", header = T, as.is = T, sep = "\t", quote = "")
genes = read.delim(file = "~/Dropbox/Kim/RNASeq_Dropbox/TCGA/stauffer.genes.entrez.attr_table", header = T, as.is = T, sep = "\t", quote = "")
gene.exp= read.table(file ="~/Dropbox/Kim/RNASeq_Dropbox/TCGA/TCGA Graphs and Files/genes.fpkm_table", header = T, as.is = T, sep = "\t", quote = "")
gene.exp1 =gene.exp[,-1]
tar = read.delim(file = "~/Dropbox/Kim/RNASeq_Dropbox/TCGA/TCGA Graphs and Files/test.tcga.target.file.mostrecent.txt", header = T, as.is = T, sep = "\t", quote = "")
missing = which(tar[,2] == "missing") 
gene.exp2=gene.exp1[,-missing]
tar1=tar[-missing,]
samples.table = read.delim(file = "~/Dropbox/Kim/RNASeq_Dropbox/TCGA/samples.table", header = T, as.is = T, sep = "\t", quote = "")
samples.table1=samples.table[-missing,]
low=which(samples.table1[,4]<.35)
gene.exp3=gene.exp2[,-low]
tar2=tar1[-low,]
samples.table2=samples.table1[-low,]
m=apply(gene.exp3, 1, mean)
high=which(log2(m+.5)>1.5)
genes1 = genes[high,] #so we are only taking the genes that are actually expressed here
gene.exp4 = gene.exp3[high,]
gene.exp5 = data.matrix(gene.exp4)
gene.exp6 = log2(gene.exp5+.5)
mut_miss = which(tar2$MLL3_Bin=="Na")
tar3=tar2[-mut_miss,]
gene.exp7 = gene.exp6[,-mut_miss]
males = which(tar3[,14] == "MALE") 
gene.exp8=gene.exp7[,-males]
tar4=tar3[-males,]
ERpos=which(tar4$ER.Status=="Positive")
head(ERpos)
tar4ER=tar4[ERpos,]
gene.exp8ER=gene.exp8[,ERpos]
dim(tar4ER) #make sure the dimensions match for all these
dim(gene.exp8ER)
LumA=which(tar4ER$Call=="LumA")
LumB=which(tar4ER$Call=="LumB")
Lum_all=union(LumA,LumB)
tar4Lum=tar4ER[Lum_all,]  
gene.exp8Lum=gene.exp8ER[,Lum_all] 
dim(tar4Lum)
dim(gene.exp8Lum) 

mod = model.matrix(~ as.factor(tar4Lum$MLL3_Bin) + as.factor(tar4Lum$Call) + as.factor(tar4Lum$hist_type_general))
mod0 = model.matrix(~1, data=tar4Lum) 
svaobj = sva(gene.exp8Lum,mod,mod0)
save(svaobj, file = "stauffer.ERLum.SVA")
#svaobj1 = sva(iso.exp8ER,mod,mod0)

cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

fcleanERLumMLL3 = cleanY(gene.exp8Lum,mod,svaobj$sv)
save(fcleanERLumMLL3, file = "stauffer.ERLumMLL3.fclean")

matfcleanERLumMLL3=c()
betfcleanERLumMLL3 = c()
residualsfcleanERLumMLL3=c()
pvfcleanERLumMLL3 =c()
for (i in 1:dim(fcleanERLumMLL3)[1]) {	##for each sample
  ref=as.character(genes1[i,4])
  print(ref)
  y=gene.exp8Lum[i,]
  yv=as.vector(unlist(y))
  ### CHECK ALL OF THESE COLUMNS
  ### mll3 mutation status, luminal status, ER status, menopause_status, number_of_lymphnodes_positive_by_he, proliferation_score, age_at_initial_pathologic_diagnosis
  ### could also add in the tar3ER$expmutstatus or create an interaction term between mut status and expression level
  ### figure out how to add those variables
  s=as.factor(tar4Lum$PIK3CA_bin) #PIK3CA mutation status
  w=as.factor(tar4Lum$MLL3_Bin) # MLL3 mutation status
  t=as.factor(tar4Lum$hist_type_general) #histological subtype (Ductal, Lobular, or Mixed)
  u=as.factor(tar4Lum$race) #race . . . to do SVA with this need to make all the NA samples either white or take them out
  v=as.factor(tar4Lum$age_at_initial_pathologic_diagnosis) # age
  x=as.factor(tar4Lum$Call) # Call
  z=as.factor(tar4Lum$ASCOM_bin) ##any member of ASCOM mutated (MLL4 not included)
  g=lm(yv~w+x+t) # also do this without the mll3 and can use that to check the residuals 
  anfcleanERLumMLL3=anova(g)
  summfcleanERLumMLL3=summary(g)
  residfcleanERLumMLL3 = resid(g)
  v=c(ref,anfcleanERLumMLL3[,5])
  #beta = c(ref,an[,4])
  residuals=c(i,residfcleanERLumMLL3)
  beta=c(i,summfcleanERLumMLL3$coefficients[,1])
  pval=c(i,summfcleanERLumMLL3$coefficients[,4])
  #tot.beta = c(ref2, beta)
  tot.pval = c(ref, pval)
  #residualsfcleanERLumMLL3=rbind(residualsfcleanERLumMLL3,residuals)
  matfcleanERLumMLL3=rbind(matfcleanERLumMLL3,v) ###this is the p-val associated with the F-value for each variable (MutStatus,Hist,Call) in my model
  betfcleanERLumMLL3=rbind(betfcleanERLumMLL3, beta) ###this is the estimate (so what I want to pull for my up or downregulations)
  pvfcleanERLumMLL3=rbind(pvfcleanERLumMLL3, tot.pval) ### this is the p-value for each variable compared against all other categories in that variable
  ##playing around here with capturing my output
  capture.output(anfcleanERLumMLL3, file = "anovafcleanERLumMLL3.32318.txt", append = TRUE)
  capture.output(summfcleanERLumMLL3, file = "anovasummaryfcleanERLumMLL3.32318.txt", append = TRUE)
  #capture.output(an, file="rms.v24.final.pvalue.isoform.anova", append=TRUE)
  #capture.output(an, file="rms.v24.final.pvalue.isoform.anova", append=TRUE)
  #capture.output(name, file="rms.v24.final.beta.isoform", append=TRUE)
  #capture.output(summ, file="rms.v24.final.beta.isoform", append=TRUE)
}


#residualsfcleanERLumMLL3.2 = cbind(refgene, residualsfcleanERLumMLL3)
#write.table(residualsfcleanERLumMLL3.2, file="residualsfcleanERLumMLL3_8182016.txt", col.name = F, row.names = F, quote = F, sep = "\t")



genetable = read.delim(file = "allgenesMLL3MutStatusfcleanERLumMLL32216.txt", header = T, as.is = T, sep = "\t", quote = "")
qpv_otherMLL3=qvalue(as.numeric(genetable[,2]))
summary(qpv_otherMLL3)
sigresisMLL3 = which(qpv_otherMLL3$qvalue<0.05)  
length(sigresisMLL3) 
p05genetable = genetable[sigresisMLL3,]
write.table(p05genetable, file= "q05MLL3MutStatusfcleanERLumMLL32616.txt", col.name = F, row.names = F, quote = F, sep = "\t")
sigresisMLL301 = which(qpv_otherMLL3$qvalue<0.01)  
length(sigresisMLL301) 
p01genetable = genetable[sigresisMLL301,]
write.table(p01genetable, file= "q01MLL3MutStatusfcleanERLumMLL32616.txt", col.name = F, row.names = F, quote = F, sep = "\t")




qpv_otherMLL3=qvalue(as.numeric(matfcleanERLumMLL3[,2]))
qpv_otherCall=qvalue(as.numeric(matfcleanERLumMLL3[,3]))
qpv_otherHist=qvalue(as.numeric(matfcleanERLumMLL3[,4]))
summary(qpv_otherMLL3)
summary(qpv_otherCall)
summary(qpv_otherHist)
sigresisMLL3 = which(qpv_otherMLL3$qvalue<0.05)  
length(sigresisMLL3) 
sigresisCall = which(qpv_otherCall$qvalue<0.05)  
length(sigresisCall) 
sigresisHist = which(qpv_otherHist$qvalue<0.05)  
length(sigresisHist) 
## make sure that we are not overfitting, run without the interaction term
### check q values on cross, intrinsic subtype, histologic subtype

matfcleanERLumMLL3Significant = matfcleanERLumMLL3[sigresisMLL3,]
betfcleanERLumMLL3Significant = betfcleanERLumMLL3[sigresisMLL3,]
pvfcleanERLumMLL3Significant = pvfcleanERLumMLL3[sigresisMLL3,]

##pvfclean[,4] is the column for p-value for MutStatus in DEGs
###betfclean[,3] is the column for estimate
###matfclean[,2] is the column for p-value for that category overall in assigning variance
#MutStatusSig = which(matfERLumASCOMPIKX[,3]<0.05)
##length(MutStatusSig) ###1262
#MutStatusSigPV = which(pvfERLumASCOMPIKX[,5]<0.05)
#length(MutStatusSigPV) ###1163
###use genes1 for your ref as first column
##refgene= genes1[,4] BRING THIS BACK UP IN A MINUTE!
refgene = genes1[,9]
MatfcleanERLumMLL3.1 = matfcleanERLumMLL3[,2] ###### NEED TO CHECK ALL THIS STUFF. THE Q-VALUES, AND THE COLUMN NUMBERS FOR THIS ANALYSIS, this should be correct 2-2-16
BetfcleanERLumMLL3.1 = betfcleanERLumMLL3[,3]
PVfcleanERLumMLL3.1 = pvfcleanERLumMLL3[,4] 
MLL3MutStatusANOVAtablefcleanERLumMLL3 = c()
MLL3MutStatusANOVAtablefcleanERLumMLL3 = cbind(refgene, MatfcleanERLumMLL3.1, BetfcleanERLumMLL3.1, PVfcleanERLumMLL3.1) 
MLL3MutStatusANOVAtablefcleanERLumMLL3.1 = as.table(MLL3MutStatusANOVAtablefcleanERLumMLL3)
#length(which(MLL3MutStatusANOVAtablefcleanERLumMLL3.1[,2]<0.05))
### going to keep only the genes that were differentially expressed based on MLL3 Mut Status
#KeepforSigfcleanERLumMLL3 = which(MLL3MutStatusANOVAtablefcleanERLumMLL3.1[,2]<0.05)
MLL3MutStatusANOVAtablefcleanERLumMLL3.2 = MLL3MutStatusANOVAtablefcleanERLumMLL3.1[sigresisMLL3,]
dim(MLL3MutStatusANOVAtablefcleanERLumMLL3.1)
dim(MLL3MutStatusANOVAtablefcleanERLumMLL3.2)
### now going to divide the genes into those that are upregulated vs downregulated
length(which(MLL3MutStatusANOVAtablefcleanERLumMLL3.2[,3]<0)) ###693
length(which(MLL3MutStatusANOVAtablefcleanERLumMLL3.2[,3]>0)) ###688 
MLL3MutUpfcleanERLumMLL3 = which(MLL3MutStatusANOVAtablefcleanERLumMLL3.2[,3]>0)
MLL3MutDownfcleanERLumMLL3 = which(MLL3MutStatusANOVAtablefcleanERLumMLL3.2[,3]<0)
length(MLL3MutDownfcleanERLumMLL3)
length(MLL3MutUpfcleanERLumMLL3)
MLL3MutStatusANOVAtableUpfcleanERLumMLL3.2 = MLL3MutStatusANOVAtablefcleanERLumMLL3.2[MLL3MutUpfcleanERLumMLL3,]
MLL3MutStatusANOVAtableDownfcleanERLumMLL3.2 = MLL3MutStatusANOVAtablefcleanERLumMLL3.2[MLL3MutDownfcleanERLumMLL3,]
write.table(MLL3MutStatusANOVAtableUpfcleanERLumMLL3.2, file="sigUPMLL3MutStatusfcleanERLumMLL32216.txt", col.name = F, row.names = F, quote = F, sep = "\t")
write.table(MLL3MutStatusANOVAtableDownfcleanERLumMLL3.2, file="sigDOWNMLL3MutStatusfcleanERLumMLL32216.txt", col.name = F, row.names = F, quote = F, sep = "\t")
write.table(MLL3MutStatusANOVAtablefcleanERLumMLL3.1, file= "allgenesMLL3MutStatusfcleanERLumMLL32216.txt", col.name = F, row.names = F, quote = F, sep = "\t")

