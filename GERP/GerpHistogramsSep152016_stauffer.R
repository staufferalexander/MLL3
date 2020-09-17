library(plyr)

# PTEN
# Mean	Number
# ERNeg	5.19	3
# ERPos	4.578125	16
# Model	4.715384615	13
# All	9.779444444	19
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/PTENERneg_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/PTENERpos_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandomGERP = read.delim(file = "~/Desktop/Stricker Lab/GERP/GERP/GERPfiles/GERPlatest/PTENModel_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/PTENAll_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandGERP = as.data.frame(RandomGERP)
#MutMean = 5.19	
#MutMean = 4.578125	
MutMean = 4.715384615	
#MutMean = 9.779444444	
minX = (min(RandGERP) - 1)
if (max(RandGERP) > MutMean) {
  maxX = (max(RandGERP) + 1)
} else {
  maxX = (MutMean + 1)
}
Freq = count(RandGERP, 'V1')
HG = hist(RandGERP[,1], breaks = 100000)
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "PTEN BrCa ER- Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "PTEN BrCa ER+ Gerp-Set Means")
HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "PTEN BrCa Model-Samples Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "PTEN BrCa All Samples Gerp-Set Means")
MeanGerp = abline(v= MutMean, col="red", lwd = 3, lty= 22)
MeanGerp
Larger = which(RandGERP[,1] > MutMean)
LengthLarger = length(Larger)
Significance = LengthLarger/10000
Significance

# PIK3CA
# Mean	Number
# ERNeg	5.65952381	21
# ERPos	5.668267857	168
# Model	5.707795276	127
# All	5.667296296	189
# 
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/PIK3CAERneg_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/PIK3CAERpos_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandomGERP = read.delim(file = "~/Desktop/Stricker Lab/GERP/GERP/GERPfiles/GERPlatest/PIK3CAModel_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/PIK3CAERAll_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandGERP = as.data.frame(RandomGERP)
#MutMean =5.65952381	
#MutMean =	5.668267857
MutMean =	5.707795276
#MutMean =	5.667296296
minX = (min(RandGERP) - 1)
if (max(RandGERP) > MutMean) {
  maxX = (max(RandGERP) + 1)
} else {
  maxX = (MutMean + 1)
}
Freq = count(RandGERP, 'V1')
HG = hist(RandGERP[,1], breaks = 100000)
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "PIK3CA BrCa ER- Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "PIK3CA BrCa ER+ Gerp-Set Means")
HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "PIK3CA BrCa Model-Samples Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "PIK3CA BrCa All Samples Gerp-Set Means")
MeanGerp = abline(v= MutMean, col="red", lwd = 3, lty= 22)
MeanGerp
Larger = which(RandGERP[,1] > MutMean)
LengthLarger = length(Larger)
Significance = LengthLarger/10000
Significance


# TTN
# Mean	Number
# ERNeg	2.867738235	68
# ERPos	3.408059057	227
# Model	3.40649028	164
# All	3.283510529	295
# 
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/TTNERneg_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/TTNERpos_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandomGERP = read.delim(file = "~/Desktop/Stricker Lab/GERP/GERP/GERPfiles/GERPlatest/TTNModel_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/TTNERpos_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandGERP = as.data.frame(RandomGERP)
#MutMean = 2.867738235
#MutMean = 	3.408059057
MutMean = 	3.40649028
#MutMean = 3.283510529	
minX = (min(RandGERP) - 1)
if (max(RandGERP) > MutMean) {
  maxX = (max(RandGERP) + 1)
} else {
  maxX = (MutMean + 1)
}
Freq = count(RandGERP, 'V1')
HG = hist(RandGERP[,1], breaks = 100000)
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "TTN BrCa ER- Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "TTN BrCa ER+ Gerp-Set Means")
HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "TTN BrCa Model-Samples Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "TTN BrCa All Samples Gerp-Set Means")
MeanGerp = abline(v= MutMean, col="red", lwd = 3, lty= 22)
MeanGerp
Larger = which(RandGERP[,1] > MutMean)
LengthLarger = length(Larger)
Significance = LengthLarger/10000
Significance


# TP53
# Mean	Number
# ERNeg	4.095168421	95
# ERPos	4.101234043	94
# Model	4.165342105	76
# All	4.098185185	189
# 
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/TP53ERneg_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/TP53ERpos_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/TP53Model_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandomGERP = read.delim(file = "~/Desktop/GERPfiles/TP53All_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandGERP = as.data.frame(RandomGERP)
#MutMean =	4.095168421
#MutMean =	4.101234043
#MutMean =	4.165342105
MutMean =	4.098185185
minX = (min(RandGERP) - 1)
if (max(RandGERP) > MutMean) {
  maxX = (max(RandGERP) + 1)
} else {
  maxX = (MutMean + 1)
}
Freq = count(RandGERP, 'V1')
HG = hist(RandGERP[,1], breaks = 100000)
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "TP53 BrCa ER- Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "TP53 BrCa ER+ Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "TP53 BrCa Model-Samples Gerp-Set Means")
HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "TP53 BrCa All Samples Gerp-Set Means")
MeanGerp = abline(v= MutMean, col="red", lwd = 3, lty= 22)
MeanGerp
Larger = which(RandGERP[,1] > MutMean)
LengthLarger = length(Larger)
Significance = LengthLarger/10000
Significance

# KMT2D
# Mean	Number
# ERNeg	3.884555556	9
# ERPos	2.58	17
# Model	2.206363636	11
# All	3.031576923	26
# 
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2DERneg_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2DERpos_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2DModel_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2DAll_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandGERP = as.data.frame(RandomGERP)
#MutMean =3.884555556
#MutMean =2.58
#MutMean =2.206363636
MutMean =3.031576923
minX = (min(RandGERP) - 1)
if (max(RandGERP) > MutMean) {
  maxX = (max(RandGERP) + 1)
} else {
  maxX = (MutMean + 1)
}
Freq = count(RandGERP, 'V1')
HG = hist(RandGERP[,1], breaks = 100000)
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2D BrCa ER- Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2D BrCa ER+ Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2D BrCa Model-Samples Gerp-Set Means")
HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2D BrCa All Samples Gerp-Set Means")
MeanGerp = abline(v= MutMean, col="red", lwd = 3, lty= 22)
MeanGerp
Larger = which(RandGERP[,1] > MutMean)
LengthLarger = length(Larger)
Significance = LengthLarger/10000
Significance


# KDM6A
# Mean	Number
# ERNeg	5.48	2
# ERPos	5.14	5
# Model	5.03	4
# All	5.237142857	7
# 
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KDM6AERneg_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KDM6AERpos_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KDM6AModel_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KDM6AAll_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandGERP = as.data.frame(RandomGERP)
#MutMean =	5.48
#MutMean =	5.14
#MutMean =	5.03
MutMean =	5.237142857
minX = (min(RandGERP) - 1)
if (max(RandGERP) > MutMean) {
  maxX = (max(RandGERP) + 1)
} else {
  maxX = (MutMean + 1)
}
Freq = count(RandGERP, 'V1')
HG = hist(RandGERP[,1], breaks = 100000)
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KDM6A BrCa ER- Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KDM6A BrCa ER+ Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KDM6A BrCa Model-Samples Gerp-Set Means")
HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KDM6A BrCa All Samples Gerp-Set Means")
MeanGerp = abline(v= MutMean, col="red", lwd = 3, lty= 22)
MeanGerp
Larger = which(RandGERP[,1] > MutMean)
LengthLarger = length(Larger)
Significance = LengthLarger/10000
Significance


# KMT2A
# Mean	Number
# ERNeg	5.52	2
# ERPos	2.231	20
# Model	2.74	10
# All	2.53	22
#
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2AERneg_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2AERpos_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2AModel_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2AAll_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandGERP = as.data.frame(RandomGERP)
#MutMean =5.52	
#MutMean =2.231
#MutMean =2.74
MutMean =2.53
minX = (min(RandGERP) - 1)
if (max(RandGERP) > MutMean) {
  maxX = (max(RandGERP) + 1)
} else {
  maxX = (MutMean + 1)
}
Freq = count(RandGERP, 'V1')
HG = hist(RandGERP[,1], breaks = 100000)
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2A BrCa ER- Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2A BrCa ER+ Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2A BrCa Model-Samples Gerp-Set Means")
HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2A BrCa All Samples Gerp-Set Means")
MeanGerp = abline(v= MutMean, col="red", lwd = 3, lty= 22)
MeanGerp
Larger = which(RandGERP[,1] > MutMean)
LengthLarger = length(Larger)
Significance = LengthLarger/10000
Significance


# KMT2B
# Mean	Number
# ERNeg	5.453333333	3
# ERPos	2.103	7
# Model	2.3842	5
# All	3.1081	10
# 
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2BERneg_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2BERpos_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2BModel_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2BAll_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandGERP = as.data.frame(RandomGERP)
#MutMean =	5.453333333
#MutMean =	2.103
#MutMean =	2.3842
MutMean =	3.1081
minX = (min(RandGERP) - 1)
if (max(RandGERP) > MutMean) {
  maxX = (max(RandGERP) + 1)
} else {
  maxX = (MutMean + 1)
}
Freq = count(RandGERP, 'V1')
HG = hist(RandGERP[,1], breaks = 100000)
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2B BrCa ER- Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2B BrCa ER+ Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2B BrCa Model-Samples Gerp-Set Means")
HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2B BrCa All Samples Gerp-Set Means")
MeanGerp = abline(v= MutMean, col="red", lwd = 3, lty= 22)
MeanGerp
Larger = which(RandGERP[,1] > MutMean)
LengthLarger = length(Larger)
Significance = LengthLarger/10000
Significance

# All = ERPos + ERNeg, does not include the "ER status Not available" samples
# KMT2C
# Mean	Number Muts
# ERNeg	4.645888889	9
# ERPos	3.474159091	44
# Model	4.612592593	27
# All	3.673132075	53
#
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2CERneg_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2CERpos_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandomGERP = read.delim(file = "~/Desktop/Stricker Lab/GERP/GERP/GERPfiles/GERPlatest/KMT2CModel_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
#RandomGERP = read.delim(file = "~/Desktop/GERPfiles/KMT2CAll_RandomGerpMeans.txt", header = F, as.is = T, sep = "\t", quote = "")
RandGERP = as.data.frame(RandomGERP)
#MutMean =	4.645888889	
#MutMean =	3.474159091
MutMean =	4.612592593
#MutMean =	3.673132075
minX = (min(RandGERP) - 1)
if (max(RandGERP) > MutMean) {
  maxX = (max(RandGERP) + 1)
} else {
  maxX = (MutMean + 1)
}
Freq = count(RandGERP, 'V1')
HG = hist(RandGERP[,1], breaks = 100000)
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2C BrCa ER- Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2C BrCa ER+ Gerp-Set Means")
HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "MLL3 BrCa Model-Samples Gerp-Set Means")
#HG = hist(RandGERP[,1],breaks=10000, xlim = c(minX, maxX), ylim = c(0,(max(HG$counts)+5)), xlab = "GERP Scores", main = "KMT2C BrCa All Samples Gerp-Set Means")
MeanGerp = abline(v= MutMean, col="red", lwd = 3, lty= 22)
MeanGerp
Larger = which(RandGERP[,1] > MutMean)
LengthLarger = length(Larger)
Significance = LengthLarger/10000
Significance
