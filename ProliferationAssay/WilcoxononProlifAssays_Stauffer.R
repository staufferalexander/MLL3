#Finding significance between shLucif and shMLL3 growth rate in response to Tam Fulv BYL
setwd("~/Desktop/")
require(RColorBrewer)
require(pcaMethods)
require(qvalue)
require(limma)
require(marray)
require(sva)
require(Gviz)
require(genefu)
library(ggplot2)
require(gplots)
library(gplots)
require(reshape2)
library(reshape2)
library(plyr)
library(dplyr)
#Cats will be the individual data points
Cats = read.table("/path.to/your/data/", header = T, quote = "", as.is = T, sep = "\t")

Cats$replicate = as.character(Cats$replicate) 
Cats$time = as.character(Cats$time)
Cats$concentration = as.character(Cats$concentration)


Fulv = which((Cats$agent == "Fulv"))
Cats.Fulv = Cats[Fulv,]
Tam = which((Cats$agent =="Tam"))
Cats.Tam = Cats[Tam,]

Dose.F.0 = which((Cats.Fulv$concentration == "0"))
Dose.F.1 = which((Cats.Fulv$concentration == "1"))
Dose.T.0 = which((Cats.Tam$concentration == "0"))
Dose.T.1 = which((Cats.Tam$concentration == "1"))
Cats.Fulv.0 = Cats.Fulv[Dose.F.0,]
Cats.Fulv.1 = Cats.Fulv[Dose.F.1,]
Cats.Tam.0 = Cats.Tam[Dose.T.0,]
Cats.Tam.1 = Cats.Tam[Dose.T.0,]

day8 = which((Cats.Fulv$time == "8"))
Cats.Fulv.8 = Cats.Fulv[day8,]
day10 = which((Cats.Fulv$time == "10"))
Cats.Fulv.10 = Cats.Fulv[day10,]
day14 = which((Cats.Fulv$time == "14"))
Cats.Fulv.14 = Cats.Fulv[day14,]
day16 = which((Cats.Fulv$time == "16"))
Cats.Fulv.16 = Cats.Fulv[day16,]


Dose.8.0.1 =  which((Cats.Fulv.8$concentration == "0.1"))
Dose.10.0.1 =  which((Cats.Fulv.10$concentration == "0.1"))
Dose.14.0.1 =  which((Cats.Fulv.14$concentration == "0.1"))
Dose.16.0.1 =  which((Cats.Fulv.16$concentration == "0.1"))
Dose.14.1 =  which((Cats.Fulv.14$concentration == "1"))
Dose.16.1 =  which((Cats.Fulv.16$concentration == "1"))
Dose.8.0 =  which((Cats.Fulv.8$concentration == "0"))
Dose.10.0 =  which((Cats.Fulv.10$concentration == "0"))
Dose.14.0 =  which((Cats.Fulv.14$concentration == "0"))
Dose.16.0 =  which((Cats.Fulv.16$concentration == "0"))
Cats.Fulv.8.0.1 = Cats.Fulv.8[Dose.8.0.1,]
Cats.Fulv.10.0.1 = Cats.Fulv.10[Dose.10.0.1,]
Cats.Fulv.14.0.1 = Cats.Fulv.14[Dose.14.0.1,]
Cats.Fulv.16.0.1 = Cats.Fulv.16[Dose.16.0.1,]
Cats.Fulv.14.1 = Cats.Fulv.14[Dose.14.1,]
Cats.Fulv.16.1 = Cats.Fulv.16[Dose.16.1,]
Cats.Fulv.8.0 = Cats.Fulv.8[Dose.8.0,]
Cats.Fulv.10.0 = Cats.Fulv.10[Dose.10.0,]
Cats.Fulv.14.0 = Cats.Fulv.14[Dose.14.0,]
Cats.Fulv.16.0 = Cats.Fulv.16[Dose.16.0,]

###	Wilcoxon rank sum test

wilcox.test(Cats.Tam.0$GRvalue ~ Cats.Tam.0$cell_line, data = Cats.Tam.0,
            alternative = c("two.sided"),
            mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95)
