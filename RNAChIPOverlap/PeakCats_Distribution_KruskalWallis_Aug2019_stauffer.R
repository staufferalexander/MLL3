setwd("~/Desktop/Kim/")
require(RColorBrewer)
require(pcaMethods)
require(qvalue)
require(limma)
require(marray)
require(sva)
require(Gviz)
require(genefu)
library(ggplot2)
library(scales)
require(gplots)
require(reshape2)
library(plyr)
library(dplyr)

V = read.table("/Users/staufferkm/Desktop/PeakAssocDEG_MEME_Sept2019/PeakLevelAnalysisTable_meltfixed_Aug2019.txt", header = T, quote = "", as.is = T, sep = "\t")
V$CatAm = NA
for (i in 1:nrow(V)) {
  V[i,5] = paste(V[i,2], V[i,3], sep = ".")
}

E.Loss.Up = which((V$Cat == "ERa loss") & (V$DOE =="Up"))
E.Loss.Down = which((V$Cat == "ERa loss") & (V$DOE =="Down"))
E.Gain.Up = which((V$Cat == "ERa gain") & (V$DOE =="Up"))
E.Gain.Down = which((V$Cat == "ERa gain") & (V$DOE =="Down"))
E.Main.Up = which((V$Cat == "ERa maintain") & (V$DOE =="Up"))
E.Main.Down = which((V$Cat == "ERa maintain") & (V$DOE =="Down"))

VE.Loss.Up = V[E.Loss.Up,]
VE.Loss.Down = V[E.Loss.Down,]
VE.Gain.Up = V[E.Gain.Up,]
VE.Gain.Down = V[E.Gain.Down,]
VE.Main.Up = V[E.Main.Up,]
VE.Main.Down = V[E.Main.Down,]

hist(VE.Loss.Up$Peak)
V$CatAm = as.factor(V$CatAm)

what = which(V$Peak < (-1000000))
WHAT = V[what,]

#V$CatAm <- ordered(V$CatAm, levels = c("ctrl", "trt1", "trt2"))

group_by(V, V$CatAm) %>%
  summarise(
    count = n(),
    mean = mean(Peak, na.rm = TRUE),
    sd = sd(Peak, na.rm = TRUE),
    median = median(Peak, na.rm = TRUE),
    IQR = IQR(Peak, na.rm = TRUE)
  )

# A tibble: 6 x 6
# `V$CatAm`         count    mean      sd median     IQR
# <fct>             <int>   <dbl>   <dbl>  <dbl>   <dbl>
#   1 Down.ERa gain      1360  24271. 442012.   9799 566575 
# 2 Down.ERa loss      3789 -13810. 433852.      0 483823 
# 3 Down.ERa maintain  1196  16616. 402121.      0 369807.
# 4 Up.ERa gain         439 -68331. 409784.      0 475504.
# 5 Up.ERa loss        1273 -12472. 418422.      0 456183 
# 6 Up.ERa maintain     449 -19116. 326865.      0 184897 
###Kruskal-Wallis test by rank is a non-parametric alternative to one-way ANOVA test, which extends the two-samples Wilcoxon test in the situation where there are more than two groups. It’s recommended when the assumptions of one-way ANOVA test are not met.
#kruskal.test(y~A)
kruskal.test(V$Peak ~ V$Cat, data = V)

# Kruskal-Wallis rank sum test
# 
# data:  V$Peak by V$CatAm
# Kruskal-Wallis chi-squared = 19.394, df = 5, p-value = 0.001623

##As the p-value is less than the significance level 0.05, we can conclude that there are significant differences between the treatment groups.
##From the output of the Kruskal-Wallis test, we know that there is a significant difference between groups, but we don’t know which pairs of groups are different.

pairwise.wilcox.test(V$Peak, V$Cat, p.adjust.method = "BH")


# 
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  V$Peak and V$CatAm 
# 
# Down.ERa gain Down.ERa loss Down.ERa maintain Up.ERa gain Up.ERa loss
# Down.ERa loss     0.0121        -             -                 -           -          
#   Down.ERa maintain 0.3124        0.1186        -                 -           -          
#   Up.ERa gain       0.0055        0.1256        0.0376            -           -          
#   Up.ERa loss       0.0613        0.7000        0.3124            0.1186      -          
#   Up.ERa maintain   0.0376        0.6608        0.1186            0.5163      0.5163     
# 
# P value adjustment method: BH 
VE.Loss = rbind(VE.Loss.Up, VE.Loss.Down)
VE.Gain = rbind(VE.Gain.Up, VE.Gain.Down)
VE.Maintain = rbind (VE.Main.Up, VE.Main.Down)

ks.test(VE.Loss$Peak, VE.Gain$Peak, p.adjust.method = "BH")
# Two-sample Kolmogorov-Smirnov test
# 
# data:  VE.Loss$Peak and VE.Gain$Peak
# D = 0.054309, p-value = 0.0007956
# alternative hypothesis: two-sided

ks.test(VE.Gain$Peak, VE.Maintain$Peak, p.adjust.method = "BH")
# Two-sample Kolmogorov-Smirnov test
# 
# data:  VE.Gain$Peak and VE.Maintain$Peak
# D = 0.13458, p-value = 6.062e-14
# alternative hypothesis: two-sided

ks.test(VE.Loss$Peak, VE.Maintain$Peak, p.adjust.method = "BH")
# Two-sample Kolmogorov-Smirnov test
# 
# data:  VE.Loss$Peak and VE.Maintain$Peak
# D = 0.091348, p-value = 2.006e-09
# alternative hypothesis: two-sided

ks.test(VE.Loss.Up$Peak, VE.Loss.Down$Peak, p.adjust.method = "BH")
# Two-sample Kolmogorov-Smirnov test
# 
# data:  VE.Loss.Up$Peak and VE.Loss.Down$Peak
# D = 0.025962, p-value = 0.5419
# alternative hypothesis: two-sided

ks.test(VE.Gain.Up$Peak, VE.Gain.Down$Peak, p.adjust.method = "BH")
# Two-sample Kolmogorov-Smirnov test
# 
# data:  VE.Gain.Up$Peak and VE.Gain.Down$Peak
# D = 0.1164, p-value = 0.0002486
# alternative hypothesis: two-sided

ks.test(VE.Main.Up$Peak, VE.Main.Down$Peak, p.adjust.method = "BH")
# Two-sample Kolmogorov-Smirnov test
# 
# data:  VE.Main.Up$Peak and VE.Main.Down$Peak
# D = 0.10491, p-value = 0.001515
# alternative hypothesis: two-sided


ks.test()


### Make some overlaid histograms http://www.cookbook-r.com/Graphs/Plotting_distributions_(ggplot2)/
#Down DEG - ERa Loss v ERa gain v ERa Maintain --> should have diff dist
#Up DEG - ERa Loss v ERa gain v ERa Maintain --> loss and gain should have similar dist
#ERa Loss - Up v Down --> should have similar dist
#ERa Gain - Up v Down --> should have similar dist
#ERa Maintain - Up v Down --> should have different dist

d = which(V$DOE == "Down")
D = V[d,]
u = which(V$DOE == "Up")
U = V[u,]
l = which(V$Cat == "ERa loss")
L = V[l,]
g = which(V$Cat == "ERa gain")
G = V[g,]
m = which(V$Cat =="ERa maintain")
M = V[m,]
#Down
ggplot(D, aes(x=D$Peak, fill=D$CatAm)) +
  geom_histogram(binwidth=80000, position="identity", alpha = 0.5) +
  scale_fill_manual(values=c("#42adf5", "#073b59", "#2145d9")) +
  #xlim(-1000000,1000000) +
  labs(title="ERa Peaks in ZR751 Associated with Downregulated Genes", x="Peak Distance to TSS of DEG", y = "Number of Peaks")
#Up
ggplot(U, aes(x=U$Peak, fill=U$CatAm)) +
  geom_histogram(binwidth=80000, position="identity", alpha = 0.7) +
  scale_fill_manual(values=c("#c27d06", "#c2bc11", "#ffe62b")) +
  #xlim(-2000000,1000000) +
  labs(title="ERa Peaks in ZR751 Associated with Upregulated Genes", x="Peak Distance to TSS of DEG", y = "Number of Peaks")
#Loss
ggplot(L, aes(x=L$Peak, fill=L$CatAm)) +
  geom_histogram(binwidth=80000, position="identity", alpha = 0.5) +
  scale_fill_manual(values=c("#0000FF", "#FFFF00")) +
  #xlim(-2000000,1000000) +
  labs(title="ERa Peaks in ZR751 Lost Upon MLL3 KD\nAssociated with DEG", x="Peak Distance to TSS of DEG", y = "Number of Peaks")
#Gain
ggplot(G, aes(x=G$Peak, fill=G$CatAm)) +
  geom_histogram(binwidth=80000, position="identity", alpha = 0.5) +
  scale_fill_manual(values=c("#0000FF", "#FFFF00")) +
  #xlim(-2000000,1000000) +
  labs(title="ERa Peaks in ZR751 Gained Upon MLL3 KD\nAssociated with DEG", x="Peak Distance to TSS of DEG", y = "Number of Peaks")
#Maintain
ggplot(M, aes(x=M$Peak, fill=M$CatAm)) +
  geom_histogram(binwidth=80000, position="identity", alpha = 0.5) +
  scale_fill_manual(values=c("#0000FF", "#FFFF00")) +
  #xlim(-2000000,1000000) +
  labs(title="ERa Peaks in ZR751 Maintained Upon MLL3 KD\nAssociated with DEG", x="Peak Distance to TSS of DEG", y = "Number of Peaks")
