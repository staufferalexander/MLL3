##### GRMetrics
# inputData - the name of an input data frame with the following columns as well as other grouping columns
# concentration - column with concentration values (not log transformed) of the perturbagen on which dose-response curves will be evaluated
# cell_count - column with the measure of cell number (or a surrogate of cell number) after treatment
# cell_count__time0 - column with initial (Time 0) cell counts - the measure of cell number in untreated wells grown in parallel until the time of treatment
# cell_count__ctrl - column with the Control cell count: the measure of cell number in control (e.g. untreated or DMSO-treated) wells from the same plate
# All other columns will be treated as additional keys on which the data will be grouped (e.g. cell_line, drug, time, replicate)

# install GRMetrics
#source("http://bioconductor.org/biocLite.R")
#biocLite("GRmetrics")
#install.packages("devtools")
#devtools::install_github('hadley/ggplot2')
#install.packages("GRmetrics")
library("GRmetrics")

##load data
input = read.delim("/path/to/GRmetrics/data", header = T, sep = "\t", quote = "",  as.is = T)
input1 = as.data.frame(input)
which(sapply(input1, function(x) { length(unique(x)) }) == 1)
which(sapply(input1, function(x) (is.character(x) | is.factor(x)) & length(unique(x))<2))
length(unique(cell_line))
length(unique(time))

 ###################################################################
 #Calculate GR values and solve for GR metrics parameters (i.e. fit curves)
drc_output = GRmetrics::GRfit(input1, groupingVariables = c('cell_line', 'time'))

#See overview of output data (SummarizedExperiment object)
drc_output

#Review output table of GR metrics parameters
head(GRmetrics::GRgetMetrics(drc_output))

#View descriptions of each GR metric (or goodness of fit measure)
View(GRmetrics::RgetDefs(drc_output))

#Review output table of GR values
head(GRmetrics::GRgetValues(drc_output))

#View grouping variables used for calculation
GRmetrics::GRgetGroupVars(drc_output)

# export your results
# Write GR metrics parameter table to tab-separated text file
write.table(GRmetrics::GRgetMetrics(drc_output), file = "/path/to/your/output/", quote = FALSE, sep = "\t", row.names = FALSE)
# Write original data plus GR values to comma-separated file
write.table(GRmetrics::GRgetValues(drc_output), file = "/path/to/you/output/", quote = FALSE, sep = "\t", row.names = FALSE)

# draw GR dose-response curves with plotly or with ggplot2. You can also specify the range of the graph.
# Draw dose-response curves
GRmetrics::GRdrawDRC(drc_output, metric = "GR", experiments = "all", min = "auto", max = "auto", points = TRUE, curves = TRUE, plotly = FALSE)
day4_Tam_dose_response = GRdrawDRC(drc_output, metric = "GR", experiments = c('ZR751shMLL2 Tam 4 days', 'ZR751shMLL3 Tam 4 days','ZR751shKDM6A Tam 4 days', 'ZR751shLucif Tam 4 days'), min = "auto", max = "auto", points = TRUE, curves = TRUE, plotly = FALSE)
day8_Tam_dose_response =  GRmetrics::GRdrawDRC(drc_output, metric = "GR", experiments = c('ZR751shMLL3 Tamoxifen 8', 'ZR751shLucif Tamoxifen 8'), min = "auto", max = "auto", points = TRUE, curves = TRUE, plotly = FALSE)

#GRmetrics::GRdrawDRC(day8_Fulv_dose_response, metric = "GR", experiments = "all", min = "auto", max = "auto", points = TRUE, curves = TRUE, plotly = FALSE)

# draw scatterplots and boxplots of GR metrics with plotly or ggplot2
output1 = GRfit(input1, groupingVariables = c('cell_line', 'time'), case = "A")

## Draw scatterplots
#scatter_drugs = GRscatter(output1, 'GR50', 'agent', c('Fulv','Tam'), 'BYL')
#scatter_drugs_noplotly = GRscatter(output1, 'GR50', 'agent', c('Fulv','Tam'), 'BYL', plotly = FALSE)

# Draw boxplots
boxplot_MLL3Lucif_noplotly = GRbox(output1, metric ='GRinf', groupVariable = c('cell_line'), pointColor = 'cell_line', factors = c('ZR751shMLL3', 'ZR751shLucif'), plotly = FALSE)

############## GR Metrics that are calculated
# GR50 -> the concentration at which the effect reaches a GR value of 0.5 based on interpolation of the fitted curve.
# GRmax -> the effect at the highest tested concentration. Note that GRmax can differ from GRinf if the dose-response does not reach its plateau value.
# GR_AOC -> the area over the dose-response curve, which is the integral of 1-GR(c) over the range of concentrations tested, normalized by the range of concentration.
# GEC50 -> the drug concentration at half-maximal effect, which reflects the potency of the drug.
# GRinf -> GR(c->inf), which reflects asymptotic drug efficacy.
# h_GR -> The Hill coefficient of the fitted (GR) curve, which reflects how steep the dose response curve is
# r2_GR -> The coefficient of determination - essentially how well the (GR) curve fits to the data points
# pval_GR -> The p-value of the F-test comparing the fit of the (GR) curve to a horizontal line fit
# flat_fit_GR -> For data that doesn’t significantly fit better to a curve than a horizontal line fit (p > 0.05), the y value (GR) of the flat line
# IC50 -> The concentration at which relative cell count = 0.5
# Emax -> The maximal effect of the drug (minimal relative cell count value)
# AUC -> The ‘Area Under the Curve’ - The area below the fitted (traditional) dose response curve
# EC50 -> The concentration at half-maximal effect (not growth rate normalized)
# Einf -> The asymptotic effect of the drug (not growth rate normalized)
# h -> The Hill coefficient of the fitted (traditional) dose response curve, which - reflects how steep the dose response curve is
# r2_IC -> The coefficient of determination - essentially how well the (traditional) curve fits to the data points
# pval_IC -> The p-value of the F-test comparing the fit of the (traditional) curve to a horizontal line fit
# flat_fit_IC -> For data that doesn’t significantly fit better to a curve than a horizontal line fit (p > 0.05), the y value (relative cell count) of the flat line
# In addition to the metrics, the scripts report the r-squared of the fit and evaluate the significance of the sigmoidal fit based on an F-test. If the fit is not significant (p > 0.05), the sigmoidal fit is replaced by a constant value (flat fit). This can be circumvented by using the “force” option in the GRfit function.

