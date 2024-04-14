######## meta analysis code using Hedge's D transformation

rm(list=ls())

# loading packages
library(metafor)

# importing data
meta_data<-read.csv("hedges_sheet.csv")
# effect sizes 1 = urban habitats and 2 = rural habitats


#############
# Identify significant means
head(meta_data)
sigmean <- which(meta_data$p.value < 0.05)
par(mfrow = c(1, 2))

# Plot 1: Hedge's  D vs. Sample size
plot(meta_data$d, meta_data$n, xlab = "Hedge's D", ylab = "Sample size")
points(meta_data$d[sigmean], meta_data$n[sigmean], pch = 16, col = "purple")
abline(v = 0.54, lty = 2)

# Plot 2: Hedge's D vs. Precision (vd)
plot(meta_data$d, 1/meta_data$vd, xlab = "Hedge's D", ylab = "Precision (1/vd)")
points(meta_data$d[sigmean], 1/meta_data$vd[sigmean], pch = 16, col = "purple")
abline(v = 0.54, lty = 2)
#############

#trying the above but removing outliers from Li et al. 
# Exclude rows 6-8 -> massive outliers for sample size
meta_data_subset <- meta_data[-(6:8), ]
plot(meta_data_subset$d, meta_data_subset$n, xlab = "Hedge's D", ylab = "Sample size")
points(meta_data_subset$d[sigmean], meta_data_subset$n[sigmean], pch = 16, col = "purple")
abline(v = 0.54, lty = 2)

# Plot 2: Hedge's D vs. Precision (vd)
plot(meta_data_subset$d, 1/meta_data_subset$vd, xlab = "Hedge's D", ylab = "Precision (1/vd)")
points(meta_data_subset$d[sigmean], 1/meta_data_subset$vd[sigmean], pch = 16, col = "purple")
abline(v = 0.54, lty = 2)

########################

# running a basic random effects meta-analysis (no random effects or moderators)
par(mfrow = c(1, 1))
model <- rma(d, vd, data=meta_data)
summary(model)
# creating a forest plot to visualise the results
forest(model, 
       cex.lab = 0.8,
       cex.axis = 0.8,
       addfit = TRUE,
       shade = "zebra",
       slab = paste(meta_data$reference, sep = ", "),
       cex = 0.8,
       order = "obs",
       xlab = "Effect Size")

# funnel plot for above model
funnel(model, col="#9933ff", pch = 16, xlab = "Effect Size", ylab = "vd")

?influence
###########################

# running a meta-analysis with random effects 
## a meta-analysis with random terms 
meta_data$id <- 1:18
#using taxon as a random effect
meta2<-rma.mv(yi=d,V=vd,random=~1|taxon/id,data = meta_data)
meta2

forest(meta2, 
       cex.lab = 0.8,
       cex.axis = 0.8,
       addfit = TRUE,
       shade = "zebra",
       slab = paste(meta_data$reference, sep = ", "),
       cex = 0.8,
       order = "obs",
       xlab = "Effect Size (d)")
summary(meta2)

funnel(meta2,  col="#9933ff", pch = 16, xlab = "Effect Size (d)", ylab = "Sampling Variance (vd)")

#################################

## meta-analysis with study as the random effect
# don't want to include both at once at risk of model-overfitting
meta3<-rma.mv(yi=d,V=vd,random=~1|reference/id,data = meta_data)
meta3

forest(meta3, 
       cex.lab = 0.8,
       cex.axis = 0.8,
       addfit = TRUE,
       shade = "zebra",
       slab = paste(meta_data$reference, sep = ", "),
       cex = 0.8,
       order = "obs",
       xlab = "Effect Size")
summary(meta3)
funnel(meta3,  col="#9933ff", pch = 16, xlab = "Effect Size", ylab = "vd" )
# May be more biologically meaningful to include taxon as the random effect. 

#testing for asymmetry (i.e. bias in the literature)
regtest(model)

#running a weighted analysis to account for publication bias (based on p-values)
install.packages("weightr")
library(weightr)
wf.meta <- weightfunct(meta_data$d, meta_data$vd, table = T)
wf.meta
#indicates there is no publication bias