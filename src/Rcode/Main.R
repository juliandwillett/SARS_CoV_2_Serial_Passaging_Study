setwd("~/ThesisWork/ViralCellMouseStudies")
library(VariantAnnotation) ; library(glue) ; library(magrittr)
library(tidyverse) ; library(vroom) ; library(ggpubr)
library(ggrepel) ; library(tidyr)
source('Functions.R') ; source('Gene_Functions.R')

`%notin%` <- Negate(`%in%`)

######################

#Get data for merged data
all.organs.merged.batch = checkVariantsLate(3,print.vep.for.vep=F,include.primers=F)
lung.samples = list(all.organs.merged.batch$Beta$Lung,all.organs.merged.batch$Delta$Lung)
p0.p10.samples = checkVariantsEarly(print.vep.for.vep=F,include.primers=F)
maf.changes = getSingleSampleMAFChange(lung.samples,p0.p10.samples)

#####################
# Fig 2: Comparing MAF by passage/sample
makeVariantMAFPlotPassage(all.organs.merged.batch,p0.p10.samples,4)
makeVariantMAFPlotSample(all.organs.merged.batch,4)

####################
# Fig 3: Compare gene and VEP proportions
gene = plotVectorChangeSameOrgan(all.organs.merged.batch,p0.p10.samples,
                                 4,'Gene',lims=c(0,0.5)) #nonsig
severity = plotVectorChangeSameOrgan(all.organs.merged.batch,p0.p10.samples,
                                 4,'VEPSeverity',lims=c(0,1)) #nonsig

#####################
# Fig 4: single variant change across passages
# Identify variants that change in frequency during the study
changing.variants = getChangingVariants(maf.changes)

#now plot variant allele changes
plotSingleSampleMAFChange(maf.changes,c('A4981G','C22674T','G23593T'))