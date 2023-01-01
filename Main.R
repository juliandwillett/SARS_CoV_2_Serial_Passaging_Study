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
# Identify variants that change in frequency during the study
changing.variants = getChangingVariants(maf.changes)

#now plot variant allele changes
plotSingleSampleMAFChange(maf.changes,c('A4981G','A21651C','C22674T'))
plotSingleSampleMAFChange(maf.changes,c('C23606T','T27638C','C28253T'))


# variant.depths = makeVariantDepthPlotGeneral(NA,NA,all.organs.merged.batch,no.early=T,4)


allele.changes = getAlleleChanges(all.organs.merged.batch)
# allele.changes.withprimers = getAlleleChanges(all.organs.merged.batch.includingprimers) #to check for dropped variants
allele.changes.all.organs = makeDFForExcel(allele.changes,p0.p10.samples)

table.with.reads = makeTableWithReadsByVariant(all.organs.merged.batch) 

# makeVariantPlot(allele.changes.all.organs,c('A4981G','GTCTGGTTTTA11287GA','C22674T'))
# makeVariantPlot(allele.changes.all.organs,c('G23012A','A23063T','G28368A'))
# makeVariantPlot(allele.changes.all.organs,c('C21614T','A21801C','A22206G'))
# makeVariantPlot(allele.changes.all.organs,c('G23012A','A23063T','C24237T'))

#Variants to highlight: C10809T

#plots for Louis
plotSingleSampleMAFChange(maf.changes,c('G210T','C751T','C5184T','C6449T'))
plotSingleSampleMAFChange(maf.changes,c('C9891T','C10809T','T11418C','C11514T'))
plotSingleSampleMAFChange(maf.changes,c('A12555C','A13192G','G15451A','C16338T'))
plotSingleSampleMAFChange(maf.changes,c('C16466T','C21618G','C22227T','T22917G'))
plotSingleSampleMAFChange(maf.changes,c('C22995A','G23012A','A23063T','C23604G'))
plotSingleSampleMAFChange(maf.changes,c('C25469T','C26250T','T26767C','T27638C'))
plotSingleSampleMAFChange(maf.changes,c('C27752T','G28368A','G28881T','G29402T'))

plotLocationAllelesThatChange(allele.changes)

#focus on lung, track changes in each sample vs mouse pop

plotSingleSampleMAFChange(maf.changes,c('A4981G','GTCTGGTTTTA11287GA','C22674T'))

p13.gene = plotVectorChangeSamePassage(all.organs.merged.batch,3,'Gene',lims=c(0,0.5))

#################################
######Gets variant info 

p0.p10.variants.firstbatch = checkVariantsEarly(1,print.vep.for.vep=F)
all.organs.first.batch = checkVariantsLate(1,print.vep.for.vep = F)
all.organs.p0.removed.first.batch = removeOriginalVariants(p0.p10.variants.firstbatch[[1]],all.organs.first.batch)

p0.p10.variants.secondbatch = checkVariantsEarly(2,print.vep.for.vep=F)
all.organs.second.batch = checkVariantsLate(2,print.vep.for.vep = F)
all.organs.p0.removed.second.batch = removeOriginalVariants(p0.p10.variants.firstbatch[[1]],all.organs.second.batch)

#Regarding comparing batch 1 and batch 2 samples, this was done using
#bcftools stats

#make plot of variant allele depth before removing p0: Use first batch p0 for second batch as not seq
variant.depths.general.firstbatch = makeVariantDepthPlotGeneral(p0.p10.variants.firstbatch[[1]],p0.p10.variants.firstbatch[[2]],all.organs.first.batch)
variant.depths.general.secondbatch = makeVariantDepthPlotGeneral(p0.p10.variants.firstbatch[[1]],p0.p10.variants.secondbatch[[2]],all.organs.second.batch)

variant.depths.lung.firstbatch = makeVariantDepthPlotOrgan(all.organs.first.batch,4)
variant.depths.lung.secondbatch = makeVariantDepthPlotOrgan(all.organs.second.batch,4)

variant.depths.lung = makeVariantDepthPlotOrgan(all.organs.second.batch,4,30000)

#plot of variant allele depth after removing p0
variant.depth.lung.p0rem.firstbatch = makeVariantDepthPlotOrgan(all.organs.p0.removed.second.batch,4)
variant.depth.lung.p0rem.secondbatch = makeVariantDepthPlotOrgan(all.organs.p0.removed.second.batch,4)

#next plot how variables change by passage
#studying this at population level vs single sample changes
#as looking at how this changes as a whole
l.shift.firstbatch = plotVectorChangeSameOrgan(all.organs.p0.removed.first.batch,4,'VcfForm',lims=c(0,0.3)) #nonsig
l.gene.firstbatch = plotVectorChangeSameOrgan(all.organs.p0.removed.first.batch,4,'Gene',lims=c(0,0.75)) #nonsig
l.severity.firstbatch = plotVectorChangeSameOrgan(all.organs.p0.removed.first.batch,4,'VEPSeverity',lims=c(0,1)) 

l.shift.secondbatch = plotVectorChangeSameOrgan(all.organs.p0.removed.second.batch,4,'VcfForm',lims=c(0,0.3)) #nonsig
l.gene.secondbatch = plotVectorChangeSameOrgan(all.organs.p0.removed.second.batch,4,'Gene',lims=c(0,0.75)) #nonsig
l.severity.secondbatch = plotVectorChangeSameOrgan(all.organs.p0.removed.second.batch,4,'VEPSeverity',lims=c(0,1)) 

#Next get the novel variants appearing each passage
novel.variants = isolateNewVariantAlleleByPassage(all.organs.p0.removed)

#pull human data
human.mutation.bases = readRDS('~/ThesisWork/ClosedProjects/Sars_Cov_2_GenomesProject/final.mutation.bases.RDS')
human.mutation.bases.all = human.mutation.bases[[1]] %>% add_row(human.mutation.bases[[2]]) %>%
  add_row(human.mutation.bases[[3]]) %>% add_row(human.mutation.bases[[4]]) %>% 
  dplyr::filter(Lineage == 'B.1.617.2' | Lineage == 'B.1.351')
summary.file = vroom('~/ThesisWork/ClosedProjects/Sars_Cov_2_GenomesProject/LATEST_REPORT.csv')
