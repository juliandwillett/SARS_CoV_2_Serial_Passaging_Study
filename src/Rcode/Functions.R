makeAwkLinesRemovePrimers = function() {
  loc = read.table('Primers_Positions.bed')
  lines = character()
  for (item in 1:nrow(loc)) {
    r1 = loc$V2[[item]] ; r2 = loc$V3[[item]]
    if (item %% 2 == 1)
      lines[[item]] = paste0("awk '{if ($2 < ",r1," || $2 > ",r2,") print $0}' tmp > tmp2")
    else
      lines[[item]] = paste0("awk '{if ($2 < ",r1," || $2 > ",r2,") print $0}' tmp2 > tmp")
  }
  write.table(lines,'awklines.txt',quote=F,row.names=F,col.names=F)
}
checkVariantsEarly = function(print.vep.for.vep,include.primers) {
  variants = list()
  
  for (pass in 1:2) {
    variants[[pass]] = list()
    for (lin in 1:2) {
      #beta is SA, Delta is IN
      if (lin == 1) lin.label = 'SA'
      else lin.label = 'IN'
      
      if (pass == 1) variants[[pass]][[lin]] = getVariants(batch=4,glue('Sample_{lin.label}_SARS-COV2'),
                                            print.vep.for.vep,include.primers)
      else variants[[pass]][[lin]] = getVariants(batch=4,glue('Sample_{lin.label}_P10'),
                                                 print.vep.for.vep,include.primers)
    }
    names(variants[[pass]]) = c('Beta','Delta')
  }
  names(variants) = c('P0','P10')
  
  return(variants)
}
checkVariantsLate = function(batch,print.vep.for.vep,include.primers) {
  out.list = list()
  
  for (lineage in 1:2) {
    out.list[[lineage]] = list()
    curr.lineage = c('SA','DEL')[[lineage]]
    for (organ in 1:5) {
      out.list[[lineage]][[organ]] = list()
      curr.organ = c('B','H','K','L','S')[[organ]]
      for (passage in 1:3) {
        out.list[[lineage]][[organ]][[passage]] = list()
        curr.passage = c('13','17','20')[[passage]]
        for (sample in 1:3) {
          if (batch==1)
            variants = getVariants(batch=1,glue('Sample_P{curr.passage}-{curr.lineage}-{sample}-{curr.organ}_2'),print.vep.for.vep,include.primers)
          else if (batch == 2)
            variants = getVariants(batch=2,glue('P{curr.passage}-{curr.lineage}-{sample}-{curr.organ}_2'),print.vep.for.vep,include.primers)
          else if (batch == 3)
            variants = getVariants(batch=3,glue('P{curr.passage}-{curr.lineage}-{sample}-{curr.organ}_2'),print.vep.for.vep,include.primers)
          out.list[[lineage]][[organ]][[passage]][[sample]] = variants
        }
        names(out.list[[lineage]][[organ]][[passage]]) = c('Sample1','Sample2','Sample3')
      }
      names(out.list[[lineage]][[organ]]) = c('P13','P17','P20')
    }
    names(out.list[[lineage]]) = c('Brain','Heart','Kidney','Lung','Spleen')
  }
  names(out.list) = c('Beta','Delta')
  return(out.list)
}
removeOriginalVariants = function(p0,samples) { #not considering p10 b/c unclear sampling
  out.samples = samples
  
  for (lin in 1:2) {
    for (organ in 1:5) {
      curr.p0 = p0[[lin]]
      p0.vcf.form = paste(curr.p0$Ref,curr.p0$Pos,curr.p0$Alt)
      for (pass in 1:3) {
        print(glue('On organ {organ} lin {lin} passage {pass}'))
        for (s in 1:3) {
          cs = samples[[lin]][[organ]][[pass]][[s]] %>%
            dplyr::mutate(VcfForm=paste(Ref,Pos,Alt)) %>%
            dplyr::filter(VcfForm %notin% p0.vcf.form)
          cs.vep = getVEPData(cs,print.vep.for.vep=F)
          if (length(which(is.na(cs.vep$VEPConsequence)))>0) stop()
          out.samples[[lin]][[organ]][[pass]][[s]] = cs.vep
        }
      }
    }
  }
  return(out.samples)
}
getVariants = function(batch,sample,print.vep.for.vep,include.primers) {
  #print.vep.for.vep outputs missing vep lines
  ###if T, can be pasted into VEP. Otherwise, prints line to add to
  ###df in the "Gene_Functions.R" df that has most consequential
  ###change from VEP logged
  df = data.frame(Sample=character(),Pos=numeric(),Ref=character(),
                  Alt=character(),Depth=numeric(),MAF=numeric(),
                  ReadsSupporting=numeric(),
                  Gene=character(),AAChange=character(),
                  VEPConsequence=character(),VEPSeverity=character(),
                  VcfForm=character())
  if (batch == 1) f = list.files('Vcfs_FirstRun_PrimersFiltered',sample,full.names=T)
  else if (batch == 2) f = list.files('Vcfs_SecondRun_PrimersFiltered',sample,full.names=T)
  else if (batch == 3 & !include.primers) f = list.files('Vcfs_Merged_PrimersFiltered_wHeader',sample,full.names=T)
  else if (batch == 3 & include.primers) f = list.files('Vcfs_Merged',sample,full.names=T)
  else if (batch == 4) f = list.files('Vcfs_Merged_P0_P10_PrimersFiltered_wHeader',sample,full.names=T)
  
  print(glue('Reading file: {f}'))
  
  data = readVcf(f)
  refs = as.character(data@fixed@listData[["REF"]])
  alts = getAlts(data)
  positions = data@rowRanges@ranges@start
  depths = data@info@listData[["DP"]]
  mafs = data@info@listData[["AF"]]@unlistData
  
  match = which(str_detect(alts,','))
  if (length(match) > 0) { #some loci have multiple alleles, this deals with them
    indices.to.drop = numeric()
    maf.indices.to.drop = numeric()
    for (m in match) { #for cases of multiple matches, index in vector
      spl = str_split(alts[[m]],',')[[1]]
      if (length(spl)>2) stop('More than 2 alleles')
      refs %<>% append(refs[[m]]) %>% append(refs[[m]])
      alts %<>% append(spl[[1]]) %>% append(spl[[2]])
      positions %<>% append(positions[[m]]) %>% append(positions[[m]])
      depths %<>% append(depths[[m]]) %>% append(depths[[m]])
      mafs %<>% append(mafs[[m]]) %>% append(mafs[[m+1]])
      indices.to.drop %<>% append(m)
      maf.indices.to.drop %<>% append(m) %>% append(m+1)
    }
    refs = refs[-indices.to.drop]
    alts = alts[-indices.to.drop]
    positions = positions[-indices.to.drop]
    depths = depths[-indices.to.drop]
    mafs = mafs[-maf.indices.to.drop]
  }

  df %<>% add_row(Sample=sample,Ref=refs,Pos=positions,Alt=alts,
                  Depth=depths,MAF=mafs,ReadsSupporting=depths*mafs,
                  Gene=NA,AAChange=NA,
                  VEPConsequence=NA,VcfForm=glue("{refs}{positions}{alts}"))
  if (print.vep.for.vep) return(getVEPData(df,print.vep.for.vep))
  else {
    df = getVEPData(df,print.vep.for.vep=F)
    #if (length(which(is.na(df$VEPConsequence)))>0) stop()
    return(df %>% dplyr::mutate(Shift=glue('{Ref}>{Alt}')))
  }
}
getAlts = function(data) {
  alts = character()
  len = length(data@fixed@listData[["ALT"]])
  for (item in 1:len) {
    curr.alt = as.character(data@fixed@listData[["ALT"]][[item]])
    if (length(curr.alt)==1) alts[[item]] = curr.alt
    else alts[[item]] = glue('{curr.alt[[1]]},{curr.alt[[2]]}')
  }
  return(alts)
}
makeVariantMAFPlotPassage = function(data,pre,organ) {
  #ggplot for signif: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  plts = list()
  stats.out = list()
  
  for (lineage in 1:2) { #make figure for each lineage
    plt.df = data.frame(Passage=character(),Depth=numeric(),MAF=numeric())
    
    for (p in 1:3) { #go through each lung passage (13/17/20)
      for (s in 1:3) { #go through each sample in lung
        if (p==1) passage='P13' 
        else if (p==2) passage='P17'
        else passage='P20'
        
        tissue = data[[lineage]][[organ]]
        
        plt.df %<>% add_row(Passage=passage,Depth=tissue[[p]][[s]]$Depth,
                            MAF = tissue[[p]][[s]]$MAF)
      }
    }
    #Now add P0/P10 data
    plt.df %<>% add_row(Passage='P0',Depth=pre[[1]][[lineage]]$Depth,
                        MAF=pre[[1]][[lineage]]$MAF) %>%
      add_row(Passage='P10',Depth=pre[[2]][[lineage]]$Depth,
              MAF=pre[[2]][[lineage]]$MAF)
    
    plts[[lineage]] = plt.df
    comparisons = list(c('P0','P10'),c('P10','P13'),c('P10','P17'),
                       c('P10','P20'),c('P13','P17'),c('P13','P20'),
                       c('P17','P20'))

    comp.df = data.frame(Comparison=character(),PVal=numeric(),AdjPVal=numeric())
    # comp.df.depth = data.frame(Comparison=character(),PVal=numeric(),AdjPVal=numeric())
    for (comp in comparisons) {
      comp.df %<>% add_row(Comparison=paste(comp,collapse=' '),
                           PVal=t.test(plt.df$MAF[which(plt.df$Passage==comp[[1]])],
                                       plt.df$MAF[which(plt.df$Passage==comp[[2]])])$p.value,
                           AdjPVal=NA)
      # comp.df.depth %<>% add_row(Comparison=paste(comp,collapse=' '),
      #                      PVal=t.test(plt.df$Depth[which(plt.df$Passage==comp[[1]])],
      #                                  plt.df$Depth[which(plt.df$Passage==comp[[2]])])$p.value,
      #                      AdjPVal=NA)
    }
    comp.df$AdjPVal = round(p.adjust(comp.df$PVal,method='BH'),2)
    # comp.df %<>% dplyr::mutate(StarSignif=as.character(AdjPVal)) %>%
    #   dplyr::mutate(StarSignif=ifelse(as.numeric(StarSignif)>0.05,'n.s.',StarSignif)) %>%
    #   dplyr::mutate(StarSignif=ifelse(as.numeric(StarSignif)>0.01,'*',StarSignif)) %>%
    #   dplyr::mutate(StarSignif=ifelse(as.numeric(StarSignif)<=0.01,'**',StarSignif))
    
    # comp.df.depth$AdjPVal = round(p.adjust(comp.df.depth$PVal,method='BH'),2)
    if (lineage==1) print('Beta') else print('Delta')
    print('MAF stats')
    print(comp.df)
    # print('Depth stats')
    # print(comp.df.depth)
    # stats.out[[lineage]] = list(comp.df,comp.df.depth)
    stats.out[[lineage]] = comp.df
    num.maf.one = c(nrow(plt.df %>% dplyr::filter(Passage=='P0',MAF==1)),
                    nrow(plt.df %>% dplyr::filter(Passage=='P10',MAF==1)),
                    nrow(plt.df %>% dplyr::filter(Passage=='P13',MAF==1)),
                    nrow(plt.df %>% dplyr::filter(Passage=='P17',MAF==1)),
                    nrow(plt.df %>% dplyr::filter(Passage=='P20',MAF==1)))
    num.maf.half = c(nrow(plt.df %>% dplyr::filter(Passage=='P0',MAF==0.5)),
                    nrow(plt.df %>% dplyr::filter(Passage=='P10',MAF==0.5)),
                    nrow(plt.df %>% dplyr::filter(Passage=='P13',MAF==0.5)),
                    nrow(plt.df %>% dplyr::filter(Passage=='P17',MAF==0.5)),
                    nrow(plt.df %>% dplyr::filter(Passage=='P20',MAF==0.5)))
    
    plt = ggplot(plt.df,aes(x=Passage,y=MAF)) + geom_boxplot(outlier.shape=NA) +
            geom_point(aes(color=Depth),size=5) + stat_boxplot(geom='errorbar') +
            theme_bw() + ylim(c(0.5,1.5)) +
            stat_compare_means(comparisons=comparisons,size=6,label.y=c(1.03,1.05,1.12,1.19,1.26,1.33,1.4)) +
            # stat_compare_means(label.y=1.05,size=6) +
            annotate('text',x=1:length(table(plt.df$Passage)),y=0.75,size=6,label=glue('n={table(plt.df$Passage)}')) +
            annotate('text',x=1:length(table(plt.df$Passage)),y=0.95,size=6,label=glue('n={num.maf.one}')) +
            annotate('text',x=1:length(table(plt.df$Passage)),y=0.55,size=6,label=glue('n={num.maf.half}')) +
    guides(fill='none') + theme(text = element_text(size=20))
    plt = ggplot_build(plt)
    plt$data[[4]]$annotation[1:3] = round(comp.df$AdjPVal[1],3)
    plt$data[[4]]$annotation[4:6] = round(comp.df$AdjPVal[2],3)
    plt$data[[4]]$annotation[7:9] = round(comp.df$AdjPVal[3],3)
    plt$data[[4]]$annotation[10:12] = round(comp.df$AdjPVal[4],3)
    plt$data[[4]]$annotation[13:15] = round(comp.df$AdjPVal[5],3)
    plt$data[[4]]$annotation[16:18] = round(comp.df$AdjPVal[6],3)
    plt$data[[4]]$annotation[19:21] = round(comp.df$AdjPVal[7],3)
    
    print(plt)
  }
  return(list(plts,stats.out))
}
makeVariantMAFPlotSample = function(data,organ.num) {
  plts = list()
  for (lineage in 1:2) {
    plt.df = data.frame(Passage=character(),Sample=character(),Depth=numeric(),
                        MAF=numeric())
    
    for (p in 1:3) { #go through each passage (13/17/20)
      for (s in 1:3) { #go through each sample in lung
        if (p==1) passage='P13' 
        else if (p==2) passage='P17'
        else passage='P20'
        
        plt.df %<>% add_row(Passage=passage,Sample=as.character(s),
                            Depth=data[[lineage]][[organ.num]][[p]][[s]]$Depth,
                            MAF=data[[lineage]][[organ.num]][[p]][[s]]$MAF)
      }
    }
    
    plts[[lineage]] = plt.df
    P13.anova.p = summary(aov(MAF~Sample,data = plt.df %>% dplyr::filter(Passage=='P13')))[[1]][['Pr(>F)']][[1]]
    P17.anova.p = summary(aov(MAF~Sample,data = plt.df %>% dplyr::filter(Passage=='P17')))[[1]][['Pr(>F)']][[1]]
    P20.anova.p = summary(aov(MAF~Sample,data = plt.df %>% dplyr::filter(Passage=='P20')))[[1]][['Pr(>F)']][[1]]
    adj.p.vals = round(p.adjust(c(P13.anova.p,P17.anova.p,P20.anova.p),
                          method='BH'),3)
    print(lineage) ; print('MAF') ; print(adj.p.vals)
    
    P13.anova.p.depth = summary(aov(Depth~Sample,data = plt.df %>% dplyr::filter(Passage=='P13')))[[1]][['Pr(>F)']][[1]]
    P17.anova.p.depth = summary(aov(Depth~Sample,data = plt.df %>% dplyr::filter(Passage=='P17')))[[1]][['Pr(>F)']][[1]]
    P20.anova.p.depth = summary(aov(Depth~Sample,data = plt.df %>% dplyr::filter(Passage=='P20')))[[1]][['Pr(>F)']][[1]]
    adj.p.vals.depth = round(p.adjust(c(P13.anova.p.depth,P17.anova.p.depth,P20.anova.p.depth),
                                method='BH'),3)
    print(lineage) ; print('Depth') ; print(adj.p.vals.depth)
    
    print(ggplot(plt.df,aes(x=Passage,y=MAF,fill=Sample)) +
            geom_boxplot(outlier.shape=NA) + stat_boxplot(geom='errorbar') +
            geom_point(position=position_dodge(width=0.8),aes(color=Depth),size=5) + theme_bw() +
            # stat_compare_means(aes(group=Sample),size=6,label='p.format',label.y=0.3) +
            theme(text = element_text(size=20)) + ylim(c(0,1)) +
            # annotate('text',x=1:length(table(plt.df$Passage)),size=6,y=0.75,label=glue('n={table(plt.df$Passage)}')) +
            # annotate('text',x=1:length(table(plt.df$Passage)),y=0.95,size=6,label=glue('n={table(plt.df$Passage[which(plt.df$MAF == 1)])}')) +
            # annotate('text',x=1:length(table(plt.df$Passage)),y=0.55,size=6,label=glue('n={table(plt.df$Passage[which(plt.df$MAF == 0.5)])}')) +
            annotate('text',x=1:length(table(plt.df$Passage)),size=6,y=0.3,label=glue('p={adj.p.vals}')) +
            guides(fill='none') + geom_hline(yintercept = 0.35))
  }
  return(plts)
}
plotVectorChangeSameOrgan = function(data,pre,organ.num,col,lims) {
  plts = stats = list()
  for (lin in 1:2) {
    plt.df = data.frame(Passage=character(),Sample=character(),Entry=character(),
                        PassageProportion=numeric(),Position=numeric())
    for (pass in 1:3) {
      for (s in 1:3) {
        entries = data[[lin]][[organ.num]][[pass]][[s]][[col]]
        curr.pass = c('P13','P17','P20')[[pass]]
        plt.df %<>% add_row(Passage=curr.pass,Sample=as.character(s),
                            Entry=entries,PassageProportion=NA,
                            Position=data[[lin]][[organ.num]][[pass]][[s]]$Pos)
      }
    }
    # now add p0 and p10
    plt.df %<>% add_row(Passage='P0',Sample=NA,Entry=pre[[1]][[lin]][[col]],
                        PassageProportion=NA,Position=pre[[1]][[lin]]$Pos) %>%
      add_row(Passage='P10',Sample=NA,Entry=pre[[2]][[lin]][[col]],
              PassageProportion=NA,Position=pre[[2]][[lin]]$Pos)
    
    #now get the proportions by passage
    for (pass in c('P0','P10','P13','P17','P20')) {
      curr.rows = which(plt.df$Passage == pass)
      curr.entries = plt.df$Entry[curr.rows]
      props = numeric()
      
      for (e in 1:length(curr.entries)) {
        prop = length(which(curr.entries[[e]] == curr.entries)) / length(curr.entries)
        props %<>% append(prop)
      }
      plt.df[curr.rows,'PassageProportion'] = props
    }
    plt.df$Passage = factor(plt.df$Passage,
                            levels=c('P0','P10','P13','P17','P20'))
    
    if (col == 'VEPSeverity') plt.df$Entry = factor(plt.df$Entry,levels=c('Modifier','Low','Moderate','High'))

    if (col == 'VcfForm') {
      long.entries = which(nchar(plt.df$Entry) > 15)
      plt.df$Entry[long.entries] = glue('X {plt.df$Position[long.entries]} X')
      plt.df$Position = as.numeric(str_sub(plt.df$Entry,3,nchar(plt.df$Entry)-2))
      plt.df %<>% dplyr::arrange(Position)
      vcf.order = unique(plt.df$Entry)
    }

    plts[[lin]] = plt.df
    
    #Now report the stats:
    stats[[lin]] = compareVectorsSameOrgan(data[[lin]],pre[[1]][[lin]],
                                           pre[[2]][[lin]],organ.num,col)
  }
  facet.plt = (plts[[1]] %>% dplyr::mutate(Lineage='Beta')) %>%
    add_row(plts[[2]] %>% dplyr::mutate(Lineage='Delta'))
  
  plt = ggplot(facet.plt,aes(x=PassageProportion,y=Entry,fill=Passage)) + 
    geom_bar(stat='identity',position=position_dodge()) + 
    xlab('Proportion') + xlim(c(lims[[1]],lims[[2]])) +
    theme_bw() + theme(text = element_text(size=20)) + ylab('') +
    facet_grid(cols=vars(Lineage))
  if (col == 'VcfForm') plt = plt + scale_y_discrete(limits = rev(vcf.order))
  else plt = plt + scale_y_discrete(limits=rev)
  print(plt)
  
  return(list(plts,stats))
}
compareVectorsSameOrgan = function(df,p0,p10,organ.num,col) {
  #takes single lineage data
  out.df = data.frame(Comparison=character(),PVal=numeric(),AdjPVal=numeric())
  test.df = data.frame(Entry=character(),Passage=character())
  for (pass in c('P13','P17','P20')) { #go by passage
    for (s in 1:3) { #go by sample
      test.df %<>% add_row(Entry=as.character(df[[organ.num]][[pass]][[s]][[col]]),
                           Passage=pass)
    }
  }
  test.df %<>% add_row(Entry=p0[[col]],Passage='P0') %>%
    add_row(Entry=p10[[col]],Passage='P10')
  
  for (comp in c('P0vP10','P10vP13','P10vP17','P10vP20','P13vP17',
                 'P13vP20','P17vP20')) {
    curr.passages = str_split(comp,'v')[[1]]
    curr.df = test.df %>% dplyr::filter(Passage %in% curr.passages)
    pval = chisq.test(table(curr.df$Passage,curr.df$Entry))$p.value
    out.df %<>% add_row(Comparison=comp,PVal=pval,AdjPVal=NA)
  }
  
  out.df$AdjPVal = p.adjust(out.df$PVal,method='BH')
  if (length(which(out.df$AdjPVal <= 0.05)) > 0) print('Significant association')
  else print('No statistically significant associations')
  
  return(out.df)
}
compareVectorsSameOrganDifferentLineages = function(df,organ.num,col) {
  out.df = data.frame(Comparison=character(),PVal=numeric(),AdjPVal=numeric())
  test.df = data.frame(Lineage=character(),Entry=character(),Passage=numeric())
  for (lin in 1:2) {
    for (pass in 1:3) { #go by passage
      for (s in 1:3) { #go by sample
        test.df %<>% add_row(Lineage=as.character(lin),
                             Entry=as.character(df[[lin]][[organ.num]][[pass]][[s]][[col]]),
                             Passage=pass)
      }
    }
  }
  for (item in 1:3) { #compare each set
    curr.df = test.df %>% dplyr::filter(Passage == item)
    beta.df = curr.df %>% dplyr::filter(Lineage == "1")
    delta.df = curr.df %>% dplyr::filter(Lineage == "2")
    
    beta.table = table(beta.df$Passage,beta.df$Entry)
    delta.table = table(delta.df$Passage,delta.df$Entry)
    if (ncol(beta.table) < ncol(delta.table)) beta.table %<>% cbind(High=0)
    
    pval = chisq.test(beta.table[1,],delta.table[1,])$p.value
    out.df %<>% add_row(Comparison=as.character(item),PVal=pval,AdjPVal=NA)
  }
  out.df$AdjPVal = p.adjust(out.df$PVal,method='BH')
  if (length(which(out.df$AdjPVal <= 0.05)) > 0) print('Significant association')
  else print('No statistically significant associations')
  
  return(out.df)
}
getAlleleChanges = function(data,p0) {
  prelim.df = data.frame(Variant=character(),Position=numeric(),
                              Gene=character(),AA=character(),Lineage=character(),
                              Organ=character(),Passage=numeric())
  df = data.frame(Variant=character(),Position=numeric(),
                  Gene=character(),AA=character(),Lineage=character(),
                  Organ=character(),P13.MAF=numeric(),P17.MAF=numeric(),
                  P20.MAF=numeric())
  
  #first get all the variants
  for (lin in 1:length(data)) {
    for (organ in 1:length(data[[lin]])) {
      for (pass in 1:length(data[[lin]][[organ]])) {
        for (s in 1:length(data[[lin]][[organ]][[pass]])) {
          vcf.forms = data[[lin]][[organ]][[pass]][[s]][['VcfForm']]
          prelim.df %<>% add_row(Variant = vcf.forms,
                        Position=data[[lin]][[organ]][[pass]][[s]][['Pos']],
                                 Gene=data[[lin]][[organ]][[pass]][[s]][['Gene']],
                        AA=data[[lin]][[organ]][[pass]][[s]][['AAChange']],
                          Lineage=as.character(lin),Organ=as.character(organ),
                        Passage=pass)
        }
      }
    }
  }
  unique.variants = unique(prelim.df$Variant)
  
  #now get the AF
  for (var in 1:length(unique.variants)) {
    print(glue('On var {var} of {length(unique.variants)}'))
    for (lin in 1:2) {
      for (organ in 1:5) {
        maf.vec = numeric()
        for (pass in 1:3) {
          tmp = prelim.df %>% dplyr::filter(Variant==unique.variants[[var]],
                                            Lineage==lin,Organ==organ,Passage==pass)
          maf.vec %<>% append(nrow(tmp)/3)
        }
        df %<>% add_row(Variant=unique.variants[[var]],Lineage=as.character(lin),
                        Organ=as.character(organ),P13.MAF=maf.vec[[1]],
                        P17.MAF=maf.vec[[2]],P20.MAF=maf.vec[[3]])
      }
    }
  }
  
  #now do stats on proportion change
  df %<>% dplyr::mutate(PVal=NA,AdjPVal=NA)
  for (row in 1:nrow(df)) {
    df$PVal[[row]] = prop.test(c(df$P13.MAF[[row]],df$P17.MAF[[row]],df$P20.MAF[[row]]),
                               c(1,1,1))$p.value
    if (is.na(df$PVal[[row]])) df$PVal[[row]] = 1 #ie no change
  }
  df$AdjPVal = p.adjust(df$PVal,method='BH')
  df$Lineage = str_replace(df$Lineage,'1','Beta')
  df$Lineage = str_replace(df$Lineage,'2','Delta')
  df$Organ = str_replace(df$Organ,'1','Brain')
  df$Organ = str_replace(df$Organ,'2','Heart')
  df$Organ = str_replace(df$Organ,'3','Kidney')
  df$Organ = str_replace(df$Organ,'4','Lung')
  df$Organ = str_replace(df$Organ,'5','Spleen')
  
  #now get the gene and AA change
  print('Getting gene and AA changes')
  for (row in 1:nrow(df)) {
    df$Position[[row]] = prelim.df[which(prelim.df$Variant==df$Variant[[row]]),'Position'][1]
    df$Gene[[row]] = prelim.df[which(prelim.df$Variant==df$Variant[[row]]),'Gene'][1]
    df$AA[[row]] = prelim.df[which(prelim.df$Variant==df$Variant[[row]]),'AA'][1]
  }
  
  return(df %>% dplyr::arrange(Lineage,Organ))
}
isolateNewVariantAlleleByPassage = function(data) {
  #Sample 1 was the same across passages, so going by sample vs aggregating
  new.variants = list()
  for (organ in 1:5) {
    new.variants[[organ]] = list()
    for (lineage in 1:2) {
      new.variants[[organ]][[lineage]] = list()
      for (sample in 1:3) {
        p13.vars = data[[organ]][[lineage]][[1]][[sample]][['VcfForm']]
        p17.vars = data[[organ]][[lineage]][[2]][[sample]][['VcfForm']]
        new.variants[[organ]][[lineage]][[1]] = data[[organ]][[lineage]][[1]][[sample]] %>% 
          dplyr::filter(VcfForm %notin% p13.vars)
        new.variants[[organ]][[lineage]][[2]] = data[[organ]][[lineage]][[2]][[sample]] %>% 
          dplyr::filter(VcfForm %notin% p13.vars)
        new.variants[[organ]][[lineage]][[3]] = data[[organ]][[lineage]][[3]][[sample]] %>% 
          dplyr::filter(VcfForm %notin% p13.vars,VcfForm %notin% p17.vars)
      }
    }
  }
  return(new.variants)
}
getVariantByAA = function(data,lin,gene,aa) {
  variants = character()
  for (organ in 1:length(data[[lin]])) {
    for (pass in 1:length(data[[lin]][[organ]])) {
      for (s in 1:length(data[[lin]][[organ]][[pass]])) {
        match = data[[lin]][[organ]][[pass]][[s]] %>% dplyr::filter(Gene==gene,AAChange==aa)
        if (nrow(match)==0) next
        else variants %<>% append(match$VcfForm[[1]])
      }
    }
  }
  print(variants)
}
makeDFForExcel = function(a.c,pre) { #take allele changes
  df.list = list()
  for (lin in c('Beta','Delta')) {
    curr.df = a.c %>% dplyr::filter(Lineage == lin)
    out.df = data.frame(Variant=character(),Gene=character(),AA=character(),
                        Lineage=character(),B.P13.MAF=numeric(),B.P17.MAF=numeric(),
                        B.P20.MAF=numeric(),H.P13.MAF=numeric(),H.P17.MAF=numeric(),
                        H.P20.MAF=numeric(),K.P13.MAF=numeric(),K.P17.MAF=numeric(),
                        K.P20.MAF=numeric(),L.P13.MAF=numeric(),L.P17.MAF=numeric(),
                        L.P20.MAF=numeric(),S.P13.MAF=numeric(),S.P17.MAF=numeric(),
                        S.P20.MAF=numeric())
    for (var in unique(curr.df$Variant)) {
      b.data = curr.df %>% dplyr::filter(Variant==var,Organ=='Brain')
      h.data = curr.df %>% dplyr::filter(Variant==var,Organ=='Heart')
      k.data = curr.df %>% dplyr::filter(Variant==var,Organ=='Kidney')
      l.data = curr.df %>% dplyr::filter(Variant==var,Organ=='Lung')
      s.data = curr.df %>% dplyr::filter(Variant==var,Organ=='Spleen')
      
      out.df %<>% add_row(Variant=var,Gene=b.data$Gene,AA=b.data$AA,Lineage=lin,
                          B.P13.MAF=b.data$P13.MAF,B.P17.MAF=b.data$P17.MAF,
                          B.P20.MAF=b.data$P20.MAF,H.P13.MAF=h.data$P13.MAF,
                          H.P17.MAF=h.data$P13.MAF,H.P20.MAF=h.data$P20.MAF,
                          K.P13.MAF=k.data$P13.MAF,K.P17.MAF=k.data$P17.MAF,
                          K.P20.MAF=k.data$P20.MAF,L.P13.MAF=l.data$P13.MAF,
                          L.P17.MAF=l.data$P17.MAF,L.P20.MAF=l.data$P20.MAF,
                          S.P13.MAF=s.data$P13.MAF,S.P17.MAF=s.data$P17.MAF,
                          S.P20.MAF=s.data$P20.MAF)
    }
    for (col in 5:ncol(out.df)) {
      out.df[[col]] = round(out.df[[col]],1)
    }
    out.df %<>% dplyr::mutate(Position=as.numeric(gsub("[^0-9.-]", "", Variant))) %>%
      dplyr::arrange(Position)
    
    out.df %<>% dplyr::mutate(PresentP0 = out.df$Variant %in% pre[[1]][[lin]][['VcfForm']],
                              PresentP10 = out.df$Variant %in% pre[[1]][[lin]][['VcfForm']])
    df.list[[lin]] = out.df
  }
  return(df.list)
}
makeTableWithReadsByVariant = function(df) {
  out.df = data.frame(Variant=character(),Lin=character(),Organ=character(),
                      Passage=character(),Sample=character(),
                      NumReads=numeric(),Depth=numeric(),MAF=numeric(),
                      Position=numeric())
  for (lin in 1:length(df)) {
    for (org in 1:length(df[[lin]])) {
      for (pass in 1:length(df[[lin]][[org]])) {
        for (s in 1:length(df[[lin]][[org]][[pass]])) {
          curr.df = df[[lin]][[org]][[pass]][[s]]
          out.df %<>% add_row(Variant=curr.df$VcfForm,Lin=names(df)[[lin]],
                              Organ=names(df[[lin]])[[org]],
                              Passage=names(df[[lin]][[org]])[[pass]],
                              Sample=names(df[[lin]][[org]][[pass]])[[s]],
                              NumReads=curr.df$ReadsSupporting,
                              Depth=curr.df$Depth,MAF=curr.df$MAF,
                              Position=curr.df$Pos)
        }
      }
    }
  }
  return(out.df %>% dplyr::arrange(Position,Lin,Organ,Passage,Sample))
}
makeVariantPlot = function(df,vars) {
  fig.df = data.frame(Variant=character(),Lineage=character(),
                      Passage=character(),Freq=numeric(),
                      Organ=character(),PresentP0=logical(),PresentP10=logical())
  for (lin in 1:length(df)) {
    lin.name = names(df)[[lin]]
    for (var in vars) {
      new.df = df[[lin]] %>% dplyr::filter(Variant==var)
      fig.df %<>% add_row(Variant=var,Lineage=lin.name,Passage='P13',Freq=new.df$B.P13.MAF,Organ='Brain',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P17',Freq=new.df$B.P17.MAF,Organ='Brain',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P20',Freq=new.df$B.P20.MAF,Organ='Brain',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P13',Freq=new.df$H.P13.MAF,Organ='Heart',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P17',Freq=new.df$H.P17.MAF,Organ='Heart',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P20',Freq=new.df$H.P20.MAF,Organ='Heart',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P13',Freq=new.df$K.P13.MAF,Organ='Kidney',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P17',Freq=new.df$K.P17.MAF,Organ='Kidney',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P20',Freq=new.df$K.P20.MAF,Organ='Kidney',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P13',Freq=new.df$L.P13.MAF,Organ='Lung',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P17',Freq=new.df$L.P17.MAF,Organ='Lung',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P20',Freq=new.df$L.P20.MAF,Organ='Lung',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P13',Freq=new.df$S.P13.MAF,Organ='Spleen',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P17',Freq=new.df$S.P17.MAF,Organ='Spleen',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10) %>%
        add_row(Variant=var,Lineage=lin.name,Passage='P20',Freq=new.df$S.P20.MAF,Organ='Spleen',PresentP0=new.df$PresentP0,PresentP10=new.df$PresentP10)
    }
  }
  
  fig.df %<>% dplyr::mutate(PresentP0=ifelse(PresentP0 == TRUE,'Present P0','Not Present P0')) %>%
    dplyr::mutate(PresentP10=ifelse(PresentP10 == TRUE,'Present P10','Not Present P10')) %>%
    dplyr::mutate(Variant=factor(Variant,levels=vars))
  
  ggplot(fig.df,aes(x=Passage,y=Freq,color=Organ,label=glue("{PresentP0}\n{PresentP10}"))) + 
    geom_point(size=3,aes(shape=Organ)) + theme_bw() + geom_line(aes(group=Organ)) +
    theme(text=element_text(size=16)) + ylim(c(0,1)) +
    geom_text(aes(x='P17',y=0.5),col='black') +
    scale_color_brewer(palette='Dark2') +
    scale_shape_manual(values=c(16,2,9,15,7)) +
    facet_grid(rows=vars(Lineage),cols=vars(Variant))
}
plotLocationAllelesThatChange = function(df) {
  df %<>% dplyr::filter(abs(P13.MAF - P17.MAF) > 0.5 | abs(P13.MAF - P20.MAF) > 0.5 |
                          abs(P17.MAF != P20.MAF) > 0.5)
  fig.df = data.frame(Position=df$Position,Lineage=df$Lineage,
                      Organ=df$Organ,Passage='P13',Freq=df$P13.MAF) %>%
    add_row(Position=df$Position,Lineage=df$Lineage,
            Organ=df$Organ,Passage='P17',Freq=df$P17.MAF) %>%
    add_row(Position=df$Position,Lineage=df$Lineage,
            Organ=df$Organ,Passage='P20',Freq=df$P20.MAF)
  fig = ggplot(fig.df,aes(x=Position,y=Freq,color=Passage,shape=Passage)) +
    geom_point(size=2) + theme_bw() + facet_grid(rows=vars(Organ),cols=vars(Lineage))
  print(fig)
  return(fig.df)
}
getSingleSampleMAFChange = function(df,pre) {
  sample.changes = list()
  all.variants = character()
  for (l in 1:2) { #by lineage
    print(glue('On lin {l} of 2'))
    sample.changes[[l]] = list()
    for (s in 1:3) { #by sample
      s.df = data.frame(Variant=character(),Gene=character(),
                        AAChange=character(),Passage=numeric(),
                        MAF=numeric(),Pos=numeric(),MAF_Change=logical())
      for (p in 1:3) { #by passage
        s.df %<>% add_row(Variant=df[[l]][[p]][[s]]$VcfForm,
                          Gene=df[[l]][[p]][[s]]$Gene,
                          AAChange=df[[l]][[p]][[s]]$AAChange,
                          Passage=p,MAF=df[[l]][[p]][[s]]$MAF,
                          Pos=df[[l]][[p]][[s]]$Pos)
      }
      #now add p0
      s.df %<>% add_row(Variant=pre[[1]][[l]]$VcfForm,
                        Gene=pre[[1]][[l]]$Gene,
                        AAChange=pre[[1]][[l]]$AAChange,
                        Passage=0,MAF=pre[[1]][[l]]$MAF,Pos=pre[[1]][[l]]$Pos)
      s.df %<>% add_row(Variant=pre[[2]][[l]]$VcfForm,
                        Gene = pre[[2]][[l]]$Gene,
                        AAChange=pre[[2]][[l]]$AAChange,
                        Passage=0.5,MAF=pre[[2]][[l]]$MAF,Pos=pre[[2]][[l]]$Pos)
      
      s.df %<>% dplyr::mutate(Passage=as.character(Passage)) %>%
        dplyr::mutate(Passage=ifelse(Passage=='0','P0',Passage)) %>%
        dplyr::mutate(Passage=ifelse(Passage=='0.5','P10',Passage)) %>%
        dplyr::mutate(Passage=ifelse(Passage=='1','P13',Passage)) %>%
        dplyr::mutate(Passage=ifelse(Passage=='2','P17',Passage)) %>%
        dplyr::mutate(Passage=ifelse(Passage=='3','P20',Passage)) 
      
      #put in MAF=0 entries for missing variants.
      for (v in unique(s.df$Variant)) {
        all.variants %<>% append(v)
        rows = s.df %>% dplyr::filter(Variant == v)
        if (nrow(rows) < 5) {
          if ('P0' %notin% rows$Passage)
            s.df %<>% add_row(Variant=v,Gene=rows$Gene[[1]],AAChange=rows$AAChange[[1]],Passage='P0',MAF=0,Pos=rows$Pos[[1]])
          if ('P10' %notin% rows$Passage)
            s.df %<>% add_row(Variant=v,Gene=rows$Gene[[1]],AAChange=rows$AAChange[[1]],Passage='P10',MAF=0,Pos=rows$Pos[[1]])
          if ('P13' %notin% rows$Passage)
            s.df %<>% add_row(Variant=v,Gene=rows$Gene[[1]],AAChange=rows$AAChange[[1]],Passage='P13',MAF=0,Pos=rows$Pos[[1]])
          if ('P17' %notin% rows$Passage)
            s.df %<>% add_row(Variant=v,Gene=rows$Gene[[1]],AAChange=rows$AAChange[[1]],Passage='P17',MAF=0,Pos=rows$Pos[[1]])
          if ('P20' %notin% rows$Passage)
            s.df %<>% add_row(Variant=v,Gene=rows$Gene[[1]],AAChange=rows$AAChange[[1]],Passage='P20',MAF=0,Pos=rows$Pos[[1]])
        }
      }
      #now add logical to keep track of variants changing in freq
      for (v in unique(s.df$Variant)) {
        rows = s.df %>% dplyr::filter(Variant == v)
        if (max(rows$MAF)-min(rows$MAF) > 0.4)
          s.df$MAF_Change[which(s.df$Variant==v)] = TRUE
        else s.df$MAF_Change[which(s.df$Variant==v)] = FALSE
      }
      sample.changes[[l]][[s]] = s.df %>% dplyr::arrange(Pos,Passage)
    }
    names(sample.changes[[l]]) = c('S1','S2','S3')
  }
  names(sample.changes) = c('Beta','Delta')
  
  #now add variants missing altogether in either lineage
  for (lin in 1:2) {
    for (s in 1:3) {
      for (var in unique(all.variants)) {
        if (var %notin% sample.changes[[lin]][[s]]$Variant) {
          pos = parse_number(var)
          
          sample.changes[[lin]][[s]] %<>% 
            add_row(Variant=var,Gene=NA,AAChange=NA,Passage='P0',MAF=0,Pos=pos) %>%
            add_row(Variant=var,Gene=NA,AAChange=NA,Passage='P10',MAF=0,Pos=pos) %>%
            add_row(Variant=var,Gene=NA,AAChange=NA,Passage='P13',MAF=0,Pos=pos) %>%
            add_row(Variant=var,Gene=NA,AAChange=NA,Passage='P17',MAF=0,Pos=pos) %>%
            add_row(Variant=var,Gene=NA,AAChange=NA,Passage='P20',MAF=0,Pos=pos) 
        }
      }
      sample.changes[[lin]][[s]] %<>% dplyr::arrange(Pos,Passage)
    }
  }
  return(sample.changes)
}
plotSingleSampleMAFChange = function(df,vars) {
  fig.df = data.frame(Variant=character(),Lineage=numeric(),
                      Passage=character(),Sample=numeric(),
                      MAF=numeric())
  for (l in 1:2) {
    for (s in 1:3) {
      fig.df %<>% add_row(Variant=df[[l]][[s]]$Variant,Lineage=l,
                          Passage=df[[l]][[s]]$Passage,Sample=s,
                          MAF=df[[l]][[s]]$MAF)
    }
  }
  fig.df %<>% dplyr::filter(Variant %in% vars) %>%
    dplyr::mutate(Sample=as.character(Sample)) %>%
    dplyr::mutate(Lineage=ifelse(Lineage=='1','Beta','Delta')) %>%
    dplyr::mutate(LabelNumber=NA)
  
  #remove busy labels when not needed
  for (lin in c('Beta','Delta')) {
    for (pass in c('P0','P10','P13','P17','P20')) {
      for (var in vars) {
        tmp = fig.df %>% dplyr::filter(Lineage==lin,Passage==pass,Variant==var)
        if (length(unique(tmp$MAF)) > 1)
          fig.df[which(fig.df$Lineage==lin & fig.df$Passage==pass & fig.df$Variant==var),'AllSameMAFPass'] = tmp$Sample
      }
    }
  }
  fig.df$Variant = factor(fig.df$Variant,levels=vars)
  
  print(ggplot(fig.df,aes(x=Passage,y=MAF,color=Sample,shape=Sample,label=AllSameMAFPass,group=Sample)) + 
    geom_point(size=3) + theme_bw() + 
    theme(text=element_text(size=16)) + ylim(c(0,1)) +
    scale_color_brewer(palette='Dark2') +
    geom_label_repel(show.legend=F) + geom_line() +
      
    facet_grid(rows=vars(Lineage),cols=vars(Variant)))
  return(fig.df)
}
getChangingVariants = function(data) {
  out.df = data.frame(Variant=character(),Gene=character(),
                      AAChange=character(),Position=numeric())
  for (lin in names(data)) {
    for (s in names(data[[1]])) {
      df = data[[lin]][[s]] %>% dplyr::filter(MAF_Change == TRUE)
      out.df %<>% add_row(Variant=df$Variant,Gene=df$Gene,AAChange=df$AAChange,Position=df$Pos) %>% distinct()
    }
  }
  return(out.df %>% dplyr::arrange(Position))
}
