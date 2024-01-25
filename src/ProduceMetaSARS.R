library(writexl) ; library(readxl)
library(tidyverse)
library(magrittr)
library(glue)

write_biosample_attributes_illumina = function() {
  setwd("~/Documents/PhD Work/SARS CoV 2 Evolution/ZenodoFiles/")
  files_run1 = list.files("Run1")
  names = gsub(".fastq.gz", "", files_run1)
  df = data.frame(sample_name=names,sample_title=c(1:96),Organism="SARS-CoV-2",
                  collected_by="Univ of Laval",collection_date="2023",
                  geographic_location="Canada",host="mus musculus",
                  host_disease="COVID-19",isolate="Lung",
                  isolation_source="BSL3",tmp=c(1:96),run="illumina_1")
  files_run2 = list.files("Run2")
  names = gsub(".fastq.gz", "", files_run2)
  df %<>% add_row(sample_name=names,sample_title=c(97:184),Organism="SARS-CoV-2",
                  collected_by="Univ of Laval",collection_date="2023",
                  geographic_location="Canada",host="mus musculus",
                  host_disease="COVID-19",isolate="Lung",
                  isolation_source="BSL3",tmp=c(97:184),run="illumina_2")
  
  # Now ensure 4 files for each sample per run are clearly connected
  df_comp = df %>% mutate(File1=NA,File2=NA,File3=NA,File4=NA)
  subs <- sub("_L00.*", "", df_comp$sample_name) %>% unique()
  for (sub in subs) {
    rows = which(str_detect(df_comp$sample_name,sub))
    s_names = df_comp$sample_name[rows]

    df_comp$File1[rows] = s_names[[1]]
    df_comp$File2[rows] = s_names[[2]]
    if (length(rows)>2) df_comp$File3[rows] = s_names[[3]]
    if (length(rows)>2) df_comp$File4[rows] = s_names[[4]]
  }
  
  write_xlsx(df_comp,"biosample_attributes_Illumina.xlsx")
}
write_meta_data_illumina = function() {
  template_file = read_xlsx("biosample_attributes_Illumina.xlsx")
  meta_df = data.frame(sample_name = template_file$sample_name,
                       library_ID = glue("Nextera DNA Flex {template_file$tmp}"),
                       title = template_file$sample_name,
                       library_strategy = "AMPLICON",
                       library_source = "VIRAL RNA",
                       library_selection = "RT-PCR",
                       library_layout = "paired",
                       platform = "ILLUMINA",
                       instrument_model = "Illumina MiSeq",
                       design_description = "Step a: dx.doi.org/10.17504/protocols.io.bjgekjte. Step b: dx.doi.org/10.17504/protocols.io.ewov18e4ygr2/v2. Step c: dx.doi.org/10.17504/protocols.io.bjgnkjve.",
                       filetype = "fastq",
                       filename = paste0(template_file$sample_name,".fastq.gz"),
                       file1 = template_file$File1,
                       file2 = template_file$File2,
                       file3 = template_file$File3,
                       file4 = template_file$File4)
  write_xlsx(list(SRA_data = meta_df),"sra_metadata_Illumina.xlsx")
}

############
write_biosample_attributes_nanopore = function() {
  setwd("~/Documents/PhD Work/SARS CoV 2 Evolution/ZenodoFiles/")
  files_run1 = list.files("NanoporeData/Nanopore_Run1")
  names = gsub(".fastq.gz", "", files_run1)
  df = data.frame(sample_name=names,sample_title=c(1:22),Organism="SARS-CoV-2",
                  collected_by="Univ of Laval",collection_date="2023",
                  geographic_location="Canada",host="mus musculus",
                  host_disease="COVID-19",isolate="Lung",
                  isolation_source="BSL3",tmp=c(1:22),run="nanopore_1")
  files_run2 = list.files("NanoporeData/Nanopore_Run2")
  names = gsub(".fastq.gz", "", files_run2)
  df %<>% add_row(sample_name=names,sample_title=c(23:44),Organism="SARS-CoV-2",
                  collected_by="Univ of Laval",collection_date="2023",
                  geographic_location="Canada",host="mus musculus",
                  host_disease="COVID-19",isolate="Lung",
                  isolation_source="BSL3",tmp=c(23:44),run="nanopore_2")
  
  # Now ensure 4 files for each sample per run are clearly connected
  df_comp = df 
  
  write_xlsx(df_comp,"biosample_attributes_nanopore.xlsx")
}
write_meta_data_nanopore = function() {
  template_file = read_xlsx("biosample_attributes_nanopore.xlsx")
  meta_df = data.frame(sample_name = template_file$sample_name,
                       library_ID = glue("Libprep {template_file$tmp}"),
                       title = template_file$sample_name,
                       library_strategy = "AMPLICON",
                       library_source = "VIRAL RNA",
                       library_selection = "RT-PCR",
                       library_layout = "single",
                       platform = "OXFORD_NANOPORE",
                       instrument_model = "PromethION",
                       design_description = "Nanopore library preparation was made following Reiling et al. 2020 and libraries were sequenced on the PromethION 24 sequencer with PromethION Flow Cells V.9.4.1 for a total of 10 million reads",
                       filetype = "fastq",
                       filename = paste0(template_file$sample_name,".gz"))
  write_xlsx(list(SRA_data = meta_df),"sra_metadata_nanopore.xlsx")
}