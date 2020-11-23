

message("Loading Libraries")

suppressMessages(library(argparser))
p = arg_parser("utility to count amplicons for SwabSeq")
p=add_argument(p,"--rundir",  default=".", help="path containing SampleSheet")
p=add_argument(p,"--basespaceID",  default=NA, help="BaseSpace Run ID")
p=add_argument(p,"--threads", default=1, help="number of threads for bcl2fastq & amatch")
p <- add_argument(p, "--git", default = TRUE, help = "commit to GitHub")
args=parse_args(p) 

#load required packages
suppressMessages(library(tidyverse))
suppressMessages(library(Rqc))
suppressMessages(library(savR))
suppressMessages(library(Biostrings))
suppressMessages(library(xlsx))

rundir=args$rundir
basespaceID=args$basespaceID
threads = args$threads



# setwd(rundir)
if (file.exists(rundir)){
  setwd(file.path(rundir))
} else {
  dir.create(file.path(rundir))
  setwd(file.path(rundir))
  
}


#-----------------------------------------------------------------------------------------------------

# if fastqs don't exist grab them from basespace
fastqR1  <- paste0(rundir, 'bcls/out/Undetermined_S0_R1_001.fastq.gz')
if(!file.exists(fastqR1)) {
  #Pull BCLs from basespace [skip this section if you already placed bcls in rundir/bcls/] ------------
  #if running miseq then paste run id here
  #if miseq run then grab from basespace, otherwise place bcls here and skip lines 12-23
  
  if(is.na(basespaceID)){
    basespaceID <- tail(strsplit(rundir,"/")[[1]],1)
    system(paste("bs download run --name", basespaceID, "-o bcls/"))
  } else{
    system(paste("bs download run --id", basespaceID, "-o bcls/"))
  }
  
  # run bcl2fastq to generate fastq.gz files (no demux is happening here)
  setwd(paste0(rundir,'bcls/'))
  #note this is using 64 threads and running on a workstation, reduce threads if necessary
  system(paste("bcl2fastq --runfolder-dir . --output-dir out/ --create-fastq-for-index-reads  --ignore-missing-bcl --use-bases-mask=Y26,I10,I10 --processing-threads", threads, "--no-lane-splitting --sample-sheet /dev/null"))
  #for sharing, make tar of bcls directory
  # system(paste("tar -cvf", paste0(rundir,'bcls.tar'), "../bcls/"))
  #-----------------------------------------------------------------------------------------------------
}

setwd(file.path(rundir))

# Align
system(paste("python3.8 ../../code/dict_align.py --rundir ./ --dictdir ../../hash_tables/"))


###################
# Reformat output #
###################

# Munginging sample sheet-------------------------------------------------------------------
ss=read.delim(paste0('../../misc/SampleSheet.csv'), stringsAsFactors=F, skip=14, sep=',') %>% 
  mutate(mergedIndex = paste0(index, index2))


results <- read_csv("results.csv") %>% 
  dplyr::rename(amplicon = amps,
                Count = `0`,
                index = i1,
                index2 = i2) %>% 
  # left_join(bc_map, by = "sequence") %>% 
  mutate(mergedIndex = paste0(index, index2)) %>% 
  left_join(ss, by = c("mergedIndex","index","index2"))


# Check direction of index 1 before merging w/ 384 plate map
ind1 <- read_tsv("../../hash_tables/ind1.txt", col_names = FALSE) %>% pull(X1)
ind2 <- read_tsv("../../hash_tables/ind2.txt", col_names = FALSE) %>% pull(X1)
pm384 <- read_csv("../../misc/384_plate_map.csv")

if (sum(unique(results$index) %in% pm384$index) <= 10) { # 10 of the indices are RC of each other
  pm384$index <- as.character(reverseComplement(DNAStringSet(pm384$index)))
}
if (sum(unique(results$index2) %in% pm384$index2) <= 12) { # 10 of the indices are RC of each other
  pm384$index2 <- as.character(reverseComplement(DNAStringSet(pm384$index2)))
}

results <- results %>%
  left_join(pm384) %>%
  separate(pm, into = c("pm_384","row_384","col_384")) %>%
  mutate(row_384 = as.numeric(row_384),
         col_384 = as.numeric(col_384))

ss <- ss %>%
  left_join(pm384) %>%
  separate(pm, into = c("pm_384","row_384","col_384")) %>%
  mutate(row_384 = as.numeric(row_384),
         col_384 = as.numeric(col_384))


# Check RC for indexes
ind1_rc <- DNAStringSet(ind1) %>% reverseComplement() %>% as.character()
ind2_rc <- DNAStringSet(ind2) %>% reverseComplement() %>% as.character()
if(sum(results$index %in% ind1) < sum(results$index %in% ind1_rc)){
  ind1 <- ind1_rc
}
if(sum(results$index2 %in% ind2) < sum(results$index2 %in% ind2_rc)){
  ind2 <- ind2_rc
}

results <- results %>% 
  mutate(index = factor(index, levels = ind1),
         index2 = factor(index2, levels = ind2))

write_csv(results, paste0(rundir, 'countTable.csv')) 
saveRDS(results, file=paste0(rundir, 'countTable.RDS'),version=2)

##################
# Save QC Report #
##################


##################
# Save QC Report #
##################

classification <- results %>%
  filter(!is.na(Plate_ID)) %>% 
  right_join(ss) %>% 
  group_by_at(names(.)[!names(.) %in% c("Count", "amplicon")]) %>% 
  summarise(S2_spike = sum(Count[grepl("S2_spike_0",amplicon)], na.rm = TRUE),
            S2 = sum(Count[amplicon == "S2"], na.rm = TRUE),
            RPP30 = sum(Count[amplicon == "RPP30"], na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(S2_spike = ifelse(is.na(S2_spike), 0, S2_spike),
         S2 = ifelse(is.na(S2), 0, S2),
         RPP30 = ifelse(is.na(RPP30), 0, RPP30)) %>%
  mutate(s2_vs_spike = ((S2 + 1) / (S2_spike + 1)),
         classification = NA,
         classification = ifelse(S2 + S2_spike < 500 & RPP30 < 10,
                                 "failed: low S2 & RPP30",
                                 ifelse(S2 + S2_spike < 500 & RPP30 >= 10,
                                        "failed: low S2",
                                        ifelse(S2 + S2_spike >= 500 & RPP30 < 10,
                                               "failed: low RPP30",
                                               ifelse(s2_vs_spike > 0.003 & RPP30 >= 10,
                                                      "COVID_pos",
                                                      ifelse(s2_vs_spike < 0.003 & RPP30 >= 10,
                                                             "COVID_neg",
                                                             classification)))))) %>% 
  dplyr::select(index, index2, pm_384, row_384, col_384, Plate_ID, Sample_Well, S2_spike, S2, RPP30, s2_vs_spike, classification)


amp.match.summary.df <- results %>% 
  group_by(amplicon) %>% 
  summarise(sum = sum(Count)) %>% 
  mutate(amplicon = ifelse(is.na(amplicon),
                           "no_align",
                           amplicon))

amp.match.summary <- amp.match.summary.df$sum
names(amp.match.summary) <- amp.match.summary.df$amplicon

# Illumina stats
sav=savR(rundir)
tMet=tileMetrics(sav)
phiX=mean(tMet$value[tMet$code=='300'])
clusterPF=mean(tMet$value[tMet$code=='103']/tMet$value[tMet$code=='102'], na.rm=T)
clusterDensity=mean(tMet$value[tMet$code=='100']/1000)
clusterDensity_perLane=sapply(split(tMet, tMet$lane), function(x) mean(x$value[x$code=='100']/1000))    
seq.metrics=data.frame("totReads"=format(sum(amp.match.summary),  big.mark=','),
                       "totReadsPassedQC"=format(sum(amp.match.summary[!(names(amp.match.summary) %in% 'no_align')]), big.mark=','),
                       "phiX"=paste(round(phiX,2), "%"), "clusterPF"=paste(round(clusterPF*100,1), "%"),
                       "tot_phiX" = format(round(phiX * sum(amp.match.summary)), big.mark = ','),
                       "clustDensity"=paste(round(clusterDensity,1), 'K/mm^2'), 
                       "clustDensity_perLane"=paste(sapply(clusterDensity_perLane, round,1),collapse=' '))

# Read stats
fastq_dir <- paste0(rundir,"out/")

qcRes <- rqc(path = fastq_dir, pattern = ".fastq.gz", openBrowser=FALSE, workers = 6)
read_quality <- rqcCycleQualityBoxPlot(qcRes) + ylim(0,NA)
seq_cont_per_cycle <- rqcCycleBaseCallsLinePlot(qcRes)
read_freq_plot <- rqcReadFrequencyPlot(qcRes)
base_calls_plot <- rqcCycleBaseCallsLinePlot(qcRes)


params <- list(
  experiment = strsplit(rundir,"/") %>% unlist() %>% tail(1),
  amp.match.summary = amp.match.summary,
  results = results,
  seq.metrics = seq.metrics,
  classification = classification,
  # ind_na = ind_na,
  # qcRes = qcRes,
  read_quality = read_quality,
  seq_cont_per_cycle = seq_cont_per_cycle,
  read_freq_plot = read_freq_plot,
  base_calls_plot = base_calls_plot
)

rmarkdown::render(
  input = "../../code/qc_report.Rmd",
  output_file = paste0(params$experiment,".pdf"),
  output_dir = rundir,
  params = params,
  envir = new.env(parent = globalenv())
)


##################
# Push to GitHub #
##################

if(args$git){
  # Commit countTable.csv
  exp_name <- strsplit(rundir,"/") %>% unlist() %>% tail(1)
  
  # Add countTable.csv to git
  system(paste0("git add ", 
                # rundir, 
                "countTable.csv"))
  
  # Commit countTable.csv to git
  system(paste0("git commit ", 
                # rundir, 
                "countTable.csv", 
                " -m '",
                exp_name,
                " has finished'"))
  
  
  # Add pdf to git
  pdf_name <- paste0(exp_name,".html")
  system(paste0("git add ", 
                # rundir, 
                pdf_name))
  
  
  system(paste0("git commit ", 
                # rundir, 
                pdf_name, 
                " -m '",
                exp_name,
                " results summary'"))
  
  # Add spreadsheets to git
  system(paste0("git add ", 
                "SwabSeq.xlsx"))
  
  system(paste0("git commit ", 
                "SwabSeq.xlsx", 
                " -m '",
                exp_name,
                " SwabSeq.xlsx'"))
  
  system(paste0("git add ", 
                "SampleSheet.csv"))
  
  
  system(paste0("git commit ", 
                "SampleSheet.csv", 
                " -m '",
                exp_name,
                " SampleSheet.csv'"))
  
  # Add Analysis.Rmd to git
  system(paste0("git add ", 
                "Analysis.Rmd"))
  
  system(paste0("git commit ", 
                "Analysis.Rmd", 
                " -m '",
                exp_name,
                " Analysis.Rmd'"))
  
  system("git push")
}