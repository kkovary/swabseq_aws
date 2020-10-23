#usage: countAmplicons.R [--] [--help] [--opts OPTS] [--rundir RUNDIR]
#       [--basespaceID BASESPACEID] [--countsOnly COUNTSONLY]
#       [--extendedAmplicons EXTENDEDAMPLICONS]
#
#utility to count amplicons for SwabSeq
#
#flags:
#  -h, --help               show this help message and exit
#
#optional arguments:
#  -x, --opts               RDS file containing argument values
#  -r, --rundir             path containing SampleSheet [default: .]
#  -b, --basespaceID        BaseSpace Run ID [default: 1000]
#  -c, --countsOnly         only output table of counts [default: TRUE]
#  -e, --extendedAmplicons  additional swabseq amplicons [default:
#                           FALSE]


#parse command line arguments and exit fast if run with -h
library(argparser)
p = arg_parser("utility to count amplicons for SwabSeq")
p=add_argument(p,"--rundir",  default=".", help="path containing SampleSheet")
# p=add_argument(p,"--outdir", default = ".", help = "path to data output")
p=add_argument(p,"--basespaceID",  default=NULL, help="BaseSpace Run ID")
#p=add_argument(p,"--countsOnly",  default=TRUE, help="only output table of counts")
p=add_argument(p,"--extendedAmplicons", default=F, help="additional swabseq amplicons")
p=add_argument(p,"--bcSwap", default=F, help="Analysis of index mis-assignment")
p=add_argument(p,"--lbuffer", default=30000001, help="how many reads to load into ram at a time")
p=add_argument(p,"--threads", default=1, help="number of threads for bcl2fastq & amatch")
p <- add_argument(p, "--allUDI", default = TRUE, help = "measure all UDIs")
# p <- add_argument(p, "--qc", default = FALSE, help = "output qc report")
p <- add_argument(p, "--git", default = TRUE, help = "commit to GitHub")
args=parse_args(p) 

#load required packages
library(ShortRead)
library(stringdist)
library(tidyverse)
library(Rqc)
library(savR)
# library(openxlsx)
#for example ... 
#rundir='/data/Covid/swabseq/runs/SwabSeqV17/'
#basespaceID="195853687"
rundir=args$rundir
basespaceID=args$basespaceID
outputCountsOnly=args$countsOnly
extendedAmplicons=args$extendedAmplicons
runSwap = args$bcSwap
threads = args$threads
nbuffer = as.numeric(args$lbuffer)
#requires BCLs (or basespaceID), python3.7 in path, bcl2fastq in path and SampleSheet.csv (or xls)



# setwd(rundir)
if (file.exists(rundir)){
  setwd(file.path(mainDir, rundir))
} else {
  dir.create(file.path(rundir))
  setwd(file.path(rundir))
  
}

#-----------------------------------------------------------------------------------------------------

# if fastqs don't exist grab them from basespace
fastqR1  <- paste0(rundir, 'out/Undetermined_S0_R1_001.fastq.gz')
if(!file.exists(fastqR1)) {
  #Pull BCLs from basespace [skip this section if you already placed bcls in rundir/bcls/] ------------
  #if running miseq then paste run id here
  #if miseq run then grab from basespace, otherwise place bcls here and skip lines 12-23
  
  if(is.null(basespaceID)){
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

# read in n reads at a time (reduce if RAM limited)
# nbuffer=3e7
# nbuffer=as.numeric(args$lbuffer)
#print(nbuffer)


# error correct the indices and count amplicons
errorCorrectIdxAndCountAmplicons=function(rid, count.table, ind1,ind2,e=1){
  # get set of unique expected index1 and index2 sequences
  index1=unique(count.table$index)
  index2=unique(count.table$index2)
  # for subset of reads where matching amplicon of interest (rid)
  # match observed index sequences to expected sequences allowing for e hamming distance
  i1m=amatch(ind1[rid],index1, method='hamming', maxDist=e, matchNA=F, nthread = 6)
  i2m=amatch(ind2[rid],index2, method='hamming', maxDist=e, matchNA=F, nthread = 6)
  # combine the error corrected indices together per read
  idm=paste0(index1[i1m], index2[i2m])
  #match error corrected indices to lookup table and count
  tS2=table(match(idm, count.table$mergedIndex))
  #get index in lookup table for indices with at least one observed count
  tbix=match(as.numeric(names(tS2)), 1:nrow(count.table))
  #increment these samples by count per sample
  count.table$Count[tbix]=as.vector(tS2)+count.table$Count[tbix]
  
  
  # Index1 not aligned
  i1_na <- table(ind1[rid][is.na(i1m)])
  
  # Index 2 not aligned
  i2_na <- table(ind2[rid][is.na(i2m)])
  
  # Combine results into a list
  ind_na <- list(i1_na = i1_na,
                 i2_na = i2_na)
  
  return(list(count.table = count.table,
              ind_na = ind_na))
}

# Barcode swap function
bcSwap=function(rid, count.table, ind_df, ind1,ind2,e=1){
  # get set of unique expected index1 and index2 sequences
  index1=unique(count.table$index)
  index2=unique(count.table$index2)
  # for subset of reads where matching amplicon of interest (rid)
  # match observed index sequences to expected sequences allowing for e hamming distance
  i1m=amatch(ind1[rid],index1, method='hamming', maxDist=e, matchNA=F, nthread = threads)
  i2m=amatch(ind2[rid],index2, method='hamming', maxDist=e, matchNA=F, nthread = threads)
  # combine the error corrected indices together per read
  
  
  x <- tibble(ind1 = factor(index1[i1m], levels = index1), 
              ind2 = factor(index2[i2m], levels = index2), 
              idm = paste0(ind1, ind2)) %>% 
    filter(!is.na(ind1),
           !is.na(ind2)) %>% 
    group_by(ind1, ind2) %>% 
    summarise(n = n())
  
  return(bind_rows(ind_df,x))
}

# Amplicon List-----------------------------------------------------------------------
#expected amplicons, note will update for RPP30 spike-in 
#RPP
#CGCAGA gccttcaggtcagaacccgc
#RPP spike
#GCGTCA gccttcaggtcagaacccgc
# amplicons=list(
#     S2=      'TATCTTCAACCTAGGACTTTTCTATT',
#     S2_spike='ATAGAACAACCTAGGACTTTTCTATT',
#     RPP30='CGCAGAGCCTTCAGGTCAGAACCCGC'
# )
# if(extendedAmplicons) {
#     amplicons=list(
#         S2=      'TATCTTCAACCTAGGACTTTTCTATT',
#         S2_spike='ATAGAACAACCTAGGACTTTTCTATT',
#         RPP30='CGCAGAGCCTTCAGGTCAGAACCCGC',
#         RPP30_spike='GCGTCAGCCTTCAGGTCAGAACCCGC'
#     )
# }

# Munginging sample sheet-------------------------------------------------------------------
ss=read.delim(paste0('../../misc/SampleSheet.csv'), stringsAsFactors=F, skip=14, sep=',')
ss$mergedIndex=paste0(ss$index, ss$index2)

if(args$allUDI){
  all_udi <- read_csv("../../misc/bc-map-all.csv")
  
  # filter for merged index
  i5_dir <- all_udi$ins[which(ss$mergedIndex == all_udi$mergedIndex)][1]
  all_udi <- all_udi %>% 
    filter(ins == i5_dir,
           !mergedIndex %in% ss$mergedIndex) %>% 
    select(-ins)
  
  ss <- full_join(ss, all_udi)
}

bc_map <- read_csv("../../misc/barcode-map.csv") %>% 
  filter(bc_set == ss$bc_set[1])
amplicons <- split(bc_map$sequence, bc_map$amplicon)

# this code would be obviated if indices designate wells, for most analyses here there are different indices for s2/s2spike and rpp30
# subset of indices for S2/S2 spike
if(sum(grepl('-1$', ss$Sample_ID))==0){
  ssS=ss
  ssR=ss
} else {
  ssS=ss[grep('-1$', ss$Sample_ID),]
  #subset of indices for RPP30
  ssR=ss[grep('-2$', ss$Sample_ID),]
}
#initalize output count tables ------------------------------------------------------------
temp.table.names <- c()
count.tables <- list()
for (i in seq_along(amplicons)) {
  if (grepl("RPP30", names(amplicons)[i])) {
    temp.table <- ssR
  } else{
    temp.table <- ssS
  }
  temp.table$Count <- 0
  temp.table$amplicon <- names(amplicons)[i]
  
  temp.table.name <- paste0(names(amplicons)[i], ".table")
  temp.table.names[i] <- temp.table.name
  
  count.tables[[i]] <- temp.table
  
  # assign(temp.table.name, temp.table)
}
names(count.tables) <- temp.table.names

# S2.table=ssS; S2_spike.table=ssS; RPP30.table=ssR; RPP30_spike.table=ssR;
# S2.table$Count=0; S2_spike.table$Count=0; RPP30.table$Count=0; RPP30_spike.table$Count=0
# S2.table$amplicon='S2'; S2_spike.table$amplicon='S2_spike'; RPP30.table$amplicon='RPP30'; RPP30_spike.table$amplicon='RPP30_spike'

S2.swap <- tibble(ind1 = factor(), ind2 = factor(), n = numeric())
S2_spike.swap <- tibble(ind1 = factor(), ind2 = factor(), n = numeric())
RPP30.swap <- tibble(ind1 = factor(), ind2 = factor(), n = numeric())
#------------------------------------------------------------------------------------------

fastq_dir  <- paste0(rundir, 'out/')
in.fileI1  <- paste0(fastq_dir, 'Undetermined_S0_I1_001.fastq.gz')
in.fileI2  <- paste0(fastq_dir, 'Undetermined_S0_I2_001.fastq.gz')
in.fileR1  <- paste0(fastq_dir, 'Undetermined_S0_R1_001.fastq.gz')

################################
# Alignment and quantification #
################################

i1 <- FastqStreamer(in.fileI1, nbuffer, readerBlockSize = 1e9, verbose = T)
i2 <- FastqStreamer(in.fileI2, nbuffer, readerBlockSize = 1e9, verbose = T)
r1 <- FastqStreamer(in.fileR1, nbuffer, readerBlockSize = 1e9, verbose = T)

repeat{
  rfq1 <- yield(i1) 
  if(length(rfq1) == 0 ) { break }
  rfq2 <- yield(i2) 
  rfq3 <- yield(r1) 
  ind1 <- sread(rfq1)
  ind2 <- sread(rfq2)
  rd1  <- sread(rfq3)
  
  #match amplicon
  amp.match=lapply(amplicons, function(x) {amatch(rd1, x, method='hamming', maxDist=1, matchNA=F, nthread = threads)})
  
  #summary
  amp.match.summary <- sapply(amp.match, function(x) sum(!is.na(x)))
  
  # no_align <- sum(table(.Internal(unlist(lapply(amp.match, function(x) which(is.na(x))), FALSE, FALSE))) == length(amp.match.summary))
  # amp.match.summary <- c(amp.match.summary, no_align)
  # names(amp.match.summary) <- c(names(amp.match.summary[-length(amp.match.summary)]),"no_align")
  
  # Speed fix from Josh
  amt=do.call('cbind', amp.match)
  amt[is.na(amt)]=0
  amt=rowSums(amt)
  no_align=sum(amt==0)
  amp.match.summary <- c(amp.match.summary, no_align)
  names(amp.match.summary) <- c(names(amp.match.summary[-length(amp.match.summary)]),"no_align")
  
  print(amp.match.summary)
  
  #convert to indices
  am.mat=do.call('cbind', amp.match)
  per.amplicon.row.index=apply(!is.na(am.mat), 2, function(x) which(x==TRUE))
  
  #for each amplicon of interest count up reads where indices match expected samples
  ind_na <- list()#initiate unaligned index counts
  for(i in seq_along(count.tables)){
    count.table <- errorCorrectIdxAndCountAmplicons(per.amplicon.row.index[[i]], 
                                                    count.tables[[i]], 
                                                    ind1,
                                                    ind2, 
                                                    1)
    count.tables[[i]] <- count.table$count.table
    ind_na[[i]] <- count.table$ind_na
  }
  names(ind_na) <- names(count.tables)
  
  
  # if(runSwap){
  #     S2.swap       = bcSwap(per.amplicon.row.index$S2, S2.table, S2.swap, ind1,ind2, 1)
  #     S2_spike.swap = bcSwap(per.amplicon.row.index$S2_spike, S2_spike.table, S2_spike.swap, ind1,ind2,1)
  #     RPP30.swap    = bcSwap(per.amplicon.row.index$RPP30, RPP30.table, RPP30.swap, ind1,ind2,1)
  # }
  
}
close(i1); close(i2); close(r1);
results <- count.tables
# results=list(S2.table=S2.table, S2_spike.table=S2_spike.table, RPP30.table=RPP30.table)

if(extendedAmplicons){ results=list(S2.table=S2.table, S2_spike.table=S2_spike.table, RPP30.table=RPP30.table,RPP30_spike.table=RPP30_spike.table) }

do.call('rbind', results) %>% write_csv(paste0(rundir, 'countTable.csv')) 
saveRDS(results, file=paste0(rundir, 'countTable.RDS'),version=2)

if(runSwap){
  swap <- list(S2.swap = S2.swap, 
               S2_spike.swap = S2_spike.swap,
               RPP30.swap = RPP30.swap)
  saveRDS(swap, file = paste0(rundir, 'swap.RDS'),version=2)
}

##################
# Save QC Report #
##################


# Illumina stats
sav=savR(paste0(rundir))
tMet=tileMetrics(sav)
phiX=mean(tMet$value[tMet$code=='300'])
clusterPF=mean(tMet$value[tMet$code=='103']/tMet$value[tMet$code=='102'], na.rm=T)
clusterDensity=mean(tMet$value[tMet$code=='100']/1000)
clusterDensity_perLane=sapply(split(tMet, tMet$lane), function(x) mean(x$value[x$code=='100']/1000))    
seq.metrics=data.frame("totalReads"=format(sum(amp.match.summary),  big.mark=','),
                       "totalReadsPassedQC"=format(sum(amp.match.summary[!(names(amp.match.summary) %in% 'no_align')]), big.mark=','),
                       "phiX"=paste(round(phiX,2), "%"), "clusterPF"=paste(round(clusterPF*100,1), "%"),
                       "tot_phiX" = format(round(phiX * sum(amp.match.summary)), big.mark = ','),
                       "clusterDensity"=paste(round(clusterDensity,1), 'K/mm^2'), 
                       "clusterDensity_perLane"=paste(sapply(clusterDensity_perLane, round,1),collapse=' '))

# Read stats
qcRes <- rqc(path = fastq_dir, pattern = ".fastq.gz", openBrowser=FALSE)
read_quality <- rqcCycleQualityBoxPlot(qcRes) + ylim(0,NA)
seq_cont_per_cycle <- rqcCycleBaseCallsLinePlot(qcRes)
read_freq_plot <- rqcReadFrequencyPlot(qcRes)
base_calls_plot <- rqcCycleBaseCallsLinePlot(qcRes)

params <- list(
  experiment = strsplit(rundir,"/") %>% unlist() %>% tail(1),
  amp.match.summary = amp.match.summary,
  results = results,
  seq.metrics = seq.metrics,
  ind_na = ind_na,
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

######################################
# Add analysis template to directory #
######################################
if (!file.exists("Analysis.Rmd")) {
  system("cp ../../code/Analysis.Rmd .")
}

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
  pdf_name <- paste0(exp_name,".pdf")
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

# this is moved to mungeTables() in helper_functions.R
#if(outputCountsOnly==FALSE) {
#    results=readRDS(paste0(rundir, 'countTable.RDS'))
#    attach(results)
#    # re-formatting ... output count table 
#    df=rbind(S2.table,S2_spike.table,RPP30.table)
#    df$virus_copy=as.factor(df$virus_copy) 
#    df$Col=as.factor(gsub('^.', '', df$Sample_Well))
#    df$Row=factor(gsub('..$', '', df$Sample_Well), levels=rev(toupper(letters[1:8])))
#    df$Sample=paste0(df$Plate_ID, '-' ,df$Sample_Well)
#    df$Plate_ID=as.factor(df$Plate_ID)
#    df$Plate_ID=factor(df$Plate_ID, levels(df$Plate_ID)[order(as.numeric(gsub('Plate', '', levels(df$Plate_ID))))])  
#    df$Plate_384=as.factor(df$Plate_384)
#    #to do update this for extended amplicon set
#    df$amplicon=factor(df$amplicon, level=c('S2', 'S2_spike', 'RPP30'))
#    df %>%
#      write_csv(paste0(rundir, 'annotated_df.csv')) 
#
#
#    #plate visualization 
#    df %>%
#      ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
#      geom_raster() +
#      coord_equal() +
#      facet_grid(amplicon~Plate_384+Plate_ID) +
#      scale_fill_viridis_c(option = 'plasma')
#    ggsave(paste0(rundir,'plateVis.png'))
#
#
#
#    #assay results
#    dfs= df %>%filter(amplicon=='S2') %>%  
#      count(Sample_Well, wt=Count, name='S2_total_across_all_wells') %>%
#      right_join(df)
#    dfs= df %>%filter(amplicon=='S2'|amplicon=='S2_spike') %>%  
#      count(Sample, wt=Count, name='Stotal') %>%
#      right_join(dfs)
#    dfs= dfs %>% count(Sample, wt=Count, name='well_total') %>%
#      right_join(dfs) %>% 
#      select(-mergedIndex, -Sample_ID, -index, -index2 ) %>% 
#      spread(amplicon, Count) %>% 
#      mutate(S2_normalized_to_S2_spike=(S2+1)/(S2_spike+1))%>%
#      mutate(RPP30_Detected=RPP30>10) %>%  
#      mutate(SARS_COV_2_Detected=S2_normalized_to_S2_spike>.003)
#    dfs$SARS_COV_2_Detected[!dfs$RPP30_Detected]='Inconclusive'
#    dfs$SARS_COV_2_Detected[dfs$Stotal<2000]='Inconclusive'
#    dfs %>%
#      write_csv(paste0(rundir, 'Calls_per_sample.csv')) 
#}
