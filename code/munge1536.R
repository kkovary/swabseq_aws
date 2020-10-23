library(gdata)
library(argparser)
library(xlsx)
# library(openxlsx)

p <- arg_parser("script to output table of barcodes in 96 well format")
p <- add_argument(p,"--rundir",  default=".", help="path to file directory")
# p <- add_argument(p,"--plate",  default=".", help="plate that primers came from")
# p <- add_argument(p,"--quadrant",  default=".", help="which quadrant")

args <- parse_args(p) 
# setwd(args$rundir)

sampleXLS <- list.files()[grepl("SwabSeq.*xlsx",list.files())][1]
SwabSeq <- read.xls(sampleXLS,
                    sheet = 1)

SwabSeq$Plate <- as.character(SwabSeq$Plate)
SwabSeq$Plate_384_Quadrant <- as.character(SwabSeq$Plate_384_Quadrant)

primers <- list(
  's2_r'=read.xls('../../misc/idt_1536.xlsx', sheet=1, header=F, stringsAsFactors=F),
  's2_f'=read.xls('../../misc/idt_1536.xlsx', sheet=2, header=F, stringsAsFactors=F),
  'rpp30_r'=read.xls('../../misc/idt_1536.xlsx', sheet=3, header=F, stringsAsFactors=F),
  'rpp30_f'=read.xls('../../misc/idt_1536.xlsx', sheet=4, header=F, stringsAsFactors=F)
  
  # 's2_r'=openxlsx::read.xlsx('../../misc/idt_1536.xlsx', sheet=1, colNames = FALSE),
  # 's2_f'=openxlsx::read.xlsx('../../misc/idt_1536.xlsx', sheet=2, colNames = FALSE),
  # 'rpp30_r'=openxlsx::read.xlsx('../../misc/idt_1536.xlsx', sheet=3, colNames = FALSE),
  # 'rpp30_f'=openxlsx::read.xlsx('../../misc/idt_1536.xlsx', sheet=4, colNames = FALSE)
  )


plates=c(rep(1,384),rep(2,384), rep(3,384), rep(4,384))

sonly=lapply(primers, function(x) gsub('^.*_', '', x[,1]))

sonlyp=lapply(sonly, function(x) split(x, plates))


asgrids=lapply(sonlyp, function(x) lapply(x, function(y) t(matrix(y, 24,16))))

asgrids2=lapply(asgrids, function(y) lapply(y ,function(x) cbind(c('',toupper(letters[1:16])), rbind(seq(1,24), x))))

asgrids3=lapply(asgrids2, function(x) do.call('rbind', x))

# for(x in names(asgrids3)){
#   write.table(asgrids3[[x]], paste0('~/Downloads/', x, '.tsv'), 
#               sep='\t', quote=F, row.names=F, col.names=F)
# }




# Function to go from 384 to 96 well format
bcPlate_map <- function( 
  amp = "s2", 
  dir = "f", 
  plate = 1, 
  quadrant = "A"){
  
  positions <- switch(
    quadrant,
    "A" = list(cols = seq(from = 1, to = 24, by = 2),
               rows = seq(from = 1, to = 16, by = 2)),
    "B" = list(cols = seq(from = 2, to = 24, by = 2),
               rows = seq(from = 1, to = 16, by = 2)),
    "C" = list(cols = seq(from = 1, to = 24, by = 2),
               rows = seq(from = 2, to = 16, by = 2)),
    "D" = list(cols = seq(from = 2, to = 24, by = 2),
               rows = seq(from = 2, to = 16, by = 2))
  )
  
  positions <- expand.grid(positions$rows, positions$cols)
  colnames(positions) <- c("row","col")
  
  plate_384 <- asgrids[paste0(amp,"_",dir)][[1]][[plate]]
  
  plate_96 <- matrix(data = sapply(1:96, function(x) plate_384[positions[x,1], positions[x,2]]),
                     nrow = 8,
                     ncol = 12)
  return(plate_96)
}


######################
# List of i7 Primers #
######################
i7 <- list()
for(i in 1:nrow(SwabSeq)){
  # i7[i] <- SwabSeq$Plate[i]
  i7[[i]] <- bcPlate_map(dir = "r",
                         plate = SwabSeq$Plate_384_Number[i],
                         quadrant = SwabSeq$Plate_384_Quadrant[i])
}


for(i in seq_along(i7)){
  if(i == 1){
    i7_out <- rbind(as.character(c(1:12)), i7[[i]], as.character(rep(NA, 12)))
    i7_out <- cbind(c(SwabSeq$Plate[i], LETTERS[1:8], NA),
                    i7_out)
  } else{
    i7_temp <- rbind(as.character(c(1:12)), i7[[i]], as.character(rep(NA, 12)))
    i7_temp <- cbind(c(SwabSeq$Plate[i], LETTERS[1:8], NA),
                     i7_temp)
    i7_out <- rbind(i7_out, i7_temp)
  }
}

######################
# List of i5 Primers #
######################
i5 <- list()
for(i in 1:nrow(SwabSeq)){
  # i7[i] <- SwabSeq$Plate[i]
  i5[[i]] <- bcPlate_map(dir = "f",
                         plate = SwabSeq$Plate_384_Number[i],
                         quadrant = SwabSeq$Plate_384_Quadrant[i])
}


for(i in seq_along(i5)){
  if(i == 1){
    i5_out <- rbind(as.character(c(1:12)), i5[[i]], as.character(rep(NA, 12)))
    i5_out <- cbind(c(SwabSeq$Plate[i], LETTERS[1:8], NA),
                    i5_out)
  } else{
    i5_temp <- rbind(as.character(c(1:12)), i5[[i]], as.character(rep(NA, 12)))
    i5_temp <- cbind(c(SwabSeq$Plate[i], LETTERS[1:8], NA),
                     i5_temp)
    i5_out <- rbind(i5_out, i5_temp)
  }
}

#######################################
# Append primer lists to SwabSeq.xlsx #
#######################################
# wb <- openxlsx::loadWorkbook("SwabSeq.xlsx")
# 
# openxlsx::addWorksheet(wb, sheetName="i7")
# openxlsx::addWorksheet(wb, sheetName="i5")
# 
# openxlsx::writeData(wb,"i7",colNames = FALSE,i7_out)
# openxlsx::writeData(wb,"i5",colNames = FALSE,i5_out)
# 
# openxlsx::saveWorkbook(wb, "SwabSeq.xlsx", overwrite = TRUE)

# openxlsx::write.xlsx(i5_out,
#            "SwabSeq.xlsx",
#            sheetName="i5",
#            col.names = FALSE,
#            row.names = FALSE,
#            showNA = FALSE,
#            append=TRUE
# )

xlsx::write.xlsx(i7_out, 
           sampleXLS,
           sheetName="i7", 
           col.names = FALSE,
           row.names = FALSE,
           showNA = FALSE,
           append=TRUE
)

xlsx::write.xlsx(i5_out, 
           sampleXLS,
           sheetName="i5", 
           col.names = FALSE,
           row.names = FALSE,
           showNA = FALSE,
           append=TRUE
)

# 
# # i7 / index1
# p1_i7 <- bcPlate_map(dir = "r", plate = 2, quadrant = "q1")
# p2_i7 <- bcPlate_map(dir = "r", plate = 2, quadrant = "q2")
# p3_i7 <- bcPlate_map(dir = "r", plate = 2, quadrant = "q3")
# 
# # i5 / index2
# p1_i5 <- bcPlate_map(dir = "f", plate = 2, quadrant = "q1")
# p2_i5 <- bcPlate_map(dir = "f", plate = 2, quadrant = "q2")
# p3_i5 <- bcPlate_map(dir = "f", plate = 2, quadrant = "q3")
# 
# write_csv(as.data.frame(p3_i7), "~/Downloads/x.csv")
