library(tidyverse)

# bring into environment
L_GAM <- read.table("./data/EN51_N123primersonly_GAMKL_unique_final.txt")
FR_GAM <- read.table("./data/EN50_N123_FR1-N123primers_GAMKL_unique_final.txt")
LP_GAM <- read.table("./data/EN52_N123_LP_GAMKL_unique_final.txt")
L_GA <- read.table("./data/EN54_N123primersonly_GAKL_unique_final.txt")
FR_GA <- read.table("./data/EN53_N123_FR1-N123primers_GAKL_unique_final.txt")
LP_GA <- read.table("./data/EN55_N123_LP_GAKL_unique_final.txt")

# save as list
N123_libraries <- list('L_GAM' = L_GAM, 'FR_GAM' = FR_GAM,
                       'LP_GAM' = LP_GAM, 'L_GA' = L_GA,
                       'FR_GA' = FR_GA, 'LP_GA' = LP_GA)

rm(list(L_GAM, FR_GAM, LP_GAM, L_GA, FR_GA, LP_GA))

# function to rename columns
RenameCols <- function(x){
  colnames(x) <- c("Reads", "H3nts", "L3nts", "H3aas", "L3aas", "HV",
                   "HJ", "HD", "H3len", "LV", "LJ", "L3len", "Hiso", "Liso")
  return(as.data.frame(x)) # return data frame ensures variable types intact
}

# apply column renaming across all libraries
N123_libraries <- lapply(N123_libraries, RenameCols)


#################### ANALYZE HV #####################

# function to sort HV count by descending frequency 
AnalyzeHV <- function(x){
  out <- x %>%
    group_by(HV) %>%
    summarise(
      count = length(HV))%>%
    arrange(desc(count))
  return(out)
}

# apply HV count and sort across all libraries and save
HV_analysis <- lapply(N123_libraries, AnalyzeHV)

# you can now access HV counts like this:
HV_analysis$L_GAM

# you can get the number of IGHV1-2*02 in each library:
IGHV1_2_counts <- as.numeric(sapply(HV_analysis,
                                     function(x) x[x$HV == 'IGHV1-2*02', 2]))

# you can get the proportion of IGHV1-2*02 (as a %) like this:
IGHV1_2_props <- IGHV1_2_counts /
  as.numeric(sapply(HV_analysis, function(x) sum(x$count))) * 100

# plot these proportions
plot(IGHV1_2_props, ylim = c(0,10), xaxt = 'n',
     xlab = '', type = 'h', ylab = 'Percentage',
     main = 'Proportion of IGHV1-2*02 found in libraries')
axis(1, at = seq(1,6), labels = names(N123_libraries),las=2)


############### ANALYZE H3LEN #####################

# check H3len for IGHV1-2*02
AnalyzeH3len <- function(x){
  out <- x %>%
    filter(HV == 'IGHV1-2*02') %>%
    arrange(desc(H3len)) %>%
    select(H3len)
  return(as.data.frame(out))}

# apply to get H3len for all 
H3_len_values <- lapply(N123_libraries, AnalyzeH3len)

# change plot window to 2x3
pdf(file = 'h3len.pdf', width = 6)
par(mfrow = c(2, 3))

# plot histograms of H3len values for all 6 libraries
for(i in seq(length(H3_len_values))){
  hist(H3_len_values[[i]]$H3len, breaks = 15, xlab = '', ylab = '',
       main = names(H3_len_values)[i], xlim = c(0,30))
  # pause before each plot
  Sys.sleep(1)}
dev.off()

# return graph window back to single plot
par(mfrow = c(1, 1))

# what are the average H3len values for IGHV1-2*02 across libraries?
avg_H3len <- sapply(H3_len_values, function(x) mean(x$H3len))

# you can plot it like this but they're all very similar
plot(avg_H3len, main = 'Average H3len values for IGHV1-2*02',
     xaxt = 'n', xlab = '', ylab = 'Avg H3len',
     ylim = c(0,18))
axis(1, at = seq(1,6), labels = names(N123_libraries),las=2)

############# TOP HVs FROM LP_GAM ###################

# get names of top 10 HVs from LP_GAM
# this is really messy; look for a simpler way
top_HVs <- as.character(unlist(as.data.frame(HV_analysis$LP_GAM[1:15,1])))

# write function to extract frequency of all HVs in a list from libraries
ExtractHVFreqs <- function(x, HV_list){

  HVFreq <- function(x, HV_name){
    out1 <- as.data.frame(filter(x, HV == HV_name))
    return(nrow(out1))}
  
  out <- sapply(HV_list, HVFreq, x= x)
  return(out)
}

# apply across libraries and bind into data frame
top_HV_freqs <- as.data.frame(do.call(cbind,
                        lapply(N123_libraries,
                               ExtractHVFreqs, HV_list = top_HVs)))

# divide by size of library to get proportion
top_HV_props <- top_HV_freqs /
  sapply(N123_libraries, function(x) nrow(x)) * 100

# generate plot
plot(top_HV_props$L_GAM, xaxt = 'n', ylab = 'Percentage', xlab = '',
     cex.lab = .7, cex.axis = .7, type = 'b',
     pch = 0, ylim = c(0, max(top_HV_props))) # square
axis(side = 1, at = seq(1, 15), labels = top_HVs, cex.axis = .52,las=2)
lines(top_HV_props$FR_GAM, type = 'b', pch = 1) # circle
lines(top_HV_props$LP_GAM, type = 'b', pch = 2) # triangle
lines(top_HV_props$L_GA, type = 'b', pch = 3) # plus
lines(top_HV_props$FR_GA, type = 'b', pch = 4) # x
lines(top_HV_props$LP_GA, type = 'b', pch = 5) # down triangle
legend("topright", inset=.05, title="Library",
       names(N123_libraries), pch = seq(0,5), horiz=FALSE,
       cex = 0.7)
