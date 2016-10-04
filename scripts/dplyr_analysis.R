library(tidyverse)

# bring into environment
L_GAM <- read.table(file = "~/../Desktop/N123_data/EN51_N123primersonly_GAMKL_unique_final.txt")
FR_GAM <- read.table("~/../Desktop/N123_data/EN50_N123_FR1-N123primers_GAMKL_unique_final.txt")
LP_GAM <- read.table("~/../Desktop/N123_data/EN52_N123_LP_GAMKL_unique_final.txt")
L_GA <- read.table("~/../Desktop/N123_data/EN53_N123_FR1-N123primers_GAKL_unique_final.txt")
FR_GA <- read.table("~/../Desktop/N123_data/EN54_N123primersonly_GAKL_unique_final.txt")
LP_GA <- read.table("~/../Desktop/N123_data/EN55_N123_LP_GAKL_unique_final.txt")

# save as list
N123_libraries <- list('L_GAM' = L_GAM, 'FR_GAM' = FR_GAM,
                       'LP_GAM' = LP_GAM, 'L_GA' = L_GA,
                       'FR_GA' = FR_GA, 'LP_GA' = LP_GA)


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

# you can get the number of IGHV1-18*01 in each library:
IGHV1_18_counts <- as.numeric(sapply(HV_analysis,
                                     function(x) x[x$HV == 'IGHV1-18*01', 2]))

# you can get the proportion of IGHV1-18*01 (as a %) like this:
IGHV1_18_props <- IGHV1_18_counts /
  as.numeric(sapply(HV_analysis, function(x) sum(x$count))) * 100

# plot these proportions
plot(IGHV1_18_props, ylim = c(0,10), xaxt = 'n',
     xlab = '', type = 'h', ylab = 'Percentage',
     main = 'Proportion of IGHV1-18*01 found in libraries')
axis(1, at = seq(1,6), labels = names(N123_libraries))


############### ANALYZE H3LEN #####################

# check H3len for IGHV1-18*01
AnalyzeH3len <- function(x){
  out <- x %>%
    filter(HV == 'IGHV1-18*01') %>%
    arrange(desc(H3len)) %>%
    select(H3len)
  return(as.data.frame(out))}

# apply to get H3len for all 
H3_len_values <- lapply(N123_libraries, AnalyzeH3len)

# change plot window to 2x3
par(mfrow = c(2, 3))

# plot histograms of H3len values for all 6 libraries
for(i in seq(length(H3_len_values))){
  hist(H3_len_values[[i]]$H3len, breaks = 15, xlab = '', ylab = '',
       main = names(H3_len_values)[i], xlim = c(0,30))
  # pause before each plot
  Sys.sleep(1)}

# return graph window back to single plot
par(mfrow = c(1, 1))

# what are the average H3len values for IGHV1-18*01 across libraries?
avg_H3len <- sapply(H3_len_values, function(x) mean(x$H3len))

# you can plot it like this but they're all very similar
plot(avg_H3len, main = 'Average H3len values for IGHV1-18*01',
     xaxt = 'n', xlab = '', ylab = 'Avg H3len',
     ylim = c(0,18))
axis(1, at = seq(1,6), labels = names(N123_libraries))

############# TOP HVs FROM LP_GAM ###################

# get names of top 10 HVs from LP_GAM
top_HVs <- as.data.frame(HV_analysis$LP_GAM[1:10,1])

# function to extract one HV from data set
ExtractFreq <- function(HV, data){
  freq <- data[data$HV == HV, 1]
  return(freq)
}

# function to extract all HVs from data set
extract_all_freqs <- function(dat){
  sapply(top_HVs, extract_freq, data = dat)
}

# extract all top HVs from all data sets
top_HV_all <- lapply(HV_values, extract_all_freqs)
