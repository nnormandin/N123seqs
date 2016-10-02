# analysis script

# bring into environment
L_GAM <- read.table("~/N123seqs/data/EN51_N123primersonly_GAMKL_unique_final.txt")
FR_GAM <- read.table("~/N123seqs/data/EN50_N123_FR1-N123primers_GAMKL_unique_final.txt")
LP_GAM <- read.table("~/N123seqs/data/EN52_N123_LP_GAMKL_unique_final.txt")
L_GA <- read.table("~/N123seqs/data/EN53_N123_FR1-N123primers_GAKL_unique_final.txt")
FR_GA <- read.table("~/N123seqs/data/EN54_N123primersonly_GAKL_unique_final.txt")
LP_GA <- read.table("~/N123seqs/data/EN55_N123_LP_GAKL_unique_final.txt")


# examine size
dim(L_GAM)

# look at it
head(L_GAM, 10)

# rename columns
N123_libraries <- list(L_GAM, FR_GAM, LP_GAM, L_GA, FR_GA, LP_GA)

RenameCols <- function(x){
  colnames(x) <- c("Reads", "H3nts", "L3nts", "H3aas", "L3aas", "HV",
                          "HJ", "HD", "H3len", "LV", "LJ", "L3len", "Hiso", "Liso")
}

N123_libraries_test <- lapply(N123_libraries, RenameCols)

# check column types
sapply(L_GAM, class)

# plot HV
plot(L_GAM$HV)

# make HV a data table
HV <- L_GAM$HV

# check counts of HV and sort by descending
HV_levels <- summary(HV)
HV_sorted <- sort(HV_levels, decreasing = TRUE)
plot(HV_sorted, type = 'h')

# check H3len for highest HV level
highest_HV <- L_GAM[L_GAM$HV == 'IGHV1-18*01',]
hist(highest_HV$H3len, breaks = 45)

# check values of H3len
unique(highest_HV$H3len)


# function to analyze HVs
Analyze_HV <- function(x){
  HV_levels <- summary(x$HV)
  HV_sorted <- as.data.frame(sort(HV_levels, decreasing = TRUE))
  HV_sorted$HV <- rownames(HV_sorted)
  rownames(HV_sorted) <- NULL
  colnames(HV_sorted) <- c('Number', 'HV')
  HV_sorted$Number <- as.numeric(HV_sorted$Number)
  total <- sum(HV_sorted$Number)
  HV_sorted$Number <- (HV_sorted$Number/total)*100
  return(HV_sorted)
}

# apply to all data sets
HV_values <- lapply(N123_libraries, Analyze_HV)

merged <- do.call(merge.data.frame, HV_values)


# analyze top HVs across data sets
top_HVs <- HV_values[[3]][1:10,2]

# function to extract one HV from data set
extract_freq <- function(HV, data){
  freq <- data[data$HV == HV, 1]
  return(freq)
}

# function to extract all HVs from data set
extract_all_freqs <- function(dat){
  sapply(top_HVs, extract_freq, data = dat)
}

# extract all top HVs from all data sets
top_HV_all <- lapply(HV_values, extract_all_freqs)

# bind results
top_HV_all <- do.call(cbind, top_HV_all)

# remove NA value
top_HV_all[4,1] <- 0

