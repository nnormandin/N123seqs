# analysis script

# bring into environment
L_GAM <- read.table("~/N123seqs/data/EN51_N123primersonly_GAMKL_unique_final.txt")

# examine size
dim(L_GAM)

# look at it
head(L_GAM, 10)

# rename columns
colnames(L_GAM) <-c("Reads","H3nts","L3nts","H3aas","L3aas","HV",
                    "HJ","HD","H3len","LV","LJ","L3len","Hiso","Liso")

# check column types
sapply(L_GAM, class)

# plot HV
plot(L_GAM$HV)

# make HV a data table
HV <- L_GAM$HV

# check counts of HV and sort by descending
HV_levels <- summary(HV)
sort(HV_levels, decreasing = TRUE)


