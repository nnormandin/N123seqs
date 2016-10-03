# dplyr method of analysis script
library(tidyverse)

# bring into environment
L_GAM <- read.table("~/../Desktop/N123_data/EN50_N123_FR1-N123primers_GAMKL_unique_final.txt")

# rename columns
colnames(L_GAM) <-c("Reads","H3nts","L3nts","H3aas","L3aas","HV",
                    "HJ","HD","H3len","LV","LJ","L3len","Hiso","Liso")

# dplyr method
HV_Count <- L_GAM %>%
  group_by(HV) %>%
  summarise(
    count = length(HV))%>%
  arrange(desc(count))

# check H3len for highest HV
H3 <- L_GAM %>%
  filter(HV == 'IGHV3-23*01') %>%
  arrange(desc(H3len)) %>%
  select(HV, H3len)

H3_Count <- H3 %>%
  group_by(H3len) %>%
  summarise(
    count = length(H3len)) %>%
  arrange(-desc(count))

# plot frequency of H3len values
plot(x = H3_Count$H3len, y = H3_Count$count)

# table method
HV <- L_GAM$HV
table(HV)

# check values of H3len
unique(H3$H3len)


# method to aggregate HV counts
double_HV <- rbind(HV_Count, HV_Count)
double_HV %>%
  group_by(HV) %>%
  summarise(
    total = sum(count))

# plot?
newplot <- ggplot(as.data.frame(HV_Count), aes(HV, count))
newplot + geom_point()
