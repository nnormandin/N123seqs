---
title: "N123 Analysis"
author: "Nick & Erica"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r, include = FALSE}
library(dplyr)
```

### Pre-Processing

The first step is to bring the six libraries in to the environment using `read.table`.

```{r, echo = FALSE}
L_GAM <- read.table(file = "~/../Desktop/N123_data/EN51_N123primersonly_GAMKL_unique_final.txt")
FR_GAM <- read.table("~/../Desktop/N123_data/EN50_N123_FR1-N123primers_GAMKL_unique_final.txt")
LP_GAM <- read.table("~/../Desktop/N123_data/EN52_N123_LP_GAMKL_unique_final.txt")
L_GA <- read.table("~/../Desktop/N123_data/EN53_N123_FR1-N123primers_GAKL_unique_final.txt")
FR_GA <- read.table("~/../Desktop/N123_data/EN54_N123primersonly_GAKL_unique_final.txt")
LP_GA <- read.table("~/../Desktop/N123_data/EN55_N123_LP_GAKL_unique_final.txt")
```

We can then save them as a list so that we can iterate through the libraries using custom functions.

```{r}
N123_libraries <- list('L_GAM' = L_GAM, 'FR_GAM' = FR_GAM,
                       'LP_GAM' = LP_GAM, 'L_GA' = L_GA,
                       'FR_GA' = FR_GA, 'LP_GA' = LP_GA)
```

Since our data didn't come with attribute names embedded in the files, we will have to write a function and apply it across the libraries to manually add variable names.

```{r}
# function to rename columns
RenameCols <- function(x){
  colnames(x) <- c("Reads", "H3nts", "L3nts", "H3aas", "L3aas", "HV",
                   "HJ", "HD", "H3len", "LV", "LJ", "L3len", "Hiso", "Liso")
  return(as.data.frame(x)) # return data frame ensures variable types intact
}

# apply column renaming across all libraries
N123_libraries <- lapply(N123_libraries, RenameCols)
```


### Analyzing HV






