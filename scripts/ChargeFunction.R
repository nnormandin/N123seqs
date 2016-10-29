FindCharge <- function(column){
  
  Dcount <- str_count(column, 'D')
  Ecount <- str_count(column, 'E')
  Hcount <- str_count(column, 'H')
  Kcount <- str_count(column, 'K')
  Rcount <- str_count(column, 'R')
  
  charge <- (-1*Dcount) + (-1*Ecount)+Hcount+Kcount+Rcount
  return(charge)
}

AddCharge <- function(lib, column){
  
  colidx <- which.max(sapply(colnames(lib),
                             function(x) grepl(as.character(column), x)))
  
  charge <- FindCharge(lib[,colidx])
  
  lib <- dplyr::mutate(lib, Charge = charge)
  
  return(lib)
  
}