
FindSigInColumn <- function(lib, column1, text1,
                     column2 = 0, text2 = 0,
                     Union = FALSE){
  
  # figure out which column number corresponds to column1 name  
  colidx <- which.max(sapply(colnames(lib),
                             function(x) grepl(as.character(column1), x)))
  
  if(text2 == 0){
    idx <- grep(as.character(text1), lib[,colidx])
    freq <- length(idx)/nrow(lib) * 100
    cat(paste('\n', freq, '% have that string'))
    return(lib[idx,])
  }
  
  else{
    
    if(column2 != 0){
      
      colidx2 <- which.max(sapply(colnames(lib),
                                  function(x) grepl(as.character(column2), x)))
      
    }
    else{
      colidx2 = colidx
    }
    
    idx1 <- grep(as.character(text1), lib[,colidx])
    idx2 <- grep(as.character(text2), lib[,colidx2])
    
    if(returnAll == TRUE){
      idx <- union(idx1, idx2)
    }
    else{
      idx <- intersect(idx1, idx2)
    }
    
    freq <- length(idx)/nrow(lib) * 100
    cat(paste('\n', freq, '% have that string'))
    return(lib[idx,])
    
  }
  
}