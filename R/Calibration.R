`Calibration` <- function(GD = NULL,testing = NULL, num = 100){
  # Objects: split the individuals to different groups based on genotype
  # input: GD, big matrix, m by n, including training and testing group
  #       testing: indicator for testing group
  #       num: how many individuals will be selected as training population
  # Output: cali, matrix, num by length(testing)
  n = ncol(GD)
  m = nrow(GD)
  n.test = length(testing)
  cali = matrix(NA,num,n.test)
  if ((n - n.test) < num){
    cali = matrix(1,(n - n.test),n.test)
    cali[1:nrow(cali),] = setdiff(seq(1:n),testing)
    return(cali)
  }else{
    cali.ind = setdiff(seq(1:n),testing)
    for ( i in 1:n.test){
        cali[,i] = CDcali(GD1,GD2,num)
    }
    return(cali)
  }
}

`CDcali` <- function(GD1,GD2,num){
  ## Objects: calculate CD
}
