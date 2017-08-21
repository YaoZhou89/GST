# GST
Genomic Selection Tools

# Installation

library(devtools)

install_github("YaoZhou89/GST")

# Running the tests
source('../calcAIvar.R')

source('../iGS.eps.R')

source('../gBLUP.eps.R')

myY = read.table("phenotype.txt",head=T)

n = nrow(myY)

##kinship calculated by pepis, uncompress firstly
#http://bioinfo.noble.org/PolyGenic_QTL/

f = file("Ka.txt","r")

ka = matrix(NA,n,n)

for (i in 1:n){

  line = readLines(f,n=1)
  
  line.split = strsplit(line,"\t")
  
  line.split = unlist(line.split)
  
  len.line = length(line.split)
  
  ka[i,1:len.line] = line.split
  
  ka[1:len.line,i] = line.split
  
}

close(f)

##KD kinship

f = file("Kd.txt","r")

kd = matrix(NA,n,n)

for (i in 1:n){

  line = readLines(f,n=1)
  
  line.split = strsplit(line,"\t")
  
  line.split = unlist(line.split)
  
  len.line = length(line.split)
  
  kd[i,1:len.line] = line.split
  
  kd[1:len.line,i] = line.split
  
}

close(f)

ka = matrix(sapply(ka,as.numeric),n,n)

kd = matrix(sapply(kd,as.numeric),n,n)

K = list(ka=ka,kd = kd)

myGS = iGS.eps(Y= myY,K = K,alg="fs")

BLUP =  myGS$BLUP[,1]

phe = myGS$phe[,1]

# Authors
Yao Zhou (yao_wsu.zhou@wsu.edu)
