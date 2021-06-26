J=read.csv("datasets/full coded data.csv",header=T)

source("fonctions/fctauxini.R")
source("fonctions/gibbsparam.R")
source("fonctions/gibbstps.R")
source("fonctions/gibbstranfter.R")
source("fonctions/initialisation.R")
source("fonctions/lkldpruningbis.R")
source("fonctions/simudata.R")
source("fonctions/Gibbspartieldata.R")
source("fonctions/gibbsrho.R")
source("fonctions/modiftopo.R")
source("fonctions/SMCforet.R")
source("fonctions/SMCforet_aux.R")
source("fonctions/testclades.R")

library("parallel")

levels(J$handshape)=c(levels(J$handshape),"L","M","Autre")

J$handshape[J$handshape=="1_b"]="1"
J$handshape[J$handshape=="1_nb"]="1"
J$handshape[J$handshape=="4_nb"]="4"
J$handshape[J$handshape=="5_b"]="5"
J$handshape[J$handshape=="5_nb"]="5"
J$handshape[J$handshape=="B_b"]="B"
J$handshape[J$handshape=="B_nb"]="B"
J$handshape[J$handshape=="E_b"]="E"
J$handshape[J$handshape=="E_bnb"]="E"
J$handshape[J$handshape=="E_nb"]="E"
J$handshape[J$handshape=="L_b"]="L"
J$handshape[J$handshape=="L_bnb"]="L"
J$handshape[J$handshape=="L_eb"]="L"
J$handshape[J$handshape=="L_ebnb"]="L"
J$handshape[J$handshape=="L_enb"]="L"
J$handshape[J$handshape=="L_nb"]="L"
J$handshape[J$handshape=="M_enb"]="M"
J$handshape[J$handshape=="M_eb"]="M"
J$handshape[J$handshape=="S_b"]="S"
J$handshape[J$handshape=="V_nb"]="V"

J$handshape[J$handshape=="3_b"]="Autre"
J$handshape[J$handshape=="3"]="Autre"
J$handshape[J$handshape=="4"]="Autre"
J$handshape[J$handshape=="6"]="Autre"
J$handshape[J$handshape=="7_enb"]="Autre"
J$handshape[J$handshape=="K"]="Autre"
J$handshape[J$handshape=="Q"]="Autre"
J$handshape[J$handshape=="R"]="Autre"
J$handshape[J$handshape=="U"]="Autre"
J$handshape[J$handshape=="U_b"]="Autre"
J$handshape[J$handshape=="U_nb"]="Autre"
J$handshape[J$handshape=="W_nb"]="Autre"
J$handshape[J$handshape=="W"]="Autre"
J$handshape[J$handshape=="X"]="Autre"

J$handshape[J$handshape=="M"]="Autre"
J$handshape[J$handshape=="Y"]="Autre"
J$handshape[J$handshape=="I"]="Autre"
J$handshape[J$handshape=="H"]="Autre"

J$handshape=droplevels(J$handshape)

qui=c(3,5,17,19)
#qui=c(3,5)
#quelleslangues=levels(J[,2])[-c(6,3,12)]
#nlanguages=length(levels(J[,2])[-c(6,3,12)])
quelleslangues=levels(J[,2])
nlangues=length(levels(J[,2]))

rdirichlet <- function (n, alpha) 
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}

k=c()
for (i in 1:100){
  if (sum(as.integer(J[which(J[,2]%in%quelleslangues),1])==i)==nlangues & !is.na(sum(sapply(1:4,function(x){sum(as.integer(J[which(as.integer(J[,1])==i),qui[x]]))})))){
    k=c(k,i)
  }
}

quelmeanings=k

nch=length(qui)

Dtot=list()
for (i in 1:nch){
  Dtot[[i]]=matrix(ncol=nlangues,nrow=length(k))
  for (j in 1:length(quelleslangues)){
    Dtot[[i]][,j]=J[which(J[,2]==quelleslangues[j]& as.integer(J[,1])%in%k),qui[i]]
  }
}

quelmeanings=k

dist=matrix(ncol=24,nrow=24)

for (i in 1:24){
  for (j in 1:24){
    dist[i,j]=sum(sapply(1:length(Dtot),function(x){Dtot[[x]][,i]!=Dtot[[x]][,j]}))
  }
}

library(corrplot)
library(phytools)

foo2=dist[-c(10,4,13,22,23),-c(10,4,13,22,23)]

foo=foo2/max(foo2)
diag(foo)=NA
foo = foo - min(foo, na.rm=T)
foo = foo - max(foo, na.rm=T)/2
foo = foo/max(foo, na.rm=T)
diag(foo) = 0
foo = -foo

corrplot(foo)


ttt=nj(dist)
ttt=quelleslangues
plot(ttt)
  