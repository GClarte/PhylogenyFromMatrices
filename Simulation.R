
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



rdirichlet <- function (n, alpha) 
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}


set.seed(541)
tr =rcoal(15)[[1]]
tr$edge.length=300/hauteur(tr,16)*tr$edge.length

set.seed(NULL)

#S=matrix(c(1,2,2,1,2,3,3,2,3,4,4,3,4,1,1,4,4,2,2,4),nrow=2)
S = matrix(c(1, 2, 2, 1, 2, 3, 3, 2, 3, 4, 4, 3, 4, 5, 5, 4, 5, 1, 1, 5, 4, 1, 1, 4), nrow =
             2)
S2 = matrix(c(1, 2, 2, 1),nrow=2)
tauxtransf = 1
TrpossA = lapply(split(S, col(S)), function(x) {
  A = diag(rep(1, 5))
  A[x[1], x[2]] = tauxtransf
  A[x[1], x[1]] = 1 - tauxtransf
  return(A)
})
TrpossB = lapply(split(S2, col(S2)), function(x) {
  A = diag(rep(1, 2))
  A[x[1], x[2]] = tauxtransf
  A[x[1], x[1]] = 1 - tauxtransf
  return(A)
})
Trpossb = list(TrpossA, TrpossA, TrpossB, TrpossB)
Xreel = list(
  transf(tr, Trpossb[[1]], 0.01, rep(1 / 12, 12)),
  transf(tr, Trpossb[[2]], 0.01, rep(1 / 12, 12)),
  transf(tr, Trpossb[[3]], 0.001, rep(1 / 2, 2)),
  transf(tr, Trpossb[[4]], 0.001, rep(1 / 2, 2))
)

Lareel=c(.01,.01,.001,.001)

Tpsreel = lapply(Xreel, function(x) {
  lapply(x, function(z) {
    sort(runif(length(z), 0, 1))
  })
})

Rhoreel=.001

Lreelnb = matrix(rpois(100 * 28, .001 * tr$edge.length),
                 ncol = 28,
                 byrow = T)
Lreel = Lreelnb
Lreel[Lreelnb != 0] = rbeta(length(which(Lreelnb != 0)), Lreel[Lreelnb !=
                                                                 0], 1)
D1 = evolarbre(tr, Xreel[[1]], rep(1 / 5, 5), 100, 16, 29, Trpossb[[1]], 5, .001, Lreel)
Dred1 = D1[, 1:15]
D2 = evolarbre(tr, Xreel[[2]], rep(1 / 5, 5), 100, 16, 29, Trpossb[[2]], 5, .001, Lreel)
Dred2 = D2[, 1:15]
D3 = evolarbre(tr, Xreel[[3]], rep(1 / 2, 2), 100, 16, 29, Trpossb[[3]], 2, .001, Lreel)
Dred3 = D3[, 1:15]
D4 = evolarbre(tr, Xreel[[4]], rep(1 / 2, 2), 100, 16, 29, Trpossb[[4]], 2, .001, Lreel)
Dred4 = D4[, 1:15]
Dtot = list(Dred1, Dred2, Dred3, Dred4)

Bruitreel=.001

passages=c(lapply(1:2,function(x){split(S,col(S))}),lapply(1:2,function(x){split(S2,col(S2))}))


#paramètres


nch=4
Dat=Dtot
hyperpbini=lapply(passages,function(x){rep(1,length(x))})
Pri=list(c(1,100,0,.1),c(1,100,0,.1),c(1,500,0,.1),c(1,500,0,.1))
Prirho=list(function(x){dgamma(x,1,8000,log=T)
},list(
  function(){rgamma(1,1,800)},
  function(x){dgamma(x,1,800,log=T)}
))
loiini=list(rep(1/5,5),rep(1/5,5),rep(1/2,2),rep(1/2,2))
bruit=.0001

priortree1=function(x,agemax){
  if(length(x$root)==1){
    return(0)
  } else {
    WWW=0
    for (y in x$root){    
      AA=sousarbre(x,y)
      if (length(AA)>0){
        WWW=WWW+dbeta(hauteur(x,y)/agemax,1,length(AA))
      }
    }
      return(sum(WWW))
  }
}

priortree2=function(x,y){
  dbeta((sum(x$edge.length)-2*y)/((length(x$tip.label)-2)*y),1,5)
}

prioryule <- function(x,y){
  tips=Prior$tipprior
  if (!all(sapply(tips,function(tip){is.monophyletic.perso(x,tip)}))){
    return(-Inf)
  } else {
    if(length(x$root[[2]]>0)){
      K=sapply(x$root[[2]],function(z){y-hauteur(x,z)})
      return(-20*length(x$tip.label)*log(sum(x$edge.length)+sum(K)))
    } else {
      return(-20*length(x$tip.label)*log(sum(x$edge.length)))
    }}
}
pribruit=list(list(function(){runif(1,10^-6,10^-2)},function(x){
  if(x>10^-6 & x<10^-2){
    return(-log(x))
  } else {
    return(-Inf)
  }
}),
list(function(x){rnorm(1,x,.0001)},
     function(x,y){dnorm(x,y,.0001,log=T)}))
#list(function(x){exp(log(x)+rnorm(1,0,.00001))},
# function(x,y){-(log(y)-log(x))/2*.00001-log(x)}))

priorunif <- function(tr,agemax){
  tips=Prior$tipprior
  cladeage=Param$Cladeage
  if (!all(sapply(tips,function(tip){is.monophyletic.perso(tr,tip)})) ||
      !all(sapply(cladeage,function(clade){ageconstraints(tr,clade)}))){
    return(-Inf)
  } else {
    hhh=hauteur(tr,tr$root[[2]])
    if (!(hhh>200&hhh<1000)){
      return(-Inf)
    } else {
      if (length(cladeage)==0){
        return(-length(tr$tip.label)*log(hauteur(tr,tr$root[[2]])))
      } else {
        return(-(length(tr$tip.label)-sum(sapply(cladeage,function(x){length(x[[1]])-1})))*log(hauteur(tr,tr$root[[2]])))
      }
    }
  }
}

priorgamma <- function(tr,agemax){
  tips=Prior$tipprior
  cladeage=Param$Cladeage
  if (!all(sapply(tips,function(tip){is.monophyletic.perso(tr,tip)})) ||
      !all(sapply(cladeage,function(clade){ageconstraints(tr,clade)}))){
    return(-Inf)
  } else {
    hhh=hauteur(tr,tr$root[[2]])
    if (!(hhh>200&hhh<1000)){
      return(-Inf)
    } else {
      if (length(cladeage)==0){
        return(-length(tr$tip.label)*log(hauteur(tr,tr$root[[2]]))+
                 dgamma(sum(tr$edge.length),10,1/90,log=T))
      } else {
        return(-(length(tr$tip.label)-sum(sapply(cladeage,function(x){length(x[[1]])-1})))*log(hauteur(tr,tr$root[[2]]))+
                 dgamma(sum(tr$edge.length),10,1/90,log=T))
      }
    }
  }
}

Prior=list(Prirho=c(1,1000),
           hyperpbini=hyperpbini,
           Prila=Pri,
           priortree=priorunif,
           #priortree=priorcoal,
           #pribeta=c(.9,.9,.9,.9),
           pribeta=matrix(rep(c(10,1),nch),ncol=nch,nrow=2),
             bruittemp=c(.1,.05,.03,.02,.01,.009,.008,.007,.006,.005,.001),
           pribruit=pribruit,#prior sur le bruit avec tout défini
           tipprior=list(c(10,11,9),c(8,12,13,15)))
#tipprior=list())
Param=list(Trpossliste=passages,
           nch=nch,nph=c(5,5,2,2),
           npart=2000,
           npas=1000,
           npasfin=10000,
           Nmin=1000,
           loiini=loiini,
           Prob=c(1,.01,.1,.5,1,0,.01),
           Probfin=c(1,.01,.1,.5,1,.1,.01),
           agemax=c(10,1/30),
           tiplabel=tr$tip.label,
           batches=c(lapply(1:(ncol(Dat[[1]])),function(x){nrow(Dat[[1]])})),
           ReconstructClade=list(c(10,9),c(6,2)),
           Cladeage=list(list( c(10,9),c(10,200)),list( c(8,15),c(20,300)))
)

VV2=SMCbruit(Dat,Param,Prior,1,T)

save.image(paste("SMCforetsimu3", date(), ".RData"), version = 2)
