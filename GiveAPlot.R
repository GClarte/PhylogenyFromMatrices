#we use ggplot for some nice plots
#and phytools for the consensus trees
library(ggplot2)
library(phytools)

#This file uses the result of the function SMCbruit stored in VV2.

npart=length(VV2[[1]])

#name of you experiment, will be used for recording the outputs
name="EABNZ_24_17h"

simu=F

#===========================================================================
#                         PARAMETERS PLOTS
#===========================================================================
#for the lambdas
lat=t(sapply(VV2[[1]],function(x){x[[2]]}))
lat=data.frame(value=c(lat),character=c(t(matrix(as.character(1:length(VV2[[1]][[1]][[2]])),nrow=length(VV2[[1]][[1]][[2]]),ncol=npart))))

h = ggplot(data=lat) + geom_density(aes(x=value,colour=character)) + theme(axis.title.y=element_blank(),
                                                                           axis.text.y=element_blank(),
                                                                           axis.ticks.y=element_blank()) +
  labs(x=expression(lambda))

if(simu){
  h=h+geom_vline(xintercept = Lareel[1])}
h
ggsave(paste("lambdas_",name,".pdf",sep=""),height=5,width=6.5,unit="cm")

lat=t(sapply(VV2[[1]],function(x){x[[2]][c(1,2)]}))
lat=data.frame(value=c(lat),character=c(t(matrix(as.character(1:2),nrow=2,ncol=npart))))

h = ggplot(data=lat) + geom_density(aes(x=value,colour=character)) + theme(axis.title.y=element_blank(),
                                                                           axis.text.y=element_blank(),
                                                                           axis.ticks.y=element_blank())+ 
  labs(x=expression(lambda))

if(simu){
  h=h+geom_vline(xintercept = Lareel[1])}
h
ggsave(paste("lambdas_1_",name,".pdf",sep=""),height=5,width=7,unit="cm")

lat=t(sapply(VV2[[1]],function(x){x[[2]][c(3,4)]}))
lat=data.frame(value=c(lat),character=c(t(matrix(as.character(3:4),nrow=2,ncol=npart))))
h = ggplot(data=lat) + geom_density(aes(x=value,colour=character)) + theme(axis.title.y=element_blank(),
                                                                           axis.text.y=element_blank(),
                                                                           axis.ticks.y=element_blank()) + 
  labs(x=expression(lambda))

if(simu){
  h=h+geom_vline(xintercept = Lareel[3])}
h
ggsave(paste("lambdas_2_",name,".pdf",sep=""),height=5,width=7,unit="cm")


#for Rho
h=ggplot(data=data.frame(value=sapply(VV2[[1]],function(x){x[[4]]})),aes(x=value))+geom_density()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + labs(x=expression(rho))

if(simu){
  h=h+geom_vline(xintercept = Rhoreel)}
h
ggsave(paste("rho_",name,".pdf",sep=""),height=5,width=7,unit="cm")


#for the noise
h=ggplot(data=data.frame(value=sapply(VV2[[1]],function(x){x[[8]]})),aes(x=value))+geom_density()+  theme(axis.title.y=element_blank(),
                                                                                                                       axis.text.y=element_blank(),
                                                                                                                       axis.ticks.y=element_blank())+
  labs(x=expression(nu)) 
if(simu){
  h=h+geom_vline(xintercept =Bruitreel)}
h
ggsave(paste("bruit_",name,".pdf",sep=""),height=5,width=7,unit="cm")

#for the betas
betat=t(sapply(VV2[[1]],function(x){x[[7]]}))
betat=data.frame(value=c(betat),character=c(t(matrix(as.character(1:length(VV2[[1]][[1]][[2]])),nrow=length(VV2[[1]][[1]][[2]]),ncol=npart))))

h = ggplot(data=betat) + geom_density(aes(x=value,colour=character)) + theme(axis.title.y=element_blank(),
                                                                             axis.text.y=element_blank(),
                                                                             axis.ticks.y=element_blank()) + 
  labs(x=expression(beta))
if (simu) {
  h=h+geom_vline(xintercept = c(1,1))
}
h
ggsave(paste("beta_",name,".pdf",sep=""),height=5,width=7,unit="cm")

#transformation probabilities
for (ch in 1:nch){
  
  Prt=sapply(VV2[[1]],function(x){x$P[[ch]]})
  
  qq=numeric()
  qq2x=numeric()
  qq2y=numeric()
  qq3x=numeric()
  qq3y=numeric()
  if (simu){
    nomsval=as.character(1:(nph[ch]))
  } else {
    nomsval=levels(J[,qui[ch]])
  }
  for (j in 1:nrow(Prt)){
    qq=c(qq,Prt[j,])
    qq2x=c(qq2x,rep(nomsval[passages[[ch]][[j]][1]],ncol(Prt)))
    qq2y=c(qq2y,rep(nomsval[passages[[ch]][[j]][2]],ncol(Prt)))
    qq3x=c(qq3x,nomsval[passages[[ch]][[j]][1]])
    qq3y=c(qq3y,nomsval[passages[[ch]][[j]][2]])
    if(simu){
      true=c(1,0,1,1,
             1,1,0,0,
             0,1,1,0,
             1,0,1,1,
             1,0,0,1)
    }
  }
  if(simu){
    dataPrt <-data.frame(Freq=qq,Transition_x=qq2x,Transition_y=qq2y)
    dataPrtTrue=data.frame(Truth=true,Transition_x=qq3x,Transition_y=qq3y)
  } else{
    dataPrt <-data.frame(Freq=qq,Transition_x=qq2x,Transition_y=qq2y)
  }
  
  h=ggplot(dataPrt,aes(Transition_y,Transition_x,fill=Freq,col=Freq)) +
    geom_tile()+guides(color = FALSE)+
    theme_minimal()+
    theme(legend.position = "left",axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 45, hjust=1))+
    scale_fill_gradient2(low="white",mid="yellow",high="red")+
    scale_size_manual(values=c(dot=NA,no_dot=1),guide="none")
  #scale_color_gradient2(low="white",mid="yellow",high="red")
  if(simu){
    h=h+
      geom_point(data=dataPrtTrue,aes(size=ifelse(Truth,"no_dot","dot"),x=qq3y,y=qq3x),inherit.aes = F)
  }
  
  ggsave(paste("transformations_ch",as.character(ch),name,".pdf",sep=""), height=6,width=6,unit="cm")
}
    #==========================================================================
#                       ROOT AGE
#==========================================================================
h=ggplot(data=data.frame(value=sapply(VV2[[1]],function(x){hauteur(x[[1]],x[[1]]$root)})),aes(x=value))+geom_density()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(x=expression(t[r])) #+geom_vline(xintercept =0.0001*Param$agemax)

if (simu) {
  h=h+geom_vline(xintercept = hauteur(tr,tr$root[[2]]))
}
h
ggsave(paste("root_age_",name,".pdf",sep=""),height=5,width=7,unit="cm")

#==========================================================================
#                       AGE OF AN INTERNAL NODE
#==========================================================================

tipsplot=c(16,17)
nom1=Param$tiplabel[tipsplot[1]]
nom2=Param$tiplabel[tipsplot[2]]

h=ggplot(data=data.frame(value=sapply(VV2[[1]],function(x){hauteur(x[[1]],nearestcommonancestor(x[[1]],tipsplot))})),aes(x=value))+geom_density()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(x=paste("Posterior age NCA(",nom1,",",nom2,")"))

if (simu) {
  h=h+geom_vline(xintercept = hauteur(tr,nearestcommonancestor(tr,tipsplot)))
}
h
ggsave(paste("nca_",as.character(tipsplot[1]),"_",as.character(tipsplot[2]),"_age",name,".pdf",sep=""),height=5,width=7,unit="cm",device=NULL)



#==========================================================================
#                       CONSENSUS TREES
#==========================================================================


#we include a last resampling, if the weights are not all the same
#the format of the tree we use is not exactly the one from phytools, if you try plotting one of the tree, there is a high chance you get an infinite loop somewhere in phytool code
npart=ncol(VV2[[2]])
pds=sapply(VV2[[1]],function(x){x$weight})
JJJ=1:npart
if (length(table(pds))!=1){
  JJJ=sample(1:npart,replace = T,prob = sapply(pds,length))
}

multiphylo=lapply(JJJ,function(x){VV2[[1]][[x]][[1]]})
class(multiphylo) <- "multiPhylo"

writeNexus(multiphylo,paste(name,".nex",sep=""))
multiphylo=read.nexus(paste(name,".nex",sep=""))
writeNexus(consensus.edges(multiphylo),paste(name,"_consensus.nex",sep=""))#can take a few seconds
plot(consensus.edges(multiphylo))

#true tree
writeNexus(tr,paste(name,"_vrai.nex",sep=""))


#In order to get a consensus tree with all interesting value, we used TREEANNOTATOR, the command you can use is :
#treeannotator -target name_consensus.nex name.nex name_annote.nex
#the resulting nexus file (name_annote.nex) can be opened with Figtree

#=======================================================================
#                             COGNACY
#=======================================================================

#pour les cognacy

trouvercheminaux <- function(tr,a,b,chemin){
  #trouve le chemin dans tr de a à b
  if (a==b){
    return(chemin)
  } else if (a == tr$root){
    w=which(tr$edge[,2]==b)
    return(trouvercheminaux(tr,a,tr$edge[w,1],c(chemin,w)))
  } else if(b==tr$root){
    w=which(tr$edge[,2]==a)
    return(trouvercheminaux(tr,tr$edge[w,1],b,c(chemin,w)))
  } else {
    w1=which(tr$edge[,2]==a)
    w2=which(tr$edge[,2]==b)
    return(trouvercheminaux(tr,tr$edge[w1,1],tr$edge[w2,1],c(chemin,w1,w2)))
  }
}
#which meanings apear in the dataset
quelmeanings=1:100
#two indices of the leaves for which you want to check cognacy status
Lg1=4
Lg2=6

tr$root=tr$root[[2]]

if(simu){
  vrai=rowSums(Lreelnb[,trouvercheminaux(tr,Lg1,Lg2,numeric())])==0 #true cognacy status if available
}
UU=sapply(1:npart, function(x){rowSums(VV2[[1]][[x]][[9]][,trouvercheminaux(VV2[[1]][[x]][[1]],Lg1,Lg2,numeric())])==0})
datacogn <- data.frame(Posterior=rowMeans(UU),Nb=1:nrow(Dat[[1]]),Meaning=quelmeanings)
datacogn$Meaning=factor(quelmeanings,levels=as.character(quelmeanings))
h=ggplot(data=datacogn,aes(x=Nb,y=Posterior))+geom_point() +geom_text(aes(label=Meaning),hjust=0, vjust=0)
if(simu){
  h=ggplot(data=datacogn,aes(x=Nb,y=Posterior,color=vrai))+geom_point() +geom_text(aes(label=Meaning),hjust=0, vjust=0)
}
h
ggsave(paste("cognats_",name,"_",as.character(Lg1),"_",as.character(Lg2),".pdf",sep=""),width=17,height=10,unit="cm")


#================================================================================
#                                    INTERNAL STATES
#================================================================================

#beware that the following code must be adapted when changing the number of characters, the values, etc.
library(tidyr)

loiini=Param$loiini
ncogn=nrow(Dat[[1]])

O=lapply(VV2[[1]],function(x){evolloitout(x,Param)})#may take a few minutes

if(simu){
  
  O1=lapply(O,function(x){t(t(x[[2]][[1]])/colSums(x[[2]][[1]]))})
  O2=lapply(O,function(x){t(t(x[[2]][[2]])/colSums(x[[2]][[2]]))})
  O3=lapply(O,function(x){t(t(x[[2]][[3]])/colSums(x[[2]][[3]]))})
  O4=lapply(O,function(x){t(t(x[[2]][[4]])/colSums(x[[2]][[4]]))})
  
  O1mean=Reduce('+',O1)/npart
  O2mean=Reduce('+',O2)/npart
  O3mean=Reduce('+',O3)/npart
  O4mean=Reduce('+',O4)/npart
  
} else {
  
  O1=lapply(O,function(x){t(t(x[[1]][[1]])/colSums(x[[1]][[1]]))})
  O2=lapply(O,function(x){t(t(x[[1]][[2]])/colSums(x[[1]][[2]]))})
  O3=lapply(O,function(x){t(t(x[[1]][[3]])/colSums(x[[1]][[3]]))})
  
  O1mean=Reduce('+',O1)/npart
  O2mean=Reduce('+',O2)/npart
  O3mean=Reduce('+',O3)/npart
  
}


#If you have your true data, here D1[,17]
if(simu){
  trueO1=matrix(0,ncol=ncogn,nrow=length(loiini[[1]]))
  trueO4=matrix(0,ncol=ncogn,nrow=length(loiini[[4]]))
  for (i in 1:ncogn){
    trueO1[D1[i,17],i]=1
    trueO4[D4[i,17],i]=1
  }
}

if(simu){
  dataO1=data.frame("Posterior"=as.vector(O1mean),
                    "Value"=rep(1:length(loiini[[1]]),ncogn),
                    "Meaning"=as.vector(sapply(1:ncogn,function(x){rep(x,length(loiini[[1]]))})),
                    "Truth"=as.vector(trueO1))
  dataO4=data.frame("Posterior"=as.vector(O4mean),
                    "Value"=rep(1:length(loiini[[4]]),ncogn),
                    "Meaning"=as.vector(sapply(1:ncogn,function(x){rep(x,length(loiini[[4]]))})),
                    "Truth"=as.vector(trueO4))
} else {
  dataO1=data.frame("Posterior"=as.vector(O1mean),
                    "Value"=rep(levels(J[,qui[1]]),ncogn),
                    "Meaning"=as.vector(sapply(1:ncogn,function(x){rep(x,length(loiini[[1]]))})))
  dataO2=data.frame("Posterior"=as.vector(O2mean),
                    "Value"=rep(levels(J[,qui[2]]),ncogn),
                    "Meaning"=as.vector(sapply(1:ncogn,function(x){rep(x,length(loiini[[2]]))})))
  dataO3=data.frame("Posterior"=as.vector(O3mean),
                    "Value"=rep(levels(J[,qui[3]]),ncogn),
                    "Meaning"=as.vector(sapply(1:ncogn,function(x){rep(x,length(loiini[[3]]))})))
}

h=ggplot(dataO1,aes(Value,Meaning,fill=Posterior)) +
  geom_tile()+
  theme_minimal()+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank())+
  scale_fill_gradient2(low="white",mid="yellow",high="red")
  #scale_color_gradient2(low="white",mid="yellow",high="red")
if(simu){
  h=h+
    geom_point(aes(size=ifelse(Truth,"no_dot","dot"),x=Value,y=Meaning),inherit.aes = F)+
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
    scale_size_manual(values=c(dot=NA,no_dot=1),guide="none")
}
h

ggsave(paste("internal_value_6_2_ch1",name,".pdf",sep=""),width=5,height=15,unit="cm")

h=ggplot(dataO4,aes(Value,Meaning,fill=Posterior)) +
  geom_tile()+
  theme_minimal()+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank())+
  scale_fill_gradient2(low="white",mid="yellow",high="red")
  #scale_size_manual(values=c(dot=NA,no_dot=1),guide="none")
#scale_color_gradient2(low="white",mid="yellow",high="red")
if(simu){
  h=h+
    geom_point(aes(size=ifelse(Truth,"no_dot","dot"),x=Value,y=Meaning),inherit.aes = F)+
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
    scale_size_manual(values=c(dot=NA,no_dot=1),guide="none")
}
h

ggsave(paste("internal_value_6_2_ch4",name,".pdf",sep=""),width=3.5,height=15,unit="cm")



#================================================================================
#                                    OTHER CHECKS
#================================================================================

#You can compute the likelihood of the particles to check
lkldtraitement <- function(v){
  nch=Param$nch
  Trposs = lapply(1:nch, function(y) {
    lapply(Param$Trpossliste[[y]], function(x) {
      A = diag(rep(1, Param$nph[y]))
      A[x[1], x[2]] = v[[7]][y]
      A[x[1], x[1]] = 1 - v[[7]][y]
      return(A)
    })
  })
  L = v[[9]]
  L[v[[9]] != 0] = rbeta(length(which(v[[9]] != 0)), L[v[[9]] != 0], 1)
  loiapplistetot = lapply(1:nch, function(x) {
    lapply(1:length(v[[1]]$edge.length), function(y) {
      loiapp(v[[5]][[x]][[y]], L[L[, y] > 0, y], loiini[[x]], Trposs[[x]], v[[8]], v[[1]]$edge.length[y], v[[6]][[x]][[y]])
    })
  })
  Mt = lapply(1:nch, function(z) {
    lapply(1:length(v[[1]]$edge.length), function(y) {
      transitionmatrix(v[[5]][[z]][[y]], Trposs[[z]], Param$nph[z], v[[8]], v[[1]]$edge.length[y])
    })
  })
  loiapplistetot = lapply(1:nch, function(x) {
    lapply(1:length(v[[1]]$edge.length), function(y) {
      loiapp(v[[5]][[x]][[y]], L[L[, y] > 0, y], Param$loiini[[x]], Trposs[[x]], v[[8]], v[[1]]$edge.length[y], v[[6]][[x]][[y]])
    })
  })
  lkldinter = lapply(1:nch, function(x) {
    pruninglkldlistini(v[[1]],
                       v[[5]][[x]],
                       Mt[[x]],
                       L,
                       Dat[[x]],
                       Param$nph[x],
                       1+length(Mt[[1]]),
                       v[[1]]$root,
                       Trposs[[x]],
                       Param$loiini[[x]],
                       bruit,
                       loiapplistetot[[x]])
  })
  tr=v[[1]]
  tr$root=list(c(1,2),v[[1]]$root)
  return(lkldtopotout(tr,v[[5]],v[[2]],v[[9]],v[[4]],L))
}

Lfin=sapply(VV2[[1]],lkldtraitement)
hist(Lfin)

#and plot the desperately degenerating genealogy of the particles

npart=ncol(VV2[[2]])
x1=numeric()
y1=numeric()
xend=numeric()
yend=numeric()
color=numeric()
actif=1:npart
for (i in ncol(VV2[[3]]):2){
  x1=c(x1,rep(i,npart))
  y1=c(y1,1:npart)
  xend=c(xend,rep(i-1,npart))
  yend=c(yend,VV2[[3]][,i])
  colortemp=rep(20/npart,npart)
  colortemp[actif]=1
  color=c(color,colortemp)
  actif=unique(VV2[[3]][actif,i])
}
geneal=data.frame(x1=x1,xend=xend,y1=y1,yend=yend,color=color)

b <- ggplot(geneal) + geom_segment(aes(x = x1, y = y1, xend = xend, yend = yend),alpha=color, data = geneal)+xlab("generation")+ylab("particle")
#b = b  + geom_vline(data=data.frame(pp=cumsum(sapply(Param$batches,length))),aes(xintercept=pp,color="red"))
b
ggsave(paste("geneal_",name,'.pdf',sep=""),height=5,width=10, unit="cm")

