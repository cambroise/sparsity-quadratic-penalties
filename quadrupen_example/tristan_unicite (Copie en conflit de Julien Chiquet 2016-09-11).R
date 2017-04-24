###################################################################
####                                                           ####
####        Genomic selection for the CornFed dataset          ####
####                                                           ####
###################################################################


rm(list=ls())
setwd('D:/Enseignement/Exterieur/MasterX/Exemples/DonneesCF_Dent_codees/')
library(quadrupen)
library(glmnet)


        #### Fonctions personnalisees


debut <- function(x){x[1:10,1:10]}
fin <- function(x){x[(nrow(x)-10):nrow(x),(ncol(x)-10):ncol(x)]}


        #### Parametres


NomGeno <- 'Genotypes.csv'
NomPheno <- 'Phenotypes.csv'
NomCarte <- 'carte.csv'
CV <- 5
PhenoInt <- "Tasseling_GDUb6"


        #### Importation donnees


geno <- t(read.table(file=NomGeno,header=T))
pheno <- read.table(file=NomPheno,header=T)
carte <-read.table(file=NomCarte,header=T)

debut(geno)
dim(geno)
debut(pheno)
dim(pheno)


        #### Controle qualite


#Variables constantes
Variances <- apply(geno,2,var)
VarNulles <- which(Variances==0)
geno <- geno[,-VarNulles]
carte <- carte[-VarNulles,]
dim(geno) 

#Heterozygotie residuelle            
unique(apply(geno,1,function(x) length(which((x==0)|(x==2)))))
dim(geno)

#Doublons            
length(unique(carte$Pos))
Doublons <- which(diff(carte$Pos)==0)
band(cor(geno[,sort(c(Doublons,Doublons+1))]),-1,1)
##M26264-M26265, M27119-M27120

#Corrections et definitions
Retirer <- which(colnames(geno) %in% c(paste('M',Doublons+1,sep=''),'M26264','M27119'))
geno <- geno[,-Retirer]
geno <- scale(geno)            
carte <- carte[-Retirer,]
Nind <- nrow(geno)
Nsnp <- ncol(geno)
y <- pheno[[PhenoInt]]
yc <- scale(y)


        #### Modele oligogénique


#Premiers pas...
fit.l1 <- glmnet(geno,yc, standardize=F)
plot(fit.l1,xvar="lambda")

#Choix du lambda par VC
Vfold.l1 <- cv.glmnet(geno,yc, foldid=repart, standardize=F)
plot(Vfold.l1)
Vfold.l1$lambda.1se
Vfold.l1$cvm[which(Vfold.l1$lambda==Vfold.l1$lambda.min)]

#Etude du Lasso optimal
fit.opt.l1 <- glmnet(geno,yc,lambda=Vfold.l1$lambda.min, standardize=F)
fit.opt.l1$df
length(which(fit.opt.l1$beta != 0))

#Unicite
Crit <- abs(crossprod(geno,yc-geno%*%fit.opt.l1$beta))
Crit[which(fit.opt.l1$beta != 0)]     
Vfold.l1$lambda.min*Nind            

