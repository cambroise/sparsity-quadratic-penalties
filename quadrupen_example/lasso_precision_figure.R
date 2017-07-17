###################################################################
####                                                           ####
####        Genomic selection for the CornFed dataset          ####
####                                                           ####
###################################################################


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

##Premiers pas...
fit.l1 <- glmnet(geno,yc, standardize=F)
plot(fit.l1,xvar="lambda")

#Choix du lambda par VC
repart <- sample(1:CV,Nind,replace=T)
Vfold.l1 <- cv.glmnet(geno,yc, foldid=repart, standardize=F)
plot(Vfold.l1)
Vfold.l1$lambda.1se
Vfold.l1$cvm[which(Vfold.l1$lambda==Vfold.l1$lambda.min)]

#Etude du Lasso optimal
fit.opt.l1 <- glmnet(geno,yc,lambda=Vfold.l1$lambda.min, standardize=F)
BHL1 <- fit.opt.l1$beta
NonNullCoeff1 <- which(BHL1 != 0)

#Etude du lasso optimal via quadrupen
fit.opt.l1.qpn <- elastic.net(geno, yc, lambda1 = Vfold.l1$lambda.min*Nind, lambda2 = 0, normalize = F)
BHL2 <- as.vector(fit.opt.l1.qpn@coefficients)
NonNullCoeff2 <- which(BHL2 != 0)

#Comparaison des deux algos
plot(BHL1,BHL2)
#plot(fit.opt.l1$beta,BHL2)
abline(a=0,b=1,col=2)
#Mais...
length(NonNullCoeff1)
length(NonNullCoeff2)
min(abs(BHL1[NonNullCoeff1]))

#Consequences : unicite
Crit1 <- abs(crossprod(geno,yc-geno%*%BHL1))
Crit1[NonNullCoeff1]     
Vfold.l1$lambda.min*Nind            

Crit2 <- abs(crossprod(geno,yc-geno%*%BHL2))
Crit2[NonNullCoeff2]   

lambda <- (Vfold.l1$lambda.min*Nind)
length(which(round(Crit2,10)==round(Vfold.l1$lambda.min*Nind,10)))

library(ggplot2)

abs.gradient <- data.frame(gradient = c(Crit2[NonNullCoeff2],Crit1[NonNullCoeff1])-lambda, method = rep(c("quadrupen", "glmnet"),sapply(list(Crit2[NonNullCoeff2],Crit1[NonNullCoeff1]),length)))

p1 <- ggplot(abs.gradient, aes(gradient,fill=method, group=method)) 
p1 <- p1 + theme_bw() + geom_histogram() + scale_x_log10()
p1 <- p1 + xlab("") + ylab("") + theme(text = element_text(size = 30)) + scale_fill_grey()

p2 <- ggplot(abs.gradient, aes(y=gradient,x=method, fill=method)) 
p2 <- p2 + geom_boxplot(outlier.colour="black", outlier.shape=16,  outlier.size=2,  notch = TRUE, notchwidth = 0.5)

ggsave("../figures/gradient.pdf", p1, width = 15, height=10)
