library(ggplot2)
library(viridis)
library(fields)

rm(list=ls())

set.seed(12345)

########################################
### DATA
########################################
NOM <- "Boucholeurs"
load(paste0("./",NOM,"/",NOM,"_zoom.RData"))
n1 <- dim(HE)[1]
n2 <- dim(HE)[2]
N <- dim(HE)[3]
d <- ncol(doe)
if (is.null(Smax)==FALSE) Area <- Smax

x1 <- 1:n1
x2 <- 1:n2
xx <- expand.grid(x1,x2)

image.plot(x1,x2,(matrix(PP,n1,n2)),col=c(rgb(0,0,0,0),tim.colors(10)))#,xlim=c(50,100),ylim=c(50,100))
ou.spa <- which(as.vector(PP)>0.20)
points(xx[ou.spa,1],xx[ou.spa,2],pch=".")

### HE where Proba>seuil
HE0 <- matrix(0,N,length(ou.spa))
for (k in 1:length(ou.spa)){
	HE0[,k]<-HE[xx[ou.spa[k],1],xx[ou.spa[k],2],]
}
HE0[which(is.na(HE0))] <- 0

### Cases with flooding
ff1 <- which(Area >0)
HE1 <- HE0[ff1,]
doe1 <- doe[ff1,]
Area1 <- Area[ff1]
