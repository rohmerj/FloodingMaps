library(ggplot2)
library(viridis)
library(fields)
library(DiceKriging)
library(hetGP)

rm(list=ls())

set.seed(12345)

########################################
### FUNCTIONS
########################################
Q2<-function(y,yhat){
	return(1-sum((y-yhat)^2)/sum((y-mean(y))^2))
}

RMSE<-function(y,yhat){
	return(sqrt(mean((y-yhat)^2)))
}

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

########################################
### KM classique
########################################
LL <- RRall <- NULL
#####
BB <- 5
hh <- hist(Area[ff1],breaks=BB)
id <- cut(Area[ff1],hh$breaks)
OU <- NULL
ccc <- 0
NB <- c(25,10,5,5,2,1)
for (i in 1:(BB+1)){
	ou <- which(id==levels(id)[i])
	ss <- sample(1:length(ou),NB[i],replace=F)
	for (ii in 1:NB[i]){
		ccc <- ccc + 1
		if (ss[ii]>0) OU[ccc] <- ou[ss[ii]]
	}
}
print(Area[ff1][OU])
OU <- na.omit(OU)
#OU <- 1:length(ff1)
#####
plt <- p.hat <- p.true <- df0 <- list()

for (cc in OU){

X.tr <- doe1[-cc,]
n.tr <- nrow(X.tr)
X.te <- doe1[cc,]
n.te <- nrow(X.te)
HE.te <- HE1[cc,]
HE.tr <- HE1[-cc,]

pca <- prcomp(HE.tr, center = TRUE)
inertiey<-cumsum(pca$sdev^2/sum(pca$sdev^2))
ny <- 5#which.min((inertiey-0.999)^2)
BETA<-pca$x[,1:ny]
BASIS<-t(pca$rotation[,1:ny])
Y.tr <- BETA

ypred <- NULL
for (pc in 1:ny){

	KM<-km(formula=~., design=X.tr,response=BETA[,pc],nugget.estim=TRUE)
	df.te <- data.frame(matrix(unlist(X.te),ncol=d))
	names(df.te) <- names(X.tr)
	ypred[pc] <- predict(KM,newdata=df.te[1,],type="UK")$mean

}

truncy <- ypred %*% BASIS
Yt.pca <- t(truncy)+pca$center

#plot(leaveOneOut.km(KM, type="UK")$mean,BETA[,1])
#plot(Yt.pca,HE.te,main=X.te)
#abline(0,1)


df0[[cc]] <- df <- data.frame(x1=xx[ou.spa,1],x2=xx[ou.spa,2],He.hat=Yt.pca,He=HE.te,err=abs(HE.te-Yt.pca))
df$err <- cut(df$err,breaks=c(0,0.05,0.1,0.25,0.5,0.75,1,2))

p.hat[[cc]] <- ggplot(df,aes(x=x1,y=x2,color=abs(He.hat)))+geom_point()+ 
  scale_color_viridis(limits=c(0,max(max(df$He.hat),max(df$He))),option="D",direction=-1)
p.true[[cc]] <- ggplot(df,aes(x=x1,y=x2,color=abs(He)))+geom_point()+ 
  scale_color_viridis(limits=c(0,max(max(df$He.hat),max(df$He))),direction=-1)
plt[[cc]] <- ggplot(df,aes(x=x1,y=x2,color=err))+geom_point()+ 
  scale_color_viridis(discrete=T)

fff <- which(df$He>0)
(RRall[cc] <- RMSE(df$He[fff],df$He.hat[fff]))
(LL[cc] <- length(which(df$He.hat>0.01 & df$He==0))/length(ou.spa))

}## cc

########################################
### KM PC all
########################################
LL.PCall <- RRall.PCall <- NULL

#####
plt.PCall <- p.hat.PCall <- p.true.PCall <- df0.PCall <- list()

for (cc in OU){

X.tr <- doe1[-cc,]
n.tr <- nrow(X.tr)
X.te <- doe1[cc,]
n.te <- nrow(X.te)
HE.te <- HE1[cc,]
HE.tr <- HE1[-cc,]

pca <- prcomp(HE.tr, center = TRUE)
inertiey<-cumsum(pca$sdev^2/sum(pca$sdev^2))
ny <- 5#which.min((inertiey-0.999)^2)
BETA<-pca$x[,1:ny]
BASIS<-t(pca$rotation[,1:ny])
Y.tr <- BETA


ypred <- NULL
KM<-km(formula=~., design=X.tr,response=BETA[,1],nugget.estim=TRUE,scaling=T)
df.te <- data.frame(matrix(unlist(X.te),ncol=d))
names(df.te) <- names(X.tr)
ypred[1] <- predict(KM,newdata=df.te[1,],type="UK")$mean

for (pc in 2:ny){
	#X.tr1 <- data.frame(X.tr,PC=BETA[,1:(pc-1)])
	#X.te1 <- c(X.te,ypred[1:(pc-1)])
	X.tr1 <- data.frame(X.tr,PC=BETA[,(pc-1)])
	X.te1 <- c(X.te,ypred[(pc-1)])
	KM<-km(formula=~., design=X.tr1,response=BETA[,pc],nugget.estim=TRUE,scaling=T)
	df.te <- data.frame(matrix(unlist(X.te1),ncol=length(X.te1)))
	names(df.te) <- names(X.tr1)
	ypred[pc] <- predict(KM,newdata=df.te[1,],type="UK")$mean

}

truncy <- ypred %*% BASIS
Yt.pca <- t(truncy)+pca$center

#plot(Yt.pca,HE.te,main=X.te)
#abline(0,1)


df0.PCall[[cc]] <- df <- data.frame(x1=xx[ou.spa,1],x2=xx[ou.spa,2],He.hat=Yt.pca,He=HE.te,err=abs(HE.te-Yt.pca))
df$err <- cut(df$err,breaks=c(0,0.05,0.1,0.25,0.5,0.75,1,2))

p.hat.PCall[[cc]] <- ggplot(df,aes(x=x1,y=x2,color=abs(He.hat)))+geom_point()+ 
  scale_color_viridis(limits=c(0,max(max(df$He.hat),max(df$He))),option="D",direction=-1)
p.true.PCall[[cc]] <- ggplot(df,aes(x=x1,y=x2,color=abs(He)))+geom_point()+ 
  scale_color_viridis(limits=c(0,max(max(df$He.hat),max(df$He))),direction=-1)
plt.PCall[[cc]] <- ggplot(df,aes(x=x1,y=x2,color=err))+geom_point()+ 
  scale_color_viridis(discrete=T)

fff <- which(df$He>0)
(RRall.PCall[cc] <- RMSE(df$He[fff],df$He.hat[fff]))
(LL.PCall[cc] <- length(which(df$He.hat>0.01 & df$He==0))/length(ou.spa))

}## cc


cc <- OU[39]
summary(abs(plt[[cc]]$data[,"He"]-plt[[cc]]$data[,"He.hat"]))
summary(abs(plt.PCall[[cc]]$data[,"He"]-plt.PCall[[cc]]$data[,"He.hat"]))

plot(ecdf(abs(plt[[cc]]$data[,"He"]-plt[[cc]]$data[,"He.hat"])),main=Area1[cc])
lines(ecdf(abs(plt.PCall[[cc]]$data[,"He"]-plt.PCall[[cc]]$data[,"He.hat"])),col="green")

x11()
plot(plt[[cc]]$data[,"He"],plt[[cc]]$data[,"He.hat"],ylim=c(0,1),xlim=c(0,1))
abline(0,1)
points(plt.PCall[[cc]]$data[,"He"],plt.PCall[[cc]]$data[,"He.hat"],col=2)

x11()
#grid.arrange(p.true[[cc]],p.hat[[cc]],p.hat.PCall[[cc]])
plt.PCall[[cc]]

plot(1:7,as.numeric(table(cut(df0[[cc]]$err,breaks=c(0,0.05,0.1,0.25,0.5,0.75,1,2)))/length(ou.spa)),ylab="",ylim=c(0,1))
points(1:7+0.2,table(cut(df0.PCall[[cc]]$err,breaks=c(0,0.05,0.1,0.25,0.5,0.75,1,2)))/length(ou.spa),col=2)

save(df0,df0.PCall,plt,plt.PCall,file="./resu/Boucholeurs_test2.RData")

E <- E.PC <- RR <- RR.PC <- NULL
for (i in 1:length(OU)){
	cc <- OU[i]
	E[i] <- as.numeric(table(cut(df0[[cc]]$err,breaks=c(0,0.1,2)))/length(ou.spa))[1]
	E.PC[i] <- as.numeric(table(cut(df0.PCall[[cc]]$err,breaks=c(0,0.1,2)))/length(ou.spa))[1]
	RR[i] <- RMSE(plt[[cc]]$data[,"He"],plt[[cc]]$data[,"He.hat"])
	RR.PC[i] <- RMSE(plt.PCall[[cc]]$data[,"He"],plt.PCall[[cc]]$data[,"He.hat"])
}

plot(ecdf(E))
lines(ecdf(E.PC),col=2)


plot(ecdf((RR)))
lines(ecdf((RR.PC)),col=2)

