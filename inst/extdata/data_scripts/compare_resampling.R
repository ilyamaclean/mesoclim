


var<-"Tmin - bilinear"
r<-mesoclimate$tmin
lyrstat<-global(r,mean,na.rm=TRUE)
qtls<-quantile(lyrstat[,1],prob=p)
sel<-c()
for(q in qtls) sel<-c(sel,which.min(abs(lyrstat[,1] - q)))
out_r1<-r[[sel]]
names(out_r1)<-paste(var,terra::time(r[[sel]]))

var<-"Tmin - cubic"
r<-mesoclimate3$tmin
lyrstat<-global(r,mean,na.rm=TRUE)
qtls<-quantile(lyrstat[,1],prob=p)
sel<-c()
for(q in qtls) sel<-c(sel,which.min(abs(lyrstat[,1] - q)))
out_r2<-r[[sel]]
names(out_r2)<-paste(var,terra::time(r[[sel]]))

var<-"Tmin - cubicspline"
r<-mesoclimate2$tmin
lyrstat<-global(r,mean,na.rm=TRUE)
qtls<-quantile(lyrstat[,1],prob=p)
sel<-c()
for(q in qtls) sel<-c(sel,which.min(abs(lyrstat[,1] - q)))
out_r3<-r[[sel]]
names(out_r3)<-paste(var,terra::time(r[[sel]]))

out_r<-c(out_r1,out_r2,out_r3)
p=c(0, 0.5, 1)
par(mfrow=c(3,3), mai=c(1,0.1,0.1,0.1))
plot(out_r,main=names(out_r),font.main=1, cex.main=1, nc=length(p))



var<-"Tmax - bilinear"
r<-mesoclimate$tmax
lyrstat<-global(r,mean,na.rm=TRUE)
qtls<-quantile(lyrstat[,1],prob=p)
sel<-c()
for(q in qtls) sel<-c(sel,which.min(abs(lyrstat[,1] - q)))
out_r1<-r[[sel]]
names(out_r1)<-paste(var,terra::time(r[[sel]]))

var<-"Tmax - cubic"
r<-mesoclimate3$tmax
lyrstat<-global(r,mean,na.rm=TRUE)
qtls<-quantile(lyrstat[,1],prob=p)
sel<-c()
for(q in qtls) sel<-c(sel,which.min(abs(lyrstat[,1] - q)))
out_r2<-r[[sel]]
names(out_r2)<-paste(var,terra::time(r[[sel]]))

var<-"Tmax - cubicspline"
r<-mesoclimate2$tmax
lyrstat<-global(r,mean,na.rm=TRUE)
qtls<-quantile(lyrstat[,1],prob=p)
sel<-c()
for(q in qtls) sel<-c(sel,which.min(abs(lyrstat[,1] - q)))
out_r3<-r[[sel]]
names(out_r3)<-paste(var,terra::time(r[[sel]]))

out_r<-c(out_r1,out_r2,out_r3)
p=c(0, 0.5, 1)
par(mfrow=c(3,3), mai=c(1,0.1,0.1,0.1))
plot(out_r,main=names(out_r),font.main=1, cex.main=1, nc=length(p))


var<-"Diurnal temp range - bilinear"
r<-mesoclimate$tmax-mesoclimate$tmin
lyrstat<-global(r,mean,na.rm=TRUE)
qtls<-quantile(lyrstat[,1],prob=p)
sel<-c()
for(q in qtls) sel<-c(sel,which.min(abs(lyrstat[,1] - q)))
out_r1<-r[[sel]]
names(out_r1)<-paste(var,terra::time(r[[sel]]))

var<-"Diurnal temp range - cubic"
r<-mesoclimate3$tmax-mesoclimate3$tmin
lyrstat<-global(r,mean,na.rm=TRUE)
qtls<-quantile(lyrstat[,1],prob=p)
sel<-c()
for(q in qtls) sel<-c(sel,which.min(abs(lyrstat[,1] - q)))
out_r2<-r[[sel]]
names(out_r2)<-paste(var,terra::time(r[[sel]]))

var<-"Diurnal temp range - cubicspline"
r<-mesoclimate2$tmax-mesoclimate2$tmin
lyrstat<-global(r,mean,na.rm=TRUE)
qtls<-quantile(lyrstat[,1],prob=p)
sel<-c()
for(q in qtls) sel<-c(sel,which.min(abs(lyrstat[,1] - q)))
out_r3<-r[[sel]]
names(out_r3)<-paste(var,terra::time(r[[sel]]))

out_r<-c(out_r1,out_r2,out_r3)
p=c(0, 0.5, 1)
par(mfrow=c(3,3), mai=c(1,0.1,0.1,0.1))
plot(out_r,main=names(out_r),font.main=1, cex.main=1, nc=length(p))
