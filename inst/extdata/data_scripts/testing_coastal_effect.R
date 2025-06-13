tc
sstf
u2
wdir
ndri<-8
smooth<-5

# Downscale tmax - no cold air drainage
sst<-sstdata
thgto<-2
whgti<-10
rh<-climdata_m$relhum
pk<-climdata_m$pres
tmax<-climdata_m$tmax
thgti<-climdata_m$tempheight_m
whgti<-climdata_m$windheight_m
tme<-climdata_m$tme
u2<-.windhgt(climdata_m$windspeed, whgti, thgto) # wind at downscale temperature height above ground (for cad)
dtmc<-climdata_m$dtm
windspeed<-climdata_m$windspeed
winddir<-climdata_m$winddir

t<-tiles[[5]]
dtmt<-crop(dtmf,t)
plot(dtmt)
if(!is.logical(basins)) basins_tile<-crop(basins,t) else basins_tile<-NA
if(!is.logical(wca)) wca_tile<-.is(crop(.rast(wca,dtmf),t)) else wca_tile<-NA
if(!is.logical(skyview)) sky_tile<-crop(skyview,t) else sky_tile<-NA
if(!is.logical(horizon)) hor_tile<-crop(horizon,t) else hor_tile<-NA

# Calc windspeed at output height if required
uzf<-winddownscale(windspeed,winddir,dtmt,dtmm,dtmc,wca_tile,whgti,zo=2)



# Check lapserate!!! - grid effect!!
lrc<-lapserate(.is(tmax), .is(rh), .is(pk))
lrc<-.rast(lrc,dtmc)
if (crs(lrc) != crs(dtmt)) {
  lrcp<-project(lrc,crs(dtmt))
  lrf<-.resample(lrcp,dtmt)
} else lrf<-.resample(lrc,dtmt)
tmaxf<-.tempelev(tmax,lrf,lrc,dtmt+thgto,ifel(is.na(dtmc),0,dtmc+thgti))
plot(tmaxf[[1]])

# if coastal
# Interpolate sst to all landcells and across full timeseries so no NA
if (any(global(sst,anyNA))) sstinterp<-.spatinterp(sst) else sstinterp<-sst
sstinterp<-.tmeinterp(sstinterp,NA,tme)
if (crs(sst) != crs(dtmf)) sstinterp<-project(sstinterp,crs(dtmt))
sstf<-.resample(sstinterp,dtmt,method="cubic")
plot(sstf[[1]])




#tmaxf<-.tempcoastal(tc=tmaxf,sst=sstf,u2=uzf,wdir=climdata$winddir,dtmf,dtmm,dtmc)
tc<-tmaxf

u2<-uzf
wdir<-climdata_m$winddir
ndir<-8
smooth<-5

# Calculate coastal exposure for each wind direction
landsea<-ifel(is.na(dtmm),NA,1)
lsr<-array(NA,dim=c(dim(dtmt)[1:2],ndir))
for (i in 0:(ndir-1)) {
  r<-coastalexposure(landsea, e=ext(dtmt), i%%ndir*(360/ndir)) # seems to return a raster 1 row and col too big!
  r<-.correctcoastal(r)
  lsr[,,i+1]<-.is(r)
}
# smooth
lsr2<-lsr
for (i in 0:(ndir-1)) lsr2[,,i+1]<-0.25*lsr[,,(i-1)%%ndir+1]+0.5*lsr[,,i%%ndir+1]+0.25*lsr[,,(i+1)%%ndir+1]
lsr2<-.is( focal(.rast(lsr2,dtmt),w=smooth,fun="mean",na.policy="omit",na.rm=TRUE) )
plot(.rast(lsr2,dtmt))

lsm<-apply(lsr,c(1,2),mean)
tst<-min(lsm,na.rm=T)
# slot in wind speeds
if(!inherits(wdir,"SpatRaster")) wdr<-.rast(wdir,dtmc) else wdr<-wdir
ew<-ext(wdr)
xy<-data.frame(x=(ew$xmin+ew$xmax)/2,y=(ew$ymin+ew$ymax)/2)
wdir<-as.numeric(extract(wdr,xy))[-1]
if (is.na(wdir[1])) wdir<-apply(.is(wdr),3,median,na.rm=TRUE)
# Calculate array land-sea ratios for every hour
i<-round(wdir/(360/ndir))%%ndir
lsr<-lsr2[,,i+1]
# Calculate sstf weighting upwind
# derive power scaling coefficient from wind speed
p2<-0.10420*sqrt(.is(u2))-0.47852
# calculate logit lsr and lsm
llsr<-log(lsr/(1-lsr))
llsm<-.rta(log(lsm/(1-lsm)),dim(lsr)[3])
# Calculate mins and maxes
s<-which(is.na(lsr) == FALSE & lsr > 0 & lsr < 1)
llsr[lsr==0]<-log(min(lsr[s])/(1-min(lsr[s])))
llsr[lsr==1]<-log(max(lsr[s])/(1-max(lsr[s])))
s<-which(is.na(lsm) == FALSE & lsm > 0 & lsm < 1)
llsm[lsm==0]<-log(min(lsm[s])/(1-min(lsm[s])))
llsm[lsm==1]<-log(max(lsm[s])/(1-max(lsm[s])))
# predict logit swgt
lswgt<- -0.1095761+p2*(llsr+3.401197)-0.1553487*llsm
swgt<-.rast(1/(1+exp(-lswgt)),tc)
tcp<-swgt*sstf+(1-swgt)*tc


# calculate aggregation factor - DOESN'T WORK VERY WELL - perhaps focal smoothing would be better???
af<-res(dtmc)[1]/res(dtmt)[1]
tcc<-resample(aggregate(tcp,af,na.rm=TRUE),tcp)
final<-tc+(tcp-tcc)
# final3<-tc+(tcp-tcc)*swgt

out<-c(tc[[1]],tcp[[1]],final[[1]], final3[[1]])
names(out)<-c(paste("input",round(global(out[[1]],min,na.rm=T),0),round(global(out[[1]],max,na.rm=T),0)),
              paste("no correction",round(global(out[[2]],min,na.rm=T),0),round(global(out[[2]],max,na.rm=T),0)),
              paste("correction",round(global(out[[3]],min,na.rm=T),0),round(global(out[[3]],max,na.rm=T),0)),
              paste("correct_opt",round(global(out[[4]],min,na.rm=T),0),round(global(out[[4]],max,na.rm=T),0))
              )
panel(out)
dif<-c(out[[3]]-out[[1]],out[[2]]-out[[1]],out[[4]]-out[[1]])
names(dif)<-c("Original function","Orig excl last step","New v")
panel(dif)

