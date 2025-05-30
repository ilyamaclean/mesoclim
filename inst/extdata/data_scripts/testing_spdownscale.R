# Set up libraries files and directories for mac
setup_file<-"inst/extdata/data_scripts/setup_mesoclim_mac.R"
source(setup_file)

member<-"01"
startyear<-"2018"
endyear<-"2018"
startdate<-as.POSIXlt(paste0(startyear,'/05/01'),tz="UTC")
enddate<-as.POSIXlt(paste0(endyear,'/05/31'),tz="UTC")

dtmuk<-terra::rast(ukdtm_file) # 50m resolution dtm of all UK
coast_v<-terra::project(terra::vect(coast_file),dtmuk)

dtmf<-rast(climdata$dtm)

# input data
climdata<-ukcpinput
dtmc<-rast(climdata$dtm)
dtmm<-crop(dtmuk,dtmc) %>% aggregate(2)
dtmm<-mask(dtmm,coast_v)

dir_sst<-"inst/extdata/sst"

# Windshelter coef
t0<-now()
wca<-calculate_windcoeffs(dtmc,aggregate(dtmm,5),dtmf,zo=2)
plot(.rast(wca,dtmf))

# Cold air drainage basins - as above and ONLY if using coastal correction - can take several minutes for large areas
basins<-basindelin(dtmf, boundary = 2)
plot(basins)
print(now()-t0)

# Sea surface data
sstdata<-create_ukcpsst_data(dir_sst,startdate,enddate,dtmc,member)
plot(sstdata[[1]])

## spatial downscaling
wca<-calculate_windcoeffs(ukcpinput$dtm,dtmm,dtmf,zo=2)
basins<-basindelin(dtmf, boundary = 2)

# one step
mesoc<-spatialdownscale(climdata, sstdata, dtmf, dtmm, basins , wca, cad = TRUE,
                        coastal = TRUE, thgto =2, whgto=2,
                        rhmin = 20, pksealevel = TRUE, patchsim = TRUE, terrainshade = TRUE,
                        precipmethod = "Elev",fast = TRUE, noraincut = 0, toArrays=TRUE)


plothgto =2
whgto=2
rhmin = 20
pksealevel = TRUE
patchsim = FALSE
coastal<-TRUE
terrainshade = TRUE
precipmethod = "Elev"
fast = TRUE
noraincut = 0.01
#Test spatialdownscale function
tme<-as.POSIXlt(ukcpinput$tme,tz="UTC")

# Find out whether daily or hourly
tint<-as.numeric(tme[2])-as.numeric(tme[1])
hourly<-TRUE
if (abs(tint-86400) < 5) hourly=FALSE

# Extract variables - calculate tmean if daily
if(class(climdata$dtm)[1]=='PackedSpatRaster') dtmc<-rast(ukcpinput$dtm) else dtmc<-ukcpinput$dtm
if (hourly) {
  tc<-ukcpinput$temp
} else {
  tmin<-ukcpinput$tmin
  tmax<-ukcpinput$tmax
  tmean<-.hourtoday(.is(temp_dailytohourly(.rast(tmin,dtmc), .rast(tmax,dtmc), ukcpinput$tme)),mean) # required for relhum downscaling
}
rh<-ukcpinput$relhum
pk<-ukcpinput$pres
wspeed<-ukcpinput$windspeed
wdir<-ukcpinput$winddir
swrad<-ukcpinput$swrad
lwrad<-ukcpinput$lwrad
prec<-ukcpinput$prec
whgti<-ukcpinput$windheight_m
thgti<-ukcpinput$tempheight_m



### WIND - AGGREGATES dtmm - ADD CONDITIONAL ON IT BEINF G SAME res AS dtmf
uzf<-winddownscale(climdata$windspeed,climdata$winddir,dtmf,aggregate(dtmm,5),dtmc,wca,zi=climdata$tempheight_m,zo=thgto)
plot(uzf[[1]])

#### TEMPERATURE
tempvar<-"tmax"
rh<-climdata$relhum
pk<-climdata$pres
tc<-climdata[[tempvar]]
thgti<-climdata$tempheight_m
whgti<-climdata$windheight_m
tme<-climdata$tme

# Test cf of tmin / tmax
tminf<-tempdownscale(climdata,sstdata,dtmf,dtmm,basins,uzf,cad=TRUE,coastal=TRUE,'tmin',thgto,whgto)
tmaxf<-tempdownscale(climdata,sstdata,dtmf,dtmm,basins,uzf,cad=FALSE,coastal=TRUE,'tmax',thgto,whgto)




# Check tax>tmin (coastal effects can switch) - NEEDS FASTER METHOD!!!
t0<-now()
diurnal<-(tmaxf-tminf<0)
if(unique(all(diurnal))==0){
  temp_tmin<-ifel(diurnal,tmaxf,tminf)
  tmaxf<-ifel(diurnal,tminf,tmaxf)
  tminf<-temp_tmin
  remove(temp_tmin)
  remove(diurnal)
}
now()-t0






# CHECK BELOW
u2<-.windhgt(climdata$windspeed, whgti, thgto) # wind at downscale temperature height above ground (for cad)
swrad<-climdata$swrad
lwrad<-climdata$lwrad
dtmc<-rast(climdata$dtm)

## ELEVATION
# Calculate lapse rates (coarse and fine scale)
n<-dim(tc)[3]
lrc<-lapserate(.is(tc), .is(rh), .is(pk))
lrc<-.rast(lrc,dtmc)
if (crs(lrc) != crs(dtmf)) {
lrcp<-project(lrc,crs(dtmf))
lrf<-.resample(lrcp,dtmf)
} else lrf<-.resample(lrc,dtmf)

# Calculate elevation effects (assume NA in dtmc =sea = height 0) - add input / output heights to elevations
tcf<-.tempelev(tc,lrf,lrc,dtmf+thgto,ifel(is.na(dtmc),0,dtmc)+thgti)
plot(tcf[[1]])

## CAD
mu<-.cad_multiplier(dtmf, basins, refhgt=thgti)
plot(mu,main='mu')

st<-.cad_conditions(tc,u2,swrad=0,lwrad,dtmc,dtmf,refhgt=thgti) #swrad=0 as we assume tmin in night
tcad<-.apply_cad(lrf,mu,st)
plot(max(tcad),main="Max tcad")

tcf<-tcf+tcad
plot(tcf[[1]])


### COASTAL

# sst data for coastal
sstdata<-.spatinterp(sstdata)
sstdata<-.tmeinterp(sstdata,time(sstdata),climdata$tme)
sstf<-.resample(sstdata,dtmf,method="cubic")

# next gets errorError in lsr[, , i + 1] <- .is(r) : number of items to replace is not a multiple of replacement length
tcf2<-.tempcoastal(tc=tcf,sst=sstf,u2=uzf,wdir=climdata$winddir,dtmf,dtmm,dtmc)

### COASTAL microclim detail

########### MAKE SURE wind and temp passed as parameter are already downscaled !!!!!!!
u2<-uzf
wdir<-climdata$winddir
tc<-tcf

landsea<-ifel(is.na(dtmm),NA,1) # dtmm same resolution as dtmf
lsr<-array(NA,dim=c(dim(dtmf)[1:2],8))  ##
for (i in 0:7) {
  r<-coastalexposure(landsea, e=ext(dtmf), i%%8*45) # seems to return a raster 1 row and col too big!
  r<-.correctcoastal(r)
  lsr[,,i+1]<-.is(r) ######CHANGED THIS!!!!!
}# results looks OK
plot(.rast(lsr,dtmf))

# smooth
lsr2<-lsr
for (i in 0:7) lsr2[,,i+1]<-0.25*lsr[,,(i-1)%%8+1]+0.5*lsr[,,i%%8+1]+0.25*lsr[,,(i+1)%%8+1]
lsm<-apply(lsr,c(1,2),mean)
tst<-min(lsm,na.rm=T)
#. plot(.rast(lsm,dtmf)) = mean as per microclima lsa

  # slot in wind speeds
  wdr<-.rast(wdir,dtmc)
  ew<-ext(wdr)
  xy<-data.frame(x=(ew$xmin+ew$xmax)/2,y=(ew$ymin+ew$ymax)/2)
  wdir<-as.numeric(extract(wdr,xy))[-1]
  if (is.na(wdir[1])) wdir<-apply(.is(wdr),3,median,na.rm=TRUE)

  # Calculate array land-sea ratios for every hour
  i<-round(wdir/45)%%8
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
  plot(tcp[[1]],main='before aggreg')

  # calculate aggregation factor - DOESN'T WORK VERY WELL - perhaps focal smoothing would be better???
  af<-res(dtmc)[1]/res(dtmf)[1]
  tcc<-resample(aggregate(tcp,af,na.rm=TRUE),tcp)
  tcp<-tc+(tcp-tcc)
  plot((tc+(tcp-tcc))[[1]],main='after aggreg')


  testr<-focal(tcp,fun='mean',na.policy='omit')
  plot(testr,main='focal')
