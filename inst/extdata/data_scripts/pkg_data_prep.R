########################## UKCP18 inputs for the southwest May 2018  ##########################
ukcp_aoi<-ext(120000, 200000, -10000, 100000 )

dtmc_uk<-rast(system.file('extdata/ukcp18rcm/orog_land-rcm_uk_12km_osgb.nc',package='mesoclim'))
dtmc<-crop(dtmc_uk,ukcp_aoi)
dtmc<-project(dtmc,"EPSG:27700")
plot(dtmc)

# Preprocess UKCP18 data using constant albedo land / sea values
dir_ukcp<-"/Users/jonathanmosedale/Data/mesoclim_inputs"

collection<-'land-rcm'
domain<-'uk'
member<-'01'
rcp<-'rcp85'
startdate<-as.POSIXlt('2018/05/01')
enddate<-as.POSIXlt('2018/05/31')

# Processes using already downloaded ukcp18rcm files in dir_data
t0<-now()
ukcpinput<-ukcp18toclimarray(dir_ukcp, dtmc,  startdate, enddate,
                             collection, domain, member)
print(now()-t0)
crs(ukcpinput$dtm)<-"EPSG:27700"
ukcpinput$dtm<-wrap(ukcpinput$dtm)
usethis::use_data(ukcpinput,overwrite=TRUE)

#write_climdata(ukcpinput,"data/ukcpinput.rda",overwrite=TRUE)

########################## Sea Surface temperature data - to match ukcpinput$dtm ##########################
dir_sst<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/SST"
sst<-create_ukcpsst_data(dir_sst,as.POSIXlt('2018/05/01'),as.POSIXlt('2018/05/31'),dtmc=ukcpinput$dtm, member="01")
plot(c(sst,ukcpinput$dtm))


########################## UKCP18 inputs for southwest for the future May 2030  ##########################
startdate<-as.POSIXlt('2030/05/01')
enddate<-as.POSIXlt('2030/05/31')
ukcpfuture<-ukcp18toclimarray(dir_ukcp, dtmc,  startdate, enddate,
                             collection, domain, member)
print(now()-t0)
crs(ukcpfuture$dtm)<-"EPSG:27700"
ukcpfuture$dtm<-wrap(ukcpfuture$dtm)
usethis::use_data(ukcpfuture,overwrite=TRUE)

futuredata<-read_climdata(mesoclim::ukcpfuture)


########################## UKCP18 inputs for INLAND location for May 2018  ##########################
ukcp_aoi<-ext(400000, 480000, 250000, 350000 )
dtmc_uk<-rast(system.file('extdata/ukcp18rcm/orog_land-rcm_uk_12km_osgb.nc',package='mesoclim'))
dtmc_inland<-crop(dtmc_uk,ukcp_aoi)
dtmc_inland<-project(dtmc_inland,"EPSG:27700")

####  Create dtmf and dtmm for INLAND area ####
dir_terrain50<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50"
dtmuk<-rast(file.path(dir_terrain50,"uk_dtm.tif"))
dtmfe<-ext(435000, 440000, 300000, 305000)
dtmf_inland<-crop(dtmuk,dtmfe)
dtmm_inland<-crop(dtmuk,dtmc_inland)

plot(dtmc_uk)
plot(vect(ext(dtmc_inland)),add=TRUE)

plot(dtmc_inland)
plot(dtmm_inland,add=T)
plot(vect(ext(dtmf_inland)),col='red',add=TRUE)
plot(dtmf_inland)

writeRaster(dtmf,"inst/extdata/dtms/dtmf_inland.tif",overwrite=TRUE)
writeRaster(dtmm,"inst/extdata/dtms/dtmm_inland.tif",overwrite=TRUE)

startdate<-as.POSIXlt('2018/05/01')
enddate<-as.POSIXlt('2018/05/31')

ukcpinland<-ukcp18toclimarray(dir_ukcp, dtmc_inland,  startdate, enddate,
                              collection, domain, member)
print(now()-t0)
crs(ukcpinland$dtm)<-"EPSG:27700"
ukcpinland$dtm<-wrap(ukcpinland$dtm)
usethis::use_data(ukcpinland,overwrite=TRUE)

ukcpinland<-read_climdata(mesoclim::ukcpinland)


##########################  DTMs  ##########################
#### Create 50m fine scale dtm of Lizard = lizard50m.tif
dir_terrain50<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50"
dtmuk<-rast(file.path(dir_terrain50,"uk_dtm.tif"))
e<-ext(160000, 182000, 10000, 30000)
dtm<-crop(dtmuk,e)
lsmask<-vect("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Boundaries/CTRY_DEC_2023_UK_BGC.shp") %>% project(dtm)
dtm<-mask(dtm,lsmask)
dtm<-project(dtm,"EPSG:27700")
plot(dtm)
writeRaster(dtm,"inst/extdata/dtms/lizard50m.tif",overwrite=TRUE)

####  Create smaller extent dtmf near Porthleven = dtmf.tif
e<-ext(162000, 167000, 23000, 29000)
dtm<-crop(dtmuk,e)
dtm<-mask(dtm,lsmask)
plot(dtm)
writeRaster(dtm,"inst/extdata/dtms/dtmf.tif",overwrite=TRUE)

#### Create corresponding medium dtm - coarser scale and wider extent than dtmf = dtmm.tif
dir_terrain50<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50"
dtmuk<-rast(file.path(dir_terrain50,"uk_dtm.tif"))
lsmask<-vect("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Boundaries/CTRY_DEC_2023_UK_BGC.shp") %>% project("EPSG:27700")
e<-ext(130000, 200000, 10000, 70000)
dtm<-crop(dtmuk,e)
plot(dtm)
dtm<-project(dtm,"EPSG:27700")
dtm<-mask(dtm,lsmask)
dtmm<-mask(terra::aggregate(dtm,20,  na.rm=TRUE), lsmask)
plot(dtmm)
writeRaster(dtmm,system.file("extdata/dtms/dtmm.tif",package='mesoclim'),overwrite=TRUE)

#### Crop daily mesoclimate outputs to smaller area - NOT USED NOW
daily100m<-lapply(mesoclimate,function(x) if(class(x)[1]=="SpatRaster") wrap(crop(x,e)) else x)
lapply(daily100m,class)
usethis::use_data(daily100m,overwrite=TRUE)




########################## Parcels shape file ##########################



########################## Bias correction data ##########################

#### 1km Observational data for May 2018
climdata<-read_climdata(ukcpinput)
dtmc<-climdata$dtm
plot(dtmc)

dir_haduk<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/haduk_1km_monthly"

tasmin<-trim(crop(rast(file.path(dir_haduk,"tasmin_hadukgrid_uk_1km_day_20180501-20180531.nc")),dtmc))
tasmax<-trim(crop(rast(file.path(dir_haduk,"tasmax_hadukgrid_uk_1km_day_20180501-20180531.nc")),dtmc))
rain<-trim(crop(rast(file.path(dir_haduk,"rainfall_hadukgrid_uk_1km_day_20180501-20180531.nc")),dtmc))

crs(tasmin)<-'EPSG:27700'
crs(tasmax)<-'EPSG:27700'
crs(rain)<-'EPSG:27700'

writeRaster(tasmin,"inst/extdata/haduk/tasmin1km.tif",overwrite=TRUE)
writeRaster(tasmax,"inst/extdata/haduk/tasmax1km.tif",overwrite=TRUE)
writeRaster(rain,"inst/extdata/haduk/rainfall1km.tif",overwrite=TRUE)


### Mask of Land proportions in 12km UKCP grid cells
# See create_ukcp18rcm_lsm.R
lsm12km.r<-rast("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50/uk_seamask_12km.tif")
plot(lsm12km.r)
writeRaster(lsm12km.r,"inst/extdata/biascorrect/uk_seamask_12km.tif",overwrite=TRUE)

### Mask of land proportions in 1km HadUK cells
lsm1km.r<-rast("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50/uk_seamask_1km.tif")
plot(lsm1km.r)
writeRaster(lsm1km.r,"inst/extdata/biascorrect/uk_seamask_1km.tif",overwrite=TRUE)


##### ERA5 data for bias corrections - southwest


