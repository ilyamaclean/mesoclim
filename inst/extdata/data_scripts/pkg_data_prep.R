### UKCP18 inputs for the southwest May 2018

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
enddate<-as.POSIXlt('2018/05/31') ### CHANGE FUNCTION SO THIS IS END DATE NOT ENDDATE +1

# Processes using already downloaded ukcp18rcm files in dir_data
t0<-now()
ukcpinput<-ukcp18toclimarray(dir_ukcp, dtmc,  startdate, enddate,
                             collection, domain, member)
print(now()-t0)
crs(ukcpinput$dtm)<-"EPSG:27700"
ukcpinput$dtm<-wrap(ukcpinput$dtm)
usethis::use_data(ukcpinput,overwrite=TRUE)

#write_climdata(ukcpinput,"data/ukcpinput.rda",overwrite=TRUE)

# UKCP land sea mask


#### Sea Surface temperature data - to match ukcpinput$dtm
dir_sst<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/SST"
sst<-create_ukcpsst_data(dir_sst,as.POSIXlt('2018/05/01'),as.POSIXlt('2018/05/31'),dtmc=ukcpinput$dtm, member="01")
plot(c(sst,ukcpinput$dtm))


#### Parcels shape file



#### Create 50m fine scale dtm of Lizard
dir_terrain50<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50"
dtmuk<-rast(file.path(dir_terrain50,"uk_dtm.tif"))
e<-ext(160000, 182000, 10000, 30000)
dtm<-crop(dtmuk,e)
lsmask<-vect("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Boundaries/CTRY_DEC_2023_UK_BGC.shp") %>% project(dtm)
dtm<-mask(dtm,lsmask)
dtm<-project(dtm,"EPSG:27700")
plot(dtm)
writeRaster(dtm,"inst/extdata/dtms/lizard50m.tif",overwrite=TRUE)

# Create smaller extent dtmf near Porthleven
e<-ext(162000, 167000, 23000, 29000)
dtm<-crop(dtmuk,e)
dtm<-mask(dtm,lsmask)
plot(dtm)
writeRaster(dtm,"inst/extdata/dtms/dtmf.tif",overwrite=TRUE)


#### Create corresponding medium dtm - coarser scale and wider extent than dtmf
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

#### Crop mesoclimate outputs to smaller area
daily100m<-lapply(mesoclimate,function(x) if(class(x)[1]=="SpatRaster") wrap(crop(x,e)) else x)
lapply(daily100m,class)
usethis::use_data(daily100m,overwrite=TRUE)

##### HadUK data for bias correction
dtm<-rast(system.file("extdata/dtms/dtmf.tif",package='mesoclim'))
plot(dtm)

dir_haduk<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Add_Trees/HadUK"
tasmin<-crop(rast(file.path(dir_haduk,"tasmin_hadukgrid_uk_1km_day_20180501-20180531.nc")),dtm)
tasmax<-crop(rast(file.path(dir_haduk,"tasmax_hadukgrid_uk_1km_day_20180501-20180531.nc")),dtm)
rain<-crop(rast(file.path(dir_haduk,"rainfall_hadukgrid_uk_1km_day_20180501-20180531.nc")),dtm)

crs(tasmin)<-'EPSG:27700'
crs(tasmax)<-'EPSG:27700'
crs(rain)<-'EPSG:27700'

writeRaster(tasmin,"inst/extdata/haduk/tasmin1km.tif",overwrite=TRUE)
writeRaster(tasmax,"inst/extdata/haduk/tasmax1km.tif",overwrite=TRUE)
writeRaster(rain,"inst/extdata/haduk/rainfall1km.tif",overwrite=TRUE)

##### ERA5 data for bias corrections - southwest


