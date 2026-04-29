########################## UKCP18 inputs for the southwest May 2018  ##########################
# UK 100m data
dtm100m<-rast("/Users/jonathanmosedale/Data/mesoclim_inputs/gb_100m_masked_dtm.tif") # 100 m raster masked for coast
# Coastal mask vector
coastal.v<-vect("/Users/jonathanmosedale/Data/mesoclim_inputs/CTRY_DEC_2023_UK_BFC.shp")

# Get dtmf
aoi<-ext(c(165000,170000,17000,22000))
dtmf<-crop(dtm100m,aoi)
plot(dtmf)

# Get dtmc by applying 10 km buffer
ukdtmc<-rast(system.file('extdata/ukcp18rcm/orog_land-rcm_uk_12km_osgb.nc',package='mesoclim'))
crs(ukdtmc)<-"EPSG:27700"
dtmf.bbox<-vect(ext(dtmf),crs="EPSG:27700")
wideraoi<-buffer(dtmf.bbox,10000)

dtmc<-terra::crop(ukdtmc,wideraoi,snap="out")
plot(dtmc)


# Get dtmm by aggreegating 
agg.f<-(res(dtmc)[1]/res(dtmf)[1])/10
dtmm<-crop(aggregate(dtm100m,agg.f,na.rm=TRUE),dtmc)
#dtmm<-terra::mask(dtmm,coastal.v,touches=TRUE)
#plot(crop(coastal.v,dtmm),col="red")
plot(dtmc)
plot(dtmm,add=T)
plot(dtmf,add=T)

writeRaster(dtmf,"inst/extdata/dtms/dtmf.tif",overwrite=TRUE)
writeRaster(dtmm,"inst/extdata/dtms/dtmm.tif",overwrite=TRUE)
writeRaster(dtmc,"inst/extdata/dtms/dtmf.tif",overwrite=TRUE)


########################## Preprocess UKCP18 data using constant albedo land / sea values
# TO DO - DROP cloud cover var??
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
# Save as arrays and packed spatraster
ukcpinput$dtm<-wrap(ukcpinput$dtm)

usethis::use_data(ukcpinput,overwrite=TRUE)

# write_climdata(ukcpinput,"data/ukcpinput.rda",overwrite=TRUE)

########################## Sea Surface temperature data - to match ukcpinput$dtm ##########################
dir_sst<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/SST"
ukcpsst<-create_ukcpsst_data(dir_sst,as.POSIXlt('2018/05/01'),as.POSIXlt('2018/05/31'),dtmc=dtmc, member="01")
plot(c(ukcpsst,dtmc))
usethis::use_data(ukcpsst,overwrite=TRUE)


########################## Preprocess ERA5 data
dir_era5data<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/era5"

# Get ERA5 dtmc by applying 10 km buffer
ancillary<-rast(file.path(dir_era5data,"era5_ancillary.nc"))
era5lsm<-ancillary$lsm %>%project("EPSG:4326")
if(max(values(era5lsm))!=1){
  mx<-as.numeric(global(era5lsm,max))
  mn<-as.numeric(global(era5lsm,min))
  era5lsm<- (era5lsm-mn) / (mx-mn)
}
plot(era5lsm)
# Calculate elevation from geopotential
RadEarth = 6371229 
gravity = 9.80665
era5elev<-(ancillary$z * RadEarth)/(gravity * RadEarth - ancillary$z)

latlonaoi<-project(wideraoi,crs(era5elev))
era5dtmc<-terra::crop(project(era5elev,"EPSG:4326"),latlonaoi,snap="out")

plot(era5dtmc)
plot(project(dtmf,"EPSG:4326"),add=TRUE)
writeRaster(era5dtmc,"inst/extdata/dtms/era5dtmc.tif",overwrite=TRUE)
writeRaster(era5lsm,"inst/extdata/dtms/era5lsm.tif",overwrite=TRUE)

ncfile<-file.path(dir_era5data,"era5_surface_ukeire_2018.nc")
# Processes using already downloaded ukcp18rcm files in dir_data
t0<-now()
era5input<-era5toclimarray(ncfile, dtmc=era5dtmc, lsm=era5lsm, aoi=latlonaoi, startdate= startdate, enddate=enddate)
now()-t0
dim(era5input$temp)
plot(era5input$dtm)
era5input$dtm<-wrap(era5input$dtm)
usethis::use_data(era5input,overwrite=TRUE)



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

##### Scotland HIghlands 20x20 1km res dtmf
# UKCP 2015
e<-ext(c( 250000, 270000, 920000, 940000))
highland_dtmf<-crop(aggregate(dtmuk,20),e)
writeRaster(highland_dtmf,"inst/extdata/dtms/altnaharra_1km.tif",overwrite=TRUE)
usethis::use_data(ukcphighland,overwrite=TRUE)
#Haduk 2015
writeRaster(hadukhland,"inst/extdata/haduk/altnaharra_rainfall.tif",overwrite=TRUE)

# climate data ukcp native res
sel<-which(lapply(ukcphighland,class)=="SpatRaster")
ukcphighland[sel]<-lapply(ukcphighland[sel],wrap)
usethis::use_data(ukcphighland,overwrite=TRUE)
ukcphighland<-read_climdata(mesoclim::ukcphighland)
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


##### Altnaharra daily climate data and dtmm
altnaharra_era5daily_2020
usethis::use_data(altnaharra_era5daily_2020,overwrite=TRUE)
altnaharra_5km<-aggregate(dtmuk,100) %>% crop(vect(st_buffer(locos.sf,50000)))
plot(altnaharra_5km)
plot(dtm1km,add=T)

sel<-which(lapply(shap_era5daily_2020,class)=="SpatRaster")
shap_era5daily_2020[sel]<-lapply(shap_era5daily_2020[sel],wrap)
usethis::use_data(shap_era5daily_2020,overwrite=TRUE)

