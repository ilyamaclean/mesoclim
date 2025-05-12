### UKCP18 inputs for the southwest May 2018

ukcp_aoi<-ext(24000, 336000, -60000, 204000 )
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
enddate<-as.POSIXlt('2018/06/01') ### CHANGE FUNCTION SO THIS IS END DATE NOT ENDDATE +1

# Processes using already downloaded ukcp18rcm files in dir_data
t0<-now()
ukcpinput<-ukcp18toclimarray(dir_ukcp, dtmc,  startdate, enddate,
                             collection, domain, member)
print(now()-t0)

ukcpinput$dtm<-wrap(ukcpinput$dtm)
usethis::use_data(ukcpinput,overwrite=TRUE)
#write_climdata(ukcpinput,"data/ukcpinput.rda",overwrite=TRUE)

#### dtm
dit_terrain50<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50"
dtmuk<-rast(file.path(dit_terrain50,"uk_dtm.tif"))
e<-ext(160000, 181000, 11000, 30000)
dtm<-crop(dtmuk,e)
lsmask<-vect("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Boundaries/CTRY_DEC_2023_UK_BGC.shp") %>% project(dtm)
dtm<-mask(dtm,lsmask)
dtm<-project(dtm,"EPSG:27700")
plot(dtm)
writeRaster(dtm,system.file("extdata/dtms/dtmf.tif",package='mesoclim'),overwrite=TRUE)

#### Crop mesoclimate outputs to smaller area
daily100m<-lapply(mesoclimate,function(x) if(class(x)[1]=="SpatRaster") wrap(crop(x,e)) else x)
lapply(daily100m,class)
usethis::use_data(daily100m)

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


