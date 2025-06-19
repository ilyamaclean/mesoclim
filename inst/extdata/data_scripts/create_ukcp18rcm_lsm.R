library(terra)
library(sf)


#### Get digital terrain at 50m and coastal vectors
dir_terrain50<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50"

dtmuk<-rast(file.path(dir_terrain50,"uk_dtm.tif"))
plot(dtmuk)

#bdry_gpkg<- "/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/bdline_gpkg_gb/Data/bdline_gb.gpkg"
#st_layers(bdry_gpkg)

lsmask<-st_read("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Boundaries/CTRY_DEC_2023_UK_BGC.shp")
lsmask<-st_union(lsmask[which(lsmask$CTRY23NM %in% c("Scotland","England","Wales")),]) %>% st_sf()

# Mask dtm with coastal vector - make sure same crs
crs(dtmuk)<-crs(lsmask)
dtmuk<-mask(dtmuk,vect(lsmask))
writeRaster(dtmuk,"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50/uk_seamask_dtm.tif",overwrite=TRUE)

# To include IoM at much lower res
lsmask2<-st_read("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Boundaries/digimap-gb-outlines/greatbritain.shp")
area<-st_area(lsmask2)
lsmask2<-lsmask2[order(area,decreasing = TRUE),]
iom<-lsmask2[8,]
st_crs(iom)<-st_crs(lsmask)
plot(iom)
lsmask<-rbind(lsmask,iom)
plot(lsmask$geometry)


# Calculate proportion of each ukcp cell that is land as defined by Terrain 50 and coastal mask - make sure same crs - restrict to UK wothout NI

# Grid of UK land cells in ukcp18rcm
ukcpgrid<-vect(file.path("/Users/jonathanmosedale/Data/ukcp-spatial-files-master/spatial-files/ukcp18-uk-land-12km/ukcp18-uk-land-12km.shp"))

# Get ukcp18 rcm dtm
ukcp18dtm<-rast(system.file("extdata/ukcp18rcm/orog_land-rcm_uk_12km_osgb.nc",package="mesoclim"))
crs(ukcp18dtm)<-crs(lsmask)
plot(ukcp18dtm); plot(vect(lsmask),add=T); plot(ukcpgrid,add=T)


ukcp18dtm<-crop(ukcp18dtm,lsm12km)
plot(ukcp18dtm)
plot(vect(ukcpgrid),add=T)


# Align dtmuk with ukcp18dtm
dtmuk<-extend(dtmuk,ukcp18dtm)
ext(dtmuk)<-ext(ukcp18dtm)
dtmuk<-crop(dtmuk,ukcp18dtm)
ext(dtmuk)<-ext(ukcp18dtm)

# Calculate 12km land sea mask by aggregation
lsm50m<-ifel(is.na(dtmuk),0,1)
plot(lsm50m)
aggfact<-res(ukcp18dtm)[1]/res(dtmuk)[1]
lsm12km.r<-aggregate(lsm50m,fact=aggfact,fun=sum)/aggfact^2
range(lsm12km.r)
plot(c(ukcp18dtm,lsm12km.r))
# Trim to GB land only with one cell padding
lsm12km.r<-trim(lsm12km.r,value=0,padding=1)
plot(lsm12km.r)

writeRaster(lsm12km.r,"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50/uk_seamask_12km.tif",overwrite=TRUE)

# Compare with HadUK 12km grid - only data where 12km ukcp18 elevation data!!
dir_had12km<-"/Users/jonathanmosedale/Data/haduk_12km"
tasmax<-rast(file.path(dir_had12km,"tasmax_hadukgrid_uk_12km_day_20180501-20180531.nc"))
plot(mask(lsm12km.r,crop(tasmax[[1]],lsm12km.r)))
range(mask(lsm12km.r,crop(tasmax[[1]],lsm12km.r)))





# Could plot mean/max/min difference between 12km Haduk obs and UKCP18 for same days
# Plot mean differences with proportion land/sea of each cell

######### Create HadUK 1km lsm mask ######
dtmuk<-rast("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50/uk_seamask_dtm.tif")
dir_haduk1km<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/haduk_1km_monthly"
r<-rast(file.path(dir_haduk1km,"tasmax_hadukgrid_uk_1km_day_20210901-20210930.nc"))[[1:10]]
crs(r)<-"EPSG:27700"
r<-crop(r,lsm12km.r)

lsm1km.r<-disagg(lsm12km.r,12)
writeRaster(lsm1km.r,"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Data/Terrain50/uk_seamask_1km.tif",overwrite=TRUE)

lsm<-ifel(lsm1km.r==0,NA,lsm1km.r)
r<-mask(r,lsm)

ri<-mask(.spatinterp(r),lsm)
plot(ri)

r<-ifel(is.na(r),NA,1)
r<-crop(extend(r,lsm12km.r),lsm12km.r)
r12km<-aggregate(r,12,fun=sum,na.rm=TRUE)/12^2
msk<-ifel(lsm12km.r==0,NA,1)
missing12km<-mask(msk,r12km,inverse=TRUE)
plot(missing12km,main="Cells without HadUK data")

missing1km<-disagg(missing12km,12)
plot(missing1km)
missdtm<-mask(extend(dtmuk,missing50m),missing50m)
plot(missdtm)




lsm50m<-ifel(is.na(dtmuk),0,1)
plot(lsm50m)
aggfact<-res(r)[1]/res(dtmuk)[1]
lsm1km.r<-aggregate(lsm50m,fact=aggfact,fun=sum)/aggfact^2
range(lsm1km.r)
plot(lsm1km.r)

rc<-crop(r,lsm1km.r)
plot(rc)

missingcoast<-mask(ifel(lsm1km.r==0,NA,lsm1km.r),rc,inverse=TRUE)
plot(missingcoast)

plot(c(r,lsm1km.r))
# Trim to GB land only with one cell padding
lsm12km.r<-trim(lsm12km.r,value=0,padding=1)
plot(lsm12km.r)


