dir_root<-"D:/mesoclim_inputs"

# UKCP data
dir_ukcp<- file.path(dir_root,'ukcp18rcm')
dir_ukcpsst<- file.path(dir_root,'ukcp18sst')

# dir for writinng models and data for corretion
dir_bmodels<-file.path(dir_root,'bias_correct_models')
dir_bcdata<-'D:/mesoclim_outputs/bc_data'

# DTM and land/sea mask files
ukcpdtm_file<-file.path(dir_ukcp,'orog_land-rcm_uk_12km_osgb.nc')

ukdtm_file<-file.path(dir_root,'dtm','uk_dtm.tif') # At output resolution

# UKCP18rcm elevation
ukcpdtm<-rast(ukcpdtm_file)

ftr_startdate<-as.POSIXlt(paste0(2026,"/01/01"))
ftr_enddate<-as.POSIXlt(paste0(2026,"/12/31"))

# Get original ukcp18 data
collection<-'land-rcm'
domain<-'uk'
member<-'01'
rcp<-'rcp85'
datares<-'12km'

ftr_data<-ukcp18toclimarray(dir_ukcp, ukcpdtm,  ftr_startdate, ftr_enddate,
                            collection, domain, member,wind_hgt=10,toArrays=FALSE)


# Get corrected from applying bias correction models - 5year run
cordata<-read_climdata(file.path(dir_bcdata,"ukcp18rcm_2026_bc_2017_2020_models.Rds"))

# Plot difference
par(mfrow=c(1,3))
for (v in names(cordata)[5:12]){
  print(v)
  dif<-cordata[[v]]-ftr_data[[v]]
  plot(mean(dif),main=paste(v,'- mean difference'))
  plot(max(dif),main=paste(v,'- max difference'))
  plot(min(dif),main=paste(v,'- min difference'))

}

# As above using 4 extreme years to correct
cordata2<-read_climdata(file.path(dir_bcdata,"ukcp18rcm_2026_bc_4yr_models.Rds"))

# Plot difference
par(mfrow=c(1,3))
for (v in names(cordata2)[5:12]){
  print(v)
  dif<-cordata2[[v]]-ftr_data[[v]]
  plot(mean(dif),main=paste('Ext models',v,'- mean difference'))
  plot(max(dif),main=paste('Ext models',v,'- max difference'))
  plot(min(dif),main=paste('Ext models',v,'- min difference'))
}
