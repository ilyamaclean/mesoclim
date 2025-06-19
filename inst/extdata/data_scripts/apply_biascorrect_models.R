#####################################################################
# Applies previously calculated bias correction models to ukcp18 derived
# data (as output by ukcp18toclimarray)
#####################################################################

library(devtools)
load_all()

library(terra)
library(sf)
library(mesoclim)
library(lubridate)
library(mesoclimAddTrees)
library(ncdf4)

#####################################################################
# Setup parameters etc
#####################################################################
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


###################################################################################
# Load bias correct models and ukcp weather data
###################################################################################

model_list<-readRDS(file.path(dir_bmodels,'bias_correct_models_4yrs.Rds'))
model_list<-readRDS(file.path(dir_bmodels,'bias_correct_models_2017_2021'))


collection<-'land-rcm'
domain<-'uk'
member<-'01'
rcp<-'rcp85'
datares<-'12km'

ftr_data<-ukcp18toclimarray(dir_ukcp, ukcpdtm,  ftr_startdate, ftr_enddate,
                            collection, domain, member,wind_hgt=10,toArrays=FALSE)


########################## Apply models to UKCP dataset ###############################################
#wsc<-biascorrectukcpone(ws1, ws2, ws3, rangelims = 1.1)
vars<-c("relhum","pres","lwrad","swrad","tmax","tmin","windspeed","prec")
vars %in% names(ftr_data)
vars %in% names(model_list)
#ukcp_vars<-c('hurs','psl','rls')

ftr_corrected<-list()
ftr_corrected$dtm<-ftr_data$dtm
ftr_corrected$tme<-ftr_data$tme
ftr_corrected$windheight_m<-ftr_data$windheight_m
ftr_corrected$tempheight_m<-ftr_data$tempheight_m

t0<-now()
for (n in 1:length(vars)) {
  v<-vars[n]
  #ukcpv<-ukcp_vars[n]
  print(v)
  if(v!='prec') ftr_corrected[[v]]<-biascorrect_apply(ftr_data[[v]], model_list[[v]], rangelims = 1.1) else  ftr_corrected[[v]]<-precipcorrect_apply(ftr_data[[v]],model_list[[v]])
  terra::time(ftr_corrected[[v]])<-ftr_corrected$tme
  names(ftr_corrected[[v]])<-ftr_corrected$tme
  # Write corrected ukcp data
  filename<-paste0(paste(v,rcp,collection, domain, datares, format(ftr_startdate,format="%Y%m%d"),sep="_"),"-",format(ftr_enddate,format="%Y%m%d"),".nc")
  print(filename)
  #writeCDF(ukcp_corrected,fileout=file.path(dir_ukcpcor,filename),overwrite=TRUE,compression=9,varname=ukcpv) # ,longname=longname,unit=unit
  print(now()-t0)
}
# Correction of 1 yr data took nearly 5mins
msk<-ifel(is.na(ftr_corrected$tmax)[[1]],NA,1)
plot(msk)

#lapply(ftr_data,global,fun=summary)
#lapply(ftr_corrected,summary)
plot_q_layers(mask(ftr_data$tmax,msk),vtext ='orig')
plot_q_layers(ftr_corrected$tmax,vtext ='corrected')
plot_q_layers(mask(ftr_data$tmin,msk),vtext ='orig')
plot_q_layers(ftr_corrected$tmin,vtext ='corrected')
plot_q_layers(mask(ftr_data$prec,msk),vtext ='orig')
plot_q_layers(ftr_corrected$prec,vtext ='corrected')


## Save corrected data
filename<-'ukcp18rcm_2026_bc_2017_2020_models.Rds'
write_climdata(ftr_corrected,filepath=file.path(dir_bcdata,filename),overwrite=TRUE)




