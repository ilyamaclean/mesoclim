
# Crop mesoclimate outputs
daily100m<-lapply(mesoclimate,function(x) if(class(x)[1]=="SpatRaster") wrap(crop(x,e)) else x)
lapply(daily100m,class)
usethis::use_data(daily100m)

