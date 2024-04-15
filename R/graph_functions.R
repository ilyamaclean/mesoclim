#' Plot timeseries plots of spatial mean, min and max values by various timesteps
#'
#' @param r stack of spatrasters with time values
#' @param var variable name used in plot titles
#' @param idx time index indicating whether summary plots by year, month, week, day of year etc
#' @param lgd whether to include legend
#'
#' @return NA plots series of timeploat
#' @export
#' @import magrittr
#'
#' @examples
#' plot_timestats_r(ukcp18rcm$climarray$tmax,'tmax',idx='doy',lgd=FALSE)
plot_timestats_r<-function(r,var,idx=c('years', 'months', 'week',  'doy', 'yearmonths', 'yearweeks', '7days','hour'),lgd=FALSE){
  if(idx=='hour'){
    time_mean<- global(r,"mean", na.rm=TRUE)
    time_max<-global(r, "max", na.rm=TRUE)
    time_min<-global(r, "min", na.rm=TRUE)
    plot_df<-as.data.frame(cbind(tstep=terra::time(r),mean=time_mean$mean,max=time_max$max,min=time_min$min))
    plot_df<-plot_df[order(plot_df$tstep),]
  } else{
  time_mean<-tapp(r,index=idx,fun=mean) %>% global("mean", na.rm=TRUE)
  time_max<-tapp(r,index=idx,fun=max) %>% global("max", na.rm=TRUE)
  time_min<-tapp(r,index=idx,fun=min) %>% global("min", na.rm=TRUE)
  plot_df<-as.data.frame(cbind(tstep=as.numeric(sapply(strsplit(rownames(time_mean),'_'),tail,1)),mean=time_mean$mean,max=time_max$max,min=time_min$min))
  plot_df<-plot_df[order(plot_df$tstep),]
  }
  matplot(plot_df$tstep, plot_df[,2:4], type = "l", lty = 1,
          col = c("green", "red", "blue"), xlab = idx, ylab = var, font.main = 1,
          tck = 0.02, cex.main=1, cex.axis=0.7, main = paste(var,'by',idx), cex.main=1)
  if(lgd==TRUE) legend("topright", legend = c("Mean", "Max", "Min"), cex=0.5,
         col = c("green", "red", "blue"),
         lty = 1)
}

#' Plots chosen layers from raster stack corresponding to selected quantiles of spatial mean
#'
#' @param r Spatraster stack (time and name values used in plot titles)
#' @param p vector of quantile probabilities to plot (0:1)
#' @param fun = mean or other spatial function used to summaries layers for which quantiles dtermined
#'
#' @return plots series of raster plots equal to the length of p
#' @export
#'
#' @examples
plot_q_layers<-function(r,p=c(0, 0.5, 1),fun='mean', common_scale=FALSE){
  par(mfrow=c(1,3), mai=c(1,0.1,0.1,0.1))
  # Find corresponding layers
  lyrstat<-global(r,fun,na.rm=TRUE)
  qtls<-quantile(lyrstat[,1],prob=p)
  sel<-c()
  for(q in qtls) sel<-c(sel,which.min(abs(lyrstat[,1] - q)))
  # Plot
  out_r<-r[[sel]]
  terra::time(out_r)<-terra::time(r[[sel]])
  names(out_r)<-paste(names(r[[sel]]),'(q=',names(qtls),')',as.POSIXlt(terra::time(out_r)))
  if(common_scale){
    rng=range(global(out_r,'range',na.rm=TRUE))
    plot(out_r,main=names(out_r),font.main=1, buffer=FALSE,nc=3,nr=1, cex.main=1, range=rng, nc=length(p))
  } else {
    plot(out_r,main=names(out_r),font.main=1, cex.main=1, nc=length(p))
  }
}

.dir_to_cardinal <- function(x) {
  upper <- seq(from = 11.25, by = 22.5, length.out = 17)
  cardinals <- c('N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW', 'N')
  ifelse(x>360 | x<0,NA,cardinals[findInterval(x,upper,rightmost.closed = T)+1])
}

#' Plot 'star' or 'radar' chart of wind cardinal directions
#'
#' @param winddir wind direction 3D Spatraster
#' @param windspeed windspeed 3D Spatraster of same geometry as winddir]
#'
#' @return plots star charts showing frequency of wind directions and of
#' average and maximum windspeeds by direction
#' @export
#' @import fmsb
#'
#' @examples
plot_wind<-function(winddir,windspeed){
  # Classify by cardinal group
  dgrp<- .dir_to_cardinal(c(winddir))
  dgrp_mean<-tapply(c(windspeed),dgrp,mean)
  dgrp_max<-tapply(c(windspeed),dgrp,max)
  dgrp_freq<-as.array(table(dgrp))

  # Mean/Max windspeed plot by dir
  cardorder<- c('N', 'NNW', 'NW', 'WNW', 'W', 'WSW', 'SW', 'SSW', 'S', 'SSE', 'SE', 'ESE', 'E', 'ENE', 'NE', 'NNE')
  df<-as.data.frame(matrix(0, 4, 16) )
  names(df)<-cardorder
  rownames(df)<-c('Max','Min','Mean windspeed','Max windspeed')
  df[1,]<-rep(max(dgrp_max),length(cardorder)) # max speed across  all directions
  df[2,]<-rep(0,length(cardorder)) # min speed across all dir
  df[3,names(dgrp_mean)]<-dgrp_mean
  df[4,names(dgrp_max)]<-dgrp_max

  par(mar=c(1,1,1,1))
  par(mfrow=c(1,2))
  fmsb::radarchart(df, title='Mean and max windspeeds', cex=0.75)
  legend(
    x = "bottomleft", legend = c('Mean','Maximum'), horiz = TRUE,
    bty = "n", pch = 20 , col = c("black", "red"),
    text.col = "black", cex = 1, pt.cex = 1
  )

  # Freq plot by dir
  df<-as.data.frame(matrix(0, 3, 16) )
  names(df)<-cardorder
  rownames(df)<-c('Max','Min','Frequency')
  df[1,]<-rep(max(dgrp_freq),length(cardorder)) # max speed across  all directions
  df[2,]<-rep(0,length(cardorder)) # min speed across all dir
  df[3,names(dgrp_freq)]<-dgrp_freq

  fmsb::radarchart(df, title='Wind direction frequency',pcol=c('green'), cex=0.75)

  legend(
    x = "bottomleft", legend = c('Frequency'), horiz = TRUE,
    bty = "n", pch = 20 , col = c("green"),
    text.col = "black", cex = 1, pt.cex = 1
  )

}




