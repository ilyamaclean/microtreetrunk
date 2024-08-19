#' @title Generates plant area index profile
#' @description Generates a vector of length `n` of plausible plant area index values
#' @param PAI Total plant area index of canopy
#' @param skew number between 0 and 10 indicating the degree of skew towards top of
#' canopy in canopy foliage (see details)
#' @param spread positive non-zero number less than 100 indicating the degree of spread in
#' canopy foliage (see details)
#' @param n Number of plant area index values to generate. Default: 100 (see details)
#' @return a vector of length `n` of plant area index values the sum of which equals `PAI`
#' @details when specifying `skew`, lower numbers indicate greater skew towards top of
#' canopy (5 = symmetrical). In specifying `spread` a value of one indicates almost
#' all the foliage in concentrated in one canopy layer, whereas a value of 100 indicates
#' completely evenly spread.
#' @examples
#' pai <- PAIgeometry(3, 7, 70)
#' plot(pai, type = "l")
#' @rdname PAIgeometry
#' @export
PAIgeometry <- function(PAI, skew, spread, n = 100) {
  skew<-10-skew
  # Plant area index of canopy layer
  shape1<-100/spread
  x<-c(1:n)/(n+1)
  if (skew>5) {
    shape2<-(10-skew)/5+1
    shape2<-shape2/2*shape1
    y<-rev(dbeta(x,shape1,shape2))
  } else {
    shape2<-(skew+5)/5
    shape2<-shape2/2*shape1
    y<-dbeta(x,shape1,shape2)
  }
  y<-PAI/sum(y)*y
  y
}
#' @title Calculates solar position
#' @description Calculates the solar zenith and azimuth angles
#' @param tme POSIXlt object of times in UTC
#' @param lat latitude (decimal degrees)
#' @param long (decimal degrees, -ve west of Greenwich meridian)
#' @return a list of zenith and azimuth angles (decimal degrees)
#' @examples
#' # At noon at Porthleven, Cornwall on 21st Jun
#' tme <- as.POSIXlt(0, origin = "2023-06-21 12:00", tz = "UTC")
#' solarposition(tme, 50.08, -5.31)
#' # Hourly at Porthleven, Cornwall on 21st Jun
#' tme <- as.POSIXlt(c(0:23) * 3600, origin = "2023-06-21 00:00", tz = "UTC")
#' sp <- solarposition(tme, 50.08, -5.31)
#' sa <- 90 - sp$zenith
#' plot(sa ~ as.POSIXct(tme), pch = 19, cex = 4, col = "yellow",
#'      xlab = "Time", ylab = "Solar altitude")
#' @rdname solarposition
#' @export
solarposition<-function(tme,lat,long) {
  # Calculate Astronomical Julian day
  dd<-with(tme,mday+(hour+(min+sec/60)/60)/24)
  ma<-with(tme,mon+1+((mon+1)<3)*12)
  ya<-with(tme,year+1900+((mon+1)<3)*-1)
  jd<-trunc(365.25*(ya+4716))+trunc(30.6001*(ma+1))+dd-1524.5
  B<-(2-trunc(ya/100)+trunc(trunc(ya/100)/4))
  jd<-jd+(jd>2299160)*B
  # Calculate solar time
  lt<-with(tme,hour+(min+sec/60)/60)
  m<-6.24004077+0.01720197*(jd-2451545)
  eot<- -7.659*sin(m)+9.863*sin(2*m+3.5932)
  st<-lt+(4*long+eot)/60
  # Calculate solar zenith
  lat<-lat*pi/180
  tt<-0.261799*(st-12)
  dec<-(pi*23.5/180)*cos(2*pi*((jd-159.5)/365.25))
  coh<-sin(dec)*sin(lat)+cos(dec)*cos(lat)*cos(tt)
  z<-acos(coh)*(180/pi)
  # Calculate solar azimuth
  sh<-sin(dec)*sin(lat)+cos(dec)*cos(lat)*cos(tt)
  hh<-(atan(sh/sqrt(1-sh^2)))
  sazi<-cos(dec)*sin(tt)/cos(hh)
  cazi<-(sin(lat)*cos(dec)*cos(tt)-cos(lat)*sin(dec))/
    sqrt((cos(dec)*sin(tt))^2+(sin(lat)*cos(dec)*cos(tt)-cos(lat)*sin(dec))^2)
  sqt<-1-sazi^2
  sqt[sqt<0]<-0
  azi<-180+(180*atan(sazi/sqrt(sqt)))/pi
  azi[cazi<0 & sazi<0]<-180-azi[cazi<0 & sazi<0]
  azi[cazi<0 & sazi>=0]<-540-azi[cazi<0 & sazi>=0]
  solar<-list(zenith=z,azimuth=azi)
  return(solar)
}
#' @title Two-stream radiation
#' @description Calculates upward and downward shortwave direct and diffuse radiation
#' @param microin an object of class `microin` as returned by [modelin()]
#' @param i time step for wich radiation is required
#' @return a list of downward and upward diffuse and direct radiation streams.
#' @details Applies the variant of the Dickenson-Sellers two-stream model
#' described by Yuan et al (2017). A demonstration of its application is given
#' in the accompanying vignette.
#' @rdname canopyrad
#' @export
canopyrad<-function(microin,i) {
  # Calculate radiation parameters
  dfo<-microin@weather[i,]
  tp<-.canopyradp(microin@vegp,dfo,microin@groundp$gref)
  paia<-microin@vegp$paia
  pait<-microin@vegp$pait
  # Calculate radiation streams
  swrad<-microin@weather$swrad[i]
  if (swrad > 0) {
    difrad<-microin@weather$difrad[i]
    dirr<-swrad-difrad
    # Downward direct
    Rbdown<-with(tp,((1-Fb)*exp(-kd*paia)+Fb)*dirr)
    # Direct contribution to downward diffuse
    Rdbm<-with(tp,((1-Fd)*((p8/sig)*exp(-kd*paia)+p9*exp(-h*paia)+p10*exp(h*paia))+Fd)*dirr)
    Rdbm[Rdbm<0]<-0
    # Diffuse contribution to downward diffuse
    Rddm<-with(tp,((1-Fd)*(p3*exp(-h*paia)+p4*exp(h*pait))+Fd)*difrad)
    Rddown<-Rdbm+Rddm
    # Direct contribution to upward
    Rdbm<-with(tp,((1-Fd)*((p5/sig)*exp(-kd*paia)+p6*exp(-h*paia)+p7*exp(h*paia))+Fd)*dirr)
    # Diffuse contribution to upward
    Rddm<-with(tp,((1-Fd)*(p1*exp(-h*paia)+p2*exp(h*pait))+Fd)*difrad)
    Rdup<-Rdbm+Rddm
    Rbdown[Rbdown<0]<-0
    Rddown[Rddown<0]<-0
    Rdup[Rdup<0]<-0
  } else {
    Rbdown<-rep(0,length(paia))
    Rddown<-rep(0,length(paia))
    Rdup<-rep(0,length(paia))
  }
  return(list(Rbdown=Rbdown,Rddown=Rddown,Rdup=Rdup))
}
#' @title Wind profile coefficient
#' @description Creates a vector by which wind speed at the top of the canopy
#' is multiplied to derive the vertical wind-height profile.
#' @param microin an object of class `microin` as returned by [modelin()]
#' @return a vector of wind speed coefficients (0-1)
#' @details See accompanying vignette for application of this function.
#' @rdname canwind
#' @export
canwind <- function(microin) {
  # Calculate mixing length and attenuation coefficient
  paii<-microin@vegp$paii
  d<-with(microin@vegp,0.65*veghgt)
  zm<-with(microin@vegp,0.1*veghgt)
  l_m<-with(microin@vegp,(0.32*(veghgt-d))/log((veghgt-d)/zm)) # mixing length
  wpai<-with(microin@vegp,(paii/mean(paii))*sum(paii))
  at<-with(microin@vegp,(0.2*wpai*veghgt/l_m)^0.5)
  ws<-1
  zh<- -1/length(paii)
  for (i in 2:length(microin@vegp$paii)) ws[i]<-ws[i-1]*exp(at[i]*zh)
  ws<-rev(ws)
  ws
}
#' @title Calculates ground surface temperature
#' @description Applies NicheMapR model to derive ground surface temperatures
#' @param microin an object of class `microin` as returned by [modelin()]
#' @param soilparams a data.frame of soil parameters (see inbuilt dataset `soilparams`)
#' @return Adds soil surface temperatures to [modelin()]
#' @details See accompanying vignette for application of this function. Requires
#' NicheMapR (github.com/mrke/NicheMapR)
#' @rdname groundtemp
#' @import NicheMapR
#' @export
groundtemp<-function(microin,soilparams) {
  twostreamp<-.canopyradp(microin@vegp,microin@weather,microin@groundp$gref)
  micro<-.NicheMapRin(microin,soilparams,twostreamp)
  microut<-microclimate(micro)
  soil<-as.data.frame(microut$soil)
  dfo<-microin@weather
  dfo$soiltemp<-soil$D0cm
  microin@weather<-dfo
  microin@twostreamp<-twostreamp
  return(microin)
}
#' @title Creates microclimate model input
#' @description Creates an S4 object of class `microin` used as an input to
#' the microclimate model
#' @param weather a data,frame of hourly weather (see details).
#' @param prec a vector of daily precipitation (see details).
#' @param vegp a list of vegetation parameters (see details).
#' @param groundp a list of soil parameters (see details).
#' @param lat latitude (decimal degrees)
#' @param long longitude (decimal degrees, -ve west of Greenwich meridian)
#' @param reqhgt height above ground for which model should be run (m)
#' @param n number of vertical layers over which to run model (default 100)
#' @return an S4 object of class `microin` used as an input to the microclimate model
#' @details See accompanying vignette for application of this function. The inbuilt
#' datasets `weather`, `prec`, `vegp` and `groundp` give examples of the components
#' used for building the model input. User-created inputs should follow the same format
#' and units as in these inbuilt datasets.
#' @rdname modelin
#' @export
modelin<-function(weather,prec,vegp,groundp,lat,long,reqhgt,n = 100) {
  # Calculate pai in each layer
  vegp$paii<-with(vegp,PAIgeometry(pait,skew,spread,n=n))
  paia<-rev(cumsum(rev(vegp$paii)))
  z<-(c(1:n)/n)*vegp$veghgt
  sel<-which.min(abs(z-reqhgt))[1]
  vegp$paia<-paia[sel]
  vegp$z<-z
  # Calculate solar position
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  sp<-solarposition(tme,lat,long)
  weather$zenith<-sp$zenith
  weather$azimuth<-sp$azimuth
  weather$si<-.solarcoef(sp,groundp)
  # Create input datasets
  #setClass("microin",slots=c(weather="data.frame",prec="numeric",vegp="list",groundp="list",
  #                           loc="numeric",twostreamp="list",reqhgt="numeric"))
  microin<-new("microin",weather=weather,prec=prec,vegp=vegp,groundp=groundp,
               loc=c(long,lat),twostreamp=list(),reqhgt=reqhgt)
  return(microin)
}
#' @title spline interpolate weather inputs
#' @description Spline interpolates hourly weather inputs to shorter time intervals.
#' @param microin an object of class `microin` as returned by [modelin()]
#' @param secs time interval in seconds of spline-interpolated data. Default: 300
#' @return an object of class `microin` contained a data.frame of spline interpolated weather
#' variables and ground surface tmeperatures.
#' @details See accompanying vignette for application of this function.
#' @rdname weatherspline
#' @export
weatherspline<-function(microin,secs=300) {
  # Obs Time
  d<-dim(microin@weather)[1]
  n<-(d-1)*3600/secs+1
  tme<-as.numeric(as.POSIXlt(microin@weather$obs_time,tz="UTC"))
  tme<-spline(tme,n=n)$y
  tme<-as.POSIXlt(tme,origin="1970-01-01 00:00",tz="UTC")
  # Temperature
  temp<-with(microin@weather,spline(temp,n=n)$y)
  # Relative humidity
  ea<-with(microin@weather,.satvp(temp)*relhum/100)
  ea<-spline(ea,n=n)$y
  es<-.satvp(temp)
  rh<-(ea/es)*100
  rh[rh>100]<-100
  # Pressure
  pres<-with(microin@weather,spline(pres,n=n)$y)
  # Shortwave radiation
  csr<-.clearskyrad(microin@weather$zenith)
  csf<-microin@weather$swrad/csr
  csf[1]<-mean(csf[1:24],na.rm=T)
  csf[length(csf)]<-mean(csf[(d-23):d],na.rm=T)
  csf<-na.approx(csf)
  csf[csf<0]<-0
  csf[csf>1]<-1
  csf<-spline(csf,n=n)$y
  csf[csf<0]<-0
  csf[csf>1]<-1
  z<-with(microin@weather,spline(zenith,n=n)$y)
  csr<-.clearskyrad(z)
  swrad<-csf*csr
  # Diffuse radiation
  dfr<-with(microin@weather,difrad/swrad)
  dfr[is.na(dfr)]<-1
  dfr<-spline(dfr,n=n)$y
  difrad<-swrad*dfr
  # Downward longwave radiation
  skyem<-with(microin@weather,spline(skyem,n=n)$y)
  skyem<-.setlims(microin@weather$skyem,skyem)
  temp<-.setlims(microin@weather$temp,temp)
  lwdown<-skyem*5.67*10^-8*(temp+273.15)^4
  lwdown[lwdown<0]<-0
  # Wind speed
  windspeed<-with(microin@weather,spline(windspeed,n=n)$y)
  # azi
  azi<-with(microin@weather,spline(azimuth,n=n)$y)
  # si
  si<-with(microin@weather,spline(si,n=n)$y)
  # soiltemp
  soiltemp<-with(microin@weather,spline(soiltemp,n=n)$y)
  # Set limits
  rh<-.setlims(microin@weather$relhum,rh)
  pres<-.setlims(microin@weather$pres,pres)
  swrad<-.setlims(microin@weather$swrad,swrad)
  difrad<-.setlims(microin@weather$difrad,difrad)
  windspeed<-.setlims(microin@weather$windspeed,windspeed)
  soiltemp<-.setlims(microin@weather$soiltemp,soiltemp)
  z<-.setlims(microin@weather$z,z)
  azi<-.setlims(microin@weather$azi,azi)
  si<-.setlims(microin@weather$si,si)
  si[si<0]<-0
  dfo<-data.frame(obs_time=tme,temp=temp,relhum=rh,pres=pres,swrad=swrad,
                  difrad=difrad,lwdown=lwdown,windspeed=windspeed,soiltemp=soiltemp,
                  zenith=z,azimuth=azi,si=si)
  microin@weather<-dfo
  return(microin)
}
#' @title Create model inputs for first time-step
#' @description The microclimate model is run sequentially in time-increments, using
#' as an input, the output from the previous time-step. This function creates the
#' input for the first time-step
#' @param microin an object of class `microin` as returned by [modelin()]
#' @return a list containing vectors of air and leaf temperatures and
#' relative humidifies (see details)
#' @details Three vectors of length `n` are returned, the `ith` element of
#' which corresponds to a height (i / n) x canopy height above ground. An example
#' of the application of this function is shown in the accompanying vignette.
#' @rdname setinit
#' @export
setinit<-function(microin) {
  n<-length(microin@vegp$z)
  ta<-spline(c(microin@weather$soiltemp[1],microin@weather$temp[1]),n=n)$y
  miconeout<-list(tair=ta,tleaf=rep(microin@weather$temp[1],n),
                  rh=rep(microin@weather$relhum[1],n))

  return(miconeout)
}
#' @title Runs one step of microclimate model
#' @description This function runs one step of the microclimate model
#' @param microin an object of class `microin` as returned by [modelin()].
#' @param miconeout a list of model outputs fromthe previous time-step or as
#' returned by [setinit()]
#' @param ws a vector of wind coeffients as returned by [canwind()].
#' @param i which time-increment of the model to run. Default 1 (the first time-step).
#' @param iter Number of iterations over which Langrangian microclimate model
#' is run witin each time-step to derive stable values. Default: 5
#' @param plotout optional logical indicating whether to plot outputs. Default: TRUE
#' @param plotvar One of `tair` (air temperature, the default), `tleaf` (leaf
#' temperature) or `rh` (relative humidity). Ignored if `plotout = FALSE`.
#' @return a list containing vectors of air and leaf temperatures and
#' relative humidifies (see details)
#' @details Three vectors of length `n` are returned, the `ith` element of
#' which corresponds to a height (i / n) x canopy height above ground. An example
#' of the application of this function is shown in the accompanying vignette.
#' @rdname runmicroone
#' @export
runmicroone<-function(microin,miconeout,ws,i=1,iter=5,plotout=TRUE,plotvar="tair") {
  # Calculate H
  HT<-.calcH(microin,miconeout,ws,i)
  # Calculate top of canopy temperature and humidity
  h<-.CalcTh(microin,HT,i)
  # Run Langrangian model
  miconeout<-.LangrangianR(HT,microin,miconeout,h,R=iter)
  if (plotout) {
    z<-microin@vegp$z
    if (plotvar=="tair") {
      x<-miconeout$tair
      plot(z~x,type="l",xlab="Air temperature",ylab="Height (m)")
    }
    if (plotvar=="tleaf") {
      x<-miconeout$tleaf
      plot(z~x,type="l",xlab="Leaf temperature",ylab="Height (m)")
    }
    if (plotvar=="rh") {
      x<-miconeout$rh
      plot(z~x,type="l",xlab="Relative humidity",ylab="Height (m)")
    }
  }
  return(list(miconeout=miconeout,HT=HT))
}
#' @title Runs microclimate model
#' @description This function runs the full microclimate model.
#' @param microin an object of class `microin` as returned by [modelin()]
#' @param interval time increment in seconds over which to run model
#' @param iter Number of iterations over which Langrangian microclimate model
#' is run witin each time-step to derive stable values. Default: 5
#' @param tomax time step to which to run model. If `NA` model is run for entire
#' duration
#' @return an S4 object of class `microout` containing the folowing slots:
#' @return `Microclimate` a list containing vectors of air and foliage temperatures (deg C)
#' and relative humidifies (percentage) for each time increment.
#' @return `Radiation`  a list containing vectors of downard direct, downward diffuse,
#' upward diffuse and upward and downward longwave radiation (W/m^2)
#' @return `Inputweather` a data.frame of pline-interpolated input weather
#' variables used to run the model.
#' @details See accompanying vignette for applictaion of the model and a detailed
#' description of model outputs.
#' @rdname runmicro
#' @export
runmicro<-function(microin,interval=300,iter=5,tomax=NA) {
  cat("Computing ground surface temperatures using NicheMapR\n")
  # Calculate soil temperature and add to microin
  microin<-groundtemp(microin,soilparams)
  # Interpolate to correct interval (default 300 seconds / 5 mins)
  cat("Interpolating model inputs to correct interval\n")
  microin<-weatherspline(microin,secs=interval)
  if (is.na(tomax)) tomax<-dim(microin@weather)[1]
  # Get additional inputs for running model
  ws<-canwind(microin) # wind shelter coefficients
  miconeout<-setinit(microin) # initial values
  # sort out paia
  microin@vegp$paia<-rev(cumsum(rev(microin@vegp$paii)))
  # Create variables for return
  Rbdown<-0
  Rddown<-0
  Rdup<-0
  Lwdown<-0
  Lwup<-0
  Tair<-0
  RH<-0
  tleaf<-0
  # find out which entry to extract
  z<-microin@vegp$z
  dif<-abs(z-microin@reqhgt)
  sel<-which.min(dif)[1]
  cat("Initialising model\n")
  # Run model for 100 timesteps to initialise
  for (i in 1:100) {
    runm<-runmicroone(microin,miconeout,ws,i,iter,plotout=FALSE)
    miconeout<-runm$miconeout
  }
  # Run model
  cat("Running model\n")
  for (i in 1:tomax) {
    runm<-runmicroone(microin,miconeout,ws,i,iter,plotout=FALSE)
    miconeout<-runm$miconeout
    HT<-runm$HT
    Rbdown[i]<-HT$Rbdown[sel]
    Rddown[i]<-HT$Rddown[sel]
    Rdup[i]<-HT$Rdup[sel]
    Lwdown[i]<-HT$Lwdown[sel]
    Lwup[i]<-HT$Lwup[sel]
    Tair[i]<-miconeout$tair[sel]
    RH[i]<-miconeout$rh[sel]
    tleaf[i]<-mean(miconeout$tleaf)
  }
  # Create input datasets
  wsp<-ws[sel]*microin@weather$windspeed[1:tomax]
  Microclimate<-list(Tair=Tair,Tfoliage=tleaf,Relhum=RH, Windspeed=wsp)
  Radiation<-list(Rdirdown=Rbdown,Rdifdown=Rddown,Rdifup=Rdup,Rlwdown=Lwdown,Rlwup=Lwup)
  Inputdata<-microin@weather[1:tomax,]
  #setClass("microout",slots=c(Microclimate="list",Radiation="list",InputWeather="data.frame"))
  microout<-new("microout",Microclimate=Microclimate,Radiation=Radiation,InputWeather=Inputdata)
  return(microout)
}
#' @title Plot model results
#' @description This function is used plot time-series of model outputs
#' @param microout an object of class `microout` as returned by [runmicro()]
#' @param var one of `Tair` (air temperature), `Tfoliage` (foliage temperature),
#' `Relhum` (relative humidity), `Windspeed` (wind speed) `Rdirdown` (Downward direct radiation), `Rdifdown`
#' (downward diffuse radiation), `Rdifup` (upward shortwave radiation, assumed
#' entirely diffuse), `Rlwdown` (downward longwave radiation), Rlwup (upward
#' longwave radiation). Default: 'Tair'
#' @details See accompanying vignette for application of this function.
#' @rdname plotmicro
#' @export
plotmicro<-function(microout,var="Tair") {
  tme<-as.POSIXlt(microout@InputWeather$obs_time,tz="UTC")
  tme<-as.POSIXct(tme)
  if (var=="Tair") x<-microout@Microclimate$Tair
  if (var=="Tfoliage") x<-microout@Microclimate$Tfoliage
  if (var=="Relhum") x<-microout@Microclimate$Relhum
  if (var=="Windspeed") x<-microout@Microclimate$Windspeed
  if (var=="Rdirdown") x<-microout@Radiation$Rdirdown
  if (var=="Rdifdown") x<-microout@Radiation$Rdifdown
  if (var=="Rdifup") x<-microout@Radiation$Rdifup
  if (var=="Rlwdown") x<-microout@Radiation$Rlwdown
  if (var=="Rlwup") x<-microout@Radiation$Rlwup
  plot(x~tme,type="l",xlab="Time",ylab=var)
}
