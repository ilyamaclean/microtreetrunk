#' internal function for setting limits when splining weather data
.setlims<-function(xs,xl) {
  mx<-max(xs)
  mn<-min(xs)
  xl[xl>mx]<-mx
  xl[xl<mn]<-mn
  xl
}
#' Internal function for calaculating time-step
.tstep<-function(InputWeather,i) {
  if (i > 1) {
    x<-with(InputWeather,difftime(obs_time[i],obs_time[i-1],units="secs"))
  } else x<-with(InputWeather,difftime(obs_time[i+1],obs_time[i],units="secs"))
  as.numeric(x)
}
#' Internal function = temperatures derived using Penman-Monteith equation
.penman<-function(modin, tsu, fwet) {
  # Find which segment is closest to the solar azimuth
  sazi<-modin@InputWeather$azimuth%%360
  asps<-modin@FixedParams$asps
  asps2<-c(asps,360)
  aspsm<-matrix(rep(asps2,each=length(sazi)),ncol=length(asps2),nrow=length(sazi))
  sazim<-matrix(rep(sazi,length(asps2)),ncol=length(asps2),nrow=length(sazi))
  xx<-abs(aspsm-sazim)
  seg<-apply(xx,1,which.min)
  seg[seg>length(asps)]<-1
  # Absorbed radiation
  salt<-90-modin@InputWeather$zenith
  dni<-modin@Radiation$Rdirdown/cos(modin@InputWeather$zenith*pi/180)
  dif<-modin@Radiation$Rdifdown+modin@Radiation$Rdifup
  index<-cos(salt*pi/180)*cos((sazi-asps[seg])*pi/180)
  index[index<0]<-0
  dirr<-index*dni
  dirr[dirr<0]<-0
  dirr[dirr>1352]<-1352
  swabs<-(1-modin@Treeparams$tref)*(dirr+dif)
  lwabs<-(1-modin@Treeparams$tem)+with(modin@Radiation,Rlwdown+Rlwup)
  lwabs[lwabs<0]<-0
  Rabs<-(swabs+lwabs)/2
  # forced convection
  d<-with(modin@Treeparams,2*radius)
  gHfo<-0.135*sqrt(modin@Microclimate$Windspeed/d)
  # free convection
  dT<-tsu-modin@Microclimate$Tair
  gHfr<-0.033*(abs(dT)/d)^0.25
  gHa<-pmax(gHfo,gHfr)
  # Radiative conductance
  Tt<-(modin@Microclimate$Tair+tsu)/2
  gr<-4*modin@Treeparams$tem*5.67*10^-8*(Tt+273.15)^3
  gHr<-gHa+gr
  # Slope of saturated vapour pressure curve
  De<-(4098*.satvp(Tt))/(Tt+237.3)^2
  DD<-with(modin@Microclimate,.satvp(Tair)*(1-Relhum/100))
  s<-De/modin@InputWeather$pres
  # Penman-Monteith equation
  lambda<-44526
  Rnet<-Rabs-modin@Treeparams$tem*5.67*10^-8*(modin@Microclimate$Tair+273.15)^4
  top<-Rnet-fwet*gHa*(DD/modin@InputWeather$pres)
  btm<-29.3*gHr+44526*fwet*s*gHa
  dT<-top/btm
  Ts<-modin@Microclimate$Tair+dT
  Ts
}
#' Internal function for calculating radiation absorbed
.trunkradabs<-function(modin,i) {
  # Extract variables
  salt<-90-modin@InputWeather$zenith[i]
  sazi<-modin@InputWeather$azimuth[i]
  asps<-modin@FixedParams$asps
  dni<-modin@Radiation$Rdirdown[i]/cos(modin@InputWeather$zenith[i]*pi/180)
  dif<-modin@Radiation$Rdifdown[i]+modin@Radiation$Rdifup[i]
  # Direct radiation absorbption
  index<-cos(salt*pi/180)*cos((sazi-asps)*pi/180)
  index[index<0]<-0
  # Absorbed shortwave radiation
  swabs<-(1-modin@Treeparams$tref)*(index*dni+dif) # original
  # swabs<-(1-modin@Treeparams$tref)*(index*dni+ 0.5*dif) # changed
  lwabs<-(1-modin@Treeparams$tem)+with(modin@Radiation,Rlwdown[i]+Rlwdown[i])
  # lwabs<-(1-modin@Treeparams$tem)+with(modin@Radiation,0.5*(Rlwdown[i]+Rlwdown[i])) # changed
  rabs<-swabs+lwabs
  rabs
}
#' Internal function for calculating radiation emitted
.treeradem<-function(modin,i) {
  Rem<-modin@Treeparams$tem*5.67*10^-8*(modin@treetemps[1,]+273.15)^4
  Rem
}
#' Internal function for calculating sensible heat
.treesensible<-function(modin,i) {
  # forced convection
  d<-with(modin@Treeparams,2*radius)
  gHfo<-rep(0.135*sqrt(modin@Microclimate$Windspeed[i]/d),modin@Treeparams$nsegs)
  # free convection
  dT<-modin@treetemps[1,]-modin@Microclimate$Tair[i]
  gHfr<-0.033*(abs(dT)/d)^0.25
  gHa<-pmax(gHfo,gHfr)
  # Sensible heat per unit area
  H<-gHa*29.3*(dT)
  H
}
#' Internal function for calculating latent heat
.treelatent<-function(modin,i,fwet) {
  if (fwet>0) {
    # forced convection
    d<-with(modin@Treeparams,2*radius)
    gHfo<-rep(0.135*sqrt(modin@Microclimate$Windspeed[i]/d),modin@Treeparams$nsegs)
    # free convection
    dT<-modin@treetemps[1,]-modin@Microclimate$Tair[i]
    gHfr<-0.033*(abs(dT)/d)^0.25
    gHa<-pmax(gHfo,gHfr)
    # Latent heat per unit area
    ea<-with(modin@Microclimate,0.6108*exp((17.27*Tair[i])/(Tair[i]+237.3))*(Relhum[i]/100))
    es<-0.6108*exp((17.27*modin@treetemps[1,])/(modin@treetemps[1,]+237.3))
    L<-44398*gHa/modin@InputWeather$pres[i]*(es-ea)*fwet
    L[L<0]<-0
  } else L<-rep(0,modin@Treeparams$nsegs)
  L
}
#' Internal function for calculating temperature of outer layer
.treeouter<-function(modin,EB,i) {
  # Adjust to Joules
  sarea<-with(modin@Treeparams,(2*pi*radius)/nsegs)
  tstep<-.tstep(modin@InputWeather,i)
  k<-EB/(abs(modin@Microclimate$Tair[i]-modin@treetemps[1,]))
  EB<-sarea*tstep*EB
  # Ensure covergence in temperatures if time-step is too long
  mu<-exp(-k*tstep/modin@FixedParams$scph[1,])
  mu<-1
  # Temperature due to energy fluxes at surface
  to<-modin@treetemps[1,]+(EB*mu)/modin@FixedParams$scph[1,]
  to
}
#' Internal function for calculating layer heat exchange
.layerexchange<-function(modin,i) {
  tstep<-.tstep(modin@InputWeather,i)
  tm<-modin@treetemps
  Kl<-modin@FixedParams$Kl
  scph<-modin@FixedParams$scph
  for (lyr in 2:modin@Treeparams$nlayers) {
    mu<-exp(-(Kl[lyr-1]*tstep)/(scph[lyr-1,]+scph[lyr,])/2)
    Hl<-Kl[lyr-1]*tstep*(tm[lyr,]-tm[lyr-1,])*mu # Flux from lyr to outer (-ve if lyr colder than outer)
    tm[lyr-1,]<-tm[lyr-1,]+Hl/scph[lyr-1,]  # Temperature adjust outer
    tm[lyr,]<-tm[lyr,]-Hl/scph[lyr,] # Temperature adjust lyr for H1
  }
  tm
}
#' Internal function for calculating segment heat exchange
.segmentexchange<-function(modin,i) {
  tstep<-.tstep(modin@InputWeather,i)
  tm<-modin@treetemps
  Ks<-modin@FixedParams$Ks
  scph<-modin@FixedParams$scph
  for (lyr in 1:modin@Treeparams$nlayers) {
    to<-tm[lyr,]
    tp<-c(to[2:length(to)],to[1])
    tn<-c(to[length(to)],to[1:(length(to)-1)])
    mu<-exp(-Ks[lyr]*tstep/scph[lyr,])
    # HS<-Ks[lyr]*(to-tp)+Ks[lyr]*(to-tn)*mu # temp diverges?
    HS<-Ks[lyr]*(tp-to)+Ks[lyr]*(tn-to)*mu # makes temp converge between adjacent segs
    tm[lyr,]<-to+(HS*tstep)/scph[lyr,]
  }
  tm
}
#' internal function for smoothing segment heats to avoid temperature yo-yo
.smoothsegments<-function(modin,i) {
  tstep<-.tstep(modin@InputWeather,i)
  tm<-modin@treetemps
  for (lyr in 3:modin@Treeparams$nlayers) {
    to<-tm[lyr,]
    tp<-c(to[2:length(to)],to[1])
    tn<-c(to[length(to)],to[1:(length(to)-1)])
    tm[lyr,] <-(to+tp+tn)/3
  }
  tm
}
#' Internal function for setting limits to tree trunk temperatures
.limtm<-function(tm,tmin,tmax) {
  tm[tm<tmin]<-tmin
  tm[tm>tmax]<-tmax
  tm
}
