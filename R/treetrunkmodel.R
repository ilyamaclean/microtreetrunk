#' @title Creates tree trunk model input
#' @description Creates an S4 object of class `treein` used as an input to
#' the tree trunk model
#' @param microout an object of class `microout` as returned by [runmicro()].
#' @param treeparams a list of tree parameters required for running the tree trunk model (see details).
#' @return an S4 object of class `treein` used as an input to the tree trunk model
#' @details See accompanying vignette for application of this function. The inbuilt
#' dataset `treeparams` give an examples of the components used for building the
#' model input. User-created inputs should follow the same format and units as in these inbuilt datasets.
#' @rdname modelprep
#' @export
modelprep<-function(microout,treeparams) {
  # ======================================= #
  # ~~~~~~~~~~~ Check inputs ~~~~~~~~~~~~~~ #
  # ======================================= #
  nlayers<-treeparams$nlayers
  nsegs<-treeparams$nsegs
  if (nlayers < 3) stop("Three or more layers needed")
  if (nsegs < 3) stop("Three or more segments needed")
  with(treeparams, if (length(K) == 1) K<-rep(K,nlayers))
  with(treeparams, if (length(cph) == 1) cph<-rep(cph,nlayers))
  with(treeparams, if (length(K) != nlayers) stop("K must be a single value or have length nlayers"))
  with(treeparams, if (length(cph) != nlayers) stop("cph] must be a single value or have length nlayers"))
  tz<-microout@Microclimate$Tair[1]
  tm<-matrix(tz,ncol=nsegs,nrow=nlayers)
  # ======================================= #
  # ~~~~ Time invariant calculations ~~~~~~ #
  # ======================================= #
  # Aspect of each segment
  asps<-c(0:(nsegs-1))*360/nsegs
  # Actual specific heat capacity of each segment and layer
  nds<-with(treeparams,radius-c(0:nlayers)^2*radius/nlayers^2)  # nodes distance from centre (m)
  carea<-pi*nds^2 # circle area (m^2)
  larea<-carea[1:nlayers]-carea[2:length(carea)] # layer volume per unit vertical area (m^3)
  lcph<-larea*treeparams$cph*10^6  # layer specific heat capacity (J/K)
  scph<-matrix(rep(lcph/nsegs,nsegs),ncol=nsegs) # segment specific heat capacity (J/K)
  # Actual thermal conductivities between each layer (in each segment)
  ndsc<-(nds[1:nlayers]+nds[2:length(nds)])/2  # distance form centre og each layer centre
  # distance from centre to edge
  z1<-ndsc-nds[2:length(nds)] # distance centre to edge
  z1<-z1[-nlayers]
  z2<-nds[2:nlayers]-ndsc[2:length(ndsc)]  # distance from edge to centre
  # Condictivity between layers per m2 cross-sectional area (W/K):
  Kl<-with(treeparams,K[1:(nlayers-1)]/z1+K[2:nlayers]/z2)
  circ<-2*pi*nds  # circumferences of circles
  csa<-circ[2:length(circ)]/nsegs # cross sectional area between circles
  Kl<-Kl*csa[1:(nlayers-1)]  # # Condictivity between layers (W/K)
  # Actual thermal conductivities between each segment (in each layer)
  circ<-2*pi*ndsc  # circumferences of cirles at node cenre (m)
  Ks<-with(treeparams,K/(circ/nsegs)) # Thermal conductivity between each segment per m2 cross-sectional area (W/K)
  lt<-nds[1:nlayers]-nds[2:length(nds)]
  Ks<-Ks*lt # Thermal conductivity between each segment (W/K)
  FixedParams<-list(asps=asps,scph=scph,Kl=Kl,Ks=Ks)
  setClass("treein",slots=c(Microclimate="list",Radiation="list",InputWeather="data.frame",
                            Treeparams="list",FixedParams="list",treetemps="matrix"))
  treein<-new("treein",Microclimate=microout@Microclimate,Radiation=microout@Radiation,
              InputWeather=microout@InputWeather,Treeparams=treeparams,
              FixedParams=FixedParams,treetemps=tm)
  return(treein)
}
#' @title Spline interpolates microclimate input
#' @description Spline interpolates microclimate input to tree trunk model to shorter time intervals.
#' @param modin an object of class `treein` as returned by [modelprep()].
#' @param secs time-interval (seconds) of required output data
#' @return an S4 object of class `treein` with higher temperal resolution data
#' @details See accompanying vignette for application of this function.
#' @rdname microspline
#' @export
microspline<-function(modin,secs=60) {
  # Get n
  tst<-.tstep(modin@InputWeather,2)
  d<-dim(modin@InputWeather)[1]
  n<-(d-1)*tst/secs+1
  # Spine Microclimate
  Tair<-spline(modin@Microclimate$Tair,n=n)$y
  Tfoliage<-spline(modin@Microclimate$Tfoliage,n=n)$y
  Relhum<-spline(modin@Microclimate$Relhum,n=n)$y
  Windspeed<-spline(modin@Microclimate$Windspeed,n=n)$y
  # Set limits
  Tair<-.setlims(modin@Microclimate$Tair,Tair)
  Tfoliage<-.setlims(modin@Microclimate$Tfoliage,Tfoliage)
  Relhum<-.setlims(modin@Microclimate$Relhum,Relhum)
  Windspeed<-.setlims(modin@Microclimate$Windspeed,Windspeed)
  Microclimate<-list(Tair=Tair,Tfoliage=Tfoliage,Relhum=Relhum,Windspeed=Windspeed)
  # Spine Radiation
  Rdirdown<-spline(modin@Radiation$Rdirdown,n=n)$y
  Rdifdown<-spline(modin@Radiation$Rdifdown,n=n)$y
  Rdifup<-spline(modin@Radiation$Rdifup,n=n)$y
  Rlwdown<-spline(modin@Radiation$Rlwdown,n=n)$y
  Rlwup<-spline(modin@Radiation$Rlwup,n=n)$y
  # Set limits
  Rdirdown<-.setlims(modin@Radiation$Rdirdown,Rdirdown)
  Rdifdown<-.setlims(modin@Radiation$Rdifdown,Rdifdown)
  Rdifup<-.setlims(modin@Radiation$Rdifup,Rdifup)
  Rlwdown<-.setlims(modin@Radiation$Rlwdown,Rlwdown)
  Rlwup<-.setlims(modin@Radiation$Rlwup,Rlwup)
  Radiation<-list(Rdirdown=Rdirdown,Rdifdown=Rdifdown,Rdifup=Rdifup,
                  Rlwdown=Rlwdown,Rlwup=Rlwup)
  # InputWeather
  modin@InputWeather$si[modin@InputWeather$si<0]<-0
  InputWeather<-matrix(0,ncol=dim(modin@InputWeather)[2],nrow=n)
  tme<-as.numeric(as.POSIXlt(modin@InputWeather$obs_time))
  tme<-spline(tme,n=n)$y
  tme<-as.POSIXlt(tme,origin="1970-01-01 00:00",tz="UTC")
  for (i in 2:dim(modin@InputWeather)[2]) {
    InputWeather[,i]<-spline(modin@InputWeather[,i],n=n)$y
  }
  # Set Limits
  for (i in 2:dim(modin@InputWeather)[2]) {
    InputWeather[,i]<-.setlims(modin@InputWeather[,i],InputWeather[,i])
  }
  # Convert to data.frame
  InputWeather<-as.data.frame(InputWeather)
  names(InputWeather)<-names(modin@InputWeather)
  InputWeather$obs_time<-tme
  # Assign to modin
  modin@Microclimate<-Microclimate
  modin@Radiation<-Radiation
  modin@InputWeather<-InputWeather
  return(modin)
}
#' @title Runs one step of tree trunk model
#' @description Runs a single step of the tree trunk model
#' @param modin an object of class `treein` as returned by [modelprep()].
#' @param i time-step of model to run.
#' @param fwet fraction of tree trunk surface acting like a saturated water surface.
#' @return an S4 object of class `treein` used as an input to the next iteration of the tree trunk model
#' @details See accompanying vignette for application of this function.
#' @rdname treeonestep
#' @export
treeonestep<-function(modin,i,fwet=0.1) {
  # Calculate radiation absorbed by outer layer
  Rabs<-.trunkradabs(modin,i)*0.5
  # Calculate radiation emitted by outer layer
  Rem<-.treeradem(modin,i)
  # Calculate sensible heat
  H<-.treesensible(modin,i)
  # Calculate latent heat
  L<-.treelatent(modin,i,fwet)
  # Calculate outer layer temperature
  EB<-Rabs-Rem-H-L
  modin@treetemps[1,]<-.treeouter(modin,EB,i)
  # Calculate layer exchange
  modin@treetemps<-.layerexchange(modin,i)
  # Calculate segment exchange
  modin@treetemps<-.segmentexchange(modin,i)
  #modin@treetemps<-.smoothsegments(modin,i)
  return(modin)
}
#' @title Runs tree trunk model
#' @description Runs a single step of the tree trunk model
#' @param modin an object of class `treein` as returned by [modelprep()].
#' @param fwet fraction of tree trunk surface acting like a saturated water surface.
#' @param tomax time step to which to run model. If `NA` model is run for entire
#' duration
#' @param daysinit number of days to use as a burn-in when starting the model (see details).
#' @return an object of class `treemodout` - a list of the following: (i) an array of temperatures for each layer,
#' segment and time-interval; (ii) a POSIXlt object of the corresponding times
#' for each time-interval.
#' @details this function calls [treeonestep()] repeatedly for each time-increment.
#' To ensure plausible starting values, temperatures in each segment and layer are
#' set to air temperature and the model run for `daysinit` days to introduce a
#' realistic temperature profile. The model is then run from the beginning. See accompanying
#' vignette for a demonstration of this function.
#' @rdname runtreetrunk
#' @export
runtreetrunk<-function(modin,fwet=0.1,tomax=NA,daysinit=1) {
  n<-length(modin@Microclimate$Tair)
  if (is.na(tomax)) tomax<-n
  # Calculate temperature limits
  TPM<-.penman(modin,modin@Microclimate$Tair,fwet=0.1)
  Tmx<-TPM+5
  ea<-with(modin@Microclimate,.satvp(Tair)*Relhum/100)
  TD<-.dewpoint(ea,modin@Microclimate$Tair)
  Tmn<-TD-5
  # Initialise model
  id<-(24*3600*daysinit)/.tstep(modin@InputWeather,2)# time-steps in first day
  for (i in 1:id) {
    modin<-treeonestep(modin,i,fwet)
    modin@treetemps<-.limtm(modin@treetemps,Tmn[i],Tmx[i])
  }
  Tout<-array(0,dim=c(dim(modin@treetemps),tomax))
  for (i in 1:tomax) {
    modin<-treeonestep(modin,i,fwet)
    modin@treetemps<-.limtm(modin@treetemps,Tmn[i],Tmx[i])
    Tout[,,i]<-modin@treetemps
  }
  treemodout<-list(Tout=Tout,
                   tme=as.POSIXlt(modin@InputWeather$obs_time[1:tomax],tz="UTC"))
  class(treemodout)<-"treemodout"
  return(treemodout)
}
#' @title plots a time-series of tree-trunk temperatures
#' @description plots a time-series of temperatures of the outer layer of the tree-trunk
#' @param treemodout an object of class `treemodout` as returned by [runtreetrunk()].
#' @param fun one of min, max or mean (see details).
#' @param col color of plot line
#' @details If `fun = mean`, within each time-step, the mean temperature of the outer ring
#' is plotted. If `fun = min` or `fun = max`, within each time-step, the temperatures in
#' the coldest or warmest segment are selected for plotting. See accompanying vignette for
#' a demonstration of this function.
#' @rdname timeseriesplot
#' @export
timeseriesplot<-function(treemodout,fun=max,col="black") {
  # Select outer ring
  Tm<-treemodout$Tout[1,,]
  tx<-apply(Tm,2,fun)
  par(mar=c(6,6,3,3))
  plot(tx~as.POSIXct(treemodout$tme),type="l",cex.axis=2,cex.lab=2,lwd=1,
       xlab = "Time", ylab = expression("Temperature ("*~degree*C*")"),
       col=col)
}
#' @title plots a cross-section of tree trunk temperatures
#' @description for a specified time period, plots a raster representing the
#' cross-sectional temperatures of the tree-trunk.
#' time-series of temperatures of the outer layer of the tree-trunk
#' @param modin an object of class `treein` as returned by [modelprep()].
#' @param modout an object of class `treemodout` as returned by [runtreetrunk()].
#' @param i the time-increment for which to plot a cross-section of temperatures.
#' @details Requires the packages `fields` and `raster` to be installed as
#' temperatures within each segement and layer, as returned by [runtreetrunk()]
#' are interpolated to a 100 x 100 pixel resolution using the `fields` prior to
#' plotting. See accompanying vignette for a demonstration of this function..
#' @rdname crossectionplot
#' @export
crossectionplot<-function(modin,modout,i,zlim=c(NA,NA)) {
  asps<-modin@FixedParams$asps
  tm<-modout$Tout[,,i]
  radius<-modin@Treeparams$radius
  nlayers<-modin@Treeparams$nlayers
  nds<-radius-c(0:nlayers)^2*radius/nlayers^2  # nodes distance from centre (m)
  ndsc<-(nds[1:nlayers]+nds[2:length(nds)])/2  # dis
  ys<-tm*0
  xs<-ys
  for (s in 1:dim(tm)[2]) {
    xs[,s]<-ndsc*sin(asps[s]*pi/180)
    ys[,s]<-ndsc*cos(asps[s]*pi/180)
  }
  # Create temperature matrix
  xs<-as.vector(xs)
  ys<-as.vector(ys)
  tc<-as.vector(tm)
  xy<-data.frame(x=xs,y=ys)
  tps <- suppressWarnings(Tps(xy, tc))
  m<-matrix(NA,ncol=100,nrow=100)
  r<-raster(m)
  extent(r)<-c(-radius,radius,-radius,radius)
  r<-suppressWarnings(interpolate(r,tps))
  # Create mask
  xys<-as.data.frame(xyFromCell(r, c(1:10000)))
  xys$h<-sqrt(xys$x^2+xys$y^2)
  xys$z<-1
  xys$z[xys$h>radius]<-NA
  xys$h<-NULL
  msk<-rasterFromXYZ(xys)
  # Mask and plot
  r<-mask(r,msk)
  v<-getValues(r)
  if (is.na(zlim)[1]) zlim=c(min(v,na.rm=T),max(v,na.rm=T))
  mypal <- colorRampPalette(c("darkblue", "blue", "green", "yellow", "orange",
                              "red"))(255)
  tme<-modout$tme[i]
  plot(r,col=mypal,zlim=zlim,main=tme)
}
