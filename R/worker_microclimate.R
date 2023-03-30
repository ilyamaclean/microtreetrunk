#' Internal function for calculating solar coefficient
.solarcoef<-function(solar,groundp) {
  z<-solar$zenith*pi/180
  a<-solar$azimuth*pi/180
  sl<-groundp$slope*pi/180
  ap<-groundp$aspect*pi/180
  index<-cos(z)*cos(sl)+sin(z)*sin(sl)*cos(a-ap)
  index[index<0]<-0
  index
}
#' Internal function for calculating canopy extinction coefficient
.cank<-function(z,x) {
  z[z>90]<-90
  Z<-z*pi/180
  if (x==1) {
    k<-1/(2*cos(Z))
  } else if (x==0) {
    k<-2/(pi*tan(0.5*pi-Z))
  } else if (is.infinite(x)) {
    k<-1
  } else {
    k<-sqrt((x^2+(tan(Z)^2)))/(x+1.774*(x+1.182)^(-0.733))
  }
  k
}
#' Internal function for calculating clear sky radiation
.clearskyrad <- function(zenith) {
  z<-zenith*pi/180
  # Calaculate clear sky radiation
  m<-35*cos(z)*((1224*cos(z)^2+1)^(-0.5))
  TrTpg<-1.021-0.084*(m*0.961337+0.051)^0.5
  Tw<-1-0.077*(0.4192448*m)^0.3
  od<-TrTpg*Tw*0.935*m
  Ic<-1352.778*(cos(z))*od
  Ic[Ic>1352.778]<-1352.778
  Ic[is.na(Ic)]<-0
  Ic
}
#' Internal function for calculating saturated vapor pressure
.satvp<-function(tc) {
  0.6108*exp((17.27*tc)/(tc+237.3))
}
#' Internal function for calculating dewpoint
.dewpoint <- function(ea, tc = 11) {
  e0 <- 611.2/1000
  L <- (2.501*10^6) - (2340 * tc)
  T0 <- 273.15
  Rv <- 461.5
  it <- 1/T0 - (Rv/L) * log(ea/e0)
  Tdew <- 1/it - 273.15
  Tdew
}
#' Internal function for calculate parameters of two stream canopy radiation model
.canopyradp<-function(vegp,weather,gref) {
  # Calculate k and kd
  k<-.cank(weather$zenith,vegp$vegx)
  z<-weather$zenith*pi/180
  kd<-k*cos(z)/weather$si
  kd[weather$si==0]<-1
  # Adjust paramaters for gap fraction and inclined surface
  pai<-with(vegp,pait/(1-clump))
  pai_a<-with(vegp,paia/(1-clump))
  z[z>pi/2]<-pi/2
  sk<-weather$si/cos(z)
  gref2<-gref/sk
  k0<-.cank(weather$zenith,0)
  Fd<-vegp$clump^(1/k0)
  Fb<-vegp$clump^(k/k0)
  # Calculate first set of two-stream base parameters
  om<-with(vegp,lref+ltra)
  a<-1-om
  del<-with(vegp,lref-ltra)
  mla<-(9.65*(3+vegp$vegx)^(-1.65))
  mla[mla>pi/2]<-pi/2
  J<-cos(mla)^2
  gma<-0.5*(om+J*del)
  s<-0.5*(om+J*del/kd)*kd
  sstr<-om*kd-s
  # Calculate second set of two-stream base parameters
  h<-sqrt(a^2+2*a*gma)
  sig<-kd^2+gma^2-(a+gma)^2
  S1<-exp(-h*pai)
  S2<-exp(-kd*pai)
  u1<-a+gma*(1-1/gref)
  u2<-a+gma*(1-gref)
  D1<-(a+gma+h)*(u1-h)*1/S1-(a+gma-h)*(u1+h)*S1
  D2<-(u2+h)*1/S1-(u2-h)*S1
  # Calculate Diffuse radiation parameters
  p1<-(gma/(D1*S1))*(u1-h)
  p2<-(-gma*S1/D1)*(u1+h)
  p3<-(1/(D2*S1))*(u2+h)
  p4<-(-S1/D2)*(u2- h)
  # Calculate Direct radiation parameters
  u1<-a+gma*(1-1/gref2)
  u2<-a+gma*(1-gref2)
  D1<-(a+gma+h)*(u1-h)*1/S1-(a+gma-h)*(u1+h)*S1
  D2<-(u2+h)*1/S1-(u2-h)*S1
  p5<- -s*(a+gma-kd)-gma*sstr
  v1<-s-(p5*(a+gma+kd))/sig
  v2<-s-gma-(p5/sig)*(u1+kd)
  p6<-(1/D1)*((v1/S1)*(u1-h)-(a+gma-h)*S2*v2)
  p7<-(-1/D1)*((v1*S1)*(u1+h)-(a+gma+h)*S2*v2)
  p8<-sstr*(a+gma+kd)-gma*s
  v3<-(sstr+gma*gref2-(p8/sig)*(u2-kd))*S2
  p9<-(-1/D2)*((p8/(sig*S1))*(u2+h)+v3)
  p10<-(1/D2)*(((p8*S1)/sig)*(u2-h)+v3)
  twostreamp<-list(p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6,p7=p7,p8=p8,p9=p9,p10=p10,
                   Fd=Fd,Fb=Fb,sig=sig,h=h,kd=kd)
  return(twostreamp)
}
#' Internal function for calculating effective shade factor
.ShadeF<-function(microin,twostreamp) {
  pai<-with(microin@vegp,pait/(1-clump))
  pai_a<-with(microin@vegp,paia/(1-clump))
  dirr<-with(microin@weather,swrad-difrad)
  gRbdown<-with(twostreamp,((1-Fb)*exp(-kd*pai)+Fb)*dirr)
  gRdbm<-with(twostreamp,((1-Fd)*((p8/sig)*exp(-kd*pai)+p9*exp(-h*pai)+p10*exp(h*pai))+Fd)*dirr)
  gRdbm[gRdbm<0]<-0
  gRddm<-with(twostreamp,((1-Fd)*(p3*exp(-h*pai)+p4*exp(h*pai))+Fd)*microin@weather$difrad)
  gRddown<-gRdbm+gRddm
  gdown<-gRbdown+gRddown
  Shade<-1-gdown/microin@weather$swrad
  Shade[is.na(Shade)]<-mean(Shade,na.rm=T)
  return(Shade)
}
#' Internal function for creating NIcheMapR input
.NicheMapRin<-function(microin,soilparams,twostreamp) {
  loc<-microin@loc
  tmehr<-as.POSIXlt(microin@weather$obs_time,tz="UTC")
  nyears<-length(unique(tmehr$year))
  fail<-nyears*24*365
  ystart<-tmehr$year[1]+1900
  yfinish<-tmehr$year[length(tmehr)]+1900
  yearlist<-seq(ystart,(ystart+(nyears-1)),1)
  doy <- unique(as.numeric(strftime(tmehr,format="%j")))
  ndays <- unique(paste(as.numeric(strftime(tmehr,format="%j")),
                        as.numeric(strftime(tmehr,format="%y"))))
  ndays <- length(ndays)
  doynum<-ndays
  ida<-ndays
  microdaily<-1
  daystart<-1
  # Set root properties
  L<-c(0,0,8.2,8,7.8,7.4,7.1,6.4,5.8,4.8,4,1.8,0.9,0.6,0.8,0.4,0.4,0,0)*10000
  R1<-0.001
  RW<-2.5e+10
  RL<-2e+06
  PC<- -1500
  SP<-10
  IM<-1e-06
  # Set snow properties
  snowtemp<-1.5
  snowdens<-0.375
  densfun<-c(0.5979,0.2178,0.001,0.0038)
  snowmelt<-1
  undercatch<-1
  rainmelt<-0.0125
  grasshade<-ifelse(microin@vegp$veghgt<0.5,1,0) #VT
  ### LAI etc
  PAI<-rep(microin@vegp$pai,ndays)
  MAXSHADES<-rep(100,ndays)
  # Calculate shade factor
  Shade<-.ShadeF(microin,twostreamp)
  MINSHADES<-matrix(Shade,ncol=24)
  MINSHADES<-apply(MINSHADES,1,mean)*100
  intercept<-mean(MINSHADES)/100*0.3 # snow interception
  x<-t(as.matrix(as.numeric(c(loc[1],loc[2]))))
  ALREF<-abs(trunc(x[1]))
  HEMIS<-ifelse(x[2]<0,2,1)
  ALAT<-abs(trunc(x[2]))
  AMINUT<-(abs(x[2])-ALAT)*60
  ALONG<- abs(trunc(x[1]))
  ALMINT<-(abs(x[1])-ALONG)*60
  azmuth<-groundp$aspect
  lat<-as.numeric(loc[2])
  long<-as.numeric(loc[1])
  Density<-2.56
  Thcond<-2.5
  SpecHeat<-870
  if (is.na(microin@groundp$soiltype) == FALSE) {
    sel<-which(soilparams$Soil.type == groundp$soiltype)
    if (length(sel) == 0) stop("Erroneous soil type specified")
    PE<-rep(soilparams$psi_e[sel],19)
    BB<-rep(soilparams$b[sel],19)
    BD<-rep(soilparams$rho[sel],19)
    KS<-rep(CampNormTbl9_1$Ks[sel],19)
  }
  BulkDensity<-BD[seq(1,19,2)]
  #########################
  ZENhr<-microin@weather$zenith
  ZENhr[ZENhr>90]<-90
  TAIRhr<-microin@weather$temp
  SOLRhr<-microin@weather$swrad
  sb<-5.67*10^-8
  IRDhr<-microin@weather$skyem*sb*(microin@weather$temp+273.15)^4
  RHhr<-microin@weather$relhum
  RHhr[RHhr>100]<-100
  RHhr[RHhr<0]<-0
  e0<-.satvp(TAIRhr)
  ea<-e0*(RHhr/100)
  eo<-1.24*(10*ea/(TAIRhr+273.15))^(1/7)
  CLDhr<-((microin@weather$skyem-eo)/(1-eo))*100
  CLDhr[CLDhr<0]<-0
  CLDhr[CLDhr>100]<-100
  WNhr<-microin@weather$windspeed
  WNhr[is.na(WNhr)]<-0.1
  PRESShr<-microin@weather$pres*1000
  RAINFALL<-microin@prec
  RAINFALL[RAINFALL<0.1]<-0
  ZENhr2<-ZENhr
  ZENhr2[ZENhr2!=90]<-0
  dmaxmin<-function(x,fun) {
    dx <- t(matrix(x, nrow = 24))
    apply(dx, 1, fun)
  }
  TMAXX<-dmaxmin(TAIRhr,max)
  TMINN<-dmaxmin(TAIRhr,min)
  CCMAXX<-dmaxmin(CLDhr,max)
  CCMINN<-dmaxmin(CLDhr,min)
  RHMAXX<-dmaxmin(RHhr,max)
  RHMINN<-dmaxmin(RHhr,min)
  WNMAXX<-dmaxmin(WNhr,max)
  WNMINN<-dmaxmin(WNhr,min)
  PRESS<-dmaxmin(PRESShr,min)
  ###
  slope <- 0
  azmuth <- 0
  relhum <- 1
  optdep.summer<-as.data.frame(rungads(loc[2],loc[1],relhum, 0))
  optdep.winter<-as.data.frame(rungads(loc[2],loc[1],relhum, 1))
  optdep<-cbind(optdep.winter[,1],rowMeans(cbind(optdep.summer[,2],optdep.winter[,2])))
  optdep<-as.data.frame(optdep)
  colnames(optdep)<-c("LAMBDA","OPTDEPTH")
  a<-lm(OPTDEPTH~poly(LAMBDA,6,raw=TRUE),data=optdep)
  LAMBDA<-c(290,295,300,305,310,315,320,330,340,350,360,370,380,390,400,420,440,460,480,
            500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,
            880,900,920,940,960,980,1000,1020,1080,1100,1120,1140,1160,1180,1200,1220,
            1240,1260,1280,1300,1320,1380,1400,1420,1440,1460,1480,1500,1540,1580,1600,
            1620,1640,1660,1700,1720,1780,1800,1860,1900,1950,2000,2020,2050,2100,2120,
            2150,2200,2260,2300,2320,2350,2380,2400,2420,2450,2490,2500,2600,2700,2800,
            2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000)
  TAI<-predict(a,data.frame(LAMBDA))
  RAINFALL<-RAINFALL
  ALLMINTEMPS<-TMINN
  ALLMAXTEMPS<-TMAXX
  ALLTEMPS<-cbind(ALLMAXTEMPS,ALLMINTEMPS)
  WNMAXX<-WNMAXX
  WNMINN<-WNMINN
  WNhr<-WNhr
  REFLS<-rep(groundp$gref,ndays)
  PCTWET <-rep(0,ndays)
  soilwet<-RAINFALL
  soilwet[soilwet<=1.5]<-0
  soilwet[soilwet>0]<-90
  if (ndays < 1) PCTWET<-pmax(soilwet,PCTWET)
  Intrvls<-rep(0,ndays)
  Intrvls[1]<-1
  Numtyps<-10
  Nodes<-matrix(data=0,nrow=10,ncol=ndays)
  Nodes[1:10,]<-c(1:10)
  ALREF<-abs(trunc(x[1]))
  HEMIS<-ifelse(x[2]<0,2,1)
  ALAT<-abs(trunc(x[2]))
  AMINUT<-(abs(x[2])-ALAT)*60
  ALONG<-abs(trunc(x[1]))
  ALMINT<-(abs(x[1])-ALONG)*60
  avetemp<-(sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
  soilinit<-rep(avetemp,20)
  tannul<-mean(unlist(ALLTEMPS))
  deepsoil<-rep(mean(TAIRhr),ndays)
  SLES<-matrix(nrow=ndays,data=0)
  SLES<-SLES+0.95
  moists2<-matrix(nrow=10,ncol=ndays,data=0)
  moists2[1:10,]<-c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3)
  moists<-moists2
  soilprops<-matrix(data=0,nrow=10,ncol=5)
  soilprops[,1]<-BulkDensity
  soilprops[,2]<-min(0.26,1-BulkDensity/Density)
  soilprops[,3]<-Thcond
  soilprops[,4]<-SpecHeat
  soilprops[,5]<-Density
  soilprops[1:2,3] <- 0.2
  soilprops[1:2,4] <- 1920
  hourly<-1
  if (length(microin@prec) == length(TAIRhr)) {
    rainhourly<-1
    RAINhr<-prec
    RAINhr[RAINhr<0.1]<-0
    raintest<-RAINhr
  } else if (length(microin@prec) == length(TAIRhr)/24) {
    rainhourly<-0
    RAINhr<-rep(microin@prec/24,each=24)
  } else stop("Precipitation must be daily or hourly")
  snowmodel <- 0
  RUF<-0.1*microin@vegp$veghgt
  D0<-0.65*microin@vegp$veghgt
  Refhyt<-ifelse(microin@vegp$veghgt>2,microin@vegp$veghgt,2)
  microinput<-c(ndays,RUF,1.5,0.05,Refhyt,Numtyps,0,0,0,0,1,ida,
                HEMIS,ALAT,AMINUT,ALONG,ALMINT,ALREF,groundp$slope,azmuth,0,1,
                microdaily,tannul,0.0167238,1,snowtemp,snowdens,snowmelt,undercatch,
                rainmult=1,runshade=0,1,maxpool=1000,0,snowmodel,rainmelt,
                0,densfun,hourly,rainhourly,0,0,RW,PC,RL,SP,R1,IM,
                500,0,0,fail,0,intercept,grasshade,0,0,D0)
  ###
  doy1<-matrix(data=0,nrow=ndays,ncol=1)
  SLES1<-matrix(data=0,nrow=ndays,ncol=1)
  MAXSHADES1<-matrix(data=0,nrow=ndays,ncol=1)
  MINSHADES1<-matrix(data=0,nrow=ndays,ncol=1)
  TMAXX1<-matrix(data=0,nrow=ndays,ncol=1)
  TMINN1<-matrix(data=0,nrow=ndays,ncol=1)
  CCMAXX1<-matrix(data=0,nrow=ndays,ncol=1)
  CCMINN1<-matrix(data=0,nrow=ndays,ncol=1)
  RHMAXX1<-matrix(data=0,nrow=ndays,ncol=1)
  RHMINN1<-matrix(data=0,nrow=ndays,ncol=1)
  WNMAXX1<-matrix(data=0,nrow=ndays,ncol=1)
  WNMINN1<-matrix(data=0,nrow=ndays,ncol=1)
  REFLS1<-matrix(data= 0,nrow=ndays,ncol=1)
  PCTWET1<-matrix(data=0,nrow=ndays, ncol=1)
  RAINFALL1<-matrix(data=0,nrow=ndays,ncol=1)
  tannul1<-matrix(data=0,nrow=ndays,ncol=1)
  moists1<-matrix(data=0,nrow=10,ncol=ndays)
  SLES1[1:ndays]<-SLES
  MAXSHADES1[1:ndays]<-MAXSHADES
  MINSHADES1[1:ndays]<-MINSHADES
  TMAXX1[1:ndays]<-TMAXX
  TMINN1[1:ndays]<-TMINN
  CCMAXX1[1:ndays]<-CCMAXX
  CCMINN1[1:ndays]<-CCMINN
  RHMAXX1[1:ndays]<-RHMAXX
  RHMINN1[1:ndays]<-RHMINN
  WNMAXX1[1:ndays]<-WNMAXX
  WNMINN1[1:ndays]<-WNMINN
  REFLS1[1:ndays]<-REFLS
  PCTWET1[1:ndays]<-PCTWET
  raind<-matrix(RAINhr,ncol=24,byrow=TRUE)
  raind<-apply(raind,1,sum)
  RAINFALL1[1:ndays]<-raind
  tannul1[1:ndays]<-tannul
  moists1[1:10, 1:ndays] <- moists
  tides<-matrix(data=0,nrow=24*ndays,ncol=3)
  TIMAXS<-c(1,1,0,0)
  TIMINS<-c(0,0,1,1)
  LAI<-0.8*PAI
  DEP = c(0,2.5,5,10,15,20,30,50,100,200)
  hori = rep(0,36)
  DD=rep(2.65,19)
  micro<-list(tides=tides,microinput=microinput,doy=doy,SLES=SLES1,DEP=DEP,Nodes=Nodes,
              MAXSHADES=MAXSHADES,MINSHADES=MINSHADES,TIMAXS=TIMAXS,TIMINS=TIMINS,TMAXX=TMAXX1,
              TMINN=TMINN1,RHMAXX=RHMAXX1,RHMINN=RHMINN1,CCMAXX=CCMAXX1,CCMINN=CCMINN1,
              WNMAXX=WNMAXX1,WNMINN=WNMINN1,TAIRhr=TAIRhr,RHhr=RHhr,WNhr=WNhr,CLDhr=CLDhr,
              SOLRhr=SOLRhr,RAINhr=RAINhr,ZENhr=ZENhr,IRDhr=IRDhr,REFLS=REFLS1,PCTWET=PCTWET1,
              soilinit=soilinit,hori=hori,TAI=TAI,soilprops=soilprops,moists=moists1,
              RAINFALL=RAINFALL1,tannulrun=deepsoil,PE=PE,KS=KS,BB=BB,BD=BD,DD=DD,L=L,LAI=LAI)
  return(micro)
}
#' Internal function for calculating stomatal conductance
.stomcond <- function(Rsw, gsmax, q50 = 100) {
  rpar <- Rsw * 4.6
  gs <- (gsmax * rpar) / (rpar + q50)
  gs
}
#' Internal function for calculating sensible heat below canopy
.calcH<-function(microin,miconeout,ws,i) {
  z<-microin@vegp$z
  dfo<-microin@weather[i,]
  # Calculate radiation
  rads<-canopyrad(microin,i)
  # Calculate absorbed shortwave (two sides)
  k<-.cank(dfo$zenith,microin@vegp$vegx)
  swabs<-(1-microin@vegp$lref)*with(rads,Rddown+Rdup+k*Rbdown)
  # Calculate absorbed longwave (two sides)
  tra<-exp(-microin@vegp$paia) # Transmission from top
  trb<-with(microin@vegp,exp(-(pait-paia))) # Transmission from btm
  lwfol<-microin@vegp$em*5.67*10^-8*(mean(miconeout$tleaf)+273.15)^4
  lwsoil<-0.97*5.67*10^-8*(dfo$soiltemp+273.15)^4
  lwdn<-tra*dfo$lwdown+(1-tra)*lwfol
  lwup<-trb*lwsoil+(1-trb)*lwfol
  lwabs<-microin@vegp$em*(lwdn+lwup)
  # Calculate absorbed and emitted radiation (two-sides)
  Rabs<-swabs+lwabs
  Rem<-2*microin@vegp$em*5.67*10^-8*(miconeout$tleaf+273.15)^4
  # Calculate stomatal conductance
  gsd<-.stomcond(rads$Rddown+rads$Rbdown,microin@vegp$gsmax)
  gsu<-.stomcond(rads$Rdup,microin@vegp$gsmax)
  gs<-0.5*gsd+0.5*gsu
  # Adjust wind to canopy top using standard profile
  vz<-microin@vegp$veghgt+2
  uh2<-dfo$windspeed*log(67.8*vz-5.42)/4.87 # 2m above canopy
  uh2[uh2<0.5]<-0.5
  d<-0.65*microin@vegp$veghgt
  zm<-0.1*microin@vegp$veghgt
  uf<-(0.4*uh2)/log((vz-d)/zm)
  uh<-(uf/0.4)*log((microin@vegp$veghgt-d)/zm)
  uz<-uh*ws
  # Calculate boundary layer conductance
  gHa<-0.19*sqrt(uz/microin@vegp$leafd)
  gHf<-0.0525*with(miconeout,abs(tleaf-tair))/microin@vegp$leafd^0.25
  gHa<-pmax(gHa,gHf)
  # Calculate vapour conductance
  gv<-1/(1/gHa+1/gs)
  # Calculate Latent heat exchange
  es<-.satvp(miconeout$tleaf)
  ea<-with(miconeout,.satvp(tair)*(rh/100))
  L<-2*44355.38*gv/dfo$pres*(es-ea)
  L[L<0]<-0
  # Calculate one-sided H
  H<-0.5*(Rabs-Rem-L)
  # Calculate leaf temperature
  tleaf<-miconeout$tair+H/(29.3*gHa)
  # Recalaculate L
  es<-.satvp(tleaf)
  L<-44355.38*gv/dfo$pres*(es-ea)
  return(list(H=H,L=L,tleaf=tleaf,gv=gv,gHa=gHa,z=z,dfo=dfo,uf=uf,
              Rbdown=rads$Rbdown,Rddown=rads$Rddown,Rdup=rads$Rdup,Lwdown=lwdn,Lwup=lwup))
}
#' Internal function for calculating temperature and humidity and top of canopy
.CalcTh<-function(microin,HT,i) {
  H<-sum(HT$H*microin@vegp$paii)
  vz<-microin@vegp$veghgt+2
  d<-0.65*microin@vegp$veghgt
  zh<-0.02*microin@vegp$veghgt
  lr<-log((microin@vegp$veghgt-d)/zh)/log((vz-d)/zh)
  TH<-microin@weather$temp[i]+(H/(503.96*HT$uf))*log((vz-d)/zh)
  eA<-with(microin@weather,.satvp(temp[i])*relhum[i]/100)
  dp<-.dewpoint(eA,microin@weather$temp[i])
  TH<-ifelse(TH<dp,dp,TH)
  Th<-TH-(TH-microin@weather$temp[i])*lr
  # Compute vapour pressure at top of canopy
  esa<-.satvp(TH)
  eh<-esa-(esa-eA)*lr
  return(list(Th=Th,eh=eh))
}
#' Internal function for running L-Nf-T Langrangian model
.Langrangian<-function(HT,microin,miconeout,h) {
  z<-microin@vegp$z
  n<-length(z)
  # Compute thermal diffusivity within canopy
  a<-sqrt((0.1*HT$uf*microin@vegp$veghgt)/(1.25^2))/0.3
  TL<-0.3
  ow<-a*(0.75+0.5*cos(pi*(1-z/microin@vegp$veghgt)))
  Kc<-ow^2*TL
  # Compute source concentration
  leafdens<-with(microin@vegp,paii/veghgt)
  ST<-leafdens*HT$H
  SL<-leafdens*HT$L
  # Compute heat flux from ground
  ea<-with(miconeout,.satvp(tair)*rh/100)
  dfo<-HT$dfo
  dCT<-(dfo$soiltemp-miconeout$tair)*29.3*43
  eG<-.satvp(dfo$soiltemp) # soil effective vapour pressure
  dCL<-(eG-ea)*44526*43/dfo$pres
  GT<-0
  GL<-0
  for (i in 1:n) {
    GT[i]<-(sum(Kc[1:i])/i)*dCT[i]/(microin@vegp$veghgt*i/n)
    GL[i]<-(sum(Kc[1:i])/i)*dCL[i]/(microin@vegp$veghgt*i/n)
  }
  H<-HT$H+GT
  L<-HT$L+GL
  # Compute reference height near-field and far-field
  Zeta<-abs((microin@vegp$veghgt-z)/(ow*TL))
  kn<- -0.39894*log(1-exp(-Zeta))-0.15623*exp(-Zeta)
  CnzrT<-with(microin@vegp,sum((ST/ow)*(kn*((veghgt-z)/(ow*TL))+kn*((veghgt+z)/(ow*TL))),na.rm=T))
  CnzrL<-with(microin@vegp,sum((SL/ow)*(kn*((veghgt-z)/(ow*TL))+kn*((veghgt+z)/(ow*TL))),na.rm=T))
  CT<-0
  CL<-0
  for (i in 1:length(z)) {
    # Near field
    Zeta<-abs((z[i]-z)/(ow*TL))
    kn<- -0.39894*log(1-exp(-Zeta))-0.15623*exp(-Zeta)
    CnT<-sum((ST/ow)*(kn*((z[i]-z)/(ow*TL))+kn*((z[i]+z)/(ow*TL))),na.rm=T)
    CnL<-sum((SL/ow)*(kn*((z[i]-z)/(ow*TL))+kn*((z[i]+z)/(ow*TL))),na.rm=T)
    # Far field
    CfsT<-sum(H[i:n]/Kc[i:n])*(microin@vegp$veghgt/n)
    CfsL<-sum(L[i:n]/Kc[i:n])*(microin@vegp$veghgt/n)
    CfT<-(h$Th*29.3*43)-CnzrT+CfsT
    CfL<-(h$eh*44526*(43/dfo$pres))-CnzrL+CfsL
    # Both
    CT[i]<-CfT+CnT
    CL[i]<-CfL+CnL
  }
  ta<-CT/(29.3*43)
  # Limits
  mn<-min(dfo$soiltemp,h$Th)-5
  mx<-max(dfo$soiltemp,h$Th)+5
  ta[ta<mn]<-mn
  ta[ta>mx]<-mx
  ea<-(CL*dfo$pres)/(44526*43)
  ea[ea<0.01]<-0.01
  es<-.satvp(ta)
  rh<-(ea/es)*100
  rh[rh>100]<-100
  rh[rh<10]<-10
  # Calculate tleaf
  tleaf<-ta+HT$H/(29.3*HT$gHa)
  # Limits
  tleaf[tleaf<mn]<-mn
  mx<-max(dfo$soiltemp,h$Th)+10
  tleaf[tleaf>mx]<-mx
  # Set at dewpoint
  dp<-.dewpoint(ea,ta)
  tleaf[tleaf<dp]<-dp[tleaf<dp]
  return(list(tair=ta,tleaf=tleaf,rh=rh))
}
#' #' Internal function for iterating L-Nf-T Langrangian model
.LangrangianR<-function(HT,microin,miconeout,h,R=5) {
  for (i in 1:R) {
    Lout<-.Langrangian(HT,microin,miconeout,h)
    miconeout<-data.frame(tair=(Lout$tair+miconeout$tair)/2,
                          tleaf=(Lout$tleaf+miconeout$tleaf)/2,
                          rh=(Lout$rh+miconeout$rh)/2)
  }
  return(miconeout)
}
