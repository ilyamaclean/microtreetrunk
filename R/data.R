#' A list of soil properties by soil type used in NicheMapR
#'
#' @format a data.frame with the following columns:
#' \describe{
#'  \item{Texture}{soil type (see [soilparams()])}
#'  \item{Silt}{proportion of silt (0-1)}
#'  \item{Clay}{proportion of clay (0-1)}
#'  \item{airentry}{}
#'  \item{b}{}
#'  \item{Ks}{}
#'  \item{field}{}
#'  \item{wilting}{}
#' }
"CampNormTbl9_1"
#'
#'#'
#' A list of soil and ground parameters required for running the model
#'
#' @format a list with the following elements:
#' \describe{
#'  \item{slope}{slope of gorund surface (decimal degrees)}
#'  \item{aspect}{aspect of ground surface (decimal degrees from N)}
#'  \item{gref}{ground refklectance (0-1)}
#'  \item{soil type}{soil type (see [soilparams()])}
#' }
"groundp"
#'
#' Microclimate model outputs at 5 m above
#' ground at Caerthillian Cove, Lizard, Cornwall (49.96807N, 5.215668W) for the
#' period 16th-30th June 2017, assuming a 10 m tall canopy and with the following
#' slots.
#' @format an object of class microout with the following slots:
#' \describe{
#'   \item{Microclimate}{A list of air and foliage temperatures and relative humidities}
#'   \item{Radiation}{A list of upward and downward radiative fluxes}
#'   \item{InputWeather}{A data.frame of weather variables provided as inputs}
#' }
"microout"
#'
#' A vector of daily rainfall at Caerthillian Cove, Lizard, Cornwall
#' (49.96807N, 5.215668W) for the period 16th-30th June 2017
#' @format A vector of daily rainfall (mm/day)
"prec"
#'
#' A table of soil parameters for different soil types
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Soil.type}{description of soil type}
#'   \item{Smax}{Volumetric water content at saturation (m^3 / m^3)}
#'   \item{Smin}{Residual water content (m^3 / m^3)}
#'   \item{Ksat}{Saturated hydraulic conductivity (kg / m^3 / day)}
#'   \item{Vq}{Volumetric quartz content of soil}
#'   \item{Vm}{Volumetric mineral content of soil}
#'   \item{Vo}{Volumetric organic content of soil}
#'   \item{Mc}{Mass fraction of clay}
#'   \item{rho}{Soil bulk density (Mg / m^3)}
#'   \item{b}{Shape parameter for NicheMapR soil moisture model (dimensionless, > 1)}
#' }
#' @source: \url{https://onlinelibrary.wiley.com/doi/full/10.1002/ird.1751}
"soilparams"
#' A list of vegetation parameters required for running the microclimate model
#'
#' @format a list with the following elements:
#' \describe{
#'  \item{pait}{total plant index value for the canopy}
#'  \item{vegx}{ratio of vertical to horizontal projections of leaf foliage}
#'  \item{clump}{fraction of radiation passing through larger gaps in the canopy (0-1)}
#'  \item{lref}{leaf reflectance (0-1)}
#'  \item{ltra}{leaf transmittance (0-1)}
#'  \item{veghgt}{canopy height (m)}
#'  \item{em}{leaf emissivity (0-1)}
#'  \item{gsmax}{maximum stomatal conductance of leaves (mol / m^2 / s)}
#'  \item{leafd}{leaf diameter (m)}
#'  \item{paii}{individual leaf area index values for each canopy layer, here set to NA as skew and spread are used to estimate this when running the model}
#'  \item{skew}{degree of skew towards top of canopy in canopy foliage (0-10, see [PAIgeometry())}
#'  \item{spread}{degree of spread in canopy foliage (0-100, see [PAIgeometry())}
#' }
"vegp"
#'
#' A list of tree parameters required for running the tree trunk model
#'
#' @format a list with the following elements:
#' \describe{
#'  \item{nlayers}{total plant index value for the canopy}
#'  \item{nsegs}{ratio of vertical to horizontal projections of leaf foliage}
#'  \item{tref}{tree trunk reflectance (0-1)}
#'  \item{tem}{tree trunk emissivity (0-1}
#'  \item{radius}{tree trunk radius (m)}
#'  \item{K}{a vector of thermal conductivities of each layer (W/m/K)}
#'  \item{cph}{a vector of volumeric specific heat capacities of each layer (MJ/m^3)}
#' }
"treeparams"
#'
#' A data frame of hourly weather at Caerthillian Cove, Lizard, Cornwall
#' (49.96807N, 5.215668W) for the period 16th-30th June 2017
#'
#' @format a data frame with the following elements:
#' \describe{
#'  \item{obs_time}{POSIXlt object of dates and times}
#'  \item{temp}{temperature (degrees C)}
#'  \item{relhum}{relative humidity (percentage)}
#'  \item{pres}{atmospheric press (kPa)}
#'  \item{swrad}{Total incoming shortwave radiation (W / m^2)}
#'  \item{difrad}{Diffuse radiation (W / m^2)}
#'  \item{skyem}{Sky emissivity (0-1)}
#'  \item{windspeed}{Wind speed (m/s)}
#'  \item{winddir}{Wind direction (decimal degrees)}
#' }
"weather"
#'
#'
#'

