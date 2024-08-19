
#' @title treein S4 class
#' input class from microspline or tree trunk modelling
#' @slot Microclimate list.
#' @slot Radiation list.
#' @slot InputWeather data.frame.
#' @slot Treeparams list.
#' @slot FixedParams list.
#' @slot treetemps matrix.
#' @name treein
#' @export
#' @importFrom methods new
setClass(
  "treein",
  slots=c(
    Microclimate="list",
    Radiation="list",
    InputWeather="data.frame",
    Treeparams="list",
    FixedParams="list",
    treetemps="matrix"
  )
)

#' @title microin S4 class
#' eg used in modelin function
#' @slot weather data.frame.
#' @slot prec numeric.
#' @slot vegp list.
#' @slot groundp list.
#' @slot loc numeric.
#' @slot twostreamp list.
#' @slot reqhgt numeric.
#' @name microin
#' @export
#' @importFrom methods new
setClass(
  "microin",
  slots=c(
    weather="data.frame",
    prec="numeric",
    vegp="list",
    groundp="list",
    loc="numeric",
    twostreamp="list",
    reqhgt="numeric"
  )
)

#' @title microout S4 class
#' Output by runmicro function
#' @slot Microclimate list.
#' @slot Radiation list.
#' @slot InputWeather data.frame.
#' @name microout
#' @export
#' @importFrom methods new
setClass(
  "microout",
  slots=c(
    Microclimate="list",
    Radiation="list",
    InputWeather="data.frame"
  )
)
