
#' @title Compute Intuitively Meaningful Curve Parameters
#' @description From the raw fit vector of either a sigmoidal or double‐sigmoidal
#'   model, calculates intuitively interpretable quantities (e.g. peak times,
#'   slopes at midpoints, plateau heights) and appends them to the parameter list.
#' @param parameterVector A named list or data frame produced by
#'   \code{\link{sigFit}} or \code{\link{sigFit2}}, containing at minimum:
#'   \describe{
#'     \item{\code{model}}{Either \dQuote{sigmoidal} or \dQuote{doublesigmoidal}.}
#'     \item{\code{*_Estimate}}{Raw parameter estimates (renormalized back to the original scale).}
#'     \item{\code{dataScalingParameters.timeRange},
#'           \code{dataScalingParameters.intensityRange},
#'           \code{dataScalingParameters.intensityMin}}{
#'       Scaling parameters used to normalize/de‐normalize the data.}
#'   }
#' @param stepSize Numeric; finite‐difference increment for slope calculation. Default is \code{1e-05}.
#' @return The input \code{parameterVector}, augmented with new fields:
#'   for sigmoidal: \code{midPoint_x}, \code{midPoint_y}, \code{slope}, etc.;
#'   for double‐sigmoidal: \code{maximum_x}, \code{midPoint1_x}, \code{slope1},
#'   \code{finalAsymptoteIntensity}, etc.
#' @export
parameterCalculation_h0 <- function (parameterVector, stepSize = 1e-05)
{
  if (parameterVector$model == "sigmoidal") { # MT making slight adjustments 8/26
    parameterVector$maximum_x <- NA
    parameterVector$maximum_y <- parameterVector$maximum_Estimate

    parameterVector$midPoint_x <- parameterVector$midPoint_Estimate

    parameterVector$h0_y <- parameterVector$h0_Estimate # changed

    # changing midpoint_y calculation (MT changing 8/26):
    amp <- parameterVector$maximum_y - parameterVector$h0_y
    parameterVector$midPoint_y <- parameterVector$h0_y + amp/2

    # parameterVector$midPoint_y <- parameterVector$maximum_y/2


    # subtracting h0 from maximum here (MT changing 8/26)
    parameterVector$slope <- parameterVector$slopeParam_Estimate *
      (parameterVector$maximum_y - parameterVector$h0_y) * (1/4)

    # changing computation of increment time (MT changing 8/26)
    if (parameterVector$slopeParam_Estimate <= 0) {
      warning("slopeParam_Estimate should be > 0 for increasing sigmoids.")
      inc_time <- Inf
    } else {
      inc_time <- 4 / parameterVector$slopeParam_Estimate
    }
    parameterVector$incrementTime <- inc_time
    # parameterVector$incrementTime <- parameterVector$maximum_y / parameterVector$slope (old version)

    parameterVector$startPoint_x <- parameterVector$midPoint_x -
      (parameterVector$incrementTime / 2)

    parameterVector$startPoint_y <- parameterVector$h0_y # changed, potentially set this to 0 to see what happens

    parameterVector$reachMaximum_x <- parameterVector$midPoint_x +
      (parameterVector$incrementTime / 2)
    parameterVector$reachMaximum_y <- parameterVector$maximum_y

    parameterVector$additionalParameters <- TRUE
  }

  if (parameterVector$model == "doublesigmoidal") {
    parameterVector$maximum_x <- f_argmax_doublesigmoidal_h0(parameterVector)
    parameterVector$maximum_y <- parameterVector$maximum_Estimate

    parameterVector$midPoint1_x <- f_mid1_doublesigmoidal_h0(parameterVector)
    parameterVector$midPoint1_y <- doubleSigmoidalFitFormula_h0(x = parameterVector$midPoint1_x,
                                                                finalAsymptoteIntensityRatio = parameterVector$finalAsymptoteIntensityRatio_Estimate,
                                                                maximum = parameterVector$maximum_y,
                                                                slope1Param = parameterVector$slope1Param_Estimate,
                                                                midPoint1Param = parameterVector$midPoint1Param_Estimate,
                                                                slope2Param = parameterVector$slope2Param_Estimate,
                                                                midPointDistanceParam = parameterVector$midPointDistanceParam_Estimate,
                                                                h0 = parameterVector$h0_Estimate)

    parameterVector$midPoint2_x <- f_mid2_doublesigmoidal_h0(parameterVector)
    parameterVector$midPoint2_y <- doubleSigmoidalFitFormula_h0(x = parameterVector$midPoint2_x,
                                                                finalAsymptoteIntensityRatio = parameterVector$finalAsymptoteIntensityRatio_Estimate,
                                                                maximum = parameterVector$maximum_y,
                                                                slope1Param = parameterVector$slope1Param_Estimate,
                                                                midPoint1Param = parameterVector$midPoint1Param_Estimate,
                                                                slope2Param = parameterVector$slope2Param_Estimate,
                                                                midPointDistanceParam = parameterVector$midPointDistanceParam_Estimate,
                                                                h0 = parameterVector$h0_Estimate)

    parameterVector$slope1 <- f_slope_doublesigmoidal_h0(parameterVector$midPoint1_x, # changed
                                                         parameterVector, timeStep = stepSize)

    parameterVector$slope2 <- f_slope_doublesigmoidal_h0(parameterVector$midPoint2_x, # changed
                                                         parameterVector, timeStep = stepSize)

    parameterVector$finalAsymptoteIntensity <- parameterVector$finalAsymptoteIntensityRatio_Estimate * parameterVector$maximum_y

    parameterVector$incrementTime <- parameterVector$maximum_y/parameterVector$slope1

    parameterVector$startPoint_x <- parameterVector$midPoint1_x -
      (parameterVector$incrementTime/2)

    parameterVector$startPoint_y <- parameterVector$h0_Estimate # changed, not 0 anymore

    parameterVector$reachMaximum_x <- parameterVector$midPoint1_x +
      (parameterVector$incrementTime/2)

    if (parameterVector$reachMaximum_x > parameterVector$maximum_x) {
      parameterVector$reachMaximum_x <- parameterVector$maximum_x
      parameterVector$warning.reachMaximum_cor = TRUE
    }

    parameterVector$reachMaximum_y <- parameterVector$maximum_y

    parameterVector$decrementTime <- (parameterVector$maximum_y -
                                        parameterVector$finalAsymptoteIntensity)/(-parameterVector$slope2)

    parameterVector$startDeclinePoint_x <- parameterVector$midPoint2_x -
      (parameterVector$decrementTime/2)

    #MT change 8/29 -- sig function has maximum_x set to NA automatically, need to guard against that

    if (isTRUE(parameterVector$reachMaximum_x > parameterVector$maximum_x)) { # new version
      parameterVector$reachMaximum_x <- parameterVector$maximum_x
      parameterVector$warning.reachMaximum_cor <- TRUE
    }

    #old version (8/29)
    # if (parameterVector$startDeclinePoint_x < parameterVector$maximum_x) {
    #   parameterVector$startDeclinePoint_x <- parameterVector$maximum_x
    #   parameterVector$warning.startDeclinePoint_cor = TRUE
    # }

    parameterVector$startDeclinePoint_y <- parameterVector$maximum_y
    parameterVector$endDeclinePoint_x <- parameterVector$midPoint2_x +
      (parameterVector$decrementTime/2)
    parameterVector$endDeclinePoint_y <- parameterVector$finalAsymptoteIntensity
    parameterVector$additionalParameters = TRUE
  }
  return(parameterVector)
}
