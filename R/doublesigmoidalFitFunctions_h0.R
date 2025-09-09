#' Double‐sigmoidal NLS fit (with baseline h0)
#'
#' Fits `doubleSigmoidalFitFormula_h0()` to a `data.frame(time, intensity)`
#' via `minpack.lm::nlsLM()`, appending h0 and AIC/BIC.
#'
#' @param dataInput      data.frame or list with `$timeIntensityData`.
#' @param tryCounter     Integer: which random‐start iteration this is.
#' @param startList      Named list of starting values (including h0).
#' @param lowerBounds    Named vector of lower bounds.
#' @param upperBounds    Named vector of upper bounds.
#' @param min_Factor     Passed to `nlsLM()`.
#' @param n_iterations   Max iterations for `nlsLM()`.
#' @return A one‐row data.frame of estimates, fit‐statistics, and
#'   `h0_Estimate`.
#' @export
#' @importFrom minpack.lm nlsLM
sigFit2 <- function (
    dataInput,
    tryCounter,
    startList = list(
      finalAsymptoteIntensityRatio = 0.5,
      maximum                      = 1,
      slope1Param                  = 1,
      midPoint1Param               = 0.33,
      slope2Param                  = 1,
      midPointDistanceParam        = 0.29,
      h0                           = 0.1     # ← added
    ),

    lowerBounds = c(
      finalAsymptoteIntensityRatio = 0,
      maximum                      = 0.3,
      slope1Param                  = 0.01,
      midPoint1Param               = -0.52,
      slope2Param                  = 0.01,
      midPointDistanceParam        = 0.04,
      h0                           = 0    # ← added
    ),

    upperBounds = c(
      finalAsymptoteIntensityRatio = 1.5,
      maximum                      = 1.5,
      slope1Param                  = 180,
      midPoint1Param               = 1.15,
      slope2Param                  = 180,
      midPointDistanceParam        = 0.63,
      h0                           = 0.3    # ← added
    ),
    min_Factor   = 1/2^20,
    n_iterations = 1000
)
{
  isalist     <- (is.list(dataInput) & !is.data.frame(dataInput))
  if (isalist) {
    dataFrameInput <- dataInput$timeIntensityData
  }
  isadataframe <- is.data.frame(dataInput)
  if (isadataframe) {
    dataFrameInput <- dataInput
  }
  if (tryCounter == 1) {
    counterDependentStartList <- startList
  }
  if (tryCounter != 1) {
    randomVector <- stats::runif(length(startList), 0, 1)
    names(randomVector) <- names(startList)    # now includes "h0"
    counterDependentStartVector <- randomVector * (upperBounds - lowerBounds) + lowerBounds
    counterDependentStartList <- as.list(counterDependentStartVector)
  }

  theFitResult <- try(
    minpack.lm::nlsLM(
      intensity ~ doubleSigmoidalFitFormula_h0(
        time,
        finalAsymptoteIntensityRatio,
        maximum,
        slope1Param,
        midPoint1Param,
        slope2Param,
        midPointDistanceParam,
        h0
      ),
      data    = dataFrameInput,
      start   = counterDependentStartList,
      control = list(maxiter = n_iterations, minFactor = min_Factor),
      lower   = lowerBounds,
      upper   = upperBounds,
      trace   = FALSE,
    ),
    silent = TRUE
  )

  if (!inherits(theFitResult, "try-error")) {
    parameterMatrix <- summary(theFitResult)$parameters
    colnames(parameterMatrix) <- c("Estimate", "Std_Error", "t_value", "Pr_t")
    parameterVector <- c(t(parameterMatrix))
    # append 4 h0 names at the end
    names(parameterVector) <- c(
      "finalAsymptoteIntensityRatio_N_Estimate", "finalAsymptoteIntensityRatio_Std_Error",
      "finalAsymptoteIntensityRatio_t_value",    "finalAsymptoteIntensityRatio_Pr_t",
      "maximum_N_Estimate", "maximum_Std_Error", "maximum_t_value", "maximum_Pr_t",
      "slope1Param_N_Estimate","slope1Param_Std_Error","slope1Param_t_value","slope1Param_Pr_t",
      "midPoint1Param_N_Estimate","midPoint1Param_Std_Error","midPoint1Param_t_value","midPoint1Param_Pr_t",
      "slope2Param_N_Estimate","slope2Param_Std_Error","slope2Param_t_value","slope2Param_Pr_t",
      "midPointDistanceParam_N_Estimate","midPointDistanceParam_Std_Error","midPointDistanceParam_t_value","midPointDistanceParam_Pr_t",
      "h0_N_Estimate","h0_Std_Error","h0_t_value","h0_Pr_t"   # ← appended
    )
    parameterVector <- c(
      parameterVector,
      residual_Sum_of_Squares = sum((as.vector(stats::resid(theFitResult)))^2)[1],
      log_likelihood         = as.vector(stats::logLik(theFitResult))[1],
      AIC_value              = as.vector(stats::AIC(theFitResult))[1],
      BIC_value              = as.vector(stats::BIC(theFitResult))[1]
    )
    parameterList <- as.list(parameterVector)
    parameterList$isThisaFit <- TRUE
    parameterList$startVector <- counterDependentStartList
    if (isalist) parameterList$dataScalingParameters <- as.list(dataInput$dataScalingParameters)
    parameterList$model <- "doublesigmoidal"
    parameterList$additionalParameters <- FALSE
    parameterDf <- as.data.frame(parameterList)
    parameterDf <- renormalizeParameters_h0(parameterDf, isalist)
  } else {
    parameterVector <- rep(NA, 24 + 4)
    names(parameterVector) <- c(
      "finalAsymptoteIntensityRatio_N_Estimate","finalAsymptoteIntensityRatio_Std_Error",
      "finalAsymptoteIntensityRatio_t_value","finalAsymptoteIntensityRatio_Pr_t",
      "maximum_N_Estimate","maximum_Std_Error","maximum_t_value","maximum_Pr_t",
      "slope1Param_N_Estimate","slope1Param_Std_Error","slope1Param_t_value","slope1Param_Pr_t",
      "midPoint1Param_N_Estimate","midPoint1Param_Std_Error","midPoint1Param_t_value","midPoint1Param_Pr_t",
      "slope2Param_N_Estimate","slope2Param_Std_Error","slope2Param_t_value","slope2Param_Pr_t",
      "midPointDistanceParam_N_Estimate","midPointDistanceParam_Std_Error","midPointDistanceParam_t_value","midPointDistanceParam_Pr_t",
      "h0_N_Estimate","h0_Std_Error","h0_t_value","h0_Pr_t"   # ← appended
    )
    parameterVector <- c(
      parameterVector,
      residual_Sum_of_Squares = Inf,
      log_likelihood         = NA,
      AIC_value              = NA,
      BIC_value              = NA
    )
    parameterList <- as.list(parameterVector)
    parameterList$isThisaFit <- FALSE
    parameterList$startVector     <- counterDependentStartList
    if (isalist) parameterList$dataScalingParameters <- as.list(dataInput$dataScalingParameters)
    parameterList$model <- "doublesigmoidal"
    parameterDf <- as.data.frame(parameterList)
    parameterDf <- renormalizeParameters_h0(parameterDf, isalist)
  }

  return(parameterDf)  # unchanged
}


#' Double‐sigmoidal fit formula with baseline h0
#'
#' Evaluates the full 5‐parameter “rise+fall” model at `x`.
#'
#' @param x                         Numeric vector.
#' @param finalAsymptoteIntensityRatio Ratio of second plateau to first.
#' @param maximum                   Maximum plateau height.
#' @param slope1Param               First slope (>0).
#' @param midPoint1Param            First midpoint (>0).
#' @param slope2Param               Second slope (>0).
#' @param midPointDistanceParam     Distance to second midpoint (>0).
#' @param h0                        Baseline when x is negative infinity.
#' @return Numeric vector of same length as `x`.
#' @export
doubleSigmoidalFitFormula_h0 <- function (x, finalAsymptoteIntensityRatio, maximum, slope1Param,
                                          midPoint1Param, slope2Param, midPointDistanceParam, h0)
{
  if (slope1Param < 0) {
    stop("slope1Param should be a positive number")
  }
  if (slope2Param < 0) {
    stop("slope2Param should be a positive number. It is the absolute value of the second slopeParam")
  }
  if (midPointDistanceParam < 0) {
    stop("midPointDistanceParam should be a positive number. It is the distance between two steppest points of exponential phase and lysis")
  }
  #if (finalAsymptoteIntensityRatio < 0 | finalAsymptoteIntensityRatio >
  #  1) {
  #  stop("finalAsymptoteIntensityRatio should be a number between 0 and 1")
  #}
  #if (maximum < 0) {
  #  stop("maximum should be a positive number")
  #}
  optimizeIntervalMin <- midPoint1Param - 2 * midPointDistanceParam
  optimizeIntervalMax <- midPoint1Param + 3 * midPointDistanceParam
  xmax <- stats::optimize(
    f1,
    c(optimizeIntervalMin, optimizeIntervalMax),
    tol = 1e-4,
    B1 = slope1Param,
    M1 = midPoint1Param,
    B2 = slope2Param,
    L  = midPointDistanceParam,
    maximum = TRUE
  )
  argumentt <- xmax$maximum
  constt <- f0(argumentt, slope1Param, midPoint1Param, slope2Param,
               midPointDistanceParam)
  y <- f2_h0(x, finalAsymptoteIntensityRatio, maximum, slope1Param,
             midPoint1Param, slope2Param, midPointDistanceParam,
             constt, argumentt, h0)
  return(y)
}



#' Double‐sigmoidal base function
#'
#' A helper “base” function used internally by the double‐sigmoidal fit formula.
#'
#' @param x Numeric vector of time points.
#' @param B1 Positive slope parameter for the first rise.
#' @param M1 Midpoint of the first rise.
#' @param B2 Positive slope parameter for the second fall.
#' @param L  Distance between the two midpoints.
#' @return Numeric vector of the same length as `x`.
#' @keywords internal
f0 <- function (x, B1, M1, B2, L) #unchanged
{
  (1/((1 + exp(-B1 * (x - M1))) * (1 + exp(B2 * (x - (M1 + L))))))
}


# parameterCalculation_h0 Double Sigmoidal Helper Functions
#' @title Log‐transformed Double‐Sigmoidal Core Function
#' @description Computes the log of the core double‐sigmoidal form used for
#'   optimizing the peak time.  This function is internal to the fitAndCategorize parameter
#'   calculation.
#' @param x Numeric vector of time points.
#' @param B1 Numeric; slope of the rising phase.
#' @param M1 Numeric; midpoint of the rising phase.
#' @param B2 Numeric; absolute slope of the declining phase.
#' @param L  Numeric; distance between the two steepest points.
#' @return Numeric vector of \code{log(1/((1+e^{-B1(x-M1)})*(1+e^{B2(x-(M1+L))})))}.
#' @keywords internal
f1 <- function (x, B1, M1, B2, L) # untouched
{
  log((1/((1 + exp(-B1 * (x - M1))) * (1 + exp(B2 * (x - (M1 +
                                                            L)))))))
}

#' Heaviside Piecewise double‐sigmoidal
#'
#' Internal helper that stitches together two scaled copies of `f0()` at a
#' switch point.
#'
#' @param x        Numeric vector of time points.
#' @param A2       Final asymptote ratio.
#' @param Ka       Maximum plateau height.
#' @param B1       First slope (see `f0()`).
#' @param M1       First midpoint.
#' @param B2       Second slope.
#' @param L        Midpoint distance.
#' @param const    Normalizing constant (value of `f0()` at its peak).
#' @param argument Switch point on the x‐axis.
#' @param h0       Baseline intensity.
#' @return Numeric vector of same length as `x`.
#' @keywords internal
f2_h0 <- function (x, A2, Ka, B1, M1, B2, L, const, argument, h0)
{
  fBasics::Heaviside(x - argument) * (f0(x, B1, M1, B2, L) * ((Ka - A2 * Ka)/(const)) + A2 * Ka) +
    (1 - fBasics::Heaviside(x - argument)) * (f0(x, B1, M1, B2, L) * ((Ka - h0)/(const)) + h0)
}

#' Half-Maximum Midpoint Function for Double-Sigmoidal Curve
#'
#' Given the parameters of the double-sigmoidal backbone \(f_0(x)\) and the peak
#' value `const`, returns \(f_0(x) - \tfrac{1}{2}\,const\).  This is used with
#' `uniroot()` to find the time points where the curve crosses half its maximum.
#'
#' @param x Numeric vector of time points at which to evaluate.
#' @param B1 Numeric; the first slope parameter of the double-sigmoidal backbone.
#' @param M1 Numeric; the first midpoint parameter (rise center).
#' @param B2 Numeric; the second slope parameter of the double-sigmoidal backbone.
#' @param L Numeric; the distance between the two midpoints.
#' @param const Numeric; the peak (maximum) value of \(f_0(x)\).
#' @return Numeric vector of the same length as `x`, giving \(\,f_0(x) - \tfrac{1}{2}\,const\).
#' @noRd
f0mid <- function (x, B1, M1, B2, L,const) {
  -const / 2 + 1 / ( ( 1 + exp(-B1 * (x-M1)) ) * ( 1 + exp(B2*(x-(M1+L))) ) )
}

# Argmax
#' @title Argmax for Double‐Sigmoidal Peak Time
#' @description Finds the time (within the scaled range) at which the fitted
#'   double‐sigmoidal curve reaches its maximum.
#' @param parameterVector A named list or data frame containing:
#'   \describe{
#'     \item{\code{slope1Param_N_Estimate}}{Rising slope (raw).}
#'     \item{\code{slope2Param_N_Estimate}}{Declining slope (raw).}
#'     \item{\code{midPoint1Param_N_Estimate}}{Rising midpoint (raw).}
#'     \item{\code{midPointDistanceParam_N_Estimate}}{Distance between midpoints (raw).}
#'     \item{\code{dataScalingParameters.timeRange}}{Total time range scaling factor.}
#'   }
#' @return Numeric scalar: the time (in original units) of the curve’s peak.
#' @keywords internal
f_argmax_doublesigmoidal_h0 <- function (parameterVector)
{
  slope1Param <- parameterVector$slope1Param_N_Estimate
  slope2Param <- parameterVector$slope2Param_N_Estimate
  midPointDistanceParam <- parameterVector$midPointDistanceParam_N_Estimate
  finalAsymptoteIntensityRatio <- parameterVector$finalAsymptoteIntensityRatio_N_Estimate
  maximum <- parameterVector$maximum_N_Estimate
  midPoint1Param <- parameterVector$midPoint1Param_N_Estimate
  timeRange <- parameterVector$dataScalingParameters.timeRange
  if (slope1Param < 0) {
    stop("slope1Param should be a positive number")
  }
  if (slope2Param < 0) {
    stop("slope2Param should be a positive number. It is the absolute value of the slope2Param")
  }
  if (midPointDistanceParam < 0) {
    stop("midPointDistanceParam should be a positive number. It is the distance between two steppest points of exponential phase and lysis")
  }
  #if (finalAsymptoteIntensityRatio < 0 | finalAsymptoteIntensityRatio >
  #  1) {
  #  stop("finalAsymptoteIntensityRatio should be a number between 0 and 1")
  #}
  #if (maximum < 0) {
  #  stop("maximum should be a positive number")
  #}
  optimizeIntervalMin <- midPoint1Param - 2 * midPointDistanceParam
  optimizeIntervalMax <- midPoint1Param + 3 * midPointDistanceParam
  xmax <- stats::optimize(f1, c(optimizeIntervalMin, optimizeIntervalMax),
                          tol = 1e-04, B1 = slope1Param, M1 = midPoint1Param,
                          B2 = slope2Param, L = midPointDistanceParam, maximum = TRUE)
  argumentt <- xmax$maximum * timeRange
  return(argumentt)
}

# Mid1
#' Find First Half-Maximum Point of Double-Sigmoidal Fit
#'
#' Given a parameter data frame from a double-sigmoidal fit, finds the
#' first time point (x) at which the model reaches half of its maximum
#' response by solving for the root of the mid-point equation.
#'
#' @param parameterDf A data.frame or list containing at least:
#'   - `dataScalingParameters.timeRange`
#'   - `slope1Param_Estimate`
#'   - `midPoint1Param_Estimate`
#'   - `slope2Param_Estimate`
#'   - `midPointDistanceParam_Estimate`
#' @return A numeric scalar giving the time of the first half-maximum point.
#' @importFrom stats optimize uniroot
#' @export
f_mid1_doublesigmoidal_h0 <- function (parameterDf)
{
  max_x <- parameterDf$dataScalingParameters.timeRange
  xmax <- stats::optimize(f1, interval = c(-1.125 * max_x,
                                           max_x * 3), tol = 1e-04, B1 = parameterDf$slope1Param_Estimate,
                          M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
                          L = parameterDf$midPointDistanceParam_Estimate, maximum = TRUE)
  argumentt <- xmax$maximum
  constt <- f0(argumentt, B1 = parameterDf$slope1Param_Estimate,
               M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
               L = parameterDf$midPointDistanceParam_Estimate)
  mid1x <- try(expr = stats::uniroot(f0mid, interval = c(-1.125 *
                                                           max_x, argumentt), tol = 1e-04, B1 = parameterDf$slope1Param_Estimate,
                                     M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
                                     L = parameterDf$midPointDistanceParam_Estimate, const = constt),
               silent = TRUE)
  if (!inherits(mid1x, "try-error"))  {
    rangeList <- seq(-1, 30)
    rangeListCounter <- -1
    while (class(mid1x) == "try-error" & rangeListCounter <=
           30) {
      rangeListSelection <- rangeList[rangeListCounter +
                                        2]
      mid1x <- try(expr = stats::uniroot(f0mid, interval = c((-1) *
                                                               (2^rangeListSelection) * max_x, argumentt),
                                         tol = 1e-04, B1 = parameterDf$slope1Param_Estimate,
                                         M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
                                         L = parameterDf$midPointDistanceParam_Estimate,
                                         const = constt), silent = TRUE)
      rangeListCounter <- rangeListCounter + 1
    }
  }
  if (inherits(mid1x, "try-error")) {
    stop("f_mid1_doublesigmoidal_h0(): unable to find first half-max root with uniroot()")
  }
  return(mid1x$root)
}

# Mid2
#' Find Second Half-Maximum Point of Double-Sigmoidal Fit
#'
#' For a double-sigmoidal fit, locates the time at which the response
#' declines back to half-maximum, by root-finding on the mid-point function
#' over the interval from the peak to beyond the original time range.
#'
#' @param parameterDf A data.frame or list containing the same fields as
#'   \code{f_mid1_doublesigmoidal_h0()}, including \code{dataScalingParameters.timeRange}.
#' @return A numeric scalar giving the time of the second half-maximum point.
#' @importFrom stats optimize uniroot
#' @export
f_mid2_doublesigmoidal_h0 <- function (parameterDf)
{
  max_x <- parameterDf$dataScalingParameters.timeRange
  xmax <- stats::optimize(f1, interval = c(-1.125 * max_x,
                                           max_x * 3), tol = 1e-04, B1 = parameterDf$slope1Param_Estimate,
                          M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
                          L = parameterDf$midPointDistanceParam_Estimate, maximum = TRUE)
  argumentt <- xmax$maximum
  constt <- f0(argumentt, B1 = parameterDf$slope1Param_Estimate,
               M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
               L = parameterDf$midPointDistanceParam_Estimate)
  mid2x <- try(stats::uniroot(f0mid, interval = c(argumentt,
                                                  max_x * (1.25)), tol = 1e-04, B1 = parameterDf$slope1Param_Estimate,
                              M1 = parameterDf$midPoint1Param_Estimate, B2 = parameterDf$slope2Param_Estimate,
                              L = parameterDf$midPointDistanceParam_Estimate, const = constt),
               silent = TRUE)
  #if (!inherits(mid2x, "try-error")) { (MT changing statement 8/26)
  if (inherits(mid2x, "try-error")) {
    rangeList <- seq(-1, 30)
    rangeListCounter <- -1
    # while (class(mid2x) == "try-error" & rangeListCounter <=
    #        30) { (MT changing statement 8/26)
    while (inherits(mid2x, "try-error") && rangeListCounter <= 30) {
      rangeListSelection <- rangeList[rangeListCounter +
                                        2]
      mid2x <- try(stats::uniroot(f0mid, interval = c(argumentt,
                                                      max_x * (2^rangeListSelection)), tol = 1e-04,
                                  B1 = parameterDf$slope1Param_Estimate, M1 = parameterDf$midPoint1Param_Estimate,
                                  B2 = parameterDf$slope2Param_Estimate, L = parameterDf$midPointDistanceParam_Estimate,
                                  const = constt), silent = TRUE)
      rangeListCounter <- rangeListCounter + 1
    }
  }

  if (inherits(mid2x, "try-error")) {
    stop("f_mid2_doublesigmoidal_h0(): unable to find second half-max root with uniroot()")
  }

  return(mid2x$root)
}


# Slope
#' @title Numerical Slope for Double‐Sigmoidal at a Given Time
#' @description Approximates the first derivative of the fitted
#'   double‐sigmoidal curve at time \code{x} using a 5‐point stencil.
#' @param x Numeric scalar; time point at which to compute the slope.
#' @param parameterDf A named list or data frame of fitted parameters, including:
#'   \describe{
#'     \item{\code{finalAsymptoteIntensityRatio_Estimate}}{Final ratio.}
#'     \item{\code{maximum_Estimate}}{Maximum intensity.}
#'     \item{\code{slope1Param_Estimate}}{Rising slope parameter.}
#'     \item{\code{midPoint1Param_Estimate}}{Rising midpoint.}
#'     \item{\code{slope2Param_Estimate}}{Declining slope parameter.}
#'     \item{\code{midPointDistanceParam_Estimate}}{Distance between midpoints.}
#'     \item{\code{h0_Estimate}}{Baseline intensity.}
#'   }
#' @param timeStep Numeric; small increment used for finite differences. Default 1e-05.
#' @return Numeric scalar: approximate slope (dY/dX) at \code{x}.
#' @keywords internal
f_slope_doublesigmoidal_h0 <- function (x, parameterDf, timeStep = 1e-05)
{
  fxp2h <- doubleSigmoidalFitFormula_h0(x = x + 2 * timeStep,
                                        finalAsymptoteIntensityRatio = parameterDf$finalAsymptoteIntensityRatio_Estimate,
                                        maximum = parameterDf$maximum_Estimate, slope1Param = parameterDf$slope1Param_Estimate,
                                        midPoint1Param = parameterDf$midPoint1Param_Estimate,
                                        slope2Param = parameterDf$slope2Param_Estimate, midPointDistanceParam = parameterDf$midPointDistanceParam_Estimate, h0 = parameterDf$h0_Estimate)
  fxph <- doubleSigmoidalFitFormula_h0(x = x + timeStep, finalAsymptoteIntensityRatio = parameterDf$finalAsymptoteIntensityRatio_Estimate,
                                       maximum = parameterDf$maximum_Estimate, slope1Param = parameterDf$slope1Param_Estimate,
                                       midPoint1Param = parameterDf$midPoint1Param_Estimate,
                                       slope2Param = parameterDf$slope2Param_Estimate, midPointDistanceParam = parameterDf$midPointDistanceParam_Estimate, h0 = parameterDf$h0_Estimate)
  fxmh <- doubleSigmoidalFitFormula_h0(x = x - timeStep, finalAsymptoteIntensityRatio = parameterDf$finalAsymptoteIntensityRatio_Estimate,
                                       maximum = parameterDf$maximum_Estimate, slope1Param = parameterDf$slope1Param_Estimate,
                                       midPoint1Param = parameterDf$midPoint1Param_Estimate,
                                       slope2Param = parameterDf$slope2Param_Estimate, midPointDistanceParam = parameterDf$midPointDistanceParam_Estimate, h0 = parameterDf$h0_Estimate)
  fxm2h <- doubleSigmoidalFitFormula_h0(x = x - 2 * timeStep,
                                        finalAsymptoteIntensityRatio = parameterDf$finalAsymptoteIntensityRatio_Estimate,
                                        maximum = parameterDf$maximum_Estimate, slope1Param = parameterDf$slope1Param_Estimate,
                                        midPoint1Param = parameterDf$midPoint1Param_Estimate,
                                        slope2Param = parameterDf$slope2Param_Estimate, midPointDistanceParam = parameterDf$midPointDistanceParam_Estimate, h0 = parameterDf$h0_Estimate)
  der <- (-1 * fxp2h + 8 * fxph - 8 * fxmh + 1 * fxm2h)/(12 *
                                                           timeStep)
  return(der)
}

#' Convert normalized min.pack estimated parameters back to raw scale (with h0)
#'
#' Takes a one‐row parameter data.frame and “unscales” its
#' time‐normalized and intensity‐normalized columns.
#'
#' @param parameterDF Data.frame with columns `*_N_Estimate` and
#'   `dataScalingParameters.*`.
#' @param isalist     Logical; if `TRUE` input came from a list,
#'   otherwise already raw.
#' @return A one‐row data.frame with `*_Estimate` columns on the raw
#'   scale (including `h0_Estimate`).
#' @export
#' @examples
#' df <- data.frame(
#'   model="doublesigmoidal",
#'   maximum_N_Estimate=1, finalAsymptoteIntensityRatio_N_Estimate = .5,
#'   slope1Param_N_Estimate=1, midPoint1Param_N_Estimate = .3,
#'   slope2Param_N_Estimate=1, midPointDistanceParam_N_Estimate = .2,
#'   h0_N_Estimate = .1,
#'   dataScalingParameters.timeRange = 100,
#'   dataScalingParameters.intensityRange = 10,
#'   dataScalingParameters.intensityMin = 0
#' )
#' renormalizeParameters_h0(df, isalist=TRUE)
renormalizeParameters_h0 <- function(parameterDF, isalist) {
  model         <- parameterDF$model
  dataInputName <- parameterDF$dataInputName
  if (isalist) {
    # raw plateau heights
    h1_raw <- parameterDF$maximum_N_Estimate  * parameterDF$dataScalingParameters.intensityRange + parameterDF$dataScalingParameters.intensityMin
    h2_raw <- (fAIR_n <- parameterDF$finalAsymptoteIntensityRatio_N_Estimate * parameterDF$maximum_N_Estimate) * parameterDF$dataScalingParameters.intensityRange + parameterDF$dataScalingParameters.intensityMin
    parameterDF$finalAsymptoteIntensityRatio_Estimate <-
      h2_raw / h1_raw



    # parameterDF$finalAsymptoteIntensityRatio_Estimate <-
    #  parameterDF$finalAsymptoteIntensityRatio_N_Estimate #old version


    parameterDF$maximum_Estimate <-
      parameterDF$maximum_N_Estimate *
      parameterDF$dataScalingParameters.intensityRange +
      parameterDF$dataScalingParameters.intensityMin
    parameterDF$slope1Param_Estimate <-
      parameterDF$slope1Param_N_Estimate /
      parameterDF$dataScalingParameters.timeRange
    parameterDF$midPoint1Param_Estimate <-
      parameterDF$midPoint1Param_N_Estimate *
      parameterDF$dataScalingParameters.timeRange
    parameterDF$slope2Param_Estimate <-
      parameterDF$slope2Param_N_Estimate /
      parameterDF$dataScalingParameters.timeRange
    parameterDF$midPointDistanceParam_Estimate <-
      parameterDF$midPointDistanceParam_N_Estimate *
      parameterDF$dataScalingParameters.timeRange
    parameterDF$h0_Estimate <-
      parameterDF$h0_N_Estimate *
      parameterDF$dataScalingParameters.intensityRange +
      parameterDF$dataScalingParameters.intensityMin
  }
  if (!isalist) {
    # if already raw, just copy over
    parameterDF$finalAsymptoteIntensityRatio_Estimate <-
      parameterDF$finalAsymptoteIntensityRatio_N_Estimate
    parameterDF$maximum_Estimate <-
      parameterDF$maximum_N_Estimate
    parameterDF$slope1Param_Estimate <-
      parameterDF$slope1Param_N_Estimate
    parameterDF$midPoint1Param_Estimate <-
      parameterDF$midPoint1Param_N_Estimate
    parameterDF$slope2Param_Estimate <-
      parameterDF$slope2Param_N_Estimate
    parameterDF$midPointDistanceParam_Estimate <-
      parameterDF$midPointDistanceParam_N_Estimate
    parameterDF$h0_Estimate <-
      parameterDF$h0_N_Estimate
  }
  parameterDF$model         <- model
  parameterDF$dataInputName <- dataInputName
  return(parameterDF)
}
