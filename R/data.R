## Dataset documentation.
##
## Deliberately written as `@name` + `@docType data` + NULL rather than as a
## bare "itaload" string literal. The string-literal idiom makes roxygen2 look
## the object up in the loaded package, via `::` semantics in some versions, so
## documenting the package fails with
##
##     'itaload' is not an exported object from 'namespace:doobr'
##
## whenever data/itaload.rda is absent or has not been registered as lazy data
## -- a fresh clone that did not fetch the binary, a build before
## `LazyData: true` took effect, or a roxygen2/pkgload version that resolves the
## topic through the namespace exports. The form below needs no object at all,
## so `devtools::document()` works on the sources alone.

#' Italian daily electricity load
#'
#' Daily electricity load for Italy, aggregated from hourly observations, with
#' the daily mean temperature and the calendar variables used in the worked
#' example. Covers 1 January 2023 to 30 September 2025.
#'
#' The series has the features that exercise the method: strong deterministic
#' seasonality at the weekly and annual frequencies, genuine conditional
#' heteroskedasticity, and a horizon long enough for the sampling error of a
#' path ensemble to matter.
#'
#' @name itaload
#' @docType data
#' @keywords datasets
#' @usage data(itaload)
#' @format A data frame with 1004 rows and 7 columns:
#' \describe{
#'   \item{date}{Date.}
#'   \item{load}{Daily total load, MWh. Mean 426 GWh, coefficient of variation 0.15.}
#'   \item{temperature}{Daily mean temperature, degrees Celsius. Note that this
#'     is a \emph{stochastic} regressor: see the caveat in the electricity
#'     vignette before using it as if it were known in advance.}
#'   \item{holiday}{Italian national holiday not falling at a weekend.}
#'   \item{bridge}{Weekday sandwiched between two non-working days.}
#'   \item{dow}{Day of the week, a factor with levels \code{Mon} to \code{Sun}.}
#'   \item{month}{Month, a factor with levels \code{Jan} to \code{Dec}.}
#' }
#' @source Aggregated from hourly Italian system load and temperature records.
#'   The construction script is in \code{data-raw/make_itaload.R}.
#' @examples
#' data(itaload)
#' plot(itaload$date, itaload$load / 1e3, type = "l",
#'      xlab = "", ylab = "GWh", main = "Italian daily electricity load")
NULL
