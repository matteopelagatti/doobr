## Builds the `itaload` dataset shipped with the package.
##
## Source: hourly Italian electricity load and temperature, aggregated to daily
## totals (load) and daily means (temperature).  The calendar columns are
## precomputed so that the vignette needs no extra dependency: `holiday` flags
## Italian national holidays not falling at a weekend, `bridge` flags weekdays
## sandwiched between two non-working days.
raw <- utils::read.csv("data-raw/itaload.csv", stringsAsFactors = FALSE)
itaload <- data.frame(
  date        = as.Date(raw$date),
  load        = as.numeric(raw$load),
  temperature = as.numeric(raw$temperature),
  holiday     = as.logical(raw$holiday),
  bridge      = as.logical(raw$bridge),
  dow         = factor(raw$dow, levels = 1:7,
                       labels = c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")),
  month       = factor(raw$month, levels = 1:12, labels = month.abb)
)
stopifnot(nrow(itaload) == 1004L, !anyNA(itaload$load))
save(itaload, file = "data/itaload.rda", compress = "xz", version = 2)
cat(sprintf("itaload: %d rows, %s to %s, %.1f KB\n", nrow(itaload),
            min(itaload$date), max(itaload$date),
            file.size("data/itaload.rda") / 1024))
