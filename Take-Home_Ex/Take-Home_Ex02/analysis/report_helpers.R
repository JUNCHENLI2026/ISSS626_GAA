# Locale-independent English month labels for both reporting documents.
format_month_en <- function(x) {
  dates <- as.Date(x)
  labels <- paste(month.name[as.integer(format(dates, "%m"))],
                  format(dates, "%Y"))
  labels[is.na(dates)] <- NA_character_
  labels
}
