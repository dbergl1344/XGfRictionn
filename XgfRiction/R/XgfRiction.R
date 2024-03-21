library(roxygen2)
library(devtools)
create_package("~/Documents/GitHub/XGfRictionn/XgfRiction")
library(REddyProc)
library(caret)
library(Metrics)
library(stats)
library(tidymodels)
library(ggplot2)


#' This class represents a processor for XGFriction data.
#' XGFrictionProcessor Class
#'
#' This class represents a processor for XGFriction data.
#'
#' @field site The name of the site associated with the processor.
#' @import REddyProc
#' @export
XGFrictionProcessor <- setRefClass("XGFrictionProcessor",
                                   fields = list(site = "character"),
                                   methods = list(
                                     initialize = function(site) {
                                       site <<- site
                                     },
                                     initialize_xgfriction_processing = function(dataframe, LatDeg, LongDeg, TimeZoneHour) {
                                       require(REddyProc)
                                       # Check if the data is in 30-minute timestamp format
                                       if (any(diff(dataframe$DateTime) != 30*60)) {
                                         mins <- 15*60
                                         print("Adding 15 minutes to DateTime")
                                         print(dataframe$DateTime + mins)
                                         dataframe$DateTime <- (dataframe$DateTime + mins)
                                       }

                                       # Initialize the processing
                                       xgfriction_proc <- REddyProc::sEddyProc$new(site, dataframe$DateTime,
                                                                                   c('NEE', 'Rg', 'Tair', 'VPD', 'Ustar'))
                                       xgfriction_proc$sSetLocationInfo(LatDeg = LatDeg, LongDeg = LongDeg, TimeZoneHour = TimeZoneHour)
                                       return(xgfriction_proc)
                                     }
                                   )
)


load_all()

