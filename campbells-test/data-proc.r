# Process data for model runs
# 
#

##### Read the data #####
campdat <- read.csv("UPC-level-Exponential-3.csv", nrows = 3000, stringsAsFactors = FALSE) 

cc <- sapply(campdat, class)

cc["UPC"] <- "character"
cc["SCV_MAG_SPEND"] <- "numeric"
cc["SCV_NTL_NEWSP_SPND"] <- "numeric"


campdat <- read.csv("UPC-level-Exponential-3.csv", colClasses = cc)
campdat <- data.table(campdat)
setkey(campdat, UPC)


############# Create a sequential id mapped to UPCs ######
idmap <- campdat[,list(UPC = unique(UPC))]
setkey(idmap, UPC)
idmap[, id := seq(.N)]
campdat <- idmap[campdat]

############# Add variables for max ACV, ACV range #######
acvrange <- campdat[, list(acv.max = max(ACV.MultiOutlet), acv.range = diff(range(ACV.MultiOutlet))), keyby = UPC]
campdat <- campdat[acvrange]
numtriers <- campdat[, list(numtriers = unique(NO_TRIERS)), keyby = UPC]

############ Drop variables from previous model runs #########
campdat[,`:=`(p0 = NULL, lambda = NULL, PredTrialPct = NULL, MAPE = NULL)]

############ Filter out cases with 30 or fewer triers and with max ACV < 5 ###########
campdat2 <- campdat[NO_TRIERS > 30 & !is.na(acv.max) & acv.max > 5]

############ Split data by sequential id
campdat2lst <- split(campdat2, campdat2$id)

save(campdat2lst, file = "campbells-model-data.rdata")
