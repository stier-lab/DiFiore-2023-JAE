##########################################
## Libraries
##########################################

library(tidybayes)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(here)
library(lme4)
library(lmerTest)
library(purrr)
library(rethinking)
library(rstanarm)

# size of experimental tanks (sides and bottom -- i.e. areal NOT volumetric)
tsize <- 137/2 * 76 /10000

set.seed(112618)

#--------------------------------------------
## Read in all data directly from EDI
#--------------------------------------------

# Package ID: knb-lter-sbc.156.1 Cataloging System:https://pasta.edirepository.org.
# Data set title: SBC LTER: Reef: The size-dependent functional response of lobster foraging on purple urchin.
# Data set creator:    - Santa Barbara Coastal LTER 
# Data set creator:  Bartholomew DiFiore -  
# Data set creator:  Adrian Stier -  
# Contact:    - Information Manager, Santa Barbara Coastal LTER   - sbclter@msi.ucsb.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/156/1/809131efb9b03a7bf5306dc3d343d5ca" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "id",     
                 "size",     
                 "treatment",     
                 "initial",     
                 "killed",     
                 "udiam",     
                 "mr",     
                 "mc"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$id)!="factor") dt1$id<- as.factor(dt1$id)
if (class(dt1$size)=="factor") dt1$size <-as.numeric(levels(dt1$size))[as.integer(dt1$size) ]               
if (class(dt1$size)=="character") dt1$size <-as.numeric(dt1$size)
if (class(dt1$treatment)!="factor") dt1$treatment<- as.factor(dt1$treatment)
if (class(dt1$initial)=="factor") dt1$initial <-as.numeric(levels(dt1$initial))[as.integer(dt1$initial) ]               
if (class(dt1$initial)=="character") dt1$initial <-as.numeric(dt1$initial)
if (class(dt1$killed)=="factor") dt1$killed <-as.numeric(levels(dt1$killed))[as.integer(dt1$killed) ]               
if (class(dt1$killed)=="character") dt1$killed <-as.numeric(dt1$killed)
if (class(dt1$udiam)=="factor") dt1$udiam <-as.numeric(levels(dt1$udiam))[as.integer(dt1$udiam) ]               
if (class(dt1$udiam)=="character") dt1$udiam <-as.numeric(dt1$udiam)
if (class(dt1$mr)=="factor") dt1$mr <-as.numeric(levels(dt1$mr))[as.integer(dt1$mr) ]               
if (class(dt1$mr)=="character") dt1$mr <-as.numeric(dt1$mr)
if (class(dt1$mc)=="factor") dt1$mc <-as.numeric(levels(dt1$mc))[as.integer(dt1$mc) ]               
if (class(dt1$mc)=="character") dt1$mc <-as.numeric(dt1$mc)

write.csv(dt1, "data/cleaned/loburc_cleaned.csv", row.names = F)



# Package ID: knb-lter-sbc.26.21 Cataloging System:https://pasta.edirepository.org.
# Data set title: SBC LTER: Reef: Long-term experiment: Kelp removal: Urchin size frequency distribution.
# Data set creator:    - Santa Barbara Coastal LTER 
# Data set creator:  Daniel C Reed -  
# Data set creator:  Robert J Miller -  
# Contact:    - Information Manager, Santa Barbara Coastal LTER   - sbclter@msi.ucsb.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/26/21/ee9c73aa15fb459f8f2fbee2d1f33044" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "YEAR",     
                 "MONTH",     
                 "DATE",     
                 "SITE",     
                 "TRANSECT",     
                 "TREATMENT",     
                 "SP_CODE",     
                 "SIZE",     
                 "COUNT",     
                 "AREA",     
                 "NOTES",     
                 "SCIENTIFIC_NAME",     
                 "COMMON_NAME",     
                 "TAXON_KINGDOM",     
                 "TAXON_PHYLUM",     
                 "TAXON_CLASS",     
                 "TAXON_ORDER",     
                 "TAXON_FAMILY",     
                 "TAXON_GENUS",     
                 "GROUP",     
                 "SURVEY",     
                 "MOBILITY",     
                 "GROWTH_MORPH"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$MONTH)!="factor") dt1$MONTH<- as.factor(dt1$MONTH)                                   
# attempting to convert dt1$DATE dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1DATE<-as.Date(dt1$DATE,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1DATE) == length(tmp1DATE[!is.na(tmp1DATE)])){dt1$DATE <- tmp1DATE } else {print("Date conversion failed for dt1$DATE. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1DATE) 
if (class(dt1$SITE)!="factor") dt1$SITE<- as.factor(dt1$SITE)
if (class(dt1$TRANSECT)!="factor") dt1$TRANSECT<- as.factor(dt1$TRANSECT)
if (class(dt1$TREATMENT)!="factor") dt1$TREATMENT<- as.factor(dt1$TREATMENT)
if (class(dt1$SP_CODE)!="factor") dt1$SP_CODE<- as.factor(dt1$SP_CODE)
if (class(dt1$SIZE)=="factor") dt1$SIZE <-as.numeric(levels(dt1$SIZE))[as.integer(dt1$SIZE) ]               
if (class(dt1$SIZE)=="character") dt1$SIZE <-as.numeric(dt1$SIZE)
if (class(dt1$COUNT)=="factor") dt1$COUNT <-as.numeric(levels(dt1$COUNT))[as.integer(dt1$COUNT) ]               
if (class(dt1$COUNT)=="character") dt1$COUNT <-as.numeric(dt1$COUNT)
if (class(dt1$AREA)=="factor") dt1$AREA <-as.numeric(levels(dt1$AREA))[as.integer(dt1$AREA) ]               
if (class(dt1$AREA)=="character") dt1$AREA <-as.numeric(dt1$AREA)
if (class(dt1$NOTES)!="factor") dt1$NOTES<- as.factor(dt1$NOTES)
if (class(dt1$SCIENTIFIC_NAME)!="factor") dt1$SCIENTIFIC_NAME<- as.factor(dt1$SCIENTIFIC_NAME)
if (class(dt1$COMMON_NAME)!="factor") dt1$COMMON_NAME<- as.factor(dt1$COMMON_NAME)
if (class(dt1$TAXON_KINGDOM)!="factor") dt1$TAXON_KINGDOM<- as.factor(dt1$TAXON_KINGDOM)
if (class(dt1$TAXON_PHYLUM)!="factor") dt1$TAXON_PHYLUM<- as.factor(dt1$TAXON_PHYLUM)
if (class(dt1$TAXON_CLASS)!="factor") dt1$TAXON_CLASS<- as.factor(dt1$TAXON_CLASS)
if (class(dt1$TAXON_ORDER)!="factor") dt1$TAXON_ORDER<- as.factor(dt1$TAXON_ORDER)
if (class(dt1$TAXON_FAMILY)!="factor") dt1$TAXON_FAMILY<- as.factor(dt1$TAXON_FAMILY)
if (class(dt1$TAXON_GENUS)!="factor") dt1$TAXON_GENUS<- as.factor(dt1$TAXON_GENUS)
if (class(dt1$GROUP)!="factor") dt1$GROUP<- as.factor(dt1$GROUP)
if (class(dt1$SURVEY)!="factor") dt1$SURVEY<- as.factor(dt1$SURVEY)
if (class(dt1$MOBILITY)!="factor") dt1$MOBILITY<- as.factor(dt1$MOBILITY)
if (class(dt1$GROWTH_MORPH)!="factor") dt1$GROWTH_MORPH<- as.factor(dt1$GROWTH_MORPH)

# Convert Missing Values to NA for non-dates

dt1$SIZE <- ifelse((trimws(as.character(dt1$SIZE))==trimws("-99999")),NA,dt1$SIZE)               
suppressWarnings(dt1$SIZE <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$SIZE))==as.character(as.numeric("-99999"))),NA,dt1$SIZE))
dt1$COUNT <- ifelse((trimws(as.character(dt1$COUNT))==trimws("-99999")),NA,dt1$COUNT)               
suppressWarnings(dt1$COUNT <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$COUNT))==as.character(as.numeric("-99999"))),NA,dt1$COUNT))
dt1$AREA <- ifelse((trimws(as.character(dt1$AREA))==trimws("-99999")),NA,dt1$AREA)               
suppressWarnings(dt1$AREA <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$AREA))==as.character(as.numeric("-99999"))),NA,dt1$AREA))
dt1$SCIENTIFIC_NAME <- as.factor(ifelse((trimws(as.character(dt1$SCIENTIFIC_NAME))==trimws("-99999")),NA,as.character(dt1$SCIENTIFIC_NAME)))

write.csv(dt1, "data/LTER/LTE_Urchin_All_Years_20210209.csv", row.names = F)





# Package ID: knb-lter-sbc.77.6 Cataloging System:https://pasta.edirepository.org.
# Data set title: SBC LTER: Reef: Abundance, size and fishing effort for California Spiny Lobster (Panulirus interruptus), ongoing since 2012.
# Data set creator:    - Santa Barbara Coastal LTER 
# Data set creator:  Daniel C Reed -  
# Data set creator:  Robert J Miller -  
# Contact:    - Information Manager, Santa Barbara Coastal LTER   - sbclter@msi.ucsb.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/77/6/f32823fba432f58f66c06b589b7efac6" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "YEAR",     
                 "MONTH",     
                 "DATE",     
                 "SITE",     
                 "TRANSECT",     
                 "REPLICATE",     
                 "SIZE_MM",     
                 "COUNT",     
                 "NUM_AO",     
                 "AREA"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$MONTH)!="factor") dt1$MONTH<- as.factor(dt1$MONTH)                                   
# attempting to convert dt1$DATE dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1DATE<-as.Date(dt1$DATE,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1DATE) == length(tmp1DATE[!is.na(tmp1DATE)])){dt1$DATE <- tmp1DATE } else {print("Date conversion failed for dt1$DATE. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1DATE) 
if (class(dt1$SITE)!="factor") dt1$SITE<- as.factor(dt1$SITE)
if (class(dt1$TRANSECT)!="factor") dt1$TRANSECT<- as.factor(dt1$TRANSECT)
if (class(dt1$REPLICATE)!="factor") dt1$REPLICATE<- as.factor(dt1$REPLICATE)
if (class(dt1$SIZE_MM)=="factor") dt1$SIZE_MM <-as.numeric(levels(dt1$SIZE_MM))[as.integer(dt1$SIZE_MM) ]               
if (class(dt1$SIZE_MM)=="character") dt1$SIZE_MM <-as.numeric(dt1$SIZE_MM)
if (class(dt1$COUNT)=="factor") dt1$COUNT <-as.numeric(levels(dt1$COUNT))[as.integer(dt1$COUNT) ]               
if (class(dt1$COUNT)=="character") dt1$COUNT <-as.numeric(dt1$COUNT)
if (class(dt1$NUM_AO)=="factor") dt1$NUM_AO <-as.numeric(levels(dt1$NUM_AO))[as.integer(dt1$NUM_AO) ]               
if (class(dt1$NUM_AO)=="character") dt1$NUM_AO <-as.numeric(dt1$NUM_AO)
if (class(dt1$AREA)=="factor") dt1$AREA <-as.numeric(levels(dt1$AREA))[as.integer(dt1$AREA) ]               
if (class(dt1$AREA)=="character") dt1$AREA <-as.numeric(dt1$AREA)

# Convert Missing Values to NA for non-dates

dt1$SIZE_MM <- ifelse((trimws(as.character(dt1$SIZE_MM))==trimws("-99999")),NA,dt1$SIZE_MM)               
suppressWarnings(dt1$SIZE_MM <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$SIZE_MM))==as.character(as.numeric("-99999"))),NA,dt1$SIZE_MM))
dt1$COUNT <- ifelse((trimws(as.character(dt1$COUNT))==trimws("-99999")),NA,dt1$COUNT)               
suppressWarnings(dt1$COUNT <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$COUNT))==as.character(as.numeric("-99999"))),NA,dt1$COUNT))
dt1$NUM_AO <- ifelse((trimws(as.character(dt1$NUM_AO))==trimws("-99999")),NA,dt1$NUM_AO)               
suppressWarnings(dt1$NUM_AO <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$NUM_AO))==as.character(as.numeric("-99999"))),NA,dt1$NUM_AO))
dt1$AREA <- ifelse((trimws(as.character(dt1$AREA))==trimws("-99999")),NA,dt1$AREA)               
suppressWarnings(dt1$AREA <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$AREA))==as.character(as.numeric("-99999"))),NA,dt1$AREA))

write.csv(dt1, "data/LTER/Lobster_Abundance_All_Years_20210412.csv", row.names = F)





# Package ID: knb-lter-sbc.50.10 Cataloging System:https://pasta.edirepository.org.
# Data set title: SBC LTER: Reef: Annual time series of biomass for kelp forest species, ongoing since 2000.
# Data set creator:    - Santa Barbara Coastal LTER 
# Data set creator:  Daniel C Reed -  
# Data set creator:  Robert J Miller -  
# Contact:    - Information Manager, Santa Barbara Coastal LTER   - sbclter@msi.ucsb.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/50/10/24d18d9ebe4f6e8b94e222840096963c" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "YEAR",     
                 "MONTH",     
                 "DATE",     
                 "SITE",     
                 "TRANSECT",     
                 "SP_CODE",     
                 "PERCENT_COVER",     
                 "DENSITY",     
                 "WM_GM2",     
                 "DM_GM2",     
                 "SFDM",     
                 "AFDM",     
                 "SCIENTIFIC_NAME",     
                 "COMMON_NAME",     
                 "TAXON_KINGDOM",     
                 "TAXON_PHYLUM",     
                 "TAXON_CLASS",     
                 "TAXON_ORDER",     
                 "TAXON_FAMILY",     
                 "TAXON_GENUS",     
                 "GROUP",     
                 "MOBILITY",     
                 "GROWTH_MORPH",     
                 "COARSE_GROUPING"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

# attempting to convert dt1$DATE dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1DATE<-as.Date(dt1$DATE,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1DATE) == length(tmp1DATE[!is.na(tmp1DATE)])){dt1$DATE <- tmp1DATE } else {print("Date conversion failed for dt1$DATE. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1DATE) 
if (class(dt1$SITE)!="factor") dt1$SITE<- as.factor(dt1$SITE)
if (class(dt1$TRANSECT)!="factor") dt1$TRANSECT<- as.factor(dt1$TRANSECT)
if (class(dt1$SP_CODE)!="factor") dt1$SP_CODE<- as.factor(dt1$SP_CODE)
if (class(dt1$PERCENT_COVER)=="factor") dt1$PERCENT_COVER <-as.numeric(levels(dt1$PERCENT_COVER))[as.integer(dt1$PERCENT_COVER) ]               
if (class(dt1$PERCENT_COVER)=="character") dt1$PERCENT_COVER <-as.numeric(dt1$PERCENT_COVER)
if (class(dt1$DENSITY)=="factor") dt1$DENSITY <-as.numeric(levels(dt1$DENSITY))[as.integer(dt1$DENSITY) ]               
if (class(dt1$DENSITY)=="character") dt1$DENSITY <-as.numeric(dt1$DENSITY)
if (class(dt1$WM_GM2)=="factor") dt1$WM_GM2 <-as.numeric(levels(dt1$WM_GM2))[as.integer(dt1$WM_GM2) ]               
if (class(dt1$WM_GM2)=="character") dt1$WM_GM2 <-as.numeric(dt1$WM_GM2)
if (class(dt1$DM_GM2)=="factor") dt1$DM_GM2 <-as.numeric(levels(dt1$DM_GM2))[as.integer(dt1$DM_GM2) ]               
if (class(dt1$DM_GM2)=="character") dt1$DM_GM2 <-as.numeric(dt1$DM_GM2)
if (class(dt1$SFDM)=="factor") dt1$SFDM <-as.numeric(levels(dt1$SFDM))[as.integer(dt1$SFDM) ]               
if (class(dt1$SFDM)=="character") dt1$SFDM <-as.numeric(dt1$SFDM)
if (class(dt1$AFDM)=="factor") dt1$AFDM <-as.numeric(levels(dt1$AFDM))[as.integer(dt1$AFDM) ]               
if (class(dt1$AFDM)=="character") dt1$AFDM <-as.numeric(dt1$AFDM)
if (class(dt1$SCIENTIFIC_NAME)!="factor") dt1$SCIENTIFIC_NAME<- as.factor(dt1$SCIENTIFIC_NAME)
if (class(dt1$COMMON_NAME)!="factor") dt1$COMMON_NAME<- as.factor(dt1$COMMON_NAME)
if (class(dt1$TAXON_KINGDOM)!="factor") dt1$TAXON_KINGDOM<- as.factor(dt1$TAXON_KINGDOM)
if (class(dt1$TAXON_PHYLUM)!="factor") dt1$TAXON_PHYLUM<- as.factor(dt1$TAXON_PHYLUM)
if (class(dt1$TAXON_CLASS)!="factor") dt1$TAXON_CLASS<- as.factor(dt1$TAXON_CLASS)
if (class(dt1$TAXON_ORDER)!="factor") dt1$TAXON_ORDER<- as.factor(dt1$TAXON_ORDER)
if (class(dt1$TAXON_FAMILY)!="factor") dt1$TAXON_FAMILY<- as.factor(dt1$TAXON_FAMILY)
if (class(dt1$TAXON_GENUS)!="factor") dt1$TAXON_GENUS<- as.factor(dt1$TAXON_GENUS)
if (class(dt1$GROUP)!="factor") dt1$GROUP<- as.factor(dt1$GROUP)
if (class(dt1$MOBILITY)!="factor") dt1$MOBILITY<- as.factor(dt1$MOBILITY)
if (class(dt1$GROWTH_MORPH)!="factor") dt1$GROWTH_MORPH<- as.factor(dt1$GROWTH_MORPH)
if (class(dt1$COARSE_GROUPING)!="factor") dt1$COARSE_GROUPING<- as.factor(dt1$COARSE_GROUPING)

# Convert Missing Values to NA for non-dates

dt1$PERCENT_COVER <- ifelse((trimws(as.character(dt1$PERCENT_COVER))==trimws("-99999")),NA,dt1$PERCENT_COVER)               
suppressWarnings(dt1$PERCENT_COVER <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$PERCENT_COVER))==as.character(as.numeric("-99999"))),NA,dt1$PERCENT_COVER))
dt1$DENSITY <- ifelse((trimws(as.character(dt1$DENSITY))==trimws("-99999")),NA,dt1$DENSITY)               
suppressWarnings(dt1$DENSITY <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$DENSITY))==as.character(as.numeric("-99999"))),NA,dt1$DENSITY))
dt1$WM_GM2 <- ifelse((trimws(as.character(dt1$WM_GM2))==trimws("-99999")),NA,dt1$WM_GM2)               
suppressWarnings(dt1$WM_GM2 <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$WM_GM2))==as.character(as.numeric("-99999"))),NA,dt1$WM_GM2))
dt1$DM_GM2 <- ifelse((trimws(as.character(dt1$DM_GM2))==trimws("-99999")),NA,dt1$DM_GM2)               
suppressWarnings(dt1$DM_GM2 <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$DM_GM2))==as.character(as.numeric("-99999"))),NA,dt1$DM_GM2))
dt1$SFDM <- ifelse((trimws(as.character(dt1$SFDM))==trimws("-99999")),NA,dt1$SFDM)               
suppressWarnings(dt1$SFDM <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$SFDM))==as.character(as.numeric("-99999"))),NA,dt1$SFDM))
dt1$AFDM <- ifelse((trimws(as.character(dt1$AFDM))==trimws("-99999")),NA,dt1$AFDM)               
suppressWarnings(dt1$AFDM <- ifelse(!is.na(as.numeric("-99999")) & (trimws(as.character(dt1$AFDM))==as.character(as.numeric("-99999"))),NA,dt1$AFDM))
dt1$SCIENTIFIC_NAME <- as.factor(ifelse((trimws(as.character(dt1$SCIENTIFIC_NAME))==trimws("-99999")),NA,as.character(dt1$SCIENTIFIC_NAME)))
dt1$TAXON_PHYLUM <- as.factor(ifelse((trimws(as.character(dt1$TAXON_PHYLUM))==trimws("-99999")),NA,as.character(dt1$TAXON_PHYLUM)))
dt1$TAXON_CLASS <- as.factor(ifelse((trimws(as.character(dt1$TAXON_CLASS))==trimws("-99999")),NA,as.character(dt1$TAXON_CLASS)))
dt1$TAXON_ORDER <- as.factor(ifelse((trimws(as.character(dt1$TAXON_ORDER))==trimws("-99999")),NA,as.character(dt1$TAXON_ORDER)))
dt1$TAXON_FAMILY <- as.factor(ifelse((trimws(as.character(dt1$TAXON_FAMILY))==trimws("-99999")),NA,as.character(dt1$TAXON_FAMILY)))

write.csv(dt1, "data/LTER/Annual_All_Species_Biomass_at_transect_20210108.csv", row.names = F)






























