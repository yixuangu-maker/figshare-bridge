#The objective of this script is to extract structural parameters 
#for PEGDA hydrogel formulations.
#This data will be used in an upcoming publication on hydrogel swelling 
#properties.
#The end result will be one csv file of the most important summary data

#Make sure that the working directory has the original csv and this 
#script in it.
#Example: setwd("C:/Users/nrich/Box/Spring 2020/Swelling Ratio Manuscript/Data/TidyData-RScripts")

#call the tidyverse library so we can use dplyr to analyze and develop the data. 
#"tidyverse" package should already be installed.
library(tidyverse)

#import the tidy summary data of PEGDA Swelling
swellingdata <- read_csv("PEGDA-Summary.csv")

#Import the tidy network parameters data
paramdata <- read_csv("Network-Param-Tidy.csv")

#Import the PEGDA DOF Data
dofdata <- read_csv("PEGDA-FCED-Tidy.csv")
dofdata %>%
  rename(
    fced = FCED
  ) -> dofdata

#All we need from each formulation is the mean reference polymer volume fraction
#and the mean swollen polymer volume fraction
swellingdata %>%
  select(
    pvf_0,
    ci,
    mean_rel_volfrac,
    mean_swo_volfrac,
    mean_ref_ratio
  ) ->
  swellingdata

#further simplify naming conventions
swellingdata %>%
  rename(
    rel = mean_rel_volfrac,
    swol = mean_swo_volfrac,
    ratio = mean_ref_ratio
  ) ->
  swellingdata

#merge each formulation's swelling data with the PEGDA Polymer parameters
#need to round PVF-0 in swellingdata and dofdata for a successful left-join
swellingdata <- mutate(
  swellingdata,
  Polymer = "PEGDA",
  pvf_0 = round(pvf_0, 3)
)

dofdata <- mutate(
  dofdata,
  pvf_0 = round(pvf_0, 3)
)


structuredata <- left_join(swellingdata, paramdata)
structuredata <- left_join(structuredata, dofdata)

#calculate Mc using the equilibrium swelling equation the assumption of
#zero chain-end defects
#and an unknow junction functionality (but it is probably very high)
get_Mc_PEGDA <- function(swol, Chi, V1, dens, rel, fced) {
  a <- log(1-swol) + swol + Chi * swol ^ 2
  a <- a/(V1*dens*rel^(2/3)*swol^(1/3))
  a <- (-1)*a*(1-fced)
  a <- 1/a
  a
}
#add mc to table for each formulation.
structuredata <- mutate(
  structuredata,
  Mc = get_Mc_PEGDA(
    swol = swol,
    Chi = Chi,
    V1 = V1,
    dens = dens,
    rel = rel,
    fced = fced
  )
)

#calculate the mesh size
get_ms_nofun <- function(swol, len, Crat, atom, Mc, Mr) {
  a = swol^(-1/3)*(len^2*Crat*atom*Mc/Mr)^(1/2)
  a
}

#add mesh size (in nanometers) to table for each formulation
structuredata <- mutate(
  structuredata,
  ms = get_ms_nofun(
    swol = swol,
    len = len,
    Crat = Crat,
    atom = atom,
    Mc = Mc,
    Mr = Mr
  )
)

#calculate the shear modulus
get_shear_nofun <- function(Rgas = 8314, Temp = 298, dens, Mc, rel, swol, fced) {
  a = Rgas*Temp*dens/Mc*rel^(2/3)*swol^(1/3)*(1-fced)
  a
}

#add shear modulus (in kPa) to the table for each formulation
structuredata <- mutate(
  structuredata,
  shear = get_shear_nofun(
    dens = dens,
    Mc = Mc,
    rel = rel,
    swol = swol,
    fced = fced
  )
)

#organize all the relevant formulation-dependent data
#(remove all the standard network parameters)
structuredata %>%
  select(
    Polymer,
    pvf_0,
    ci,
    rel,
    swol,
    ratio,
    Mc,
    fced,
    ms,
    shear
  ) ->
  structuredata

#export everything to a csv.
write.csv(structuredata, "PEGDA-Structure.csv")