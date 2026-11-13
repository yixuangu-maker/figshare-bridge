#The objective of this script is to extract structural parameters for PVA hydrogel formulations.
#This data will be used in an upcoming publication on hydrogel swelling properties.
#The end result will be one csv file of the most important summary data

#Make sure that the working directory has the original csv and this script in it.
#Example: setwd("C:/Users/nrich/Box/Spring 2020/Swelling Ratio Manuscript/Data/TidyData-RScripts")

#call the tidyverse library so we can use dplyr to analyze and develop the data. 
#"tidyverse" package should already be installed.
library(tidyverse)

#import the tidy summary data of PVA Swelling
swellingdata <- read_csv("PVA-Summary.csv")

#Import the tidy network parameters data
paramdata <- read_csv("Network-Param-Tidy.csv")

#All we need from each formulation is the mean relaxed polymer volume fraction
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

#merge each formulation's swelling data with the PVA Polymer parameters
swellingdata <- mutate(
  swellingdata,
  Polymer = "PVA"
)
structuredata <- left_join(swellingdata, paramdata)

#calculate Mc using the equilibrium swelling equation with molecular weight-based chain end defects
#and a know junction functionality
get_Mc_PVA <- function(fun, Mn, swol, Chi, V1, dens, rel) {
  a <- log(1-swol) + swol + Chi * swol ^ 2
  a <- a/((1-2/fun)*V1*dens*rel^(2/3)*swol^(1/3))
  a <- fun/((fun-2)*Mn) - a
  a <- 1/a
  a
}
#add mc to table for each formulation.
structuredata <- mutate(
  structuredata,
  Mc = get_Mc_PVA(
    fun = fun,
    Mn = Mn,
    swol = swol,
    Chi = Chi,
    V1 = V1,
    dens = dens,
    rel = rel
  )
)

#calculate the frequency of chain-end defects (fced)
get_fced <- function(fun, Mc, Mn) {
  a = fun*Mc/((fun-2)*Mn)
  a
}

#add fced to table for each formulation
structuredata <- mutate(
  structuredata,
  fced = get_fced(
    fun = fun,
    Mc = Mc,
    Mn = Mn
  )
)

#calculate the mesh size
get_ms <- function(swol, fun, len, Crat, atom, Mc, Mr) {
  a = swol^(-1/3)*((1-2/fun)*len^2*Crat*atom*Mc/Mr)^(1/2)
  a
}

#add mesh size (in nanometers) to table for each formulation
structuredata <- mutate(
  structuredata,
  ms = get_ms(
    swol = swol,
    fun = fun,
    len = len,
    Crat = Crat,
    atom = atom,
    Mc = Mc,
    Mr = Mr
  )
)

#calculate the shear modulus
get_shear <- function(Rgas = 8314, Temp = 298, fun, fced, dens, Mc, rel, swol) {
  a = Rgas*Temp*(1-2/fun)*(1-fced)*dens/Mc*rel^(2/3)*swol^(1/3)
  a
}

#add shear modulus (in kPa) to the table for each formulation
structuredata <- mutate(
  structuredata,
  shear = get_shear(
    fun = fun,
    fced = fced,
    dens = dens,
    Mc = Mc,
    rel = rel,
    swol = swol
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
write.csv(structuredata, "PVA-Structure.csv")