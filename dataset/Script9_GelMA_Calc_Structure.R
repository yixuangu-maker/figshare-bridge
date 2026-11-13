#The objective of this script is to extract structural parameters 
#for GelMA hydrogel formulations.
#This data will be used in an upcoming publication on hydrogel swelling 
#properties.
#The end result will be one csv file of the most important summary data

#Make sure that the working directory has the original csv and this 
#script in it.
#Example: setwd("C:/Users/nrich/Box/Spring 2020/Swelling Ratio Manuscript/Data/TidyData-RScripts")

#call the tidyverse library so we can use dplyr to analyze and develop the data. 
#"tidyverse" package should already be installed.
library(tidyverse)

#import the tidy summary data of GelMA Swelling
swellingdata <- read_csv("GelMA-Calc-Summary.csv")

#Import the tidy network parameters data
paramdata <- read_csv("Network-Param-Tidy.csv")

#simplify naming conventions
swellingdata %>%
  rename(
    pvf_0 = "PVF-0",
    ci = "CI",
    rel = "PVF-R",
    ratio = "RefRatio"
  ) ->
  swellingdata

#merge each formulation's swelling data with the GelMA Polymer parameters
#AND add swollen polymer volume fraction based on reference polymer volume fraction and reference ratio
swellingdata <- mutate(
  swellingdata,
  Polymer = "GelMA",
  swol = rel / ratio
)
structuredata <- left_join(swellingdata, paramdata)

#calculate Mc using the equilibrium swelling equation with molecular 
#weight-based chain end defects
#and an unknow junction functionality (but it is probably very high)
get_Mc_GelMA <- function(fun, Mn, swol, Chi, V1, dens, rel) {
  a <- log(1-swol) + swol + Chi * swol ^ 2
  a <- a/((1-2/fun)*V1*dens*rel^(2/3)*swol^(1/3))
  a <- fun/((fun-2)*Mn) - a
  a <- 1/a
  a
}
#add mc to table for each formulation.
structuredata <- mutate(
  structuredata,
  Mc = get_Mc_GelMA(
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
write.csv(structuredata, "GelMA-Calc-Structure.csv")