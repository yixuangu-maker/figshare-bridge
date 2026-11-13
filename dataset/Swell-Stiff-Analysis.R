#The objective of this script is to predict intermediate structural parameters and resulting physical properties
#for PVA hydrogel formulations based on swollen polymer network theory.
#The input data will be the csv file containing experimental data for all PVA hydrogel formulations
#the output will add predicted structural parameters
# 1. Degree of Polymerization Between Junctions
# 2. Junction Functionality
# 3. Initial Polymer Volume Fraction
# 4. Frequency of Chain-End Defects
# 5. Structure-Predicted Shear Modulus

#call the tidyverse library so we can use dplyr to analyze and develop the data. 
#"tidyverse" package should already be installed.
library(tidyverse)

#set working directory
homedir <- "C:/Users/nr23379/Box/PVA Hydrogel Data/PVA Stiffness Paper"
setwd(homedir)

#retrieve summary tables from each experimental group
sweldata <- read_csv("Swelling/Swel-summary-data.csv")
tensdata <- read_csv("Stiffness/Tensile/Tens-summary-data.csv")
rheodata <- read_csv("Stiffness/Rheological/Rheo-summary-data.csv")
compdata <- read_csv("Stiffness/Compressive/Comp-summary-data.csv")
mpendata <- read_csv("Stiffness/MacroPenetration/MPen-summary-data.csv")

#remove index column (...1) for each table 
sweldata %>% select(-...1) -> sweldata
tensdata %>% select(-...1) -> tensdata
rheodata %>% select(-...1) -> rheodata
compdata %>% select(-...1) -> compdata
mpendata %>% select(-...1) -> mpendata

#further trim the compdata
compdata %>% select(
  polymer, MWrange, conc, interval, quickform,
  dry_n_comp = dry_n_dma,
  dry_g = dry_g_s1,
  dry_g_sd = dry_g_s1_sd,
  dry_r2 = dry_r2_s1,
  wet_n_comp = wet_n_dma,
  wet_g = wet_g_s1,
  wet_g_sd = wet_g_s1_sd,
  wet_r2 = wet_r2_s1
) -> compdata

#join the data for swelling & stiffness
expdata <- left_join(sweldata, tensdata, by = c( "polymer", "MWrange", "conc", "interval", "quickform"))
expdata <- left_join(expdata, rheodata, by = c("polymer", "MWrange", "conc", "interval", "quickform"))
expdata <- left_join(expdata, compdata, by = c("polymer", "MWrange", "conc", "interval", "quickform"))
expdata <- left_join(expdata, mpendata, by = c("polymer", "MWrange", "conc", "interval"))

#Import the tidy network parameters data
paramdata <- read_csv("PVA-M-Param-Tidy.csv")

#join parameters with appropriate formulations
structuredata <- left_join(expdata, paramdata)

#add the initial polymer volume fraction and degree of polymerization between junctions (j) for each formulation
structuredata <- mutate(
  structuredata,
  ini_volfrac = conc / 1000,
  Nj = 2*interval / fun
)

#calculate the frequency of chain-end defects (fced)
get_fced <- function(fun, Nj, Mr, Mn) {
  a = fun*Nj*Mr/((fun-2)*Mn)
  a
}

#add fced to table for each formulation
structuredata <- mutate(
  structuredata,
  fced = get_fced(
    fun = fun,
    Nj = Nj,
    Mr = Mr,
    Mn = Mn
  )
)

#calculate the shear modulus (kPa)
#note: This assumes no temperature change between relaxed, swollen state, and deformed states
get_shear <- function(Rgas = 8314, Temp = 298, fun, fced, dens, Nj, Mr, ini, ref_rat) {
  a = Rgas*Temp*(1-2/fun)*(1-fced)*dens/(Mr*Nj)*ini*ref_rat^(-1/3)
  a
}

#add shear modulus (in kPa) to the table for each formulation
structuredata <- mutate(
  structuredata,
  pred_shear = get_shear(
    fun = fun,
    fced = fced,
    dens = dry_dens,
    Nj = Nj,
    Mr = Mr,
    ini = ini_volfrac,
    ref_rat = ref_ratio
  )
)

#export everything to csvs
write.csv(structuredata, "swell-stiff-summary.csv")

