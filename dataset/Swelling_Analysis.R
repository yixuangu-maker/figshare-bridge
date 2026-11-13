#The objective of this script is to extract formulation-dependent and batch-dependent swelling parameters and dry densities in experimental PVA hydrogels
#This data will be used in an upcoming publication on hydrogel swelling properties.
#The end result will be two csv files of the most important summary data, for use in GraphPad Prism to generate nice figures of the data.

#Make sure that the working directory has the original csv and this script in it.

#call the tidyverse library so we can use dplyr to analyze and develop the data. 
#"tidyverse" package should already be installed.
library(tidyverse)

homedir <- "C:/Users/nrich/Box/PVA Hydrogel Data/PVA Stiffness Paper/Swelling"
setwd(homedir)

#import the tidy data of PVA Swelling
rawsubdir <- "Tidy Swelling Data"

file_list_justnames <- list.files(path = rawsubdir)
file_list <- paste( rawsubdir, file_list_justnames, sep = "/")

#bind all swelling analysis files into one master table
if(length(file_list) > 1) {
  sweldata <- read_csv(file_list[1])
  
  for (i in 2:length(file_list)) {
    sweldata <- bind_rows(sweldata, read_csv(file_list[i]))
  }
  
}else if (length(file_list) == 1){
  sweldata <- read_csv(file_list[1])
}
#not gonna code a result for if file list is not at least one... This whole script is useless then and will quickly throw errors.

#function to get volume from buoyancy weighing measurements
#default heptane density is based on 5 measurements of the mass of 1 mL of heptane.
get_vol_from_air_hep <- function(air_mass, hep_mass, hep_dens=0.6635) {
  (air_mass - hep_mass)/hep_dens
}

#adds volume columns to swelling data table for each state
sweldata <- mutate(
  sweldata, 
  rel_vol = get_vol_from_air_hep(
    rel_air, 
    rel_hep
  ),
  swo_vol = get_vol_from_air_hep(
    swo_air, 
    swo_hep
  ),
  dry_vol = get_vol_from_air_hep(
    dry_air, 
    dry_hep
  )
)

#add dry density data for each observation
sweldata <- mutate(
  sweldata, 
  dry_dens = dry_air / dry_vol
)

#convert each swollen mass and volume to mass and volume fractions based on dry data
sweldata <- mutate(
  sweldata, 
  rel_massfrac = dry_air / rel_air,
  swo_massfrac = dry_air / swo_air,
  rel_volfrac = dry_vol / rel_vol,
  swo_volfrac = dry_vol / swo_vol
)

#additionally want the reference ratio for the swollen states
sweldata %>% mutate(
  ref_ratio = rel_volfrac / swo_volfrac,
  quickform = paste(polymer, MWrange, conc, interval, sep = "-") #added in quickform here
) -> sweldata

sweldata %>%
  group_by(polymer, MWrange, conc, interval, quickform) %>%
  summarize(
    n_swel = n(),
    dry_dens_avg = mean(dry_dens),
    rel_volfrac_avg = mean(rel_volfrac),
    swo_volfrac_avg = mean(swo_volfrac),
    ref_ratio_avg = mean(ref_ratio),
    sd_dry_dens = sd(dry_dens),
    sd_rel_volfrac = sd(rel_volfrac),
    sd_swo_volfrac = sd(swo_volfrac),
    sd_ref_ratio = sd(ref_ratio)
  ) -> summary_data


summary_data %>% select(
  polymer, MWrange, conc, interval, quickform,
  n_swel,
  dry_dens = dry_dens_avg,
  sd_dry_dens,
  rel_volfrac = rel_volfrac_avg,
  sd_rel_volfrac,
  swo_volfrac = swo_volfrac_avg,
  sd_swo_volfrac,
  ref_ratio = ref_ratio_avg, 
  sd_ref_ratio
) -> summary_data

#save to statistically compare and further graph in GraphPad Prism
write.csv(summary_data, "Swel-summary-data.csv")
