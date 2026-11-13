#The objective of this script is to extract formulation-dependent and batch-dependent swelling parameters and dry densities in experimental PVA hydrogels
#This data will be used in an upcoming publication on hydrogel swelling properties.
#The end result will be two csv files of the most important summary data, for use in GraphPad Prism to generate nice figures of the data.

#Make sure that the working directory has the original csv and this script in it.

#call the tidyverse library so we can use dplyr to analyze and develop the data. 
#"tidyverse" package should already be installed.
library(tidyverse)

#import the tidy data of PVA Swelling
swellingdata <- read_csv("PVA-Swelling-Tidy.csv")

#cut out unnecessary date data (for personal record-keeping)
swellingdata %>%
  select(
    -"Date",
    -"Replicate"
  ) ->
  swellingdata


#Replace each column header with an abbreviation for my own sanity
swellingdata %>%
  rename(
    set = "Set",
    pvf_0 = "Initial Polymer Volume Fraction",
    ci = "Crosslinking Interval",
    rel_mass = "Relaxed Air",
    rel_hept = "Relaxed Heptane",
    swo_mass = "Swollen Air",
    swo_hept = "Swollen Heptane",
    dry_mass = "Dried Air",
    dry_hept = "Dried Heptane"
  ) ->
  swellingdata

#remove incomplete sets
swellingdata %>%
  drop_na() ->
  swellingdata

#function to get volume from buoyancy weighing measurements
#default heptane density is based on 5 measurements of the mass of 1 mL of heptane.
get_vol_from_air_hept <- function(air_mass, hept_mass, hept_dens=0.6635) {
  (air_mass - hept_mass)/hept_dens
}

#adds volume columns to swelling data table for each state
swellingdata <- mutate(
  swellingdata, 
  rel_vol = get_vol_from_air_hept(
    rel_mass, 
    rel_hept
  ),
  swo_vol = get_vol_from_air_hept(
    swo_mass, 
    swo_hept
  ),
  dry_vol = get_vol_from_air_hept(
    dry_mass, 
    dry_hept
  )
)

#remove apparent masses in heptane since those are no longer relevant
swellingdata <- select(
  swellingdata,
  -rel_hept,
  -swo_hept,
  -dry_hept
)

#add dry density data for each observation
swellingdata <- mutate(
  swellingdata, 
  dry_dens = dry_mass / dry_vol
)  

#convert each swollen mass and volume to mass and volume fractions based on dry data
swellingdata <- mutate(
  swellingdata, 
  rel_massfrac = dry_mass / rel_mass,
  swo_massfrac = dry_mass / swo_mass,
  rel_volfrac = dry_vol / rel_vol,
  swo_volfrac = dry_vol / swo_vol
)

#remove remaining absolute data points
swellingdata <- select(
  swellingdata,
  -rel_mass,
  -swo_mass,
  -dry_mass,
  -rel_vol,
  -swo_vol,
  -dry_vol
)

#additionally want the reference ratio for the swollen states
swellingdata <- mutate(
  swellingdata, 
  ref_ratio = rel_volfrac / swo_volfrac
)

#now to summarize everything by batch and formulation
batchdata <- data.frame()

swellingdata %>%
  group_by(set, pvf_0, ci) %>%
  summarize(
    n = n(),
    mean_rel_volfrac = mean(rel_volfrac),
    sd_rel_volfrac = sd(rel_volfrac),
    mean_swo_volfrac = mean(swo_volfrac),
    sd_swo_volfrac = sd(swo_volfrac),
    mean_rel_massfrac = mean(rel_massfrac),
    sd_rel_massfrac = sd(rel_massfrac),
    mean_swo_massfrac = mean(swo_massfrac),
    sd_swo_massfrac = sd(swo_massfrac),
    mean_dry_dens = mean(dry_dens),
    sd_dry_dens = sd(dry_dens),
    mean_ref_ratio = mean(ref_ratio),
    sd_ref_ratio = sd(ref_ratio)
  ) ->
  batchdata

#same summary, but the two batches are lumped together
formdata <- data.frame()

swellingdata %>%
  group_by(pvf_0, ci) %>%
  summarize(
    n = n(),
    mean_rel_volfrac = mean(rel_volfrac),
    sd_rel_volfrac = sd(rel_volfrac),
    mean_swo_volfrac = mean(swo_volfrac),
    sd_swo_volfrac = sd(swo_volfrac),
    mean_rel_massfrac = mean(rel_massfrac),
    sd_rel_massfrac = sd(rel_massfrac),
    mean_swo_massfrac = mean(swo_massfrac),
    sd_swo_massfrac = sd(swo_massfrac),
    mean_dry_dens = mean(dry_dens),
    sd_dry_dens = sd(dry_dens),
    mean_ref_ratio = mean(ref_ratio),
    sd_ref_ratio = sd(ref_ratio)
  ) ->
  formdata

#calculate the ideal relationship between mass fraction and volume fraction
get_ideal_volfrac_from_massfrac <- function(mean_massfrac, polymer_density, solvent_density = 1.00) {
  1/((1/mean_massfrac-1)*(polymer_density/solvent_density)+1)
}

#for the no-batch data, include columns for the ideal reference mass fraction and the ideal swollen mass fraction
formdata <- mutate(
  formdata, 
  ideal_rel_volfrac = get_ideal_volfrac_from_massfrac(
    mean_massfrac = mean_rel_massfrac,
    polymer_density = 1.27 #Should be consistent with the Network-Param-Tidy "dens" value for PVA
  ),
  ideal_swo_volfrac = get_ideal_volfrac_from_massfrac(
    mean_massfrac = mean_swo_massfrac,
    polymer_density = 1.27 #Should be consistent with the Network-Param-Tidy "dens" value for PVA
  )
)

#export batchdata and formdata for visual/statistical analysis in excel then GraphPad Prism
write.csv(batchdata, "PVA-Batch-Summary.csv")
write.csv(formdata, "PVA-Summary.csv")
