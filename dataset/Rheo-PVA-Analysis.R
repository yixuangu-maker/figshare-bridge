#Rheo-Analysis-Script
#Specifically for analyzing the Sandblasted wet and dry series of PVA-SML-S1

library(tidyverse)

#set working directory
homedir <- "C:/Users/nrich/Box/PVA Hydrogel Data/PVA Stiffness Paper/Stiffness/Rheological"
setwd(homedir)

rawsubdir <- "HalfN-PVA-M-Rheo"

#extract name of each file in a raw data subdirectory
#and a second list with just the names of each file to make for slightly easier name assignment later.
file_list_justnames <- list.files(path = rawsubdir)
file_list <- paste( rawsubdir, file_list_justnames, sep = "/")

#generate summary tibble with designated column names but no values yet
#will fill this in the following for loop
summarycolnames <- c(
  "exp_id",
  "quickform",
  "stormod",
  "stormod_sd",
  "lossmod",
  "lossmod_sd",
  "stor_loss_frac"
)
summary_data <- read_csv("\n", col_names = summarycolnames, col_types = "ccnnnnn")

for(i in 1:length(file_list)){
  
  #temporary storage for each file to extract relevant data
  temp_data <- read_csv(file = file_list[i], skip = 10)
  
  #just keep the storage modulus, loss modulus, and strain data from the original file (and convert column header for clarity)
  temp_data %>% select(
    stormod = Pa, #storage modulus is in pascals
    lossmod = Pa_1, #loss modulus is in pascals
    strain = '%' #strain is in percent
  ) -> temp_data 
  
  #convert storage modulus and loss modulus to kPa for consistency with other measurement methods
  temp_data %>% mutate(
    stormod = stormod / 1000,
    lossmod = lossmod / 1000
  ) -> temp_data
  
  #extract first five real storage modulus values
  temp_data %>%
    select(stormod) %>%
    slice(n=1:5) ->
    temp_stormod_data
  
  #extract average statistics for the first five values
  stormod_mean <- mean(temp_stormod_data$stormod)
  stormod_sd <- sd(temp_stormod_data$stormod)
  
  #extract first five real loss modulus values
  #note, may differ in location from storage modulus values due to unrealistic negative values
  temp_data %>%
    select(lossmod) %>%
    filter(lossmod > 0) %>%
    slice(n=1:5) ->
    temp_lossmod_data
  
  #extract average statistics for the first five loss modulus values
  lossmod_mean <- mean(temp_lossmod_data$lossmod)
  lossmod_sd <- sd(temp_lossmod_data$lossmod)
  
  #extract ratio of loss modulus to storage modulus
  stor_loss_frac = lossmod_mean / stormod_mean
  
  #Make the file name without the .csv file-type
  experiment_name <- substr(file_list_justnames[i], start = 7, stop = nchar(file_list_justnames[i])-4)
  
  #Add in quick formulation ID for grouped analysis in graphpad!
  quickform <- substr(experiment_name, start = 8, stop = nchar(experiment_name)-2)
  
  #append linear regression statistics from each experiment to summary data table
  summary_data %>% add_row(
    exp_id = experiment_name,
    quickform = quickform,
    stormod = stormod_mean,
    stormod_sd = stormod_sd,
    lossmod = lossmod_mean,
    lossmod_sd = lossmod_sd,
    stor_loss_frac = stor_loss_frac
  ) -> summary_data
  
}

#Save summary data to a csv in the home directory
write.csv(summary_data, "perexperimentdata.csv")


#separate the second half of the experiment id
summary_data %>% separate(
  col = exp_id,
  into = c("drywet", "set","polymer","MWrange","conc","interval","replicate"),
  sep = "-",
  convert = TRUE
) -> expanded_data

#summarize the data and extract formulation-averaged shear moduli, loss moduli, and the ratio between the two
expanded_data %>%
  group_by(drywet, polymer, MWrange, conc, interval, quickform) %>%
  summarize(
    n_rheo = n(),
    mean_stormod = mean(stormod),
    sd_stormod = sd(stormod),
    mean_lossmod = mean(lossmod),
    sd_lossmod = sd(lossmod),
    mean_ve_frac = mean(stor_loss_frac),
    sd_ve_frac = sd(stor_loss_frac)
  ) -> avgdata


avgdata %>% filter(drywet == "Dry") %>% ungroup() %>% select(
  polymer, MWrange, conc, interval, quickform,
  dry_n_rheo = n_rheo,
  dry_stormod = mean_stormod,
  dry_stormod_sd = sd_stormod,
  dry_lossmod = mean_lossmod,
  dry_lossmod_sd = sd_lossmod,
  dry_vefrac = mean_ve_frac,
  dry_vefrac_sd = sd_ve_frac
)-> drydata

avgdata %>% filter(drywet == "Wet") %>% ungroup %>% select(
  polymer, MWrange, conc, interval, quickform,
  wet_n_rheo = n_rheo,
  wet_stormod = mean_stormod,
  wet_stormod_sd = sd_stormod,
  wet_lossmod = mean_lossmod,
  wet_lossmod_sd = sd_lossmod,
  wet_vefrac = mean_ve_frac,
  wet_vefrac_sd = sd_ve_frac
)-> wetdata

left_join (drydata, wetdata) -> finaldata

#save to statistically compare in GraphPad Prism
write.csv(finaldata, "Rheo-summary-data.csv")
