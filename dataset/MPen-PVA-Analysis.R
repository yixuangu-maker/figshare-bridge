#The goal of this script is to import a bunch of tab-separated text files that represent dynamic mechanical analyzer (DMA) experiments
#and analyze those data to extract shear moduli values
#as well as compare the effects of a dry versus wet environment and second vs first runs\
#and do a little analysis or at least visualization of the difference between the 

library(tidyverse)
library(cowplot)

#set working directory
homedir <- "C:/Users/nrich/Box/PVA Hydrogel Data/PVA Stiffness Paper/Stiffness/MacroPenetration"
setwd(homedir)

#store subdirectory names (subdirectories already made in windows)
rawsubdir <- "MPen-Raw-Data"
prosubdir <- "MPen-Processed-Data"
figsubdir <- "MPen-Figs"

#extract name of each file in each raw data subdirectory into an overarching file list
#and a second list with just the names of each file to make for slightly easier name assignment later.
file_list_justnames <- list.files(path = rawsubdir)
file_list <- paste( rawsubdir, file_list_justnames, sep = "/")
#generate tibble with no designated column names but no values yet
summarycolnames <- c(
  "exp_id",
  "quickform",
  "fit_r2",
  "fit_int",
  "fit_int_sterr",
  "fit_slope",
  "fit_slope_sterr"
)
summary_data <- read_csv("\n", col_names = summarycolnames, col_types = "ccnnnnn")

#for loop used to analyze each file
for(i in 1:length(file_list)){
  
  #this should work with original text files in UTF-16
  temp_data <- read.csv(file =  file_list[i], skip = 55, fileEncoding = "UTF-16", sep = "\t", header = FALSE)
  
  #convert to tibble format (for tidyverse analysis)
  temp_data <- as_tibble(temp_data)
  
  #just keep the time, force, and length data from the original file (and convert column header for clarity)
  temp_data %>% select(
    time = V1, #time is in minutes
    force = V4, #force is in newtons
    length = V6 #length is in millimeters
    ) -> temp_data 
  
  #subset the data to the first 0.6 minutes, approximately encapsulating the penetration phase
  temp_data <- filter(temp_data, time < 0.6)
  
  #assign row IDs
  temp_data <- rowid_to_column(temp_data, "ID")
  
  #filter out data points before the lowest force point
  temp_datapt <- filter(temp_data, force == min(force))
  temp_data <- filter(temp_data, ID >= temp_datapt$ID[1])
  
  #filter out data points after the highest force point
  temp_datapt <- filter(temp_data, force == max(force))
  temp_data <- filter(temp_data, ID <= temp_datapt$ID[1])
  
  temp_data %>% mutate(
    #add column for applied force
    #in milliNewtons to make sure final G (slope) will be in kPa
    normforce = (force - first(force))*1000,
    #add column for applied displacement
    normdisp = first(length) - length,
    #add column for the full Hertzian term
    #assumes tip radius = 1.5 mm and Poisson's ratio = 0.5
    #also converts so final G (slope) will be in kPa
    hertz = 1000*1000*16* .0015^(1/2) * (0.001*normdisp) ^ (3/2) / 3
    ) -> temp_data
  
  #filter out any possible erroneous less-than-zero data
  temp_data <- filter(temp_data, normforce >= 0 && normdisp >= 0)
  
  #Make the file name without the .txt file-type
  experiment_name <- substr(file_list_justnames[i], start = 1, stop = nchar(file_list_justnames[i])-4)
 
  #save the processed data for each experiment as a csv in the processed data folder
  write.csv(temp_data, paste(prosubdir, "/", experiment_name, ".csv", sep = ""))
  
  #create linear regression model
  temp_model <- lm(normforce ~ hertz, data = temp_data)
  
  #create figure with actual data and linear regression comparison
  fit_fig <- ggplot(temp_data, aes(hertz, normforce)) +
    geom_point() +
    stat_smooth(method = lm) +
    theme_cowplot() +
    labs(x = "Hertzian Indentation (mm^2)", y = "Force (mN)") +
    scale_y_continuous(
      # don't expand y scale at the lower end
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  
  #save figure
  ggsave(
    filename = paste(figsubdir, "/", experiment_name, ".tiff", sep = ""),
    plot = fit_fig,
    width = 20,
    height = 16,
    units = "cm",
    dpi = 300
  )
  
  #append linear regression statistics from each experiment to summary data table
  fit_int <- summary(temp_model)$coefficients[1,1]
  fit_int_sterr <- summary(temp_model)$coefficients[1,2]
  fit_slope <- summary(temp_model)$coefficients[2,1]
  fit_slope_sterr <- summary(temp_model)$coefficients[2,2]
  fit_r2 <- summary(temp_model)$r.squared
  
  #Add in quick formulation ID for grouped analysis in graphpad!
  quickform <- substr(experiment_name, start = 4, stop = nchar(experiment_name)-2)
  
  summary_data %>% add_row(
    exp_id = experiment_name,
    quickform = quickform,
    fit_r2 = fit_r2,
    fit_int = fit_int,
    fit_int_sterr = fit_int_sterr,
    fit_slope = fit_slope,
    fit_slope_sterr = fit_slope_sterr
  ) -> summary_data
  
}

#Save summary data to a csv in the home directory
write.csv(summary_data, "perexperimentdata.csv")

#separate the experiment id
summary_data %>% separate(
  col = exp_id,
  into = c("set","polymer","MWrange", "conc","interval","replicate"),
  sep = "-",
  convert = TRUE
) -> expanded_data


#summarize the data and extract shear moduli from compressive data assuming a Neo-Hookean material model
#NOTE: Neo-Hookean material model assumption is clearly false based on stress/strain figures...
#need to look back into material models to find something that is consistent with both the tensile and compressive data
expanded_data %>%
  group_by(polymer, MWrange, conc, interval, quickform) %>%
  summarize(
    n_mpen = n(),
    mean_g = mean(fit_slope),
    sd_g = sd(fit_slope),
    mean_r2 = mean(fit_r2)
  ) -> avgdata

#save to statistically compare and further graph in GraphPad Prism
write.csv(avgdata, "MPen-summary-data.csv")

