#The goal of this script is to import a bunch of tab-separated text files that represent dynamic mechanical analyzer (DMA) experiments
#and analyze those data to extract shear moduli values
#as well as compare the effects of a dry versus wet environment and second vs first runs\
#and do a little analysis or at least visualization of the difference between the 

library(tidyverse)
library(cowplot)

#set working directory
homedir <- "C:/Users/nrich/Box/PVA Hydrogel Data/PVA Stiffness Paper/Stiffness/Compressive"
setwd(homedir)

#store subdirectory names (subdirectories already made in windows)
rawsubdir <- "Comp Raw Data"
prosubdir <- "Comp Processed Data"
wcfigsubdir <- "Comp Figures/WholeCurves"
p1figsubdir <- "Comp Figures/Phase1Figs"
s1figsubdir <- "Comp Figures/Subset1Figs"

#extract name of each file in each raw data subdirectory into an overarching file list
#and a second list with just the names of each file to make for slightly easier name assignment later.
datafolders_list  <- list.dirs(path = rawsubdir, full.names = FALSE, recursive = FALSE)
file_list_justnames <- vector()
file_list <- vector()
for (j in 1:length(datafolders_list)) {
  temp_file_list_justnames <- list.files(path = paste(rawsubdir, datafolders_list[j], sep = "/"))
  temp_file_list <- paste( rawsubdir, datafolders_list[j], temp_file_list_justnames, sep = "/")
  file_list_justnames <- c(file_list_justnames, temp_file_list_justnames)
  file_list <- c(file_list, temp_file_list)
}

#generate tibble with no designated column names but no values yet
summarycolnames <- c(
  "exp_id",
  "quickform",
  "p1_r2",
  "p1_int",
  "p1_int_sterr",
  "p1_slope",
  "p1_slope_sterr",
  "s1_r2",
  "s1_int",
  "s1_int_sterr",
  "s1_slope",
  "s1_slope_sterr"
)
summary_data <- read_csv("\n", col_names = summarycolnames, col_types = "ccnnnnnnnnnn")

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
  
  temp_data %>% mutate(
    #add column for stretch
    stretch = length / first(length),
    #add column for area (assuming starting diameter is 10 mm) (in meters)
    area = pi * 0.005 ^ 2 * 1/stretch,
    #add column for current (Cauchy) stress (kPa)
    stress = force / area / 1000,
    #add column for Neo-Hookean imcompressible solid "strain" stretch^2-1/stretch
    compstrain = -1*(stretch ^ 2 - 1/stretch)
    ) -> temp_data
  
  #Make the file name without the .txt file-type
  experiment_name <- substr(file_list_justnames[i], start = 1, stop = nchar(file_list_justnames[i])-4)
 
  #save the processed data for each experiment as a csv in the processed data folder
  write.csv(temp_data, paste(prosubdir, "/", experiment_name, ".csv", sep = ""))
  
  #create a ggplot figure showing the full stress-strain loop over the entire experiment, with the line changing color gradually over time
  #to help visualize the effect of conditioning
  g <- ggplot(temp_data, aes(x=compstrain, y = stress, colour = time)) +
    geom_path(size = 1.2) +
    theme_cowplot() +
    scale_color_gradient2(
      low = "blue",
      mid = "grey",
      high = "red",
      midpoint = 2
    ) +
    labs(x = "L^2-1/L", y = "Stress (kPa)") +
    scale_y_continuous(
      # don't expand y scale at the lower end
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  
  #save figure to whole curve figure folder
  ggsave(
    filename = paste(wcfigsubdir, "/", experiment_name, ".tiff", sep = ""),
    plot = g,
    width = 20,
    height = 16,
    units = "cm",
    dpi = 300
  )
  
  #separate data subsets based on time (0-1 min) for the compression phass for subsequent linear analysis
  #these represent the first compression phase (excludes decompression phases and second compression phase)
  phase1data <- filter(temp_data, time < 1)
  
  #create linear regression model
  phase1_model <- lm(stress ~ compstrain, data = phase1data)
  
  #create figure with actual data and linear regression comparison
  phase1_fig <- ggplot(phase1data, aes(compstrain, stress)) +
    geom_point() +
    stat_smooth(method = lm) +
    theme_cowplot() +
    labs(x = "L^2-1/L", y = "Stress (kPa)") +
    scale_y_continuous(
      # don't expand y scale at the lower end
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  
  #save figure
  ggsave(
    filename = paste(p1figsubdir, "/", experiment_name, ".tiff", sep = ""),
    plot = phase1_fig,
    width = 20,
    height = 16,
    units = "cm",
    dpi = 300
  )
  
  #append linear regression statistics from each experiment to summary data table
  p1_int <- summary(phase1_model)$coefficients[1,1]
  p1_int_sterr <- summary(phase1_model)$coefficients[1,2]
  p1_slope <- summary(phase1_model)$coefficients[2,1]
  p1_slope_sterr <- summary(phase1_model)$coefficients[2,2]
  p1_r2 <- summary(phase1_model)$r.squared
  
  #Create a subset of the upper half of strains in phase 1
  sub1data <- filter(phase1data, compstrain > 0.5*max(phase1data$compstrain))
  
  #create linear regression model for the subset dataset
  sub1_model <- lm(stress ~ compstrain, data = sub1data)
  
  #create figure with actual data and linear regression comparison
  sub1_fig <- ggplot(sub1data, aes(compstrain, stress)) +
    geom_point() +
    stat_smooth(method = lm) +
    theme_cowplot() +
    labs(x = "L^2-1/L", y = "Stress (kPa)") +
    scale_y_continuous(
      # don't expand y scale at the lower end
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  
 
  #save figure
  ggsave(
    filename = paste(s1figsubdir, "/", experiment_name, ".tiff", sep = ""),
    plot = sub1_fig,
    width = 20,
    height = 16,
    units = "cm",
    dpi = 300
  )
  
  #append linear regression statistics from each experiment and phase to summary data table
  s1_int <- summary(sub1_model)$coefficients[1,1]
  s1_int_sterr <- summary(sub1_model)$coefficients[1,2]
  s1_slope <- summary(sub1_model)$coefficients[2,1]
  s1_slope_sterr <- summary(sub1_model)$coefficients[2,2]
  s1_r2 <- summary(sub1_model)$r.squared
  
  #Add in quick formulation ID for grouped analysis in graphpad!
  quickform <- substr(experiment_name, start = 8, stop = nchar(experiment_name)-2)
  
  summary_data %>% add_row(
    exp_id = experiment_name,
    quickform = quickform,
    p1_r2 = p1_r2,
    p1_int = p1_int,
    p1_int_sterr = p1_int_sterr,
    p1_slope = p1_slope,
    p1_slope_sterr = p1_slope_sterr,
    s1_r2 = s1_r2,
    s1_int = s1_int,
    s1_int_sterr = s1_int_sterr,
    s1_slope = s1_slope,
    s1_slope_sterr = s1_slope_sterr
  ) -> summary_data
  
}

#Save summary data to a csv in the home directory
write.csv(summary_data, "perexperimentdata.csv")

#separate the experiment id
summary_data %>% separate(
  col = exp_id,
  into = c("drywet","set","polymer","MWrange", "conc","interval","replicate"),
  sep = "-",
  convert = TRUE
) -> expanded_data


#summarize the data and extract shear moduli from compressive data assuming a Neo-Hookean material model
#NOTE: Neo-Hookean material model assumption is clearly false based on stress/strain figures...
#need to look back into material models to find something that is consistent with both the tensile and compressive data
expanded_data %>%
  group_by(drywet, polymer, MWrange, conc, interval, quickform) %>%
  summarize(
    n_dma = n(),
    mean_g1 = mean(p1_slope),
    sd_g1 = sd(p1_slope),
    mean_r21 = mean(p1_r2),
    mean_gs1 = mean(s1_slope),
    sd_gs1 = sd(s1_slope),
    mean_r2s1 = mean(s1_r2)
  ) -> avgdata

avgdata %>% filter(drywet == "DRY") %>% ungroup() %>% select(
  polymer, MWrange, conc, interval, quickform,
  dry_n_dma = n_dma,
  dry_g_s1 = mean_gs1,
  dry_g_s1_sd = sd_gs1,
  dry_g_p1 = mean_g1,
  dry_g_p1_sd = sd_g1,
  dry_r2_s1 = mean_r2s1,
  dry_r2_p1 = mean_r21
)-> drydata

avgdata %>% filter(drywet == "WET") %>% ungroup %>% select(
  polymer, MWrange, conc, interval, quickform,
  wet_n_dma = n_dma,
  wet_g_s1 = mean_gs1,
  wet_g_s1_sd = sd_gs1,
  wet_g_p1 = mean_g1,
  wet_g_p1_sd = sd_g1,
  wet_r2_s1 = mean_r2s1,
  wet_r2_p1 = mean_r21
)-> wetdata

left_join (drydata, wetdata) -> finaldata

#save to statistically compare and further graph in GraphPad Prism
write.csv(finaldata, "Comp-summary-data.csv")

