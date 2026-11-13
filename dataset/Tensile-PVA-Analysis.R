#This script is for the analysis of hydrogel uniaxial tensile data
#For each data set, inputs include
# 1) force, time, and displacement data from the CellScale program
# 2) tracking coordinates vs time data from the CellScale program
# 3) thickness of the sample, measured independently using calipers

#The warning for "missing column names" are not important, just means some extra white-space is imported on some of the Poisson's
#ratio tables because of some quirk of reading the csv files, but that does not affect calculations.

#call the tidyverse library so we can use dplyr to analyze and develop the data. 
#"tidyverse" package should already be installed.
library(tidyverse)
library(cowplot) #for making decently themed plots

#setwd
homedir <- "C:/Users/nrich/Box/PVA Hydrogel Data/PVA Stiffness Paper/Stiffness/Tensile"
setwd(homedir)
#main subdirectories with data
set1subdir <- "Tens-S1-PVA-M"
set2subdir <- "Tens-S2-PVA-M"
#listing all sub-subdirectories within the data subdirectories. Each contains data from a tensile experiment.
set1subdir_list_justnames  <- list.dirs(path = set1subdir, full.names = FALSE, recursive = FALSE)
set1subdir_list <- paste(set1subdir, set1subdir_list_justnames, sep = "/")
set2subdir_list_justnames  <- list.dirs(path = set2subdir, full.names = FALSE, recursive = FALSE)
set2subdir_list <- paste(set2subdir, set2subdir_list_justnames, sep = "/")

#make one overarching list of data-containing folders that can be accessed from the homedir
datadir_list <- c(set1subdir_list, set2subdir_list)
datadir_list_justnames <- c(set1subdir_list_justnames, set2subdir_list_justnames)

#Poisson figures folder, will save all Poisson's Ratio figures here
poissonfigsubdir <- "Poisson Figs"
#image-based analysis figures folder, will save all image-based analysis plots here
imagefigsubdir <- "Image-Based Figs"
#displacement-based analysis figures folder, will save all displacement-based analysis plots here
dispfigsubdir <- "Disp-Based Figs"
#thickness raw data source
set1thick <- read_csv("Thickness/S1-PVA-M-Thick-Tidy.csv")
set2thick <- read_csv("Thickness/S2-PVA-M-Thick-Tidy.csv")

#creates a simple "thickness" table that can be left-joined onto each dataset when needed.
thick <- bind_rows(set1thick, set2thick)
thick %>% group_by(set, polym, MWrange, conc, interval, sample) %>%
  summarize(
    count = n(),
    avgthick = mean(`thickness (mm)`)
  ) -> thick

thick %>% mutate(
  exp_id = paste(set, polym, MWrange, conc, interval, sample, sep = "-")
) -> thick
thick %>% ungroup() %>% select(exp_id, avgthick) -> thick

#generate tibble with no designated column names but no values yet
summarycolnames <- c(
  "exp_id",
  "poisson_mean",
  "numpoints", #this is the number of datapoints that are actually used in the linear model
  "avgthickness",
  "image_r2",
  "image_int",
  "image_int_sterr",
  "image_slope",
  "image_slope_sterr",
  "disp_r2",
  "disp_int",
  "disp_int_sterr",
  "disp_slope",
  "disp_slope_sterr"
)
perexperiment_summarydata <- read_csv("\n", col_names = summarycolnames, col_types = "cnnnnnnnnnnnnn")

#loops over each experiment folder, analyzing data and generating figures for each experiment.
#Key data from each experiment is added to a list of per-experiment observations
#which are eventually further summarized into a per-formulation analysis of tensile results.
for (i in 1:length(datadir_list)) {
  #the experiment name is taken from the experiment folder's name
  exp_id <- datadir_list_justnames[i]
  #import the relevant force-displacement csv
  forceraw <- read_csv(paste(datadir_list[i],"F.csv", sep = "/"))
  #trim to the 1-Stretch cycle (this is the cycle that all the image-based analysis was done on)
  forceraw %>% filter(Cycle == "1-Stretch", SetName == "Experimental") -> forceraw
  
  #import the relevant tracking csv
  #for tracking, we want one dataset for the principal strain 1 (main stretch)
  #and a second dataset for the principal strain 2 (lateral compression)
  #to calculate Poisson's ratio
  #note, the skips and n-max values are hard-coded, but all datasets should be the same in this regard. Will see when testing the full loop.
  strain1raw <- read_csv(paste(datadir_list[i],"C.csv", sep = "/"), skip = 96, n_max = 91)
  strain2raw <- read_csv(paste(datadir_list[i],"C.csv", sep = "/"), skip = 191, n_max = 91)
  
  #We calculate Poisson's ratio from the two strain raw datasets
  left_join(strain1raw, strain2raw, by = "Time(S)") -> prepoisson
  for ( j in 1:16) {
    col1 <- paste("Strain1_",j, sep = "")
    col2 <- paste("Strain2_",j, sep = "")
    col3 <- paste("Poisson_",j, sep = "")
    prepoisson %>% mutate(
      !!as.symbol(col3) := -1*!!as.symbol(col2)/ !!as.symbol(col1)  #I do not know why, but the := made this work. 
      #Something to do w/ "Glue" formatting for dplyr (borrowed from Stackoverflow)
    ) -> prepoisson
  }
  #add an explicit row number (image number) for later use
  prepoisson %>% mutate(image_number = row_number()) -> prepoisson
  #get rid of first row which should have a bunch of NaN values from dividing by zero.
  prepoisson <- prepoisson[-1,]
  #average over poisson values for each row
  prepoisson %>% mutate(Poisson_mean = rowMeans(select(., starts_with("Poisson")))) -> prepoisson
  
  #Problem here! A large outlier can completely throw off the dataset!
  #created another column with filtered Poisson's ratios (0 < x < 1), 
  #Note these limits are a little arbitrary. That's why I made it a separate column so we can see which data were omitted...
  prepoisson %>% mutate(
    no_outliers = replace(Poisson_mean, Poisson_mean > 1 | Poisson_mean < 0, NA)
  ) -> prepoisson
  
  #average over column (over time)
  poissonmean <- mean(prepoisson$no_outliers, na.rm = TRUE) #na.rm since no_outliers might have a bunch of NA values
  #add this to prepoisson datatable for reference
  prepoisson %>% mutate(avgpoissonratio = poissonmean) -> prepoisson
  #Generate a plot of how the average Poisson's ratio changes over time
  g <- ggplot(prepoisson, aes(x=image_number , y = Poisson_mean)) +
    geom_hline(yintercept=0.5, linetype="dashed", color = "red", size = 1.5) + #red line is at ideal 0.5
    geom_hline(yintercept = poissonmean, linetype="dashed", color = "blue", size = 1.5) + #blue line represents average of all data
    geom_path(size = 1.2) +
    theme_cowplot() +
    labs(x = "Image Number", y = "Poisson's Ratio") +
    scale_y_continuous(
      # don't expand y scale at the lower end
      expand = expansion(mult = c(0, 0.05)),
      limits = c(0,1)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  
  #save figure to Poisson folder
  ggsave(
    filename = paste(poissonfigsubdir, "/",exp_id,".tiff", sep = ""),
    plot = g,
    width = 20,
    height = 16,
    units = "cm",
    dpi = 300
  )
  
  #save Poisson calc file to data folder
  write.csv(prepoisson, paste(datadir_list[i], "poissoncalc.csv", sep = "/"))
  
  ##Analyzing strain-force data
  #average over the first principal strain
  strain1raw %>% mutate(strain1mean = rowMeans(select(., starts_with("Strain1")))+1) -> avgstrain
  #select only the time ID and the mean strain
  avgstrain %>% select(Time_S = `Time(S)`, strain1mean) -> avgstrain
  #join to the force data
  left_join(forceraw, avgstrain, by = "Time_S") -> forceraw
  #clear data with negative forces
  forceraw %>% filter(Force_N >= 0) -> forceraw
  
  #add exp_id to forceraw table
  forceraw %>% mutate(exp_id = exp_id) -> forceraw
  #left-join thickness to forceraw
  left_join(forceraw, thick) -> forceraw

  
  #Normalize the stretch to the lowest value with non-negative force
  #also mutate in all the other relevant terms!
  forceraw %>% mutate(
    cs_area = avgthick * 10, #assuming gauge width of 10 mm based on the dimensions of the dogbone die cutter
    stretch = strain1mean / strain1mean[1],
    area_mm = cs_area / stretch,
    stress_kPa = Force_N/area_mm *1000,
    eng_stress = Force_N/cs_area * 1000,
    disp_stretch = (30+Displacement_mm)/30,
    NHTerm = stretch^2 - 1/stretch,
    NHTerm_disp = disp_stretch^2 - 1/disp_stretch
  ) -> forceraw
  
  #create linear models for image-based and displacement-based data
  image_model <- lm(stress_kPa ~ NHTerm, data = forceraw)
  disp_model <- lm(eng_stress ~ NHTerm_disp, data = forceraw)
  
  #create figure with actual data and linear regression comparison for each model
  image_fig <- ggplot(forceraw, aes(NHTerm, stress_kPa)) +
    geom_point() +
    stat_smooth(method = lm) +
    theme_cowplot() +
    labs(x = "L^2-1/L", y = "Stress (kPa)") +
    scale_y_continuous(
      # don't expand y scale at the lower end
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  
  disp_fig <- ggplot(forceraw, aes(NHTerm_disp, eng_stress)) +
    geom_point() +
    stat_smooth(method = lm) +
    theme_cowplot() +
    labs(x = "L^2-1/L", y = "Engineering Stress (kPa)") +
    scale_y_continuous(
      # don't expand y scale at the lower end
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  
  
  #save figures for each model
  ggsave(
    filename = paste(imagefigsubdir, "/",exp_id,".tiff", sep = ""),
    plot = image_fig,
    width = 20,
    height = 16,
    units = "cm",
    dpi = 300
  )
  
  ggsave(
    filename = paste(dispfigsubdir, "/",exp_id,".tiff", sep = ""),
    plot = disp_fig,
    width = 20,
    height = 16,
    units = "cm",
    dpi = 300
  )
  
  #save the developed forceraw data to the folder
  write.csv(forceraw, paste(datadir_list[i], "stressstraincalc.csv", sep = "/"))
  
  #append linear regression statistics from each method
  image_int <- summary(image_model)$coefficients[1,1]
  image_int_sterr <- summary(image_model)$coefficients[1,2]
  image_slope <- summary(image_model)$coefficients[2,1]
  image_slope_sterr <- summary(image_model)$coefficients[2,2]
  image_r2 <- summary(image_model)$r.squared
  
  disp_int <- summary(disp_model)$coefficients[1,1]
  disp_int_sterr <- summary(disp_model)$coefficients[1,2]
  disp_slope <- summary(disp_model)$coefficients[2,1]
  disp_slope_sterr <- summary(disp_model)$coefficients[2,2]
  disp_r2 <- summary(disp_model)$r.squared
  
  #also the number of time points used in the analysis
  numpoints <- nrow(forceraw)
  
  #need to add everything into a row of the overall per-experiment analysis table
  perexperiment_summarydata %>% add_row(
    exp_id = exp_id,
    poisson_mean = poissonmean,
    numpoints = numpoints,
    avgthickness = mean(forceraw$avgthick),
    image_r2 = image_r2,
    image_int = image_int,
    image_int_sterr = image_int_sterr,
    image_slope = image_slope,
    image_slope_sterr = image_slope_sterr,
    disp_r2 = disp_r2,
    disp_int = disp_int,
    disp_int_sterr = disp_int_sterr,
    disp_slope = disp_slope,
    disp_slope_sterr = disp_slope_sterr
  ) -> perexperiment_summarydata
}

#save data in this per-experiment state
write.csv(perexperiment_summarydata, "perexperimentdata.csv")

#break down experiment IDs for further grouping
perexperiment_summarydata %>% separate(
  col = exp_id,
  into = c("set","polymer","MWrange", "conc","interval","sample"),
  sep = "-",
  convert = TRUE
) -> expanded_data

#create a row for quick formulation reference later
expanded_data %>% mutate(
  quickform = paste(polymer, MWrange, conc, interval, sep = "-")
) -> expanded_data

expanded_data %>%
  group_by(polymer, MWrange, conc, interval, quickform) %>%
  summarize(
    n_tens = n(),
    poissonratio = mean(poisson_mean),
    sd_poissonratio = sd(poisson_mean),
    thickness = mean(avgthickness),
    image_shear = mean(image_slope),
    sd_image_shear = sd(image_slope),
    disp_shear = mean(disp_slope),
    sd_disp_slope = sd(disp_slope)
  ) -> tensilesummarydata

#save to statistically compare and further graph in GraphPad Prism
write.csv(tensilesummarydata, "Tens-summary-data.csv")
