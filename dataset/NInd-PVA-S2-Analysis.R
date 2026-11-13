#The goal of this script is to import a bunch of tab-separated text files that represent dynamic mechanical analyzer (DMA) experiments
#and analyze those data to extract shear moduli values

library(tidyverse)

#set working directory
homedir <- "C:/Users/nrich/Box/PVA Hydrogel Data/PVA Stiffness Paper/Stiffness/NanoIndentation/NR-PVA-S2"
setwd(homedir)

#store subdirectory names (subdirectories already made in windows)
rawsubdir <- "Raw Data"

#extract name of each file in each raw data subdirectory into an overarching file list
#and a second list with just the names of each file to make for slightly easier name assignment later.
dir_list <- list.dirs(path = rawsubdir, recursive = FALSE)

#generate tibble with no designated column names but no values yet
summarycolnames <- c(
  "exp_id",
  "quickform",
  "n_NInd",
  "g_avg",
  "g_sd"
)
summary_data <- read_csv("\n", col_names = summarycolnames, col_types = "ccnnn")

#for loop used to analyze each file
for(i in 1:length(dir_list)){
  
  #this should work with original text files in UTF-16
  temp_data <- read.csv(file = paste(
    dir_list[i], 
    "NR-PVA S-1 E-eff vs XY position.txt", 
    sep = "/"), sep = "\t", 
    skip = 2, 
    header = FALSE, 
    col.names = c("A", "B", "C")
    )
  
  #convert to tibble format (for tidyverse analysis)
  temp_data <- as_tibble(temp_data)
  
  #pivot longer to simplify mutations and averaging
  temp_data %>% pivot_longer(
    cols = everything(),
    names_to = "x",
    values_to = "E"
  ) -> temp_data
  
  temp_data %>% mutate(
    #add column for g (keeping assumption that Poisson's Ratio = 0.5, also converting to kPa)
    g = E / 3000
    ) -> temp_data
  
  #summarize the data from this folder/matrix scan
  temp_data %>% summarize(
    n_NInd = n(),
    g_avg = mean(g),
    g_sd = sd(g)
  ) -> temp_data
  
  #define experiment ID from folder name
  exp_id <- substr(dir_list[i], start = 10, stop = nchar(dir_list[i]))
  quickform <- substr(exp_id, start = 1, stop = nchar(exp_id) - 7)
  
  summary_data %>% add_row(
    exp_id = exp_id,
    quickform = quickform,
    n_NInd = temp_data$n_NInd[1],
    g_avg = temp_data$g_avg[1],
    g_sd = temp_data$g_sd[1]
  ) -> summary_data
  
  
}

#separate the experiment id
summary_data %>% separate(
  col = exp_id,
  into = c("polymer","FCED", "conc","MW", "set","replicate", "side"),
  sep = "-",
  convert = TRUE, 
  remove = FALSE
) -> summary_data


#save to statistically compare and further graph in GraphPad Prism
write.csv(summary_data, "NInd-PVA-S2-Per-Experiment-Results.csv")

# Mutate to get fractional sd, then remove data where sd > 30% of values
summary_data %>% mutate(
  frac_sd = g_sd / g_avg
) -> summary_data

summary_data %>% filter(
  frac_sd < 0.3
) -> summary_data

#removed so the pivot wider would work well
summary_data %>% select(-exp_id) -> summary_data

#The next few operations convert "Side A vs B" to "Hi vs Lo Side"
summary_data %>% pivot_wider(
  names_from = side,
  names_glue = "{side}_{.value}",
  values_from = c(n_NInd, g_avg, g_sd, frac_sd)
) ->wide_data

wide_data %>% mutate(
  Hi_avg = pmax(A_g_avg, B_g_avg, na.rm = TRUE),
  Lo_avg = pmin(A_g_avg, B_g_avg, na.rm = TRUE)
) -> wide_data

wide_data %>% mutate(
  Hi_sd = ifelse(is.na(A_g_avg), B_g_sd, ifelse(Hi_avg == A_g_avg, A_g_sd, B_g_sd)),
  Lo_sd = ifelse(is.na(A_g_avg), B_g_sd, ifelse(Lo_avg == A_g_avg, A_g_sd, B_g_sd)),
  Hi_n = ifelse(is.na(A_g_avg), B_n_NInd, ifelse(Hi_avg == A_g_avg, A_n_NInd, B_n_NInd)),
  Lo_n = ifelse(Hi_avg == Lo_avg, NA, ifelse(Hi_avg == A_g_avg, A_n_NInd, B_n_NInd))
) -> wide_data

wide_data %>% select(
  -starts_with("A_"),
  -starts_with("B_")
) -> wide_data

#Here we actually calculate differences and averages over the sides, propagating errors.
wide_data %>% mutate(
  diff_avg = ifelse(Hi_avg - Lo_avg == 0, NA, Hi_avg - Lo_avg),
  diff_err = ifelse(Hi_avg - Lo_avg == 0, NA, sqrt(Hi_sd ^ 2 + Lo_sd ^ 2)),
  tot_avg = (Hi_avg + Lo_avg)/2,
  tot_n = ifelse(is.na(Lo_n), 9, 18),
  tot_err = ifelse(is.na(Lo_n), Hi_sd, tot_avg * sqrt((Hi_sd/Hi_avg)^2 + (Lo_sd/Lo_avg)^2))
) -> wide_data

#removing now-obsolete side data
wide_data %>% select(
  -starts_with("Hi_"),
  -starts_with("Lo_")
) -> wide_data

#export a spreadsheet here for creating a figure of sample side differences.
write.csv(wide_data, "NInd-PVA-S2-Side-Difference-Results.csv")
#for the figure, make tot_avg stiffness the x and difference the y with std dev error bars. 
#See if there is any trends based on stiffness or formulation

#finally to summarize NInd-PVA-S2 data, grouping side differences by formulation
wide_data %>% group_by(polymer, FCED, conc, MW, quickform) %>% summarize(
  g_avg = mean(tot_avg),
  g_sd =  sd(tot_avg),
  sidediff = mean(diff_avg, na.rm = TRUE),
  side_sd = sd(diff_avg, na.rm = TRUE),
  n_NInd = n()
) -> final_data

#publish the data
write.csv(final_data, "NInd-PVA-S2-Summary.csv")