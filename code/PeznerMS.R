# - - - - - - - - - - - - - - - - - - - - - - - -
# Code for Pezner et al. MS - Coral reef hypoxia analysis
# Last updated, January 10, 2023
# Written by: Ariel Pezner (apezner@ucsd.edu)
# - - - - - - - - - - - - - - - - - - - - - - - - 

### Load packages and data, etc. ----

## Set working directory (will vary by computer)
setwd("~/Documents/SIO/O2 Variability")

## Load necessary packages
library(ggplot2);library(scales);library(ggridges);library(dplyr);library(ggthemes);library(patchwork);library(gsw);library(tidyverse);
library(readr);library(DescTools);library(ggmap);library(pastecs);library(forcats);library(reshape2);library(lubridate);library(paletteer)

## Load color palettes
library(RColorBrewer); library(wesanderson)

## Load functions to calculate DO saturation concentration (2 ways)
source('Code/R/calcDOatsat.R') 

## Load necessary data
load("RDataFiles/DO_all_Pezner.Rda") # all DO data for all sites
cmip6_sst = read_csv("Data/cmip6_sst_increase.csv") # CMIP6 predictions of mean T change by location by 2100

### Calculating the Impact of Warming on Solubility ----

## Add CMIP temperature projection data to new DO_all dataframe called "DO_all_proj", matching by location
DO_all_proj=merge(DO_all,cmip6_sst,by="loc")

## Add a +6ºC heatwave temperature scenario
DO_all_proj = DO_all_proj %>%
  mutate(heatwave = 6)

## Calculate new solubility if temperature increased by different projection amounts (measured T + CMIP projection)
DO_all_proj = DO_all_proj %>%
  mutate(DOsol_ssp126 = calcDOatsat_GG(temp + ssp126, sal),
         DOsol_ssp245 = calcDOatsat_GG(temp + ssp245, sal),
         DOsol_ssp370 = calcDOatsat_GG(temp + ssp370, sal),
         DOsol_ssp585 = calcDOatsat_GG(temp + ssp585, sal),
         DOsol_heatwave = calcDOatsat_GG(temp + heatwave, sal))

## Calculate ∆DOsol (change in DO due to changes in solubility)
DO_all_proj = DO_all_proj %>%
  mutate(deltaDOsol_ssp126 = oxyatsat-DOsol_ssp126,
         deltaDOsol_ssp245 = oxyatsat-DOsol_ssp245,
         deltaDOsol_ssp370 = oxyatsat-DOsol_ssp370,
         deltaDOsol_ssp585 = oxyatsat-DOsol_ssp585,
         deltaDOsol_heatwave = oxyatsat-DOsol_heatwave)


### Calculate some daily statistics needed later on ----      

## Calculate daily DO minimum 
DO_all_proj_daily = DO_all_proj %>%
  mutate(day = day(datetime), month = month(datetime)) %>%
  group_by(site,month,day) %>%
  summarize(dailymin_oxy = min(oxy, na.rm = TRUE)) # calculate daily min of present day DO (will need this later)

## Remove first and last days because we only want to use full days (not partials)
DO_all_proj_daily_sub = DO_all_proj_daily %>% 
  group_by(site) %>% 
  slice(2:n()) # removes first day (first row for that site)

DO_all_proj_daily_sub = DO_all_proj_daily_sub %>% 
  group_by(site) %>%  
  slice(1:(n()-1)) # removes last day (last row for that site)

## Now average these values across all days for each site for entire deployment (without first and last day)
daily_means = DO_all_proj_daily_sub %>%
  group_by(site) %>%
  summarize(mean_dailymin_oxy = mean(dailymin_oxy,na.rm = TRUE))

## Also,calculate mean DO across entire data set for each site:
mean_oxy = DO_all%>% 
            group_by(site) %>% 
            summarize(mean_oxy = mean(oxy, na.rm = TRUE))

## Merge dfs:
  DO_off_calc=merge(mean_oxy, daily_means, by="site")

### Calculating the Impact of Warming on Respiration ----

## Using Q10 equation, assuming Q10 = 2 (respiration rate doubles with 10ºC increase in T):
#
#    R2 = R1 * Q10^[(T2-T1)/10]
#    which is equivalent to:
#    R2 = R1 * Q10^(Trise/10)
#
# also assume R1 = 200 (this value doesn't matter as long as use the same for all calcs since only ratio of R2/R1 matters in the end),
# T1 is present day (measured T) and T2 is projected increased T (T1 + Trise)

  # Make function to calculate R2
  calcR2 <- function(R1,Q10,Trise) {
    
    # Input:
    #   R1 = present day respiration rate
    #   Q10 = change in R for 10º increase in T (e.g. if Q10 = 2, then R doubles for every 10º increase)
    #   Trise = increase in temp from present day (T2-T1 = Trise)
    # 
    # Output:
    #  R2 = projected respiration rate given increase in T
    #
    
    # Calculate R2
    R2 = R1 * (Q10 ^(Trise/10))
    
  }
  
  # Add R1 column to dataframe
  DO_all_proj = DO_all_proj %>%
    mutate(R1 = 200)
  
  # Calculate R2s with different temperature increases
  DO_all_proj = DO_all_proj %>%
    mutate(R2_ssp126 = calcR2(R1,2,ssp126),
           R2_ssp245 = calcR2(R1,2,ssp245),
           R2_ssp370 = calcR2(R1,2,ssp370),
           R2_ssp585 = calcR2(R1,2,ssp585),
           R2_heatwave = calcR2(R1,2,heatwave))

## Calculate ∆DO (amount to shift DO concentration curve down as a result of change in R)

  # Basic principle
  #   
  #  DOoff = meanDO - meanDOmin         [Typical offset from mean DO and mean daily min (respiration signal)]
  #  DOoff_Q10 = R2/R1 * DOoff          [change in DO due to warming, which is proportional to changes in R as a result of warming]
  #  deltaDO_Q10 = DOoff_Q10 - DOoff    [shift applied to entire time series]
  
  # Calculate DOoff
  DO_off_calc = DO_off_calc %>%
    mutate(DOoff = mean_oxy - mean_dailymin_oxy)
  
  # Add DOoff to DO_all_proj
  DO_all_proj=merge(DO_all_proj, DO_off_calc, by="site")
  
  # Calculate DOoff_Q10
  DO_all_proj = DO_all_proj %>%
    mutate(DOoff_Q10_ssp126 = (R2_ssp126/R1)*DOoff,
           DOoff_Q10_ssp245 = (R2_ssp245/R1)*DOoff,
           DOoff_Q10_ssp370 = (R2_ssp370/R1)*DOoff,
           DOoff_Q10_ssp585 = (R2_ssp585/R1)*DOoff,
           DOoff_Q10_heatwave = (R2_heatwave/R1)*DOoff)
  
  # Calculate ∆DO_Q10  
  DO_all_proj = DO_all_proj %>%
    mutate(deltaDO_Q10_ssp126 = DOoff_Q10_ssp126-DOoff,
           deltaDO_Q10_ssp245 = DOoff_Q10_ssp245-DOoff,
           deltaDO_Q10_ssp370 = DOoff_Q10_ssp370-DOoff,
           deltaDO_Q10_ssp585 = DOoff_Q10_ssp585-DOoff,
           deltaDO_Q10_heatwave = DOoff_Q10_heatwave-DOoff)

### Combining the Impacts of Warming on Solubility and Respiration ----

## Calculate ∆DO_T_Q10
DO_all_proj = DO_all_proj %>%
  mutate(deltaDO_T_Q10_ssp126 = deltaDO_Q10_ssp126 + deltaDOsol_ssp126,
         deltaDO_T_Q10_ssp245 = deltaDO_Q10_ssp245 + deltaDOsol_ssp245,
         deltaDO_T_Q10_ssp370 = deltaDO_Q10_ssp370 + deltaDOsol_ssp370,
         deltaDO_T_Q10_ssp585 = deltaDO_Q10_ssp585 + deltaDOsol_ssp585,
         deltaDO_T_Q10_heatwave = deltaDO_Q10_heatwave + deltaDOsol_heatwave)

## Calculate new projected DO
DO_all_proj = DO_all_proj %>%
  mutate(projDO_ssp126 = oxy - deltaDO_T_Q10_ssp126,
         projDO_ssp245 = oxy - deltaDO_T_Q10_ssp245,
         projDO_ssp370 = oxy - deltaDO_T_Q10_ssp370,
         projDO_ssp585 = oxy - deltaDO_T_Q10_ssp585,
         projDO_heatwave = oxy - deltaDO_T_Q10_heatwave)

## Replace negative DO values with zeroes (assume 0 is lowest possible concentration)
DO_all_proj$projDO_ssp126[DO_all_proj$projDO_ssp126 < 0] <- 0
DO_all_proj$projDO_ssp245[DO_all_proj$projDO_ssp245 < 0] <- 0
DO_all_proj$projDO_ssp370[DO_all_proj$projDO_ssp370 < 0] <- 0
DO_all_proj$projDO_ssp585[DO_all_proj$projDO_ssp585 < 0] <- 0
DO_all_proj$projDO_heatwave[DO_all_proj$projDO_heatwave < 0] <- 0


# Calculate n (number of observations)
n_vals_all = DO_all_proj %>% group_by(site) %>% tally()

### Create Figure 2A-----

## Rearrange data and subset for ease of plotting:
  
  # Create new df with just a subset of the bigger df (just need projected DO and other identifiers)
  DO_all_proj_subset = DO_all_proj %>% 
    select(site,loc,ocean,ocean_num,duration,recfreq,depth,datetime,oxy,projDO_ssp126,projDO_ssp245,projDO_ssp370,
           projDO_ssp585,projDO_heatwave)
  
  # Rename columns for alphabetical order when plotting later
  DO_all_proj_subset = DO_all_proj_subset %>% 
    rename(present = oxy,
           ssp126 = projDO_ssp126,
           ssp245 = projDO_ssp245,
           ssp370 = projDO_ssp370,
           ssp585 = projDO_ssp585,
           heatwave = projDO_heatwave) 
  
  # Combine all DO projections + present day into one oxy_proj column
  DO_all_proj_melt = melt(DO_all_proj_subset,id.vars = c("site","loc","ocean","ocean_num","duration","recfreq","depth","datetime"), 
                          measure.vars = c("heatwave","ssp585","ssp370","ssp245","ssp126","present"),
                          variable.name = "projection",value.name = "oxy_proj")

## Calculate daily stats
  # Calculate daily min/max oxygen: present day, under each temp scenario in a melted (vertical) format that will be easier to plot
  DO_all_proj_melt_daily = DO_all_proj_melt %>% 
    mutate(day = day(datetime),month = month(datetime)) %>%
    group_by(projection,site,month,day, ocean_num) %>%
    summarize(dailymin_oxy = min(oxy_proj,na.rm = TRUE),
              dailymax_oxy = max(oxy_proj,na.rm = TRUE))
  
  # Remove first and last days because we only want to use full days (not partials)
  DO_all_proj_melt_daily = DO_all_proj_melt_daily %>% 
    group_by(projection,site) %>% 
    slice(2:n()) # removes first day (first row for that site)
  
  DO_all_proj_melt_daily = DO_all_proj_melt_daily %>% 
    group_by(projection,site) %>% 
    slice(1:(n()-1)) # removes last day (last row for that site)
  
  # Now average these values across all days for each site for entire deployment (minus first and last day still)
  DO_melt_daily = DO_all_proj_melt_daily %>% 
    group_by(projection,site,ocean_num) %>% 
    summarize(mean_min = mean(dailymin_oxy,na.rm = TRUE),
              mean_max = mean(dailymax_oxy,na.rm = TRUE),
              min_min = min(dailymin_oxy,na.rm = TRUE),
              max_max = max(dailymax_oxy,na.rm = TRUE))
  
## Make a bar plot with mean and extreme daily DO ranges under each scenario for each site [FIGURE 2a]

  # Rearrange order of projections to plot from present to heatwave left to right, then create line-dot plot
  pf1 = DO_melt_daily %>%
    mutate(projection = fct_relevel(projection,"present", "ssp126", "ssp245", "ssp370","ssp585","heatwave")) %>%
    ggplot() +
    
    # Add horizontal lines to denote thresholds
    geom_hline(yintercept = 153, linetype = "dashed", color = "grey", alpha = 1)+
    geom_hline(yintercept = 122, linetype = "dashed", color = "grey", alpha = 1)+
    geom_hline(yintercept = 92, linetype = "dashed", color = "grey", alpha = 1)+
    geom_hline(yintercept = 61, linetype = "dashed", color = "grey", alpha = 1)+
    
    # Plot mean min and mean max
    geom_linerange(aes(x=fct_reorder(site,ocean_num), ymin = mean_min, ymax = mean_max, color=projection, group = projection),
                   position=position_dodge(width=0.8), size = 1.4)+
    
    #Plot max max and min min
    geom_linerange(aes(x=fct_reorder(site,ocean_num), ymin = min_min, ymax = max_max, color=projection, group = projection),
                   position=position_dodge(width=0.8), size = 1.4, alpha = 0.3)+
    
    # Fix various details
    scale_y_continuous(expand = c(0, 0))+
    scale_color_brewer(palette = "RdYlBu", direction = -1, labels=c("Present Day","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5","Heatwave"))+
    labs(x="", y= expression(paste("Dissolved oxygen (µmol ",O[2], " ",kg^-1,")")), color = "Projection") +
    theme_bw() + 
    theme(legend.position=c(.7,0.9),legend.direction = "horizontal")+
    theme(axis.text.x = element_text(angle = 90))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))


### Sub-sampling all data to 30 min recording frequency------
  
### First we need to sub-sample the data for all sites to ensure all have a 30 min sampling frequency
  # There are 6 different sampling frequencies: 1, 0.5, 0.25, 0.167, 0.083, and 0.003 hr
  # so datasets can be grouped according to this somewhat. We do not need to do this for
  # datasets that are already recording at 0.5 hr-1 freq.
  
  ## Figure out which sites have which recording frequency:
  unique(DO_all_proj$site[DO_all_proj$recfreq == 1]) # Crocker 1a, 1b, 1c; Palmyra 1b
  unique(DO_all_proj$site[DO_all_proj$recfreq == 0.5]) # Bocas 1a,1b,1c,1d,2a,2b; Heron 1,2,3; Kaneohe 1a,1b,2; Okinawa 1; Palmyra 1a,2
  unique(DO_all_proj$site[DO_all_proj$recfreq == 0.25]) # Dongsha 1, Dongsha 2, Hog 1a, Hog 1b
  unique(DO_all_proj$site[DO_all_proj$recfreq == 0.167]) # Dongsha 3, Okinawa 2, Okinawa 3, Taiping 1, Taiping 2
  unique(DO_all_proj$site[DO_all_proj$recfreq == 0.083]) # Tutuila
  unique(DO_all_proj$site[DO_all_proj$recfreq == 0.003]) # Baker, Jarvis, Palmyra 3
  
  ## Create a list to store all of this data in for tidyness
  DO_all_proj_30min = list()
  
  ## Sub-sample: rec freq = 0.003 (every 180 rows) for Baker, Jarvis, Palmyra 3
  DO_all_proj_30min[[1]] = DO_all_proj_melt %>%
    filter(site == "Baker" | site == "Jarvis" | site == "Palmyra 3") %>% 
    group_by(site,projection) %>% 
    arrange(datetime) %>%
    slice(which(row_number() %% 180 == 1))
  
  ## Sub-sample: rec freq = 0.083 (every 6 rows) for Tutuila
  DO_all_proj_30min[[2]] = DO_all_proj_melt %>%
    filter(site == "Tutuila") %>% 
    group_by(projection) %>% 
    arrange(datetime) %>%
    slice(which(row_number() %% 6 == 1))
  
  ## Sub-sample: rec freq = 0.167 (every 3 rows) for Dongsha 3, Okinawa 2, Okinawa 3, Taiping 1, Taiping 2
  DO_all_proj_30min[[3]] = DO_all_proj_melt %>%
    filter(site == "Dongsha 3" |site == "Okinawa 2" | site == "Okinawa 3" | site == "Taiping 1" | 
             site == "Taiping 2") %>% 
    group_by(site, projection) %>% 
    arrange(datetime) %>%
    slice(which(row_number() %% 3 == 1))  
  
  ## Sub-sample: rec freq = 0.25 (every 2 rows) for Dongsha 1, Dongsha 2, Hog 1a, Hog 1b
  DO_all_proj_30min[[4]] = DO_all_proj_melt %>%
    filter(site == "Dongsha 1" |site == "Dongsha 2" | site == "Hog 1a" | site == "Hog 1b") %>% 
    group_by(site, projection) %>% 
    arrange(datetime) %>% 
    slice(which(row_number() %% 2 == 1))
  
  ## Select site data for those with 0.5 hr-1 rec freq: Bocas 1a,1b,1c,1d,2a,2b; Heron 1,2,3; Kaneohe 1a,1b,2; Okinawa 1; Palmyra 1a,2
  DO_all_proj_30min[[5]] = DO_all_proj_melt %>%
    filter(loc == "Bocas del Toro" | loc == "Kaneohe Bay" | loc == "Heron Island" | site == "Okinawa 1" | 
             site == "Palmyra 1a" | site == "Palmyra 2")
  
  ## Interpolate data with rec freq = 1 : Crocker 1a, 1b, 1c; Palmyra 1b
  
  ## Make a function to filter data, interpolate, make new df
  interp_dat <- function (sitename, proj) {
    # Create temporary df from DO_all_proj_melt of desired site and projection, arrange by datetime so it's in order
    temp = DO_all_proj_melt %>% filter(site == sitename & projection == proj) %>%  arrange(datetime)
    
    # Interpolate using approx, where xout is datetimes with 30 min intervals instead of 1 hr, get datetime and interpolated oxy_proj
    datetime = approx(temp$datetime,temp$oxy_proj,xout = seq(min(temp$datetime),max(temp$datetime),by='30 min'))$x
    oxy_proj = approx(temp$datetime,temp$oxy_proj,xout = seq(min(temp$datetime),max(temp$datetime),by='30 min'))$y
    
    # Create a new df with the interpolated datetime and oxy_proj and then add in all identifiers
    new = data.frame(datetime,oxy_proj)
    new = new %>% 
      mutate(site = rep(temp$site[1], length(new$oxy_proj)),
             loc = rep(temp$loc[1], length(new$oxy_proj)),
             ocean = rep(temp$ocean[1], length(new$oxy_proj)),
             ocean_num = rep(temp$ocean_num[1], length(new$oxy_proj)),
             duration = rep(temp$duration[1], length(new$oxy_proj)),
             recfreq = rep(temp$recfreq[1], length(new$oxy_proj)),
             depth = rep(temp$depth[1], length(new$oxy_proj)),
             projection = rep(temp$projection[1], length(new$oxy_proj)))
    new <- new[, c("site", "loc", "ocean","ocean_num","duration","recfreq","depth","datetime","projection","oxy_proj")] #reorder columns to match
    
    # print(new)
  }
  
  ## Perform this interpolation function on each dataset and for each projection scenario
  DO_all_proj_30min[[6]] = interp_dat("Crocker 1a","present")
  DO_all_proj_30min[[7]] = interp_dat("Crocker 1a","ssp126")
  DO_all_proj_30min[[8]] = interp_dat("Crocker 1a","ssp245")
  DO_all_proj_30min[[9]] = interp_dat("Crocker 1a","ssp370")
  DO_all_proj_30min[[10]] = interp_dat("Crocker 1a","ssp585")
  DO_all_proj_30min[[11]] = interp_dat("Crocker 1a","heatwave")
  
  DO_all_proj_30min[[12]] = interp_dat("Crocker 1b","present")
  DO_all_proj_30min[[13]] = interp_dat("Crocker 1b","ssp126")
  DO_all_proj_30min[[14]] = interp_dat("Crocker 1b","ssp245")
  DO_all_proj_30min[[15]] = interp_dat("Crocker 1b","ssp370")
  DO_all_proj_30min[[16]] = interp_dat("Crocker 1b","ssp585")
  DO_all_proj_30min[[17]] = interp_dat("Crocker 1b","heatwave")
  
  DO_all_proj_30min[[18]] = interp_dat("Crocker 1c","present")
  DO_all_proj_30min[[19]] = interp_dat("Crocker 1c","ssp126")
  DO_all_proj_30min[[20]] = interp_dat("Crocker 1c","ssp245")
  DO_all_proj_30min[[21]] = interp_dat("Crocker 1c","ssp370")
  DO_all_proj_30min[[22]] = interp_dat("Crocker 1c","ssp585")
  DO_all_proj_30min[[23]] = interp_dat("Crocker 1c","heatwave")
  
  DO_all_proj_30min[[24]] = interp_dat("Palmyra 1b","present")
  DO_all_proj_30min[[25]] = interp_dat("Palmyra 1b","ssp126")
  DO_all_proj_30min[[26]] = interp_dat("Palmyra 1b","ssp245")
  DO_all_proj_30min[[27]] = interp_dat("Palmyra 1b","ssp370")
  DO_all_proj_30min[[28]] = interp_dat("Palmyra 1b","ssp585")
  DO_all_proj_30min[[29]] = interp_dat("Palmyra 1b","heatwave")
  
  ## Now, combine all of the newly formatted data frames together into one
  DO_all_proj_30min_all = bind_rows(DO_all_proj_30min)
  
  ## Calculate n (number of observations)
  n_vals_all_30min = DO_all_proj_30min_all %>% filter(projection == "present") %>% group_by(site) %>% tally()
  
### Figure 2B-E and F: Intensity, Duration, and Severity Calculation [Figs. S5-7] -------   
  
### Theory of Intensity and Severity Calculations ---
  # After Hauri et al., 2013 (https://doi.org/10.1002/grl.50618): 
  #
  # Severity = Intensity * Duration
  #
  # where:
  # Intensity = threshold - mean during event
  #
  # so, 
  # Intensity = hyp_threshold - mean(DO)_hyp_event
  #
  # and
  # Duration = duration of hypoxic event of certain threshold
  
  ## Add ID columns for 4 thresholds (<153, <122, <92, and <61 µmol/kg) 
  DO_all_proj_30min_all = DO_all_proj_30min_all %>% 
    arrange(loc, site, projection, datetime) %>% 
    mutate(hyp_153_ID = ifelse(oxy_proj <=153, 1, 0), #same as hyp_all_ID
           hyp_122_ID = ifelse(oxy_proj <=122, 1, 0),
           hyp_92_ID = ifelse(oxy_proj <=92 , 1, 0),
           hyp_61_ID = ifelse(oxy_proj <=61, 1, 0)) #same as hyp_sev_ID 
  
### WORKING WITH INDIVIDUAL OBSERVATIONS  ---
  
  ## Where oxy_proj < threshold, calculate intensity below threshold, otherwise put NA:
  intens_hyp_thresholds = DO_all_proj_30min_all %>% 
    mutate(intens_153 = ifelse(oxy_proj <=153, 153 - oxy_proj, NA),
           intens_122 = ifelse(oxy_proj <=122, 122 - oxy_proj, NA),
           intens_92  = ifelse(oxy_proj <=92, 92 - oxy_proj, NA),
           intens_61  = ifelse(oxy_proj <=61, 61 - oxy_proj, NA))
  
  ## Calculate number of data points below each threshold for each projection:    
  num_intens = intens_hyp_thresholds %>% 
    group_by(projection) %>% 
    summarize("153" = length(which(intens_153 > 0)),
              "122" = length(which(intens_122 > 0)),
              "92" = length(which(intens_92 > 0)),
              "61" = length(which(intens_61 > 0)))
  
  num_intens_melt = melt(num_intens, variable.name = "threshold", value.name = "num_pts")
  
  ## And plot:
  # select colors
  mycols = brewer.pal(n = 9, "YlGnBu")[4:7]
  
  # Number of events by threshold vs projection
  figS6a = num_intens_melt %>% 
    mutate(projection = fct_relevel(projection,"present", "ssp126", "ssp245", "ssp370","ssp585","heatwave")) %>% 
    ggplot()+
    geom_line(aes(x = projection, y = num_pts, color = threshold, group = threshold), size = 1.5) +
    geom_point(aes(x = projection, y = num_pts, color = threshold, group = threshold), size = 1.5) +
    theme_classic() + labs(x = "", y = "Number of observations", color = expression(paste("Threshold (µmol ",O[2], " ",kg^-1,")"))) +
    scale_color_manual(values = mycols) + scale_x_discrete(labels=c("Present Day","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5","Heatwave"))+
    theme(legend.position = c(0.4,0.95),
          legend.direction = "horizontal",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))
  
  # Percent change in number of observations by threshold vs projection
  pchange_intens = num_intens
  
  pchange_intens$`153` = ((pchange_intens$`153`-pchange_intens$`153`[pchange_intens$projection == "present"])/
                            pchange_intens$`153`[pchange_intens$projection == "present"])*100
  pchange_intens$`122` = ((pchange_intens$`122`-pchange_intens$`122`[pchange_intens$projection == "present"])/
                            pchange_intens$`122`[pchange_intens$projection == "present"])*100
  pchange_intens$`92` = ((pchange_intens$`92`-pchange_intens$`92`[pchange_intens$projection == "present"])/
                           pchange_intens$`92`[pchange_intens$projection == "present"])*100
  pchange_intens$`61` = ((pchange_intens$`61`-pchange_intens$`61`[pchange_intens$projection == "present"])/
                           pchange_intens$`61`[pchange_intens$projection == "present"])*100
  # remove present row
  pchange_intens = subset(pchange_intens, projection!="present")
  # melt
  pchange_intens_melt = melt(pchange_intens, variable.name = "threshold", value.name = "pchange")
  
  # Custom colors
  mycols2 = rev(brewer.pal(n = 6, "RdYlBu")[1:5])

  # Figure 2F: Plot percent change in number of hypoxic events of 4 thresholds for each projection relative to present
  pf3 = pchange_intens_melt %>% 
    mutate(projection = fct_relevel(projection,"present", "ssp126", "ssp245", "ssp370","ssp585","heatwave")) %>% 
    ggplot()+
    geom_linerange(aes(x=threshold, ymin = 0, ymax = pchange, color=projection, group = projection),
                   position=position_dodge(width=0.85), size = 5.5)+
    theme_classic() + labs(color = "Projection", y = "% change in hypoxia", x = expression(paste("Threshold (µmol ",O[2], " ",kg^-1,")"))) +
    theme(legend.position = "top") +
    scale_color_manual(values = mycols2, labels=c("SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5","Heatwave"))+
    scale_fill_manual(values = mycols2, labels=c("SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5","Heatwave"))+
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black")) +
    scale_y_continuous(expand = c(0, 0)) # force y-axis to start at 0
  
  
### WORKING WITH EVENTS (groups of data points below threshold)  ---    
  
  ## Make data frame that groups by rle(hyp_all_ID), site, and proj and then takes average of oxy_proj for each group where ID == 1
  # This will be part of intensity calculation  
  # Code source: https://stackoverflow.com/questions/59991203/how-to-group-consecutive-rows-having-same-event-and-find-average
  
  library(data.table)
  
  event_hyp_153 = DO_all_proj_30min_all %>%
    arrange(loc, site, projection, datetime) %>%  #make sure datetime is in order before grouping events
    group_by(loc,site, projection, grp = rleid(hyp_153_ID),recfreq) %>% # group by RLE of hyp_153_ID as well
    summarise(ID = first(hyp_153_ID),
              mean_oxy = mean(oxy_proj), # mean of oxy_proj for run of 1's
              dur_hrs = as.numeric(difftime(last(datetime), first(datetime),units="hours"))) %>%  
    filter(ID == 1) # only want rows where hyp_ID == 1
  
  event_hyp_122 = DO_all_proj_30min_all %>%
    arrange(loc, site, projection, datetime) %>% 
    group_by(loc,site, projection, grp = rleid(hyp_122_ID),recfreq) %>% # group by RLE of hyp_122_ID as well
    summarise(ID = first(hyp_122_ID),
              mean_oxy = mean(oxy_proj), # mean of oxy_proj for run of 1's
              dur_hrs = as.numeric(difftime(last(datetime), first(datetime),units="hours"))) %>%  
    filter(ID == 1) # only want rows where hyp_ID == 1
  
  event_hyp_92 = DO_all_proj_30min_all %>%
    arrange(loc, site, projection, datetime) %>% 
    group_by(loc,site, projection, grp = rleid(hyp_92_ID),recfreq) %>% # group by RLE of hyp_92_ID as well
    summarise(ID = first(hyp_92_ID),
              mean_oxy = mean(oxy_proj), # mean of oxy_proj for run of 1's
              dur_hrs = as.numeric(difftime(last(datetime), first(datetime),units="hours"))) %>%  
    filter(ID == 1) # only want rows where hyp_ID == 1
  
  event_hyp_61 = DO_all_proj_30min_all %>%
    arrange(loc, site, projection, datetime) %>% 
    group_by(loc, site, projection, recfreq, grp = rleid(hyp_61_ID)) %>% # group by RLE of hyp_61_ID as well
    summarise(ID = first(hyp_61_ID),
              mean_oxy = mean(oxy_proj), # mean of oxy_proj for run of 1's
              dur_hrs = as.numeric(difftime(last(datetime), first(datetime),units="hours"))) %>%  
    filter(ID == 1) # only want rows where hyp_ID == 1
  
  ## For hypoxic events only lasting 1 data point (1 row), the duration will be 0. Change 0's to 1x recording frequency instead (0.5hr), as this is the
  # theoretical mean duration of that "event". For example, if recording freq is every 30 min, then one data point of hypoxia could 
  # represent 30 min of hypoxia prior to that sampling time or also 30 min of hypoxia until the next sampling event, so 30 min is mean
  
  event_hyp_153$dur_hrs <- with(event_hyp_153, replace(dur_hrs, dur_hrs == 0, 0.5))
  
  event_hyp_122$dur_hrs <- with(event_hyp_122, replace(dur_hrs, dur_hrs == 0, 0.5))
  
  event_hyp_92$dur_hrs <- with(event_hyp_92, replace(dur_hrs, dur_hrs == 0, 0.5))
  
  event_hyp_61$dur_hrs <- with(event_hyp_61, replace(dur_hrs, dur_hrs == 0, 0.5))
  
  ## Now, calculate intensity, then combine with duration information to calculate severity for each event for each site and projection     
  
  event_hyp_153 =  event_hyp_153 %>% 
    mutate(intensity = 153 - mean_oxy, # 5mg/L
           severity = intensity * dur_hrs,
           threshold = as.factor("153"))
  
  event_hyp_122 =  event_hyp_122 %>% 
    mutate(intensity = 122 - mean_oxy, # 4mg/L
           severity = intensity * dur_hrs,
           threshold = as.factor("122"))
  
  event_hyp_92 =  event_hyp_92 %>% 
    mutate(intensity = 92 - mean_oxy, # 3mg/L
           severity = intensity * dur_hrs,
           threshold = as.factor("92"))
  
  event_hyp_61 =  event_hyp_61 %>% 
    mutate(intensity = 61 - mean_oxy, # 2mg/L
           severity = intensity * dur_hrs,
           threshold = as.factor("61"))     
  
  
  events_hyp_all = rbind(event_hyp_153,event_hyp_122,event_hyp_92,event_hyp_61)

  
  # Violin plot of all DURATION by proj and thresh
  figS7a = events_hyp_all %>% 
    mutate(projection = fct_relevel(projection,"present", "ssp126", "ssp245", "ssp370","ssp585","heatwave"),
           threshold = fct_relevel(threshold,"61", "92", "122", "153")) %>%
    ggplot()+
    geom_violin(aes(x = threshold, y = dur_hrs, color = projection),alpha = 0.3)+
    geom_point(aes(x = threshold, y = dur_hrs, color = projection), size = 0.4, alpha = 0.8,position = position_dodge(width=0.9))+
    theme_classic() + labs(x = "", y = "Duration (hr)", color = "Projection")  + 
    theme(legend.position = "top")+ 
    scale_color_brewer(palette = "RdYlBu", direction = -1, labels=c("Present Day","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5","Heatwave"))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          legend.position=c(.5,0.9),
          legend.direction = "horizontal")
    
  
  # Violin plot of all INTENSITY by proj and thresh
  figS7b = events_hyp_all %>% 
    mutate(projection = fct_relevel(projection,"present", "ssp126", "ssp245", "ssp370","ssp585","heatwave"),
           threshold = fct_relevel(threshold,"61", "92", "122", "153")) %>%
    ggplot()+
    geom_violin(aes(x = threshold, y = intensity, color = projection),alpha = 0.3)+
    geom_point(aes(x = threshold, y = intensity, color = projection), size = 0.4, alpha = 0.8,position = position_dodge(width=0.9))+
    theme_classic() + 
    labs(x = "", y = expression(paste("Intensity (µmol ",O[2], " ",kg^-1,")")), color = "Projection")  + 
    scale_color_brewer(palette = "RdYlBu", direction = -1, labels=c("Present Day","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5","Heatwave"))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          legend.position="none")
  
  # Violin plot of all SEVERITY by proj and thresh
  figS7c = events_hyp_all %>% 
    mutate(projection = fct_relevel(projection,"present", "ssp126", "ssp245", "ssp370","ssp585","heatwave"),
           threshold = fct_relevel(threshold,"61", "92", "122", "153")) %>%
    ggplot()+
    geom_violin(aes(x = threshold, y = severity, color = projection),alpha = 0.3)+
    geom_point(aes(x = threshold, y = severity, color = projection), size = 0.4, alpha = 0.8,position = position_dodge(width=0.9))+
    theme_classic() + 
    labs(x = expression(paste("Threshold (µmol ",O[2], " ",kg^-1,")")), y = expression(paste("Severity DO (µmol ",O[2], " ",kg^-1,"hr)")), color = "Projection")  + 
    scale_color_brewer(palette = "RdYlBu", direction = -1, labels=c("Present Day","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5","Heatwave"))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          legend.position="none")
  
  figS7 = figS7a/figS7b/figS7c
  
  ## Calculate cumulative duration, mean duration, cumulative intesity, mean intensity, mean severity, cumulative severity
  events_hyp_summary = events_hyp_all %>% 
    group_by(projection, threshold) %>% 
    summarize(tot_dur = sum(dur_hrs),
              mean_dur = mean(dur_hrs),
              sd_dur = sd(dur_hrs),
              
              tot_intens = sum(intensity),
              mean_intens = mean(intensity),
              sd_intens = sd(intensity),
              
              tot_sev = sum(severity),
              mean_sev = mean(severity),
              sd_sev = sd(severity))
  # Plots:
  
  # Cumulative duration
  figS8a = events_hyp_summary %>% 
    mutate(projection = fct_relevel(projection,"present", "ssp126", "ssp245", "ssp370","ssp585","heatwave")) %>%
    ggplot()+
    geom_line(aes(x = projection, y = tot_dur, color = threshold, group = threshold), size = 1.5)+
    geom_point(aes(x = projection, y = tot_dur, color = threshold, group = threshold), size = 2)+
    theme_classic() + scale_color_manual(values = mycols) + theme(legend.position = "top")+
    scale_x_discrete(labels=c("Present Day","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5","Heatwave"))+
    labs(x = "", y = "Cumulative duration (hr)", color = expression(paste("Threshold (µmol ",O[2], " ",kg^-1,")"))) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          legend.position=c(.4,0.9),
          legend.direction = "horizontal")
  
  # Cumulative intensity
  figS8b = events_hyp_summary %>% 
    mutate(projection = fct_relevel(projection,"present", "ssp126", "ssp245", "ssp370","ssp585","heatwave")) %>%
    ggplot()+
    geom_line(aes(x = projection, y = tot_intens, color = threshold, group = threshold), size = 1.5)+
    geom_point(aes(x = projection, y = tot_intens, color = threshold, group = threshold), size = 2)+
    theme_classic() + scale_color_manual(values = mycols) + theme(legend.position = "top")+
    scale_x_discrete(labels=c("Present Day","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5","Heatwave"))+
    labs(x = "", y = expression(paste("Cumulative intensity (µmol ",O[2], " ",kg^-1,")")), color = expression(paste("Threshold (µmol ",O[2], " ",kg^-1,")")))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          legend.position="none")
  
  # Cumulative severity
  figS8c = events_hyp_summary %>% 
    mutate(projection = fct_relevel(projection,"present", "ssp126", "ssp245", "ssp370","ssp585","heatwave")) %>%
    ggplot()+
    geom_line(aes(x = projection, y = tot_sev, color = threshold, group = threshold), size = 1.5)+
    geom_point(aes(x = projection, y = tot_sev, color = threshold, group = threshold), size = 2)+
    scale_x_discrete(labels=c("Present Day","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5","Heatwave"))+
    theme_classic() + scale_color_manual(values = mycols) + theme(legend.position = "top")+
    labs(x = "Projection", y = expression(paste("Cumulative severity (µmol ",O[2], " ",kg^-1,"hr)")), color = expression(paste("Threshold (µmol ",O[2], " ",kg^-1,")")))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          legend.position="none")
  
  figS8 = figS8a/figS8b/figS8c
  
  
### Calculate total number of events for each projection and threshold
  events_hyp_all_tot = events_hyp_all %>% 
    group_by(projection, threshold) %>% 
    summarize(totevents = length(intensity))  
  
  ## Plot - total number of events by proj and threshold
  figS6b = events_hyp_all_tot %>% 
    mutate(projection = fct_relevel(projection,"present", "ssp126", "ssp245", "ssp370","ssp585","heatwave")) %>%
    ggplot()+
    geom_line(aes(x = projection, y = totevents, color = threshold, group = threshold), size = 1.5)+
    geom_point(aes(x = projection, y = totevents, color = threshold, group = threshold), size = 2)+
    scale_x_discrete(labels=c("Present Day","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5","Heatwave"))+
    theme_classic() + scale_color_manual(values = mycols) + theme(legend.position = "top")+
    labs(x = "Projection", y = "Number of events", color = "")+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          legend.position="none")
  
  fig6 = figS6a/figS6b
  
  ## Create Figure 2b: Isolate events of certain duration for each threshold (also used to make Table S6)
  events_hyp_all_dur = events_hyp_all %>% 
    group_by(projection, threshold, site) %>% 
    summarize(tot_less1hr = length(which(dur_hrs <=1)),
              tot_1to6hr = length(which(dur_hrs <= 6 & dur_hrs > 1)),
              tot_6to12hr = length(which(dur_hrs <= 12 & dur_hrs > 6)),
              tot_12to24hr = length(which(dur_hrs <= 24 & dur_hrs > 12)),
              tot_more24hr = length(which(dur_hrs > 24)))
  
  events_hyp_all_dur_melt = melt(events_hyp_all_dur, id.vars = c("site", "projection","threshold"), 
                                 measure.vars = c("tot_less1hr","tot_1to6hr","tot_6to12hr","tot_12to24hr","tot_more24hr"),
                                 variable.name = "duration_cat", value.name = "tot_events")
  
  events_hyp_all_dur_sites = events_hyp_all_dur_melt %>% 
                              group_by(projection, threshold, duration_cat) %>% 
                              summarize(num_sites = length(which(tot_events>0)), 
                                        num_events = sum(tot_events))
  
  pf2 = events_hyp_all_dur_sites %>% 
    mutate(projection = fct_relevel(projection,"present", "ssp126", "ssp245", "ssp370","ssp585","heatwave")) %>% 
    ggplot() +
    # Bars
    geom_bar(aes(x = duration_cat, weight = num_sites/32*100, fill = projection), position = "dodge", alpha = 1)+
    # Facet
    facet_wrap(~threshold, scales = "free")+
    # Other plot controls
    theme_classic()+ ylim(0,88)+ labs(x = "Duration of hypoxia (hours)", y = "% of sites")+
    # Custom color palette (ColorBrewer)
    scale_fill_manual(values = c("#4575B4","#91BFDB","#E0F3F8", "#FEE090","#FC8D59","#D73027"), 
                      labels=c("Present Day","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5","Heatwave"))+
    theme(legend.position="none")+ scale_x_discrete(labels=c("< 1 hr", "1 to 6 hr", "6 to 12 hr", "12 to 24 hr", "> 24hr"))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))



### Make combined Figure 2a,b,c ------
  pf4 = (pf2 | (pf3/ plot_spacer()))  + plot_layout(heights = c(2, 1), widths = c(2,1.5))
  
  fig2 = pf1/pf4
  
  
### Make table with min, max, mean daily DO, T, DOsat (Table S2) ----    
    # Calculate by site
    table_sites = DO_all%>% 
      mutate(day = day(datetime),month = month(datetime)) %>%
      group_by(site,month,day) %>%
      summarize(daily_mean_DO = mean(oxy, na.rm = TRUE),
                daily_range_DO = max(oxy, na.rm = TRUE)-min(oxy, na.rm = TRUE),
                daily_max_DO = max(oxy, na.rm = TRUE),
                daily_min_DO = min(oxy, na.rm = TRUE),
                
                daily_mean_DOsat = mean(oxysat, na.rm = TRUE),
                daily_range_DOsat = max(oxysat, na.rm = TRUE)-min(oxysat, na.rm = TRUE),
                daily_max_DOsat = max(oxysat, na.rm = TRUE),
                daily_min_DOsat = min(oxysat, na.rm = TRUE),
                
                daily_mean_T = mean(temp, na.rm = TRUE),
                daily_range_T = max(temp, na.rm = TRUE)-min(temp, na.rm = TRUE),
                daily_max_T = max(temp, na.rm = TRUE),
                daily_min_T = min(temp, na.rm = TRUE))
    
    # Remove first and last days because we only want to use full days (not partials)
    table_sites_sub = table_sites %>% 
      group_by(site) %>% 
      slice(2:n()) # removes first day (first row for that site)
    
    table_sites_sub = table_sites_sub %>% 
      group_by(site) %>%  
      slice(1:(n()-1)) # removes last day (last row for that site)
    
    # Average across remaining deployment days for each site (Table S2)
    table_sites_means =  table_sites_sub %>% 
      group_by(site) %>% 
      summarize(mean_dailymeanDO = mean(daily_mean_DO, na.rm = TRUE),
                sd_dailymeanDO = sd(daily_mean_DO, na.rm = TRUE),
                mean_dailyrangeDO = mean(daily_range_DO, na.rm = TRUE),
                sd_dailyrangeDO = sd(daily_range_DO, na.rm = TRUE),
                mean_dailymaxDO = mean(daily_max_DO, na.rm = TRUE),
                sd_dailymaxDO = sd(daily_max_DO, na.rm = TRUE),
                mean_dailyminDO = mean(daily_min_DO, na.rm = TRUE),
                sd_dailyminDO = sd(daily_min_DO, na.rm = TRUE),
                
                mean_dailymeanDOsat = mean(daily_mean_DOsat, na.rm = TRUE),
                sd_dailymeanDOsat = sd(daily_mean_DOsat, na.rm = TRUE),
                mean_dailyrangeDOsat = mean(daily_range_DOsat, na.rm = TRUE),
                sd_dailyrangeDOsat = sd(daily_range_DOsat, na.rm = TRUE),
                mean_dailymaxDOsat = mean(daily_max_DOsat, na.rm = TRUE),
                sd_dailymaxDOsat = sd(daily_max_DOsat, na.rm = TRUE),
                mean_dailyminDOsat = mean(daily_min_DOsat, na.rm = TRUE),
                sd_dailyminDOsat = sd(daily_min_DOsat, na.rm = TRUE),
                
                mean_dailymeanT = mean(daily_mean_T, na.rm = TRUE),
                sd_dailymeanT = sd(daily_mean_T, na.rm = TRUE),
                mean_dailyrangeT = mean(daily_range_T, na.rm = TRUE),
                sd_dailyrangeT = sd(daily_range_T, na.rm = TRUE),
                mean_dailymaxT = mean(daily_max_T, na.rm = TRUE),
                sd_dailymaxT = sd(daily_max_T, na.rm = TRUE),
                mean_dailyminT = mean(daily_min_T, na.rm = TRUE),
                sd_dailyminT = sd(daily_min_T, na.rm = TRUE))
    
    # Calculate by location (lumping all sites together for given location)
    table_loc = DO_all%>% 
      mutate(day = day(datetime),month = month(datetime)) %>%
      group_by(loc,month,day) %>%
      summarize(daily_mean = mean(oxy, na.rm = TRUE),
                daily_range = max(oxy, na.rm = TRUE)-min(oxy, na.rm = TRUE),
                daily_max = max(oxy, na.rm = TRUE),
                daily_min = min(oxy, na.rm = TRUE)) %>% 
      group_by(loc) %>% 
      summarize(mean_dailymean = mean(daily_mean, na.rm = TRUE),
                sd_dailymean = sd(daily_mean, na.rm = TRUE),
                mean_dailyrange = mean(daily_range, na.rm = TRUE),
                sd_dailyrange = sd(daily_range, na.rm = TRUE),
                mean_dailymax = mean(daily_max, na.rm = TRUE),
                sd_dailymax = sd(daily_max, na.rm = TRUE),
                mean_dailymin = mean(daily_min, na.rm = TRUE),
                sd_dailymin = sd(daily_min, na.rm = TRUE))
    
    # Calculate n (number of observations) used to calculate mean of daily stats (e.g. number of days)
    n_vals = table_sites_sub %>% group_by(site) %>% tally()
    
   
### Make Figure 3 - Timing of Hypoxia ----
    
    # Use only present-day data from 30 min dataframe:
    DO_all_30min_pres = DO_all_proj_30min_all %>% filter(projection == "present")
    
    # Add hour of day column
    DO_all_30min_pres$hour = hour(DO_all_30min_pres$datetime)
    
    # Find number of data points of all levels of hypoxia by hour
    timinghyp = DO_all_30min_pres %>% 
                  group_by(hour) %>% 
                  summarize(h153 = length(which(hyp_153_ID>0)),
                            h122 = length(which(hyp_122_ID>0)),
                            h92 = length(which(hyp_92_ID>0)),
                            h61 = length(which(hyp_61_ID>0)))
    
    # Calculate total number of observations for each threshold
    timinghyp_tot = timinghyp %>%  
                    summarize(tot153 = sum(h153),
                              tot122 = sum(h122),
                              tot92 = sum(h92),
                              tot61 = sum(h61))
    
    # Convert to % of observations for given threshold
    timinghyp = timinghyp %>% 
                mutate(h153 = (h153/timinghyp_tot$tot153)*100,
                       h122 = (h122/timinghyp_tot$tot122)*100,
                       h92 = (h92/timinghyp_tot$tot92)*100,
                       h61 = (h61/timinghyp_tot$tot61)*100)
    
    # Rename columns
    timinghyp = timinghyp %>% 
      mutate('153' = h153,
             '122' = h122,
             '92' = h92,
             '61' = h61)
    
    # Melt df for plotting
    timinghyp_melt = melt(timinghyp,id.vars = c("hour"), 
                          measure.vars = c("153","122","92","61"),
                          variable.name = "threshold",value.name = "pobs")
    
    # Plot
    fig3 = timinghyp_melt %>% 
    ggplot()+ 
      annotate("rect", xmin=6, xmax=18, ymin=-Inf, ymax=Inf, alpha=0.2, fill = "#eda63a")+ # add daytime box
      geom_point(aes(x = hour, y = pobs, color = threshold)) +
      scale_x_continuous("Hour of day", breaks = c(0:23))+ 
      theme(panel.grid.minor = element_blank())+ theme_classic()+ 
      labs(y="% of Observations", color = expression(paste("Threshold (µmol ",O[2], " ",kg^-1,")")))+ 
      scale_color_manual(values = mycols) + theme(legend.position = c(0.9,0.8),
                                                  legend.direction = "vertical",
                                                  axis.text.x = element_text(colour="black"), 
                                                  axis.text.y = element_text(colour="black"))
    
### Make Fig S1 - DO and Temp time series for all sites ------------

    # Calculate deployment duration, start day, end day for all sites with processed data
    durations = DO_all %>% 
      group_by(site) %>% 
      summarize(dur = round(as.numeric(max(datetime) - min(datetime)),0),
                start = min(datetime),
                end = max(datetime))
    
    # make plotting function for DO
    supp_plot_DO <- function(sitename,breaks,dform) {
      plt = DO_all %>% 
        filter(site == sitename) %>%
        ggplot(aes(x=datetime, y=oxy)) + 
        geom_line()+
        theme_classic()+
        scale_x_datetime(date_breaks = breaks, labels = date_format(dform))+
        ylim(0,400)+ ggtitle(sitename) + theme(plot.title = element_text(vjust = - 5, hjust = 0.5, size = 10))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    }
    
    # make DO plots for each site, vary date format and days or weeks between date ticks
    he1 = supp_plot_DO("Heron 1","4 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank()) 
    he2 = supp_plot_DO("Heron 2","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    he3 = supp_plot_DO("Heron 3","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    boc1a = supp_plot_DO("Bocas 1a","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    boc1b = supp_plot_DO("Bocas 1b","5 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    boc1c = supp_plot_DO("Bocas 1c","5 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    boc1d = supp_plot_DO("Bocas 1d","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    boc2a = supp_plot_DO("Bocas 2a","4 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    boc2b = supp_plot_DO("Bocas 2b","4 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    tp1 = supp_plot_DO("Taiping 1","1 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    tp2 = supp_plot_DO("Taiping 2","1 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    oki1 = supp_plot_DO("Okinawa 1","6 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    oki2 = supp_plot_DO("Okinawa 2","4 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    oki3 = supp_plot_DO("Okinawa 3","3 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    dsha1 = supp_plot_DO("Dongsha 1","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    dsha2 = supp_plot_DO("Dongsha 2","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    dsha3 = supp_plot_DO("Dongsha 3","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    bak = supp_plot_DO("Baker","1 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    jar = supp_plot_DO("Jarvis","1 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    tut = supp_plot_DO("Tutuila","7 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    pal1a = supp_plot_DO("Palmyra 1a","2 month","%b") + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    pal1b = supp_plot_DO("Palmyra 1b","2 month","%b") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    pal2 = supp_plot_DO("Palmyra 2","2 month","%b") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    pal3 = supp_plot_DO("Palmyra 3","1 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    hog1a = supp_plot_DO("Hog 1a","1 week","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    hog1b = supp_plot_DO("Hog 1b","1 week","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    kb1a = supp_plot_DO("Kaneohe 1a","4 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    kb1b = supp_plot_DO("Kaneohe 1b","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    kb2 = supp_plot_DO("Kaneohe 2","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    cr1a = supp_plot_DO("Crocker 1a","2 week","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    cr1b = supp_plot_DO("Crocker 1b","2 week","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    cr1c = supp_plot_DO("Crocker 1c","1 month","%b") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    
    
    # plot all together (export 8x10in) 
    figs1a = tp1+tp2+bak+jar+pal3+boc1a+boc1d+he2+he3+dsha1+dsha2+dsha3+kb1b+kb2+oki3+boc2a+boc2b+oki2+he1+kb1a+boc1b+boc1c+
      oki1+hog1a+hog1b+tut+cr1b+cr1a+cr1c+pal1b+pal1a+pal2 + plot_layout(ncol = 6)
    
    # make function for plotting temp on right y-axis
    supp_plotT <- function(sitename,breaks,dform) {
      plt = DO_all %>% 
        filter(site == sitename) %>%
        ggplot(aes(x=datetime, y=temp)) + 
        geom_line(color = "blue")+
        theme_classic()+
        scale_y_continuous(position = "right", limits = c(19,33))+
        scale_x_datetime(date_breaks = breaks, labels = date_format(dform))+
        ggtitle(sitename) + theme(plot.title = element_text(vjust = - 5, hjust = 0.5, size = 10))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    }
    
    # make T plots for each site, vary date format and days or weeks between date ticks
    t_he1 = supp_plotT("Heron 1","4 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank()) 
    t_he2 = supp_plotT("Heron 2","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_he3 = supp_plotT("Heron 3","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_boc1a = supp_plotT("Bocas 1a","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    t_boc1b = supp_plotT("Bocas 1b","5 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_boc1c = supp_plotT("Bocas 1c","5 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_boc1d = supp_plotT("Bocas 1d","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_boc2a = supp_plotT("Bocas 2a","4 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_boc2b = supp_plotT("Bocas 2b","4 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_tp1 = supp_plotT("Taiping 1","1 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_tp2 = supp_plotT("Taiping 2","1 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_oki1 = supp_plotT("Okinawa 1","6 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_oki2 = supp_plotT("Okinawa 2","4 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    t_oki3 = supp_plotT("Okinawa 3","3 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_dsha1 = supp_plotT("Dongsha 1","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_dsha2 = supp_plotT("Dongsha 2","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_dsha3 = supp_plotT("Dongsha 3","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    t_bak = supp_plotT("Baker","1 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_jar = supp_plotT("Jarvis","1 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_tut = supp_plotT("Tutuila","7 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_pal1a = supp_plotT("Palmyra 1a","2 month","%b") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_pal1b = supp_plotT("Palmyra 1b","2 month","%b") + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    t_pal2 = supp_plotT("Palmyra 2","2 month","%b") + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    t_pal3 = supp_plotT("Palmyra 3","1 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_hog1a = supp_plotT("Hog 1a","1 week","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    t_hog1b = supp_plotT("Hog 1b","1 week","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_kb1a = supp_plotT("Kaneohe 1a","4 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_kb1b = supp_plotT("Kaneohe 1b","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_kb2 = supp_plotT("Kaneohe 2","2 day","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_cr1a = supp_plotT("Crocker 1a","2 week","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_cr1b = supp_plotT("Crocker 1b","2 week","%b %d") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    t_cr1c = supp_plotT("Crocker 1c","1 month","%b") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y =element_blank())
    
    
    # make t plots all together (export 8 x 10")
    figs1b = t_tp1+t_tp2+t_bak+t_jar+t_pal3+t_boc1a+t_boc1d+t_he2+t_he3+t_dsha1+t_dsha2+t_dsha3+t_kb1b+t_kb2+t_oki3+t_boc2a+t_boc2b+t_oki2+t_he1+t_kb1a+t_boc1b+t_boc1c+
      t_oki1+t_hog1a+t_hog1b+t_tut+t_cr1b+t_cr1a+t_cr1c+t_pal1b+t_pal1a+t_pal2 + plot_layout(ncol = 6)
    
### Create Table S5: % of Obs Hypoxic by site -----

# Calculate length of dataset for all sites    
deplength = DO_all_proj_30min_all %>% 
            group_by(site) %>% 
            summarize(length = length(oxy_proj))

# Create Table S5        
tableS5 = intens_hyp_thresholds %>% 
      group_by(projection, site) %>% 
      summarize(t153 = length(which(intens_153 > 0)),
                t122 = length(which(intens_122 > 0)),
                t92 = length(which(intens_92 > 0)),
                t61 = length(which(intens_61 > 0)))

tableS5 = merge(tableS5,deplength,by="site")

tableS5 = tableS5 %>%
          group_by(site,projection) %>% 
          summarize("153" = (t153/length)*100,
                 "122" = (t122/length)*100,
                 "92" = (t92/length)*100,
                 "61" = (t61/length)*100)

tableS5 = tableS5 %>% mutate(across(where(is.numeric), ~ round(., 1)))

### Calculations for abstract and MS body -----

## ABSTRACT

# Number of data points that are hypoxic for each site and projection
hyp_stats = DO_all_proj_30min_all %>%
  group_by(site,projection) %>% 
  summarize(less153 = length(which(hyp_153_ID>0)),
            less122 = length(which(hyp_122_ID>0)),
            less92 = length(which(hyp_92_ID>0)),
            less61 = length(which(hyp_61_ID>0)))

# Number of sites experiencing hypoxia (at least 1 data point)
hyp_stats_count = hyp_stats %>% 
  group_by(projection) %>% 
  summarise(less153 = length(which(less153>0)),
            less122 = length(which(less122>0)),
            less92 = length(which(less92>0)),
            less61 = length(which(less61>0)))

hyp_stats_frac = hyp_stats %>% 
  group_by(projection) %>% 
  summarise(less153 = length(which(less153>0))/length(less153),
            less122 = length(which(less122>0))/length(less122),
            less92 = length(which(less92>0))/length(less92),
            less61 = length(which(less61>0))/length(less61))

#### BODY --

## Calculate mean stats used in MS for DO concentration

# Mean +/- SD: Mean Daily Range [DO]
mean(table_sites_means$mean_dailyrangeDO)
sd(table_sites_means$mean_dailyrangeDO)

# Mean +/- SD: Mean Daily Min [DO]
mean(table_sites_means$mean_dailyminDO)
sd(table_sites_means$mean_dailyminDO)

# Mean +/- SD: Mean Daily Max [DO]
mean(table_sites_means$mean_dailymaxDO)
sd(table_sites_means$mean_dailymaxDO)

# Mean +/- SD: Mean Daily Mean [DO]
mean(table_sites_means$mean_dailymeanDO)
sd(table_sites_means$mean_dailymeanDO)

## Calculate mean stats used in MS for DO % sat

# Mean +/- SD: Mean Daily Range DO sat %
mean(table_sites_means$mean_dailyrangeDOsat)
sd(table_sites_means$mean_dailyrangeDOsat)

# Mean +/- SD: Mean Daily Min DO sat %
mean(table_sites_means$mean_dailyminDOsat)
sd(table_sites_means$mean_dailyminDOsat)

# Mean +/- SD: Mean Daily Max DO sat %
mean(table_sites_means$mean_dailymaxDOsat)
sd(table_sites_means$mean_dailymaxDOsat)

# Mean +/- SD: Mean Daily Mean DO sat %
mean(table_sites_means$mean_dailymeanDOsat)
sd(table_sites_means$mean_dailymeanDOsat)

## Calculate number of hypoxic events and other stats

# PRESENT-DAY

  # Total # of events by threshold (present-day only):
    events_hyp_all_tot %>% 
      filter(projection == "present")
    
  # Max and min durations of hypoxia at present by threshold
    events_hyp_all %>% 
      filter(projection == "present") %>% 
      group_by(threshold) %>% 
      summarize(max = max(dur_hrs),
                min = min(dur_hrs))
  
  # Look at specific sites that experience severe hypoxia and their event durations  
    events_hyp_all %>% 
      filter(projection == "present" && threshold == "61") %>% 
      group_by(threshold)
    
    
# SSP1-2.6
    
    # Total # of events by threshold (SP1-2.6 only):
    events_hyp_all_tot %>% 
      filter(projection == "ssp126")
    
    # Max and min durations of hypoxia at SP1-2.6 by threshold
    events_hyp_all %>% 
      filter(projection == "ssp126") %>% 
      group_by(threshold) %>% 
      summarize(max = max(dur_hrs),
                min = min(dur_hrs))
    
    # Look at specific sites that experience severe hypoxia and their event durations  
    events_hyp_all %>% 
      filter(projection == "ssp126" && threshold == "61") %>% 
      group_by(threshold)    

# SSP5-8.5
    
    # Total # of events by threshold (SSP5-8.5 only):
    events_hyp_all_tot %>% 
      filter(projection == "ssp585")
    
    # Max and min durations of hypoxia at SSSP5-8.5 by threshold
    events_hyp_all %>% 
      filter(projection == "ssp585") %>% 
      group_by(threshold) %>% 
      summarize(max = max(dur_hrs),
                min = min(dur_hrs))
    
    # Look at specific sites that experience severe hypoxia and their event durations  
    events_hyp_all %>% 
      filter(projection == "ssp585" && threshold == "61") %>% 
      group_by(threshold)    
    
### Depth, Current Speed, Benthic Community Analysis -----
    
## Wrangle the data for plotting:

# Create depthDO df with just site and mean daily range (with SD) and mean daily min (with SD) in DO at present-day
depthDO = table_sites_means[c("site","mean_dailyrangeDO","sd_dailyrangeDO","mean_dailyminDO","sd_dailyminDO")]

# Extract depth and location data for each site  
info = DO_all %>% group_by(site, loc) %>% summarise(depth = mean(depth))

# Combine into one depthDO df
depthDO = merge(depthDO, info, by = "site")

# Load in reef type categorization, ADCP speed data, and depth data
reef_ADCP_depth = read_csv("Data/reef_ADCP_depth.csv",show_col_types = "FALSE")

# Combine into depthDO df
depthDO = merge(depthDO, reef_ADCP_depth, by = "site")

# Re-order reef type classifications to go from offshore to inshore
depthDO = depthDO %>% mutate(reeftype = fct_relevel(reeftype,"Reef Front", "Terrace", "Outer Reef Flat", "Inner Reef Flat","Back Reef Slope","Lagoon","Shallow Lagoon","Patch Reef"))

## (1) MEAN DAILY RANGE IN DO vs DEPTH 
##...............................................

  # Create power model based on equation from MATLAB curve fitting: f(x) = a*x^b
  depthvals = seq(0,20, 0.1)
  powermod = data.frame(depthvals)
  powermod$est_DO = 127.7*depthvals^-0.4088
  
  # Where: SSE = 5.781e+04
  #  R-square: 0.3212
  #  Adjusted R-square: 0.2978
  #  RMSE: 44.65
  
  # Plot mean daily range in DO (present-day) vs depth for all sites (except Hog 1a), color by reef type
  
  dplot1 =  ggplot() +
    # Scatter and error bars
    geom_point(data = depthDO, aes(x = mean_depth, y = mean_dailyrangeDO, color = reeftype), size = 5) +
    geom_errorbar(data = depthDO, aes(x = mean_depth, ymin=mean_dailyrangeDO-sd_dailyrangeDO, 
                                      ymax=mean_dailyrangeDO+sd_dailyrangeDO, color = reeftype), width=.4)+
    geom_errorbar(data = depthDO, aes(y = mean_dailyrangeDO, xmin=mean_depth-sd_depth, 
                                      xmax=mean_depth+sd_depth, color = reeftype), width=8)+
    
    # Power model
    geom_line(data = powermod, aes(x = depthvals, y = est_DO), color = 'red')+
    
    labs(x = "Mean depth (m)", y = expression(paste("Mean daily range in DO (µmol ",O[2], " ",kg^-1,")")), color = "Reef type")+
    scale_color_paletteer_d("ggthemes::Green_Orange_Teal", direction = -1)+
    ylim(0,300) + xlim(0,20) +
    theme_classic() +  theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank(),
                             axis.text.x = element_text(colour="black"), 
                             axis.text.y = element_text(colour="black"),
                             axis.text = element_text(size = 12),
                             axis.title = element_text(size = 12),
                             legend.text = element_text(size = 12),
                             legend.title = element_text(size = 12),
                             legend.position = "none")+
    annotate("text", x = 0, y=300, label = "A", fontface = "bold", size = 7)

## (2) MEAN DAILY RANGE IN DO vs MEAN FLOW SPEED 
##...............................................
  # Power model could not be fit (very poor fit)
  
  #powermod$est_DO_flow = 71.43*depthvals^0.03977
  
  # Where: SSE = 1.587e+04
  #  R-square: 0.002156
  #  Adjusted R-square: -0.09763
  #  RMSE: 39.83
  
  # Plot mean daily range in DO (present-day) vs mean flow speed for sites with flow data, color by reef type
  dplot2 =  ggplot() +
    # Scatter and error bars
    geom_point(data = depthDO, aes(x = mean_speed, y = mean_dailyrangeDO, color = reeftype), size = 5) +
    geom_errorbar(data = depthDO, aes(x = mean_speed, ymin=mean_dailyrangeDO-sd_dailyrangeDO, 
                                      ymax=mean_dailyrangeDO+sd_dailyrangeDO, color = reeftype), width=.008)+
    geom_errorbar(data = depthDO, aes(y = mean_dailyrangeDO, xmin=mean_speed-sd_speed, 
                                      xmax=mean_speed+sd_speed, color = reeftype), width=10)+
    labs(x = expression(paste("Mean flow speed (m ",s^-1,")")), y = expression(paste("Mean daily range in DO (µmol ",O[2], " ",kg^-1,")")), color = "Reef type")+
    scale_color_paletteer_d("ggthemes::Green_Orange_Teal", direction = -1)+
    ylim(0,300) + xlim(0,0.4) +
    theme_classic() +  theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank(),
                             axis.text.x = element_text(colour="black"), 
                             axis.text.y = element_text(colour="black"),
                             axis.text = element_text(size = 12),
                             axis.title = element_text(size = 12),
                             legend.text = element_text(size = 12),
                             legend.title = element_text(size = 12),
                             legend.position = c(0.8,0.75))+
    annotate("text", x = 0, y=300, label = "B", fontface = "bold", size = 7)

## (3) MEAN DAILY MIN DO vs DEPTH
##...............................................

# Create power model based on equation from MATLAB curve fitting: f(x) = a*x^b
powermod$estminDO_depth = 117.9*depthvals^0.1006

# Where: SSE = 4.191e+04
#  R-square: 0.09983
#  Adjusted R-square: 0.06879
#  RMSE: 38.02

# Plot mean daily min DO (present-day) vs depth for all sites (except Hog 1a), color by reef type
dplot3 = ggplot() +
  # Scatter and error bars
  geom_point(data = depthDO, aes(x = mean_depth, y = mean_dailyminDO, color = reeftype), size = 5) +
  geom_errorbar(data = depthDO, aes(x = mean_depth, ymin=mean_dailyminDO-sd_dailyminDO, 
                                    ymax=mean_dailyminDO+sd_dailyminDO, color = reeftype), width=.4)+
  geom_errorbar(data = depthDO, aes(y = mean_dailyminDO, xmin=mean_depth-sd_depth, 
                                    xmax=mean_depth+sd_depth, color = reeftype), width=6)+
  
  # Power model
  geom_line(data = powermod, aes(x = depthvals, y = estminDO_depth), color = 'red')+
  
  labs(x = "Mean depth (m)", y = expression(paste("Mean daily min DO (µmol ",O[2], " ",kg^-1,")")), color = "Reef type")+
  scale_color_paletteer_d("ggthemes::Green_Orange_Teal", direction = -1)+
  ylim(0,220) + xlim(0,20) +
  theme_classic() +  theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           axis.text.x = element_text(colour="black"), 
                           axis.text.y = element_text(colour="black"),
                           axis.text = element_text(size = 12),
                           axis.title = element_text(size = 12),
                           legend.text = element_text(size = 12),
                           legend.title = element_text(size = 12),
                           legend.position = "none")+
  annotate("text", x = 0, y=220, label = "C", fontface = "bold", size = 7)

## (4) MEAN DAILY MIN DO vs MEAN FLOW SPEED
##...............................................

# Create power model based on equation from MATLAB curve fitting: f(x) = a*x^b
powermod$speedvals = seq(0,0.4,0.002)
powermod$estminDO_flow = 246.2*powermod$speedvals^0.2065

# Where: SSE: 6863
#  R-square: 0.4043
#  Adjusted R-square: 0.3447
#  RMSE: 26.2

# Plot mean daily min DO (present-day) vs mean flow speed for sites with flow data, color by reef type
dplot4 =  ggplot() +
  # Scatter and error bars
  geom_point(data = depthDO, aes(x = mean_speed, y = mean_dailyminDO, color = reeftype), size = 5) +
  geom_errorbar(data = depthDO, aes(x = mean_speed, ymin=mean_dailyminDO-sd_dailyminDO, 
                                    ymax=mean_dailyminDO+sd_dailyminDO, color = reeftype), width=.010)+
  geom_errorbar(data = depthDO, aes(y = mean_dailyminDO, xmin=mean_speed-sd_speed, 
                                    xmax=mean_speed+sd_speed, color = reeftype), width=6)+
  # Power model
  geom_line(data = powermod, aes(x = speedvals, y = estminDO_flow), color = 'red')+
  
  labs(x = expression(paste("Mean flow speed (m ",s^-1,")")), y = expression(paste("Mean daily min DO (µmol ",O[2], " ",kg^-1,")")), color = "Reef type")+
  scale_color_paletteer_d("ggthemes::Green_Orange_Teal", direction = -1)+
  ylim(0,220) + xlim(0,0.4) +
  theme_classic() +  theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           axis.text.x = element_text(colour="black"), 
                           axis.text.y = element_text(colour="black"),
                           axis.text = element_text(size = 12),
                           axis.title = element_text(size = 12),
                           legend.text = element_text(size = 12),
                           legend.title = element_text(size = 12),
                           legend.position = "none")+
  annotate("text", x = 0, y=220, label = "D", fontface = "bold", size = 7)

(dplot1|dplot2)/(dplot3|dplot4)  
ggsave("~/Documents/SIO/O2 Variability/Figures/Supplement/Depth Analysis/4Panel.pdf", width = 7, height = 7, units = "in", encoding = "MacRoman")

