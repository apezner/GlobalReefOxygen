# - - - - - - - - - - - - - - - - - - - - - - - -
# Code for Pezner et al. MS - Create Figure 1 (map, climatology, distribution)
# Last updated, January 10, 2023
# Written by: Ariel Pezner (apezner@ucsd.edu)
# - - - - - - - - - - - - - - - - - - - - - - - - 

# Set working directory
setwd("~/Documents/SIO/O2 Variability")

# Load necessary packages
library(ggplot2);library(scales);library(ggridges);library(dplyr);library(ggthemes);library(patchwork);library(gsw);library(RColorBrewer)
library(tidyverse);library(readr);library(DescTools);library(ggmap);library(pastecs);library(forcats);library(reshape2);library(lubridate);
library(tidyr);library(rworldmap);library(readxl); library(mapproj)

# Load color palette
library(wesanderson)

# Load data
load("RDataFiles/DO_all_Pezner.Rda") # all DO data for all sites

## Expand existing color palettes to fit this data
nb.cols <- 12
mycolors <- colorRampPalette(rev(wes_palette("Zissou1", 100, type = "continuous")))(nb.cols)

## Make Ridgeline
    # DO []
    ridge1 <- DO_all %>%
      mutate(loc = fct_relevel(loc, levels = "Bocas del Toro","Crocker Reef","Hog Reef","Heron Island",
                               "Dongsha","Okinawa","Taiping","Baker","Tutuila","Jarvis","Palmyra","Kaneohe Bay")) %>% 
      ggplot() +
      # Hypoxia rectangles
      #annotate("rect", xmin=0, xmax=61, ymin=-Inf, ymax=Inf, alpha=0.5, fill = "#FF9999")+
      #annotate("rect", xmin=0, xmax=92, ymin=-Inf, ymax=Inf, alpha=0.4, fill = "#FF9999")+
      #annotate("rect", xmin=0, xmax=122, ymin=-Inf, ymax=Inf, alpha=0.3, fill = "#FF9999")+
      #annotate("rect", xmin=0, xmax=153, ymin=-Inf, ymax=Inf, alpha=0.2, fill = "#FF9999")+
      # Add horizontal lines to denote thresholds
      geom_vline(xintercept = 153, linetype = "dashed", color = "black", alpha = 0.6)+
      geom_vline(xintercept = 122, linetype = "dashed", color = "black", alpha = 0.6)+
      geom_vline(xintercept = 92, linetype = "dashed", color = "black", alpha = 0.6)+
      geom_vline(xintercept = 61, linetype = "dashed", color = "black", alpha = 0.6)+
      
      # Ridges
      geom_density_ridges(aes(x=oxy, y=loc, fill=loc, color = loc, group = site),scale=1.5,quantile_lines=FALSE, alpha=0.8) +
      theme_ridges() + 
      xlim(0,475)+
      theme(legend.position = "none")+
      labs(x="DO (µmol/kg)", y="Location") +
      scale_fill_manual(values  = mycolors)+
      scale_color_manual(values  = mycolors)+
      theme(axis.text.x = element_text(colour="black"), 
            axis.text.y = element_text(colour="black"))
    
    # DOsat% 
    ridge2 <-  DO_all %>%
      mutate(loc = fct_relevel(loc, levels = "Bocas del Toro","Crocker Reef","Hog Reef","Heron Island",
                               "Dongsha","Okinawa","Taiping","Baker","Tutuila","Jarvis","Palmyra","Kaneohe Bay")) %>% 
      ggplot() +
      geom_density_ridges(aes(x=oxysat, y=loc, fill=loc, color = loc, group = site),scale=1.5,quantile_lines=FALSE, alpha=0.8) +
      theme_ridges() + 
      xlim(0,255)+
      theme(legend.position = "none",axis.title.y=element_blank(),axis.text.y=element_blank())+
      labs(x="DO Saturation (%)", y="") +
      scale_fill_manual(values  = mycolors)+
      scale_color_manual(values  = mycolors)+
      theme(axis.text.x = element_text(colour="black"),
            axis.text.y=element_blank())
    
    ridge = (ridge1|ridge2)
    
## Make Map    
    # Load in site names
    locations <- read_excel("Data/SiteLocations.xlsx")
    
    # Fix locations longitude if pacific centered map
    locations$lon2 <- ifelse(locations$Lon < -25, locations$Lon + 360, locations$Lon)
    
    # Load in longitude and latitude for labels (add or subtract to adjust placement of label relative to point plotted)
    labs <- data.frame(loc=locations$Location,
                       long=c(locations$lon2[1]-28, #Bocas
                              locations$lon2[2]+24, #Crocker
                              locations$lon2[3]+19, #Hog
                              locations$lon2[4]+24, #Heron
                              locations$lon2[5]+24, #Dsha
                              locations$lon2[6]+20, #OKI
                              locations$lon2[7]+24, #TP
                              locations$lon2[8]-13, #Baker
                              locations$lon2[9]+15, #Tutuila
                              locations$lon2[10]+14,#Jarvis
                              locations$lon2[11],   #pal
                              locations$lon2[12]),#Kbay
                       lat=c(locations$Lat[1],
                             locations$Lat[2],
                             locations$Lat[3],
                             locations$Lat[4],
                             locations$Lat[5],
                             locations$Lat[6],
                             locations$Lat[7],
                             locations$Lat[8],
                             locations$Lat[9],
                             locations$Lat[10],
                             locations$Lat[11]+6,
                             locations$Lat[12]+6))
    
    
    # Create plot (pacific centered)
    #world <- map_data("world")
    mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))
    
    map = ggplot() + 
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = 'grey') + 
      coord_fixed(1.3) +
      xlab('Longitude')+
      ylab('Latitude')+
      geom_text(data=labs, aes(long, lat, label = loc),color=mycolors,size=4)+
      geom_point(data=locations, aes(lon2, Lat),color=mycolors,size=2)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.text.y = element_blank(), axis.line =  element_blank(),
            legend.key = element_blank(), axis.text.x = element_blank(),
            axis.title.x =element_blank(),axis.title.y = element_blank(),axis.ticks=element_blank())+
      coord_quickmap(xlim = c(100,335), expand = FALSE)
    

    
## Now make climatology plot
    
    DO_all$hour = hour(DO_all$datetime)
    
    # Using ALL of the data (even for long time series)
    climatology_sum_all = DO_all %>%
      group_by(loc,site,hour) %>% #using that data, group it first by site and then each site by hour of the day
      summarize(oxy_mean=MeanCI(oxy,conf.level=0.95,na.rm=TRUE)[1], #calculate the mean value for the hourly data from each site [DO]
                oxy_lwr=MeanCI(oxy,conf.level=0.95,na.rm=TRUE)[2], #calculate the lower 95% confidence interval value for the hourly data from each site [DO]
                oxy_upr=MeanCI(oxy,conf.level=0.95,na.rm=TRUE)[3]) #calculate the upper 95% confidence interval value for the hourly data from each site [DO]

    # Plot climatologies, all on one
    clim= climatology_sum_all %>%
      mutate(loc = fct_relevel(loc, levels = "Bocas del Toro","Crocker Reef","Hog Reef","Heron Island",
                               "Dongsha","Okinawa","Taiping","Baker","Tutuila","Jarvis","Kaneohe Bay",
                               "Palmyra")) %>% 
      ggplot() + 
      # Add horizontal lines to denote thresholds
      geom_hline(yintercept = 153, linetype = "dashed", color = "black", alpha = 0.6)+
      geom_hline(yintercept = 122, linetype = "dashed", color = "black", alpha = 0.6)+
      geom_hline(yintercept = 92, linetype = "dashed", color = "black", alpha = 0.6)+
      geom_hline(yintercept = 61, linetype = "dashed", color = "black", alpha = 0.6)+
      
      # Add data
      geom_ribbon(aes(x=hour,ymin=oxy_lwr, ymax=oxy_upr, group=site,fill=loc), alpha=0.3)+ #adds the confidence intervals grouped station filled by the station id for month
      geom_line(aes(x=hour, y=oxy_mean, group=site,color=loc),size=1, alpha=0.8) + #adds a line for the mean monthly temperatures
      
      theme_classic()+
      theme(legend.position="none")+ 
      labs(x="Hour of the Day",y="DO (µmol/kg)") +
      #ylim(20,300)+
      #xlim(0,24)+
      scale_x_continuous(expand = c(0, 0), breaks= c(0,6,12,18,23)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      scale_fill_manual(values  = mycolors)+
      scale_color_manual(values  = mycolors)+
      theme(axis.line = element_line(colour = "black"),
            axis.text.x = element_text(colour="black"), 
            axis.text.y = element_text(colour="black"),
            axis.ticks = element_line(colour = "black"))

## Figure 1 components in two parts
  map 
  (clim|ridge)
  