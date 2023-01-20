# - - - - - - - - - - - - - - - - - - - - - - - -
# Code for Pezner et al. MS - extracting CMIP6 data
# Last updated, January 10, 2023
# Written by: Ariel Pezner (apezner@ucsd.edu), Travis Courtney
# - - - - - - - - - - - - - - - - - - - - - - - - 

#load required packages
library(tidyverse)
library(tidync)
library(readxl)
library(patchwork)

# Set working directory (will vary by computer)
setwd("~/Documents/SIO/O2 Variability/CMIP")

#import dissolved oxygen locations and convert to degrees east format for use with CMIP6 output data
locs=read_xlsx("SiteLocations.xlsx")
locs$lat=locs$Lat
locs1=subset(locs,Lon>0)
locs1$lon=locs1$Lon
locs2=subset(locs,Lon<0)
locs2$lon=locs2$Lon+360
locations=bind_cols(loc=bind_rows(locs1,locs2)$loc,latitude=bind_rows(locs1,locs2)$lat,longitude=bind_rows(locs1,locs2)$lon)

#import CESM simulations
ssp_126_1 <- tidync("tos_Omon_CESM2-WACCM_ssp126_r1i1p1f1_gr_201501-210012.nc")
ssp_245_1 <- tidync("tos_Omon_CESM2-WACCM_ssp245_r1i1p1f1_gr_201501-206412.nc")
ssp_245_2 <- tidync("tos_Omon_CESM2-WACCM_ssp245_r1i1p1f1_gr_206501-210012.nc")
ssp_370_1 <- tidync("tos_Omon_CESM2-WACCM_ssp370_r1i1p1f1_gr_201501-206412.nc")
ssp_370_2 <- tidync("tos_Omon_CESM2-WACCM_ssp370_r1i1p1f1_gr_206501-210012.nc")
ssp_585_1 <- tidync("tos_Omon_CESM2-WACCM_ssp585_r1i1p1f1_gr_201501-210012.nc")

#establish null values prior to looping through locs
ssp_126_warming=NULL
ssp_245_warming=NULL
ssp_370_warming=NULL
ssp_585_warming=NULL

#make empty plot list
modelplots = list()

#loop through the locations to find coordinates in CESM model outputs that most closely match locs
for (i in 1:nrow(locations))
{
ssp_126=ssp_126_1 %>%  hyper_filter(lon = index == which.min(abs(lon - locations$longitude[i])), lat = index == which.min(abs(lat - locations$latitude[i]))) %>% hyper_tibble() #create tibble of Hog Reef SST projections
ssp_126$decy=ssp_126$time/365+1 #converting from days since 0001-01-01 00:00:00 with no leap years to decimal years
ssp_126$year=as.numeric(substr(ssp_126$decy, 1, 4)) #create annual bin to use for annual means
ssp_126_annual=ssp_126 %>% group_by(year) %>% summarize(sst = mean(tos)) #projected annual SST means

ssp_245_1_sst=ssp_245_1 %>%  hyper_filter(lon = index == which.min(abs(lon - locations$longitude[i])), lat = index == which.min(abs(lat - locations$latitude[i]))) %>% hyper_tibble() #create tibble of Hog Reef SST projections
ssp_245_2_sst=ssp_245_2 %>%  hyper_filter(lon = index == which.min(abs(lon - locations$longitude[i])), lat = index == which.min(abs(lat - locations$latitude[i]))) %>% hyper_tibble() #create tibble of Hog Reef SST projections
ssp_245=bind_rows(ssp_245_1_sst,ssp_245_2_sst) #bind beginning and end of century data frames
ssp_245$decy=ssp_245$time/365+1 #converting from days since 0001-01-01 00:00:00 with no leap years to decimal years
ssp_245$year=as.numeric(substr(ssp_245$decy, 1, 4)) #create annual bin to use for annual means
ssp_245_annual=ssp_245 %>% group_by(year) %>% summarize(sst = mean(tos)) #projected annual SST means

ssp_370_1_sst=ssp_370_1 %>%  hyper_filter(lon = index == which.min(abs(lon - locations$longitude[i])), lat = index == which.min(abs(lat - locations$latitude[i]))) %>% hyper_tibble() #create tibble of Hog Reef SST projections
ssp_370_2_sst=ssp_370_2 %>%  hyper_filter(lon = index == which.min(abs(lon - locations$longitude[i])), lat = index == which.min(abs(lat - locations$latitude[i]))) %>% hyper_tibble() #create tibble of Hog Reef SST projections
ssp_370=bind_rows(ssp_370_1_sst,ssp_370_2_sst) #bind beginning and end of century data frames
ssp_370$decy=ssp_370$time/365+1 #converting from days since 0001-01-01 00:00:00 with no leap years to decimal years
ssp_370$year=as.numeric(substr(ssp_370$decy, 1, 4)) #create annual bin to use for annual means
ssp_370_annual=ssp_370 %>% group_by(year) %>% summarize(sst = mean(tos)) #projected annual SST means

ssp_585=ssp_585_1 %>%  hyper_filter(lon = index == which.min(abs(lon - locations$longitude[i])), lat = index == which.min(abs(lat - locations$latitude[i]))) %>% hyper_tibble() #create tibble of Hog Reef SST projections
ssp_585$decy=ssp_585$time/365+1 #converting from days since 0001-01-01 00:00:00 with no leap years to decimal years
ssp_585$year=as.numeric(substr(ssp_585$decy, 1, 4)) #create annual bin to use for annual means
ssp_585_annual=ssp_585 %>% group_by(year) %>% summarize(sst = mean(tos)) #projected annual SST means

ssp_126_warming[i]=mean(subset(ssp_126,year>2089)$tos)-mean(subset(ssp_126,year<2021)$tos) #projected delta annual temperatures between present and end of century
ssp_245_warming[i]=mean(subset(ssp_245,year>2089)$tos)-mean(subset(ssp_245,year<2021)$tos) #projected delta annual temperatures between present and end of century
ssp_370_warming[i]=mean(subset(ssp_370,year>2089)$tos)-mean(subset(ssp_370,year<2021)$tos) #projected delta annual temperatures between present and end of century
ssp_585_warming[i]=mean(subset(ssp_585,year>2089)$tos)-mean(subset(ssp_585,year<2021)$tos) #projected delta annual temperatures between present and end of century

# Plot temp rise over time for each loc
# make plot
plot = ggplot()+
        geom_line(data=ssp_126_annual, aes(year, sst),size=1,color="#2c7bb6")+
        geom_line(data=ssp_245_annual, aes(year, sst),size=1,color="#abd9e9")+
        geom_line(data=ssp_370_annual, aes(year, sst),size=1,color="#fdae61")+
        geom_line(data=ssp_585_annual, aes(year, sst),size=1,color="#d7191c")+
        ylab("SST (°C)")+xlab("Year")+ggtitle(locations$loc[i])+
        theme_classic()+theme(text = element_text(size = 7)) + ylim(24,35)

# save to list and rename as loc
modelplots[[i]] = plot
names(modelplots)[[i]] <- locations$loc[i]

}

# create dataframe with loc info and warming
locs_warming=bind_cols(loc=locations$loc,ssp126=ssp_126_warming,ssp245=ssp_245_warming,ssp370=ssp_370_warming,ssp585=ssp_585_warming)

# Save plots to folder
 for (i in 1:nrow(locations)) {
   file_name = paste("Plots/",locations$loc[i], ".pdf", sep="")
   pdf(file_name, width = 6, height = 4)
   print(modelplots[[i]])
   dev.off()
 }

# make big plot
modelplots[[1]] + modelplots[[2]] + modelplots[[3]] + modelplots[[4]] + modelplots[[5]] + modelplots[[6]] + modelplots[[7]] + 
modelplots[[8]] + modelplots[[9]] + modelplots[[10]] + modelplots[[11]] + modelplots[[12]] +  plot_layout(ncol = 4)

# save data and export
#write.csv(locs_warming, "~/Documents/SIO/O2 Variability/Data/cmip6_sst_increase.csv", row.names=FALSE)

