rm(list=ls())

# load packages
library(sf)
library(mapplots)
library(scatterpie)
library(maptools)
library(plyr) 
library(dplyr) 
library(tidyr)
library(ggplot2)
library(ggmap)
library(mapproj)
library(ggnewscale)
library(ggspatial)
library(ggrepel)

# set wd
homewd = "/Users/flg9/Desktop/Developer/brook_lab/R/madamap/Mada-GIS/" 
setwd(paste0(homewd))

# activate API key 
register_google(key = "AIzaSyC5NF4xccTOw2UvjV4PRIElzrX8vr5Yv6U")

# create base mada map
madamap <- get_map(location = 'Madagascar',
                   zoom = 5, scale = "auto", maptype = "satellite") 

ggmap(madamap, padding = 0.02) + xlab("Longitude") + ylab("Latitude") 

# load in prevalence data 
dat <- read.csv(file = paste0(homewd,"all_NGS_distribute.csv"), 
                header = T, stringsAsFactors = F )
head(dat)
names(dat)

# subset for fecal data 
dat_fec = subset(dat, sample_type=="feces")

# add age class
# clean class
unique(dat$bat_age_class)

#and rank by rough age
unique(dat_fec$young_of_year)
dat_fec$age_class <- dat_fec$bat_age_class
dat_fec$age_class[dat_fec$age_class=="P" | dat_fec$age_class=="L"] <- "A"
dat_fec$age_class[dat_fec$age_class=="NL" | dat_fec$young_of_year=="no"] <- "A"
dat_fec$age_class[dat_fec$young_of_year=="yes"] <- "J"

# subset to only inlcude columns of interest
dat_fec <- dplyr::select(dat_fec,roost_site,latitude_s, longitude_e,
                         collection_date, age_class, bat_sex,
                         species, sampleid, KoV)
head(dat)
unique(dat$roost_site)

# get sites and assign bat species to each
coordinate <- ddply(dat_fec, .(roost_site), summarise, latitude_s=unique(latitude_s), longitude_e=unique(longitude_e))
coordinate <-subset(coordinate, 
                    roost_site=="Ambakoana" | 
                      roost_site=="AngavoKely" |
                      roost_site=="Maromizaha")
coordinate$species <- c("Pteropus rufus", "Eidolon dupreanum", 
                        "Rousettus madagascariensis")
head(coordinate)

# take a glance with new coordinate data 
p2 <- ggmap(madamap) + 
  geom_point(data=coordinate,aes(x = longitude_e, y = latitude_s, color = "white"), 
             color = 'white', size = 1.2) + xlab("Longitude") + ylab("Latitude") +
  annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
                         pad_x = unit(1, "cm"), 
                         pad_y = unit(1, "cm"),        
                         style = north_arrow_nautical(text_col = "white", 
                                                      text_size = 22))
p2

# add bat species labels to map
coordinate$label <- coordinate$species
coordinate$label[coordinate$label=="Rousettus madagascariensis"] <- "Rousettus\nmadagascariensis"
coordinate$label[coordinate$label=="Eidolon dupreanum"] <- "Eidolon\ndupreanum"
coordinate$label[coordinate$label=="Pteropus rufus"] <- "Pteropus\nrufus"

# load GPS point and label
p3 <- p2+ geom_point(aes(x=longitude_e, y=latitude_s),color="white",size=1,data=coordinate) + 
  geom_text(data= coordinate,                       #### Labeling
            aes(x=longitude_e, y=latitude_s, label=label),
            fontface="italic",
            color = "white", size=6,
            nudge_x = c(-2,-7, 8),
            nudge_y = c(4.5,-2,3.5),
            check_overlap = F) #removing overlap ensures all 3 species appear

p3

# may not be necessary 
p2b<- p2+ geom_point(aes(x=longitude_e, y=latitude_s),color="black",size=1,data=coordinate)+
  geom_text(data= coordinate,                       #### Labeling
            aes(x=longitude_e, y=latitude_s, label=label),
            fontface="italic",
            color = "#1B262C", size=3,
            nudge_x = c(-2,-3.6,8),
            nudge_y = c(3,-1.1,-.3),
            check_overlap = T) +
  annotation_scale(location = "bl", width_hint = 0.05) +    #scale
  annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
                         pad_x = unit(0.03, "cm"), 
                         pad_y = unit(0.2, "cm"),        
                         style = north_arrow_fancy_orienteering)+
  geom_text_repel(segment.colour="black")+
  theme_bw() +theme(panel.grid = element_blank(), 
                    plot.title = element_text(color="black", size=12, face="bold"),
                    plot.margin = unit(c(-1,.5,-1.5,.1),"cm"),
                    axis.title.x = element_text(color="black", size=12),
                    axis.title.y = element_text(color="black", size=12),
                    legend.position=c(.26,.90),
                    legend.margin = margin(),
                    legend.title=element_blank(),
                    legend.background = element_rect(color="gray",size = .1),
                    legend.text = element_text(size = 9,face = "italic"))

p2b

# plot one site per species according to where positives are found
dat_fec$roost_site[dat_fec$species=="Pteropus rufus"] <- "Ambakoana"
dat_fec$roost_site[dat_fec$species=="Eidolon dupreanum"] <- "AngavoKely"
dat_fec$roost_site[dat_fec$species=="Rousettus madagascariensis"] <- "Maromizaha"

# summarize longitude data for each site
dat_fec$longitude_e[dat_fec$roost_site=="Ambakoana"] <- coordinate$longitude_e[coordinate$roost_site=="Ambakoana"]
dat_fec$longitude_e[dat_fec$roost_site=="AngavoKely"] <- coordinate$longitude_e[coordinate$roost_site=="AngavoKely"]
dat_fec$longitude_e[dat_fec$roost_site=="Maromizaha"] <- coordinate$longitude_e[coordinate$roost_site=="Maromizaha"]

# summarize longitude data for each site
dat_fec$latitude_s[dat_fec$roost_site=="Ambakoana"] <- coordinate$latitude_s[coordinate$roost_site=="Ambakoana"]
dat_fec$latitude_s[dat_fec$roost_site=="AngavoKely"] <- coordinate$latitude_s[coordinate$roost_site=="AngavoKely"]
dat_fec$latitude_s[dat_fec$roost_site=="Maromizaha"] <- coordinate$latitude_s[coordinate$roost_site=="Maromizaha"]

# i tried this initially, but there are only two positives so not really relevant
dat_fec$plot_class <- NA
dat_fec$plot_class[dat_fec$age_class=="J" & dat_fec$KoV==1] <- "juvenile: KoV pos"
dat_fec$plot_class[dat_fec$age_class=="J" & dat_fec$KoV==0] <- "juvenile: KoV neg"
dat_fec$plot_class[dat_fec$age_class=="A" & dat_fec$KoV==1] <- "adult: KoV pos"
dat_fec$plot_class[dat_fec$age_class=="A" & dat_fec$KoV==0] <- "adult: KoV neg"


# final grouping for scatter pie
dat_fec$plot_class <- NA
dat_fec$plot_class[dat_fec$KoV==1] <- "KoV Pos"
dat_fec$plot_class[dat_fec$KoV==0] <- "KoV Neg"

pies <- ddply(dat_fec, .(species, roost_site, latitude_s, longitude_e, plot_class), summarise, value=length(sampleid))

tot_sum = ddply(pies,.(species), summarise,N=sum(value))

pies <- merge(pies, tot_sum, by=c("species"), all.x=T)

pies$plot_class <- factor(pies$plot_class, levels=c( "KoV Pos", "KoV Neg"))

# add colors 
colz = c('KoV Pos' ="#FF6B6B", 
         'KoV Neg' ="#46B4AF")

# now let's make our pies

#tot_coord <- join(coordinate, tot_sum, by = "species", match = "all", 
# type = "right") # df containing total counts for each pie

p4<-ggplot() + 
  geom_scatterpie(aes(x=longitude_e, y=latitude_s, r=(N/1000)), 
                  data = pies, cols="plot_class", long_format=TRUE) +
  geom_scatterpie(aes(x=longitude_e, y=latitude_s, r=(N/1000)), 
                  data = pies, cols="plot_class", long_format=TRUE) +
  scale_fill_manual(values=colz)

p4 # glance at pies


# # copy of latitude (x.) and longitude (y.)
pies$x2 <- pies$longitude_e
pies$y2 <- pies$latitude_s

# #manually move the pie chart in case there is an overlap (change x and y)
# 
pies$x2[pies$species== "Pteropus rufus"] <- pies$longitude_e[pies$species== "Pteropus rufus"] -2
pies$y2[pies$species== "Pteropus rufus"] <- pies$latitude_s[pies$species== "Pteropus rufus"] + 1


pies$x2[pies$species== "Eidolon dupreanum"] <- pies$longitude_e[pies$species== "Eidolon dupreanum"] - .5
pies$y2[pies$species== "Eidolon dupreanum"] <- pies$latitude_s[pies$species== "Eidolon dupreanum"] - 3

pies$x2[pies$species== "Rousettus madagascariensis"] <- pies$longitude_e[pies$species== "Rousettus madagascariensis"] + 4
pies$y2[pies$species== "Rousettus madagascariensis"] <- pies$latitude_s[pies$species== "Rousettus madagascariensis"] - 0

head(pies)

# this is Fig1A
p5 <- p3+
  annotate("segment", x=pies$longitude_e, xend=pies$x2,y=pies$latitude_s,yend=pies$y2,size=.7)+ # put the lines
  annotate("segment", x=pies$longitude_e, xend=pies$x2,y=pies$latitude_s,yend=pies$y2,size=.7)+ # put the lines
  geom_scatterpie(aes(x=x2, y=y2, r=(log10(N)/1.2)), 
                  data = pies, cols="plot_class", long_format=TRUE) +
  geom_scatterpie(aes(x=x2, y=y2, r=(log10(N)/1.2)), 
                  data = pies, cols="plot_class", long_format=TRUE) +
  theme_bw() +theme(panel.grid = element_blank(),
                    plot.title = element_text(color="black", size=12, face="bold"),
                    plot.margin = unit(c(-1,.5,-1.5,.1),"cm"),
                    axis.title.x = element_text(color="black", size=12),
                    axis.title.y = element_text(color="black", size=12),
                    legend.position=c(.8,.8),
                    legend.margin = margin(),
                    legend.title=element_blank(),
                    legend.text = element_text(size = 7.5)) +
  scale_fill_manual(values=colz)
# geom_scatterpie_legend(log10(c(10,100)/1.2),
# x=54.5, y=-23.5, 
#  n=2,
# labeller = function(x) paste(10^(x)*1.2,"indiv"))



p5

homewd="/Users/flg9/Desktop/"
setwd(paste0(homewd))


ggsave(file = paste0(homewd, "/final-figures/Fig1_final.pdf"),
       plot=p5,
       units="mm",  
       width=60, 
       height=60, 
       scale=3, 
       dpi = 300)











