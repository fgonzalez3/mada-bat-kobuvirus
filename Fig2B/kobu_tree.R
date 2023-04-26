rm(list=ls())

# make fig 1b - nt full genome kobuvirus

# load in necessary packages 

library(ggplot2)
library(ggtree)
library(ape)
library(plyr)
library(phytools)
library(phangorn)
library(ggnewscale)
library(phylobase)
library(stringr)
library(tidyr)
library(ggpubr)
library(grid)
library(gridExtra)
library(dplyr)
library(aptheme)


# set working directory 

homewd= '/Users/flg9/Desktop/Developer/brook_lab/kobu_tree/'

setwd(paste0(homewd))

# load fig 1b tree 
final.kobu <-  read.tree(file = paste0(homewd, "kobu.supportFBP"))

# change Node_4 to Genbank name
final.kobu$tip.label <- gsub("NODE_4", "OP287812", final.kobu$tip.label)

# root it using a rabovirus as an outgroup
final.rooted.kobu <- root(final.kobu, which(final.kobu$tip.label=='NC_026314'))


# take a quick look in base R
ggtree(final.rooted.kobu) + 
  geom_nodelab(aes(label=label), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align= FALSE, linetype="dotted", linesize = 0.1, size = 2.5) + geom_point(colour='red')


# load manual tree csv for kobuvirus 
# and fill in for NA values 

kobu.manual <- read.csv(file=paste0('kobuvirus_manual2.csv'), 
                        header=T, stringsAsFactors = F, na = "")

kobu.manual[is.na(kobu.manual)] = "NA"

# use dplyr to only include columns that will be in the final tip label

colnames(kobu.manual)

kobu.manual <- kobu.manual %>%
  dplyr::select(Accession_Number, Species, Geo_Location, host, Year)


# check unique hosts that will be used to color tip labels 
# and assign colors to them

unique(kobu.manual$host)

colz = c("Human Kobuvirus" = "red", "Bovine Kobuvirus" = "tomato", 
         "Porcine Kobuvirus" = "khaki", "Ovine Kobuvirus" = "yellow", 
         "Canine Kobuvirus" = "blue", "Rodent Kobuvirus" = "plum", 
         "Sewage Kobuvirus" = "pink", "Caprine Kobuvirus" = "darkolivegreen1",
         "Feline Kobuvirus" = "purple", "Avian Kobuvirus" = "darksalmon",  
         "Rabovirus" = "brown", "Rabbit Kobuvirus" = "darkgreen", 
         "Bat Kobuvirus" = "cyan", "Bat Picornavirus" = "black")


# pick order for the labels
C =c("Human Kobuvirus","Bovine Kobuvirus","Porcine Kobuvirus","Ovine Kobuvirus",
     "Canine Kobuvirus","Rodent Kobuvirus","Sewage Kobuvirus","Caprine Kobuvirus",
     "Feline Kobuvirus","Avian Kobuvirus","Rabovirus", "Rabbit Kobuvirus",
     "Bat Kobuvirus", "Bat Picornavirus")


# visualize again with added legend order 

ggtree(final.rooted.kobu) %<+% kobu.manual + geom_tippoint(aes(fill=host)) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  scale_fill_manual(values=colz) + theme(legend.position = c(.2,.85), legend.title = element_blank())

# add a "novel" category
# our novel virus will be highlighted on the tree

kobu.manual$novel = 0
kobu.manual$novel[kobu.manual$Accession=="OP287812"] <- 1

kobu.manual$novel <- as.factor(kobu.manual$novel)


# add a bat host category
# bat hosts will have different tip shapes 

kobu.manual$bat_host <- 0

kobu.manual$bat_host[kobu.manual$Accession=="OP287812"] <- "bat-host"
kobu.manual$bat_host[kobu.manual$Accession=="MN602325"] <- "bat-host"
kobu.manual$bat_host[kobu.manual$Accession=="MF352427"] <- "bat-host"
kobu.manual$bat_host[kobu.manual$Accession=="KJ641686"] <- "bat-host"
kobu.manual$bat_host[kobu.manual$Accession=="KJ641691"] <- "bat-host"



kobu.manual$bat_host[kobu.manual$bat_host==0] <- "non-bat-host"
kobu.manual$bat_host[kobu.manual$bat_host==1] <- "bat-host"

# give bat host category a shape and color 

kobu.manual$bat_host <- as.factor(kobu.manual$bat_host)
shapez = c("bat-host" =  24, "non-bat-host" = 21)
colz2 = c('1' =  "yellow", '0' = "white")

# clade labels
AichiA <- getMRCA(final.rooted.kobu, c("OP287812", "KJ934637"))
AichiB <- getMRCA(final.rooted.kobu, c("MN336260", "GU245693"))
AichiC <- getMRCA(final.rooted.kobu, c("NC_023422", "KY234500"))
AichiD <- getMRCA(final.rooted.kobu, c("NC_027918", "NC_027919"))
AichiEF <- getMRCA(final.rooted.kobu, c("KJ641686", "KT325852"))






# tip labels on manual csv and raxml tree do not match, so fix that 
# first make a new df with the labels numbered by root tip on the raxml tree

kobu.dat <- data.frame(Accession_Number=final.rooted.kobu$tip.label, 
                       num =1:length(final.rooted.kobu$tip.label))

# then create a 3rd df and right join the original manual and the 2nd df 
# by accession_number onto this 3rd df 

kobux <- join(kobu.manual, kobu.dat, by = "Accession_Number", match = "all", 
              type = "right")


kobux$new_label <- kobux$Accession_Number


# match tip labels from tree with csv file 
#kobux <- kobu.manual[match(kobu.dat$new_label, kobu.manual$new_label),]
#kobux <- kobux[!is.na(kobux$Accession_Number), ]


#kobux$Accession <- kobux$new_label
#final.rooted.kobu$tip.label <- kobux$Accession


# original 
# create new tip labels

#kobux$Accession_Number<- str_replace_all(kobux$Accession_Number, "NA", "")


kobux$new_label[!is.na(kobux$new_label)] <- paste(kobux$Accession_Number[!is.na(kobux$Accession_Number)], " | ", 
                                                        kobux$Species[!is.na(kobux$Species)], " | ",
                                                        kobux$host[!is.na(kobux$host)], " | ",
                                                        kobux$Geo_Location[!is.na(kobux$Geo_Location)], " | ",
                                                        kobux$Year[!is.na(kobux$Year)])

kobux$new_label <- paste(kobux$Accession_Number, " | ", 
                                                  kobux$Species, " | ",
                                                  kobux$host, " | ",
                                                  kobux$Geo_Location, " | ",
                                                  kobux$Year)

kobux$Accession_Number <- paste(kobux$Accession_Number, " | ", 
                         kobux$Species, " | ",
                         kobux$host, " | ",
                         kobux$Geo_Location, " | ",
                         kobux$Year)


kobux$Accession_Number <- kobux$new_label


#kobux$Accession_Number <- gsub("NA", "", kobux$Accession_Number)
#kobux$Accession_Number <- gsub("| |", "", kobux$Accession_Number)


final.rooted.kobu$tip.label <- kobux$Accession_Number

final.rooted.kobu$tip.label



# visualize final tree

ggtree(final.rooted.kobu) %<+% kobux + geom_tippoint(aes(color=host, fill=host, shape=bat_host)) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) +
  scale_color_manual(values=colz) + 
  scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(0,70) + ggtitle("Full Genome Kobuvirus Tree") + theme_ap(family="") + 
  theme(legend.position = "bottom", legend.title = ) + geom_treescale()

Fig2B <- ggtree(final.rooted.kobu) %<+% kobux + 
  geom_tippoint(aes(color=host, fill=host, shape=bat_host), size = 1) +
  geom_nodelab(size=1.5,nudge_x = -.5, nudge_y = .7, hjust = 0.5) +
  #scale_color_brewer("Family", palette = "Spectral") + 
  # scale_color_manual(values=colznc) + 
  #scale_fill_manual(values=colznc) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, 
               alpha=.1, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(c(0,40)) + ggtitle("Full Genome Kobuvirus Tree") + 
  theme(legend.position = "bottom", legend.title = ) + 
  geom_treescale(y = -2.5, fontsize = 5, offset = 1, color = "black") + 
geom_cladelab(node = 107, label = "Aichivirus A", align = T, offset = .001) + 
  geom_cladelab(node = 161, label = "Aichivirus B", align = T, offset = .001) + 
  geom_cladelab(node = 111, label = "Aichivirus C", align = T, offset = .001) +
  geom_cladelab(node = 173, label = "Aichivirus D", align = T, offset = .001) +
  geom_cladelab(node = 175, label = "Aichivirus ?", align = T, offset = .001)


# now export 

homewd= '/Users/flg9/Desktop/final-figures/'

setwd(paste0(homewd))

ggsave(file = paste0(homewd, "Fig2B.png"),
       units="mm",  
       width=150, 
       height=100, 
       #limitsize = F,
       scale=4)#, 


########
p1 <- 
  
  ggtree(final.rooted.kobu) %<+% kobux + geom_tippoint(aes(fill= Host, shape=bat_host), shape=21) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) + scale_fill_manual(values = colznc) + scale_shape_manual(values=shapezk) +
  new_scale_fill() + geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) + 
  scale_fill_manual(values = colz2k) + geom_tiplab(aes(fill=novel), geom="label",label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  new_scale_fill() + ggtitle("Kobuvirus Full Genome") + theme(legend.position = "bottom", legend.title = ) 



hostlegend <- get_legend(p1 + theme(legend.position = "bottom"))

p2 <- ggtree(final.rooted.kobu) %<+% kobux + geom_tippoint(aes(shape=bat_host)) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) + scale_fill_manual(values = colznc) + scale_shape_manual(values=shapezk) +
  new_scale_fill() + geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) + 
  scale_fill_manual(values = colz2k) + geom_tiplab(aes(fill=novel), geom="label",label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  new_scale_fill() + ggtitle("Kobuvirus Phylogenetic Tree") + theme(legend.position = "bottom", legend.title =) 
 

batshapelegend <- get_legend(p2 + theme(legend.position = "bottom"))

ggtree(final.rooted.kobu) %<+% kobux + geom_tippoint(aes(fill=new_class,shape=bat_host), show.legend = F) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) + scale_fill_manual(values = colznc) + scale_shape_manual(values=shapezk) +
  new_scale_fill() + geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) + 
  scale_fill_manual(values = colz2k) + geom_tiplab(aes(fill=novel), geom="label",label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  new_scale_fill() + ggtitle("Kobuvirus Phylogenetic Tree") + cowplot::get_legend(hostlegend) + cowplot::get_legend(batshapelegend) +
  theme(legend.position = "none", legend.title =)

ggtree(final.rooted.kobu) %<+% kobux + geom_tippoint(aes(fill=host,shape=bat_host), show.legend = T) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) + scale_fill_manual(values = colznc) + scale_shape_manual(values=shapezk) +
  new_scale_fill() + geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) + 
  scale_fill_manual(values = colz2k) + geom_tiplab(aes(fill=novel), geom="label",label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  new_scale_fill() + ggtitle("Kobuvirus Phylogenetic Tree") + theme(legend.position = "bottom")


+ get_legend(hostlegend +theme(legend.position = "bottom")) + 
  get_legend(batshapelegend + theme(legend.position = "bottom")) 
  




scale_color_manual(name="Host", 
                   labels = c("Human Kobuvirus", 
                              "Bovine Kobuvirus", 
                              "Porcine Kobuvirus", 
                              "Ovine Kobuvirus", 
                              "Canine Kobuvirus", 
                              "Rodent Kobuvirus", 
                              "Sewage Kobuvirus", 
                              "Caprine Kobuvirus", 
                              "Feline Kobuvirus", 
                              "Avian Kobuvirus", 
                              "Rabbit Kobuvirus", 
                              "Bat Kobuvirus"), 
                   values = c("Human Kobuvirus"="red", 
                              "Bovine Kobuvirus"="blue", 
                              "Porcine Kobuvirus"="plum", 
                              "Ovine Kobuvirus"="tomato", 
                              "Canine Kobuvirus"="yellow", 
                              "Rodent Kobuvirus"="pink", 
                              "Sewage Kobuvirus"="darkolivegreen1", 
                              "Caprine Kobuvirus"="green",
                              "Feline Kobuvirus"="purple",
                              "Avian Kobuvirus"="darksalmon", 
                              "Rabbit Kobuvirus"="khaki", 
                              "Bat Kobuvirus"="cyan"))




