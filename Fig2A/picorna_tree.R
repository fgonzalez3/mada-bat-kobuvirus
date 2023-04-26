rm(list=ls())

#time to make Fig 1a

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
library(caper)
library(dplyr)


# set wd
homewd= '/Users/flg9/Desktop/Developer/brook_lab/picorna_tree/'

setwd(paste0(homewd))

# load the Fig 1a tree

final.picorna <-  read.tree(file = paste0(homewd, "T3.raxml.supportFBP"))

# change Node_4 to Genbank name
final.picorna$tip.label <- gsub("NODE_4", "OP287812", final.picorna$tip.label)

# root it using cov

final.rooted.picorna <- root(final.picorna, which(final.picorna$tip.label=='NC_048212'))


# take a quick look
ggtree(final.rooted.picorna) + 
  geom_nodelab(aes(label=label), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align= FALSE, linetype="dotted", linesize = 0.1, size = 2.5) + geom_point(colour='red')


# load tree data for picorna nt tree
# and fill in for any NA values

picorna.manual <- read.csv(file=paste0(homewd,'picorna_manual.csv'), 
                           header=T, stringsAsFactors = F, na = "")

picorna.manual[is.na(picorna.manual)] = "NA"

# use dplyr to only include columns that will be in the final tip label
# mass package also loaded, so need to explicitly call dplyr 

colnames(picorna.manual)

picorna.manual <- picorna.manual %>%
  dplyr::select(Accession_Number, Species, Geo_Location, Family, Year)

#check unique hosts that will be used to color tip labels 

unique(picorna.manual$Family)

colz = c("Human" = "red", "Bovine" = "tomato", "Porcine" = "khaki", 
         "Canine " = "blue", "Rodent" = "plum", "Caprine" = "darkolivegreen1", 
         "Feline" = "purple", "Avian" = "darksalmon", "Bat" = "cyan", 
         "Ursid" = "cadetblue1", "Reptile" = "mediumvioletred", 
         "Erinaceidae" = "darkseagreen1", "Amphibian" = "yellow", 
         "Primate" = "pink", "Camelid" = "brown", "Equine" = "darkgreen", 
         "Fish" = "thistle", "Shrew" = "sienna", "Seal" = "wheat",
         "Marsupial" = "peru", "NA" = "black")

#pick order for the labels

C = c("Human", "Bovine", "Porcine", "Canine", "Rodent", "Caprine", "Feline", 
      "Avian", "Bat", "Ursid", "Reptile", "Erinaceidae", "Amphibian", "Primate", 
      "Camelid", "Equine", "Fish", "Shrew", "Seal", "Marsupial", "NA")

                      
#and add a "novel" category
picorna.manual$novel = 0
picorna.manual$novel[picorna.manual$Accession_Number=="OP287812"] <- 1

picorna.manual$novel <- as.factor(picorna.manual$novel)

#tip shapes for bat hosts (picorna) 
picorna.manual$bat_host <- 0

picorna.manual$bat_host[picorna.manual$Accession_Number=="OP287812"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_038313"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_038316"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_038961"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_034381"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_033820"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_030843"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_028366"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_026470"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_015934"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_015940"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_015941"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="MN602325"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="MF352419"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="MF352423"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="MF352427"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="KJ641686"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="KJ641687"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="KJ641691"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="KJ641693"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="KJ641694"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="KJ641696"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="HQ595341"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="HQ595343"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="HQ595345"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_048212"] <- "bat-host"


picorna.manual$bat_host[picorna.manual$bat_host==0] <- "non-bat-host"
picorna.manual$bat_host[picorna.manual$bat_host==1] <- "bat-host"

picorna.manual$bat_host <- as.factor(picorna.manual$bat_host)
shapez = c("bat-host" =  24, "non-bat-host" = 21)
colz2 = c('1' =  "yellow", '0' = "white")


# visualize again using these labels 

ggtree(final.rooted.picorna) %<+% picorna.manual + geom_tippoint(aes(fill=Family)) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  #scale_fill_manual(values=colz) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank())


# tip labels on manual csv and raxml tree do not match, so fix that 
# first make a new df with the labels numbered by root tip on the raxml tree

picorna.dat<- data.frame(Accession_Number=final.rooted.picorna$tip.label, 
                      num =1:length(final.rooted.picorna$tip.label))

# then create a 3rd df and right join the original manual and the 2nd df 
# by accession_number onto this 3rd df 

picornax <- join(picorna.manual, picorna.dat, by = "Accession_Number", match = "all", 
              type = "right")

# now begin to make a new label 

picornax$new_label <- picornax$Accession_Number


picorna.manual$tip.label <- NA


# change tip names


picornax$new_label[!is.na(picornax$new_label)] <- paste(picornax$Accession, " | ", 
                                                  picornax$Genus[!is.na(picornax$Genus)], " | ",
                                                  picornax$Family[!is.na(picornax$Family)], " | ",
                                                  picornax$Host[!is.na(picornax$Host)], " | ",
                                                  picornax$Geo_Location[!is.na(picornax$Geo_Location)], " | ",
                                                  picornax$Year)

picornax$new_label <- paste(picornax$Accession_Number, " | ", 
                            picornax$Species, " | ",
                            picornax$Family, " | ",
                            picornax$Geo_Location, " | ",
                            picornax$Year)

picornax$Accession_Number <- picornax$new_label


final.rooted.picorna$tip.label <- picornax$Accession_Number

final.rooted.picorna$tip.label

#picornax$new_label<- str_replace_all(picornax$new_label, "NA", "")
#picornax$new_label<- str_replace_all(picornax$new_label, "| ", "")

# visualize final tree

fig2a <- ggtree(final.rooted.picorna) %<+% picornax + 
  geom_tippoint(aes(color=Family, fill=Family, shape=bat_host), size=1) +
  geom_nodelab(size=2, nudge_x = -.1, nudge_y = .7) +
  #scale_color_manual(values=colz) + 
  #scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.1, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(c(0,10)) + ggtitle("Full Genome Picornavirus Tree") + 
  theme(legend.position = "bottom", legend.title = ) + 
  geom_treescale(y = -6, fontsize = 5, offset = 1, color = "black")


ggtree(final.rooted.picorna) %<+% picornax + 
  geom_tippoint(aes(color = Family, fill = Family, shape = bat_host), size = 1) + 
  geom_nodelab(size = 2, nudge_x = -.1, nudge_y = .7) + 
  scale_shape_manual(values = shapez) + 
  new_scale_fill() + 
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0,
              alpha = .1, size = 1.8, show.legend = F) + xlim(c(0,20))


 ggtree(rooted.tree.B) %<+% datB + 
  geom_tippoint(aes(fill=sub_group, shape=bat_host), show.legend = F, size=5) +
  geom_nodelab(size=3,nudge_x = -.03, nudge_y = .4) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=4, hjust = -.1) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "right", legend.title = element_blank()) +
  geom_treescale(fontsize=4, x=.3,y=50, linesize = .5) + 
  geom_cladelabel(node = clade.a, label = "HKU9", offset = 1, fontsize = 6.5, color="tomato", offset.text = .01) +
  geom_cladelabel(node = clade.b, label = "atop(African,italic(Eidolon))", offset = .8, fontsize = 6.5, color="tomato", parse=T) +
  geom_cladelabel(node = clade.c, label = "GCCDC1", offset = 1.05, fontsize = 6.5, color="tomato", offset.text = .01) +
  geom_cladelabel(node = clade.d, label = "BtCoV92 /\nGX2018", offset = 1.03, fontsize = 6.5, color="tomato", offset.text = .01) +
  geom_cladelabel(node = clade.e, label = "atop(Madagascar,italic(Pteropus))", offset = .7, fontsize = 6.5, color="tomato" , parse = T) +
  xlim(c(0,1.8))


 homewd="/Users/flg9/Desktop/"
 setwd(paste0(homewd))
 
 ggsave(file = paste0(homewd, "final-figures/Fig2a.png"),
        plot = fig2a,
        units="mm",  
        width=150, 
        height=100, 
        #limitsize = F,
        scale=4)#,
 


#after Gwen checks these, can later manually edit any that are "NA"

#make sure to sort in order
tree.datB <- data.frame(old_tip_label=rooted.tree.B$tip.label, num =1:length(rooted.tree.B$tip.label))
tree.datB <- merge(tree.datB, datB, by = "old_tip_label", all.x = T, sort = F)
rooted.tree.B$tip.label <- tree.datB$tip_label

datB$bat_host[datB$bat_host==0] <- "non-bat host"
datB$bat_host[datB$bat_host==1] <- "bat host"
datB$bat_host <- as.factor(datB$bat_host)

datB$novel = 0
datB$novel[datB$accession_num=="OK020086" |
             datB$accession_num=="OK067321" |
             datB$accession_num=="OK067320" | 
             datB$accession_num=="OK020089" |
             datB$accession_num=="OK067319" |
             datB$accession_num=="OK020087" |
             datB$accession_num=="OK020088" ] <- 1

datB$novel <- as.factor(datB$novel)
colz2 = c('1' =  "yellow", '0' = "white")

shapez = c("bat host" =  24, "non-bat host" = 21)


  ggtree(rooted.tree.B) %<+% datB + geom_tippoint(aes(fill=sub_group, shape=bat_host)) +
  geom_nodelab(size=1.5,nudge_x = -.02, nudge_y = .4) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) +
  new_scale_fill()+
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "right", legend.title = element_blank()) +
  xlim(c(0,2))
p2














#####now combine the two together somehow


###wrking on p1

p1 <- ggtree(rooted.tree.A) %<+% tree.dat + 
  geom_tippoint(aes(fill=sub_group, shape=bat_host), size=3, show.legend = F) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = .9) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3,  show.legend=F, size=3, hjust = -.08) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  geom_treescale(fontsize=4, x=1,y=124, linesize = .5) + 
  xlim(c(0,4))
p1

#and flip some clades
node_flip_Embeco_Merbeco1 = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "NC_019843  |  MERS  |  Homo_sapiens  |  Saudi_Arabia  |  2012" ),which(rooted.tree.A$tip.label == "NC_006213  |  HCoV_OC43  |  Homo_sapiens  |  USA  |  1960"  ))
node_flip_Merbeco_Sarbeco1 = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "NC_019843  |  MERS  |  Homo_sapiens  |  Saudi_Arabia  |  2012" ),which(rooted.tree.A$tip.label == "MZ081380  |  SARSr_CoV  |  Rhinolophus_stheno  |  China  |  2020"))
node_flip_Embeco_Sarbeco1 = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "MZ081380  |  SARSr_CoV  |  Rhinolophus_stheno  |  China  |  2020" ),which(rooted.tree.A$tip.label == "NC_006213  |  HCoV_OC43  |  Homo_sapiens  |  USA  |  1960"  ))

p1.2 <- p1 %>% ggtree::rotate(node = node_flip_Embeco_Sarbeco1 )


#p1.2 <- p1 %>% ggtree::rotate(node = node_flip_Merbeco_Sarbeco1)
#p1.3 <- p1.2 %>% ggtree::rotate(node = node_flip_Embeco_Sarbeco1)


#collapse the alpha clade (all bat CoVs)
#alpha_node = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "NC_048211 | Suncus_murinus | China | 2015" ),which(rooted.tree.A$tip.label == "NC_018871 | bat | China | 2021" ))

p1.2.leg <- ggtree(rooted.tree.A) %<+% tree.dat + 
  geom_tippoint(aes(color=sub_group, shape=bat_host), size=3) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = .9) +
  scale_fill_manual(values=colz) + 
  scale_color_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3,  show.legend=F, size=3, hjust = -.08) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size=12)) +
  geom_treescale(fontsize=4, x=1,y=124, linesize = .5) + 
  xlim(c(0,4))
p1.2.leg

#separate legend
leg.all <- cowplot::get_legend(p1.2.leg)


#new p2
p2.1 <- ggtree(rooted.tree.B) %<+% datB + 
  geom_tippoint(aes(fill=sub_group, shape=bat_host), show.legend = F, size=5) +
  geom_nodelab(size=3,nudge_x = -.03, nudge_y = .4) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=4, hjust = -.1) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "right", legend.title = element_blank()) +
  geom_treescale(fontsize=4, x=.3,y=50, linesize = .5) + 
  xlim(c(0,1.5))
p2.1 

#add lineage clade labels bars

#nodebase
clade.a <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "EF065514  |  HKU9  |  Rousettus_leschenaulti  |  China  |  2005" ),which(rooted.tree.B$tip.label == "HM211098  |  HKU9  |  Rhinolophus_sinicus  |  China  |  2005" ))
clade.b <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "OK020089  |  Rousettus_madagascariensis  |  Madagascar  |  2018" ),which(rooted.tree.B$tip.label == "MG693172  |  Eidolon_helvum  |  Cameroon  |  2013" ))
clade.c <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_030886  |  GCCDC1  |  Rousettus_leschenaulti  |  China  |  2014" ),which(rooted.tree.B$tip.label == "MT350598  |  GCCDC1  |  Eonycteris_spelaea  |  Singapore  |  2016" ))
clade.d <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "MK211379  |  GX2018  |  Rhinolophus_affinis  |  China  |  2016" ),which(rooted.tree.B$tip.label == "MK492263  |  BatCoV92  |  Cynopteris_brachyotis  |  Singapore  |  2015" ))
clade.e <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "OK020087  |  Pteropus_rufus  |  Madagascar  |  2018" ),which(rooted.tree.B$tip.label == "OK067319  |  Pteropus_rufus  |  Madagascar  |  2018" ))




p2.1 <- ggtree(rooted.tree.B) %<+% datB + 
  geom_tippoint(aes(fill=sub_group, shape=bat_host), show.legend = F, size=5) +
  geom_nodelab(size=3,nudge_x = -.03, nudge_y = .4) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=4, hjust = -.1) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "right", legend.title = element_blank()) +
  geom_treescale(fontsize=4, x=.3,y=50, linesize = .5) + 
  geom_cladelabel(node = clade.a, label = "HKU9", offset = 1, fontsize = 6.5, color="tomato", offset.text = .01) +
  geom_cladelabel(node = clade.b, label = "atop(African,italic(Eidolon))", offset = .8, fontsize = 6.5, color="tomato", parse=T) +
  geom_cladelabel(node = clade.c, label = "GCCDC1", offset = 1.05, fontsize = 6.5, color="tomato", offset.text = .01) +
  geom_cladelabel(node = clade.d, label = "BtCoV92 /\nGX2018", offset = 1.03, fontsize = 6.5, color="tomato", offset.text = .01) +
  geom_cladelabel(node = clade.e, label = "atop(Madagascar,italic(Pteropus))", offset = .7, fontsize = 6.5, color="tomato" , parse = T) +
  xlim(c(0,1.8))
p2.1 


#great, now need to flip some of the clases to match plot on the left

node_flip_Embeco_Merbeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_019843  |  MERS  |  Homo_sapiens  |  Saudi_Arabia  |  2012" ),which(rooted.tree.B$tip.label == "NC_006213  |  HCoV_OC43  |  Homo_sapiens  |  USA  |  1960"  ))
node_flip_Sarbeco_Hibeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_025217  |  Hipposideros_pratti  |  China  |  2013" ),which(rooted.tree.B$tip.label == "NC_004718  |  SARS_CoV  |  Homo_sapiens  |  Canada  |  2003" ))
node_flip_Embeco_Nobeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_019843  |  MERS  |  Homo_sapiens  |  Saudi_Arabia  |  2012" ),which(rooted.tree.B$tip.label == "EF065516  |  HKU9  |  Rousettus_leschenaulti  |  China  |  2005"   ))


p2.2 <- p2.1 %>% ggtree::rotate(node = node_flip_Embeco_Merbeco)
p2.3 <- p2.2 %>% ggtree::rotate(node = node_flip_Sarbeco_Hibeco)
p2.4 <- p2.3 %>% ggtree::rotate(node = node_flip_Embeco_Nobeco)


Fig3 <- cowplot::plot_grid(p1.2,p2.4, ncol=2, nrow=1, labels = c("(A)", "(B)"), label_size = 22, label_x = .03, label_y = .98)

Fig3all <- cowplot::plot_grid(Fig3,leg.all, ncol=1, nrow=2, rel_heights = c(1,.1))


#and save to the final figures

ggsave(file = paste0(homewd, "/final-figures/Fig3.png"),
       units="mm",  
       width=150, 
       height=100, 
       #limitsize = F,
       scale=4)#, 

