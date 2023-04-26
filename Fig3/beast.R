# now make beast tree 

rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(lubridate)
library(treeio)
library(rBt)

# make Bayesian time tree from Kobuvirus relaxed molecular clock model
# first, read in the tree

homewd= "/Users/flg9/Desktop/Developer/brook_lab/beast_tree/extra_seqs/"
setwd(paste0(homewd))

beast_tree <- read.beast(file = paste0(homewd, "bat_kobu.MCC.tree"))

# initial visualize 
ggtree(beast_tree) + 
  geom_nodelab(aes(label=label), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align= FALSE, linetype="dotted", linesize = 0.1, size = 2.5) + geom_point(colour='red')

# create df that contains accession numbers from beast tree 
# use cbind to do this, which pulls those accessions from tip label data 
# then create new column to mess with 

treedat <- cbind.data.frame(tip_name = beast_tree@phylo$tip.label)
treedat$beast_name <-treedat$tip_name

#tree <- read.annot.beast(file = paste0(homewd, "/Fig4/beast-out/AllNobeco/NobecoStrict/AvgNobecoStrictNexus.trees"))

#tree$node.label <- round(tree$posterior,2)
#treedat <- cbind.data.frame(tip_name = tree$tip.label)

treedat$accession_num <- sapply(strsplit(treedat$tip_name, "_"), function(x) x[[1]])
treedat$accession_num[treedat$accession_num=="NC"] <- c("NC_001918", "NC_011829", 
                                                        "NC_023422", "NC_027919")
#names(treedat)[names(treedat)=="tip_name"] <- "beast_name"

# load the corresponding manual tree csv data 


dat <- read.csv(file = "beast_kobuvirus_metadata_manual copy.csv", 
                header = T, stringsAsFactors = F)

# and omit any na 
dat<-na.omit(dat)

# change collection date to something more readable 

dat$collection_date <- as.Date(dat$collection_date)

# test out tree 

mrsd.dat <- max(dat$collection_date)
p1 <- ggtree(beast_tree, mrsd=mrsd.dat)  + theme_tree2()  +geom_nodelab()

# and extract tree data 

tree.dat <- p1$data

# extract node data and change labels 

node.sub <- dplyr::select(tree.dat, node, x)
names(node.sub) <-  c("node", "nodetime")

# head tree data and make new clade column  

head(dat)

dat$clade <- dat$type

# don't include?  

#dat$clade[is.na(dat$clade)] <- "African Eidolon"
dat$clade[dat$host == "Pteropus_rufus" ] <- "Madagascar Pteropus" 
dat$clade[dat$clade == "GX2018"] <- "BtCoV92 / GX2018" 
dat$clade[dat$clade == "BtCoV92"] <- "BtCoV92 / GX2018" 
#dat$clade[dat$accession_num == "KU182962"] <- "BtCoV92 / GX2018" 
#dat$clade[dat$host == "Eidolon_helvum" | dat$host == "Rousettus_madagascariensis"] <- "African Eidolon" 
#dat$clade[dat$accession_num=="MG693170"] <- "HKU9"
#dat$clade[dat$host == "Pteropus_rufus" ] <- "Madagascar Pteropus" 

# change col name in tree.dat from label to accession_num
# and merge treedat and dat by accession number 

# colnames takes df and changes name by column position within df 

colnames(tree.dat)[4] <- "accession_num"

dat.plot <- merge(treedat, dat, by="accession_num", all.x = T, sort=F)

# now make new labels for tree tips using paste function

head(dat.plot)
dat.plot$new_label = NA

dat.plot$new_label[!is.na(dat.plot$type)] <- paste(dat.plot$accession_num[!is.na(dat.plot$type)], " | ", 
                                                     dat.plot$type[!is.na(dat.plot$type)], " | ", 
                                                     dat.plot$host[!is.na(dat.plot$type)], " | ",
                                                     dat.plot$country[!is.na(dat.plot$type)], " | ",
                                                     dat.plot$collection_year[!is.na(dat.plot$type)])

#dat.plot$new_label[is.na(dat.plot$type)] <- paste(dat.plot$accession_num[is.na(dat.plot$type)], " | ", 
 #                                                   dat.plot$host[is.na(dat.plot$type)], " | ",
 #                                                   dat.plot$country[is.na(dat.plot$type)], " | ",
 #                                                   dat.plot$collection_year[is.na(dat.plot$type)])
# check to make sure labels are fine 

dat.plot$new_label

# hard code tip label on original tree file to be new label on dat.plot

beast_tree@phylo$tip.label <- dat.plot$new_label

# create new df that only contains columns we need for final tree

dat.sub <- dplyr::select(dat.plot, new_label, collection_date, country, clade)
head(dat.sub)

# change dat.sub$clade from character to factor

dat.sub$clade <- as.factor(dat.sub$clade)

p2 <-ggtree(beast_tree, mrsd=mrsd.dat) %<+% dat.sub + geom_tippoint(aes(color=clade)) +
  geom_tiplab(size=3) + geom_nodelab(size=2,nudge_x = -15, nudge_y = .7) +
  theme_tree2() +
  theme(legend.position = c(.1,.85),
        plot.margin = unit(c(2,20,2,3), "lines")) +
  coord_cartesian(clip = "off")
p2


# second node label that is the date for E. dupreanum and the original?

#orig.date <- round(node.sub$nodetime[27],0)

### left off here

nodeEidol <- MRCA(beast_tree, which(beast_tree@phylo$tip.label == "OP287812  |  unclassifed  |  eidolon_dupreanum  |  Madagascar  |  2018"),
                 which(beast_tree@phylo$tip.label == "NC_001918  |  AiVA  |  homo_sapien  |  Japan  |  1989"))

nodeRous <- MRCA(tree, which(tree@phylo$tip.label == "OK067320  |  Rousettus_madagascariensis  |  Madagascar  |  2018"),
                 which(tree@phylo$tip.label == "HM211099  |  HKU9  |  Rousettus_leschenaulti  |  China  |  2005"))

#nodeall <- MRCA(tree, which(tree$tip.label == "KP696747  |  Pteropus_rufus  |  Madagascar  |  2011"),which(tree$tip.label == "MK211379  |  GX2018  |  Rhinolophus_affinis  |  China  |  2016"))

#orig.date <- round(node.sub$nodetime[nodeall],0)
Pruf.date <- round(node.sub$nodetime[nodePruf],0) #date that P ruf branches off from everything else
recent.age <- year(mrsd.dat) + yday(mrsd.dat)/365
Pruf.mean <- round(recent.age-tree@data$height[35],0) #sane as above--date that P. ruf branches off -

#this is taking the P. ruf tip and subtracting the branch length from present to get its position
Pruf.uci <- round(recent.age-tree@data$height_0.95_HPD[35][[1]][1],0)
Pruf.lci <- round(recent.age-tree@data$height_0.95_HPD[35][[1]][2],0)

#should be close to the lengths: tree@data$length_0.95_HPD[35] (does not have it)
Pruf.date <- paste0("~", Pruf.mean, "\n[", Pruf.lci, "-", Pruf.uci, "]")#from FigTree

Rous.mean <- round(recent.age-tree@data$height[33],0)
Rous.date <- round(node.sub$nodetime[nodeRous],0)#they are the same--we're good
Rous.uci <- round(recent.age-tree@data$height_0.95_HPD[33][[1]][1],0)
Rous.lci <- round(recent.age-tree@data$height_0.95_HPD[33][[1]][2],0)
Rous.date <- paste0("~", Rous.mean, "\n[", Rous.lci, "-", Rous.uci, "]")#from FigTree

new.nodel.lab <- rep(NA, nrow(node.sub))
#new.nodel.lab[nodeall] <- paste0("~",orig.date)
new.nodel.lab[nodePruf] <- Pruf.date
new.nodel.lab[nodeRous] <- Rous.date

dat.sub$clade <- as.character(dat.sub$clade)
dat.sub$clade[dat.sub$clade=="African Eidolon"] <- "African~italic(Eidolon)"
dat.sub$clade[dat.sub$clade=="Madagascar Pteropus"] <- "Madagascar~italic(Pteropus)"

# make novel category to highlight mada kobuvirus using country 

dat.sub$novel = "no"
dat.sub$novel[dat.sub$country=="Madagascar"] <- "yes"

colz2 = c('yes' =  "yellow", 'no' = "white")

ggtree(beast_tree, mrsd = mrsd.dat) %<+% dat.sub + geom_tippoint(aes(color = clade), size = 3) + 
  geom_nodelab(size = 2.5, nudge_x = -21, nudge_y = .7) +
  geom_nodelab(size=4,nudge_x = -55, nudge_y = -1,  color="firebrick", fontface=2, geom="label", fill="white") + 
  theme_tree2() + 
  theme(legend.position = c(.2,.75), 
        plot.margin = unit(c(.2,20,2,3), "lines")) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=2, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red") + 
  scale_x_continuous(breaks=c(1700, 1800, 1900, 2000)) + ggnewscale::new_scale_fill() + 
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=4, hjust = -.05) + scale_fill_manual(values=colz2) 
  
p3 <-ggtree(beast_tree, mrsd=mrsd.dat) %<+% dat.sub + geom_tippoint(aes(color=clade), size=3) +
  #geom_tiplab(size=3, nudge_x=5) + 
  #geom_nodelab(size=2.5,nudge_x = -21, nudge_y = .7) +
  geom_nodelab(size=4,nudge_x = -21, nudge_y = -1,  color="firebrick", fontface=2, geom="label", fill="white") +
  theme_tree2() +
  #geom_treescale(fontsize=3, x=1300,y=22, linesize = .5, width=200,label="years") + 
  #scale_color_discrete(labels=c(parse(text="African~italic(Eidolon)"), "BtCoV92 / GX2018", "HKU9", parse(text="Madagascar~italic(Pteropus)"))) +
  theme(legend.position = c(.2,.75), 
        plot.margin = unit(c(.2,20,2,3), "lines")) +
  coord_cartesian(clip = "off", xlim=c(1600, 2150)) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=2, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red") +
  guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1)) +
  scale_x_continuous(breaks=c(1700, 1800, 1900, 2000))+ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=4, hjust = -.05) + scale_fill_manual(values=colz2)

p3 <-ggtree(beast_tree, mrsd=mrsd.dat) %<+% dat.sub + geom_tippoint(aes(color=clade), size=3) +
  #geom_tiplab(size=3, nudge_x=5) + 
  geom_nodelab(size=2.5,nudge_x = -21, nudge_y = .7) +
  geom_nodelab(aes(label=new.nodel.lab), size=4,nudge_x = -55, nudge_y = -1,  color="firebrick", fontface=2, geom="label", fill="white") +
  theme_tree2() +
  #geom_treescale(fontsize=3, x=1300,y=22, linesize = .5, width=200,label="years") + 
  scale_color_discrete(labels=c(parse(text="African~italic(Eidolon)"), "BtCoV92 / GX2018", "HKU9", parse(text="Madagascar~italic(Pteropus)"))) +
  theme(legend.position = c(.2,.75), 
        plot.margin = unit(c(.2,20,2,3), "lines")) +
  coord_cartesian(clip = "off", xlim=c(1600, 2150)) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=2, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red") +
  guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1)) +
  scale_x_continuous(breaks=c(1700, 1800, 1900, 2000))+ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=4, hjust = -.05) + scale_fill_manual(values=colz2)

ggtree(beast_tree, mrsd=mrsd.dat) %<+% dat.sub + geom_tippoint(aes(color=clade), size=3) +
  geom_tiplab(size=3, nudge_x=5) + 
  geom_nodelab(size=2.5,nudge_x = -21, nudge_y = .7) +
  geom_nodelab(size=4,nudge_x = -55, nudge_y = -1,  color="firebrick", fontface=2, geom="label", fill="white") +
  theme_tree2() +
  #geom_treescale(fontsize=3, x=1300,y=22, linesize = .5, width=200,label="years") + 
  #scale_color_discrete(labels=c(parse(text="African~italic(Eidolon)"), "BtCoV92 / GX2018", "HKU9", parse(text="Madagascar~italic(Pteropus)"))) +
  theme(legend.position = c(.2,.75), 
        plot.margin = unit(c(.2,20,2,3), "lines")) +
  coord_cartesian(clip = "off") + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=2, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red") +
  guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1)) +
  ggnewscale::new_scale_fill() + 
  scale_fill_manual(values=colz2) + geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
                                                alpha=.3,  show.legend=F, size=4, hjust = -.05) + 
scale_x_continuous(breaks=c(1700, 1800, 1900, 2000)) 

p3 <-ggtree(beast_tree, mrsd=mrsd.dat, size=.8) %<+% dat.sub +
  geom_tippoint(aes(color=clade), size=4) +
  #scale_color_manual(values=c('Bat AstV' = "#0CB702", 'Bovine AstV' = "#00BFC4", 'Feline AstV' = "#C77CFF", 'Human AstV' = "#ED68ED", 'Mink AstV' = "#CD9600", 'Murine AstV' = "#ABA300", 'Porcine AstV' = "#7CAE00")) +
  #scale_fill_manual(values=c('Bat AstV' = "#0CB702", 'Bovine AstV' = "#00BFC4", 'Feline AstV' = "#C77CFF", 'Human AstV' = "#ED68ED", 'Mink AstV' = "#CD9600", 'Murine AstV' = "#ABA300", 'Porcine AstV' = "#7CAE00")) +
  #geom_nodelab(size=2.5,nudge_x = -21, nudge_y = .7) +
  # geom_nodelab(aes(label=new.nodel.lab), size=4,nudge_x = -55, nudge_y = -1,  color="firebrick", fontface=2, geom="label", fill="white") +
  #geom_nodelab(size=1.8,nudge_x = -.05, nudge_y = .7) +
  #geom_treescale(fontsize=2.5) + 
  theme_tree2() +
  theme(legend.position = c(.07,.7), plot.margin = unit(c(.2,24,3,3), "lines"), legend.title = element_blank()) +
  coord_cartesian(clip = "off") + 
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=3, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red") +
  guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1)) +
  scale_x_continuous(breaks=c(-28000, -23000, -18000, -13000, -8000, -3000, 2000),
                     labels=c(30000, 25000, 20000, 15000, 10000, 5000,  0)) +
  xlab("years to MRCA") +
  ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=3) + scale_fill_manual(values=colz2) 
#xlim(c(0,6))

# now save
homewd="/Users/flg9/Desktop/"
setwd(paste0(homewd))

ggsave(file = paste0(homewd, "/final-figures/BEAST.png"),
       units="mm",  
       width=90, 
       height=60, 
       #limitsize = F,
       scale=3)#, 

