## orf aa tree
rm(list=ls())


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
library(aptheme)

# set wd
homewd= '/Users/flg9/Desktop/Developer/brook_lab/orf_tree/'

setwd(paste0(homewd))

# load the orf tree
final.orf <- read.tree(file = paste0(homewd, "T3.raxml.supportFBP"))

#final.orf <- read.tree(file = paste0(homewd, "orf.supportFBP"))

# change weird underscore after each tip label to match .csv file accessions 
name <- unique(final.orf$tip.label)
name <- str_sub(name,1,nchar(name)-1)

final.orf$tip.label <- name 

# change Node_4 to Genbank name
final.orf$tip.label <- gsub("NODE_4", "OP287812", final.orf$tip.label)

# root it
rooted.orf.tree <- root(final.orf, which(final.orf$tip.label=="NC_026314"))


# take a quick look in base R
ggtree(rooted.orf.tree) + 
  geom_nodelab(aes(label=label), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align= FALSE, linetype="dotted", linesize = 0.1, size = 2.5) + geom_point(colour='red')


# load tree data prepared from elsewhere
orf.manual <- read.csv(file=paste0(homewd, 'kobu_orf.csv'), header=T, stringsAsFactors = F)

# color tips by host Family
unique(orf.manual$Family)

colznc = c("Bovine" = "red", "Ovine" = "blue", "Rabbit" = "plum",
           "Canine" = "pink", "Hyaenidae" = "brown", "Caprine" = "darkolivegreen1", "Feline" = "green",
           "Rodent" = "purple", "Human" = "darksalmon", "Avian" = "khaki", "Bat" = "cyan", "Porcine" = "tomato")

C =c("Human","Bovine","Porcine","Ovine",
     "Canine","Rodent","Sewage","Caprine",
     "Feline","Avian","Rabovirus", "Rabbit",
     "Bat", "Bat")


# and add a "novel" category
orf.manual$novel = 0
orf.manual$novel[orf.manual$Accession=="OP287812"] <- 1
orf.manual$novel <- as.factor(orf.manual$novel)


# now add a bat category
orf.manual$bat_host <- 0

orf.manual$bat_host[orf.manual$Accession=="OP287812"] <- "bat-host"
orf.manual$bat_host[orf.manual$Accession=="MN602325"] <- "bat-host"
orf.manual$bat_host[orf.manual$Accession=="MF352427"] <- "bat-host"


# now differentiate between bat-hosts and non-bat-hosts
orf.manual$bat_host[orf.manual$bat_host==0] <- "non-bat-host"
orf.manual$bat_host[orf.manual$bat_host==1] <- "bat-host"

# finally add shapes and colors to bat-hosts/novel viruses
orf.manual$bat_host <- as.factor(orf.manual$bat_host)
shapezk = c("bat-host" =  24, "non-bat-host" = 21)
colz2k = c('1' =  "yellow", '0' = "white")


# label clades

AichiC <- getMRCA(rooted.orf.tree, c("NC_023422", "KY234500"))
AichiB <- getMRCA(rooted.orf.tree, c("MN336260", "GU245693"))
AichiD <- getMRCA(rooted.orf.tree, c("NC_027918", "NC_027919"))
AichiA <- getMRCA(rooted.orf.tree, c("OP287812", "KJ934637"))
AichiEF <- getMRCA(rooted.orf.tree, c("MN602325", "KT325852"))




# create df using the imported tree tip order 
#kobu.dat.orf <- data.frame(new_label=rooted.orf.tree$tip.label, num =1:length(rooted.orf.tree$tip.label))
#orf.manual$new_label <- orf.manual$Accession


# now match with manual .csv by accession
#kobux.orf <- orf.manual[match(kobu.dat.orf$new_label, orf.manual$new_label),]

#kobux.orf$Accession <- kobux.orf$new_label
#rooted.orf.tree$tip.label <- kobux.orf$Accession


# tip labels on manual csv and raxml tree do not match, so fix that 
# first make a new df with the labels numbered by root tip on the raxml tree

kobu.dat <- data.frame(Accession=final.orf$tip.label, 
                       num =1:length(final.orf$tip.label))

# then create a 3rd df and right join the original manual and the 2nd df 
# by accession_number onto this 3rd df 

kobux.orf <- join(orf.manual, kobu.dat, by = "Accession", match = "all", 
              type = "right")

# then create new tip labels 

kobux.orf$new_label <- kobux.orf$Accession

kobux.orf$new_label[!is.na(kobux.orf$new_label)] <- paste(kobux.orf$Accession, " | ", 
                                                  kobux.orf$Species[!is.na(kobux.orf$Species)], " | ",
                                                  kobux.orf$Family[!is.na(kobux.orf$Family)], " | ",
                                                  kobux.orf$Host[!is.na(kobux.orf$Host)], " | ",
                                                  kobux.orf$Geo_Location[!is.na(kobux.orf$Geo_Location)], " | ",
                                                  kobux.orf$Year)



kobux.orf$Accession <- kobux.orf$new_label

kobux.orf$new_label<- str_replace_all(kobux.orf$new_label, "NA", "")

rooted.orf.tree$tip.label <- kobux.orf$Accession

# now visualize 

ggtree(rooted.orf.tree) %<+% kobux.orf + geom_tippoint(aes(color=Family, fill=Family, shape=bat_host)) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) +
  #scale_color_manual(values=colznc) + 
  #scale_fill_manual(values=colznc) +
  scale_shape_manual(values=shapezk) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2k) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(c(0,25)) + ggtitle("Kobuvirus ORF Tree") + 
  theme(legend.position = "bottom", legend.title = ) + geom_treescale() + 
  geom_cladelab(node = 149, label = "Aichivirus A", align = T, offset = .4) + 
  geom_cladelab(node = 133, label = "Aichivirus B", align =T, offset = .4) + 
  geom_cladelab(node = 131, label = "Aichivirus C", align = T, offset = .4) +
  geom_cladelab(node = 177, label = "Aichivirus D", align = T, offset = .4) +
  geom_cladelab(node = 147, label = "Aichivirus ?", align = T, offset = .4)
  
ggtree(rooted.orf.tree) %<+% kobux.orf + geom_tippoint(aes(color=Family, fill=Family, shape=bat_host)) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) +
  #scale_color_brewer("Family", palette = "Spectral") + 
 #scale_color_manual(values=colznc) + 
#scale_fill_manual(values=colznc) +
  scale_shape_manual(values=shapezk) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2k) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(c(0,25)) + ggtitle("Kobuvirus ORF Tree") + 
  theme(legend.position = "bottom", legend.title = ) + 
  geom_treescale(y = -2) + theme_ap(family = "")

Fig2C <- ggtree(rooted.orf.tree) %<+% kobux.orf +
  geom_tippoint(aes(color=Family, fill=Family, shape=bat_host), size = 1) +
  geom_nodelab(size= 1.5, nudge_x = -.12, nudge_y = .7, geom = "text") +
  scale_shape_manual(values = shapezk) +
  new_scale_fill() + 
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha = .3, size = 1.8, show.legend = F) + 
  scale_fill_manual(values = colz2k) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) + 
  xlim(c(0,15)) + ggtitle("Kobuvirus ORF Tree") + 
  theme(legend.position = "bottom", legend.title = ) + 
  geom_treescale(y = -2.5, fontsize = 5, offset = 1, color = "black") + 
geom_cladelab(node = 149, label = "Aichivirus A", align = T, offset = .5) + 
  geom_cladelab(node = 133, label = "Aichivirus B", align = T, offset = .5) + 
  geom_cladelab(node = 131, label = "Aichivirus C", align = T, offset = .5) +
  geom_cladelab(node = 177, label = "Aichivirus D", align = T, offset = .5) +
  geom_cladelab(node = 147, label = "Aichivirus ?", align = T, offset = .5)

# now export 

homewd= '/Users/flg9/Desktop/final-figures/'

setwd(paste0(homewd))

ggsave(file = paste0(homewd, "Fig2C.png"),
       units="mm",  
       width=150, 
       height=100, 
       #limitsize = F,
       scale=4)#, 







