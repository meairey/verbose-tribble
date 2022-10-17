## Read in libraries ---------
library(ggplot2)
library(dplyr)
library(gridExtra)
library(mapproj)

library(tidyr)
library(tidyverse)
library(networkD3)
library(usmap)


## Read in data -------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/verbose-tribble")
oto=read.csv("Data/Otolith_Data_CSV.csv") %>%
  filter(Area != "NA")

## PCA of rim and core values 

pca_dat = oto %>% select(Core_O,Rim_O,Core_C, Rim_C, Group.C, Individual) %>% 
  rename("Sample_ID" = Individual) %>%
  mutate(individal = row_number()) %>% 
  pivot_longer(1:4, names_to = "group") %>%
  separate(group, into = c("age", "isotope")) %>%  
  pivot_wider(names_from = isotope, values_from = value) %>%
  unite("group",c(Group.C, age), remove = F) %>%
  na.omit()



results = prcomp(pca_dat %>% select(O,C))


results$x %>% as.data.frame() %>% 
  ggplot(aes(x = PC1, y = PC2, color = pca_dat$group)) +
  geom_text(aes(label = as.character(pca_dat$Sample_ID))) +
  stat_ellipse(level = .95)+
  scale_color_manual(
    values = as.numeric(as.factor(pca_dat$group))
  )

## These are the individuals that do not fit into an ellipse with 97.5% stat ellipses
list.975 = c(4,6,8,19,23,1)
list.975.samp = (pca_dat %>% filter(individal %in% list.975))$Sample_ID %>% unique
list.95 = c(4,6,8,19,23,1,35,10,14)
list.95.samp = (pca_dat %>% filter(individal %in% list.95))$Sample_ID %>% unique

## Graph of d18O for un-classifiable individuals 

pca_dat %>% filter(individal %in% list.95) %>% ggplot(aes(x = C, y = O, color = Group.C, shape = age)) + 
  geom_point()

## Rim vs. Core d18O graph ----
ggplot(data=oto, aes(x = Rim_O, y=Core_O, color=Depth, shape = Area)) + 
  geom_point(size=5)  + 
  scale_color_gradient(low='#91d0f2', high='#39658e') + 
  geom_abline(intercept=0) + 
  ylim(-1.3,1) +
  xlim(-1.3,1) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text=element_text(size=13))

## Rim vs. Core d13C graph ---- 
ggplot(data=oto, aes(x = Rim_C, y=Core_C, color=Depth, shape = Area)) +  
  geom_point(size=5)  + 
  scale_color_gradient(low='#91d0f2', high='#39658e') + 
  geom_abline(intercept=0) + 
  ylim(-7.7, -3) +
  xlim(-7.7,-3) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text=element_text(size=13)) 

## Rim vs Core d18O graph --- 
# maybe don't use this one? 
ggplot(oto, aes(Group.C, Rim_O, fill=Group.C)) +
  geom_boxplot() + 
  scale_fill_brewer(palette = "Blues") + 
  scale_x_discrete(labels=c("A"="0-15",
                            "B"="15-30",
                            "C"="30-45",
                            "D"="45-60",
                            "E"="60-75",
                            "G"="90-105")) +
  xlab("Depth (m)") +
  guides(fill=F) +
  theme(text = element_text(size=26))+
  ylab(expression({delta}^18*O~'\u2030'))

## Depth vs. 018 graph -----
oto %>%
  ggplot(aes(x = Depth, 
             y = Rim_O)) +
  geom_point(cex = 4, aes(pch = Area)) + 
  ylab(({delta}^18*O~'\u2030')) +
  xlab("Depth (m)") +
  theme(text = element_text(size=23)) + 
  geom_smooth(method = "lm")


## Depth vs. C13 graph -----
oto %>%
  ggplot(aes(y = Depth,
             x = Rim_C)) +
  geom_point(cex = 4, aes(pch = Area)) + 
  xlab(({delta}^13*C~'\u2030')) + 
  ylab("Depth (m)") +
  theme(text = element_text(size=23)) + 
  geom_smooth(method = "lm")

## Summary regressions ------

summary(lm(oto$Rim_C ~ oto$Depth))

summary(lm(oto$Rim_O ~ oto$Depth))

summary(aov(oto$Rim_C ~ oto$Depth))
summary(aov(oto$Rim_O ~ oto$Depth))


## Standard deviations of replicates -----
oto %>%
  filter(Individual %in% c("158")) %>%
  mutate(Individual.1 = parse_number(Individual.1)) %>%
  group_by(Individual.1) %>%
  summarise_at(vars(Core_C, Core_O, Rim_C, Rim_O), sd) 

%>%
  mutate_all(mean, na.rm=T) %>%
  unique()

## Other replicates , "338", "338-B", "338?", "130"
## Testing out differences and size 

oto %>% 
  mutate(dif_C = Rim_C - Core_C) %>%
  mutate(dif_O = Rim_O - Core_O) %>%
  select(Depth, dif_O, dif_C, Wt.g.) %>%
  ggplot(aes(x = dif_C, y = Wt.g.)) + 
  geom_point() + 
  geom_smooth(method = "lm")


cat = oto %>% 
  mutate(dif_C = Rim_C - Core_C) %>%
  mutate(dif_O = Rim_O - Core_O) %>%
  select(Depth, dif_O, dif_C, Wt.g.)


summary(lm(cat$dif_C ~ cat$Wt.g.))

## Sites Map ------ 

lat_mat = data.frame(lat = oto$lat, 
           long = oto$long, depth = oto$Depth) %>%
  unique()

# Creating map 
data("countryExData", envir=environment(), package="rworldmap")

mymap <- joinCountryData2Map(countryExData, 
                             joinCode = "ISO3",
                             nameJoinColumn = "ISO3V10", 
                             mapResolution = "low")
mymap <- fortify(mymap) 

mypoints <- data.frame(lat = rep(55, 3),
                       long = c(-145, -147, -149))

# Plotting Map and Data 
ggplot() + 
  coord_map(xlim = c(-85, -79.5), 
            ylim = c(24.7, 30)) +
  geom_polygon(data = mymap,
               aes(long, lat, group = group), 
               color = "grey20",
               fill = "grey15", 
               size = 0.3) +
  xlim(-90, -79) +
  ylim(23, 39) +
  geom_point(data=lat_mat, 
             aes(x=long,
                 y=lat,
                 col =2),
             size=5) + 
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = "none") +
  theme(text = element_text(size=26)) 





## Zoom in to WFS sites

ggplot() + 
  coord_map(xlim = c(-84, -80), 
            ylim = c(25, 26.8)) +
  geom_polygon(data = mymap,
               aes(long, lat, group = group), 
               color = "grey20",
               fill = "grey15", 
               size = 0.3) +
  theme_classic() + 
  xlim(-90, -79) +
  ylim(23.5, 39) +
  geom_point(data=lat_mat, 
             aes(x=long,
                 y=lat,
                 col = -depth),
             size=5) + 
  xlab("Longitude") +
  ylab("Latitude") + 
  labs(color= "Depth")


+
  theme(legend.position = "none") 

state <- map_data("state")
Florida <- subset(state, region=="florida")

ggplot() + 
  coord_fixed(1.3) + 
  geom_polygon(data = Florida, mapping = aes(x = long, y = lat, group = group), color="black", fill="gray") + 
  
  geom_polygon(color="black", fill=NA) + 
  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+ 
  geom_point(aes(x = lat_mat$long, y =lat_mat$lat, color = 2))



## Zoom in to Keys sites 
ggplot() + 
  coord_map(xlim = c(-82, -80.2), 
            ylim = c(24.5, 25.3)) +
  geom_polygon(data = mymap,
               aes(long, lat, group = group), 
               color = "grey20",
               fill = "grey15", 
               size = 0.3) +
  xlim(-90, -79) +
  ylim(23, 39) +
  theme_classic() + 
  geom_point(data=lat_mat, 
             aes(x=long,
                 y=lat,
                 col =-depth),
             size=5) + 
  xlab("Longitude") +
  ylab("Latitude") +
  labs(color = "Depth")
  
  
  
  theme(legend.position = "none") 



## Zoom out to either south east or whole USA


usa <- map_data('usa')
state <- map_data("state")
ggplot() + 
  theme_classic() +
  geom_polygon(data = state, aes(x=long, y=lat, group=group),color = "white", fill='black') + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  coord_fixed(1.3) + 
  geom_point(aes(x = lat_mat$long, y =lat_mat$lat, color = -lat_mat$depth)) + 
  theme(legend.position = "none") 
  




oto$Individual %>% unique() %>% length()



### Alex's figure comments 


## Sankey ---------------------
## https://www.r-graph-gallery.com/sankey-diagram.html - Check out this link for tutorial 
### I manually move the labels to outside the diagram so don't search that out

data = read.csv("../Data/sankey_dataframe_new.csv", header=T)



nodes = data.frame(
  name = c(as.character(data$source), 
           as.character(data$target)) %>% unique()
)
data$val = as.character(data$val)
#data$value = 1 I am not sure what the point of this line is. Sorry future me! 

data$IDsource = match(data$source, nodes$name)-1
data$IDtarget = match(data$target, nodes$name)-1
## data = data[-which(data$val == "0"),] Don't seem to need this with this new data sheet



nodes = data.frame(name = rep(c("W:<45m",
                                "W:45-90m",
                                "W:>90m",
                                "K:<10m",
                                "K:>10m"),2))




ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'

p <- sankeyNetwork(Links = data,
                   Nodes = nodes,
                   Source = "IDsource",
                   Target = "IDtarget",
                   Value = "val",
                   NodeID = "name", 
                   colourScale=ColourScal,
                   sinksRight=FALSE,
                   nodeWidth=40,
                   fontSize=30,
                   nodePadding = 20)


p


#### playing around with something that I shouldn't because this is supposed to be clean -----------------------

oto %>%
  select(Area,Rim_C, Rim_O, Core_C, Core_O, Wt.g., TL.mm., Age) %>%
  mutate(dif_C = Rim_C - Core_C) %>%
  mutate(dif_O = Rim_O - Core_O) %>%
  ggplot(aes(x = TL.mm., y = dif_C, color = Area)) + geom_point()



