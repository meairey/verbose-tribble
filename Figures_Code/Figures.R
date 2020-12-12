## Read in libraries ---------
library(ggplot2)
library(dplyr)
library(gridExtra)
library(mapproj)
library(rworldmap)


## Read in data -------
oto=read.csv("CSVs/Otolith_Data_CSV.csv") %>%
  filter(Area != "NA")


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
ggplot(data=oto, aes(x = Rim_C, y=Core_C, color=Depth, shape = Area)) +  geom_point(size=5)  + 
  scale_color_gradient(low='#91d0f2', high='#39658e') + 
  geom_abline(intercept=0) + 
  ylim(-7.7, -3) +
  xlim(-7.7,-3) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text=element_text(size=13)) 

## Rim vs Core d18O graph --- 
# maybe don't use this one? 
ggplot(oto, aes(Group.C, Rim_O, fill=Group.C)) + geom_boxplot() +scale_fill_brewer(palette = "Blues") + scale_x_discrete(labels=c("A"="0-15","B"="15-30","C"="30-45","D"="45-60", "E"="60-75","G"="90-105")) + xlab("Depth (m)")+  guides(fill=F) +theme(text = element_text(size=26))+ylab(expression({delta}^18*O~'\u2030'))

## Depth vs. 018 graph -----
oto %>%
  ggplot(aes(x = Depth, 
             y = Rim_O, 
             col = Area)) +
  geom_point() + 
  ylab(({delta}^18*O~'\u2030'))

## Depth vs. C13 graph 
oto %>%
  ggplot(aes(x = Depth,
             y = Rim_C, 
             col = Area)) +
  geom_point() + 
  ylab(({delta}^13*C~'\u2030'))


## Sites Map ------ 

lat_mat = data.frame(lat = oto$lat, 
           long = oto$long) %>%
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

## Sankey 
cat = rep(letters[1:6],2)
data = read.csv("Data/sankey_dataframe.csv",header=T)


