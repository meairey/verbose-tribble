## Read in libraries ---------
library(ggplot2)
library(dplyr)
library(gridExtra)

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

