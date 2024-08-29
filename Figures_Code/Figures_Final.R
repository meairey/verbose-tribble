## Read in libraries ---------
library(ggplot2)
library(dplyr)
library(gridExtra)
library(mapproj)

library(tidyr)
library(tidyverse)
library(networkD3)
library(usmap)
library(SIBER)

library("RColorBrewer")
## Jags is installing in a strange place and I dont have time to figure it out. So use the line below to tell the computer where to find Rjags and it should load the package in fine 
Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.0")
library(rjags)




#install.packages("rjags")

## Read in data -------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/verbose-tribble")
oto=read.csv("Data/Otolith_Data_CSV.csv") %>%
  filter(Area != "NA")

## PCA of rim and core values 

pca_dat = oto %>% 
  filter(Individual != "338-B",  
         Individual != "338?") %>%
  select(Core_O,Rim_O,Core_C,
         Rim_C, Group, Individual) %>% 
  rename("Sample_ID" = Individual) %>%
  group_by(Group, Sample_ID) %>%
  #summarise_at(c("Core_O", "Rim_O", "Core_C", "Rim_C"), mean, na.rm = TRUE) %>%
  summarise(across(c(Core_O, Rim_O, Core_C, Rim_C), list(mean = ~mean(.)))) %>%
  ungroup() %>% 
  mutate(individal = row_number()) %>% 
  pivot_longer(3:6, names_to = "group") %>%
  #pivot_longer(1:4, names_to = "group") %>%
  #filter(individal %nin% c(9,10)) %>%
  separate(group, into = c("age", "isotope", "mean")) %>%  
  pivot_wider(names_from = isotope, values_from = value) %>%
  unite("group",c(Group, age), remove = F) %>%
  na.omit() %>%
  select(-mean)



results = prcomp(pca_dat %>% select(O,C))


results$x %>% as.data.frame() %>% 
  ggplot(aes(x = PC1, y = PC2, color = pca_dat$group)) +
  geom_text(aes(label = as.character(pca_dat$individal))) +
  stat_ellipse(level = .99)+
  scale_color_manual(
    values = as.numeric(as.factor(pca_dat$age))
  )


## Fish sizes 

oto %>% 
  filter(Individual != "338-B",  
         Individual != "338?") %>% select(Group, Individual,TL.mm.) %>% 
  unique() %>% group_by(Group) %>% summarize(mean = median(TL.mm., na.rm=T))



## Group groupings - I think this is the one you should use 

list.95 = c(22,37,16,3,14,10,23,17,8,34,35,28,4)
list.95.samp = (pca_dat %>% filter(individal %in% list.95))$Sample_ID %>% unique
list.975 = c(22,37,16,3,14,10,23,8,34,28,4)
list.975.samp = (pca_dat %>% filter(individal %in% list.975))$Sample_ID %>% unique
list.99 = c(37,16,3,10,14,8)
list.99 = c(16,3,14,8,14,37)
list.99.samp = (pca_dat %>% filter(individal %in% list.99))$Sample_ID %>% unique
# Group.C groupings
## These are the individuals that do not fit into an ellipse with 97.5% stat ellipses
list.975 = c(4,6,8,19,23,1)
list.975.samp = (pca_dat %>% filter(individal %in% list.975))$Sample_ID %>% unique
list.95 = c(4,6,8,19,23,1,35,10,14)
list.95.samp = (pca_dat %>% filter(individal %in% list.95))$Sample_ID %>% unique


# I combined the K_A and K_B groups because their isotopes looked very similar? Here are new groupings that are outside the stat ellipses

list.975 = c(1,6,16,3,12,7,5,15)
list.975.samp = (pca_dat %>% filter(individal %in% list.975))$Sample_ID %>% unique()
list.95 = c(1,6,8,16,3,12,7,5,15)
list.95.samp = (pca_dat %>% filter(individal %in% list.975))$Sample_ID %>% unique
list.99 = c(1,3,12,7,15)

## I combined K_A and W_A because their isotopes arent actually that different from eachother

list.975 = c(13,10,1,3,12,6,14,9)
list.975.samp = (pca_dat %>% filter(individal %in% list.975))$Sample_ID %>% unique()

## I combined W_B and W_C together because their isotopes are similar
list.99 = c(1,13,10,3,6,12)
list.99.samp = list.975.samp = (pca_dat %>% filter(individal %in% list.99))$Sample_ID %>% unique()
list.975 = c(1,13,10,3,12,6,14,18)
list.975.samp = (pca_dat %>% filter(individal %in% list.975))$Sample_ID %>% unique()
list.95 = c(1,13,10,3,12,6,14,18,19)
list.95.samp = (pca_dat %>% filter(individal %in% list.95))$Sample_ID %>% unique()


## Graph of d18O for un-classifiable individuals 

pca_dat %>% filter(individal %in% list.975.samp) %>% ggplot(aes(x = C, y = O, color = Group, shape = age)) + 
  geom_point()

pca_dat %>% filter(individal %in% list.95) %>% ggplot(aes(x =individal, y = O, color = Group, shape = age)) + 
  geom_point()

## Rim vs. Core d18O graph ----
ggplot(data=oto, aes(x = Rim_O, y=Core_O, color=Depth, shape = Area)) + 
  geom_point(size=5)  + 
  theme_minimal() + 
  scale_color_gradient(low='gray', high='black') + 
  geom_abline(intercept=0) + 
  ylim(-1.3,1) +
  xlim(-1.3,1) + 
  ylab((expression({delta}^18*O~'\u2030')) )+
  xlab((expression({delta}^18*O~'\u2030')) )+
  theme(axis.text=element_text(size=13))

## Rim vs. Core d13C graph ---- 
ggplot(data=oto, aes(x = Rim_C, y=Core_C, color=Depth, shape = Area)) +  
  geom_point(size=5)  + 
  scale_color_gradient(low='gray', high='black') + 
  geom_abline(intercept=0) + 
  ylim(-7.7, -3) +
  xlim(-7.7,-3) +
  theme_minimal() + 
  ylab((expression({delta}^13*C~'\u2030')) )+
  xlab((expression({delta}^13*C~'\u2030')) )+
  theme(axis.text=element_text(size=13)) 

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
  ylab(expression({delta}^18*O~'\u2030')) + 
  ggtitle("delete")

## Depth vs. 018 graph -----

summary(lm(oto$Rim_C ~ oto$Depth,na.rm = T))

oto %>%
  ggplot(aes(x = Depth, 
             y = Rim_O)) +
  theme_minimal() +
  geom_point(cex = 4) + 
  ylab(({delta}^18*O~'\u2030')) +
  xlab("Depth (m)") +
  theme(text = element_text(size=23)) + 
  geom_smooth(method = "lm", color = 1, se = F) + 
  #geom_vline(xintercept = 0, linetype = 1) +
  geom_vline(xintercept = 35, linetype = 1) + 
  ggtitle("Delete")
  #geom_vline(xintercept = 95, linetype = 3) 


## Depth vs. C13 graph -----
oto %>%
  ggplot(aes(x = Depth,
             y = Rim_C)) +
  theme_minimal() +
  geom_point(cex = 4) + 
  ylab(({delta}^13*C~'\u2030')) + 
  xlab("Depth (m)") +
  theme(text = element_text(size=23)) + 
  geom_smooth(method = "lm", color = 1, se = F) + 
  #geom_vline(xintercept = 0, linetype = 1) +
  geom_vline(xintercept = 35, linetype = 1) +
  ggtitle("delete")
  #geom_vline(xintercept = 95, linetype = 3)

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
                     long = oto$long, depth = oto$Depth, group = oto$Group) %>%
  unique()

#install.packages("rworldmap")
library(rworldmap)
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
            ylim = c(24, 26.8)) +
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
                 col = -depth, 
                 shape = group),
             size=5) + 
  xlab("") +
  ylab("") + 
  labs(color= "Depth", 
       shape = "Depth Step") +
  scale_shape_discrete(labels = c("Shallow: 0-35m", "Deep: 45 - 93m")) + theme(axis.text = element_text(size = 18))



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
            ylim = c(24, 25.3)) +
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
                 col =-depth, 
                 shape = group),
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
  geom_point(aes(x = lat_mat$long, y =lat_mat$lat, color = -lat_mat$depth, 
                 shape = lat_mat$group)) + 
  theme(legend.position = "none") 


lat_mat



oto$Individual %>% unique() %>% length()



### Alex's figure comments 


## Sankey ---------------------
## https://www.r-graph-gallery.com/sankey-diagram.html - Check out this link for tutorial 
### I manually move the labels to outside the diagram so don't search that out

data = read.csv("../Data/sankey_dataframe_new.csv", header=T)

data = read.csv("Data/sankey_dataframe_new.csv", header=T)

data = read.csv("Data/sankey_dataframe_102322.csv", header = T)


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

nodes = data.frame(name = rep(c("< 45m",
                                "> 45m"),2))




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

pdf("sankey_newboogy.pdf", height = 11, width = 11, onefile = T)
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
dev.off()


#### playing around with something that I shouldn't because this is supposed to be clean -----------------------

oto %>%
  select(Area,Rim_C, Rim_O, Core_C, Core_O, Wt.g., TL.mm., Age) %>%
  mutate(dif_C = Rim_C - Core_C) %>%
  mutate(dif_O = Rim_O - Core_O) %>%
  ggplot(aes(x = TL.mm., y = dif_C, color = Area)) + geom_point()


oto %>% group_by(Group_old) %>% select(Group_old,Group_new, Individual) %>% unique() -> cat

oto %>% ggplot(aes(x = Rim_C, y = Rim_O, col = Group.C)) + geom_point() + stat_ellipse()

sd_range = oto %>% group_by(Group.C) %>% summarize(mean_c = mean(Rim_C, na.rm = T), sd_c = sd(Rim_C, na.rm = T), mean_O = mean(Rim_O, na.rm = T), sd_O = sd(Rim_O, na.rm = T)) %>% rownames_to_column(var = "group")


## Graph with biplot and 
ggplot() + geom_point(data = sd_range,
                      aes(x = mean_c,y = mean_O, col = group))+
  geom_pointrange(data = sd_range, aes(x = mean_c, y = mean_O,ymin = mean_O - sd_O, ymax = mean_O +sd_O, col = group)) +
  geom_pointrange(data = sd_range, aes(x = mean_c, y = mean_O,xmax = mean_c + sd_c, xmin = mean_c - sd_c, col = group)) + 
  theme_minimal() + 
  geom_point(data = oto, aes(x = oto$Rim_C,
                             y = oto$Rim_O, 
                             col = as.factor(as.numeric(as.factor(oto$Group.C)))
                             ))


oto %>% ggplot(aes(x = Depth, y = water.temp..F...Surface.Bottom.)) + geom_point()


## May 21 2023 - I'm not quite sure what i was doing above. A reviewer asked me to make sure that there were significant differences between chosen groups. I'm using the column Group.C for this analysis which should have three groups KA, WA, WB

WA = oto %>% filter(Group.C == "W_A")
KA = oto %>% filter(Group.C == "K_A")
WB = oto %>% filter(Group.C == "W_B")
KA_WA = oto %>% filter(Group.C %in% c("W_A", "K_A"))

## Shallow are different for C but not different for O
wilcox.test(WA$Rim_O, KA$Rim_O) ## Not significant - no difference 
wilcox.test(WA$Rim_C, KA$Rim_C) ## significant - difference 

## Deep WFS is different from shallow Keys and shallow WFS
wilcox.test(WB$Rim_C, KA$Rim_C) ## sig
wilcox.test(WB$Rim_C, WA$Rim_C) ## sig

wilcox.test(WB$Rim_O, KA$Rim_O) ## sig
wilcox.test(WB$Rim_O, WA$Rim_O) ## sig

## Deep vs. shallow pooled 

wilcox.test(WB$Rim_C, KA_WA$Rim_C) # sig
wilcox.test(WB$Rim_O, KA_WA$Rim_O) # sig


(oto$TL.mm.) %>% min(na.rm =  T)

## Fig 4 Setting this up to go into SIBER -----------------------
## Jags 
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3
## All 3 ellipses in one graph -----------------
siber_oto = oto %>% select(Rim_C, Rim_O, Group.C) %>% as.data.frame() %>% 
  mutate(community = 1) %>% 
  rename(iso1 = Rim_C) %>% 
  rename(iso2 = Rim_O) %>% 
  rename(group = Group.C) %>%
  na.omit()
# SIBER
palette(c("gray","gray","#444444"))
siber.example <- createSiberObject(siber_oto)
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^18*O~'\u2030'),
                x.limits = c(-7.2,-2.5),
                y.limits = c(-2,2))
legend("topright",c("Deep Stock", "Shallow Stock - WFS", "Shallow Stock - Keys"), col = c("black", "gray", "gray"), lty = c(1,1,2), cex = .25)

plot(x = siber_oto$iso1, 
     y =siber_oto$iso2,
     xlab=expression({delta}^13*C~'\u2030'),
     ylab=expression({delta}^18*O~'\u2030'),
     xlim = c(-7.2,-2.5),
     ylim= c(-2,2),
     col = as.numeric(as.factor(siber_oto$group)),
     cex = .5, 
     pch = 18)


legend("topright",c("Deep Stock", "Shallow Stock - WFS", "Shallow Stock - Keys"), col = c("black", "gray", "gray"), lty = c(1,1,2), cex = .25)


cords_ellip3 = addEllipse(siber.example$ML.mu[[1]][ , , 3],                                                        siber.example$ML.cov[[1]][ , , 3],                                              m = NULL,
                          n = 500,
                          p.interval = 0.90,
                          ci.mean = F,
                          col = "black",
                          lty = 1,
                          lwd = 2)
cords_ellip3_mean = addEllipse(siber.example$ML.mu[[1]][ , , 3],                                                        siber.example$ML.cov[[1]][ , , 3],                                     m = 500,
                               n = 500,
                               p.interval = 0.90,
                               ci.mean = T,
                               col = "black",
                               lty = 1,
                               lwd = 5)

cords_ellip1 = addEllipse(siber.example$ML.mu[[1]][ , , 1],                                                        siber.example$ML.cov[[1]][ , , 1],                                              m = NULL,
                n = 500,
                p.interval = 0.90,
                ci.mean = FALSE,
                col = "gray",
                lty = 2,
                lwd = 2)
cords_ellip1_mean = addEllipse(siber.example$ML.mu[[1]][ , , 1],                                                        siber.example$ML.cov[[1]][ , , 1],                                            m = 500,
                          n = 500,
                          p.interval = 0.90,
                          ci.mean = T,
                          col = "gray",
                          lty = 2,
                          lwd = 5)
cords_ellip2 = addEllipse(siber.example$ML.mu[[1]][ , , 2],                                                        siber.example$ML.cov[[1]][ , , 2],                                          m = NULL,
                   n = 100,
                   p.interval = 0.90,
                   ci.mean = FALSE,
                   col = "gray",
                   lty = 1,
                   lwd = 2)
cords_ellip2_mean = addEllipse(siber.example$ML.mu[[1]][ , , 2],                                                        siber.example$ML.cov[[1]][ , , 2],                                          m = 500,
                          n = 500,
                          p.interval = 0.90,
                          ci.mean = T,
                          col = "gray",
                          lty = 5,
                          lwd = 3)

## Final Graphic of all 3 groups 
ggplot() + geom_point(data = siber_oto,aes(x= iso1, y = iso2, col = group, shape = group)) + 
  geom_point(aes(x = cords_ellip1[,1], y = cords_ellip1[,2]),
             col = "#bcbcbc", size = 1)+
  geom_point(aes(x = cords_ellip1_mean[,1], y = cords_ellip1_mean[,2]), 
             col = "#bcbcbc", size = 2) +
  geom_point(aes(x = cords_ellip2[,1], y = cords_ellip2[,2]), 
             col = "#999191", size = 1) + 
  geom_point(aes(x = cords_ellip2_mean[,1], y = cords_ellip2_mean[,2]),
             col =  "#999191", size = 2) +
  geom_point(aes(x = cords_ellip3[,1], y = cords_ellip3[,2]), 
             col = "black", size = 1) + 
  geom_point(aes(x = cords_ellip3_mean[,1], y = cords_ellip3_mean[,2]),
             col = "black", size = 2) +  
  scale_color_manual(values = c("#bcbcbc", "#999191", "#5b5b5b"), labels = c("Shallow - Keys", "Shallow - WFS", "Deep - WFS")) + 
  scale_shape_manual(values = c(3,2,1),labels = c("Shallow - Keys", "Shallow - WFS", "Deep - WFS")) + 
  theme_minimal() + 
  ylab(expression({delta}^18*O~'\u2030')) + 
  xlab(expression({delta}^13*C~'\u2030')) + 
  labs(color = "Stock", shape = "Stock") + 
  theme(text = element_text(size = 17))
## Just 2 ellipses ---------------------------------------
siber_oto = oto %>% select(Rim_C, Rim_O, Group) %>% as.data.frame() %>% 
  mutate(community = 1) %>% 
  rename(iso1 = Rim_C) %>% 
  rename(iso2 = Rim_O) %>% 
  rename(group = Group) %>%
  na.omit()



palette(c("gray", "black"))
siber.example <- createSiberObject(siber_oto)
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^18*O~'\u2030'),
                cex = 0.5,
                x.limits = c(-7.2,-2.5),
                y.limits = c(-2,2))
legend("topright",c("Deep Stock", "Shallow Stock"), col = c("black", "gray"), lty = c(1,1), cex = .5)
cords_ellip3 = addEllipse(siber.example$ML.mu[[1]][ , , 2],                                                        siber.example$ML.cov[[1]][ , , 2],                                              m = NULL,
                          n = 500,
                          p.interval = 0.95,
                          ci.mean = F,
                          col = "black",
                          lty = 1,
                          lwd = 2)
cords_ellip3_mean = addEllipse(siber.example$ML.mu[[1]][ , , 2],                                                        siber.example$ML.cov[[1]][ , , 2],                                     m = 500,
                               n = 500,
                               p.interval = 0.95,
                               ci.mean = T,
                               col = "black",
                               lty = 1,
                               lwd = 5)

cords_ellip1 = addEllipse(siber.example$ML.mu[[1]][ , , 1],                                                        siber.example$ML.cov[[1]][ , , 1],                                              m = NULL,
                          n = 500,
                          p.interval = 0.95,
                          ci.mean = FALSE,
                          col = "gray",
                          lty = 2,
                          lwd = 2)
cords_ellip1_mean = addEllipse(siber.example$ML.mu[[1]][ , , 1],                                                        siber.example$ML.cov[[1]][ , , 1],                                            m = 500,
                               n = 500,
                               p.interval = 0.95,
                               ci.mean = T,
                               col = "gray",
                               lty = 2,
                               lwd = 5)

## Final Graphic of 2 final groups 
ggplot() + geom_point(data = siber_oto,aes(x= iso1, y = iso2, col = group, shape = group)) + 
  geom_point(aes(x = cords_ellip1[,1], y = cords_ellip1[,2]),
             col = "#bcbcbc", size = 1)+
  geom_point(aes(x = cords_ellip1_mean[,1], y = cords_ellip1_mean[,2]), 
             col = "#bcbcbc", size = 2) +
  geom_point(aes(x = cords_ellip3[,1], y = cords_ellip3[,2]), 
             col = "black", size = 1) + 
  geom_point(aes(x = cords_ellip3_mean[,1], y = cords_ellip3_mean[,2]),
             col = "black", size = 2) +  
  scale_color_manual(values = c("#bcbcbc", "#999191", "#5b5b5b"), labels = c("Shallow", "Deep")) + 
  scale_shape_manual(values = c(3,1),labels = c("Shallow", "Deep")) + 
  theme_minimal() + 
  ylab(expression({delta}^18*O~'\u2030')) + 
  xlab(expression({delta}^13*C~'\u2030')) + 
  labs(color = "Stock", shape = "Stock")+ 
  theme(text = element_text(size = 17))

## Overlap of ellipses on biplot ------------------

siber_oto = oto %>% select(Rim_C, Rim_O, Group.C) %>% as.data.frame() %>% 
  mutate(community = 1) %>% 
  rename(iso1 = Rim_C) %>% 
  rename(iso2 = Rim_O) %>% 
  rename(group = Group.C) %>%
  na.omit()
siber.example <- createSiberObject(siber_oto)


## WA + WB -------
ellipse1 <- "1.W_A" 
ellipse2 <- "1.W_B"
ellipses.posterior <- siberMVN(siber.example, parms, priors)
bayes95.overlap = bayesianOverlap(ellipse1, ellipse2, ellipses.posterior,
                                   draws = 1000, p.interval = 0.90, n = 500)
# Main ellipse
main = 2
## Overlapping area as percent of main ellipses's area
non_overlapping_area =  bayes95.overlap[,main] - bayes95.overlap[,3]  
proportion_non_overlap = non_overlapping_area /bayes95.overlap[,main] 
overlapping_area = 100*((bayes95.overlap[,main] - non_overlapping_area) / bayes95.overlap[,main])
#bayes_AB = bayes95.overlap
overlap_AB = overlapping_area
vector_WA.B = c(mean(overlapping_area),quantile(overlapping_area, .95),quantile(overlapping_area, .05))

bayes95.overlap = bayes_AB 

## KA + WB ----------------------
ellipse1 <- "1.K_A" 
ellipse2 <- "1.W_B"

bayes95.overlap = bayesianOverlap(ellipse1, ellipse2, ellipses.posterior,
                                  draws = 1000, p.interval = 0.90, n = 500)
# Main ellipse
main = 1
## Overlapping area as percent of main ellipses's area
non_overlapping_area =  bayes95.overlap[,main] - bayes95.overlap[,3]  
proportion_non_overlap = non_overlapping_area /bayes95.overlap[,main] 
overlapping_area = 100*((bayes95.overlap[,main] - non_overlapping_area) / bayes95.overlap[,main])
#bayes_KAB = bayes95.overlap
overlap_KAB = overlapping_area
vector_KA.B = c(mean(overlapping_area),quantile(overlapping_area, .95),quantile(overlapping_area, .05))
#bayes95.overlap = bayes_KAB

## deep vs. shallow 


siber_oto = oto %>% select(Rim_C, Rim_O, Group) %>% as.data.frame() %>% 
  mutate(community = 1) %>% 
  rename(iso1 = Rim_C) %>% 
  rename(iso2 = Rim_O) %>% 
  rename(group = Group) %>%
  na.omit()
siber.example <- createSiberObject(siber_oto)

ellipse1 <- "1.K_A" 
ellipse2 <- "1.W_B"

bayes95.overlap = bayesianOverlap(ellipse1, ellipse2, ellipses.posterior,
                                  draws = 1000, p.interval = 0.90, n = 500)
# Main ellipse
main = 2
## Overlapping area as percent of main ellipses's area
non_overlapping_area =  bayes95.overlap[,main] - bayes95.overlap[,3]  
proportion_non_overlap = non_overlapping_area /bayes95.overlap[,main] 
overlapping_area = 100*((bayes95.overlap[,main] - non_overlapping_area) / bayes95.overlap[,main])
#bayes_deepshallow = bayes95.overlap
overlap_deepshallow = overlapping_area
vector_deepshallow = c(mean(overlapping_area),quantile(overlapping_area, .95),quantile(overlapping_area, .05))
bayes95.overlap = bayes_deepshallow



##bayes90KA_WA = bayes.prop.95.over
##bayesproportion_overlapping = overlapping_area ## KA_WA
#bayes90KA_WB = bayes.prop.95.over
#bayes90KA_WB_overlap = overlapping_area 
#bayes90WA_WB_overlap= overlapping_area
#bayes90WA_WB = bayes.prop.95.over 



## I think it should be facet wrapped based on depth? and have two ellipses that represent the isotopic niches for shallow vs. deep and then have the core values overlaid onto that? 


packageVersion("r2WINBUGS")
#install.packages("R2WinBUGS")
citation("SIBER")
packageVersion("SIBER")
