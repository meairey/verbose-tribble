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



## Summary regressions ------
# Determine if there are relationships between isotopes and depth 
summary(lm(oto$Rim_C ~ oto$Depth))

summary(lm(oto$Rim_O ~ oto$Depth))

## Depict the Rim values
core = (pca_dat %>% filter(age == "Core"))
Rim = (pca_dat %>% filter(age == "Rim"))
ggplot() +
  geom_point(aes(x = Rim$C, y = Rim$O, col = Rim$Group)) + 
  stat_ellipse(aes(x = Rim$C, y = Rim$O, col = Rim$Group), level = .99) + 
  geom_point(aes(y = core$O, x = core$C))

shallow = pca_dat %>% filter(age == "Rim", Group == "K_A")
deep = pca_dat %>% filter(age == "Rim", Group == "W_B")


wilcox.test(shallow$C, deep$C)
wilcox.test(shallow$O, deep$O)


shallow_detailed = oto %>% filter(Group == "K_A") %>% select(Group, Group_old, Rim_C, Rim_O)
deep_detailed = oto %>% filter(Group == "W_B") %>% select(Group, Group_old, Rim_C, Rim_O)


## These groupings show that the shallow keys and the shallow wfs can be differentiated by d13C but not by d18O 
wilcox.test((shallow_detailed %>% filter(Group_old == "K_A"))$Rim_O,(shallow_detailed %>% filter(Group_old == "W_A"))$Rim_O)
wilcox.test((shallow_detailed %>% filter(Group_old == "K_A"))$Rim_C,(shallow_detailed %>% filter(Group_old == "W_A"))$Rim_C)
mean((shallow_detailed %>% filter(Group_old == "K_A"))$Rim_C)
mean((shallow_detailed %>% filter(Group_old == "W_A"))$Rim_C)
## These groupings show that shallow WFS is sig different from deep wfs with new groupings... and shallow wfs is sig different from deep wfs 
wilcox.test((shallow_detailed %>% filter(Group_old == "W_A"))$Rim_C,(deep_detailed %>% filter(Group_old %in% c("W_B", "W_C")))$Rim_C)

wilcox.test((shallow_detailed %>% filter(Group_old == "K_A"))$Rim_O,(deep_detailed %>% filter(Group_old %in% c("W_B", "W_C")))$Rim_O)

wilcox.test((shallow_detailed %>% filter(Group_old == "W_A"))$Rim_O,(deep_detailed %>% filter(Group_old %in% c("W_B", "W_C")))$Rim_O)


## These groupings show that shallow is sig different from deep with new groupings... 
wilcox.test((shallow_detailed %>% filter(Group == "K_A"))$Rim_C,(deep_detailed %>% filter(Group_old %in% c("W_B", "W_C")))$Rim_C)
wilcox.test((shallow_detailed %>% filter(Group == "K_A"))$Rim_O,(deep_detailed %>% filter(Group_old %in% c("W_B", "W_C")))$Rim_O)



## Clean Data ----------
## PCA of rim and core values 
## Involves getting the averages of replicates 
pca_dat = oto %>% select(Core_O,Rim_O,Core_C, Rim_C, Group, Individual) %>% 
  rename("Sample_ID" = Individual) %>%
  group_by(Group, Sample_ID) %>%
  summarise(across(c(Core_O, Rim_O, Core_C, Rim_C),
                   list(mean = ~mean(.)))) %>%
  ungroup() %>% 
  mutate(individal = row_number()) %>% 
  pivot_longer(3:6, names_to = "group") %>%
  separate(group, into = c("age", "isotope", "mean")) %>%  
  pivot_wider(names_from = isotope, values_from = value) %>%
  unite("group",c(Group, age), remove = F) %>%
  na.omit() %>%
  select(-mean)

# Run PCA
results = prcomp(pca_dat %>% select(O,C))

# Plot the PCA 
# Use this to see which datapoints fall outside of Rim stat ellipses 
# I'm going with the level .99 because I have so few datapoints I don't want to exclude a lot
results$x %>% as.data.frame() %>% 
  ggplot(aes(x = PC1, y = PC2, color = pca_dat$group)) +
  geom_text(aes(label = as.character(pca_dat$individal))) +
  stat_ellipse(level = .99)+
  scale_color_manual(
    values = as.numeric(as.factor(pca_dat$age)))

# Create the list of individuals to exclude 
list.99 = c(1,13,10,3,6,12)
list.99.samp = list.975.samp = (pca_dat %>% filter(individal %in% list.99))$Sample_ID %>% unique()

## Prepare for mixFish.-------
`%nin%` = Negate(`%in%`) # sets up a way to exclude if in a string
data = pca_dat %>%
  filter(Sample_ID %nin% list.99.samp) %>%
  
  rename("CARBON" = C) %>%
  rename("OXYGEN" = O) %>%
  rename("SEASON" = age) %>%
  rename("STOCK" = Group) %>%
  # mutate(STOCK = replace(STOCK, STOCK == "W_B", "B_W  )) %>% 
  mutate(STOCK = as.factor(STOCK)) %>% ## Has to be a factor to work 
  as.data.frame() 

data = data[-which((pca_dat$age == "Core") & (pca_dat$individal %in% list.99)),]


