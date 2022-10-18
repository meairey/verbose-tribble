library(gdata)
library(R2WinBUGS)

## Strip.names function -------------------
strip.names = function (x){
  x <- as.matrix(x)
  dimnames(x) <- list(NULL, NULL)
  x}

## otolith init function ---------------------------
otolith.inits = function (otolith.object, type = c("random", "fixed")[1], n.chains = 1, 
                          max.prob = 0.9) 
{
  out.fin <- vector("list", length = n.chains)
  for (g in 1:n.chains) {
    n.mix <- length(ymix.pos <- grep("y.mix", names(otolith.object)))
    ymix.names <- names(otolith.object)[ymix.pos]
    out1 <- vector("list", length = n.mix)
    lambda <- matrix(1, otolith.object$G, otolith.object$V)
    tau <- matrix(0, otolith.object$V, otolith.object$V)
    diag(tau) <- 1
    n.base <- length(unique(otolith.object$Train.oto))
    if (type == "random") {
      for (i in 1:n.mix) {
        assign(paste("T", i, sep = ""), sample(1:n.base, 
                                               dim(otolith.object[ymix.names][[i]])[1], replace = TRUE))
        out1[[i]] <- get(paste("T", i, sep = ""))
        names(out1)[i] <- paste("T", i, sep = "")
      }
    }
    else {
      for (i in 1:n.mix) {
        PRob <- rep(0, n.base)
        PRob[sample(1:n.base, 1)] <- max.prob
        PRob[PRob == 0] <- (1 - max.prob)/(n.base - 1)
        assign(paste("T", i, sep = ""), sample(1:n.base, 
                                               dim(otolith.object[ymix.names][[i]])[1], prob = PRob, 
                                               replace = TRUE))
        out1[[i]] <- get(paste("T", i, sep = ""))
        names(out1)[i] <- paste("T", i, sep = "")
      }
    }
    out.fin[[g]] <- as.list(c(list(lambda = lambda), out1, 
                              list(tau = tau)))
  }
  return(as.list(out.fin))
}

## Otolith data function --------------
otolith.data = function (data, groups, otolith.var, base, basetows = NULL, mix, 
                         mixtows = NULL) 
{
  if (!all(is.element(levels(base[[groups]]), levels(data[[groups]])))) {
    print("One or more base groups is not present in the data.")
  }
  else if (!all(is.element(levels(mix[[groups]]), levels(data[[groups]])))) {
    print("One or more mix groups is not present in the data.")
  }
  else if (!is.null(basetows) && !all(is.element(basetows, 
                                                 levels(data[[groups]])))) {
    print("One or more basetow groups is not present in the data.")
  }
  else if (!is.null(mixtows) && !all(is.element(mixtows, levels(data[[groups]])))) {
    print("One or more mixtow groups is not present in the data.")
  }
  if (!is.null(basetows)) {
    base <- base[is.element(as.character(base[[groups]]), 
                            basetows), ]
    base[[groups]] <- drop.levels(base[[groups]])
  }
  else {
    base[[groups]] <- drop.levels(base[[groups]])
  }
  if (!is.null(mixtows)) {
    mix <- mix[is.element(as.character(mix[[groups]]), mixtows), 
    ]
    mix[[groups]] <- drop.levels(mix[[groups]])
  }
  else {
    mix[[groups]] <- drop.levels(mix[[groups]])
  }
  n.base.groups <- length(unique(base[[groups]]))
  n.mix.groups <- length(unique(mix[[groups]]))
  basecheck <- summary(base[[groups]])
  mixcheck <- summary(mix[[groups]])
  for (i in 1:n.base.groups) {
    if (basecheck[[i]] < 3) {
      return(print(paste("Error, base group ", names(basecheck[i]), 
                         " has <3 observations")))
    }
  }
  for (m in 1:n.mix.groups) {
    if (mixcheck[[m]] < 3) {
      return(print(paste("Error, mix group ", names(mixcheck[m]), 
                         "  has <3 observations")))
    }
  }
  if (!all(is.element(otolith.var, names(data)))) {
    var.diff <- setdiff(otolith.var, intersect(otolith.var, 
                                               names(data)))
    print(paste("Variable ", var.diff, " does not exist"))
  }
  ngroups <- n.base.groups
  if (is.factor(base[[groups]])) {
    T.groups <- as.numeric(as.factor(T.group.names <- as.character(base[[groups]])))
  }
  else {
    T.groups <- as.numeric(as.factor(base[[groups]]))
  }
  N.cols <- (1:dim(base)[2])[is.element(names(base), otolith.var)]
  y.base <- base[, N.cols]
  y.base <- y.base[order(T.groups), ]
  y.base <- strip.names(y.base)
  T.groups <- sort(T.groups)
  mixtows <- sort(mixtows)
  N.mixtows <- length(mixtows)
  out1 <- out2 <- vector("list", length = N.mixtows)
  out3 <- vector("list", length = ngroups)
  for (i in 1:N.mixtows) {
    mixit <- mix[mix[groups] == mixtows[i], N.cols]
    assign(paste("y.mix", i, sep = ""), mixit)
    out1[[i]] <- get(paste("y.mix", i, sep = ""))
    names(out1)[i] <- paste("y.mix", i, sep = "")
    assign(paste("N.mix", i, sep = ""), dim(get(paste("y.mix", 
                                                      i, sep = "")))[1])
    out2[[i]] <- get(paste("N.mix", i, sep = ""))[1]
    names(out2)[i] <- paste("N.mix", i, sep = "")
  }
  out1 <- sapply(out1, strip.names, simplify = FALSE)
  num.vars <- length(N.cols)
  T.R <- diag(num.vars)
  for (i in 1:ngroups) {
    assign(paste("y.mean", i, sep = ""), as.vector(apply(y.base[T.groups == 
                                                                  i, ], 2, mean)))
    out3[[i]] <- get(paste("y.mean", i, sep = ""))
    names(out3)[i] <- paste("y.mean", i, sep = "")
  }
  return(as.list(c(list(R = T.R, NTrain = diag(table(T.groups, 
                                                     T.groups))), out2, list(Ntrain.oto = dim(base)[1], G = ngroups, 
                                                                             V = num.vars, Train.oto = T.groups, alpha = rep(1/ngroups, 
                                                                                                                             ngroups)), out3, out1, list(y.base = y.base))))
}



### Otolith model function --------------

otolith.model = 
  function (otolith.data, mixparameters = NULL) {
    n.mix <- grep("N.mix", names(otolith.data))
    G <- otolith.data$G
    lnmix <- length(n.mix)
    params <- vector("character", length = (lnmix * 2) + 2)
    for (i in 1:lnmix) {
      params[i] <- paste("T", i, sep = "")
      params[i + lnmix] <- paste("P", i, sep = "")
    }
    params[lnmix * 2 + 1] <- paste("lambda")
    params[lnmix * 2 + 2] <- paste("tau")
    x.out <- paste("model;", "\n", "{")
    cat(x.out, "\n", file = "mixmodel.txt", append = FALSE)
    for (i in 1:length(n.mix)) {
      x1 <- paste("for( i in 1 : ", names(otolith.data[n.mix[i]]), 
                  ")  {", sep = "")
      x6 <- paste("T", i, "[i] ~ dcat(P", i, "[])", sep = "")
      x2 <- paste("y.mix", i, "[i,1:V] ~ dmnorm(mu.", i, "[i,1:V], tau[,])", 
                  sep = "")
      x3 <- paste("for(n in 1:V){")
      x4 <- paste("mu.", i, "[i,n] <- lambda[T", i, "[i],n]", 
                  sep = "")
      x5 <- paste("}")
      x7 <- paste("}")
      x.out <- paste(x1, x6, x2, x3, x4, x5, x7, sep = "\n")
      cat(x.out, "\n", file = "mixmodel.txt", append = TRUE)
    }
    x1 <- paste("for( j in 1 : Ntrain.oto) {", sep = "")
    x2 <- paste("y.base[j,1:V] ~ dmnorm(mu[j,1:V], tau[1:V,1:V])")
    x3 <- paste("for(n in 1:V){")
    x4 <- paste("mu[j,n] <- lambda[Train.oto[j],n]")
    x5 <- paste("}")
    x6 <- paste("}")
    x.out <- paste(x1, x2, x3, x4, x5, x6, sep = "\n")
    cat(x.out, "\n", file = "mixmodel.txt", append = TRUE)
    for (n in 1:length(n.mix)) {
      x.out <- paste("P", n, "[1:G] ~ ddirch(alpha[])", sep = "")
      cat(x.out, "\n", file = "mixmodel.txt", append = TRUE)
    }
    x1 <- paste("tau[1:V,1:V] ~ dwish(R[,],V)")
    x2 <- paste("sigma[1:V,1:V]<-inverse(tau[ , ])")
    x3 <- paste("for(k in 1:G){")
    x4 <- paste("for(kk in 1:V) {")
    x5 <- paste("lambda[k,kk]~dnorm(0.0,1.0E-6)")
    x6 <- paste("}")
    x7 <- paste("}")
    x.out <- paste(x1, x2, x3, x4, x5, x6, x7, sep = "\n")
    cat(x.out, "\n", file = "mixmodel.txt", append = TRUE)
    x1 <- paste("for( j in 1 : Ntrain.oto) {")
    x2 <- paste("yrep[j,1:V] ~ dmnorm(mu[j,1:V], tau[1:V,1:V])")
    x3 <- ("}")
    x.out <- paste(x1, x2, x3, sep = "\n")
    cat(x.out, "\n", file = "mixmodel.txt", append = TRUE)
    x.out <- paste("for( i in 1 :V) {")
    cat(x.out, "\n", file = "mixmodel.txt", append = TRUE)
    for (m in 1:G) {
      y1 <- paste("NTrain[", m, "]", sep = "")
      if (m == 1) {
        xx <- paste("1")
        x1 <- paste("yrepmean", m, "[i]<-mean(yrep[", xx, 
                    ":", y1, ",i])", sep = "")
        yy <- y1
      }
      else {
        xx <- paste("1+", yy, sep = "")
        yy <- paste(yy, "+NTrain[", m, "]", sep = "")
        x1 <- paste("yrepmean", m, "[i]<-mean(yrep[(", xx, 
                    "):(", yy, "),i])", sep = "")
      }
      x2 <- paste("otopredcheck", m, "[i]<-step(yrepmean", 
                  m, "[i]-y.mean", m, "[i])", sep = "")
      x.out <- paste(x1, x2, sep = "\n")
      cat(x.out, "\n", file = "mixmodel.txt", append = TRUE)
    }
    x.out <- paste("}")
    cat(x.out, "\n", file = "mixmodel.txt", append = TRUE)
    x.out <- paste("}")
    cat(x.out, "\n", file = "mixmodel.txt", append = TRUE)
    if (is.null(mixparameters)) {
      return(params)
    }
  }


### Mixfish Function ----------------------

mixy = function (data, method = c("combined", "otolith", "genetic")[1], 
          groups, loci = NULL, otolith.var = NULL, mixparameters = NULL, 
          base, basetows = NULL, mix, mixtows = NULL, gmethod = c("pbayes", 
                                                                  "rmprior")[1], type = c("random", "fixed")[1], max.prob = 0.9, 
          n.chains = 2, n.iter = 200, n.burnin = floor(n.iter/2), 
          n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/n.sims)), 
          n.sims = 1000, debug = FALSE, tol = 1e-04) {
  if (method == "combined") {
    data.otolith <- otolith.data(data, groups, otolith.var, 
                                 base, basetows, mix, mixtows)-
    data.genetic <- genetic.data(data, groups, loci, base, 
                                 basetows, mix, mixtows, gmethod = gmethod, tol = tol)
    mixdata <- mergelist(data.otolith, data.genetic)
    mixinits <- otolith.inits(mixdata, type, n.chains, max.prob = max.prob)
    if (is.null(mixparameters)) {
      mixparameters <- combined.model(mixdata, loci, mixparameters)
    }
    else combined.model(mixdata, loci, mixparameters)
  }
  else if (method == "otolith") {
    mixdata <- otolith.data(data, groups, otolith.var, base, 
                            basetows, mix, mixtows)
    mixinits <- otolith.inits(mixdata, type, n.chains, max.prob = max.prob)
    if (is.null(mixparameters)) {
      mixparameters <- otolith.model(mixdata, mixparameters)
    }
    else otolith.model(mixdata, mixparameters)
  }
  else if (method == "genetic") {
    mixdata <- genetic.data(data, groups, loci, base, basetows,
                            mix, mixtows, gmethod = gmethod, tol = tol)
    mixinits <- genetic.inits(mixdata, type, n.chains, max.prob = max.prob)
    if (is.null(mixparameters)) {
      mixparameters <- genetic.model(mixdata, loci, mixparameters)
    }
    else genetic.model(mixdata, loci, mixparameters)
  }
  bayes <- bugs(mixdata, mixinits, mixparameters, model.file = "mixmodel.txt", 
                n.chains = n.chains, n.iter = n.iter, debug = debug, 
                n.burnin = n.burnin, n.thin = n.thin, n.sims = n.sims, 
                bugs.directory = "c:/Program Files/WinBUGS14/", program = c("WinBUGS", 
                                                                            "OpenBUGS", "winbugs", "openbugs"), working.directory = getwd())
  return(bayes)
}

bugs()

### Import data from project directory ----------
library(dplyr)
data=read.csv("Data/lionmix.csv", header=T)
data = read.csv("Data/lionmix_modified.csv")
data = data %>% 
  filter(Final == "Y") %>%
  select(-c(STOCK_1, Final)) %>%
  na.omit()

## Define model properties -------------------
method="otolith"
groups="STOCK"
base=data[data$SEASON=="R",]
mix=data[data$SEASON=="C",]
type="random"
mixparameters = c("P1", "P2","P3","P4","P5","P6")
mixtows=c("A","B","C","D","E","G")
basetows = c("A","B","C","D","E","G")
otolith.var = c("CARBON","OXYGEN")
n.chains = 2
n.iter = 2000
cat = otolith.data(data, groups, otolith.var, base, 
                   basetows, mix, mixtows)

## Modified Model -
# Just deep and shallow WFS

lionmix_m = mixy(data = data, 
                 method = "otolith",
                 groups = "STOCK",
                 base = data[data$SEASON == "R",],
                 mix = data[data$SEASON == "C",],
                 type = "random",
                 mixparameters = c("P1", "P2"),
                 mixtows = c("A", "B"),
                 basetows = c("A","B"),
                 otolith.var = c("CARBON","OXYGEN"),
                 n.chains = 3,
                 n.iter = 10000 )
lionmix_m
plot(lionmix_m)

## Modified Model -
# Multiple depths WFS

data = read.csv("Data/lionmix_modified.csv")
data = data %>% 
  filter(Final == "Y") %>%
  select(-c(STOCK_1, STOCK, Final)) %>%
  na.omit()



lionmix_m = mixy(data = data, 
                 method = "otolith",
                 groups = "STOCK_2",
                 base = data[data$SEASON == "R",],
                 mix = data[data$SEASON == "C",],
                 type = "random",
                 mixparameters = c("P1", "P2", "P3"),
                 mixtows = c("A", "B","C"),
                 basetows = c("A","B","C"),
                 otolith.var = c("CARBON","OXYGEN"),
                 n.chains = 5,
                 n.iter = 10000 )

lionmix_m
plot(lionmix_m)




## Modified Model -----------------------------
# Multiple WFS and Keys - USE ME

## You need to have this in the main folder i was having an issue with the directory. FIle needs to be in directory folder bc thats where it writes innits

data = read.csv("Data/lionmix_modified.csv")
data = data %>% 
  #filter(STOCK_2 != "KD") %>%
  select(-c(STOCK_1,STOCK, Final)) %>%
  mutate(STOCK_2 = as.factor(STOCK_2)) %>%
  na.omit()

lionmix_mo = mixy(data = data, 
                 method = "otolith",
                 groups = "STOCK_2",
                 base = data[data$SEASON == "R",],
                 mix = data[data$SEASON == "C",],
                 type = "random",
                 mixparameters = c("P1", "P2","P3","P4"),
                 mixtows = c("A","B","C","K"),
                 basetows = c("A","B","C","K"),
                 otolith.var = c("CARBON","OXYGEN"),
                 n.chains = 5,
                 n.iter = 5000)
lionmix_mo

testing = lionmix_m[["sims.matrix"]]
testing_1 = testing[,12]
quantile(testing_1, probs = c(.05,.95, .976))


mean = testing %>%
  as.data.frame() %>%
  summarise_all(mean)
q_5 = testing %>%
  as.data.frame() %>%
  summarise_all(.,funs(quantile(., probs = .05)))
q_95 = testing %>%
  as.data.frame() %>%
  summarise_all(.,funs(quantile(., probs = .95)))
summary = rbind(mean, q_5, q_95)
results = format(round(summary,7))
write.csv(results, "mixfish_results_2.csv")

## Modified Model 
# Just Keys
data = read.csv("Data/lionmix_modified.csv")
data = data %>% 
  #filter(STOCK_1 == "K") %>%
  select(-c(STOCK_1, Final, STOCK_2)) %>%
  mutate(STOCK = as.factor(STOCK)) %>%
  na.omit()

lionmix_mu = mixy(data = data, 
                 method = "otolith",
                 groups = "STOCK",
                 base = data[data$SEASON == "R",],
                 mix = data[data$SEASON == "C",],
                 type = "random",
                 mixparameters = c("P1", "P2"),
                 mixtows = c( "K", "KD"),
                 basetows = c("K", "KD"),
                 otolith.var = c("CARBON","OXYGEN"),
                 n.chains = 3,
                 n.iter = 10000)

lionmix_mu



## Trying with the filtered data 
`%nin%` = Negate(`%in%`) # sets up a way to exclude if in a string
data = pca_dat %>% filter(Sample_ID %nin% list.95.samp) %>%
  select(Group,age, O, C) %>%
  rename("CARBON" = C) %>%
  rename("OXYGEN" = O) %>%
  rename("SEASON" = age) %>%
  rename("STOCK" = Group) %>%
  mutate(STOCK = as.factor(STOCK)) %>% 
  as.data.frame()

#data = read.csv("Data/lionmix_modified.csv")
#data = data %>% 
 # filter(STOCK_2 != "KD") %>%
  #select(-c(STOCK_1,STOCK, Final)) %>%
  #rename("STOCK" = STOCK_2)
  #na.omit()

lionmix_my = mixy(data = data, 
                 method = "otolith",
                 groups = "STOCK",
                 base = data[data$SEASON == "Rim",],
                 mix = data[data$SEASON == "Core",],
                 type = "random",
                 mixparameters = c("P1", "P2", "P3"),
                 mixtows = c("W_A","W_B","W_C"),
                 basetows = c("W_A","W_B","W_C"),
                 otolith.var = c("CARBON","OXYGEN"),
                 n.chains = 5,
                 n.iter = 5000)
lionmix_my
