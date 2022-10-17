mix.fishs = function (data, 
                      method = c("combined", "otolith", "genetic")[1], 
                      groups,
                      loci = NULL, 
                      otolith.var = NULL,
                      mixparameters = NULL, 
                      base,
                      basetows = NULL,
                      mix,
                      mixtows = NULL,
                      gmethod = c("pbayes","rmprior")[1],
                      type = c("random", "fixed")[1],
                      max.prob = 0.9,
                      n.chains = 2, 
                      n.iter = 200, 
                      n.burnin = floor(n.iter/2),
                      n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/n.sims)),
                      n.sims = 1000,
                      debug = FALSE, 
                      tol = 1e-04) {
   if (method == "otolith") {
    mixdata <- otolith.data(data, groups, otolith.var, base, 
                            basetows, mix, mixtows)
    mixinits <- otolith.inits(mixdata, type, n.chains, max.prob = max.prob)
    if (is.null(mixparameters)) {
      mixparameters <- otolith.model(mixdata, mixparameters)
    }
    else otolith.model(mixdata, mixparameters)
  }
  bayes <- bugs(mixdata, mixinits, mixparameters, model.file = "mixmodel.txt", 
                n.chains = n.chains, n.iter = n.iter, debug = debug, 
                n.burnin = n.burnin, n.thin = n.thin, n.sims = n.sims, 
                bugs.directory = "c:/Program Files/WinBUGS14/", program = c("WinBUGS", 
                                                                            "OpenBUGS", "winbugs", "openbugs"), working.directory = getwd())
  return(bayes)
}



mix=read.csv("../Data/lionmix.csv", header=T)

data = mix
groups = "STOCK"
mixtows = c("A","B","C","D","E","G")
basetows =  c("A","B","C","D","E","G")
otolith.var = c("CARBON","OXYGEN")
max.prob = 0.9
n.chains = 2
n.iter = 2000
n.burnin = floor(n.iter/2)
n.sims = 1000
n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/n.sims))
debug = FALSE
tol = 1e-04
mixparameters = c("P1", "P2","P3","P4","P5","P6")
base = lmix[lmix$SEASON=="R",]
gmethod = c("pbayes","rmprior")[1]

mix = lmix[lmix$SEASON=="C",]
type = c("random")[1]

mixdata = mix

inits1 = list(lambda= structure(.Data= c(1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00), .Dim=c(6, 2)), T1=c(2.00000E+00, 2.00000E+00, 4.00000E+00, 5.00000E+00), T2=c(6.00000E+00, 6.00000E+00, 5.00000E+00), T3=c(4.00000E+00, 5.00000E+00, 5.00000E+00, 6.00000E+00), T4=c(6.00000E+00, 2.00000E+00, 6.00000E+00, 2.00000E+00, 5.00000E+00, 6.00000E+00, 3.00000E+00, 1.00000E+00), T5=c(3.00000E+00, 3.00000E+00, 5.00000E+00, 4.00000E+00), T6=c(4.00000E+00, 5.00000E+00, 5.00000E+00, 5.00000E+00, 1.00000E+00, 1.00000E+00, 4.00000E+00), tau= structure(.Data= c(1.00000E+00, 0.00000E+00, 0.00000E+00, 1.00000E+00), .Dim=c(2, 2)))

inits2 = list(lambda= structure(.Data= c(1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00), .Dim=c(6, 2)), T1=c(6.00000E+00, 2.00000E+00, 4.00000E+00, 1.00000E+00), T2=c(1.00000E+00, 2.00000E+00, 5.00000E+00), T3=c(4.00000E+00, 5.00000E+00, 1.00000E+00, 5.00000E+00), T4=c(3.00000E+00, 1.00000E+00, 4.00000E+00, 3.00000E+00, 2.00000E+00, 3.00000E+00, 3.00000E+00, 6.00000E+00), T5=c(2.00000E+00, 3.00000E+00, 3.00000E+00, 3.00000E+00), T6=c(2.00000E+00, 3.00000E+00, 6.00000E+00, 3.00000E+00, 5.00000E+00, 1.00000E+00, 2.00000E+00), tau= structure(.Data= c(1.00000E+00, 0.00000E+00, 0.00000E+00, 1.00000E+00), .Dim=c(2, 2)))

mixinits = array(c(inits1, inits2),2)




bayes <- bugs(mixdata, mixinits, mixparameters, model.file = "mixmodel.txt", 
              n.chains = n.chains, n.iter = n.iter, debug = TRUE, 
              n.burnin = n.burnin, n.thin = n.thin, n.sims = n.sims, 
              bugs.directory = "C:/Program Files/WinBUGS14/", program = c("WinBUGS", 
                                                                          "OpenBUGS", "winbugs", "openbugs"),)
bayes = bugs(data = mixdata, 
             inits = mixinits, 
             parameters.to.save = mixparameters,
             model.file = "mixmodel.txt", 
             n.chains = n.chains, 
             n.iter = n.iter,
             debug = debug,
             n.burnin = n.burnin, 
             n.thin = n.thin,
             n.sims = n.sims,
             bugs.directory = "C:/Program Files/WinBUGS14/", 
             program = c("WinBUGS","OpenBUGS", "winbugs", "openbugs"),
             working.directory = getwd()
             )




                    
                      
                       
  if (method == "otolith") {
    mixdata <- otolith.data(data, groups, otolith.var, base, 
                            basetows, mix, mixtows)
    mixinits <- otolith.inits(mixdata, type, n.chains, max.prob = max.prob)
    if (is.null(mixparameters)) {
      mixparameters <- otolith.model(mixdata, mixparameters)
    }
    else otolith.model(mixdata, mixparameters)
  }




#install.packages("car")
library(coda) 
library(car)
data(Angell)
attach(Angell)
res <- lm(moral ~ hetero + mobility)
summary(res)

N <- nrow(Angell)
data <- list(moral=moral, hetero=hetero, mobility=mobility, N=N)
inits <- function(){
  list(alpha=rnorm(1), beta1=rnorm(1), beta2=rnorm(1), tau=runif(1,0,2))
}
parameters <- c("alpha", "beta1", "beta2", "tau")
sims <- bugs(model.file="model_angell.txt",
             data = data,
             parameters = parameters,
             inits = inits,
             n.chains = 2,
             n.iter = 25000, n.burnin = 20000, n.thin = 5,
             bugs.directory =  "C:/Program Files/WinBUGS14/")







mixy = function(data, groups, otolith.var = NULL, mixparameters = NULL,
                base, basetows = NULL, mix, mixtows = NULL, type = c("random", "fixed")[1], 
                max.prob = .9, 
                n.chains = 2, n.iter = 200, n.burnin = 1000,
                n.thin =  max(1, floor(n.chains * (n.iter - n.burnin)/n.sims)), 
                n.sims = 1000, debug = FALSE, tol = 1e-04){
  mixdata = otolith.data(data, groups, otolith.var, 
                               base, basetows, mix, mixtows)
  mixinits = oolith.inits(mixdata, type, n.chains, max.prop = max.prob)
  if (is.null(mixparameters)) {
    mixparameters <- otolith.model(mixdata, mixparameters)
  }
  else otolith.model(mixdata, mixparameters)
  bayes <- bugs(mixdata, mixinits, mixparameters, model.file = "mixmodel.txt", 
                n.chains = n.chains, n.iter = n.iter, debug = debug, 
                n.burnin = n.burnin, n.thin = n.thin, n.sims = n.sims, 
                bugs.directory = "c:/Program Files/WinBUGS14/", program = c("WinBUGS", 
                                                                            "OpenBUGS", "winbugs", "openbugs"), working.directory = getwd())
  return(bayes)
}

lionmix=mixy(data=hmix,groups="STOCK", base=mix[mix$SEASON=="R",],mix=mix[mix$SEASON=="C",],type="random", mixparameters = c("P1"),mixtows="H",basetows = "H", otolith.var = c("CARBON","OXYGEN"),n.chains = 2, n.iter = 2000 )

mix=read.csv("../Data/lionmix.csv", header=T)

mixy(mix)
                
data.frame(mix, groups)

                

                
                
                
                
                
                
                
                
                

