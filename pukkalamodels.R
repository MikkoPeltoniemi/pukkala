## Mikko Peltoniemi, Aug 2021, updated Dec 2021
## Pukkala et al., 2021 models
## https://academic.oup.com/forestry/advance-article/doi/10.1093/forestry/cpab008/6172082

library(ggplot2)
library(data.table)
pathtoken = "C:/Users/03011290/OneDrive - Valtion/Active/Projektit/FUNPOTENTIAL/models/pukkala2021/"
paramsI = fread(paste(pathtoken, "pukkala_increment.txt", sep=""))
paramsI = paramsI[, c(1,3, 2, 4)]
paramsS = fread(paste(pathtoken, "pukkala_survival_fixed.txt", sep=""))
paramsIng = fread(paste(pathtoken, "pukkala_ingrowth.txt", sep=""))

I.f = function(X, params=paramsI) {
  ## Note: output is diameter growth vector with estimates for different species, but with same X
  fx = as.matrix(X) %*% as.matrix(params)
  return(exp(fx))
}

S.f = function(X, params=paramsS) {
  ## Note: output is mortality probability for different species in a forest described by X
  fx = as.matrix(X) %*% as.matrix(params)
  return(1/(1+exp(-fx)) )
}

Ing.f = function(X, params=paramsIng) {
  ## Note: output is ingrowth of different species in a forest described by X
  return(exp(as.matrix(X) %*% as.matrix(params)))
}



## Test, Motti Asikkala, 1st simulations, TaimiCO2
dt= fread(paste(pathtoken, "dt.txt", sep=""))
dt = dt[, .(simulointi, vuosi, sp=puulaji, dbh=LPM, Ntrees=Runkoluku, sitetype=2, landclass=1, TS=1400)]
# dtin: 1 sp, 2 dbh, 3 N, 4 sitetype, 5 landclass, 6 TS
dtin = dt[simulointi==1 & vuosi == 17, -c(1:2)]

## PREPARE INPUTS

## override can be used to force forest variables to fixed conditions
## override=c(BAL, BALp, BALs, BALsd, ln_BAplus1, ln_TS)
prepX = function(dtin, override) {
  BALS.f = function(x, dt=dtin) {
    BAL = dt[unlist(x)[2] < dbh, pi*sum(Ntrees*(0.01*dbh/2)^2)]
    BALp = dt[sp==1 & unlist(x)[2] < dbh, pi*sum(Ntrees*(0.01*dbh/2)^2)]
    BALs = dt[sp==2 & unlist(x)[2] < dbh, pi*sum(Ntrees*(0.01*dbh/2)^2)]
    BALsd = dt[(sp>= 2) & unlist(x)[2] < dbh, 
               pi*sum(Ntrees*(0.01*dbh/2)^2)]
    c(BAL=BAL, BALp=BALp, BALs=BALs, BALsd=BALsd) / sqrt(unlist(x)[2]+1)
  }
  BALS = data.table(t(apply(dtin, 1, BALS.f)))
  if (!is.null(override)) {
    BALS[, BAL:=override[1]/sqrt(dtin$dbh+1)]
    BALS[, BALp:=override[2]/sqrt(dtin$dbh+1)]
    BALS[, BALs:=override[3]/sqrt(dtin$dbh+1)]
    BALS[, BALsd:=override[4]/sqrt(dtin$dbh+1)]
  }
  ## T‰ss‰ selitt‰j‰t koottuna Increment-mallille ekaks, lopussa lis‰tty muut
  X = dtin[, .(Intercept=1, #1
               sqrt_dbh=sqrt(dbh), #2
               dbh=dbh, #3
               logBA=log(pi*sum(Ntrees*(dbh/2)^2)+1), #4
               BALLIT=BALS[, -2],#*as.numeric(c(sp==1, sp<=2, sp>=3)), #5,6,7
               logTS=log(TS), #8
               site1=as.numeric(sitetype <= 2), #9
               site2=as.numeric(sitetype == 3), #10
               site3=as.numeric(sitetype >= 4), #11
               peat=as.numeric(landclass >= 2), #12
               dbhinter=dbh*as.numeric(sp %in% 4:5),  #13
           ## Additional for surv model:
           Aspen=as.numeric(sp == 5), #14
           birch=as.numeric(sp %in% 3:4), #15
           BALp=BALS[, 2],#*as.numeric(sp>=3), #16
           ## Additional for ingrowth model
           sqrtG=sqrt(pi*sum(Ntrees*(0.01*dbh/2)^2)), #17
           sqrtGpine= sqrt(pi*sum(Ntrees*(ifelse(sp==1, 0.01*dbh, 0)/2)^2)))]#18
 # colnames(X)[1:13] = unlist(paramsI[, 1])
  if (!is.null(override)) {
    X[, ln_TS:=override[5]]
    X[, ln_BAplus1:=override[6]]
  }
  return(X)
}


dbh_inc.f = function(dtin, params=paramsI[, -1], override) {
  X = prepX(dtin, override)
  delta_dbh=I.f(X[, 1:13], params=params)
  ## filter out correct species estimate:
  ddbh = t(as.matrix(delta_dbh))[t(as.matrix(dtin[, .(sp==1, sp==2, sp>=3)]))]
  #rowMeans(as.matrix(dtin[, .(sp==1, sp==2, sp>=3)]) * 
  #                  as.matrix(delta_dbh))*3
  list(delta=ddbh, X=X)
}

coremodel.f = function(dtin, paramsin=paramsI[, -1], 
                  paramsSin=paramsS[, -1], 
                  paramsIngin=paramsIng[, -1], 
                  override) {
  
  X = prepX(dtin, override)
  
  delta_dbh=I.f(X[, 1:13], params=paramsin)
  ## Tarkista ett‰ t‰m‰n parametrit menee oikein
  surv = S.f(X[, c(1:3, 5, 16, 6, 7, 12, 14, 15)], params = paramsSin)
  
  ing = Ing.f(X[1, c(1, 8, 17, 18, 9, 10, 11)], params = paramsIngin)
  
  ## filter out correct species estimate:
  ddbh = rowMeans(as.matrix(dtin[, .(sp==1, sp==2, sp>=3)]) * 
                    as.matrix(delta_dbh))*3
  list(delta=ddbh, surv=surv, ing=ing, X=X)
}

dtin2 = copy(dtin)
dtin2[, yr:=0]
dbhs = list()
Nsteps = 40
override = NULL 
#dtin2[, dbh:=0.1*dbh]
# 
# override = c(10, 0, 0, 0, log(1400), log(25+1) )
INGROWTH = TRUE
dtin2[, Ntrees:=5]
date()
for (ii in 1:Nsteps) {
  #delta = dbh_inc.f(dtin2[, 1:6], override = override)
  comps = coremodel.f(dtin2[, 1:6], override = override)
  dtin2[, dbh:=dbh+comps$delta]
  dtin2[, ddbh:=comps$delta]
  dtin2[, treeno:=rownames(dtin2)]
  
  ## Assume decrease of Ntrees to allow 
 
  dtin2[, Nout:=Ntrees*ifelse(sp==1, (1-comps$surv[, 1]), 
                                   ifelse(sp==2, (1-comps$surv[, 2]),
                                          (1-comps$surv[, 3])))] 
  dtin2[, Ntrees:=Ntrees-Nout ]#Ntrees*(1-0.015)^5] 
  dtin2[, yr:=5*ii]
  if (INGROWTH) {
    ## Exlude exotic prediction for sp 4 for ingrowth
    dtin2 = rbind(dtin2, data.table(sp=1:3, dbh=0, 
                                    Ntrees=as.numeric(comps$ing)[1:3], 
                                    sitetype=dtin$sitetype[1], 
                                    landclass=dtin$landclass[1], 
                                    TS=dtin$TS[1], 
                                    yr=5*ii,
                                    ddbh = 0, 
                                    treeno=dtin2$treeno[length(dtin2$treeno)], 
                                    Nout = 0) )
  }

  
  dbhs[[ii]] = copy(dtin2)
}

dbhs = rbindlist(dbhs)
G = dbhs[, .(sum(Ntrees * pi*(0.01*dbh/2)^2), sum(Ntrees)) , by=.(yr, sp)]
colnames(G) = c('yr', 'sp','G', 'Ntrees')

date()

## Plottauksia: 
p1 = ggplot(dbhs[treeno %in% c(1, 11, 12, 27, 28, 32, 33)], 
           aes(x=yr, y=dbh, group=treeno, col=factor(sp))) + 
  geom_line()
p1

p2 = ggplot(dbhs[treeno %in% c(1, 10, 20, 27, 28, 32, 33)], 
           aes(x=dbh, y=ddbh, group=treeno, col=factor(sp))) + 
  geom_line()
p2

p3 = ggplot(G, 
           aes(x=yr, y=G, group=factor(sp), col=factor(sp))) + 
  geom_line()
p3

p4 = ggplot(G, 
            aes(x=yr, y=Ntrees, group=factor(sp), col=factor(sp))) + 
  geom_line()
p4

foo = dbhs[, sum(Nout), by=.(yr, sp)]
p5 = ggplot(foo, 
            aes(x=yr, y=V1, group=factor(sp), col=factor(sp))) + 
  geom_line() + ylab(label = "Mortality, N")
p5

foo = dbhs[dbh==1.3, sum(Ntrees), by=.(yr, sp)]
p6 = ggplot(foo, 
            aes(x=yr, y=V1, group=factor(sp), col=factor(sp))) + 
  geom_line() + ylab(label = "ingrowth, N")
p6


library(ggpubr)
## yhdist‰ kuvat yhteen pdf:‰‰n
multi.page <- ggarrange(p1, p2, p3, p4, p5, p6, 
                        nrow = 3, ncol = 2)
ggsave(multi.page, filename = "puk-tests.pdf", 
         hei=13, wid=9)