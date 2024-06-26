# PGLS-models
R code of PGLS models testing AMR against life-history and ecological variables

## Preprint
The following code has been developed in collaboration with **Jose Valdebenito**. The methodology underlying the use of the scripts has been detailed in the preprint _"Urbanization spreads antimicrobial resistant enteric pathogens in wild bird microbiomes"_. Please refer to citation information at the bottom of this document.

## Files Summary



## R Packages
* phytools
* caper
* geiger
* knitr

## Load tree phylogeny
```
t <- read.tree("/path/to/directory/tree")
plotTree(t, ftype="i")
```
<img width="518" alt="1 tree" src="https://github.com/evangelosmourkas/PGLS-models/assets/73548463/920cd3d8-9ffb-4cf9-b224-ff6f4f96e1a9">

## Load data
```
df <- read.csv("/path/to/directory/data.csv")
df$Weights... <- NULL
df$Common.name <- NULL
colnames(df) <- c("species", "number", "lin", "amr", "prop","b.mass","clutch_size", "msys", "migra", "habit", "diet", "social")

df$b.mass <- log(df$b.mass)
df$habit <- relevel(factor(df$habit), ref="terrestrial")
df$migra <- relevel(factor(df$migra), ref="no")
df$diet <- relevel(factor(df$diet), ref="omnivorous")
df$social <- as.factor(df$social)
df$msys <- as.factor(df$msys)

df <- df[match(t$tip.label, df$species),]

set.seed(4)
```

Then we prepare for running the phylogenetically corrected least squares (PGLS). It will be four sets of models, one for each measure of AMR:

a.Shannon diversity of AMR
b.Shannon diversity of Clonal complexes (CC)
c.The proportion of AMR genes present per host species, excluding genes against beta-lactam drugs

```
comp1 <- comparative.data(data=df, phy=t, species, vcv = T, vcv.dim = 3)
```

a.Shannon diversity of AMR
Run all models, one after the other [(AMR parameters and Ecological variables (N=30)]
```
pgls1 <- pgls(amr ~ social, comp1, lambda = "ML")
pgls2 <- pgls(amr ~ msys, comp1, lambda = "ML")
pgls3 <- pgls(amr ~ migra, comp1, lambda = "ML")
pgls4 <- pgls(amr ~ habit, comp1, lambda = "ML")
pgls5 <- pgls(amr ~ diet+b.mass, comp1, lambda = "ML")
pgls6 <- pgls(amr ~ b.mass, comp1, lambda = "ML")
pgls7 <- pgls(amr ~ clutch_size, comp1, lambda = "ML")

pgls.res <- round(rbind(
  summary(pgls1)$coefficients,
  summary(pgls2)$coefficients,
  summary(pgls3)$coefficients,
  summary(pgls4)$coefficients,
  summary(pgls5)$coefficients,
  summary(pgls6)$coefficients,
  summary(pgls7)$coefficients),3)

pgls.res
```
We continue with remaining three set of models:
b.Shannon diversity of Clonal complexes (CC)
c.The proportion of AMR genes present per host species, excluding genes against beta-lactam drugs
```
pgls1.1 <- pgls(lin ~ social, comp1, lambda = "ML")
pgls2.1 <- pgls(lin ~ msys, comp1, lambda = "ML")
pgls3.1 <- pgls(lin ~ migra, comp1, lambda = "ML")
pgls4.1 <- pgls(lin ~ habit, comp1, lambda = "ML")
pgls5.1 <- pgls(lin ~ diet+b.mass, comp1, lambda = "ML")
pgls6.1 <- pgls(lin ~ b.mass, comp1, lambda = "ML")
pgls7.1 <- pgls(lin ~ clutch_size, comp1, lambda = "ML")

pgls.res1 <- round(rbind(
  summary(pgls1.1)$coefficients,
  summary(pgls2.1)$coefficients,
  summary(pgls3.1)$coefficients,
  summary(pgls4.1)$coefficients,
  summary(pgls5.1)$coefficients,
  summary(pgls6.1)$coefficients,
  summary(pgls7.1)$coefficients),3)

pgls1.2 <- pgls(prop ~ social, comp1, lambda = "ML")
pgls2.2 <- pgls(prop ~ msys, comp1, lambda = "ML")
pgls3.2 <- pgls(prop ~ migra, comp1, lambda = "ML")
pgls4.2 <- pgls(prop ~ habit, comp1, lambda = "ML")
pgls5.2 <- pgls(prop ~ diet+b.mass, comp1, lambda = "ML")
pgls6.2 <- pgls(prop ~ b.mass, comp1, lambda = "ML")
pgls7.2 <- pgls(prop ~ clutch_size, comp1, lambda = "ML")

pgls.res2 <- round(rbind(
  summary(pgls1.2)$coefficients,
  summary(pgls2.2)$coefficients,
  summary(pgls3.2)$coefficients,
  summary(pgls4.2)$coefficients,
  summary(pgls5.2)$coefficients,
  summary(pgls6.2)$coefficients,
  summary(pgls7.2)$coefficients),3)

pgls.res1
```
Lineage diversity
```
pgls.res2
```
Proportion of AMR

Proceed to do the conduct the part of the analysis involving the proximity to humans estimates
AMR parameters and Proximity to urbanisation (N=38)
amr = AMR diversity
lin = Lineage diversity
prop = Proportion of resistant isolates (excluding β-lactam genes)

```
tt <- read.tree("/path/to/directory/tree_submit_add_branch")
#loading three 38 tips. this one has sister branches of repeated species of length <0.00001.
tt <- drop.tip(tt, tip="Larus_californicus")
tt
```
 
```
pgd <- read.csv("/path/to/directory/prox.csv", header = T)
pgd$Weights... <- NULL
colnames(pgd)[1:8] <- c("species", "name", "country", "prox", "number", "lin", "amr", "prop")
comp <- comparative.data(tt, pgd, species, vcv=T,vcv.dim = 3)

pg1 <- pgls(amr ~ prox, comp, lambda = "ML")
pg2 <- pgls(lin ~ prox, comp, lambda = "ML")
pg3 <- pgls(plogis(prop) ~ prox, comp, lambda = "ML")
summary(pg1)
```
```
summary(pg2)
```
```
summary(pg3)
```
AMR parameters and proximity to humans, subset excluding genomes with n<2 (N=20)
```
pgd0 <- pgd[pgd$number!=1,]
comp0 <- comparative.data(tt, pgd0, species, vcv=T,vcv.dim = 3)
pg10 <- pgls(amr ~ prox, comp0, lambda = "ML")
pg20 <- pgls(lin ~ prox, comp0, lambda = "ML")
pg30 <- pgls(plogis(prop) ~ prox, comp0, lambda = "ML")
summary(pg10)
```

```
summary(pg20)
```
```
summary(pg30)
```
# How to cite
Mourkas E, Valdebenito JO, Marsh H, Hitchings MD, Cooper KK, Parker CT, Székely T, Johansson H, Ellström P, Pascoe B, Waldenström J, Sheppard SK (2023) **Urbanisation spreads antimicrobial resistant enteric pathogens in wild bird microbiomes**.
bioRxiv doi: 10.1101/2023.07.11.548564
