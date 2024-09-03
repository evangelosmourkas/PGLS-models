# PGLS-models
R code of PGLS models testing AMR against life-history and ecological variables

## Publication
The following code has been developed in collaboration with **Jose Valdebenito**. The methodology underlying the use of the script has been detailed in the article "Proximity to humans is associated with antimicrobial-resistant enteric pathogens in wild bird microbiomes" published in Current Biology.

## Files Summary

## R Packages
* phytools
* caper
* geiger
* knitr

## Load tree phylogeny
```
t <- read.tree("/path/to/directory/tree_submit")
t <- drop.tip(t, tip="Larus_californicus")
plotTree(t, ftype="i")
```
<img width="518" alt="1 tree" src="https://github.com/evangelosmourkas/PGLS-models/assets/73548463/920cd3d8-9ffb-4cf9-b224-ff6f4f96e1a9">

## Load data
```
df <- read.csv("/path/to/directory/Table S7_New ecol.csv")
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
<img width="344" alt="2 pgls res" src="https://github.com/evangelosmourkas/PGLS-models/assets/73548463/8da373ae-7823-4b7d-ab92-fbe937196006">

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

<img width="339" alt="3 pgls res1" src="https://github.com/evangelosmourkas/PGLS-models/assets/73548463/ddbbd96b-093c-4bdd-86e4-3de490244726">

```
pgls.res2
```
Proportion of AMR

<img width="340" alt="4 pgls res2" src="https://github.com/evangelosmourkas/PGLS-models/assets/73548463/933d58f2-dfa4-4096-ad91-e311fdd6df8e">

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
<img width="617" alt="5 tt" src="https://github.com/evangelosmourkas/PGLS-models/assets/73548463/dc1a372e-8835-4ba8-a6c1-62f4dc7f97e2">
 
```
pgd <- read.csv("/path/to/directory/Table S7_New prox.csv", header = T)
pgd$Weights... <- NULL
colnames(pgd)[1:8] <- c("species", "name", "country", "prox", "number", "lin", "amr", "prop")
comp <- comparative.data(tt, pgd, species, vcv=T,vcv.dim = 3)

pg1 <- pgls(amr ~ prox, comp, lambda = "ML")
pg2 <- pgls(lin ~ prox, comp, lambda = "ML")
pg3 <- pgls(plogis(prop) ~ prox, comp, lambda = "ML")
summary(pg1)
```
<img width="367" alt="6 pg1" src="https://github.com/evangelosmourkas/PGLS-models/assets/73548463/9f9ac200-2c7a-4445-a956-6be2ffd09764">

```
summary(pg2)
```

<img width="360" alt="7 pg2" src="https://github.com/evangelosmourkas/PGLS-models/assets/73548463/5eaa5d6c-b230-45ca-b272-771b601e64a5">

```
summary(pg3)
```

<img width="382" alt="8 pg3" src="https://github.com/evangelosmourkas/PGLS-models/assets/73548463/2416ab89-9d68-4c68-ac1a-edac25ba5980">

AMR parameters and proximity to humans, subset excluding genomes with n<2 (N=20)

```
pgd0 <- pgd[pgd$number!=1,]
comp0 <- comparative.data(tt, pgd0, species, vcv=T,vcv.dim = 3)
pg10 <- pgls(amr ~ prox, comp0, lambda = "ML")
pg20 <- pgls(lin ~ prox, comp0, lambda = "ML")
pg30 <- pgls(plogis(prop) ~ prox, comp0, lambda = "ML")
summary(pg10)
```

<img width="363" alt="9 pg10" src="https://github.com/evangelosmourkas/PGLS-models/assets/73548463/c22504d2-a225-4cb5-be6d-7f9a19ffe6aa">

```
summary(pg20)
```

<img width="373" alt="10 pg20" src="https://github.com/evangelosmourkas/PGLS-models/assets/73548463/6e25b923-21fb-42d5-a775-62b0a9af15a7">

```
summary(pg30)
```

<img width="375" alt="11 pg30" src="https://github.com/evangelosmourkas/PGLS-models/assets/73548463/099fb0d8-85a7-403d-a026-868c7cbc4a1f">

# How to cite
Mourkas E, Valdebenito JO, Marsh H, Hitchings MD, Cooper KK, Parker CT, Székely T, Johansson H, Ellström P, Pascoe B, Waldenström J, Sheppard SK (2023) **Proximity to humans is associated with antimicrobial resistant enteric pathogens in wild bird microbiomes**.
Current Biology. doi: https://doi.org/10.1016/j.cub.2024.07.059
