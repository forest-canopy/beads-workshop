###########################################################################
####### AFEC 2020 Multivariate analysis  ##################################
###########################################################################

## First, remove all the imported data used from previous analyses
rm(list=ls())

## Set working directory
getwd()
#setwd("C:/Users/Aki/Desktop/AFEC2020/Multivariate analysis")

## installing packages

## vegan to run a lot of multivariate analyses
library(vegan)

## ggplot makes nice ordination graphs etc.
library(ggplot2)

# organizing data (aggregatind abundances per plot etc.)
# install.packages("tidyverse")
library(tidyverse) 
library(dplyr)

## This package includes vegan and runs type II tests (adomis.II) and multiple pairwise comparisons##
BiocManager::install("mixOmics")
install.packages("RVAideMemoire")
library(RVAideMemoire) 

## This package calculated beta diversity (total, turnover and nested)
install.packages("betapart")
library(betapart)

## required for ecoCopula analysis and graph generation
#devtools::install_github("gordy2x/ecoCopula")
#devtools::install_github("caijun/ggcorrplot2")
#devtools::install_github("inSileco/graphicsutils")
library(ecoCopula)
library(corrplot)
library(graphicsutils)
library(mvtnorm)
library(sna)


########################################################################################################
### Data preparation ############################
#################################################
## Arboreal ant data collected from rubber plantations and rainforest
## from three locations (bubeng, nabanhe, menglun)
## Each location contains rubber and rainforest habitats
## A total of six plots per habitat 
## Sampling was repeated in wet and dry seasons 
## (3 locations x 2 habitats x 2 seasons x 6 plots = 36 samples)
## e.g., w_Fb2 = wet season, Forest habitat, Bubeng location, 2=plot number
## In each plot, a total 10 trees were randomly selected and collected by canopy tree ant baiting
## Making three datasets (plot-based ants; tree-based ants, and and plot-based env)

## Importing plot based ant data
ants<- read.csv("data.plots.ants.s.csv", row.name=1)
names(ants)
head(ants)

## Removing rare species using specnumber (MARGIN=2 so that count of each species are returned)
## I made a very arbitral decision to define rare species > 2 occurrence (singletons and doubletons
## are removed)
ants.count<-specnumber(ants, MARGIN=2)
ants.many<-ants[,ants.count>2]


## Also importing tree based ant data
ants.trees<-read.csv("data.trees.ants.csv", row.name=1)
names(ants.trees)
head(ants.trees)

## Tree-based data only
## Adding a dummy variable for treebased ant data, as some data contain samples with no ants!
dummy<-c(rep(1,720))
ants.trees$dummy <- dummy

## Log transforming the data (adding 1 by using log1p function so 0s are not undefined)
ants.trees.log<-log1p(ants.trees)

## Importing plot based environmental data
env<- read.csv("data.plots.env1.csv", row.names = 1)
names(env)
head(env)


#### also subsetting the ant data according to the habitats ###
ants.f<-subset(ants, env$forest_type=="forest")
ants.r<-subset(ants, env$forest_type=="rubber")
head(ants.r)

env.f<-subset(env, env$forest_type=="forest")
env.r<-subset(env, env$forest_type=="rubber")
head(env.r)

## Removing species which have only zeros
ants.f<-ants.f[, colSums(ants.f != 0) > 0]
ants.r<-ants.r[, colSums(ants.r != 0) > 0]


## Removing rare species using specnumber (MARGIN=2 so that count of each species are returned)
## I made a very arbitral decision to define rare species > 2 occurrence (singletons and doubletons
## are removed)
ants.count.f<-specnumber(ants.f, MARGIN=2)
ants.many.f<-ants.f[,ants.count.f>2]
ants.count.r<-specnumber(ants.r, MARGIN=2)
ants.many.r<-ants.r[,ants.count.r>2]




################################################################################
#######NMDS ordination using Bray-Curtis function ##############################
###BUT this is nor recommended if the stress value returns too high (0.25)
NMDS.ants.many<-metaMDS(ants.many, distance="bray",autotransform =F)
plot(NMDS.ants.many)
stressplot(NMDS.ants.many)

# We can use the functions `ordiplot` and `orditorp` to add text to the 
# plot in place of points
ordiplot(NMDS.ants.many,type="n")
orditorp(NMDS.ants.many,display="species",col="red",air=0.01) #air = amount of space between labels
orditorp(NMDS.ants.many,display="sites",cex=1.25,air=0.01)

# There are some additional functions that might of interest
ordiplot(NMDS.ants.many,type="n")
ordihull(NMDS.ants.many,groups=env$forest_type,draw="polygon",col="grey90",
         label=FALSE)
orditorp(NMDS.ants.many,display="species",col="red",air=0.01)
orditorp(NMDS.ants.many,display="sites",col=c(rep("green",36),rep("blue",36)),
         air=0.01,cex=1.25)

# Use the function ordisurf to plot contour lines
ordisurf(NMDS.ants.many~elevation, env, main="drawing coutour lines",col="forestgreen", add=T)

########################################################
### Excercise: Let's run NMDS on tree-based ant data! ##






##################################################
#### NMDS plots using ggplot2 ####################
##################################################
## Preparing site scores
data.scores <- as.data.frame(scores(NMDS.ants.many))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- env$for_season  #  add the grp variable created earlier
head(data.scores) #look at the data
unique(data.scores$grp)

## Species score data (to plot the positions of the species on an ordination graph)
species.scores <- as.data.frame(scores(NMDS.ants.many, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

## Buidling convex hull.
## We use chull function to build the data.frame to build the convex hull
grp.a <- data.scores[data.scores$grp == "forest dry", ][chull(data.scores[data.scores$grp == 
                                                                                             "forest dry", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$grp == "rubber dry", ][chull(data.scores[data.scores$grp == 
                                                                                             "rubber dry", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.c <- data.scores[data.scores$grp == "forest wet", ][chull(data.scores[data.scores$grp == 
                                                                                             "forest wet", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.d <- data.scores[data.scores$grp == "rubber wet", ][chull(data.scores[data.scores$grp == 
                                                                                             "rubber wet", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
hull.data

ggplot() + 
  geom_polygon(data=hull.data, aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.30) + # add the convex hulls
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.2) +  # add the species labels
  geom_point(data=data.scores, aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=5) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("forest dry"="red", "rubber dry"="red", "forest wet"="blue", "rubber wet"="blue")) +
  scale_shape_manual(values=c("forest dry"="circle", "rubber dry"="triangle", "forest wet"="circle", "rubber wet"="triangle")) +
  coord_equal() +
  theme_bw()+ 
  theme(#axis.text.x = element_blank(),  # remove x-axis text
    #axis.text.y = element_blank(), # remove y-axis text
    #axis.ticks = element_blank(),  # remove axis ticks
    #axis.title.x = element_text(size=18), # remove x-axis labels
    #axis.title.y = element_text(size=18), # remove y-axis labels
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),  #remove major-grid labels
    panel.grid.minor = element_blank(),  #remove minor-grid labels
    plot.background = element_blank())

###############################
##### Statistical tests #######
###############################

##Defining blocks and the number of permutations (9999)
permAll <- how(nperm = 999)
setBlocks(permAll) <- with(env, location)
permAll
#Testing ALL samples (day 0, 30, 60)
MyM.ants<-adonis.II(ants.many~forest_type*season, 
                          data=env, permutations = permAll, method = "bray") 
MyM.ants

## Pairwise test (if you get interactions or the levels of treatments are more than two)
## avaiable from RVAideMemoire package
MyM.ants.pairwise<-pairwise.perm.manova(vegdist(ants.many, method="bray"), 
                                        env$for_season, nperm=99)
MyM.ants.pairwise

## simper calculates similarity percentage of species that discriminate two groups (habitats)
MyM.ants.simper<-simper(ants.many, env$forest_type, permutations = 99)
summary(MyM.ants.simper)

## bioenv tries to find the best subset of environmental variables that explain the assemblage data
## but this may take time as this attempts to analyse every sinble possible combinations of the 
## env variables
env.num<-dplyr::select_if(env, is.numeric) #selecting only numeric var.
names(env.num)
MyM.ants.bioenv<-bioenv(ants.many, env.num, index="bray")


###Checking dispersion of the data among the forests x seasons ###
MyDM.ants<-betadisper(vegdist(ants.many, method="bray"), env$for_season)
MyDM.ants
anova(MyDM.ants)
permutest(MyDM.ants, pairwise = TRUE, permutations = 999)

## from betapart package 
beta.multi(decostand(ants.many, "pa"))

###########################################
### Summary of the results ################





#############################################################################
####### Below analyses the r'ships between ants vs env using capscale #######
#############################################################################

## We will now investigate what factors are likely to influence ant assemblages 
## We will look at individual habitats separately (forest vs rubber)
## The env factors selected by the below analysis will then be summarized using PCA (PC1 and PC2 etc)
## We will use ordistep function (backward)

#############################
#############################
### Rainforest tree data ####
#############################
#############################
## Looking at env and ant data used for this section
head(ants.many.f)
names(env.f)

## Making dummy variables for location factor
env.f$location.bubeng<-ifelse(env.f$location=='bubeng',1,0)
env.f$location.nabanhe<-ifelse(env.f$location=='nabanhe',1,0)
names(env.f)
head(env.f)

env.f.num<-dplyr::select_if(env.f, is.numeric) #selecting only numeric var.
names(env.f)

####################################################
### pairwise correlation and histograms ############
####################################################

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
#pairs(data.plots.env.num.f, lower.panel= panel.smooth, upper.panel = panel.cor,gap=0)

## put histograms on the diagonal to check the presence of ourliers
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

## making a graph
pairs(env.f.num, lower.panel = panel.smooth,
      cex = 1, pch = 1, bg = "light blue", horOdd=TRUE, 
      upper.panel = panel.cor,       ## Turn on to see the correlation values
      diag.panel = panel.hist, cex.labels = 0.9, font.labels = 2, gap=0)

## Based on the correlation, tree height and trap height are correlated so I will remove trap height
## Others (dbh vs tree height) (elevation vs nabanhe) are highly correlated but retained

## Transforming some of the env variables to achieve more or less normal distributions
names(env.f.num)
env.f.num$height<-log(env.f.num$height)
env.f.num$dbh<-log(env.f.num$dbh)
env.f.num$liana_on_horisontal<-log1p(env.f.num$liana_on_horisontal)
env.f.num$liana_to_crown<-log1p(env.f.num$liana_to_crown)
env.f.num$liana_vertical_on_trunk_cv<-log1p(env.f.num$liana_vertical_on_trunk_cv)
env.f.num$climber_on_trunk<-log1p(env.f.num$climber_on_trunk)
env.f.num$bird_nest_fern<-log1p(env.f.num$bird_nest_fern)
env.f.num$elevation<-log(env.f.num$elevation)
env.f.num$tree_density_include_focal<-log(env.f.num$tree_density_include_focal)
env.f.num$tree_richness_include_focal<-log(env.f.num$tree_richness_include_focal)
env.f.num$Temperature_at_collect<-log(env.f.num$Temperature_at_collect)
env.f.num$wind<-log1p(env.f.num$wind)
env.f.num$canopy_openess<-log1p(env.f.num$canopy_openess)


## Running capscale 
head(ants.many.f)

dbRDA.f<-capscale(ants.many.f~. , env.f.num, distance = "bray", scale=T)
summary(dbRDA.f)
plot(dbRDA.f)
anova(dbRDA.f)
anova(dbRDA.f,  by="axis", perm.max=500)
anova(dbRDA.f,  by="terms", permutations =9999)
dbRDA.f0<-capscale(ants.many.f~ 1, env.f.num, distance = "bray", scale=T)
anova(dbRDA.f, dbRDA.f0)

### Ordistep function to chose the best model. We will first use "backward" to be conservative
dbRDA.f.ordstep<-ordistep(dbRDA.f) 
dbRDA.f.ordstep
plot(dbRDA.f.ordstep) ## default CAP ordination (ugly so we will use ggplot)

## Now we will try stepwise model selection
## Scope defines the upper (maximum) formula where the model selection aims to achieve
## input object (dbRDA0) is set as a lower (starting) scope
dbRDA1.f.ordstep<-ordistep(dbRDA.f0, scope=formula(dbRDA.f))
dbRDA1.f.ordstep
dbRDA1.f.ordstep$anova
plot(dbRDA1.f.ordstep) 

# plotting the results using ordiplot
ordiplot(dbRDA1.f.ordstep,type="n")
orditorp(dbRDA1.f.ordstep,display="sites",col=c(rep("green",18),rep("blue",18)),
         air=0.01,cex=1.25)
# Use the function ordisurf to plot contour lines
ordisurf(dbRDA1.f.ordstep~elevation, env.f.num, main="drawing coutour lines",col="forestgreen", add=T)


##################################################
#### CAP plots using ggplot2 ####################
##################################################

## Preparing site scores
data.scores.cap.f <- as.data.frame(scores(dbRDA1.f.ordstep)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores.cap.f$site <- rownames(data.scores.cap.f)  # create a column of site names, from the rownames of data.scores
data.scores.cap.f$location <- env.f$location  #  add the grp variable created earlier
head(data.scores.cap.f) #look at the data
unique(data.scores.cap.f$location)


## Building convex hull for location differences
## We use chull function to build the data.frame to build the convex hull

grp.m <- data.scores.cap.f[data.scores.cap.f$location == "menglun", ][chull(data.scores.cap.f[data.scores.cap.f$location == 
                                                                                                "menglun", c("CAP1", "CAP2")]), ]  # hull values for grp Menglun
grp.bb <- data.scores.cap.f[data.scores.cap.f$location == "bubeng", ][chull(data.scores.cap.f[data.scores.cap.f$location == 
                                                                                                "bubeng", c("CAP1", "CAP2")]), ]  # hull values for grp Bubeng
grp.n <- data.scores.cap.f[data.scores.cap.f$location == "nabanhe", ][chull(data.scores.cap.f[data.scores.cap.f$location == 
                                                                                                "nabanhe", c("CAP1", "CAP2")]), ]  # hull values for grp Nabanhe

hull.data.l <- rbind(grp.m, grp.bb, grp.n)  #combine grp.a and grp.b
hull.data.l

## Biplot vectors
vec.df<-as.data.frame(dbRDA1.f.ordstep$CCA$biplot)
vec.df$labels <- rownames(vec.df)  # create a column of site names, from the rownames of data.scores
head(vec.df)  #look at the data

ggplot() + 
  geom_polygon(data=hull.data.l, aes(x=CAP1, y=CAP2, fill=location, group=location), alpha=0.30) + # add the convex hulls
  geom_point(data=data.scores.cap.f, aes(x=CAP1,y=CAP2,shape=location,colour=location),size=5) + # add the point markers
  scale_colour_manual(values=c("menglun"="red", "nabanhe"="blue", "bubeng"="green")) +
  ## Adding biplot vectors of the selected env variables
  geom_segment(data = vec.df, 
               aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour="blue",
               inherit.aes = T) + 
  geom_text(data = vec.df, 
            aes(x = CAP1, y=CAP2, label = labels),
            size=4) +
  geom_text(data=data.scores.cap.f,aes(x=CAP1,y=CAP2,label=site),size=4,vjust=0) +  # add the site labels
  
  coord_equal() +
  theme_bw()+ 
  theme(#axis.text.x = element_blank(),  # remove x-axis text
    #axis.text.y = element_blank(), # remove y-axis text
    #axis.ticks = element_blank(),  # remove axis ticks
    #axis.title.x = element_text(size=18), # remove x-axis labels
    #axis.title.y = element_text(size=18), # remove y-axis labels
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),  #remove major-grid labels
    panel.grid.minor = element_blank(),  #remove minor-grid labels
    plot.background = element_blank())


######### running PCA of the selected env variables #############
dbRDA1.f.ordstep
names(env.f.num)
env.f.num.pca<-env.f.num %>%
  dplyr::select(elevation , tree_density_include_focal , ant_species ,
                  climber_on_trunk_mean , height_mean , bird_nest_fern_cv , location.bubeng , liana_to_crown_cv ,
                  height_cv)
PCA.f <-rda(env.f.num.pca, scale = T)
summary(PCA.f)
data.rda.coord <- as.data.frame(PCA.f$CA$u[,1:2]) ## extracting PCA1 and PCA2
data.rda.coord

## running capscale using PCA1 and PCA2 as explanatory variables
dbRDA.df.pca<-capscale(ants.many.f~. , data.rda.coord, distance = "bray", scale=T)
## Checking if this model explains sig variation in ant assemblage (and they do)
anova(dbRDA.df.pca)
anova(dbRDA.df.pca,  by="terms", permutations =9999)

