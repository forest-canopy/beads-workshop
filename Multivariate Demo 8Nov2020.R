###########################################################################
####### AFEC 2020 Multivariate analysis  ##################################
###########################################################################

## First, remove all the imported data used from previous analyses
rm(list=ls())

## Set working directory
getwd()
#setwd("C:/Users/Aki/Desktop/AFEC2020/Multivariate analysis")
setwd("test")

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
gr3.sp<- read.csv("./gr3/group_3.csv", row.name=1)
names(gr3.sp)
head(gr3.sp)

## Removing rare species using specnumber (MARGIN=2 so that count of each species are returned)
## I made a very arbitral decision to define rare species > 2 occurrence (singletons and doubletons
## are removed)
count<-specnumber(gr3.sp, MARGIN=2)
gr3.sp.norare<-gr3.sp[,count>2]


## Adding a dummy variable for treebased ant data, as some data contain samples with no ants!
dummy<-c(rep(1,720)) # replace 720 with the number of YOUR samples
gr3.sp$dummy <- dummy

## Log transforming the data (adding 1 by using log1p function so 0s are not undefined)
gr3.sp.log<-log1p(gr3.sp)

## Importing environmental data
gr3.env<- read.csv("./gr3/group_3(elevation)(1).csv", row.names = 1)
names(gr3.env)
head(gr3.env)


#### also subsetting the ant data according to the habitats ###
gr3.sp.f<-subset(gr3.sp, gr3.env$forest_type=="forest")
gr3.sp.s<-subset(gr3.sp, gr3.env$forest_type=="secondary")
head(gr3.sp.f)

gr3.env.f<-subset(gr3.env, gr3.env$forest_type=="forest")
gr3.env.s<-subset(gr3.env, gr3.env$forest_type=="secondary")
head(gr3.env.s)

## Removing species which have only zeros
gr3.sp.f<-gr3.sp.f[, colSums(gr3.sp.f != 0) > 0]
gr3.sp.s<-gr3.sp.s[, colSums(gr3.sp.s != 0) > 0]


## Removing rare species using specnumber (MARGIN=2 so that count of each species are returned)
## I made a very arbitral decision to define rare species > 2 occurrence (singletons and doubletons
## are removed)
count.f<-specnumber(gr3.sp.f, MARGIN=2)
gr3.sp.f<-gr3.sp.f[,count.f>2]
count.s<-specnumber(gr3.sp.s, MARGIN=2)
gr3.sp.s<-gr3.sp.s[,count.s>2]




################################################################################
#######NMDS ordination using Bray-Curtis function ##############################
###BUT this is nor recommended if the stress value returns too high (0.25)
NMDS.gr3.raw<-metaMDS(gr3.sp, distance="bray",autotransform =F)
plot(NMDS.gr3.raw)
stressplot(NMDS.gr3.raw)

# We can use the functions `ordiplot` and `orditorp` to add text to the 
# plot in place of points
ordiplot(NMDS.gr3.raw,type="n")
orditorp(NMDS.gr3.raw,display="species",col="red",air=0.01) #air = amount of space between labels
orditorp(NMDS.gr3.raw,display="sites",cex=1.25,air=0.01)

# There are some additional functions that might of interest
ordiplot(NMDS.gr3.raw,type="n")
ordihull(NMDS.gr3.raw,groups=gr3.env$forest_type,draw="polygon",col="grey90",
         label=FALSE)
orditorp(NMDS.gr3.raw,display="species",col="red",air=0.01)
orditorp(NMDS.gr3.raw,display="sites",col=c(rep("green",36),rep("blue",36)),
         air=0.01,cex=1.25)

# Use the function ordisurf to plot contour lines
ordisurf(NMDS.gr3.raw~elevation, gr3.env, main="drawing coutour lines",col="forestgreen", add=T)

########################################################
### Excercise: Let's run NMDS on tree-based ant data! ##






##################################################
#### NMDS plots using ggplot2 ####################
##################################################
## Preparing site scores
data.scores <- as.data.frame(scores(NMDS.gr3.raw))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # creat column of site names, from the rownames of data.scores
data.scores$grp <- gr3.env$forest_type  #  add the grp variable created earlier
head(data.scores) #look at the data
unique(data.scores$grp)

## Species score data (to plot the positions of the species on an ordination graph)
species.scores <- as.data.frame(scores(NMDS.gr3.raw, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

## Buidling convex hull.
## We use chull function to build the data.frame to build the convex hull
grp.a <- data.scores[data.scores$grp == "primary", ][chull(data.scores[data.scores$grp == 
                                                                         "primary", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$grp == "secondary", ][chull(data.scores[data.scores$grp == 
                                                                           "secondary", c("NMDS1", "NMDS2")]), ]  # hull values for grp B                                                                                    "rubber wet", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
hull.data

ggplot() + 
  geom_polygon(data=hull.data, aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.30) + # add the convex hulls
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.2) +  # add the species labels
  geom_point(data=data.scores, aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=5) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("primary"="blue", "secondary"="red")) +
  scale_shape_manual(values=c("primary"="circle", "secondary"="triangle")) +
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

##Defining blocks and the number of permutations (999)
permAll <- how(nperm = 999)
setBlocks(permAll) <- with(gr3.env, location)
permAll
#Testing ALL samples (day 0, 30, 60)
MyM.gr3<-adonis.II(gr3.sp~soilN*elevation*primary, 
                   data=gr3.env.num, permutations = 999, method = "bray") 
MyM.gr3

## Pairwise test (if you get interactions or the levels of treatments are more than two)
## avaiable from RVAideMemoire package
MyM.ants.pairwise<-pairwise.perm.manova(vegdist(ants.many, method="bray"), 
                                        env$for_season, nperm=99)
MyM.ants.pairwise

## simper calculates similarity percentage of species that discriminate two groups (habitats)
MyM.gr3.simper<-simper(gr3.sp, gr3.env$forest_type, permutations = 99)
summary(MyM.gr3.simper)

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
beta.multi(decostand(gr3.sp, "pa"))

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
gr3.env$primary<-ifelse(gr3.env$forest_type=='primary',1,0)
env.f$location.nabanhe<-ifelse(env.f$location=='nabanhe',1,0)
names(gr3.env)
head(gr3.env)

gr3.env.num<-dplyr::select_if(gr3.env, is.numeric) #selecting only numeric var.
names(gr3.env.num)

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
pairs(gr3.env.num, lower.panel = panel.smooth,
      cex = 1, pch = 1, bg = "light blue", horOdd=TRUE, 
      upper.panel = panel.cor,       ## Turn on to see the correlation values
      diag.panel = panel.hist, cex.labels = 0.9, font.labels = 2, gap=0)

## Based on the correlation, tree height and trap height are correlated so I will remove trap height
## Others (dbh vs tree height) (elevation vs nabanhe) are highly correlated but retained

## Transforming some of the env variables to achieve more or less normal distributions
names(env.f.num)
gr3.env.num$soilN<-log1p(gr3.env.num$soilN)
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

dbRDA.gr3<-capscale(gr3.sp~. , gr3.env.num, distance = "bray", scale=T)
summary(dbRDA.gr3)
plot(dbRDA.gr3)
anova(dbRDA.f)
anova(dbRDA.f,  by="axis", perm.max=500)
anova(dbRDA.f,  by="terms", permutations =9999)
dbRDA.gr30<-capscale(gr3.sp~ 1, gr3.env.num, distance = "bray", scale=T)
anova(dbRDA.gr3, dbRDA.gr30)

### Ordistep function to chose the best model. We will first use "backward" to be conservative
dbRDA.f.ordstep<-ordistep(dbRDA.f) 
dbRDA.f.ordstep
plot(dbRDA.f.ordstep) ## default CAP ordination (ugly so we will use ggplot)

## Now we will try stepwise model selection
## Scope defines the upper (maximum) formula where the model selection aims to achieve
## input object (dbRDA0) is set as a lower (starting) scope
dbRDA.gr3.ordstep<-ordistep(dbRDA.gr30, scope=formula(dbRDA.gr3))
dbRDA.gr3.ordstep
dbRDA.gr3.ordstep$anova
plot(dbRDA.gr3.ordstep) 

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

