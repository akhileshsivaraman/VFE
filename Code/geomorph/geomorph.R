##### geomorph #####
setwd("~/Documents/Imperial/UEBS/geomorph/")
library(geomorph) # mac users will need to download XQuartz

#### Adams & Otarola-Castillo, 2013 ####
data("plethodon")

# plot original landmark data
plotAllSpecimens(plethodon$land)

# Procrustes superimposition and plot the aligned landmarks
Y.gpa <- gpagen(plethodon$land) # generalised procrustes analysis of 2D/3D landmark data
plotAllSpecimens(Y.gpa$coords, links = plethodon$links)

# 3D data with semilandmarks
data("scallops")
Y.gpa2 <- gpagen(A = scallops$coorddata, # object containing landmark coordinates
                 curves = scallops$curvslide, # landmarks to be treated as semilandmarks on boundary curves
                 surfaces = scallops$surfslide # landmarks to be treated as semilandmarks on surfaces
)
plotAllSpecimens(Y.gpa2$coords) # 3D plot


# principal coordinate analysis
# plotTangentSpace() removed from geomorph
# gm.prcomp, plot.gm.prcomp (doesn't exist) and picknplot.shape are now used
pca <- gm.prcomp(Y.gpa$coords) # PCA analysis on Procrustes coordinates
pc <- as.data.frame(pca$x)
plot(pc$Comp1,
     pc$Comp2) # not sure about this

# in this data set specimens represent two species in two distinct environments
# we can use MANOVA to identify sig differences
y <- two.d.array(Y.gpa$coords) # convert a 3D array into a 2D matrix
lm <- procD.lm(y~plethodon$species*plethodon$site, iter=99)

# visualise mean shape difference between groups
# obtain average landmark coordinates for each group and the overall mean
# plot differences as thin plate spline transformation grids
ref <- mshape(Y.gpa$coords)
gp1.mn <- mshape(Y.gpa$coords[,,1:20])
plotRefToTarget(ref, gp1.mn, mag = 2,
                links = plethodon$links)


# multivariate patterns of allometry can also be visualised
ratgpa <- gpagen(ratland)
plotAllometry(ratgpa$coords, ratgpa$Csize, method = "CAC")
# ratgpa is an atomic vector $ is invalid but it works elsewhere

# multivariate regression
ratlm <- procD.lm(two.d.array(ratgpa$coords)~ratgpa$Csize, iter = 999)


# combining phylogenetic data with shape data to estimate the degree of phylogenetic signal in shape
# view of the shape space with the phylogeny superimposed
data("plethspecies")
plethgpa <- gpagen(plethspecies$land)
a <- physignal(A = plethgpa$coords,
          phy = plethspecies$phy,
          iter = 99)
a$phy.signal
a$pvalue

plot(a$PaCA$phy)



#### geomorph digitize ####
### digitizing fixed landmarks
data("scallopPLY")
my.ply <- scallopPLY$ply # scallop 3D image in ply format

# read a ply file with read.ply()
# digit.fixed() to digitize fixed landmarks
fixed.lms1 <- digit.fixed(spec = my.ply, # object containing 3D vertices
                          fixed = 5 # number of landmarks to digitize
                          )
# then click on the landmarks
# for mac users, change preferences: preferences => input => emulate three button mouse
# click to rotate the image
# hold command and click to select landmarks
# new .nts file in wd

clear3d() # clear rgl plot

### sampling surface semilandmarks
surf.pts1 <- buildtemplate(spec = my.ply,
                           fixed = fixed.lms1, # matrix of coordinates selected earlier
                           surface.sliders = 100)
# build a reference template - semilandmarks in the reference are used to sample homologues in other specimens
# saves a .txt file with the coordinates of the landmarks
# the .nts file is modified too to contrain the coordinates of the landmarks and semilandmarks

# can also import landmarks saved in an .nts file with readland.nts()

# sample other specimens for homologous landmarks using digitsurface()


#### geomorph assitance ####
data("larvalMorph")
larvagpa <- gpagen(larvalMorph$tailcoords,
                   curves = larvalMorph$tail.sliders,
                   ProcD = F,
                   print.progress = F)
plot(larvagpa)

# create a data frame list
gdf <- geomorph.data.frame(larvagpa,
                           treatments = larvalMorph$treatment,
                           family = larvalMorph$family)
# create models
fit.size <- procD.lm(coords ~ log(Csize),
                     data = gdf, print.progress = F) # simple allometry model
fit.family <- procD.lm(coords ~ log(Csize) * family, 
                      data = gdf, print.progress = FALSE) # unique family allometries
fit.treatment <- procD.lm(coords ~ log(Csize)*treatments/family, 
                         data = gdf, print.progress = FALSE) # unique treatment: family allometries
## anova
anova(fit.size)
anova(fit.family)
anova(fit.treatment)
anova(fit.size, fit.family, fit.treatment) # compare models

## allometry analyses
plot(fit.size,
     type = "regression", # plot multivariate dispersion against an explanatory variable
     reg.type = "PredLine", # prediction line or a regression line
     predictor = log(gdf$Csize))
plot(fit.size,
     type = "regression", 
     reg.type = "RegScore", 
     predictor = log(gdf$Csize))

plotAllometry(fit.size, size = gdf$Csize, logsz = T, method = "PredLine")
plotAllometry(fit.size, size = gdf$Csize, logsz = T, method = "RegScore")
# plots are model-based projections of shape data => change the model and the plot changes


## two-block partial least sqaures analysis
# not based on a particular model
pls <- two.b.pls(log(gdf$Csize), gdf$coords)
pls
plot(pls)
# same as a common allometric component plot
plotAllometry(fit.size, size = gdf$Csize, logsz = T, method = "CAC")


### more complex allometry models
# family and treatment have significant effects
# do family and treatment have unique allometries or common allometries
fit.unique <- procD.lm(coords ~ log(Csize) * treatments/family,
                       data = gdf)
fit.common <- procD.lm(coords ~ log(Csize) + treatments/family,
                       data = gdf)
anova(fit.common, fit.unique)
# models are not significantly different => common allometry 

# plot by colouring points according to treatment
plotAllometry(fit.common,
              size = gdf$Csize, 
              logsz = T,
              method = "PredLine",
              pch = 1,
              col = as.numeric(gdf$treatments))
plotAllometry(fit.common,
              size = gdf$Csize, 
              logsz = T,
              method = "RegScore",
              pch = 1,
              col = as.numeric(gdf$treatments))

### more anova and pairwise comparisons
# treatment is a fixed effect
# family is a random effect
anova(fit.common)
# F-value is by default MS treatment/MS residuals
# in our mixed model anova we want to see if the treatment effect is meaningful among families
anova(fit.common, error = c("Residuals", "treatments:family", "Residuals"))

# shape covaries with size in a common way for each treatment
# let's compare treatment least-squares means to see with treatments differ in shape, accounting for allometry and for family effects
# use pairwise(fit, fit.null, groups, covariate)
# need to find null model first
reveal.model.designs(fit.common) # returns every full and reduced model for model terms in a linear model
# reduced model is the proposed null model

fit.null <- procD.lm(coords ~ log(Csize) + family,
                     data = gdf) # null model could also be log(Csize) + family
pw <- pairwise(fit.common, fit.null, groups = gdf$treatments)
pw

## tests
# distance between least squares means
summary(pw, test.type = "dist")
anova(fit.null, fit.common)
summary(pw, test.type = "dist", stat.table = F) # generate pairwise tables for each statistic

# disparity (variance) among treatments
summary(pw, test.type = "var")
morphol.disparity(fit.common, groups = gdf$treatments) # morphological disparity and pairwise comparisons among groups
# same output




