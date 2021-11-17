#----------------------------------------------------------------------------#
# Title: FC_vis
#
# Description: This script generates the figures presented at the 
#              ESA 2021 conference using the Vizumap R package. For additional
#              help using the package, see the vignette at 
#              https://lydialucchesi.github.io/Vizumap/
# Author: Petra Kuhnert
# Date: 16 November 2021.
#----------------------------------------------------------------------------#

#------------------------- loading relevant libraries -----------------------#

# Installation of Vizumap
# remotes::install_github(repo = "lydialucchesi/Vizumap", build_vignettes = TRUE, force = TRUE)

library(Vizumap) # Visualising uncertainty in spatial data
library(raster)  # shapefiles
library(ggplot2)


#------------------------- load datasets -------------------------------------#

load("data/bowen_shp.rda") # Bowen shapefile
load("data/bowen_wet.rda") # Wet season data for 18/19 and 19/20 wet seasons
head(BB_1819wet)
head(BB_1920wet)

# mapping limits to ensure comparison between maps. This can be produced
# from the findNbounds function.
load("data/bound.rda")   

#------------------------- Format the data for Vizumap ----------------------#

fc_1819 <- read.uv(data = BB_1819wet, estimate = "mean", error = "std")
fc_1920 <- read.uv(data = BB_1920wet, estimate = "mean", error = "std")

#------------------------- Bivariate Map ---------------------------#

# Create a blue/yellow palette.  There are other palettes available or
# you can create your own.
fc_pal <- build_palette(name = "BlueYellow")

#--------------- 2018/19 wet season map
fcBivMap1819 <- build_bmap(data = fc_1819, geoData = bowen_shp, 
                           id = "SUBCATS", terciles = TRUE,
                           palette = fc_pal, bound = bound) 
fcBivMap1819_p <- view(fcBivMap1819)
m1819 <- fcBivMap1819_p + ggtitle("2018/19 Fractional Cover (Dec-Feb)") + coord_fixed()
m1819

# build a key
fcBivKey1819 <- build_bkey(data = fc_1819, palette = fc_pal, bound = bound)
view(fcBivKey1819)

# attach key
attach_key(fcBivMap1819, fcBivKey1819)

#--------------- 2019/20 wet season map
fcBivMap1920 <- build_bmap(data = fc_1920, geoData = bowen_shp, 
                           id = "SUBCATS", terciles = TRUE,
                           palette = fc_pal, bound = bound) 
fcBivMap1920_p <- view(fcBivMap1920)
m1920 <- fcBivMap1920_p + ggtitle("2019/20 Fractional Cover (Dec-Feb)") + 
             coord_fixed()
m1920

# build a key
fcBivKey1920 <- build_bkey(data = fc_1920, palette = fc_pal, bound = bound)
view(fcBivKey1920)

# attach key
attach_key(fcBivMap1920, fcBivKey1920)

#------------------------- Map Pixelation ---------------------------#

# pixelate the bowen catchment (this may take a few moments)
bowen_pix <- pixelate(bowen_shp, id = "region")

# matching pixels with polygons
df <- data.frame(region = sapply(slot(bowen_shp, "polygons"),
                                 function(x) slot(x, "ID")), name = unique(bowen_shp@data$SUBCATS))


#--------------- 2018/19 wet season map

# Building the pixelate map based on random draws from a distribution
# We choose uniform for illustration but others could be considered.
fc_1819$region <- df[match(fc_1819$SUBCATS, df$name), 1]
fc_1819$region <- as.character(fc_1819$region)

fcPix1819 <- build_pmap(data = fc_1819, distribution = "uniform", 
                pixelGeo = bowen_pix, id = "region", palette = "Greens",
                border = bowen_shp)
view(fcPix1819) + ggtitle("2018/19 Fractional Cover (Dec-Feb)") + coord_equal()

# animate
fcPixA1819 <- Vizumap::animate(fcPix1819, aniLength = 30)
view(fcPixA1819)


#--------------- 2019/20 wet season map

# Building the pixelate map based on random draws from a distribution
# We choose uniform for illustration but others could be considered.
fc_1920$region <- df[match(fc_1920$SUBCATS, df$name), 1]
fc_1920$region <- as.character(fc_1920$region)

fcPix1920 <- build_pmap(data = fc_1920, distribution = "uniform", 
                        pixelGeo = bowen_pix, id = "region", palette = "Greens",
                        border = bowen_shp)
view(fcPix1920) + ggtitle("2019/20 Fractional Cover (Dec-Feb)") + coord_equal()

# animate
fcPixA1920 <- Vizumap::animate(fcPix1920, aniLength = 30)
view(fcPixA1920)

#------------------------- Glyph Rotation ---------------------------#

# Need to change projection so coordinate system is in lats/longs
proj4string(bowen_shp) <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +ellps=WGS84") 
bowen_utm <- spTransform(bowen_shp, CRS("+proj=longlat +datum=WGS84"))

# find limits
glyph_estlimits <- range(fc_1819$mean, fc_1920$mean)
glyph_errmax <- max(fc_1819$std, fc_1920$std)

#--------------- 2018/19 wet season map

bowen_gm1819 <- build_gmap(data = fc_1819, geoData = bowen_utm, id = "SUBCATS", size = 1, 
                       glyph = "icone", palette = "Greens", 
                       limits = glyph_estlimits, max_error = glyph_errmax)
view(bowen_gm1819) + coord_equal() 

#build a glyph key 
fcglyphkey <- build_gkey(data = rbind(fc_1819, fc_1920) , glyph = "icone")
view(fcglyphkey)

attach_key(bowen_gm1819 , fcglyphkey) 

#--------------- 2019/20 wet season map

bowen_gm1920 <- build_gmap(data = fc_1920, geoData = bowen_utm, id = "SUBCATS", size = 1, 
                           glyph = "icone", palette = "Greens", 
                           limits = glyph_estlimits, max_error = glyph_errmax)
view(bowen_gm1920) + coord_equal() 

#build a glyph key 
fcglyphkey <- build_gkey(data = rbind(fc_1819, fc_1920) , glyph = "icone")
view(fcglyphkey)

attach_key(bowen_gm1920, fcglyphkey)

#------------------------- Exceedance Maps ---------------------------#

# define probability distribution (uniform distribution)
# Probabilities are P[X ≤ x]
pd <- quote({ punif(q, min, max, lower.tail = TRUE) })  

# redefine the parameters of a uniform from the mean and std
funa <- function(mu, sigma){2 * mu - funb(mu, sigma)} 
funb <- function(mu, sigma){0.5 * (sqrt(12*sigma^2) + 2*mu)}

# define argument listing
args <- quote({ list(min = funa(mu = estimate, sigma = error),
                     max = funb(mu = estimate, sigma = error)) })

# capture distribution and arguments in a single list
pdflist <- list(dist = pd, args = args, th = 50)

#--------------- 2018/19 wet season map

bowen_ExcMap1819 <- build_emap(data = fc_1819, pdflist = pdflist, 
                               geoData = bowen_shp, id = "SUBCATS", 
                               key_label = "Pr[X <= 50%]")
view(bowen_ExcMap1819) + coord_equal()

#--------------- 2019/20 wet season map

bowen_ExcMap1920 <- build_emap(data = fc_1920, pdflist = pdflist, 
                               geoData = bowen_shp, id = "SUBCATS", 
                               key_label = "Pr[X <= 50%]")
view(bowen_ExcMap1920) + coord_equal()





