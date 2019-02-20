# RMap.r ####
## R script for plotting point data in a map
## Adapted from valentinitnelav's scripts (https://gist.github.com/valentinitnelav)
## Source: https://seethedatablog.wordpress.com/2017/01/01/r-pacific-centered-world-map/
## Check the RMap Gudie.md document for more instructions.

# <START> ####

# <INPUT> ####
## DATA
input.data <- read.table("samples.txt", header = T, sep="\t", stringsAsFactors = F, comment.char = "")
long.col <- "Longitude"
lat.col <- "Latitude"
name.col <- "Population"

## Plot Modes
plotmode <- 4
### Mode 1: Binned Points (integer)
data.col <- "N"
bins <- c(10,20,50,100)

### Mode 2: Color Scale (numeric)
z.col <- "N"

### Mode 3: Colored Groups (categoric)
grp.col <- "Continent"

### Mode 4: Custom (add your code in the Mode 4 sections)
x.col <- ""

### Misc data

# MAP PARAMETERS ####
## Rounded corners, auto-off if zoomed in
rounded <- TRUE
## Shift map center. Atlantica: (0), Pacifica: (180).
shift <- 180+25
## Map limits (zoom)
### Recommended ratio 2:1, correct for shift value
long.lim <- c(-180,180)
lat.lim <- c(-90,90)

#OUTPUT FILE NAME
op.name <- "map_draft"

#POL
#long.lim <- c(-15,105)
#lat.lim <- c(-50,10)

#<CORE> ####
# LIBRARIES ####
library(data.table)
library(ggplot2)
library(ggrepel)
library(rgdal)
library(rgeos)
library(maps)
library(maptools)
library(raster)

# MAP FILES ####
#download.file(url = "http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_land.zip", destfile = "ne_110m_land.zip")

# unzip the shapefile in the directory mentioned with "exdir" argument
#unzip(zipfile="ne_110m_land.zip", exdir = "ne_110m_land")
# delete the zip file
#file.remove("ne_110m_land.zip")
# read the shapefile with readOGR from rgdal package
NE_continents <- readOGR(dsn = "ne_110m_land", layer = "ne_110m_land")

# MAP SHIFT ####
# shift central/prime meridian towards west - positive values only
WGS84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +towgs84=0,0,0")

split.line = SpatialLines(list(Lines(list(Line(cbind(180-shift,c(-90,90)))), ID="line")), 
                          proj4string=WGS84)

line.gInt <- gIntersection(split.line, NE_continents)

bf <- gBuffer(line.gInt, byid=TRUE, width=0.000001)  

NE_continents.split <- gDifference(NE_continents, bf, byid=TRUE)

# GRATICULES ####
if(!all((long.lim) == c(-180,180) & (lat.lim) == c(-90,90))){rounded <- FALSE}
# create a bounding box - world extent
b.box <- as(raster::extent(180, -180, -90, 90), "SpatialPolygons")
# assign CRS to box
proj4string(b.box) <- WGS84
# create graticules/grid lines from box
grid <- gridlines(b.box, 
                  easts  = seq(from=-180, to=180, by=20),
                  norths = seq(from=-90, to=90, by=10))

# create labels for graticules
grid.lbl <- labels(grid, side = 1:4)

# transform labels from SpatialPointsDataFrame to a data table that ggplot can use
grid.lbl.DT <- data.table(grid.lbl@coords, grid.lbl@data)

# prepare labels with regular expression:
# - delete unwanted labels
grid.lbl.DT[, labels := gsub(pattern="180\\*degree|90\\*degree\\*N|90\\*degree\\*S", replacement="", x=labels)]
# - replace pattern "*degree" with "°" (* needs to be escaped with \\)
grid.lbl.DT[, lbl := gsub(pattern="\\*degree", replacement="°", x=labels)]
# - delete any remaining "*"
grid.lbl.DT[, lbl := gsub(pattern="*\\*", replacement="", x=lbl)]

# adjust coordinates of labels so that they fit inside the globe
if (rounded){
  grid.lbl.DT[, long := ifelse(coords.x1 %in% c(-180,180), coords.x1*175/180, coords.x1)]
  grid.lbl.DT[, lat  := ifelse(coords.x2 %in% c(-90,90), coords.x2*82/90, coords.x2)]
} else {
  grid.lbl.DT[, long := coords.x1]
  grid.lbl.DT[, lat  := coords.x2]
}

# transform graticules from SpatialLines to a data table that ggplot can use
grid.DT <- data.table(map_data(SpatialLinesDataFrame(sl=grid, 
                                                     data=data.frame(1:length(grid)), 
                                                     match.ID = FALSE)))
# project coordinates

PROJ <- ifelse(rounded,"+proj=eck4","+proj=longlat")
PROJ <-paste(PROJ,"+lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# assign matrix of projected coordinates as two columns in data table
grid.DT[, c("X","Y") := data.table(project(cbind(long, lat), proj=PROJ))]

# project coordinates of labels
grid.lbl.DT[, c("X","Y") := data.table(project(cbind(long, lat), proj=PROJ))]

# CONTINENTS ####
# transform split country polygons in a data table that ggplot can use
Continent.DT <- data.table(map_data(as(NE_continents.split, "SpatialPolygonsDataFrame")))
# Shift coordinates
Continent.DT[, long.new := long + shift]
Continent.DT[, long.new := ifelse(long.new > 180, long.new-360, long.new)]
# project coordinates 
Continent.DT[, c("X","Y") := data.table(project(cbind(long.new, lat), proj=PROJ))]

# MODES ####
## Mode 1: Binned Points
if (plotmode == 1){
  NE_points.df     <- input.data[,c(long.col,lat.col,name.col,data.col)]
  names(NE_points.df) <- c("lon", "lat", "Name","N")
  # split population in desired intervals/classes
  bins.lbl <- c(paste("<",bins[1]),bins[1:(length(bins)-1)],paste(">",bins[length(bins)]))
}

## Mode 2: Color Scale
if (plotmode == 2){
  NE_points.df     <- input.data[,c(long.col,lat.col,name.col,z.col)]
  names(NE_points.df) <- c("lon", "lat", "Name","Z")
  NE_points.df$Z <- as.numeric(as.character(NE_points.df$Z)) * 1.0
}

## Mode 3: Colored Groups
if (plotmode == 3){
  NE_points.df     <- input.data[,c(long.col,lat.col,name.col,grp.col)]
  names(NE_points.df) <- c("lon", "lat", "Name","G")
}

## Mode 4: Custom Data
if (plotmode == 4){
  #[Add data-procesing code here]
}


## Convergence
NE_points.df$lon.new <- NE_points.df$lon + shift
NE_points.df$lon.new <- ifelse(NE_points.df$lon.new > 180, NE_points.df$lon.new-360, NE_points.df$lon.new)
prj.coord        <- project(cbind(NE_points.df$lon.new, NE_points.df$lat), proj = PROJ)
NE_points.df.prj <- cbind(prj.coord, NE_points.df)
names(NE_points.df.prj)[1:2] <- c("X.prj","Y.prj")
NE_points.dt.prj <- data.table(NE_points.df.prj)

## Mode 1: Binned Points
if (plotmode == 1){
  # split population in desired intervals/classes
  bins.lbl <- c(paste("<",bins[1]),bins[1:(length(bins)-1)],paste(">",bins[length(bins)]))
  # this will act as point size
  NE_points.dt.prj[, N.cls := cut(N,
                                  labels=bins.lbl,
                                  breaks=c(min(N)*0.999,floor(bins*0.999),max(N)),
                                  include.lowest=FALSE, right=T)]
  NE_points.dt.prj[, N.cls := as.integer(N.cls)]
  NE_points.dt.prj[, N.cls := bins.lbl[N.cls]]
}

# LIMITS ####
lx.lim <- project(cbind(long.lim,lat.lim), proj = PROJ)
long.lim <- t(lx.lim)[1,] * (ifelse(rounded,2,1))
lat.lim <- t(lx.lim)[2,]

# MAP PLOTTING ####
## Blank map
mapplot <- ggplot() + 
  geom_polygon(data = Continent.DT, 
               aes(x = X, y = Y, group = group), 
               fill = "gray70", 
               size = 0.25) +
  geom_path(data = grid.DT[(long %in% c(-180,180) & region == "NS")
                           |(long %in% c(-180,180) & lat %in% c(-90,90) & region == "EW")],
            aes(x = X, y = Y, group = group), 
            linetype = "solid", colour = "grey50", size = .3) +
  coord_fixed(ratio=ifelse(rounded,1,1.3), xlim= long.lim, ylim= lat.lim) + 
  theme_void()

## Mode 1: Binned Points
if (plotmode == 1){
  mapplot <- mapplot +
    geom_point(data = NE_points.dt.prj, 
               aes(x = X.prj, y = Y.prj, size = N, color = N.cls), 
               alpha = .8) +
    scale_color_discrete(limits= bins.lbl)+
    guides(color=guide_legend(title="Bins"), size = guide_legend(title = data.col))
}

## Mode 2: Color Scale
if (plotmode == 2){
  mapplot <- mapplot +
    geom_point(data = NE_points.dt.prj, 
               aes(x = X.prj, y = Y.prj, colour= Z), 
               size=2,  alpha = .8) +
    scale_colour_gradient(low = "blue", high = "red", guide = guide_colorbar(title = z.col))
}

## Mode 3: Colored Groups
if (plotmode == 3){
  mapplot <- mapplot +
    geom_point(data = NE_points.dt.prj, 
               aes(x = X.prj, y = Y.prj, colour= G), 
               size=2,  alpha = .8) +
    scale_color_discrete(name = grp.col)
}

## Mode 4: Custom Data
if (plotmode == 4){
  #[Add other ggplot2 functions for your data here]
}

## Extra Layers

## Labels and Legend
mapplot <- mapplot +
  geom_text_repel(data = NE_points.dt.prj, aes(x = X.prj, y = Y.prj, label = Name),
                  # label size
                  size = 2,
                  # color of the line segments
                  segment.colour = "#777777",
                  # line segment transparency
                  segment.alpha = .5,
                  # line segment thickness
                  segment.size = .25,
                  # the force of repulsion between overlapping text labels
                  force = 3,
                  # maximum number of iterations to attempt to resolve overlaps
                  max.iter = 10e3,
                  # turn labels off in the legend
                  show.legend = FALSE) +
  theme(legend.title = element_text(colour="black", size=10, face="bold"), # adjust legend title
        legend.position = c(1.02, 0.5), # relative position of legend
        plot.margin = unit(c(t=0, r=2, b=0, l=0), unit="cm")) # adjust margins

# <OUTPUT> ####
mapplot
ggsave(paste0(op.name,".pdf"), width=25, height=14, units="cm")
ggsave(paste0(op.name,".png"), width=25, height=14, units="cm", dpi=300)

# <END> ####
