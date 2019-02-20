# RMap

*By Pavel Salazar-Fernandez (epsalazarf@gmail.com)*

*Human Population and Evolutionary Genomics Lab | LANGEBIO*

## About
*Documentation: [See the data](https://seethedatablog.wordpress.com/2016/12/23/r-simple-world-map-robinson-ggplot/)*

This document explains how to use `RMap.R`, an R script for plotting different types of data on a world map.

## Procedure
The `RMap.R`script is an adaptation from various scripts developed by [valentinitnelav](https://github.com/valentinitnelav/) which plots a high resolution world map that can be longitude shifted, zoomed in, and is ggplot2 compatible.

### Preparations

**Software:**
- [R](https://cran.r-project.org/) (v3.2+)

**R packages**
- *Plotting*: [ggplot2](http://ggplot2.org/)
- *Mapping*: maps, maptools
- *Data management*: [data.table](https://github.com/Rdatatable/data.table/wiki)
- *Other*: ggrepel, raster, rgdal, rgeos

**Install all:**
`install.packages(c("data.table", "ggplot2", "ggrepel", "maps", "maptools", "raster", "rgdal", "rgeos"), dependencies= TRUE)`

### Map Data
By default, `RMap.R` uses the [1:110m Physical Land Vector](http://www.naturalearthdata.com/downloads/110m-physical-vectors/110m-land/) from [Natural Earth](http://www.naturalearthdata.com/), which includes only continental and island polygons with no country lines. Be sure that the `ne_110m_land/` directory is available at the running location.

### Input Data
The input data requires at least three columns: name, 
longitude and latitude. All columns must be named. Coordinates cannot exceed ±180° for longitude and ±90° for latitude.

Data columns can be numeric or categorical. Other columns, such as "Color", can be also be included. 

|*Name*|*Longitude*|*Latitude*|*Data*|...|
|---|---|---|---|---|
|ABC|0.0|0.0|100|...|
|DEF|-30.0|-30.0|200|...|
|GHI|120.0|60.0|300|...|
|...|...|...|...|...|

### Script Instructions
1. Go to the **INPUT** section of the script, where most of the customization is done. 
2. **DATA**: Fill out the file and columns names.
3. **MODE**: Select the mode number for your desired plot. The mode number may require different kind of data or additional information. Input for other modes will be ignored so you may skip them. The 4th mode the most open to customization if desired.
4. **MAP PARAMETERS**: . Activate `rounded` for a elliptical map, or deactivate for a rectangular map. `shift` re-centers the map to a given longitude. `long.lim`and `lat.lim` zooms in the map. A 2:1 ratio for the longitude and latitude zoom is recommended. Zooming in the map auto-deactivates the `rounded` option.
5. Change the output file name if desired.
6. If the map files are not available or want to use a different map, de-comment the downloading code in **MAP FILES**.
7. If you are customizing the 4th mode, head to the **CORE/MODES** section to add your custom code.
8. Run the script. A preview of the plotted map should show. If it doesn't, simply call the `mapplot` object. If a map is still not shown, check your data and code for any possible errors.
9. The maps generated are saved as PDF and PNG. Both versions are visually identical but may differ from the preview.
___
