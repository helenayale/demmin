
#####################
### Starting Work ###
#####################

# set working directory
setwd("D:/R/demmin")
getwd()

# load packages
library(raster)
library(rgdal)  
library(RStoolbox)
#library(sp)
library(ggplot2)  
#library(ncdf4)
#library(SDMTools)
#library(plotKML)
#library(e1071)


########################
### loading data ###
########################

## loading & cropping function ##

loading <- function(date,shapefile){
  my_path <- getwd()
  my_path <- paste(my_path,"/",date, sep = "")
  file_list <- list.files(my_path, pattern='*.tif$')
  
  for (i in 1:length(file_list)){
    if (i == 1){
      # for the first run define the final raster file ...
      file_path <- paste(my_path,"/",file_list[i], sep = "")
      current_layer <- raster(file_path)
      in_stack <- stack(current_layer)
    } else {
      # ... and fill it with each additional run with another layer
      file_path <- paste(my_path,"/",file_list[i], sep = "")
      current_layer <- raster(file_path)
      in_stack <- stack(in_stack, current_layer)
    }
  }
  band_names <- c("pixel_qa", "radsat_qa", "aerosol_qa", "Aerosol", "Blue", "Green", "Red", "NIR", "SWIR.1", "SWIR.2")
  names(in_stack) <- band_names
  # load shapefile of the study area
  ww_felder <- readOGR(shapefile)
  ww_felder <- spTransform(ww_felder, CRS(proj4string(in_stack)))
  # crop the image to the extent of study area
  study_area <- crop(in_stack, ww_felder)
  return(study_area)
}



# load data of 2017

# put dates of data acquisition in 2017
date_string <- c("20170406","20170501","20170517","20170602","20170609","20170618","20170805","20170828","20170913")
sf_string <- "WW-felder.shp"

for (j in 1:length(date_string )){
  stacks <- loading(date_string[j],sf_string)
  eval(parse(text=paste(paste('stack',date_string[j],sep=''), '<- stacks')))
}

# check the images
# remove the ones with too many clouds (those were made into comments)
plotRGB(stack20170406,  r="Red", g="Green", b="Blue", stretch="lin")
plotRGB(stack20170501,  r="Red", g="Green", b="Blue", stretch="lin")
#plotRGB(stack20170517,  r="Red", g="Green", b="Blue", stretch="lin")
plotRGB(stack20170602,  r="Red", g="Green", b="Blue", stretch="lin")
plotRGB(stack20170609,  r="Red", g="Green", b="Blue", stretch="lin")
#plotRGB(stack20170618,  r="Red", g="Green", b="Blue", stretch="lin")
#plotRGB(stack20170805,  r="Red", g="Green", b="Blue", stretch="lin")
plotRGB(stack20170828,  r="Red", g="Green", b="Blue", stretch="lin")
#plotRGB(stack20170913,  r="Red", g="Green", b="Blue", stretch="lin")


# export the data as *.tif file
writeRaster(stack20170406, filename="stack20170406", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(stack20170501, filename="stack20170501", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(stack20170602, filename="stack20170602", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(stack20170609, filename="stack20170609", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(stack20170828, filename="stack20170828", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))


# load data of 2018

date_string <- c("20180409","20180418","20180504","20180520","20180605","20180628","20180707","20180723","20180909","20180925")
sf_string <- "WW-felder.shp"


for (j in 1:length(date_string )){
  stacks <- loading(date_string[j],sf_string)
  eval(parse(text=paste(paste('stack',date_string[j],sep=''), '<- stacks')))
}

# check the images
# remove the ones with too many clouds (those were made into comments)
plotRGB(stack20180409,  r="Red", g="Green", b="Blue", stretch="lin")
plotRGB(stack20180418,  r="Red", g="Green", b="Blue", stretch="lin")
plotRGB(stack20180504,  r="Red", g="Green", b="Blue", stretch="lin")
plotRGB(stack20180520,  r="Red", g="Green", b="Blue", stretch="lin")
#plotRGB(stack20180605,  r="Red", g="Green", b="Blue", stretch="lin")
plotRGB(stack20180628,  r="Red", g="Green", b="Blue", stretch="lin")
plotRGB(stack20180707,  r="Red", g="Green", b="Blue", stretch="lin")
plotRGB(stack20180723,  r="Red", g="Green", b="Blue", stretch="lin")
plotRGB(stack20180909,  r="Red", g="Green", b="Blue", stretch="lin")
#plotRGB(stack20180925,  r="Red", g="Green", b="Blue", stretch="lin")


# export the data as *.tif file
writeRaster(stack20180409, filename="stack20180409", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(stack20180418, filename="stack20180418", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(stack20180504, filename="stack20180504", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(stack20180520, filename="stack20180520", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(stack20180628, filename="stack20180628", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(stack20180707, filename="stack20180707", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(stack20180723, filename="stack20180723", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(stack20180909, filename="stack20180909", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))



###########################
### NDVI & LAI Analysis ###
###########################

## In this section, 
## indices are calculated


## NDVI function ##

## input: layerstack
## output: NDVI layer

NDVI <- function(stack){
  # get NIR & Red bands
  nir <- stack[["NIR"]]
  red <- stack[["Red"]]
  # calculate NDVI
  NDVI <- (nir-red) / (nir+red)
  return(NDVI)
}

## LAI function ##

## input: NDVI layer
## output: LAI layer

LAI <- function(ndvi){0.57*exp(2.33*ndvi)}  # Saito et al (2001)


# calculate NDVI for 2017
date_string <- c("20170406","20170501","20170602","20170609","20170828")

for (k in 1:length(date_string)){
  eval(parse(text=paste(paste('NDVI',date_string[k],sep='') ,'<-NDVI(',paste('stack',date_string[k],sep=''),')')))
}


# calculate NDVI for 2018
date_string <- c("20180409","20180418","20180504","20180520","20180628","20180707","20180723","20180909")

for (k in 1:length(date_string)){
  eval(parse(text=paste(paste('NDVI',date_string[k],sep='') ,'<-NDVI(',paste('stack',date_string[k],sep=''),')')))
}


# calculate LAI for 2017
date_string <- c("20170406","20170501","20170602","20170609","20170828")
for (k in 1:length(date_string)){
  eval(parse(text=paste(paste('LAI',date_string[k],sep='') ,'<-LAI(',paste('NDVI',date_string[k],sep=''),')')))
}


# calculate LAI for 2018
date_string <- c("20180409","20180418","20180504","20180520","20180628","20180707","20180723","20180909")

for (k in 1:length(date_string)){
  eval(parse(text=paste(paste('LAI',date_string[k],sep='') ,'<-LAI(',paste('NDVI',date_string[k],sep=''),')')))
}

# export the indices as *.tif file
writeRaster(NDVI20170406, filename="NDVI20170406", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(NDVI20170501, filename="NDVI20170501", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(NDVI20170602, filename="NDVI20170602", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(NDVI20170609, filename="NDVI20170609", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(NDVI20170828, filename="NDVI20170828", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))


writeRaster(NDVI20180409, filename="NDVI20180409", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(NDVI20180418, filename="NDVI20180418", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(NDVI20180504, filename="NDVI20180504", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(NDVI20180520, filename="NDVI20180520", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(NDVI20180628, filename="NDVI20180628", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(NDVI20180707, filename="NDVI20180707", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(NDVI20180723, filename="NDVI20180723", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(NDVI20180909, filename="NDVI20180909", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))

writeRaster(LAI20170406, filename="LAI20170406", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(LAI20170501, filename="LAI20170501", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(LAI20170602, filename="LAI20170602", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(LAI20170609, filename="LAI20170609", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(LAI20170828, filename="LAI20170828", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))


writeRaster(LAI20180409, filename="LAI20180409", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(LAI20180418, filename="LAI20180418", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(LAI20180504, filename="LAI20180504", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(LAI20180520, filename="LAI20180520", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(LAI20180628, filename="LAI20180628", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(LAI20180707, filename="LAI20180707", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(LAI20180723, filename="LAI20180723", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(LAI20180909, filename="LAI20180909", format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))




########################
### Temperature Data ###
########################


# temperature date are from ERA Interim, Daily data
# downloaded from ECMWF database
# raw data are in format *.nc

## temperature function

## input: max - raw data of max temperature
##        min - raw data of min temperature
##        ww_felder - study area as shapefile


temp_mean <- function(max,min,ww_felder){
  # crop the data to the extent of study area
  gg<-spTransform(ww_felder,CRS(proj4string(max)))
  crop_max <- crop(max,gg)
  crop_min <- crop(min,gg)
  # minus 273.15 because the temperature is Kalvin Temperature
  crop_max_c <- crop_max - 273.15
  crop_min_c <- crop_min - 273.15
  # calculate the mean temperature
  mean_temp <-  (crop_max_c + crop_min_c)/2
  return(mean_temp)
}

temp_pro <-  function(asa,ww_felder){
  # crop the data to the extent of study area
  gg<-spTransform(ww_felder,CRS(proj4string(asa)))
  crop_asa <- crop(asa,gg)
  # minus 273.15 because the temperature is Kalvin Temperature
  crop_asa_c <- crop_asa - 273.15
  return(crop_asa_c )
}



## prepare the temperature data for each month
# import shapefile of study field
ww_felder <- readOGR("WW-felder.shp")
# 2017 April 06
max201704 <- brick("netcdf/201704max1300.nc")
min201704 <- brick("netcdf/201704min1300.nc")

# get temperature of the month
temp_max201704 <- temp_pro(max201704,ww_felder)
temp_min201704 <- temp_pro(min201704,ww_felder)
temp_mean201704 <- temp_mean(max201704,min201704,ww_felder)
# find the date needed manually
temp_max20170406 <- temp_max201704$X2017.04.06.13.00.00
temp_min20170406 <- temp_min201704$X2017.04.06.13.00.00
temp_mean20170406 <- temp_mean201704$X2017.04.06.13.00.00


# 2017 May 01
max201705 <- brick("netcdf/201705max1300.nc")
min201705 <- brick("netcdf/201705min1300.nc")
# get temperature of the month
temp_max201705 <- temp_pro(max201705,ww_felder)
temp_min201705 <- temp_pro(min201705,ww_felder)
temp_mean201705 <- temp_mean(max201705,min201705,ww_felder)

# find the date needed manually
temp_max20170501 <- temp_max201705$X2017.05.01.13.00.00
temp_min20170501 <- temp_min201705$X2017.05.01.13.00.00
temp_mean20170501 <- temp_mean201705$X2017.05.01.13.00.00


# 2017 June 02
max201706 <- brick("netcdf/201706max1300.nc")
min201706 <- brick("netcdf/201706min1300.nc")

# get temperature of the month
temp_max201706 <- temp_pro(max201706,ww_felder)
temp_min201706 <- temp_pro(min201706,ww_felder)
temp_mean201706 <- temp_mean(max201706,min201706,ww_felder)

# find the date needed manually
temp_max20170602 <- temp_max201706$X2017.06.02.13.00.00
temp_min20170602 <- temp_min201706$X2017.06.02.13.00.00
temp_mean20170602 <- temp_mean201706$X2017.06.02.13.00.00


# 2017 June 09
temp_max20170609 <- temp_max201706$X2017.06.09.13.00.00
temp_min20170609 <- temp_min201706$X2017.06.09.13.00.00
temp_mean20170609 <- temp_mean201706$X2017.06.09.13.00.00

# 2017 July
max201707 <- brick("netcdf/201707max1300.nc")
#min201707 <- brick("netcdf/201707min1300.nc")

# get temperature of the month
temp_max201707 <- temp_pro(max201707,ww_felder)



# 2017 Aug. 28
max201708 <- brick("netcdf/201708max1300.nc")
min201708 <- brick("netcdf/201708min1300.nc")

# get temperature of the month
temp_max201708 <- temp_pro(max201708,ww_felder)
temp_min201708 <- temp_pro(min201708,ww_felder)
temp_mean201708 <- temp_mean(max201708,min201708,ww_felder)

# find the date needed manually
temp_max20170828 <- temp_max201708$X2017.08.28.13.00.00
temp_min20170828 <- temp_min201708$X2017.08.28.13.00.00
temp_mean20170828 <- temp_mean201708$X2017.08.28.13.00.00

# 2018 April 09
max201804 <- brick("netcdf/201804max1300.nc")
min201804 <- brick("netcdf/201804min1300.nc")

# get temperature of the month
temp_max201804 <- temp_pro(max201804,ww_felder)
temp_min201804 <- temp_pro(min201804,ww_felder)
temp_mean201804 <- temp_mean(max201804,min201804,ww_felder)

# find the date needed manually
temp_max20180409 <- temp_max201804$X2018.04.09.13.00.00
temp_min20180409 <- temp_min201804$X2018.04.09.13.00.00
temp_mean20180409 <- temp_mean201804$X2018.04.09.13.00.00


# 2018 April 18
temp_max20180418 <- temp_max201804$X2018.04.18.13.00.00
temp_min20180418 <- temp_min201804$X2018.04.18.13.00.00
temp_mean20180418 <- temp_mean201804$X2018.04.18.13.00.00


# 2018 May 04
max201805 <- brick("netcdf/201805max1300.nc")
min201805 <- brick("netcdf/201805min1300.nc")

# get temperature of the month
temp_max201805 <- temp_pro(max201805,ww_felder)
temp_min201805 <- temp_pro(min201805,ww_felder)
temp_mean201805 <- temp_mean(max201805,min201805,ww_felder)

# find the date needed manually
temp_max20180504 <- temp_max201805$X2018.05.04.13.00.00
temp_min20180504 <- temp_min201805$X2018.05.04.13.00.00
temp_mean20180504 <- temp_mean201805$X2018.05.04.13.00.00

# 2018 May 20
temp_max20180520 <- temp_max201805$X2018.05.20.13.00.00
temp_min20180520 <- temp_min201805$X2018.05.20.13.00.00
temp_mean20180520 <- temp_mean201805$X2018.05.20.13.00.00

#2018 June 28
max201806 <- brick("netcdf/201806max1300.nc")
min201806 <- brick("netcdf/201806min1300.nc")

# get temperature of the month
temp_max201806 <- temp_pro(max201806,ww_felder)
temp_min201806 <- temp_pro(min201806,ww_felder)
temp_mean201806 <- temp_mean(max201806,min201806,ww_felder)

# find the date needed manually
temp_max20180628 <- temp_max201806$X2018.06.28.13.00.00
temp_min20180628 <- temp_min201806$X2018.06.28.13.00.00
temp_mean20180628 <- temp_mean201806$X2018.06.28.13.00.00

# 2018 July 07
max201807 <- brick("netcdf/201807max1300.nc")
min201807 <- brick("netcdf/201807min1300.nc")


# get temperature of the month
temp_max201807 <- temp_pro(max201807,ww_felder)
temp_min201807 <- temp_pro(min201807,ww_felder)
temp_mean201807 <- temp_mean(max201807,min201807,ww_felder)

# find the date needed manually
temp_max20180707 <- temp_max201807$X2018.07.07.13.00.00
temp_min20180707 <- temp_min201807$X2018.07.07.13.00.00
temp_mean20180707 <- temp_mean201807$X2018.07.07.13.00.00

# 2018 July 23
temp_max20180723 <- temp_max201807$X2018.07.23.13.00.00
temp_min20180723 <- temp_min201807$X2018.07.23.13.00.00
temp_mean20180723 <- temp_mean201807$X2018.07.23.13.00.00

# 2018 Aug.
max201808 <- brick("netcdf/201808max1300.nc")
#min2018098 <- brick("netcdf/201808min1300.nc")

# get temperature of the month
temp_max201808 <- temp_pro(max201808,ww_felder)

# 2018 Sep. 09
max201809 <- brick("netcdf/201809max1300.nc")
min201809 <- brick("netcdf/201809min1300.nc")

# get temperature of the month
temp_max201809 <- temp_pro(max201809,ww_felder)
temp_min201809 <- temp_pro(min201809,ww_felder)
temp_mean201809 <- temp_mean(max201809,min201809,ww_felder)

# find the date needed manually
temp_max20180909 <- temp_max201809$X2018.09.09.13.00.00
temp_min20180909 <- temp_min201809$X2018.09.09.13.00.00
temp_mean20180909 <- temp_mean201809$X2018.09.09.13.00.00


#############################
### Combine and Plot Data ###
#############################


# reload the layerstacks

## reloading function
reloading <- function(stack){
  in_stack <- brick(stack)
  band_names <- c("pixel_qa", "radsat_qa", "aerosol_qa", "Aerosol", "Blue", "Green", "Red", "NIR", "SWIR.1", "SWIR.2")
  names(in_stack) <- band_names
  return(in_stack)
}


# 2017
stack20170406 <- reloading("stack20170406.tif")
stack20170501 <- reloading("stack20170501.tif")
stack20170602 <- reloading("stack20170602.tif")
stack20170609 <- reloading("stack20170609.tif")
stack20170828 <- reloading("stack20170828.tif")
# 2018
stack20180409 <- reloading("stack20180409.tif")
stack20180418 <- reloading("stack20180418.tif")
stack20180504 <- reloading("stack20180504.tif")
stack20180520 <- reloading("stack20180520.tif")
stack20180628 <- reloading("stack20180628.tif")
stack20180707 <- reloading("stack20180707.tif")
stack20180723 <- reloading("stack20180723.tif")
stack20180909 <- reloading("stack20180909.tif")
stack20180925 <- reloading("stack20180925.tif")


# plot sample points on top of the images
# import sample points
ww <-readOGR("WW.shp")
ww2<-spTransform(ww,CRS(proj4string(stack20170406)))

# plot 2017
plotRGB(stack20170406,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)
plotRGB(stack20170501,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)
plotRGB(stack20170602,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)
plotRGB(stack20170609,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)
plotRGB(stack20170828,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)

# plot 2018
plotRGB(stack20180409,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)
plotRGB(stack20180418,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)
plotRGB(stack20180504,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)
plotRGB(stack20180520,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)
plotRGB(stack20180628,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)
plotRGB(stack20180707,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)
plotRGB(stack20180723,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)
plotRGB(stack20180909,  r="Red", g="Green", b="Blue", stretch="lin")
plot(ww2,col="red",add=TRUE)


# pick out the points that are always outside clouds
# import another shapefile including only these points
ww3 <- readOGR("WW3.shp")

# transfer coordination system to the same as stacks
ww3<-spTransform(ww3,CRS(proj4string(stack20180925)))
plot(ww3,col="yellow",add=TRUE)

# get the coordinates of the sample point (in the new coordination system)
cor_ww3 <- coordinates(ww3)
cor_ww3

# transfer sample points data to the same as temperature data
gg2<-spTransform(ww3,CRS(proj4string(temp_max20170406)))
# get the coordinates of the sample point (in the new coordination system)
cor_gg2 <- coordinates(gg2)
cor_gg2 <- cor_gg2[,-3]
cor_gg2

# define a dataframe to save the indices for each year
Data_ww2017 <- data.frame("year" = 0, "date_real" = 0, "date" = 0, "Pt" =0, "NDVI" = 0,
                          "LAI" = 0, "Temp_max" = 0, "Temp_min" = 0, "Temp_mean" = 0)
Data_ww2018 <- data.frame("year" = 0, "date_real" = 0, "date" = 0,"Pt" =0, "NDVI" = 0, 
                          "LAI" = 0, "Temp_max" = 0, "Temp_min" = 0, "Temp_mean" = 0)


# extract NDVI, LAI and temperature values at the sample points
# extract 2017
date_string <- c("20170406","20170501","20170602","20170609","20170828")

for (i in 1:length(date_string)){
  # extract NDVI & LAI
  a <- vector()
  b <- vector()
  eval(parse(text=paste('a <- extract(',paste('NDVI',date_string[i],sep=''),',cor_ww3)')))
  eval(parse(text=paste('b <- extract(',paste('LAI',date_string[i],sep=''),',cor_ww3)')))
  # extract temperatures
  c <- vector()
  d <- vector()
  e <- vector()
  eval(parse(text=paste('c <- extract(',paste('temp_max',date_string[i],sep=''),',cor_gg2)')))
  eval(parse(text=paste('d <- extract(',paste('temp_min',date_string[i],sep=''),',cor_gg2)')))
  eval(parse(text=paste('e <- extract(',paste('temp_mean',date_string[i],sep=''),',cor_gg2)')))
  
  for(j in 1:length(a)){
    Data_ww2017[(i-1)*6+j,"year"] <- 2017
    Data_ww2017[(i-1)*6+j,"date_real"] <- date_string[i]
    Data_ww2017[(i-1)*6+j,"date"] <- i
    Data_ww2017[(i-1)*6+j,"Pt"] <- j
    Data_ww2017[(i-1)*6+j,"NDVI"] <- a[j]
    Data_ww2017[(i-1)*6+j,"LAI"] <- b[j]
    Data_ww2017[(i-1)*6+j,"Temp_max"] <- c[j]
    Data_ww2017[(i-1)*6+j,"Temp_min"] <- d[j]
    Data_ww2017[(i-1)*6+j,"Temp_mean"] <- e[j]
  }
}
Data_ww2017

write.csv(Data_ww2017, file="Data_ww2017.csv")

# NDVI & LAI plots for separate points
ggplot(Data_ww2017,aes(x=date, y=NDVI, color=Pt))+geom_point()+
  facet_wrap(~Pt)+stat_smooth()
ggplot(Data_ww2017,aes(x=date, y=LAI, color=Pt))+geom_point()+
  facet_wrap(~Pt)+stat_smooth()

# NDVI & LAI plots binded
ggplot(Data_ww2017,aes(x=date, y=NDVI, color=Pt))+geom_point()+stat_smooth()
ggplot(Data_ww2017,aes(x=date, y=LAI, color=Pt))+geom_point()+stat_smooth()

# temperature plots
ggplot(Data_ww2017,aes(x=date, y=Temp_max, color=Pt))+geom_point()+stat_smooth()
ggplot(Data_ww2017,aes(x=date, y=Temp_min, color=Pt))+geom_point()+stat_smooth()
ggplot(Data_ww2017,aes(x=date, y=Temp_mean, color=Pt))+geom_point()+stat_smooth()


# extract 2018
date_string <- c("20180409","20180504","20180520","20180628","20180909")
#date_string <- c("20180409","20180418","20180504","20180520","20180628","20180707","20180723","20180909")

for (i in 1:length(date_string)){
  # extract NDVI & LAI
  a <- vector()
  b <- vector()
  eval(parse(text=paste('a <- extract(',paste('NDVI',date_string[i],sep=''),',cor_ww3)')))
  eval(parse(text=paste('b <- extract(',paste('LAI',date_string[i],sep=''),',cor_ww3)')))
  # extract temperatures
  c <- vector()
  d <- vector()
  e <- vector()
  eval(parse(text=paste('c <- extract(',paste('temp_max',date_string[i],sep=''),',cor_gg2)')))
  eval(parse(text=paste('d <- extract(',paste('temp_min',date_string[i],sep=''),',cor_gg2)')))
  eval(parse(text=paste('e <- extract(',paste('temp_mean',date_string[i],sep=''),',cor_gg2)')))
  
  for(j in 1:length(a)){
    Data_ww2018[(i-1)*6+j,"year"] <- 2018
    Data_ww2018[(i-1)*6+j,"date_real"] <- date_string[i]
    Data_ww2018[(i-1)*6+j,"date"] <- i
    Data_ww2018[(i-1)*6+j,"Pt"] <- j
    Data_ww2018[(i-1)*6+j,"NDVI"] <- a[j]
    Data_ww2018[(i-1)*6+j,"LAI"] <- b[j]
    Data_ww2018[(i-1)*6+j,"Temp_max"] <- c[j]
    Data_ww2018[(i-1)*6+j,"Temp_min"] <- d[j]
    Data_ww2018[(i-1)*6+j,"Temp_mean"] <- e[j]
  }
}
Data_ww2018

write.csv(Data_ww2018, file="Data_ww2018.csv")

# NDVI & LAI plots for separate points
ggplot(Data_ww2018,aes(x=date, y=NDVI, color=Pt))+geom_point()+
  facet_wrap(~Pt)+stat_smooth()
ggplot(Data_ww2018,aes(x=date, y=LAI, color=Pt))+geom_point()+
  facet_wrap(~Pt)+stat_smooth()

# NDVI & LAI plots binded
ggplot(Data_ww2018,aes(x=date, y=NDVI, color=Pt))+geom_point()+stat_smooth()
ggplot(Data_ww2018,aes(x=date, y=LAI, color=Pt))+geom_point()+stat_smooth()


# temperature plots
ggplot(Data_ww2018,aes(x=date, y=Temp_max, color=Pt))+geom_point()+stat_smooth()
ggplot(Data_ww2018,aes(x=date, y=Temp_min, color=Pt))+geom_point()+stat_smooth()
ggplot(Data_ww2018,aes(x=date, y=Temp_mean, color=Pt))+geom_point()+stat_smooth()



Data_Temp2017<- data.frame("year" = 0, "month" = 0, "date" = 0,  "Temp" = 0)
Data_Temp2018<- data.frame("year" = 0, "month" = 0, "date" = 0,  "Temp" = 0)

# import single point shapefile
ww_single <- readOGR("WW_single.shp")
# transfer sample points data to the same as temperature data
ss<-spTransform(ww_single ,CRS(proj4string(temp_max20170406)))
# get the coordinates of the sample point (in the new coordination system)
cor_ss <- coordinates(ss)
cor_ss 


Temp_ex201704 <- extract(temp_max201704, cor_ss )
Temp_ex201705 <- extract(temp_max201705, cor_ss )
Temp_ex201706 <- extract(temp_max201706, cor_ss )
Temp_ex201707 <- extract(temp_max201707, cor_ss )
Temp_ex201708 <- extract(temp_max201708, cor_ss )
Temp_ex201804 <- extract(temp_max201804, cor_ss )
Temp_ex201805 <- extract(temp_max201805, cor_ss )
Temp_ex201806 <- extract(temp_max201806, cor_ss )
Temp_ex201807 <- extract(temp_max201807, cor_ss )
Temp_ex201808 <- extract(temp_max201808, cor_ss )
Temp_ex201809 <- extract(temp_max201809, cor_ss )



date_string <- c("201704","201705","201706","201707","201708")


for (i in 1:length(date_string)){
  # extract temperatures
  t <- vector()

  eval(parse(text=paste('t <- Temp_ex',date_string[i],sep='')))
  for(j in 1:length(t)){
    Data_Temp2017[(i-1)*31+j,"year"] <- 2017
    Data_Temp2017[(i-1)*31+j,"month"]  <- i+3
    Data_Temp2017[(i-1)*31+j,"date"]  <-  j
    Data_Temp2017[(i-1)*31+j,"Temp"]  <-  t[j]
  }
}
Data_Temp2017

Data_Temp2017 <- na.omit(Data_Temp2017, invert=FALSE)

write.csv(Data_Temp2017, file="Data_Temp2017.csv")

ggplot(Data_Temp2017,aes(x=date, y=Temp, col=month))+geom_point()+facet_wrap(~month)+stat_smooth()


date_string <- c("201804","201805","201806","201807","201808")


for (i in 1:length(date_string)){
  # extract temperatures
  t <- vector()
  
  eval(parse(text=paste('t <- Temp_ex',date_string[i],sep='')))
  for(j in 1:length(t)){
    Data_Temp2018[(i-1)*31+j,"year"] <- 2017
    Data_Temp2018[(i-1)*31+j,"month"]  <- i+3
    Data_Temp2018[(i-1)*31+j,"date"]  <-  j
    Data_Temp2018[(i-1)*31+j,"Temp"]  <-  t[j]
  }
}
Data_Temp2018

Data_Temp2018 <- na.omit(Data_Temp2018, invert=FALSE)
write.csv(Data_Temp2018, file="Data_Temp2018.csv")

ggplot(Data_Temp2018,aes(x=date, y=Temp, col=month))+geom_point()+facet_wrap(~month)+stat_smooth()

