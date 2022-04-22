#######################################################
### Seasonal rainfall in south korea during JAS
### 
###
#######################################################

# load the ncdf4 package
library(ncdf4)

# set path and filename
# at company 
setwd("c:/R/20180418_JAS rainfall forecasting/data")

# home 
#setwd("D:/R/20180418_JAS rainfall forecasting/data") 

# read rainfall data
df.rainfall <- read.csv("rainfall.csv")
str(df.rainfall)

# open a netCDF file
nc.sst <- nc_open("sst.mnmean.nc")
print(nc.sst)

# get longitude and latitude
lon <- ncvar_get(nc.sst,"lon")
lat <- ncvar_get(nc.sst,"lat")

nlon <- dim(lon)
nlat <- dim(lat)

print(c(nlon,nlat))

# get time
time <- ncvar_get(nc.sst, "time")
tunits <- ncatt_get(nc.sst,"time","units")
nt <- dim(time)

# get temperature
dname <- "sst"
mat.sst <- ncvar_get(nc.sst, dname)
fillvalue <- ncatt_get(nc.sst,dname,"_FillValue")
#dlname <- ncatt_get(nc.sst,dname,"long_name")
#dunits <- ncatt_get(nc.sst,dname,"units")
dim(mat.sst)

# convert time -- split the time units string into fields
library(chron)
library(lattice)
library(RColorBrewer)

tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth <- as.integer(unlist(tdstr)[2])
tday <- as.integer(unlist(tdstr)[3])
tyear <- as.integer(unlist(tdstr)[1])
obs.date <- as.Date(chron(time,origin=c(tmonth, tday, tyear)))

# reshape the netcdf into a dataframe
df.pos <- data.frame(expand.grid(lon,lat))
lonlat <- paste0("(E",df.pos[,1],",N",df.pos[,2], ")")

df.sst <- data.frame(as.data.frame(t(matrix(as.vector(mat.sst), nrow=nlon*nlat, ncol=nt))), date = obs.date)
names(df.sst) <- c(lonlat, "date")

# check validation of data reshaping
m <- 1
tmp_slice <- mat.sst[,,m]
tmp_slice[,6]

df.sst[1,901:1080]


# subset sst from 1980 to present
df.sst.subset <- subset(df.sst, date >= as.Date("1980-01-01"))

####################################
### draw map
###

start = 6
end = 6

sst <- data.frame()

for(year in df.rainfall$YEAR)
{
  ### set start and end 
  begda <- as.Date(paste0(year,"-",start,"-01"))
  endda <- as.Date(paste0(year,"-", end, "-20"))
  
  df.temp <- subset(df.sst.subset, date >= begda & date <= endda, select = -c(date))
  df.temp <- colMeans(df.temp)
  
  sst <- rbind(sst, df.temp)
  print(year)
}
names(sst) <- c(lonlat)

### calculate correlation of each cases

rt.cor.pearson <- c()
rt.cor.spearman <- c()
for(idx in (1:(nlon*nlat)))
{
  rt.cor.pearson <- append(rt.cor.pearson, cor(df.rainfall$JAS_RF, sst[,idx], method = "pearson"))
  rt.cor.spearman <- append(rt.cor.pearson, cor(df.rainfall$JAS_RF, sst[,idx], method = "spearman"))
}

mat.pearson <- as.matrix(rt.cor.pearson, nrow=nlon, ncol=nlat)
mat.spearman <- as.matrix(rt.cor.spearman, nrow=nlon, ncol=nlat)

head(sort(rt.cor.pearson))

head(sort(rt.cor.pearson, decreasing = T))


###################
# Boundary
###################

library(raster)
library(rasterVis)

library(maps)
library(mapdata)
library(maptools)

##
myRaster <- raster(mat.pearson)

## polygon shapefile
ext <- as.vector(extent(myRaster))

boundaries <- map('worldHires', fill=TRUE,
                  wrap = c(0,360),
                  #xlim=ext[1:2], ylim=ext[3:4],
                  plot=FALSE)

## read the map2SpatialPolygons help page for details
IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
bPols <- map2SpatialPolygons(boundaries, IDs=IDs,
                             proj4string=CRS(projection(myRaster)))

##############


library(viridis)
# levelplot of correlation
grid <- expand.grid(lon=lon, lat=lat)
cutpts <- seq(-1, 1, 0.05)
levelplot(mat.pearson ~ lon * lat, data=grid, at=cutpts, cuts=length(cutpts), pretty=T, 
          main="RAINFALL ~ SST(Mean of Jun), correlation(pearson)",
          #main="RAINFALL ~ SST(Mean of March,April,May), correlation(pearson)",
          col.regions = viridis_pal(option = "C")(length(cutpts))) + layer(sp.polygons(bPols))
          #col.regions=(rev(brewer.pal(length(cutpts),"RdBu"))))

