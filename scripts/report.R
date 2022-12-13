#### R ENVIRONMENT ####

library(ggh4x)
library(flexpolyline)
library(cluster)
library(raster)
library(rgdal)
library(soiltexture)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(lubridate)
library(tabulizer)
library(readr)
library(sf)
library(ggpubr)
library(RStoolbox)
library(ggshadow)
library(shadowtext)
library(ggspatial)
library(ggnewscale)
library(paletteer)
library(caret)
library(viridis)
library(ggimage)
library(cowplot)
library(xlsx)
library(forcats)
library(ggtern)
library(dplyr)

#### INPUT DATA FOLDER ####

#data folder
data.folder <- "E:/DATA/pambamarca/data/PLOSONE/data"

#### REQUIERED FUNCTIONS ####

#summary function shapefile
summaryShapefile <- function(InputShape=suelos.shp,columName="orden",ReferenceShape=study.area){
  columName <- paste0("^",columName,"$")
  InputShape <- spTransform(InputShape,CRS(ReferenceShape@proj4string@projargs))
  InputShape$AREA_HA <- area(InputShape) *0.0001
  InputShape <- InputShape@data
  targerCol <- grep(columName,names(InputShape))
  InputShape <- aggregate(InputShape$AREA_HA,list(InputShape[,targerCol]),FUN=sum)
  names(InputShape) <- c("Class","ha")
  InputShape$perc <- (InputShape$ha * 100) / sum(InputShape$ha)
  InputShape <- InputShape[order(InputShape$ha,decreasing = T),]
  return(InputShape)
}
#summary function raster
summaryRaster <- function(InputRaster=koeppen.ras,ReferenceShape=study.area,categorical=T){
  InputRaster <- crop(InputRaster,study.area)
  if(categorical){
    if(length(InputRaster@data@attributes)!=0){
      catVal <- as.data.frame(InputRaster@data@attributes)
    }else{
      catVal <- unique(values(InputRaster)) 
    }
    areasVal <- as.data.frame(table(values(InputRaster)))
    areasVal$perc <- (areasVal$Freq * 100) / sum(areasVal$Freq)
    areasVal$cat <- NA
    if(length(InputRaster@data@attributes)!=0){
      for(i in 1:nrow(areasVal)){
        check.val <- as.numeric(as.character(areasVal[i,"Var1"]))
        check.val <- catVal$ID %in% check.val
        areasVal[i,"cat"] <- catVal[check.val,2]
      }
    }
    areasVal <- areasVal[order(areasVal$perc,decreasing = T),]
    
    return(areasVal)
  }else{
    numVal <- values(InputRaster)
    quanVal <- summary(numVal)
    numDF <- as.data.frame(t(as.matrix(as.numeric(quanVal))))
    names(numDF) <- attr(quanVal,"names")
    numDF <- numDF[,c("Min.","1st Qu.","Median","3rd Qu.","Max.")]
    numDF <- cbind(data.frame(meanVal=round(mean(numVal,na.rm=T),2),SD=round(sd(numVal,na.rm=T),2)),numDF)
    return(numDF)
  }
}
#cross shapes
crossShapefiles <- function(x=geo.shp,y=suelos.shp,cols=c("orden","MESORELIEV"),ReferenceShape=study.area){
  x <- spTransform(x,CRS(ReferenceShape@proj4string@projargs))
  y <- spTransform(y,CRS(ReferenceShape@proj4string@projargs))
  xy <- raster::union(x,y)
  xy$AREA_HA <- round(area(xy) *0.0001,2)
  xy <- xy@data[,c(cols,"AREA_HA")]
  xy <- aggregate(xy$AREA_HA,list(xy[,cols[1]],xy[,cols[2]]),FUN=sum)
  xy$perc <- (xy[,3] * 100) / sum(xy[,3])
  xy <- xy[order(xy$perc,decreasing = T),]
  names(xy)[3] <- "AREA_HA"
  return(xy)
}
#cross shapes & raster
crossShpRas <- function(x=suelos.shp,y=altitude,col.num=3,ReferenceShape=study.area,func=mean){
  x <- spTransform(x,CRS(ReferenceShape@proj4string@projargs))
  y <- raster::extract(y,x,fun=mean)
  x.data <- data.frame(x=x@data[,col.num],y=as.numeric(y))
  x.data <- stats::aggregate(x.data$y,by=list(x.data$x),FUN=func,na.rm = TRUE)
  return(x.data)
}

#### 2.1 The Pambamarca Complex in northern Ecuador ####

#get publication PDF
input_pdf <- paste0(data.folder,"/references/Connel_Pambamarca.pdf")
#read & extract table
tbl <- extract_tables(input_pdf,pages=6)
tbl <- as.data.frame(tbl[[1]])
names(tbl) <- as.character(apply(tbl[1:2,],2,function(x){paste(x,collapse=" ")}))
tbl <- tbl[3:nrow(tbl),]
tbl <- tbl[-6,]
rownames(tbl) <- NULL
tbl[5,1] <- "Pambamarca/Frances Urco" #fix name
#select columns and prepare base
tbl[,4] <- parse_number(tbl[,4])
tbl[,5] <- parse_number(tbl[,5])
tbl[,7] <- parse_number(tbl[,7])
tbl <- tbl[,c(1,6,4,5)]
names(tbl) <- c("names","code","altitude","interior") 
#recalculate areas as HA
tbl$interior <- round(tbl$interior*0.0001,1)
#derive evident and buffered areas from pucaras
evident.shp <- st_read(paste0(data.folder,"/studyArea/pucaras_pol.shp"))
buffered.shp <- st_read(paste0(data.folder,"/studyArea/pucaras_pol_buffer.shp"))
evident.shp$evident <- st_area(evident.shp)
evident.shp <- st_drop_geometry(evident.shp)
buffered.shp$buffered <- st_area(buffered.shp)
buffered.shp <- st_drop_geometry(buffered.shp)
areas.df <- merge(evident.shp,buffered.shp,by="code")
areas.df <- areas.df[,c(1,2,4,6)]
areas.df$evident <- round(as.numeric(areas.df$evident) * 0.0001,1)
areas.df$buffered <- round(as.numeric(areas.df$buffered) * 0.0001,1)
#merge df and print
tbl <- merge(tbl,areas.df,by="code",all=T)
tbl$name.x <- NULL
tbl$code[1] <- NA
tbl.order <- c("Pi10","Pi14","Pi17","Pi18","Pi19","Pi20","Pi22","Pi23","Pi11","Pi12","Pi13","Pi15","Pi16","Pi21","Pi24","Pi25","Pi26",NA)
tbl[match(tbl.order,tbl$code),]

#### 2.2 Climate and landscape overview ####

#study area
study.area <- shapefile(paste0(data.folder,"/studyArea/pambamarca_area.shp"))
#Koppen-Geiger present climate (please see Beck et al. 2018)
koppenGeiger.classes <- data.frame(val=c(9,10,12,15,16,29),
                                   class=c("Csb","Csc","Cwb","Cfb","Cfc","ET"))
koppenGeiger.present <- raster(paste0(data.folder,"/characterization/KoppenGeiger_Present.tif"))
koppenGeiger.present <- summaryRaster(koppenGeiger.present,categorical=T)
koppenGeiger.present$cat <- sapply(koppenGeiger.present$Var1,function(x){
  x <- koppenGeiger.classes[koppenGeiger.classes$val  %in% x,2]
})
koppenGeiger.present.confidence <- raster(paste0(data.folder,"/characterization/KoppenGeiger_Present_confidence.tif"))
mean(values(koppenGeiger.present.confidence))
sd(values(koppenGeiger.present.confidence))
#Koppen-Geiger future climate 
koppenGeiger.future <- raster(paste0(data.folder,"/characterization/KoppenGeiger_Future.tif"))
koppenGeiger.future <- summaryRaster(koppenGeiger.future,categorical=T)
koppenGeiger.future$cat <- sapply(koppenGeiger.future$Var1,function(x){
  x <- koppenGeiger.classes[koppenGeiger.classes$val  %in% x,2]
})
koppenGeiger.future.confidence <- raster(paste0(data.folder,"/characterization/KoppenGeiger_Future_confidence.tif"))
mean(values(koppenGeiger.future.confidence))
sd(values(koppenGeiger.future.confidence))
#measure change
koppenGeiger.present$future <- NA
for(i in 1:nrow(koppenGeiger.present)){
  i.present <- koppenGeiger.present[i,] 
  i.future <- koppenGeiger.future[koppenGeiger.future$cat %in% i.present$cat,]
  if(nrow(i.future)==0){
    koppenGeiger.present$future[i] <- 0
  }else{
    koppenGeiger.present$future[i] <- i.future$perc
  }
}
names(koppenGeiger.present) <- c("Value","Freq","PresentPerc","class","FuturePerc")
koppenGeiger.present <- koppenGeiger.present[,c(1,2,4,3,5)]
koppenGeiger.present$diff <- koppenGeiger.present$FuturePerc - koppenGeiger.present$PresentPerc
koppenGeiger.present
#landforms
geo.shp <- shapefile(paste0(data.folder,"/characterization/geoformas_MAE.shp"))
summaryShapefile(geo.shp,"RELIEVE_GE",study.area)
summaryShapefile(geo.shp,"MACRORELIE",study.area)
summaryShapefile(geo.shp,"MESORELIEV",study.area)
#land cover
magap.shp <- shapefile(paste0(data.folder,"/characterization/cobertura2014_MAGAP.shp"))
summaryShapefile(magap.shp,"NIVEL1",study.area)
sum.all <- summaryShapefile(magap.shp,"NIVEL2",study.area)
sum.intervened <- sum.all[sum.all$Class %in% c("PASTIZAL","MOSAICO AGROPECUARIO","CULTIVO ANUAL","PLANTACION FORESTAL","INFRAESTRUCTURA","AREA POBLADA","CULTIVO PERMANENTE","OTRAS TIERRAS AGRICOLAS","AREA SIN COBERTURA VEGETAL"),]
sum.intervened
sum(sum.intervened$perc) #intevened %
100-sum(sum.intervened$perc) #non-intervened %
13.9+7.3+0.22+0.06 #croplands %
#ecosystems
mae.shp <- shapefile(paste0(data.folder,"/characterization/ecosistemas_MAE.shp"))
summaryShapefile(mae.shp,"ECOSISTEMA",study.area)
sum.natural <- sum.all[!sum.all$Class %in% c("PASTIZAL","MOSAICO AGROPECUARIO","CULTIVO ANUAL","PLANTACION FORESTAL","INFRAESTRUCTURA","AREA POBLADA","CULTIVO PERMANENTE","OTRAS TIERRAS AGRICOLAS","AREA SIN COBERTURA VEGETAL"),]
sum.natural #ecosystems percentages adjusted to land cover non-intervened areas report 
21.4+2.8 #approximate Northern Semi-deciduous Forest and Shrubland of the Valleys

#### 3.1 Rainfall factor (R-factor) for historic and climate change scenarios  ####

#clasification error of climate models (See: Cannon AJ. Reductions in daily continental-scale atmospheric circulation biases between generations of global climate models: CMIP5 to CMIP6. Environ Res Lett. 2020 Jun 1;15(6):064006. )
#IPSL-CM6A-LR 0.2-0.4, 0.6-0.8, 0.6-0.8, 0.2-0.4
mean(c(0.3,0.7,0.7,0.3))
#MIROC6 0.4-0.6, 0.8-1, 0.4-0.6, 0.2-0.4
mean(c(0.5,0.9,0.5,0.3))
#CNRM-CM6-1 0.8-1, 1.0-1.2, 0.6-0.8, 0.6-0.8
mean(c(0.9,1.1,0.7,0.7))
#GFDL-ESM4 0.8-1, 1.2-1.4, 0.8-1, 0.6-0.8
mean(c(0.9,1.3,0.9,0.7))
#historical rainfall (all models have the same historical dataset)
rainfall <- raster(paste0(data.folder,"/RUSLE/CNRM-CM6-1/r_rainfall/rainfall/historical.tif"))
dem <- raster(paste0(data.folder,"/DEM/pambamarca_100mts_dem_crop.tif"))
#rainfall for altitudes higher >3k masl
demHigh <- dem 
values(demHigh) <- ifelse(values(dem)>3000,1,NA)
rainfallHigh <- rainfall
rainfallHigh <- mask(rainfallHigh,demHigh)
plot(rainfallHigh)
summaryRaster(rainfallHigh,categorical=F)
#rainfall for altitudes lower >3k masl
demLow <- dem 
values(demLow) <- ifelse(values(dem)<3000,1,NA)
rainfallLow <- rainfall
rainfallLow <- mask(rainfallLow,demLow)
plot(rainfallLow)
summaryRaster(rainfallLow,categorical=F)
#get r-quared for historical rainfall (histPrec, histR2) and for models (scenPrec, scewnR2)
fold.ls <- list.files(paste0(data.folder,"/RUSLE"),full.names=T)
rain.data <- lapply(fold.ls,function(x){
  x.hist <- list.files(paste0(x,"/r_rainfall/rainfall"),recursive=T,full.names=T,pattern="historical.tif$")
  x.hist <-  mean(values(raster(x.hist)),na.rm=T)
  x.hist.rsqr <- readRDS(list.files(paste0(x,"/r_rainfall/models"),recursive=T,full.names=T,pattern="accuracy_historical.rds$"))
  x.scen <- list.files(paste0(x,"/r_rainfall/rainfall"),recursive=T,full.names=T,pattern="ssp585_2081-2100.tif$")
  x.scen <-  mean(values(raster(x.scen)),na.rm=T)
  x.scen.rsqr <- list.files(paste0(x,"/r_rainfall/models"),recursive=T,full.names=T,pattern="_ssp585_2081-2100.rds$")
  x.scen.rsqr <- readRDS(grep("accuracy_",x.scen.rsqr,value=T))
  return(data.frame(model=basename(x),
                    histPrec=x.hist,histR2=x.hist.rsqr$R2,
                    scenPrec=x.scen,scenR2=x.scen.rsqr$R2)
  )
})
rain.data <- do.call("rbind.data.frame",rain.data)
rain.data
#get r-factor
rfac.data <- lapply(fold.ls,function(x){
  x.hist <- list.files(paste0(x,"/r_rainfall/Rfactor"),recursive=T,full.names=T,pattern="historical.tif$")
  x.hist <-  mean(values(raster(x.hist)),na.rm=T)
  x.scen <- list.files(paste0(x,"/r_rainfall/Rfactor"),recursive=T,full.names=T,pattern="ssp585_2081-2100.tif$")
  x.scen <-  mean(values(raster(x.scen)),na.rm=T)
  return(data.frame(model=basename(x),hist=x.hist,scen=x.scen))
})
rfac.data <- do.call("rbind.data.frame",rfac.data)
rfac.data
#collect all r-square from models and standard deviation (SD)
rsqr.data <- lapply(fold.ls,function(x){
  x.scen.rsqr <- list.files(paste0(x,"/r_rainfall/models"),recursive=T,full.names=T,pattern="accuracy")
  x.scen.rsqr <- as.numeric(sapply(x.scen.rsqr[-1],function(y){
    return(readRDS(y)$R2)
  }))
  return(data.frame(model=basename(x),r2=mean(x.scen.rsqr),sd=sd(x.scen.rsqr)))
})
rsqr.data <- do.call("rbind.data.frame",rsqr.data)
rsqr.data

#### 3.2 Soil erodibility factor (K-factor) ####

#soil texture fractions downscale r-squared estimations
acc.data <- list.files(paste0(data.folder,"/RUSLE/CNRM-CM6-1/k_erodability/models"),pattern="accuracy",full.names = T)
acc.data <- lapply(acc.data,function(x){
  x.name <- gsub(".rds","",unlist(strsplit(basename(x),"_"))[2])
  return(data.frame(name=x.name,r2=readRDS(x)$R2))
})
acc.data <- do.call("rbind.data.frame",acc.data)
acc.data
mean(acc.data$r2)
sd(acc.data$r2)
#report textures based in altitudes
study.area <- shapefile(paste0(data.folder,"/studyArea/pambamarca_area.shp"))
dem <- raster(paste0(data.folder,("/DEM/pambamarca_100mts_dem_crop.tif")))
tex.data <- stack(list.files(paste0(data.folder,"/RUSLE/CNRM-CM6-1/k_erodability/texture"),pattern=".tif$",full.names = T))
#evaluate texture for altitudes higher 3600 masl
demHigh <- dem 
values(demHigh) <- ifelse(values(dem)>3600,1,NA)
texHigh <- tex.data
texHigh <- mask(texHigh,demHigh)
clay.sum <- summaryRaster(texHigh$clay,categorical=F)
silt.sum <- summaryRaster(texHigh$silt,categorical=F)
sand.sum <- summaryRaster(texHigh$sand,categorical=F)
org.sum <- summaryRaster(texHigh$organic,categorical=F)
print("HIGH")
org.sum
texHigh <- data.frame(
  "CLAY"=as.numeric(clay.sum[,c(1,3:7)]),
  "SILT"=as.numeric(silt.sum[,c(1,3:7)]),
  "SAND"=as.numeric(sand.sum[,c(1,3:7)])
)
texHigh <- TT.normalise.sum(texHigh)
geo <- TT.plot(
  class.sys = "USDA-NCSS.TT",
  tri.data = texHigh,
  col="red",
  main = ""
)
#evaluate texture for altitudes lower 3600 masl
demLow <- dem 
values(demLow) <- ifelse(values(dem)<3600,1,NA)
texLow <- tex.data
texLow <- mask(texLow,demLow)
clay.sum <- summaryRaster(texLow$clay,categorical=F)
silt.sum <- summaryRaster(texLow$silt,categorical=F)
sand.sum <- summaryRaster(texLow$sand,categorical=F)
org.sum <- summaryRaster(texLow$organic,categorical=F)
print("LOW")
org.sum
texLow <- data.frame(
  "CLAY"=as.numeric(clay.sum[,c(1,3:7)]),
  "SILT"=as.numeric(silt.sum[,c(1,3:7)]),
  "SAND"=as.numeric(sand.sum[,c(1,3:7)])
)
texLow <- TT.normalise.sum(texLow)
TT.points(
  tri.data=texLow,
  geo=geo,
  col="blue",
  pch=3
)
#check k-factor statistics (all models use the same k-factor model)
kfactor <- raster(paste0(data.folder,"/RUSLE/CNRM-CM6-1/k_erodability/Kfactor/k_factor.tif"))
plot(kfactor)
mean(values(kfactor),na.rm=T)*10
sd(values(kfactor),na.rm=T)*10
range(values(kfactor),na.rm=T)*10

#### 3.3 Slope length (S-factor) and steepness factors (L-factor)  ####

study.area <- shapefile(paste0(data.folder,"/studyArea/pambamarca_area.shp"))
lsfactor <- raster(paste0(data.folder,"/RUSLE/CNRM-CM6-1/ls_slope/ls_factor.tif"))
slope <- terrain(raster(paste0(data.folder,("/DEM/pambamarca_100mts_dem_crop.tif"))),unit="degrees")
summaryRaster(lsfactor,categorical = F)
#for high slopes
slopeHigh <- slope 
values(slopeHigh) <- ifelse(values(slope)>25,1,NA)
lsfactorH <- lsfactor 
lsfactorH <- mask(lsfactorH,slopeHigh)
plot(lsfactorH)
summaryRaster(lsfactorH,categorical=F)
#for low >3k
slopeLow <- slope 
values(slopeLow) <- ifelse(values(slope)<25,1,NA)
lsfactorL <- lsfactor 
lsfactorL <- mask(lsfactorL,slopeLow)
plot(lsfactorL)
summaryRaster(lsfactorL,categorical=F)

#### 3.4 Land cover management factor (C-factor) ####

#check composite images dates
landsat.dates <- read.csv(paste0(data.folder,"/landsatComposites/compositesImagesDates.csv"))
landsat.dates <- landsat.dates[,2]
landsat.dates <- unlist(strsplit(gsub("[[]|[]]| ","",landsat.dates),","))
landsat.dates <- ymd(landsat.dates)
table(year(landsat.dates))
table(month(landsat.dates))
length(landsat.dates)
#compare C-factor with land cover
study.area <- shapefile(paste0(data.folder,"/studyArea/pambamarca_area.shp"))
landcover <-  shapefile(paste0(data.folder,"/characterization/cobertura2014_MAGAP.shp"))
landcover <- spTransform(landcover,CRSobj=CRS("+proj=longlat +datum=WGS84 +no_defs"))
cfactor <- raster(paste0(data.folder,"/RUSLE/CNRM-CM6-1/c_coverManagement/c_factor.tif"))
landras <- cfactor
values(landras) <- NA
landras <- rasterize(landcover,landras,field=landcover$COD1,fun=mean)
#forest
forest <- landras 
values(forest) <- ifelse(values(forest)==1,1,NA)
plot(forest)
cfactorF <- cfactor
cfactorF <- mask(cfactorF,forest)
plot(cfactorF)
summaryRaster(cfactorF,categorical=F)
#croplands
croplands <- landras 
values(croplands) <- ifelse(values(croplands)==2,1,NA)
plot(croplands)
cfactorC <- cfactor
cfactorC <- mask(cfactorC,croplands)
plot(cfactorC)
summaryRaster(cfactorC,categorical=F)
#shurblands and grasslands
grass <- landras 
values(grass) <- ifelse(values(grass)==3,1,NA)
plot(grass)
cfactorG <- cfactor
cfactorG <- mask(cfactorG,grass)
plot(cfactorG)
summaryRaster(cfactorG,categorical=F)

#### 3.5 Operation of RUSLE and validation of its results with the historical climate model ####

#rusle map AREAS
rusle.map <- raster(paste0(data.folder,"/RUSLE/CNRM-CM6-1/OPERATION/HISTORICAL/historical_RUSLE.tif"))
rusle.map <- projectRaster(rusle.map,res=c(100,100),crs=CRSargs(CRS("+init=epsg:32617")))
rusle.map <- as.data.frame(as(rusle.map, "SpatialPixelsDataFrame"))
names(rusle.map)[1] <- "value"
rusle.map <- reshape2::melt(rusle.map,id.vars=c("x","y"))
rusle.map$value <- cut(rusle.map$value,breaks=c(0,5,10,20,100,1025),
                       labels=c("0 - 5","5 - 10","10 - 20",
                                "20 - 100","> 100"))
rusle.map <- as.data.frame(table(rusle.map$value))
rusle.map$ha <- (rusle.map$Freq * 100 * 100) * 0.0001
rusle.map

#### 3.6 Priorization of Pucaras and identification of gully erosion ####

#prepare table
meta.data <- list.files(paste0(data.folder,"/metadata"),full.names=T,pattern="droneTable")
meta.data <- lapply(meta.data,function(x){
  x.name <- basename(x)
  x <- read.xlsx(x,1)
  x.df <- data.frame(name=x.name,
                     startTime=min(as_datetime(x$CreateDate)),
                     endTime=max(as_datetime(x$CreateDate)),
                     imgNumber=nrow(x),
                     height=paste0(round(mean(x$RelativeAltitude),1),"±",round(sd(x$RelativeAltitude),1)),
                     ISO=names(which.max(table(x$ISO))),
                     speed=names(which.max(table(x$ExposureTime))),
                     aperture=names(which.max(table(x$ApertureValue))),
                     angle=paste0(round(mean(x$GimbalPitchDegree),1),"±",round(sd(x$GimbalPitchDegree),1)))
  
})
meta.data <- do.call("rbind.data.frame",meta.data)
meta.data
sum(meta.data$imgNumber)
#absolute
mean(c(2.9,0.8, 1.1, 2.7, 0.7, 1.6))
mean(c(6.1,2.3,1.1,3.0,1.9,10.4))
#relative
mean(c(0.1,0.2,0.3,0.7,0.2,0.3))
mean(c(0.3,0.3,0.7,1.1,0.3, 0.6))

#### 4.1 Accuracy assessment, and evaluation of RUSLE with the historical climate data ####

#confusion  matrix inputs 
rusle.map <- raster(paste0(data.folder,"/RUSLE/CNRM-CM6-1/OPERATION/HISTORICAL/historical_RUSLE.tif"))
geophotos.shp <- shapefile(paste0(data.folder,"/validation/geophotos.shp"))
#get samples
geophotos.shp <- split(geophotos.shp,geophotos.shp$class)
set.seed(663)
geophotos.shp[[1]] <- geophotos.shp[[1]][sample(1:nrow(geophotos.shp[[1]]),20),]
geophotos.shp[[2]] <- geophotos.shp[[2]][sample(1:nrow(geophotos.shp[[2]]),20),]
geophotos.shp[[3]] <- geophotos.shp[[3]][sample(1:nrow(geophotos.shp[[3]]),20),]
geophotos.shp <- rbind(geophotos.shp[[1]],geophotos.shp[[2]],geophotos.shp[[3]])
#extract RUSLE map values to samples
val <- raster::extract(rusle.map,geophotos.shp)
geophotos.shp$RUSLE <- unlist(val)
geophotos.shp <- geophotos.shp[!is.na(geophotos.shp$RUSLE),]
#prepare clases
geophotos.shp$RUSLE_cat <- NA
geophotos.shp$RUSLE_cat[geophotos.shp$RUSLE<10] <- "Low"
geophotos.shp$RUSLE_cat[geophotos.shp$RUSLE>=10 & geophotos.shp$RUSLE<20] <- "Medium"
geophotos.shp$RUSLE_cat[geophotos.shp$RUSLE>=20] <- "High"
geophotos.data <- geophotos.shp@data
geophotos.data$class <- factor(geophotos.data$class,levels=c("Low","Medium","High"))
geophotos.data$RUSLE_cat <- factor(geophotos.data$RUSLE_cat,levels=c("Low","Medium","High"))
#derive confusion matrix
confusionMatrix(geophotos.data$RUSLE_cat,geophotos.data$class)
#get soil samples laboratory analysis
soil.data <- read.xlsx(paste0(data.folder,"/validation/lab_soil.xlsx"),2)
soil.data$RUSLE_soil <- soil.data$RUSLE
soil.data$RUSLE_soil <- factor(soil.data$RUSLE_soil,labels=c("0-1","1-5","5-10","10-20","20-100",">100"))
#sand fraction
sand.test <- cor.test(soil.data$sand,soil.data$RUSLE,method="spearman",exact=F)
sand.test$p.value < 0.05
sand.test$estimate
#silt fraction
silt.test <- cor.test(soil.data$silt,soil.data$RUSLE,method="spearman",exact=F)
silt.test$p.value < 0.05
silt.test$estimate
#clay fraction
clay.test <- cor.test(soil.data$clay,soil.data$RUSLE,method="spearman",exact=F)
clay.test$p.value < 0.05
clay.test$estimate
#get buffered pucaras areas
buffered.shp <- shapefile(paste0(data.folder,"/studyArea/pucaras_pol_buffer.shp"))
code.data <- buffered.shp@data
code.data$ID <- 1:nrow(code.data)
#extract
box.data <- raster::extract(rusle.map,buffered.shp,df=T)
box.data$code <- NA
box.data <- merge(box.data,code.data,by="ID")
box.data <- box.data[,c(2,5)]
names(box.data) <- c("value","code")
box.data <- split(box.data,box.data$code)
box.data <- lapply(box.data,function(x){
  x$color <- mean(x$value)
  return(x)
})
box.data <- do.call("rbind.data.frame",box.data)
median(box.data[box.data$code=="Pi25","value"])
box.data$color <- cut(box.data$color,breaks=c(4.445182e-02,-20,-10,-5,-1,1,5,10,20,1.196098e+03),
                      labels=c("< -20","-20 - -10","-10 - -5","-5 - -1","0 - 1","1 - 5","5 - 10","10 - 20","> 20"))
#prepare colors
colors.pal <- rev(colorRampPalette(brewer.pal(name="RdBu", n = 8))(10))
colors.pal[5] <- "white"
color_expression <- bquote('Location of '~italic(.("Pucaras")))
fill_expression <- expression("Mean annual soil loss (tons*ha"^-1~"*yrs"^-1~")")
#draw boxplots
q <- ggplot() +
  theme_bw() +
  theme(axis.title.x = ggtext::element_markdown(),
        legend.position = "none",
        panel.grid.minor = element_blank()) +
  geom_boxplot(data=box.data,aes(x=reorder(code,value),y=value,fill=color)) +
  scale_fill_manual(values=colors.pal[6:10]) +
  xlab("*Pucaras* (codes)") +
  ylab(fill_expression) 
q

#### 4.2 Future climate change and subsequent RUSLE predictions ####

#get RUSLE models
rusle.dir <- paste0(list.files(paste0(data.folder,"/RUSLE"),full.names=T),"/EVALUATION")
pucara.data <- lapply(rusle.dir,function(x){
  x <- readRDS(list.files(x,full.names=T)[2])
})
pucara.data <- do.call("rbind.data.frame",pucara.data)
pucara.data$model_f <- factor(pucara.data$model,levels=c("IPSL-CM6A-LR","CNRM-CM6-1","MIROC6","GFDL-ESM4"))
#short period names
pucara.data$period <- factor(as.character(pucara.data$period),labels = c("2040","2060","2080","2100"))
pucara.data$ssp <- factor(pucara.data$ssp,labels=c("1-2.6","2-4.5","3-7.0","5-8.5"))
pi25.future <- pucara.data[pucara.data$code=="Pi25",]
pi25.future <- pi25.future[pi25.future$period=="2100" & pi25.future$ssp=="5-8.5",]
pi25.future
#plot
y_expression <- expression("Mean annual soil loss increase/decrease (tons*ha"^-1~"*yrs"^-1~")")
p <- ggplot(pucara.data,aes(x=reorder(code,value),y=value,color=ssp,group=ssp)) +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5,angle=90),
        axis.title.x = ggtext::element_markdown(),
        legend.position = "bottom",
        panel.grid.minor = element_blank()) +
  geom_line(size=0.5) +
  geom_point(size=0.75) +
  geom_hline(yintercept=0, linetype="dashed", color = "black",size=0.5) +
  scale_color_manual(values=c("#7CAE00","#00BFC4","#C77CFF","#F8766D")) +
  xlab("*Pucaras* (codes)") +
  ylab(y_expression) + 
  labs(color="SSP") +
  facet_grid(model_f~period,scales="free")
p

#### 4.3 Pucara prioritization according to their rain erosion gain according to the future climate models ####

#get data
rusle.dir <- paste0(list.files(paste0(data.folder,"/RUSLE"),full.names=T),"/EVALUATION")
pucara.data <- lapply(rusle.dir,function(x){
  x <- readRDS(list.files(x,full.names=T)[2])
})
pucara.data <- do.call("rbind.data.frame",pucara.data)
pucara.data$model_f <- factor(pucara.data$model,levels=c("IPSL-CM6A-LR","CNRM-CM6-1","MIROC6","GFDL-ESM4"))
pucara.data$period <- factor(pucara.data$period,labels=c("2040","2060","2080","2100"))
#get ranking
rank.data <- split(pucara.data,pucara.data$model)
rank.data <- lapply(rank.data,function(x){
  x.ssp <- split(x,x$ssp)
  x.ssp <- lapply(x.ssp,function(y){
    y.code <- split(y,y$code)
    y.code <- lapply(y.code,function(z){
      z <- data.frame(code=unique(z$code),meanVal=mean(z$value),model=unique(z$model),ssp=unique(z$ssp),period=z[which.max(z$value),"period"])
      return(z)
    })
    y.code <- do.call("rbind.data.frame",y.code)
    y.code <- y.code[order(y.code$meanVal,decreasing =T),]
    y.code$rank <- 1:nrow(y.code)
    return(y.code)
  })
  x.ssp <- do.call("rbind.data.frame",x.ssp)
  return(x.ssp)
})
rank.data <- do.call("rbind.data.frame",rank.data)
#first three in all models and ssp by scenario averages
three.first <- rank.data[rank.data$rank %in% 1:3,]
three.first <- as.data.frame(table(three.first$code))
three.first <- three.first[order(three.first$Freq,decreasing=T),]
three.first
#get model stats
stats.first <- rank.data[rank.data$rank %in% 1:3,]
rownames(stats.first) <- NULL
stats.first <- as.data.frame.matrix(table(stats.first$code,stats.first$model))
stats.first <- stats.first[order(rowSums(stats.first),decreasing=T),]
stats.first$total <- rowSums(stats.first)
stats.first[nrow(stats.first)+1,] <- colSums(stats.first)
rownames(stats.first)[nrow(stats.first)] <- "total"
stats.first
#get ssp stats
ssp.first <- rank.data[rank.data$rank %in% 1:3,]
rownames(ssp.first) <- NULL
ssp.first <- as.data.frame.matrix(table(ssp.first$code,ssp.first$ssp))
ssp.first <- ssp.first[order(rowSums(ssp.first),decreasing=T),]
ssp.first[nrow(ssp.first)+1,] <- colSums(ssp.first)
rownames(ssp.first)[nrow(ssp.first)] <- "total"
ssp.first
#get scenario stats
scen.first <- rank.data[rank.data$rank %in% 1:3,]
rownames(scen.first) <- NULL
scen.first <- as.data.frame.matrix(table(scen.first$code,scen.first$period))
scen.first <- scen.first[order(rowSums(scen.first),decreasing=T),]
scen.first[nrow(scen.first)+1,] <- colSums(scen.first)
rownames(scen.first)[nrow(scen.first)] <- "total"
scen.first
#last three in all models and ssp by scenario averages
three.last <- rank.data[rank.data$rank %in% 13:15,]
three.last <- as.data.frame(table(three.last$code))
three.last <- three.last[order(three.last$Freq,decreasing=T),]
three.last

#### 4.4 Assessment of high spatial resolution DTMs and gully-prone areas in prioritized Pucaras ####

#may take time processing....

#tabulate terraces
#pi10
pi10.L <- diff(c(25,445,469,482,535,548)) #cut points
pi10.S <- diff(c(78,194,304,343))
pi10.tot <- c(pi10.L,pi10.S)
pi10.tot <- c(mean(pi10.tot),sd(pi10.tot),length(pi10.tot)+2)
terraces.df <- data.frame(code="pi10",dist=pi10.tot[1],sd=pi10.tot[2],number=pi10.tot[3])
#pi12
pi12.L <- diff(c(80,113,142,159))
pi12.S <- diff(c(8,47,96))
pi12.tot <- c(pi12.L,pi12.S)
pi12.tot <- c(mean(pi12.tot),sd(pi12.tot),length(pi12.tot)+2)
pi12.tot <- data.frame(code="pi12",dist=pi12.tot[1],sd=pi12.tot[2],number=pi12.tot[3])
terraces.df <- rbind(terraces.df,pi12.tot)
#pi16
pi16.L <- diff(c(15,37,61,84,125,144))
pi16.S <- diff(c(9,21,62,84,115,129))
pi16.tot <- c(pi16.L,pi16.S)
pi16.tot <- c(mean(pi16.tot),sd(pi16.tot),length(pi16.tot)+2)
pi16.tot <- data.frame(code="pi16",dist=pi16.tot[1],sd=pi16.tot[2],number=pi16.tot[3])
terraces.df <- rbind(terraces.df,pi16.tot)
#pi19
pi19.L <- diff(c(9,29,93,122,150,190))
pi19.S <- diff(c(14,36,75,98,122))
pi19.tot <- c(pi19.L,pi19.S)
pi19.tot <- c(mean(pi19.tot),sd(pi19.tot),length(pi19.tot)+2)
pi19.tot <- data.frame(code="pi19",dist=pi19.tot[1],sd=pi19.tot[2],number=pi19.tot[3])
terraces.df <- rbind(terraces.df,pi19.tot)
#pi20
pi20.L <- diff(c(16,69,85,105,128,149,161,186,207,249))
pi20.S <- diff(c(55,66,84,99,133,142,153,172))
pi20.tot <- c(pi20.L,pi20.S)
pi20.tot <- c(mean(pi20.tot),sd(pi20.tot),length(pi20.tot)+2)
pi20.tot <- data.frame(code="pi20",dist=pi20.tot[1],sd=pi20.tot[2],number=pi20.tot[3])
terraces.df <- rbind(terraces.df,pi20.tot)
#pi25
pi25.L <- diff(c(46,92,182,400,502,594))
pi25.S <- diff(c(94,129,175,321,380,420))
pi25.tot <- c(pi25.L,pi25.S)
pi25.tot <- c(mean(pi25.tot),sd(pi25.tot),length(pi25.tot)+2)
pi25.tot <- data.frame(code="pi25",dist=pi25.tot[1],sd=pi25.tot[2],number=pi25.tot[3])
terraces.df <- rbind(terraces.df,pi25.tot)
#review axis distances
axis.shp <- st_read(paste0(data.folder,"/profiles/profiles.shp"))
axis.shp <- split(axis.shp,axis.shp$code)
axis.shp <- lapply(axis.shp,function(x){
  x.df <- as.numeric(st_length(x))
  x.df <- data.frame(code=unique(x$code),short=x.df[1],long=x.df[2],total=sum(x.df))
  return(x.df)
})
axis.shp <- do.call("rbind.data.frame",axis.shp)
rownames(axis.shp) <- NULL
axis.shp
#evident area
evident.shp <- st_read(paste0(data.folder,"/studyArea/pucaras_pol_UTM.shp"))
evident.shp <- evident.shp[evident.shp$code %in% c("Pi10","Pi12","Pi16","Pi19","Pi20","Pi25"),]
evident.shp <- evident.shp[order(match(evident.shp$code,c("Pi10","Pi12","Pi16","Pi19","Pi20","Pi25"))),]
evident.shp$area <- st_area(evident.shp)
#stream power
stream.files <- list.files(paste0(data.folder,"/streamPower"),full.names=T,pattern=".tif$")
stream.files <- lapply(1:6,function(x){
  x.ras <- raster(stream.files[x])
  x.val <- values(x.ras)
  x.val <- ifelse(x.val>=30,30,x.val)
  x.val <- ifelse(x.val>=20 & x.val<30,20,x.val)
  x.val <- ifelse(x.val<20 ,NA,x.val)
  values(x.ras) <- x.val
  return(x.ras)
})
#tabulate
gully.tab <- lapply(1:6,function(x){
  evident.x <- evident.shp[x,]
  stream.x <- stream.files[[x]]
  df.x <- extract(stream.x,evident.x,df=T)
  df.x <- table(df.x)
  df.x <- as.data.frame.matrix(df.x)
  df.x$`20` <- df.x$`20`*res(stream.x)[1]*res(stream.x)[2]
  df.x$`30` <- df.x$`30`*res(stream.x)[1]*res(stream.x)[2]
  df.x$per20 <-  sum(df.x[,1])*100/  as.numeric(evident.x$area)
  df.x$per30 <-  sum(df.x[,2])*100/  as.numeric(evident.x$area)
  df.x$perTot <- sum(df.x$per20+df.x$per30)
  df.x$areaTot <- df.x$`20`+df.x$`30`
  df.x$code <- evident.x$code
  return(df.x)
})
gully.tab <- do.call("rbind.data.frame",gully.tab)
gully.tab <- gully.tab[order(gully.tab$perTot,decreasing = T),]
gully.tab
terraces.df

#### 5.3 Identified gully-prone areas and opportunities using UAV in heritage preservation ####

#compare DEMS from Pi25
options(scipen=5)
dem.pi16 <- list.files(paste0(data.folder,"/DEM"),full.names=T,pattern="streamPower.tif$")
dem.pi16 <- lapply(dem.pi16,function(x){
  x.ras <- as.data.frame(as(raster(x), "SpatialPixelsDataFrame"))
  x.sum <- c(mean(x.ras[,1]),sd(x.ras[,1]))
  x.nam <- unlist(strsplit(basename(x),"_"))[2]
  x.sum <- data.frame(name=x.nam,
                      mean=as.numeric(x.sum[1]),
                      sd=as.numeric(x.sum[2]))
  x.ras <- reshape2::melt(x.ras,id.vars=c("x","y"))
  x.ras[,3] <- x.nam
  x <- list(x.sum,x.ras)
  return(x)
})
dem.pi16.sum <- do.call("rbind.data.frame",lapply(dem.pi16,"[[",1))
dem.pi16.sum
#map (annex)
dem.pi16 <- lapply(dem.pi16,"[[",2)
dem.pi16[[3]]$variable <- "SRTM (100 mts)"
dem.pi16[[1]]$variable <- "ASTER (30 mts)"
dem.pi16[[2]]$variable <- "UAV (0.15 mts)"
#plot
plot_func <- function(df,colorPal="Blues",fillExpression="test",inverse=F,legend=T,axis.coordinates=T) {
  extreme.vals <- quantile(df$value,probs=c(0.01,0.99),na.rm=T)
  if(inverse){
    dir <- -1
  }else{
    dir <- 1
  }
  for.scale <- st_as_sf(df,coords=c("x","y"),crs=4326)
  df$value[df$value <= extreme.vals[1]] <- extreme.vals[1]
  df$value[df$value >= extreme.vals[2]] <- extreme.vals[2]
  p <- ggplot() +
    theme_bw() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
          strip.background = element_rect(color = "black", size = 0.5),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_line(size=0.5,color="black"),
          axis.ticks.length = unit(.25, "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(color = NA, fill = NA),
          aspect.ratio = 1,
          plot.margin=unit(c(0,0.2,0,0.2),"cm"),
          legend.position=ifelse(legend,"bottom","none"),
          legend.key.width=unit(1, "cm"),
          legend.box.just = "left",
          legend.key=element_rect(colour="black")) +
    #map
    new_scale("fill") +
    geom_raster(data=df, aes(x=x, y=y, fill=value)) +
    #for scale
    geom_sf(data=for.scale,alpha=0,show.legend=F) +
    #for legend
    scale_fill_distiller(name=fillExpression,
                         palette = colorPal,
                         direction=dir,
                         trans = 'log10',
                         na.value ="white") +
    #grid
    scale_x_continuous(breaks = seq(-78.1982, -78.1934, by = 0.002),expand=c(0,0)) +
    scale_y_continuous(breaks = seq(-0.0124, -0.0079, by = 0.002),expand=c(0,0)) +
    #other
    annotation_scale(location="bl",text_col="black",style="ticks",line_col="black",pad_y = unit(0.1, "cm"),) +
    guides(fill = guide_legend(title.position="top",ncol=1)) +
    facet_wrap(variable~.)
  #for coordinates axis
  if(axis.coordinates){
    p <- p + theme(
      axis.text=element_text(size=8),
      axis.text.x = element_text(angle=90))
  }
  return(p)
}

fill_expression <- ""
pi16.aster <- plot_func(df=dem.pi16[[1]],colorPal="Greys",inverse=T,fillExpression=fill_expression,legend=T)

fill_expression <- expression("Stream power index (watts*mts"^-2~")")
pi16.srtm <- plot_func(df=dem.pi16[[3]],colorPal="Greys",inverse=T,fillExpression=fill_expression,legend=T)

fill_expression <- ""
pi16.drone <- plot_func(df=dem.pi16[[2]],colorPal="Greys",inverse=T,fillExpression=fill_expression,legend=T)

q <- ggarrange(pi16.srtm,pi16.aster,pi16.drone,common.legend = F,legend="bottom",ncol=3,nrow=1,align = "hv") 
output_folder <- "E:/Indoamerica/13_Manuscritos/33_Pambamarca_paper/figures"
ggsave(paste0(output_folder,"/Annex1_v2.jpg"),q,units="cm",width=25,height=20,bg="white")




#### 5.4 Recommendations for Pucaras management in the context of climate change ####

#landcover tabulation
landcover <- shapefile(paste0(data.folder,"/studyArea/pucaras_buffer_landCover.shp"))
landcover$area <- area(landcover)
landcover$intervened <- factor(landcover$cobertura_)
landcover$intervened <- fct_collapse(landcover$intervened,
                                     NATURAL = c("BOSQUE","VEGETACION ARBUSTIVA Y HERBACEA"),
                                     INTERVENED = c("TIERRA AGROPECUARIA","ZONA ANTROPICA"))
landcover <- split(landcover@data,landcover$code)
landcover <- lapply(landcover,function(x){
  x <- x[,c("code","intervened","cobertura0","area")]
  if(all(x$intervened=="NATURAL")){
    x$Intper <- 0
  }else if(all(x$intervened=="INTERVENED")){
    x$Intper <- 100
  }else{
    x$Intper <- (sum(x$area[x$intervened=="INTERVENED"]) * 100) / sum(x$area)
  }
  print(x)
  return(x)
})
landcover <- do.call("rbind.data.frame",landcover)
rownames(landcover) <- NULL
landcover

