#' getSoilGrid
#'
#' Downloads soil data (sand, clay, OM, gravel, bulk density and depth) from soilgrids.org in raster format, calculates the weigthed average according to the layer depths
#' @param  bbox, any type of spatial object
#' @param  global_depth, global raster of depth to the R horizon, unfortunately not available anymore in the new update of soilgrids, so, download it separately fromhttps://files.isric.org/soilgrids/former/2017-03-10/data/BDRICM_M_250m_ll.tif
#' @param  ouputdir, directory where to save the files
#' @return  A RasterBrick object with the layers organized as sand(perc),clay(perc),organic matter(perc),coarse-fragments-fraction(perc), bulk density(g cm-3) and depth (m)
#' @import raster
#' @import rgdal    
#' @keywords Soilgrids, raster
#' @export
#' @examples
#' getSoilGrid()


getSoilGrid<-function(bbox,global_depth,ouputdir=getwd()){
# function to gather raster data from soilgrids.org, the function also does a   
# require(raster)
# require(rgdal)
# weigthed average according to the point soil depth to get a value for the whole soil column
# testing:
# bbox<-readOGR("C:/Water_Data/iMHEA/iMHEA_geodata/JTU/JTU_04", "iMHEA_JTU_04_catchment")
# end testing
###############################################################################
# 01. define the bounding box
###############################################################################
if(class(bbox)=="SpatialPolygonsDataFrame"){
	limits<-bbox@bbox
}else if(class(bbox)=="RasterLayer"||class(bbox)=="RasterBrick"||class(bbox)=="RasterStack"){
	limits<-as.matrix(extent(bbox))
}else{
	limits<-bbox
}

latmin<-limits[2,1]-0.01
latmax<-limits[2,2]+0.01
lonmin<-limits[1,1]-0.01
lonmax<-limits[1,2]+0.01
subset<-paste0("&FORMAT=image/tiff&SUBSET=long(",lonmin,",",lonmax,")&SUBSET=lat(",latmin,",",latmax,")")

bbox<-extent(c(lonmin,lonmax,latmin,latmax))

###############################################################################
# 02. crop the depth, convert from cm to m the bounding box
###############################################################################
depth<-crop(global_depth,bbox)/100
###############################################################################
# 03. prepare the urls
###############################################################################

url_sg1<-"https://maps.isric.org/mapserv?map=/map/"
url_sg2<-'.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID='
url_sg3<-'&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/152160'
sand_Ids<-c("sand_0-5cm_mean","sand_5-15cm_mean","sand_15-30cm_mean","sand_30-60cm_mean","sand_60-100cm_mean","sand_100-200cm_mean")	
clay_Ids<-c("clay_0-5cm_mean","clay_5-15cm_mean","clay_15-30cm_mean","clay_30-60cm_mean","clay_60-100cm_mean","clay_100-200cm_mean")
OM_Ids<-c("soc_0-5cm_mean","soc_5-15cm_mean","soc_15-30cm_mean","soc_30-60cm_mean","soc_60-100cm_mean","soc_100-200cm_mean")
gravel_Ids<-c("cfvo_0-5cm_mean","cfvo_5-15cm_mean","cfvo_15-30cm_mean","cfvo_30-60cm_mean","cfvo_60-100cm_mean","cfvo_100-200cm_mean")
bd_Ids<-c("bdod_0-5cm_mean","bdod_5-15cm_mean","bdod_15-30cm_mean","bdod_30-60cm_mean","bdod_60-100cm_mean","bdod_100-200cm_mean")
###############################################################################
# 04. define the function to compute the weigthed average
###############################################################################
 
avg.soil.prop<-function(soil_prop,depth){
	interv<-c(0.05,0.1,0.15,0.3,0.4,1)
	
	w1<-interv/sum(interv,na.rm = TRUE)
	w2<-interv/sum(interv[1:5],na.rm = TRUE)
	w2[6]<-0
	w3<-interv/sum(interv[1:4],na.rm = TRUE)
	w3[5:6]<-0
	w4<-interv/sum(interv[1:3],na.rm = TRUE)
	w4[4:6]<-0
	w5<-interv/sum(interv[1:2],na.rm = TRUE)
	w5[3:6]<-0
	w6<-interv/interv[1]
	w6[2:6]<-0
	
	z<-ifelse(depth>1,soil_prop%*%w1,ifelse(depth>0.6,soil_prop%*%w2,ifelse(depth>0.3,soil_prop%*%w3,ifelse(depth>0.15,soil_prop%*%w4,ifelse(depth>0.05,soil_prop%*%w5,ifelse(depth>0.0,soil_prop%*%w6))))))
	return(z)
}
###############################################################################
# 05. download the files
###############################################################################
setwd(dirname(rasterTmpFile()))
# get sand, 1st layer g/kg
download.file(paste0(url_sg1,'sand',url_sg2,sand_Ids[1],subset,url_sg3),paste0(sand_Ids[1],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
sand<-brick(paste0(sand_Ids[1],".tif"))
# get clay, 1st layer g/kg
download.file(paste0(url_sg1,'clay',url_sg2,clay_Ids[1],subset,url_sg3),paste0(clay_Ids[1],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
clay<-brick(paste0(clay_Ids[1],".tif"))
# get soil organic carbon, 1st layer dg/kg
download.file(paste0(url_sg1,'soc',url_sg2,OM_Ids[1],subset,url_sg3),paste0(OM_Ids[1],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
SOC<-brick(paste0(OM_Ids[1],".tif"))
# get volumetric coarse fraction, 1st layer cm3/dm3
download.file(paste0(url_sg1,'cfvo',url_sg2,gravel_Ids[1],subset,url_sg3),paste0(gravel_Ids[1],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
gravel<-brick(paste0(gravel_Ids[1],".tif"))
# get volumetric coarse fraction, 1st layer cg/dm3
download.file(paste0(url_sg1,'bdod',url_sg2,bd_Ids[1],subset,url_sg3),paste0(bd_Ids[1],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
bd<-brick(paste0(bd_Ids[1],".tif"))

for(i in 2:6){
	
	download.file(paste0(url_sg1,'sand',url_sg2,sand_Ids[i],subset,url_sg3),paste0(sand_Ids[i],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
	sand<-addLayer(sand,raster(paste0(sand_Ids[i],".tif")))
	
	download.file(paste0(url_sg1,'clay',url_sg2,clay_Ids[i],subset,url_sg3),paste0(clay_Ids[i],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
	clay<-addLayer(clay,raster(paste0(clay_Ids[i],".tif")))	
	
	download.file(paste0(url_sg1,'soc',url_sg2,OM_Ids[i],subset,url_sg3),paste0(OM_Ids[i],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
	SOC<-addLayer(SOC,raster(paste0(OM_Ids[i],".tif")))
	
	download.file(paste0(url_sg1,'cfvo',url_sg2,gravel_Ids[i],subset,url_sg3),paste0(gravel_Ids[i],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
	gravel<-addLayer(gravel,raster(paste0(gravel_Ids[i],".tif")))
	
	download.file(paste0(url_sg1,'bdod',url_sg2,bd_Ids[i],subset,url_sg3),paste0(bd_Ids[i],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
	bd<-addLayer(bd,raster(paste0(bd_Ids[i],".tif")))
	
}

###############################################################################
# 05. calc the weighted average, project to wgs and convert units
###############################################################################
sand<-overlay(projectRaster(sand,depth),depth,fun=avg.soil.prop,filename="sand_ave.grd",overwrite=TRUE)
#conversion to percentage g/100g
sand<-calc(sand,function(x){x*0.1})

clay<-overlay(projectRaster(clay,depth),depth,fun=avg.soil.prop,filename="clay_ave.grd",overwrite=TRUE)
#conversion to percentage g/100g
clay<-calc(clay,function(x){x*0.1})

SOC<-overlay(projectRaster(SOC,depth),depth,fun=avg.soil.prop,filename="soc_ave.grd",overwrite=TRUE)
#conversion to percentage g/100g
SOM<-calc(SOC,function(x){x*1.724/100})

gravel<-overlay(projectRaster(gravel,depth),depth,fun=avg.soil.prop,filename="gravel_ave.grd",overwrite=TRUE)
#conversion to percentage m3/100m3
gravel<-calc(gravel,function(x){x*0.1})

bd<-overlay(projectRaster(bd,depth),depth,fun=avg.soil.prop,filename="bd_ave.grd",overwrite=TRUE)
#conversion to percentage g/cm3
bd<-calc(bd,function(x){x*0.01})


soil_phys<-stack(sand,clay,SOM,gravel,bd,depth)


soil_phys<-writeRaster(soil_phys,paste0(ouputdir,"/","soil_data",".nc"),format="CDF",overwrite=TRUE,varname="soil_data", varunit="various", longname="soil properties:depth,sand,clay,SOM,gravel,bd", xname="lon", yname="lat")
names(soil_phys)<-c("sand","clay","OM","gravel","bd","depth")
return(soil_phys)

}
