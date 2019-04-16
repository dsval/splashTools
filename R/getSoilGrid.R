#' getSoilGrid
#'
#' Downloads soil data (sand, clay, OM, gravel, bulk density and depth) from soilgrids.org in raster format, calculates the weigthed mean according layer depths
#' @param  bbox, any type of spatial object
#' @param  ouputdir, directory where to save the files
#' @import raster
#' @import rgdal    
#' @keywords splashTools
#' @export
#' @examples
#' splash.grid()


getSoilGrid<-function(bbox,ouputdir=getwd()){
# function to gather point data from soilgrids.org, the function also does a   
# require(raster)
# require(rgdal)
# weigthed average according to the point soil depth to get a value for the whole soil column
# testing:
# bbox<-readOGR("C:/Water_Data/iMHEA/iMHEA_geodata/JTU/JTU_04", "iMHEA_JTU_04_catchment")
# end testing
# define spatial point
if(class(bbox)=="SpatialPolygonsDataFrame"){
	limits<-bbox@bbox
}else if(class(bbox)=="RasterLayer"||class(bbox)=="RasterBrick"||class(bbox)=="RasterStack"){
	limits<-as.matrix(extent(bbox))
}else{
	limits<-bbox
}


setwd(dirname(rasterTmpFile()))
# define filenames to dwld
url_sg<-"http://85.214.241.121:8080/geoserver/ows?service=WCS&version=2.0.1&request=GetCoverage&CoverageId="
depth_Id<-"BDRICM_M_250m"
sand_Ids<-c("SNDPPT_M_sl1_250m","SNDPPT_M_sl2_250m","SNDPPT_M_sl3_250m","SNDPPT_M_sl4_250m","SNDPPT_M_sl5_250m","SNDPPT_M_sl6_250m","SNDPPT_M_sl7_250m")	
clay_Ids<-c("CLYPPT_M_sl1_250m","CLYPPT_M_sl2_250m","CLYPPT_M_sl3_250m","CLYPPT_M_sl4_250m","CLYPPT_M_sl5_250m","CLYPPT_M_sl6_250m","CLYPPT_M_sl7_250m")
OM_Ids<-c("ORCDRC_M_sl1_250m","ORCDRC_M_sl2_250m","ORCDRC_M_sl3_250m","ORCDRC_M_sl4_250m","ORCDRC_M_sl5_250m","ORCDRC_M_sl6_250m","ORCDRC_M_sl7_250m")
gravel_Ids<-c("CRFVOL_M_sl1_250m","CRFVOL_M_sl2_250m","CRFVOL_M_sl3_250m","CRFVOL_M_sl4_250m","CRFVOL_M_sl5_250m","CRFVOL_M_sl6_250m","CRFVOL_M_sl7_250m")		
bd_Ids<-c("BLDFIE_M_sl1_250m","BLDFIE_M_sl2_250m","BLDFIE_M_sl3_250m","BLDFIE_M_sl4_250m","BLDFIE_M_sl5_250m","BLDFIE_M_sl6_250m","BLDFIE_M_sl7_250m")



# define raster range for dwld
latmin<-limits[2,1]-0.01
latmax<-limits[2,2]+0.01
lonmin<-limits[1,1]-0.01
lonmax<-limits[1,2]+0.01
subset<-paste0("&subset=Long(",lonmin,",",lonmax,")&subset=Lat(",latmin,",",latmax,")")
# get depth data in meters	
download.file(paste0(url_sg,depth_Id,subset),"depth.tif",method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
depth<-raster("depth.tif")/100

# function to compute the weigthed average
avg.soil.prop<-function(soil_prop,depth){
	interv<-c(0.05,0.1,0.15,0.3,0.4,1,0.5)
	w1<-interv/sum(interv,na.rm = TRUE)
	w2<-interv/sum(interv[1:6],na.rm = TRUE)
	w2[7]<-0
	w3<-interv/sum(interv[1:5],na.rm = TRUE)
	w3[6:7]<-0
	w4<-interv/sum(interv[1:4],na.rm = TRUE)
	w4[5:7]<-0
	w5<-interv/sum(interv[1:3],na.rm = TRUE)
	w5[4:7]<-0
	w6<-interv/sum(interv[1:2],na.rm = TRUE)
	w6[3:7]<-0
	w7<-interv/interv[1]
	w7[2:7]<-0
	z<-ifelse(depth>=2,soil_prop%*%w1,ifelse(depth>1,soil_prop%*%w2,ifelse(depth>0.6,soil_prop%*%w3,ifelse(depth>0.3,soil_prop%*%w4,ifelse(depth>0.15,soil_prop%*%w5,ifelse(depth>0.05,soil_prop%*%w6,ifelse(depth>0,soil_prop%*%w7,soil_prop)))))))
	return(z)
}



# get sand data
download.file(paste0(url_sg,sand_Ids[1],subset),paste0(sand_Ids[1],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
sand<-brick(paste0(sand_Ids[1],".tif"))
for(i in 2:7){
	
	download.file(paste0(url_sg,sand_Ids[i],subset),paste0(sand_Ids[i],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
	sand<-addLayer(sand,raster(paste0(sand_Ids[i],".tif")))	
}
sand<-overlay(sand,depth,fun=avg.soil.prop,filename="sand_ave.grd",overwrite=TRUE,mode = "wb")

# get clay data
download.file(paste0(url_sg,clay_Ids[1],subset),paste0(clay_Ids[1],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
clay<-brick(paste0(clay_Ids[1],".tif"))
for(i in 2:7){
	
	download.file(paste0(url_sg,clay_Ids[i],subset),paste0(clay_Ids[i],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
	clay<-addLayer(clay,raster(paste0(clay_Ids[i],".tif")))	
}

clay<-overlay(clay,depth,fun=avg.soil.prop,filename="clay_ave.grd",overwrite=TRUE)

# get OM data
download.file(paste0(url_sg,OM_Ids[1],subset),paste0(OM_Ids[1],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
OM<-brick(paste0(OM_Ids[1],".tif"))
for(i in 2:7){
	
	download.file(paste0(url_sg,OM_Ids[i],subset),paste0(OM_Ids[i],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
	OM<-addLayer(OM,raster(paste0(OM_Ids[i],".tif")))	
}

# transform from g/Kg Org C to OM in percent
OM<-overlay(OM*1.724/10,depth,fun=avg.soil.prop,filename="OM_ave.grd",overwrite=TRUE)

# get gravel data
download.file(paste0(url_sg,gravel_Ids[1],subset),paste0(gravel_Ids[1],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
gravel<-brick(paste0(gravel_Ids[1],".tif"))
for(i in 2:7){
	
	download.file(paste0(url_sg,gravel_Ids[i],subset),paste0(gravel_Ids[i],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
	gravel<-addLayer(gravel,raster(paste0(gravel_Ids[i],".tif")))	
}

gravel<-overlay(gravel,depth,fun=avg.soil.prop,filename="gravel_ave.grd",overwrite=TRUE)


# get pH data
download.file(paste0(url_sg,bd_Ids[1],subset),paste0(bd_Ids[1],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
bd<-brick(paste0(bd_Ids[1],".tif"))
for(i in 2:7){
	
	download.file(paste0(url_sg,bd_Ids[i],subset),paste0(bd_Ids[i],".tif"),method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
	bd<-addLayer(bd,raster(paste0(bd_Ids[i],".tif")))	
}
# transform to g/cm3
bd<-overlay(bd/1000,depth,fun=avg.soil.prop,filename="bd_ave.grd",overwrite=TRUE)




soil_phys<-stack(sand,clay,OM,gravel,bd,depth)


soil_phys<-writeRaster(soil_phys,paste0(ouputdir,"/","soil_data",".nc"),format="CDF",overwrite=TRUE,varname="soil_data", varunit="various", longname="soil properties:depth,sand,clay,OM,gravel,bd", xname="lon", yname="lat")
names(soil_phys)<-c("sand","clay","OM","gravel","bd","depth")
return(soil_phys)

}
