#' getSoilSite
#'
#' Get soil data from soilgrids REST API, point data is from the 250m resolution dataset, calculates the weigthed mean according layer depths
#' @param   lat, lon
#' @param  global_depth (optional), global raster of depth to the R horizon, unfortunately not available anymore in the new update of soilgrids, so, download it separately from https://files.isric.org/soilgrids/former/2017-03-10/data/BDRICM_M_250m_ll.tif
#' @return a numeric vector with percentages of sand (w/w), clay(w/w), SOM(w/w), gravel(v/v); bulk density (g/cm3) and depth (m)
#' @import httr 
#' @keywords soil texture
#' @export
#' @examples
#' getSoilSite(lat=8.7333,lon=-70.883333) 


getSoilSite<-function(lat,lon,global_depth=NULL,set='phys'){
	# require(httr)
	sg_query<-paste0("https://rest.soilgrids.org/soilgrids/v2.0/properties/query?lon=",lon,"&lat=",lat,'&property=bdod&property=cfvo&property=clay&property=sand&property=soc&depth=0-5cm&depth=0-30cm&depth=5-15cm&depth=15-30cm&depth=30-60cm&depth=60-100cm&depth=100-200cm&value=mean')

	sgquery<- GET(sg_query)
	sg_all<-content(sgquery)
	if(!is.null(global_depth)){
		loc<-SpatialPoints(coords = cbind(x=lon,y=lat),proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
		depth<-raster::extract(global_depth,loc)/100
	}else(
		depth=2.0)
	
	# function to compute the weigthed average
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

read_f_api<-function(response,factor){
	soilprop<-unlist(response)
	soilprop<-soilprop[names(soilprop)=='depths.values.mean']
	soilprop<-(as.numeric(soilprop))*factor
	soilprop
}
#read data and convert to percentages
sand<-read_f_api(sg_all$properties$layers[[4]],0.1)
sand<-avg.soil.prop(sand,depth)
clay<-read_f_api(sg_all$properties$layers[[3]],0.1)
clay<-avg.soil.prop(clay,depth)
# transform from g/Kg OC to % OM
OM<-read_f_api(sg_all$properties$layers[[5]],(1.724/100))
OM<-avg.soil.prop(OM,depth)

gravel<-read_f_api(sg_all$properties$layers[[2]],0.1)
gravel<-avg.soil.prop(gravel,depth)

bd<-read_f_api(sg_all$properties$layers[[1]],0.01)
bd<-avg.soil.prop(bd,depth)

result<-c(sand,clay,OM,gravel,bd,depth)
names(result)<-c("sand","clay","OM","gravel","bd","depth")
return(result)
	
}
