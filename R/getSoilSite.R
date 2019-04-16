#' getSoilSite
#'
#' Get soil data (sand, clay, OM, gravel, bulk density and depth) from soilgrids REST API, point data is from the 250m resolution dataset, calculates the weigthed mean according layer depths
#' @param   lat, lon
#' @import httr 
#' @keywords splash
#' @export
#' @examples
#' splash.grid()


getSoilSite<-function(lat,lon){
	# require(httr)
	sg_query<-paste0("https://rest.soilgrids.org/query?lon=",lon,"&lat=",lat)
	sgquery<- GET(sg_query)
	sg_all<-content(sgquery)
	depth<-as.numeric(sg_all$properties$BDRICM$M$BDRICM_M)/100
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
sand<-as.numeric(sg_all$properties$SNDPPT$M)
sand<-avg.soil.prop(sand,depth)
clay<-as.numeric(sg_all$properties$CLYPPT$M)
clay<-avg.soil.prop(sand,depth)
# transform from g/Kg OC to % OM
OM<-as.numeric(sg_all$properties$ORCDRC$M)*1.724/10
OM[OM<0]<-0
OM<-avg.soil.prop(OM,depth)
gravel<-as.numeric(sg_all$properties$CRFVOL$M)
gravel<-avg.soil.prop(gravel,depth)
bd<-as.numeric(sg_all$properties$BLDFIE$M)/1000
bd<-avg.soil.prop(bd,depth)
result<-c(sand,clay,OM,gravel,bd,depth)
names(result)<-c("sand","clay","OM","gravel","bd","depth")
return(result)
	
}
