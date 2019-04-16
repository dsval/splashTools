#' getSoilNRCS
#'
#' wrapper of the soilDB to retrievesoil data (sand, clay, OM, gravel, bulk density and depth) from the NRCS sites
#' @param   coords vector c(lat,lon)
#' @param   maxdepth maximum depth in  cm, to calculate the weigthed average, if maxdepth=NULL (default) the function uses all the available depths

#' @importFrom soilDB fetchKSSL
#' @keywords splash
#' @export
#' @examples
#' splash.grid()

getSoilNRCS<-function(coords,maxdepth=NULL){
	# wrapper from soildb to match snotel and scan stations to soil profiles
	# it searches in 10 Km range, ussually every stations has a pedon
	# inputs
	# coords vector c(lat,lon)
	# testing
	# coords=full_db[358,8:9]
	# coords=c(lat[6],lon[6])
	# end testing
	# require(soilDB)
	
	bbox<-c(coords[2]-0.02, coords[1]-0.02, coords[2]+0.02, coords[1]+0.02)
	if(is.null(maxdepth)){
		prof_id<- fetchKSSL(bbox=bbox)
		if(!is.null(prof_id)){
			# first pedon found is more closer so
			prof_id<-prof_id@site$pedon_id[1]
			soil.profile <- fetchKSSL(pedon_id=prof_id)
			# assumming deeper measurements are most likely to be missing
			depths<-soil.profile@horizons$hzn_bot-soil.profile@horizons$hzn_top
			weigths<-depths/sum(depths)
			sand<-weighted.mean(soil.profile@horizons$sand,weigths,na.rm = TRUE)
			clay<-weighted.mean(soil.profile@horizons$clay,weigths,na.rm = TRUE)
			OM<-weighted.mean(soil.profile@horizons$estimated_om,weigths,na.rm = TRUE)
			gravel<-weighted.mean(soil.profile@horizons$frags,weigths,na.rm = TRUE)
			bd<-weighted.mean(soil.profile@horizons$db_od,weigths,na.rm = TRUE)
			soil.dat<-c(sand,clay,OM,gravel,bd,sum(depths)/100)
			names(soil.dat)<-c("sand","clay","OM","gravel","bulk_dens","maxdepth")
		}else{
			soil.dat<-rep(NA,6)
			names(soil.dat)<-c("sand","clay","OM","gravel","bulk_dens","maxdepth")
		}	
	}else{
		prof_id<- fetchKSSL(bbox=bbox)
		if(!is.null(prof_id)){
			# first pedon found is more closer so
			prof_id<-prof_id@site$pedon_id[1]
			soil.profile <- fetchKSSL(pedon_id=prof_id)
			# assumming deeper measurements are most likely to be missing
			depths<-soil.profile@horizons$hzn_bot[soil.profile@horizons$hzn_bot<maxdepth]-soil.profile@horizons$hzn_top[soil.profile@horizons$hzn_bot<maxdepth]
			weigths<-depths/sum(depths)
			
			sand<-weighted.mean(soil.profile@horizons$sand[soil.profile@horizons$hzn_bot<maxdepth],weigths,na.rm = TRUE)
			clay<-weighted.mean(soil.profile@horizons$clay[soil.profile@horizons$hzn_bot<maxdepth],weigths,na.rm = TRUE)
			OM<-weighted.mean(soil.profile@horizons$estimated_om[ soil.profile@horizons$hzn_bot<maxdepth],weigths,na.rm = TRUE)
			gravel<-weighted.mean(soil.profile@horizons$frags[ soil.profile@horizons$hzn_bot<maxdepth],weigths,na.rm = TRUE)
			bd<-weighted.mean(soil.profile@horizons$db_od[ soil.profile@horizons$hzn_bot<maxdepth],weigths,na.rm = TRUE)
			soil.dat<-c(sand,clay,OM,gravel,bd,sum(depths)/100)
			names(soil.dat)<-c("sand","clay","OM","gravel","bulk_dens","maxdepth")
			
		}else{
			soil.dat<-rep(NA,6)
			names(soil.dat)<-c("sand","clay","OM","gravel","bulk_dens","maxdepth")
			
		}
		
	}
	
	
	return(soil.dat)
}

