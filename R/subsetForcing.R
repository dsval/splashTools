#' subsetForcing
#'
#' Subset spatially and temporally low resolution data (currently only WFDEI) and do the downscaling
#' @param  forcing_path,varnam,dem,dem_lowres,yr.start,yr.end,ouputdir
#' @import rgdal 
#' @import raster  
#' @importFrom MASS rlm
#' @keywords splash
#' @export
#' @examples
#' splash.grid()

subsetForcing<-function(forcing_path,varnam,dem,dem_lowres,yr.start,yr.end,ouputdir=getwd(), downscale=TRUE){
	
	# require(raster)
	# require(rgdal)
	# testing
	# forcing_path = "C:/Water_Data/WFDEI"
	# varnam = "Tair"
	# bbox<-readOGR("C:/Water_Data/iMHEA/iMHEA_geodata/JTU/JTU_04", "iMHEA_JTU_04_catchment")
	# dem<-raster("C:/Water_Data/iMHEA/iMHEA_geodata/JTU/JTU_04/dem_jtu_04")
	# end testing
	
	
	elevwfdei<-dem_lowres	
	setwd(dirname(rasterTmpFile()))
	###############################################################################################
	# 1. get the area from the dem to extract the forcing data
	###############################################################################################
	
	wgs<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
	ncores<-parallel::detectCores()-1
	if(as.character(crs(dem))!=wgs){
		cat("reprojecting dem to WGS84")
		beginCluster(ncores, type='SOCK')
		dem<-projectRaster(dem,crs=wgs,filename="dem.grd",overwrite=TRUE)
		endCluster()
		gc()
	}
	areatotal<-sum(area(dem)[!is.na(dem)])
	bbox<-extent(dem)
	# bbox<-rasterToPolygons(bbox,na.rm=TRUE,fun=function(x){x==1},dissolve=TRUE)
	###############################################################################################
	# 2. read the forcing data
	###############################################################################################
	
	filenames<- list.files(path=forcing_path,pattern = paste0(varnam,".*.nc$"), full.names=TRUE)
	lowres<-brick(filenames[1])
		
	if(length(filenames)>1){
		for(i in 2: length(filenames)){lowres<-addLayer(lowres,brick(filenames[i]))}
	}
	# 2.1 set the time and subset requested time
	ind<-seq(as.Date('1979-01-01'),as.Date('2016-12-31'), by="day")
	# ind<-seq(as.Date('2010-01-01'),as.Date('2010-01-31'), by="day")
	lowres<-setZ(lowres,ind)
	years<-as.numeric(format(ind,format="%Y"))
	lowres<-subset(lowres, which(getZ(lowres) > as.Date(paste(yr.start-1,12,31,sep="-"),format="%Y-%m-%d") & getZ(rad_lowres) <as.Date(paste(yr.end+1,01,01,sep="-"),format="%Y-%m-%d") ))
	ind<-getZ(lowres)
	
	if(areatotal<2500 && varnam == "SWdown"){
		lowres<-raster::extract(lowres,bbox)
			if(length(lowres[,1])>1){
				lowres<-colMeans(lowres[[1]])
			}
			centroid<-c((bbox@xmax+bbox@xmin),(bbox@ymax+bbox@ymin))/2
			resx<-(bbox@xmax-bbox@xmin)
			resy<-(bbox@ymax-bbox@ymin)
			require(raster)
			df<-data.frame(X=centroid[1],Y=centroid[2],Z=t(lowres))
			lowres<-rasterFromXYZ(xyz = df,res=c(resx,resy),crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
			lowres<-disaggregate(lowres,8)
			lowres<-setZ(lowres,ind)
			hres<-downscaleSolar(dem,lowres)
			hres<-writeRaster(rad_in,paste0(ouputdir,"sw_in",".nc"),format="CDF",overwrite=TRUE,varname="sw_in", varunit="W/m2", longname="daily solar radiation", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(years[1]-1,"-",12,"-",31)))
	}else if(areatotal<2500 && varnam == "Tair"){
		centroid<-c((bbox@xmax+bbox@xmin),(bbox@ymax+bbox@ymin))/2
		pts_train<-SpatialPoints(data.frame(x=seq(-359.5,359.5,1),y=rep(centroid[2],720)),proj4string=CRS(wgs))
		elev_train<-raster::extract(elevwfdei,pts_train)
		temp_train<-raster::extract(lowres,pts_train)-273.15
				
		hres<-list()
		for (i in 1:length(temp_train[1,])){
			grd<-rlm(as.numeric(temp_train[,i])~elev_train)
			hres[[i]]<-grd$coefficients[1]+grd1$coefficients[2]*dem
			
		}
		gc()
		hres<-brick(hres)
		
	}else if(areatotal<2500 && varnam == "Rainf"){
		lowres<-raster::extract(lowres,bbox)
		if(length(lowres[[1]][,1])>1){
			lowres<-colMeans(lowres[[1]])
		}
		centroid<-c((bbox@xmax+bbox@xmin),(bbox@ymax+bbox@ymin))/2
		resx<-(bbox@xmax-bbox@xmin)
		resy<-(bbox@ymax-bbox@ymin)
		require(raster)
		df<-data.frame(X=centroid[1],Y=centroid[2],Z=t(lowres))
		lowres<-rasterFromXYZ(xyz = df,res=c(resx,resy),crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
		lowres<-disaggregate(lowres,8)
		lowres<-setZ(lowres,ind)
		hres<-downscaleRLM(lowres,dem)$dowscaled
		
	}
	
	return(hres)	
	
	
	
	}
