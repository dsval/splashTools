#' krigeForcing
#'
#' Wrapper of the hydroTSM::hydrokrige to apply automatic kriging over time series in parallel*
#' @param  stations.sp: SpatialPointsDataFrame object, mandatory attributes : lat, lon, code of the station, min 4 points
#' @param  data.df: xts matrix, names of the columns should be the same as in stations.sp$code
#' @param  dem: Raster layer, elevation (masl)
#' @return  A RasterBrick object with z time dimension saved to outdir as a netcdf file
#' @import raster  
#' @importFrom xts xts
#' @importFrom hydroTSM hydrokrige
#' @keywords splash
#' @export
#' @examples
#' *optional run beginCluster() first, for parallel computing
#' krigeForcing()

krigeForcing<-function(stations.sp,data.df,dem,outdir=getwd()){
	###############################################################################################
	# Calling libraries
	###############################################################################################
	# require(raster)
	# require(hydroTSM)
	# require(xts)
	# testing
	# dem<- raster("C:/Water_Data/iMHEA/iMHEA_geodata/HMT/elev.2015_2016.HMT.nc")
	# stations.sp<-readOGR(dsn="C:/Water_Data/iMHEA/iMHEA_geodata/HMT", layer= "stations",pointDropZ=TRUE)
	# stations.sp$code<-stations.sp@data$Código_iM
	# data.df<-prececu2012
	# end testing
	###########################################################################
	# 00. Check if parallel computation is required by the user and if the dimensions of the raster objects match
	###########################################################################
	on.exit(endCluster())
	clcheck<-try(getCluster(), silent=TRUE)
	if(class(clcheck)=="try-error"){
		# If no cluster is initialized, assume only one core will do the calculations, beginCluster(1) saved me the time of coding serial versions of the functions
		beginCluster(1,'SOCK')
		message('Only using one core, use first beginCluster() if you want to run splash in parallel!!')
		
	}
	y<-as.numeric(unique(format(time(data.df),'%Y')))
	
	wgs<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
	###############################################################################################
	# 00. Setting up the inputs
	###############################################################################################
	###############################################################################################
	# 00. Projecting to UTM Hydrokrigde requires this!!
	###############################################################################################
		
	
	find_zone<- function(longitude, latitude) {
		
		# Special zones for Svalbard and Norway
		if (latitude >= 72.0 && latitude < 84.0 ){
			if (longitude >= 0.0  && longitude <  9.0){
				zone=31	
			} else if (longitude >= 9.0  && longitude < 21.0){
				zone=33	
			}else if (longitude >= 21.0 && longitude < 33.0){
				zone=35
			}else if (longitude >= 33.0 && longitude < 42.0) {
				zone=37	
			}
		} else{
			zone=(floor((longitude + 180) / 6) %% 60) + 1
		}
				
		hemisphere = ifelse(latitude > 0, "north", "south")
		return(list(zone=zone,hemisphere = hemisphere))
	}
	# get utm projection info
	if(as.character(crs(stations.sp))!=wgs){
		centroid<-c((stations.sp@bbox[1,1]+stations.sp@bbox[1,2]),(stations.sp@bbox[2,1]+stations.sp@bbox[2,2]))/2
		utm.info<-find_zone(centroid[1],centroid[2])
		utm_proj<-CRS(paste0("+proj=utm +zone=",utm.info$zone," +",utm.info$hemisphere," +datum=WGS84 +units=m +no_defs+ellps=WGS84 +towgs84=0,0,0"))
		stations.sp<-spTransform(stations.sp,utm_proj)
	}
	if(as.character(crs(dem))!=as.character(crs(stations.sp))){
		cat("reprojecting dem to UTM")
		dem<-projectRaster(dem,crs=crs(stations.sp),filename="dem.grd",overwrite=TRUE)
	}
	
	###############################################################################################
	# 02. assigning elevation to coordinates
	###############################################################################################
	stations.sp$elev<-raster::extract(dem,stations.sp)
	###############################################################################################
	# 03. Building dataframes with spatialinfo
	###############################################################################################
	xgis<-data.frame(stations.sp$code,stations.sp@coords[,2],stations.sp@coords[,1],stations.sp$elev)
	names(xgis)<-c("ID","POINT_Y","POINT_X", "Z")
	
	###############################################################################################
	# 04. Building dataframes with the timeseries 
	###############################################################################################
		
	xts.dfs<-vector("list", length(data.df[,1])) 
	
	for(i in 1:length(data.df[,1])){
		xts.dfs[[i]]<-as.numeric(data.df[i,])
		names(xts.dfs[[i]])<-colnames(data.df)
	}
	
	###############################################################################################
	# 05. Building the spatial dataframes
	###############################################################################################
	
	dem.sgdf<-as(dem, 'SpatialGridDataFrame')
	names(dem.sgdf@data)<-"Z"
	###############################################################################################
	# 05. Interpolating all timesteps in parallel
	###############################################################################################
	
	p4s <- crs(stations.sp)
	build.lay<-function(x,y){
		# x= datafrom xts.dfs
		# y= interpolated layer
		if (sum(x)==0){
			y <-0	
		}
		else{
			y<-hydrokrige(x.ts= x, x.gis=xgis,
				X="POINT_X", Y="POINT_Y", sname="ID", elevation="Z",
				type= "cells",formula=value~Z,
				p4s=p4s,predictors=dem.sgdf, 
				plot=FALSE)
			y
			
			
		}
	}
	# create empty array
	x.ked<-vector("list", length(data.df[,1]))
		
	cl<-getCluster()

	parallel::clusterEvalQ(cl, library("hydroTSM"))
	parallel::clusterExport(cl, c("xts.dfs","x.ked","xgis","dem.sgdf","p4s","build.lay"),envir=environment()) 
	cat("Kriging",length(data.df[,1]),'layers',"\n")
	x.ked<-parallel::clusterMap(cl = cl, fun=build.lay,xts.dfs,x.ked)	
	# identify zero values
	indtemp<-which(sapply(x.ked, FUN=function(X) class(X)=="SpatialGridDataFrame"))
	# create a spatial template with zero values
	zerotemp<-x.ked[[indtemp[1]]]
	zerotemp@data<-zerotemp@data*0
	# replace numeric zeros with the spatial zero template
	replacezero<-function(X){
		if( class(X)!="SpatialGridDataFrame"){
			X<-zerotemp
		}else{
			X
		}
	}
	x.ked<-lapply(x.ked,replacezero)
	
	gc()
	# build the rasters
	x.ked<-lapply(x.ked,raster)
	x.ked<-stack(x.ked)
	x.ked[x.ked<0.0]<-0.0
	# reproject to wgs and write

	x.ked<-projectRaster(setZ(x.ked,as.Date(time(data.df))),crs=wgs,filename=paste0(outdir,"/",y[1],"_",y[length(y)],".","pn",".","nc"),format="CDF",overwrite=TRUE,varname="pn", varunit="mm", longname="Precipitation", xname="lon", yname="lat", zname="time")


	return(x.ked)
	
	
}
