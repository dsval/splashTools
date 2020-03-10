#' downscaleRLM
#'
#' Downscale precipitation or temperature using robust linear models, wrapper from the spatialEco package
#' @param  fields on stations.sp: lat, lon, code
#' @importFrom MASS rlm
#' @import raster  
#' @keywords splashTools
#' @export
#' @examples
#' splash.grid()

downscaleRLM<-function(low_res_raster,high_res_raster){
	# ************************************************************************
	# Name:     forcing.splash
	# Inputs:   - raster object, low_res_raster
	#           - raster object, high_res_raster
	#           - region to subset, character	
	# Returns:  raster object clepped and resampled
	# Features: Downscale a forcing dataset using robust linear modeling rlm.
	# Depends:  - ke ............. eccentricity of earth's orbit, unitless
	#           - komega ......... longitude of perihelion
	#  Ref:     Berger, A. L. (1978), Long term variations of daily insolation
	#             and quaternary climatic changes, J. Atmos. Sci., 35, 2362-2367.
	# ************************************************************************
	# rasterOptions(maxmemory=3e7, timer=TRUE, progress = "text", tmptime = 24, chunksize = 1e7, overwrite=TRUE,tolerance=0.5,todisk=TRUE)
	raster.downscale <- function(x, y, p = NULL, n = NULL, filename = FALSE,scatter = FALSE, ...) {
		# ************************************************************************
		# Name:     raster.downscale 
		# Inputs:   - raster object, low_res_raster
		#           - raster object, high_res_raster
		#           - region to subset, character	
		# Returns:  raster odownscale
		# Features: Downscaling corrected function form spatialEco progress bar bug
		
		# ************************************************************************	
		# testing
		# x=high_res_raster; y=low_res_raster[[1]];scatter=FALSE
		if(!class(y) == "RasterLayer") stop( "y is not a raster object")
		if(!class(x) == "RasterLayer" & !class(x) == "RasterStack" & !class(x) == "RasterBrick")
		stop( "x is not a raster object")  		
		x <- raster::stack(x)
		if(is.null(p) & is.null(n)) {
			warning("Population is being used and may cause memory issues")
			sub.samp <- raster::rasterToPoints(y, spatial=TRUE)
		} else {  
			if(!is.null(n)) { sampSize = n }
			if(!is.null(p)) { sampSize = round((raster::ncell(y)*p),0) }
			sub.samp <- raster::sampleRandom(y, sampSize, sp=TRUE)
		}	  
		sub.samp@data <- data.frame(sub.samp@data, raster::extract(x, sub.samp) )
		names(sub.samp@data) <- c("y", names(x))
		sub.samp <- stats::na.omit(sub.samp@data)	  
		rrr <- rlm(stats::as.formula(paste(names(sub.samp)[1], ".", sep=" ~ ")), 
			data=sub.samp, scale.est="Huber", psi=MASS::psi.hampel, init="lts",maxit = 100)
		if(scatter == TRUE) { graphics::plot(sub.samp[,2], sub.samp[,1], pch=20, cex=0.50,
			xlab=names(sub.samp)[2], ylab="y") }				
	if (filename != FALSE) {
		raster::predict(x, rrr, filename=filename, na.rm=TRUE, progress='', 
			overwrite=TRUE, ...)
		print(paste("Raster written to", filename, sep=": "))	
		return(list(model = rrr, MRE = round(mean(rrr$residuals), digits=4), 
			AIC = round(AIC(rrr), digits=4)))			  
} else {
	r <- raster::predict(x, rrr, na.rm=TRUE, progress='')
	return(list(downscale = r, model = rrr, MSE = round(mean(rrr$residuals), digits=4), 
		AIC = round(AIC(rrr), digits=4)))		
}
}
# ************************************************************************
spatial.downscale<-function(low_res_raster,high_res_raster,filename=NULL){
	
	# testing
	# low_res_raster=lowres;high_res_raster=dem
	# compareRaster(low_res_raster,high_res_raster,extent=TRUE, rowcol=FALSE, crs=TRUE, res=FALSE, orig=TRUE,
	# 	rotation=TRUE, values=FALSE, tolerance=0.5, stopiffalse=FALSE, showwarning=TRUE)
	
	if (nlayers(high_res_raster)>1){high_res_raster<-high_res_raster[[1]]}
	if(nlayers(low_res_raster)==1){dowscaled<-raster.downscale(high_res_raster, low_res_raster, scatter=FALSE)
		return(dowscaled)}
	else{results<-list()
		
		dowscaledtemp<-raster.downscale(x=high_res_raster, y=low_res_raster[[1]],scatter=FALSE)
		results$dowscaled<-brick(dowscaledtemp$downscale)
		results$models<-list(dowscaledtemp$model)
		results$MSEs<-dowscaledtemp$MSE
		results$AICs<-dowscaledtemp$AIC
		for(i in 2:nlayers(low_res_raster)){
			dowscaledtemp<-raster.downscale(x=high_res_raster, y=low_res_raster[[i]], scatter=FALSE)
			results$dowscaled<-addLayer(results$dowscaled,dowscaledtemp$downscale)
			results$models[[i]]<-dowscaledtemp$model
			results$MSEs<-c(results$MSEs,dowscaledtemp$MSE)
			results$AICs<-c(results$AICs,dowscaledtemp$AIC)	
		}
		
		return(results)	
	}	
	
	
}


downscaled<-spatial.downscale(low_res_raster,high_res_raster)
return(downscaled)
}
