#' downscaleKrige
#'
#' Downscale the input applying the autokrige from hydroTSM::hydrokrige, using elevation as predictor
#' @param  elev_hres, elevation, a high resolution Raster layer object
#' @param  rad_lowres, solar radiation (w/m2), the low resolution Raster* object with z time dimension
#' @param  ouputdir, directory where to save the files
#' @return  solar radiation (w/m2), Raster* object with the resolution of elev_hres and the same z time dimension as rad_lowres
#' @import raster  
#' @keywords Downscale, solar radiation, terrain
#' @export
#' @examples
#' # *optional run beginCluster() first, for parallel computing
#' downscaleKrige()
downscaleKrige<-function(low_res_raster,high_res_raster,ouputdir=getwd()){
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
	###########################################################################
	# 01. get the timeseries df
	###########################################################################
	timeind<-getZ(low_res_raster)
	pts<-rasterToPoints(low_res_raster,spatial=TRUE)
	pts$code<-paste0('code_',1:length(pts[,1]))
	data.df<-xts(t(data.df@data),timeind)
	names(data.df)<-pts$code
	
	downsc<-krigeForcing(pts,data.df,dem=high_res_raster,outdir = outdir)
	
	return(downsc)
}
