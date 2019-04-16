#' clipPointGrid
#'
#' function to clip points(sp) in overlapped ares where raster == value
#' @param   raster,sp,value
#' @import raster 
#' @import rgdal 
#' @keywords splash
#' @export
#' @examples
#' splash.grid()

clipPointGrid<-function(raster,sp,value){
	# function to clip points in overlapped ares where raster == value
	# require(raster)
	# require(rgdal)
	subs <- which(extract(raster == value,sp) == 1)
	clipped <- sp[subs,]
	return(clipped)
}
