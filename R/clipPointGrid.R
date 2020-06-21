#' clipPointGrid
#'
#' function to clip points(sp) in overlapped ares where raster == value
#' @param   raster, Raster layer
#' @param   sp, SpatialPoints or SpatialPointsDataFrame object
#' @param   value, numeric
#' @return  A subset of sp
#' @import raster 
#' @import rgdal 
#' @keywords splash
#' @export
#' @examples
#'clipPointGrid()

clipPointGrid<-function(raster,sp,value){
	# function to clip points in overlapped ares where raster == value
	# require(raster)
	# require(rgdal)
	subs <- which(extract(raster == value,sp) == 1)
	clipped <- sp[subs,]
	return(clipped)
}
