#' clipGrid
#'
#' clips and mask a raster object, based on a polygon
#' @param   raster, Raster* object
#' @param   poly, a SpatialPolygons* object
#' @return  A Raster* object
#' @import raster 
#' @keywords splash
#' @export
#' @examples
#' clipGrid()
clipGrid<-function(raster,poly) {
	clip1 <- crop(raster, extent(poly))
	clip2 <- mask(clip1,poly)
	return(clip2)
	}	
