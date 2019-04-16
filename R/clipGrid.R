#' clipGrid
#'
#' clips and mask a raster object, based on a polgon
#' @param   raster,poly
#' @import raster 
#' @keywords splash
#' @export
#' @examples
#' splash.grid()
clipGrid<-function(raster,poly) {
	clip1 <- crop(raster, extent(poly))
	clip2 <- mask(clip1,poly)
	return(clip2)
	}	
