#' cropmem
#'
#' crops and mask a raster object, based on a polygon
#' @param   raster, Raster* object
#' @param   poly, a SpatialPolygons* object
#' @return  A Raster* object
#' @import raster 
#' @keywords splash
#' @export
#' @examples
#' cropmem()
cropmem<-function(r,p){
	cells<-cellFromPolygon(r,p)
	nls<-nlayers(r)
	ext<-extent(p)
	colstrt<-colFromX(r,ext[1])
	rowstrt<-rowFromY(r, ext[4])
	
	if(nls>1){
		newr<-brick(rasterFromCells(r, cells[[1]], values=F),values=F, nl=nls)
	}else{
		newr<-rasterFromCells(r, cells[[1]], values=F)
	}
	vals<-getValuesBlock(r,row=rowstrt,nrows=nrow(newr),col=colstrt,ncols=ncol(newr))
	newr<-setValues(newr,vals)
	return(newr)
}


calcStrmflow<-function(ro,bflow){
	### 1. get time info
	ztime<-getZ(ro)
	time.freq<-ztime[2]-ztime[1]
	#calc factor to convert time to seconds
	tfactor<-as.numeric(time.freq, units = "days")*86400
	## calc areas
	area_km2<-area(ro[[1]])
	area_m2<-area_km2*1e6
	area_km2=mask(area_km2,ro[[1]])
	area_km2<-cellStats(area_km2,sum)
	##2. calc stormflow litre/s/km2
	ro_ltrs<-ro*area_m2
	ro_ltrs<-cellStats(ro_ltrs,sum)
	# ro_ltrs<-ro_ltrs/(area_km2*tfactor)
	# ro_ltrs<-xts(ro_ltrs,ztime)
	# m3/s
	ro_m3<-(ro_ltrs/1000)/tfactor
	## 3. calc baseflow litre/s/km2
	bflow_ltrs<-bflow*area_m2
	bflow_ltrs<-cellStats(bflow_ltrs,sum)
	# bflow_ltrs<-bflow_ltrs/(area_km2*tfactor)
	# bflow_ltrs<-xts(bflow_ltrs,ztime)
	bflow_m3<-(bflow_ltrs/1000)/tfactor
	### 4. streamflow
	#strm<-ro_ltrs+bflow_ltrs
	strm<-ro_m3+bflow_m3
	result<-data.frame(Date=ztime,RO=ro_m3,BF=bflow_m3,Q=strm)
	#result<-merge(ro_ltrs,bflow_ltrs,strm)
	#names(result)<-c('RO.l/s/km2','BFLW.l/s/km2','SFLW.l/s/km2')
	#names(result)<-c('RO.l/s/km2','BFLW.l/s/km2','SFLW.l/s/km2')
	result
}



