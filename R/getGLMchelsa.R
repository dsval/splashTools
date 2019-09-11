#' getGLMchelsa
#'
#' Function to ease the download of the high resolution climate projections from the chelsa server
#' @param   model: character describing the glm, check http://chelsa-climate.org/downloads/
#' @param   scenario: character: 'rcp45', 'rcp60' or 'rcp85'
#' @param   yr: integer year of the projection, only available: 2060 or 2080
#' @import raster 
#' @import rvest 
#' @keywords splash
#' @return list with monthly temperature and precipitation from the GLM's projections, both saved as netCDF files
#' @export
#' @examples getGLMchelsa(model='bcc',scenario='rcp85',yr=2060)

getGLMchelsa<-function(model,scenario,yr,outdir=getwd()){
	if(yr==2060){
		yr_dir<-'/2041-2060/'
		yr_info<-'2041-2060'
	}else if(yr==2080){
		yr_dir<-'/2061-2080/'
		yr_info<-'2061-2080'
	}
	url_prec<-paste0('https://www.wsl.ch/lud/chelsa/data/cmip5',yr_dir,'prec/')
	url_temp<-paste0('https://www.wsl.ch/lud/chelsa/data/cmip5',yr_dir,'temp/')
	# get the links for precipitation	
	doc_prec <- read_html(url_prec)
	all_prec <- html_attr(html_nodes(doc_prec, "a"), "href")
	# model <- "bcc"
	files_prec <- all_prec[grep(model, all_prec)]
	files_prec <- files_prec[grep(scenario, files_prec)]
	
	# get the links for temperature
	doc_temp <- read_html(url_temp)
	all_temp <- html_attr(html_nodes(doc_temp, "a"), "href")
	
	files_temp <- all_temp[grep(model, all_temp)]
	files_temp <- files_temp[grep(scenario, files_temp)]
	
	setwd(dirname(rasterTmpFile()))
	# get the datasets
	options(download.file.extra = c("--no-check-certificate"))
	download.file(paste0(url_prec,files_prec[1]),files_prec[1],method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
	dataset_prec<-brick(files_prec[1])
	dataset_temp<-brick(files_temp[1])
	for(i in 2:length(files_prec)){
		
		download.file(paste0(url_prec,files_prec[i]),files_prec[i],method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
		download.file(paste0(url_temp,files_temp[i]),files_temp[i],method = "auto",quiet = TRUE,cacheOK = TRUE,mode = "wb")
		dataset_prec<-addLayer(dataset_prec,raster(files_prec[i]))
		dataset_temp<-addLayer(dataset_temp,raster(files_temp[i]))		
	}
	gc()
	temp_filename<-do.call(rbind,strsplit(files_temp,'.nc',fixed = TRUE))[1,1]
	prec_filename<-do.call(rbind,strsplit(files_prec,'.nc',fixed = TRUE))[1,1]
	
	
	
	prec<-writeRaster(dataset_prec,filename=paste0(outdir,'/',prec_filename,'_',yr_info,'.nc'), format="CDF", overwrite=TRUE,varname="prec", varunit="mm",longname='projection precipitation GLM', xname="lon", yname="lat")
	temp<-writeRaster(dataset_temp,filename=paste0(outdir,'/',temp_filename,'_',yr_info,'.nc'), format="CDF", overwrite=TRUE,varname="temp", varunit="C/10",longname='projection precipitation GLM', xname="lon", yname="lat")
	
	return(list(prec,temp))
}

