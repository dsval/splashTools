#' getModisClim
#'
#' download and gapfill MOD07 and MYD07 atmospheric profiles, calculate saturated vapour pressure using dowsncaled LST from Microwave SSMI and Modis IR, if used with the option use.clouds=TRUE, additional files from MOD06 and MYD06 cloud product will be used to infer temperature and actual vapour pressure under the clouds below the tropopause. 
#' @param   coords vector c(lat,lon)
#' @param   start, end : data range
#' @import raster
#' @import gdalUtils
#' @import rgdal
#' @import httr
#' @import xml2
#' @keywords modis
#' @export
#' @examples
#' getModisClim()

getModisClim<-function(lat,lon,start,end,outmode=list(tile=TRUE,monthly=TRUE,use.clouds=FALSE),dem,outdir=getwd(), tmpdir=dirname(rasterTmpFile()),usr='usr',pass='pass'){
	# testing
	# on.exit(traceback(2))
	# rasterOptions(todisk=FALSE)
	# end testing
########################################################################
#1.get the urls
########################################################################
	#build the query
	cat('Retrieving the urls',"\n")
 url<- "http://modwebsrv.modaps.eosdis.nasa.gov/axis2/services/MODAPSservices/"
 if (outmode$use.clouds==TRUE){
	query_par <- list(products = gsub(" ", "", toString(c('MOD07_L2','MYD07_L2','MOD06_L2','MYD06_L2'))),
		collection = '61',
		startTime = start,
		endTime = end,
		north = lat+0.1,
		south = lat-0.1,
		east = lon+0.1,
		west = lon-0.1,
		coordsOrTiles = 'coords',
		dayNightBoth = 'DB')
 }else{
	query_par <- list(products = gsub(" ", "", toString(c('MOD07_L2','MYD07_L2'))),
		collection = 61,
		startTime = start,
		endTime = end,
		north = lat+0.1,
		south = lat-0.11,
		east = lon+0.1,
		west = lon-0.1,
		coordsOrTiles = 'coords',
		dayNightBoth = 'DB')
 }

	query_out<-httr::GET(url = paste0(url, "searchForFiles"),query = query_par)
	# get the fileids
	fileIds <- httr::content(query_out, as = "text")
	# fileIds <- xml2::read_xml(fileIds)
	fileIds<- xml2::read_html(fileIds)
	fileIds <- do.call(rbind, xml2::as_list(fileIds)[[1]][[1]][[1]])
	# get the urls
	get_urls<-function(fileid){
		files_md<-httr::GET(url = paste0(url,"getFileUrls"),query = list(fileIds=fileid))
		fileurl <- httr::content(files_md, as = "text")
		# fileurl <- xml2::read_xml(fileurl)
		fileurl <- xml2::read_html(fileurl)
		fileurl <- do.call(rbind, xml2::as_list(fileurl)[[1]][[1]][[1]])
		fileurl
	}
	
	file_urls<-mapply(FUN=get_urls,fileIds)
	file_urls<-do.call(c,file_urls)
	zdates<-do.call(rbind,strsplit(file_urls,'.',fixed=T))
	# length(zdates[1,])-4
	zdates<-strptime(paste0(zdates[,length(zdates[1,])-4],zdates[,length(zdates[1,])-3]),format='A%Y%j%H%M')
	file_urls<-cbind(file_urls,as.character(zdates))
	file_urls<- file_urls[order(file_urls[,2]),]
	########################################################################
	#1.get monthly MODIS LST urls
	########################################################################
	qextent<-query_par
	qextent[[1]]<-'MOD11B3'
	qextent[[2]]<-6
	qextent[[10]]<-'D'
	format(as.Date(start),'%m')
	format(as.Date(start),'%Y')
	qextent[[3]]<-paste0(format(as.Date(start),'%Y'),'-',format(as.Date(start),'%m'),'-01')
	qextent[[4]]<-paste0(format(as.Date(end),'%Y'),'-',format(as.Date(end),'%m'),'-01')
	query_oute<-httr::GET(url = paste0(url, "searchForFiles"),query = qextent)
	fileIdse <- httr::content(query_oute, as = "text")
	fileIdse <- xml2::read_html(fileIdse)
	fileIdse <- do.call(rbind, xml2::as_list(fileIdse)[[1]][[1]][[1]])
	# urlext<-get_urls(fileIdse[,1])
	urlext<-mapply(FUN=get_urls,fileIdse)
	urlext<-do.call(c,urlext)
	hv<-do.call(rbind,strsplit(urlext,'.',fixed = T))
	hv<-hv[,length(hv[1,])-3]
	beg<-format(as.Date(start),'%Y%j')
	til<-format(as.Date(end),'%Y%j')
	urlext<-do.call(c,regmatches(urlext, gregexpr(paste0('.*.',hv[1],'.*.'),urlext)))
	# file_urls<-file_urls[!duplicated(file_urls[,2]),]
	########################################################################
	#1.get the urls for lst microwave SSM/I-SSMIS Pathfinder
	########################################################################
	dateseq<-format(seq(as.Date(start),as.Date(end),by='day'),format='%Y.%m.%d')
	urlsSSM<-paste0('https://n5eil01u.ecs.nsidc.org/PM/NSIDC-0032.002/',as.character(dateseq),'/')
	
	get_ssmurl<-function(dayurl){
		url_file_SSM<-GET(dayurl,authenticate(usr, pass))
		url_file_SSM <- httr::content(url_file_SSM, as = "parsed")
		url_file_SSM <-xml2::as_list(url_file_SSM)
		url_file_SSM<-do.call(rbind,url_file_SSM$html$body[[4]]$table)
		# get brightness temperature 37gz V polarization Holmes, et al. 2008 doi 10.1029/2008JD010257
		url_file_SSM<-regmatches(unlist(url_file_SSM[,2]), gregexpr('.*.ML.*.A.*.37V.*.gz$', unlist(url_file_SSM[,2])))
		url_file_SSM<-do.call(c,url_file_SSM)
		paste0(dayurl,url_file_SSM)
	}
	if(outmode$monthly==FALSE){
		SSM_url<-mapply(get_ssmurl,urlsSSM,SIMPLIFY = T)	
	}
	
	
	########################################################################
	#set credentials copied from https://git.earthdata.nasa.gov/projects/LPDUR/repos/daac_data_download_r/browse/DAACDataDownload.R
# 	########################################################################
# 	# ---------------------------------SET UP ENVIRONMENT--------------------------------------------- #
# 	# IMPORTANT: Update the line below if you want to download to a different directory (ex: "c:/data/")
# 	#dl_dir <- Sys.getenv("HOME")                                 # Set dir to download files to
# 	#setwd(dl_dir)                                                # Set the working dir to the dl_dir
# 	usr <- file.path(Sys.getenv("USERPROFILE"))                  # Retrieve home dir (for netrc file)
# 	if (usr == "") {usr = Sys.getenv("HOME")} 
# 	# usr<-gsub("\\","/",usr,fixed=TRUE)			# If no user profile exists, use home
# 	netrc <- file.path(usr,'.netrc', fsep = .Platform$file.sep)  # Path to netrc file
# 	
# 	# ------------------------------------CREATE .NETRC FILE------------------------------------------ #
# 	# If you already have a .netrc file with your Earthdata Login credentials stored in your home
# 	# directory, this portion will be skipped. Otherwise you will be prompted for your NASA Earthdata
# 	# Login Username/Password and a netrc file will be created to store your credentials (in home dir)
# 	if (file.exists(netrc) == FALSE || grepl("ladsweb.modaps.eosdis.nasa.gov", readLines(netrc)) == FALSE) {
# 		netrc_conn <- file(netrc)
# 		
# 		# User will be prompted for NASA Earthdata Login Username and Password below
# 		writeLines(c("machine ladsweb.modaps.eosdis.nasa.gov",
# 			sprintf("login %s", getPass(msg = "Enter NASA Earthdata Login Username \n (or create an account at urs.earthdata.nasa.gov) :")),
# 			sprintf("password %s", getPass(msg = "Enter NASA Earthdata Login Password:"))), netrc_conn)
# 	close(netrc_conn)
# }

########################################################################
#2.download the files
########################################################################
	# tmpdir<-dirname(rasterTmpFile())
	setwd(tmpdir)
	# setwd('C:/Rcalculations')
	destfiles<-basename(file_urls[,1])
	# pb <- pbCreate(length(file_urls[,1]))
	pb <- txtProgressBar(min=1,max = length(file_urls[,1]), style = 3)
	# Loop through all files
	########################################################################
	#2.download athmosphere
	########################################################################
	
	cat('downloading atmospheric profiles',"\n")
	for (i in 1:length(file_urls[,1])) {
		if(!file.exists(destfiles[i])){
			# Write file to disk (authenticating with netrc) using the current directory/filename
			response <- httr::GET(as.character(file_urls[i,1]), write_disk(destfiles[i], overwrite = TRUE), 
				authenticate(usr, pass))
			while(file.size(destfiles[i])<=200){
				response <- httr::GET(as.character(file_urls[i,1]), write_disk(destfiles[i], overwrite = TRUE), 
					authenticate(usr, pass))
			}
			
			# Check to see if file downloaded correctly
			if (response$status_code == 200) {
				cat(destfiles[i],"downloaded","\n")
			} else {
				cat("error downloading, make sure usr/password are correct","\n")
				errordwn<-response$status_code
			}
		}
	
		setTxtProgressBar(pb,i)
	}
	close(pb)
	gc()	
	# 
	########################################################################
	#2.download MODIS LST
	########################################################################
	destfileext<-basename(as.character(urlext))
		
	if(length(destfileext)>1){pb <- txtProgressBar(min=1,max = length(destfileext), style = 3)}
	
	# destfiles<-basename(file_urls[,1])
	cat(' ',"\n")
	cat('downloading MODIS LST',"\n")
	for(i in 1:length(destfileext)){
		if(!file.exists(destfileext[i])){
		respon <-httr::GET(as.character(urlext[i]), write_disk(destfileext[i], overwrite = TRUE), 
			authenticate(usr, pass))
		while(file.size(destfileext[i])<=200){
			respon <-httr::GET(as.character(urlext[i]), write_disk(destfileext[i], overwrite = TRUE), 
				authenticate(usr, pass))
		}
		# Check to see if file downloaded correctly
		if (respon$status_code == 200) {
			cat(' ',"\n")
			cat(destfileext[i],"downloaded","\n")
		} else {
			cat("error downloading, make sure usr/password are correct","\n")
		}
		if(length(destfileext)>1){setTxtProgressBar(pb,i)}
	}
	}
	########################################################################
	#2.download SSM LST
	########################################################################
	if(outmode$monthly==FALSE){
		destfiles_ssm<-basename(SSM_url)
		cat(' ',"\n")
		cat('downloading SSM/I, SSMIS LST',"\n")
		pb <- txtProgressBar(min=1,max = length(destfiles_ssm), style = 3)
		for(i in 1:length(destfiles_ssm)){
			if(!file.exists(destfiles_ssm[i])){
				respon <-httr::GET(as.character(SSM_url[i]), write_disk(destfiles_ssm[i], overwrite = TRUE), 
					authenticate(usr, pass))
				# Check to see if file downloaded correctly
				if (response$status_code == 200) {
					cat(' ',"\n")
					cat(destfiles_ssm[i],"downloaded","\n")
				} else {
					cat("error downloading, make sure usr/password are correct","\n")
				}
				if(length(destfiles_ssm)>1){setTxtProgressBar(pb,i)}
			}
		}
		filenames_SSM<-paste0(tmpdir,'/',destfiles_ssm)	
	}
	
	

########################################################################
#3.read the data, ....mapply is working faster than overay wth a rasterstack
########################################################################	
# prepare the filenames
filenames<-paste0(tmpdir,'/',destfiles)
if (outmode$use.clouds==TRUE){
	filenames_mod06<-do.call(c,regmatches(filenames, gregexpr('.*.MOD06.*.',filenames)))
	filenames_myd06<-do.call(c,regmatches(filenames, gregexpr('.*.MYD06.*.',filenames)))
	filenames_cld<-c(filenames_mod06,filenames_myd06)
	zdates_cld<-do.call(rbind,strsplit(filenames_cld,'.',fixed=T))
	# zdates_cld<-strptime(paste0(zdates_cld[,2],zdates_cld[,3]),format='A%Y%j%H%M')
	zdates_cld<-as.POSIXct(paste0(zdates_cld[,length(zdates_cld[1,])-4],zdates_cld[,length(zdates_cld[1,])-3]),format='A%Y%j%H%M', tz = "GMT")
	filenames_cld<-filenames_cld[order(zdates_cld)]
}


filenames_mod07<-do.call(c,regmatches(filenames, gregexpr('.*.MOD07.*.',filenames)))
filenames_myd07<-do.call(c,regmatches(filenames, gregexpr('.*.MYD07.*.',filenames)))
filenames_atm<-c(filenames_mod07,filenames_myd07)
zdates_atm<-do.call(rbind,strsplit(filenames_atm,'.',fixed=T))
# zdates_atm<-strptime(paste0(zdates_atm[,2],zdates_atm[,3]),format='A%Y%j%H%M')
zdates_atm<-as.POSIXct(paste0(zdates_atm[,length(zdates_atm[1,])-4],zdates_atm[,length(zdates_atm[1,])-3]),format='A%Y%j%H%M', tz = "GMT")
filenames_atm<-filenames_atm[order(zdates_atm)]
filenamlst<-paste0(tmpdir,'/',destfileext)


############################# some constants ############################
kG <- 9.80665       # gravitational acceleration, m/s^2 (Allen, 1973)
kL <- 0.0065        # adiabatic lapse rate, K/m (Cavcar, 2000)
kMa <- 0.028963     # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
kPo <- 101325       # standard atmosphere, Pa (Allen, 1973)
kR <- 8.31447       # universal gas constant, J/mol/K (Moldover et al., 1988)
kTo <- 288.15       # base temperature, K (Berberan-Santos et al., 1997)

############################# create the functions to read the data############################
# ************************************************************************
# Name:     elev2pres
# Inputs:   double (z), meters
# Returns:  double, Pa
# Features: Calculates atmospheric pressure for a given elevation
# Depends:  - kPo ............ base pressure, Pa
#           - kTo ............ base temperature, C
#           - kL ............. temperature lapse rate, K/m
#           - kR ............. universal gas constant, J/mol/K
#           - kMa ............ molecular weight of dry air, kg/mol
#           - kG ............. gravity, m/s^2
# Ref:      Allen et al. (1998)
# ************************************************************************


elev2pres <- function(z) {
	kPo*(1 - kL*z/kTo)^(kG*kMa/(kR*kL))
}
p_dryair<-function(elev,t){
	#calc dry air density [kg/m3] using pv=nRT
	(elev2pres(elev)*kMa)/(kR*(t+273.15))
}

############################# func to read modis lst############################
readlst<-function(filename){
	# filename=filenamlst[1]
	lst<-raster(get_subdatasets(filename)[1])
	lst<-projectRaster(lst,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
	lst-273.15
}
############################# func to read modis athmos############################
readMOD07<-function(filename,output='Ta'){
	# testting
	# filename=filenames[1]
	# end testing
	#read datatsets
	sds <- get_subdatasets(filename)
	md<-gdalinfo(filename)
	north<-do.call(c,regmatches(md, gregexpr(paste0('.*.','NORTHBOUNDINGCOORDINATE','.*.'),md)))
	south<-do.call(c,regmatches(md, gregexpr(paste0('.*.','SOUTHBOUNDINGCOORDINATE','.*.'),md)))
	east<-do.call(c,regmatches(md, gregexpr(paste0('.*.','EASTBOUNDINGCOORDINATE','.*.'),md)))
	west<-do.call(c,regmatches(md, gregexpr(paste0('.*.','WESTBOUNDINGCOORDINATE','.*.'),md)))
	bbox<-c(west,east,south,north)
	# bbox<-c(md[77],md[18],md[66],md[45])
	bbox<-do.call(rbind,strsplit(bbox,'='))
	bbox<-as.numeric(bbox[,2])
	wgs<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
	
	if(output=='Tsurf'){
		tsurf <- readGDAL(sds[8], as.is = TRUE)
		tsurf <-((raster(tsurf)+15000) * 0.009999999776482582)-273.15
		## set extent and coordinate reference system
		extent(tsurf)<-extent(bbox)
		projection(tsurf) <- wgs
		return(tsurf)	
	}else if(output=='Ta'){
		# mixing ratio profile(Kg/Kg)
		vw <- readGDAL(sds[17], as.is = TRUE,silent = T)
		vw<-calc(stack(vw),fun=function(x){x*(0.001000000047497451)/1000})
		# mixing ratio at surface (Kg/Kg)
		vw_surf<-max(vw,na.rm=T)
		# mask surftemp
		fmask<-function(maxvw,vw_hp){
			masked<-ifelse(vw_hp==maxvw,1,NA)
			masked
		}
		maskTa<-overlay(vw,vw_surf,fun=fmask)
		# surface Ta
		ta_press<-readGDAL(sds[15], as.is = TRUE,silent = T)
		ta_press <- calc(stack(ta_press), fun = function(x){((x+15000) * 0.009999999776482582)-273.15})
		
		# ta_press <-((stack(ta_press)+15000) * 0.009999999776482582)-273.15
		ta_press<-overlay(ta_press,maskTa,fun=function(x,y){x*y})
		ta<-mean(ta_press,na.rm=T)
		extent(ta)<-extent(bbox)
		projection(ta)<-wgs
		extent(vw_surf)<-extent(bbox)
		projection(vw_surf)<-wgs
	
		return(list(ta=ta,vw_surf=vw_surf))
	}else if(output=='a'){
		# mixing ratio profile(Kg/Kg)
		vw <- readGDAL(sds[17], as.is = TRUE)
		vw<-calc(stack(vw),fun=function(x){x*(0.001000000047497451)/1000})
		# vw <-((stack(vw)) * 0.001000000047497451)/1000
		# mixing ratio at surface (Kg/Kg)
		vw_surf<-max(vw,na.rm=T)
		extent(vw_surf)<-extent(bbox)
		projection(vw_surf)<-wgs
		return(vw_surf)
	}else if(output=='LR'){
		# reading heigth at different geopotentials 
		tropohgt <- brick(readGDAL(sds[18], as.is = TRUE,silent = T))+32500
		# subsetting readings max aronud ~8KM (lower than tropopause)
		tropohgt<-subset(tropohgt,13:20)
		# reading air temperature profile
		ta_press <- calc(brick(readGDAL(sds[15], as.is = TRUE,silent = T)), fun = function(x){((x+15000) * 0.009999999776482582)-273.15})
		# subsetting readings max aronud ~8KM (lower than tropopause)
		ta_press<-subset(ta_press,13:20)
		extent(ta_press)<-extent(bbox)
		projection(ta_press)<-wgs
		extent(tropohgt)<-extent(bbox)
		projection(tropohgt)<-wgs
		# calculate adiabatic lapse rate		
		# fr<-function(x,y) {
		# 	
		# 	frvect<-function(x,y){if(length(x[is.na(x)])==length(x)){
		# 		# cbind(rep(NA,length(x)),rep(NA,length(x)))
		# 		c(NA,NA)
		# 	}else{
		# 		coefs<-lm(y~x)
		# 		# coefs<-mapply(function(x,y)lm(y~x)$coefficients,split(x, row(x)) ,split(y, row(y)))
		# 		return(coefs$coefficients)
		# 	}
		# 	}
		# 	coefs<-t(mapply(frvect,split(x, row(x)) ,split(y, row(y))))
		# }
		# regr<-overlay(x=tropohgt,y=ta_press,fun=fr,forcefun=T)
		return(list(ta_prof=ta_press,tropohgt=tropohgt))
	}
		
}


############################# func to read modis cloud############################
readMOD06<-function(filename){
	# returns a list with 3 datasets: cth= cloud top heigth (m),cl_bs=cloud base height(m) ,tcloud=cloud top temperature (C)
	# testting
	# filename=filenames_cld[4]
	# end testing
	#read datatsets
	sds <- get_subdatasets(filename)
	md<-gdalinfo(filename)
	north<-do.call(c,regmatches(md, gregexpr(paste0('.*.','NORTHBOUNDINGCOORDINATE','.*.'),md)))
	south<-do.call(c,regmatches(md, gregexpr(paste0('.*.','SOUTHBOUNDINGCOORDINATE','.*.'),md)))
	east<-do.call(c,regmatches(md, gregexpr(paste0('.*.','EASTBOUNDINGCOORDINATE','.*.'),md)))
	west<-do.call(c,regmatches(md, gregexpr(paste0('.*.','WESTBOUNDINGCOORDINATE','.*.'),md)))
	bbox<-c(west,east,south,north)
	# bbox<-c(md[86],md[26],md[76],md[55])
	bbox<-do.call(rbind,strsplit(bbox,'='))
	bbox<-as.numeric(bbox[,2])
	wgs<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
	# read cloud tod height
	cth <- raster(readGDAL(sds[58], as.is = TRUE,silent = T))
	extent(cth)<-extent(bbox)
	projection(cth) <- wgs
	# read cloud optical thickness
	copt<- raster(readGDAL(sds[73], as.is = TRUE,silent = T))*0.009999999776482582
	extent(copt)<-extent(bbox)
	projection(copt) <- wgs
	# read cloud effective radius
	cefr<- raster(readGDAL(sds[67], as.is = TRUE,silent = T))*0.009999999776482582
	extent(cefr)<-extent(bbox)
	projection(cefr) <- wgs
	# calculate cloud base height(m)  Welch, et al. 2008 doi 10.1175/2007JAMC1668.1
	calc_bh<-function(cot,cer,clth){
		cldbh<-clth-((4/3)*cot*cer*(0.002)^-1)^0.5
		ifelse(cldbh<0,0.0,cldbh)
	}
	cl_bs<-overlay(copt,cefr,cth,fun=calc_bh)
	# read cloud top temperature
	tcloud<-stack(readGDAL(sds[104], as.is = TRUE,silent = T))
	tcloud<-((tcloud+15000)* 0.009999999776482582)-273.15
	extent(tcloud)<-extent(bbox)
	projection(tcloud) <- wgs
	tcloud<-projectRaster(tcloud,cth)
	# cloud water path(g/m2)
	cwp<- raster(readGDAL(sds[83], as.is = TRUE,silent = T))
	extent(cwp)<-extent(bbox)
	projection(cwp) <- wgs
	# cloud water/vapour content g/m3
	cl_a<-overlay(cth,cl_bs,cwp,fun=function(x,y,z){a<-z/(x-y);ifelse(is.infinite(a),0.0,a)})
	return(list(cth=cth,cl_bs=cl_bs,tcloud=tcloud,cl_a=cl_a))
	
}
############################# func to read SSM lst############################
read_ssm<-function(filename){
	# testing
	# filename=filenames_SSM[1]
	# end testing
	con <- gzfile(filename, open = "rb")
	raw <- readBin(con, what = "int",size=2,n=1e7)
	close(con)
	raw[raw==0]<-NA
	raw<- matrix(data = raw, nrow = 586, ncol = 1383,byrow = TRUE)
	r<-raster(raw, xmn=-17324563.84, xmx=17324563.84, ymn=-7338939.46, ymx=7338939.46,crs='+init=epsg:3410')
	r<-projectRaster(r,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
	calc(r,function(x){y=(x/10);(1.11*y-15.2)-273.15})
}
############################# reading the data ############################
##reading temperature and humidity profiles
cat(' ',"\n")
cat('reading temperature and humidity profiles',"\n")
atm<-mapply(FUN=readMOD07,filenames_atm,MoreArgs=list(output='Ta'),SIMPLIFY = F)
##read atm
# atmtest<-readMOD07(filenames_atm[1],output='Ta')
##read clouds
if (outmode$use.clouds==TRUE){
	cat('reading cloud info',"\n")
	clds<-mapply(FUN=readMOD06,filenames_cld,SIMPLIFY = F)
	# cat('computing adiabatic lapse rate',"\n")
	T_prof<-mapply(FUN=readMOD07,filenames_atm,MoreArgs=list(output='LR'),SIMPLIFY = F)
}

##read modis lst
cat(' ',"\n")
cat('reading modis LST',"\n")
lst_mod<-mapply(FUN=readlst,filenamlst,SIMPLIFY = F)
##read ssm lst
if(outmode$monthly==FALSE){
	cat(' ',"\n")
	cat('reading SSM LST',"\n")
	##read ssm lst
	lst_ssm<-mapply(FUN=read_ssm,filenames_SSM,SIMPLIFY = F)
	lst_ssm<-mapply(crop,x=lst_ssm,MoreArgs = list(y=lst_mod[[1]]),SIMPLIFY = F)
	lst_ssm<-approxNA(stack(lst_ssm),rule=2)
	########################################################################
	#1. downscale lst
	########################################################################
	# get nday per month
	nds<-format(as.Date(dateseq,format='%Y.%m.%d'),'%m')
	nds<-as.data.frame(table(nds))
	# get array to do mapply
	lsmpd<-unlist(mapply(replicate,nds$Freq,lst_mod,SIMPLIFY = F))
	# get funct to downs
	cat('downscaling LST',"\n")
	dowsn<-function(x,y){downscaleRLM(x,y)$downscale}
	# gapfill
	lst_ssm<-mapply(dowsn,as.list(lst_ssm),lsmpd,SIMPLIFY = F)
}

gc()
########################################################################
#3.Organize the layers
########################################################################
# stack the grids
#air temperature at surface
ta<-list()
for(i in 1:length(atm)){ta[[i]]<-atm[[i]]$ta}
#air absolute humidity kg/kg at surface
a<-list()
for(i in 1:length(atm)){a[[i]]<-atm[[i]]$vw_surf}



if (outmode$use.clouds==TRUE){
	#air temperature profile
	ta_prof<-list()
	for(i in 1:length(T_prof)){ta_prof[[i]]<-T_prof[[i]]$ta_prof}
	#elevation at temperature profile readings
	tropo_temp<-list()
	for(i in 1:length(T_prof)){tropo_temp[[i]]<-T_prof[[i]]$tropohgt}
	#cloud top height
	cl_top_hgt<-list()
	for(i in 1:length(clds)){cl_top_hgt[[i]]<-clds[[i]]$cth}
	# cloud base heigth
	cl_b_hgt<-list()
	for(i in 1:length(clds)){cl_b_hgt[[i]]<-clds[[i]]$cl_bs}
	# cloud top temperature
	cl_t<-list()
	for(i in 1:length(clds)){cl_t[[i]]<-clds[[i]]$tcloud}
	# cloud top absolut humidity g/m3
	cl_a<-list()
	for(i in 1:length(clds)){cl_a[[i]]<-clds[[i]]$cl_a}
	
}

# fix formats and extent
tiddyup<-function(rastrlist,mold){
	result<-mapply(FUN=extend,rastrlist,MoreArgs=list(y=mold),SIMPLIFY = F)
	result<-mapply(FUN=crop,result,MoreArgs=list(y=mold),SIMPLIFY = F)
	result<-mapply(FUN=projectRaster,result,MoreArgs=list(to=mold),SIMPLIFY = F)
	# result<-mapply(FUN=mask,result,MoreArgs=list(mask=mold))
	# stack(result)
	result
}
ta<-tiddyup(ta,lst_mod[[1]])
a<-tiddyup(a,lst_mod[[1]])
if (outmode$use.clouds==TRUE){
	cl_top_hgt<-tiddyup(cl_top_hgt,lst_mod[[1]])
	cl_b_hgt<-tiddyup(cl_b_hgt,lst_mod[[1]])
	cl_t<-tiddyup(cl_t,lst_mod[[1]])
	cl_a<-tiddyup(cl_a,lst_mod[[1]])
	ta_prof<-tiddyup(ta_prof,lst_mod[[1]])
	tropo_temp<-tiddyup(tropo_temp,lst_mod[[1]])
}


########################################################################
#3.gap fill daily adiabatic Lapse rate asuming meain in 3 neigthbour pixels
########################################################################
gapfill<-function(x){
	fill.na <- function(x, i=5) {
		if( is.na(x)[i] ) {
			return( mean(x, na.rm=TRUE) )
		} else {
			return(x[i])
		}
	}
	fillLR<-focal(x,w = matrix(1,3,3), fun = fill.na, pad = TRUE, na.rm = FALSE) 
	for (i in 1:2){
		fillLR<-focal(fillLR,w = matrix(1,3,3), fun = fill.na, pad = TRUE, na.rm = FALSE)
	}
	return(fillLR)
}

gc()
dem<-projectRaster(dem,lst_mod[[1]])

if (outmode$use.clouds==TRUE){
	ta<-mapply(gapfill,ta,SIMPLIFY = F)
	a<-mapply(gapfill,a)
	########################################################################
	#3.average temperature and height geopotential profiles
	########################################################################
	###stack profiles per day##
	ta_prof<-tapply(ta_prof,as.Date(unique(file_urls[,2])),stack,simplify = FALSE)
	tropo_temp<-tapply(tropo_temp,as.Date(unique(file_urls[,2])),stack,simplify = FALSE)
	calc_avg<-function(x){
		indmean<-rep(1:8,nlayers(x)/8)
		x<-setZ(x,indmean)
		zApply(x,by=indmean,fun=mean,na.rm=T)
	}
	ta_prof<-lapply(ta_prof,calc_avg)
	tropo_temp<-lapply(tropo_temp,calc_avg)
	########################################################################
	#3.calculate the daily adiabatic lapse rate
	########################################################################
	cat('computing daily adiabatic lapse rate',"\n")
	
	fr<-function(x,y) {
		
		frvect<-function(x,y){if(length(x[is.na(x)])==length(x)){
			# cbind(rep(NA,length(x)),rep(NA,length(x)))
			c(NA,NA)
			}else{
				coefs<-lm(y~x)
				# coefs<-mapply(function(x,y)lm(y~x)$coefficients,split(x, row(x)) ,split(y, row(y)))
				return(coefs$coefficients)
			}
		}
		coefs<-t(mapply(frvect,split(x, row(x)) ,split(y, row(y)),SIMPLIFY = T))
	}
	calc_lr<-function(x,y){overlay(x,y,fun=fr,forcefun=T)}
	
	LR<-mapply(calc_lr,tropo_temp,ta_prof,SIMPLIFY = F)
	alr<-list()
	for(i in 1:length(LR)){alr[[i]]<-LR[[i]][[2]]}
	alr<-mapply(gapfill,alr,SIMPLIFY = F)
	alr<-approxNA(stack(alr),rule=2)
	########################################################################
	# dem<- raster("C:/Base_de_datos_GIS/global_data/elev_1km/elevation_masked_1KMmd_GMTEDmd.tif")
	
	########################################################################
	#3.calculate temperature under the clouds and map orographic inmmersion following Nair et al., 2008 doi: 10.1175/2007JAMC1819.1
	########################################################################
	#######################get air temperature at surface under the low clouds########################
	calc_Ta<-function(tcld,dem,cldBH,cldTH,cld_elev_prof,alr){
		#calc density air at top cloud heigh (Kg/m3(barometric formula, Po,To,g,L,R,M: constants at STP):da=(PoM/(RT))*(1-(Lh/To))^((gM/RL)-1),L adiabatic lapse rate
		# mask out clouds higher than the tropopause
		cldBH<-ifelse(cldBH>cld_elev_prof[,1],NA,cldBH)
		cldTH<-ifelse(cldTH>cld_elev_prof[,1],NA,cldTH)
		# recalc get the cloud thickness
		z_cld<-cldTH-cldBH
		z_cld<-ifelse(z_cld<0,0,z_cld)
		# get altitude from the base of the cloud to surface (m)
		z<-cldBH-dem
		# temperature in kelvin
		TclK<-tcld+273.15
		# get air density at the BASE of the cloud kg/m3
		pda<-((101325*0.0289654)/(8.3144*(TclK)))*(1-((-alr*cldTH)/288.15))^(((9.8*0.0289654)/(8.3141*0.0065))-1)
		# get cloud water vap pressure [pa]
		es<-0.6108*1000*exp((17.27*tcld)/(tcld+237.3))
		# get vapour density g/m3 moist air Oke, 1996 
		pwv<-2.17*(es/TclK)
		# get mixing ratio of the cloud Kg/Kg
		mr<-(pwv/pda)/1000
		# get adiabatic lapse rate K/m http://glossary.ametsoc.org/wiki/Saturation-adiabatic_lapse_rate
		L<-9.8*((287*(TclK)^2+ 2501000*mr*(TclK))/(1003.5*287*(TclK)^2+2501000^2*mr*0.622))
		# get altitude from the top of the cloud to surface (m)
		
		#calc air temperature at the base of the cloud
		Ta_bh<-tcld+(L*z_cld)
		# calc air tmperature at surface from the base of the clouds
		Ta_s<-Ta_bh-(alr*z)
		Ta_s<-ifelse(z<=0,Ta_bh,Ta_s)
		Ta_s<-ifelse(Ta_s<=-15,NA,Ta_s)
		return(Ta_s)
	}
	calc_Ta_rast<-function(tcld,dem,cldBH,cldTH,cld_elev_prof,alr){overlay(tcld,dem,cldBH,cldTH,cld_elev_prof,alr,fun=calc_Ta)}
	########################################################################
	#apply the daily mean ADL to all cloud temperatures
	########################################################################
	################################prepare array for each cloud image
	nimag<-as.data.frame(table(as.character(as.Date(zdates_atm))))
	
	alr_list<-unlist(mapply(replicate,nimag$Freq,as.list(alr),SIMPLIFY = F))
	tropo_temp_list<-unlist(mapply(replicate,nimag$Freq,tropo_temp,SIMPLIFY = F))
	################################ calc air temperature under the clouds at surface level
	Ta_cld_s<-mapply(calc_Ta_rast,tcld=cl_t,cldBH=cl_b_hgt,cldTH=cl_top_hgt,cld_elev_prof=tropo_temp_list,alr=alr_list,MoreArgs = list(dem=dem),SIMPLIFY = F)
	################################ calc air avg temperature per image
	calc_avgTa<-function(x,y){overlay(x,y,fun=function(x,y){rowMeans(cbind(x,y),na.rm = T)})}
	
	Ta<-mapply(calc_avgTa,ta,Ta_cld_s,SIMPLIFY = F)
	Ta<-brick(Ta)
	Ta<-approxNA(Ta)
	Ta<-setZ(Ta,zdates_atm)
	Ta<-zApply(x=Ta,by=as.Date(zdates_atm),fun=mean,na.rm=T)
	
	# Ta<-stackApply(stack(Ta),as.Date(zdates_atm),fun=mean,na.rm=T)
	########################################################################
	#calc mixing ratio [kg/kg]under the clouds, theoretical at Ta under the clouds
	# ########################################################################
	calc_a_clds<-function(lr,t){overlay(lr,t,fun=function(lr,t){
		tk<-t+273.15
		a<-(287*tk^2/2501000)*((9.8+lr*1003.5)/(-lr*2501000*0.622-9.8*tk))
		a<-ifelse(a<0,0,a)
		a
	}
	 )
	 }
	 a_clds<-mapply(calc_a_clds,as.list(alr),as.list(Ta),SIMPLIFY = F)
	
	# ########################################################################
	# #3.compute ea from mixing ratio
	# ########################################################################
	calc_ea<-function(a,t,dem){
		#dry air density[kg/m3]
		pair<-p_dryair(dem,t)
		#get vapour density [kg/m3](Oke, 1996)
		pvap<-a*pair
		ea<-pvap*461.5*(273.15+t)
		es<-0.6108*1000*exp((17.27*t)/(t+237.3))
		ea<-ifelse(ea>es,es,ea)
		ea
		
	}
	ea<-overlay(stack(a),stack(ta),dem,fun=calc_ea)
	ea_clds<-overlay(stack(a_clds),Ta,dem,fun=calc_ea)
	# ea<-approxNA(ea,rule=2)
	ea<-setZ(ea,zdates_atm)
	ea<-zApply(x=ea,by=as.Date(zdates_atm),fun=mean,na.rm=T)
	# ea<-approxNA(ea,rule=2)
	ea<-mapply(calc_avgTa,as.list(ea),as.list(ea_clds),SIMPLIFY = F)
	ea<-stack(ea)
	########################################################################
}else{
	# ########################################################################
	# #3.compute daily averages
	# ########################################################################
	cat('gapfilling',"\n")
	ta<-mapply(gapfill,ta,SIMPLIFY = F)
	a<-mapply(gapfill,a,SIMPLIFY = F)
	lst_mod<-mapply(gapfill,lst_mod,SIMPLIFY = F)
	# save(ta,file=paste0(outdir,'/','Ta_day_',hv[1],'.',beg,'.',til,'.RData'))
	# save(zdates_atm,file=paste0(outdir,'/','zdates_atm',hv[1],'.',beg,'.',til,'.RData'))
	Ta<-brick(ta)
	Ta<-approxNA(Ta,rule=2)
	# Ta<-stackApply(Ta,as.Date(zdates_atm),fun=mean,na.rm=T)
	Ta<-setZ(Ta,zdates_atm)
	Ta<-zApply(Ta,as.Date(zdates_atm),fun=mean,na.rm=T)
	

	a<-brick(a)
	a<-approxNA(a,rule=2)
	# Ta<-stackApply(Ta,as.Date(zdates_atm),fun=mean,na.rm=T)
	a<-setZ(a,zdates_atm)
	a<-zApply(a,as.Date(zdates_atm),fun=mean,na.rm=T)
	# a<-stackApply(a,as.Date(zdates_atm),fun=mean,na.rm=T)
	# a<-approxNA(a,rule=2)
	# a<-approxNA(a,rule=2)
	# ########################################################################
	# #3.compute ea from mixing ratio
	# ########################################################################
	calc_ea<-function(a,t,dem){
		#dry air density[kg/m3]
		pair<-p_dryair(dem,t)
		#get vapour density [kg/m3](Oke, 1996)
		pvap<-a*pair
		ea<-pvap*461.5*(273.15+t)
		es<-0.6108*1000*exp((17.27*t)/(t+237.3))
		ea<-ifelse(ea>es,es,ea)
		ea
		
	}
	ea<-overlay(a,Ta,dem,fun=calc_ea)
	
}

if(outmode$monthly==TRUE){
	# get es* [pa]
	es<-calc(stack(lst_mod),fun=function(ts)0.6108*1000*exp((17.27*ts)/(ts+237.3)))
	monthind<-format(getZ(Ta),'%Y-%m')
	Ta<-zApply(Ta,monthind,fun=mean,na.rm=T)
}else{
	es<-calc(stack(lst_ssm),fun=function(ts)0.6108*1000*exp((17.27*ts)/(ts+237.3)))
}
# 
# ########################################################################
# #3.compute mixing ratio under the clouds
# ########################################################################
# 
# 
# #################calc moist air density clear sky pixels kg/m3 ###################
# calcpma<-function(ah){calc(ah,fun=function(ah){1.2754*((1+ah)/(1+1.609*ah))})}
# pma<-mapply(calcpma,a)
# #################calc water vapour pressure kg/m3 ###################
# calc_ea<-function(a,pma,cld_a,ta,ta_cl_s,dem){
# 	# temperature in kelvin
# 	TclK<-ta_cl_s+273.15
# 	# calc ea where cloud inmmersion happens assuming air is saturated
# 	ea_cl<-0.6108*1000*exp((17.27*ta_cl_s)/(ta_cl_s+237.3))
# 	#calc ea under the cloud, using vapour density Oke, 1996
# 	# ea_cl<-cld_a*TclK*461.5
# 	
# 	# get air density at the top of the cloud kg/m3
# 	# pda<-((101325*0.0289654)/(8.3144*(TclK)))*(1-((0.0065*dem)/288.15))^(((9.8*0.0289654)/(8.3141*0.0065))-1)
# 	# # density of moist air in cloudy pixels kg/m3
# 	# pma_cl<-pda+(cld_a/1000)
# 	# # water vapor pressure  at surface under cloudy pixels [Pa]
# 	# ea_cl<-(pma_cl*8.314*TclK)/(18.015/1000)
# 	# water vapor pressure  at surface clear pixels [Pa]
# 	eair<-(a*1000)*(1/18.015)*pma*8.314*(ta+273.15)
# 	# water vapor pressure
# 	ea<-rowMeans(cbind(eair,ea_cl),na.rm = T)
# 	return(ea)
# }
# calc_ea_rast<-function(a,pma,cld_a,ta,ta_cl_s,dem){overlay(a,pma,cld_a,ta,ta_cl_s,dem,fun=calc_ea)}
# ea<-mapply(FUN=calc_ea_rast,a=a,pma=pma,cld_a=cl_a,ta=ta,ta_cl_s=Ta_cld_s,MoreArgs = list(dem=dem))
# 
# calc_avgTa<-function(x,y){overlay(x,y,fun=function(x,y){rowMeans(cbind(x,y),na.rm = T)})}
# 
# Ta<-mapply(calc_avgTa,ta,Ta_cld_s)
# ########################################################################
# #3.interpolate using gwr
# ########################################################################
# Ta<-setZ(stack(Ta),as.Date(zdates_atm))
# Ta<-zApply(x=Ta,by=as.Date(zdates_atm),fun=mean,na.rm=T)
# 
# ea<-setZ(stack(ea),as.Date(unique(file_urls[,2])))
# ea<-zApply(x=ea,by=as.Date(unique(file_urls[,2])),fun=mean,na.rm=T)

# testsgwrd<-stack(Ta[[1]],dem)
# testsgwrd<-rasterToPoints(testsgwrd, spatial=T)
# names(testsgwrd)<-c('ta','elev')
# testsgwrd<-subset(testsgwrd,!is.na(testsgwrd$ta) & !is.na(testsgwrd$elev))
# dist.v1<-gw.dist(dp.locat=testsgwrd@coords, p=2, theta=0, longlat=FALSE)
# bands<-bw.gwr(ta~elev,data=testsgwrd,p=2, theta=0, longlat=FALSE,dMat =dist.v1)
# 
# # testtcl<-overlay(cl_t[[2]],cl_a[[2]],dem,cl_top_hgt[[2]],fun=calc_Ta)
# 
# 
# tile<-do.call(rbind,strsplit(urlext[1],'.',fixed = T))[,7]
# beg<-format(as.Date(start),'%Y%j')
# til<-format(as.Date(end),'%Y%j')
# 
# 
# 
# ta<-setZ(ta,as.Date(unique(file_urls[,2])))
# ta<-zApply(x=ta,by=as.Date(unique(file_urls[,2])),fun=mean,na.rm=T)
# windows()
# plot(ta[[1]])
# 
# 
# 
# plot(crop(dem,ta))
# 	
# 			
# calc_tavg<-function(daylist,timeoverpass){
# 	
# }
# 
# 
# lst_ssmtest<-projectRaster(lst_ssm[[1]],lst_mod[[1]])
# 
# 
# 
# extent(lst_ssmtest)<-extent(lst_ssm[[1]])
# 
# lstdown<-downscaleRLM(lst_ssmtest,lst_mod[[1]])
# 
# dem<- raster("C:/Base_de_datos_GIS/global_data/elev_1km/elevation_masked_1KMmd_GMTEDmd.tif")
# gc()
# 
# dem<-projectRaster(crop(dem,lst_mod[[1]]),lst_mod[[1]])
# 
# # downscale krigging
# 
# lst_ssmtest<-crop(stack(lst_ssm),lst_mod[[1]])
# 
# lst_ssm_sp<-rasterToPoints(lst_ssmtest,spatial=T)
# lst_ssm_sp$code<-paste0('cod_',1:length(lst_ssm_sp[,1]))
# lst_ts<-xts(t(as.matrix(lst_ssm_sp@data[,1:2])),seq(as.Date(start),as.Date(end),by='day'))
# 
# tpn<-interpolateForcing(lst_ssm_sp,lst_ts,lst_mod[[1]],tmpdir)
# 
# 
# hmt_stat<-readOGR(dsn="C:/Water_Data/iMHEA/iMHEA_geodata/HMT", layer= "stations",pointDropZ=TRUE)
# hmt_stat$code<-hmt_stat@data$Código_iM
# save(prececu2012,file = "C:/Water_Data/iMHEA/iMHEA_geodata/HMT/prec_HMT_2015_2016.RData")
# load(file="C:/Water_Data/iMHEA/iMHEA_geodata/HMT/prec_HMT_2015_2016.RData")
# tpn<-interpolateForcing(hmt_stat,prececu2012,dem,"C:/Water_Data/iMHEA/iMHEA_geodata/HMT")
# 
# 
# 
# 
# # get surface temperature C
# Tsurf_mod11<-mapply(FUN=readlst,filenamlst)						
# Tsurf<-mapply(FUN=readMOD07,filenames,MoreArgs=list(output='Tsurf'))
# # get air temperature at surface C
# Ta<-mapply(FUN=readMOD07,filenames,MoreArgs=list(output='Ta'))
# 
# # get absolute humidity at surface kg/kg
# a<-mapply(FUN=readMOD07,filenames,MoreArgs=list(output='a'))
# gc()
# # fix formats and extent
# tiddyup<-function(rastrlist,mold){
# 	result<-mapply(FUN=extend,rastrlist,MoreArgs=list(y=mold))
# 	result<-mapply(FUN=crop,result,MoreArgs=list(y=mold))
# 	result<-mapply(FUN=projectRaster,result,MoreArgs=list(to=mold))
# 	# stack(result)
# 	result
# }
# Tsurf<-tiddyup(Tsurf,Tsurf_mod11[[1]])
# Ta<-tiddyup(Ta,Tsurf_mod11[[1]])
# a<-tiddyup(a,Tsurf_mod11[[1]])
# ########################################################################
# #4. gap filling
# ########################################################################	
# ind.complete<-file_urls[,2]

# Ta<-mapply(FUN=focal,Ta,MoreArgs = list(w = matrix(1,3,3), fun = fill.na, pad = TRUE, na.rm = FALSE),SIMPLIFY = F)
# Tsurf<-mapply(FUN=focal,Tsurf,MoreArgs = list(w = matrix(1,3,3), fun = fill.na, pad = TRUE, na.rm = FALSE),SIMPLIFY = F)
# a<-mapply(FUN=focal,a,MoreArgs = list(w = matrix(1,3,3), fun = fill.na, pad = TRUE, na.rm = FALSE),SIMPLIFY = F)
# gc()
# ########################################################################
# #6.get daily averages
# ########################################################################
# 
# Tsurf<-setZ(stack(Tsurf),as.Date(file_urls[,2]))
# Tsurf<-zApply(x=Tsurf,by=as.Date(file_urls[,2]),fun=mean,na.rm=T)
# avg<-function(x,y){
# 	ifelse(is.na(x),ifelse(is.na(y),NA,y),ifelse(is.na(y),x,(x+y)/2))
# }
# Tsurf<-overlay(Tsurf,stack(Tsurf_mod11),fun=avg)
# Ta<-setZ(stack(Ta),as.Date(file_urls[,2]))
# Ta<-zApply(x=Ta,by=as.Date(file_urls[,2]),fun=mean,na.rm=T)
# a<-setZ(stack(a),as.Date(file_urls[,2]))
# a<-zApply(x=a,by=as.Date(file_urls[,2]),fun=mean,na.rm=T)
# # calculate vpd surface - air
# calc_vpd<-function(t,ts,a){
# 	# get es* [pa]
# 	es<-0.6108*1000*exp((17.27*ts)/(ts+237.3))
# 	# density moist air [kg/m3] https://www.engineeringtoolbox.com/density-air-d_680.html
# 	pma<-1.2754*((1+a)/(1+1.609*a))
# 	# water vapor pressure [Pa]
# 	ea<-(a*1000)*(1/18.015)*pma*8.314*(t+273.15)
# 	# vdp<-ifelse(ea>es,0.0,es-ea)
# 	vpd<-es-ea
# 	vpd<-ifelse(vpd<0,0.0,vpd)
# 	return(vpd)
# }
# vpd<-overlay(Ta,Tsurf,a,fun=calc_vpd,forcefun=T)
# vpd<-setZ(vpd,getZ(Ta))	
Ta<-writeRaster(Ta,paste0(outdir,'/','Ta_day_',hv[1],'.',beg,'.',til,'.nc'),format="CDF",overwrite=TRUE,varname="Ta", varunit="C", longname="air temperature", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(as.numeric(format(as.Date(start),'%Y'))-1,"-",12,"-",31)))
ea<-writeRaster(ea,paste0(outdir,'/','ea_day_',hv[1],'.',beg,'.',til,'.nc'),format="CDF",overwrite=TRUE,varname="ea", varunit="Pa", longname="actual vapor pressure", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(as.numeric(format(as.Date(start),'%Y'))-1,"-",12,"-",31)))
es<-writeRaster(es,paste0(outdir,'/','es_day_',hv[1],'.',beg,'.',til,'.nc'),format="CDF",overwrite=TRUE,varname="es", varunit="Pa", longname="saturation vapor pressure", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(as.numeric(format(as.Date(start),'%Y'))-1,"-",12,"-",31)))


return(list(Ta=Ta,ea=ea,es=es))									
					
					
}




























