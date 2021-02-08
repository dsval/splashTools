#' getForcingNRCS
#'
#' Wraper of the RNRCS and daymetr packages to get forcing and soil moisture data at NRCS sites
#' @param  fields on stations.sp: lat, lon, code
#' @import RNRCS 
#' @import xts
#' @importFrom zoo na.approx
#' @importFrom daymetr download_daymet
#' @importFrom matrixStats rowWeightedMeans
#' @keywords splashTools
#' @export
#' @examples
#' splash.grid()

getForcingNRCS<-function(site_id,lat,lon){
	# download and process snotel and scan networks data, wrapper of grabNRCS.data
	# require(RNRCS)
	# require(daymetr)
	# require(xts)
	# require(matrixStats)
	# testing
	# site_id=US_sntl_stations$site_id[318]
	# lat<-US_sntl_stations$latitude[318]
	# lon<- US_sntl_stations$longitude[318]
	###################################################################################################
	# grab metadata
	try(meta_sntl<-grabNRCS.elements(site_id = site_id))
	if(is.na(meta_sntl[[1]])){
		result<-NULL
	}else{
	try(upper_yr_sm<-grep("Soil Moisture.*.", meta_sntl[[1]][,1]))
	if(length(upper_yr_sm)==0||is.null(meta_sntl)){
		result<-NULL	
	}else{
	###################################################################################################	
	# get start day for soil moisture
	y1<-as.Date(meta_sntl[[1]][upper_yr_sm,2],format="%Y-%m-%d")
	y1<-format(y1,format="%Y")
	###################################################################################################
	# grab data
	data_sntl_error<-try(data_sntl<-grabNRCS.data(network=do.call(rbind, strsplit(site_id, ':'))[,1], site_id=do.call(rbind, strsplit(site_id, ':'))[,2], timescale="daily",DayBgn = paste0(y1,"-","01","-","01"), DayEnd = '2015-12-31'), silent=TRUE)
	
	if(class(data_sntl_error)=="try-error"){
		result<-NULL
		message('online dataset not available!!')
		return(result)
	}
	
	# if(is.null(data_sntl)){
	# 	result<-NULL
	# 	}else{
	# 		
	# 	}
	###################################################################################################
	# define temperature C
	Ta<-xts((data_sntl$Air.Temperature.Average..degF.-35)*(5/9),as.Date(data_sntl[,1],format="%Y-%m-%d"))
	Ta[Ta>40]<-NA
	Ta[Ta<=-30]<-NA
	###################################################################################################			
	# define precipitation[mm]
	if(is.null(data_sntl$Precipitation.Accumulation..in..Start.of.Day.Values)){
		P_raw=Ta*NA
		P=Ta*NA
	}else{
		P_raw<-xts(data_sntl$Precipitation.Accumulation..in..Start.of.Day.Values*25.4,as.Date(data_sntl[,1],format="%Y-%m-%d"))
		p<-tapply(P_raw,format(time(P_raw),"%j"),mean, na.rm=TRUE)
		p<-tapply(P_raw,format(time(P_raw),"%Y"),function(x){ifelse(!is.na(x),x,p)})
		p<-do.call(c,p)
		P<-na.approx(xts(p,time(P_raw)),na.rm=FALSE)	
		# data_sntl<-xts(data_sntl[,2:18],as.Date(data_sntl[,1],format="%Y-%m-%d"))
		# grab previous year data
		
		# function to get rainfal from cumsum geiven by sntl
		deacum_prec<-function(acumts){
			deacumts<-rep(NA,length(acumts))
			for(i in 1:length(acumts)){
				j<-i+1
				if(is.na(acumts[i])){
					deacumts[i]<-NA	
				}else{
					if(i==length(acumts)){
						deacumts[i]<-0
					}else{
						if(acumts[j]>acumts[i]){
							deacumts[i]<-acumts[j]-acumts[i]
						}
						else{
							deacumts[i]<-0
						}
					}
				}
				
				
			}
			
			
			return(deacumts)
		}
		# get precipitation[mm/day]
		P<-deacum_prec(as.numeric(P))
		P<-xts(P,time(P_raw))
		P[is.na(P_raw)]<-P_raw[is.na(P_raw)]
		
		P[P>100]<-NA
		
		
	}
	###################################################################################################			
	# exclude incomplete years
	if(length(P_raw)<=366){
		result<-NULL
	}else{
	###################################################################################################			
	# get soil moisture 1st layer [mm]	
	varnames<-names(data_sntl)
	sm.names<-grep(".*Moisture.*.", varnames,value=TRUE)
	if(length(sm.names)==0){
		result<-NULL
	}else{
		sm.names.loc<-grep(".*Moisture.*.", varnames)
		depths<-do.call(rbind, strsplit(sm.names, '..',fixed = TRUE))[,2]
		depths<-do.call(rbind, strsplit(depths, 'in'))
		depths<-as.numeric(depths)
		depths.from<-depths
		
		for(i in 1:length(depths)){
			if(i==1){
				depths.from[i]<-0
			}else{
				depths.from[i]<-depths[i-1]
			}
			
		}
		depth.interv<-depths-depths.from
		depth.total<-sum(depth.interv)
		
		
		# sm to from perc per in to mm
		sm<-list()
		for(i in 1:length(depths)){
			sm[[i]]<-data_sntl[,sm.names.loc[i]]
		}
		sm<-Reduce(cbind,sm)
		sm[sm<0]<-NA
		
		depth_weigths<-depth.interv/depth.total
		if(length(depth_weigths)==1){
			sm<-sm*depth.total*24.5/100
		}else{
			sm<-rowWeightedMeans(x=sm,w=depth_weigths,na.rm = F)*depth.total*24.5/100
		}
		
		sm[sm>700]<-NA
		sm[sm<0]<-NA
		
		sm<-xts(sm,time(Ta))
		###################################################################################################	
		# get snow swe[mm] if swe not available, estimate from snowdepth, assuming 250 Kg/m3 as snow density
		if(is.null(data_sntl$Snow.Water.Equivalent..in..Start.of.Day.Values)){
			swe<-xts(data_sntl$Snow.Depth..in..Start.of.Day.Values*(250/1000)*24.5,time(Ta))
		}else{
			swe<-xts(data_sntl$Snow.Water.Equivalent..in..Start.of.Day.Values*24.5,time(Ta))
		}
		if(length(swe)==0){
			swe<-Ta*NA
		}else{
			swe[swe>2000]<-NA
			swe[swe<0]<-0
		}
		###################################################################################################	
		########################## gap fill with DAYMET ################################
		###################################################################################################	
		# get mean solar radiation from daymet * daymet provides mean rad per dayligth time, so rad mean= rad_dayligth*dayligthduration/86400
		daymet_data <- daymetr::download_daymet("Daymet",
			lat = lat,
			lon = lon,
			start = y1,
			end = 2015,
			internal = TRUE)
		ind_daym<-as.Date(paste0(daymet_data$data$year,"-",daymet_data$data$yday),format="%Y-%j")
		# get radiation in W/m2
		sw_in<-xts((daymet_data$data$srad..W.m.2.*daymet_data$data$dayl..s.)/86400,ind_daym)
		###################################################################################################	
		# gap fill precipitation and temperature with daymet data
		P_daym<-xts(daymet_data$data$prcp..mm.day.,ind_daym)
		T_daym<-xts((daymet_data[[7]]$tmax..deg.c.+daymet_data[[7]]$tmin..deg.c.)/2,ind_daym)
		P[is.na(P_raw[time(P_daym)])]<-P_daym[is.na(P_raw[time(P_daym)])]
		Ta[is.na(Ta[time(T_daym)])]<-T_daym[is.na(Ta[time(T_daym)])]
		#swe.daymet<-xts(daymet_data$data$swe..kg.m.2.*4,ind_daym)
		#swe[is.na(swe[time(swe.daymet)])]<-swe.daymet[is.na(swe[time(swe.daymet)])]
		data<-merge.xts(P,Ta,sw_in,swe,sm)
		###################################################################################################	
		# remove incomplete first year observations, subset from the next complete year	
		if(sum(!is.na(data[,5][y1]))>=365){
			data<-data
		}else{
			y1<-as.character(as.numeric(y1)+1)
			data<-data[paste0(y1,"/")]
		}
		#data<-na.approx(data)
		NonNAindex <- which(!is.na(data[,1]))
		firstNonNA <- min(NonNAindex)
		ynonna<-format(time(data[,1][firstNonNA]),format="%Y")
		###################################################################################################	
		# remove incomplete first year observations, subset from the next complete year	
		if(sum(!is.na(data[,1][ynonna]))>=365){
			data<-data[paste0(y1,"/")]
		}else{
			ynonna<-as.character(as.numeric(ynonna)+1)
			data<-data[paste0(ynonna,"/")]
		}
		#data<-na.approx(data)
		###################################################################################################	
		########################## gap fill with seasonal averages ################################
		###################################################################################################	
		# get average per day, 365 values
		t<-tapply(data$Ta,format(time(data),"%j"),mean, na.rm=TRUE)
		p<-tapply(data$P,format(time(data),"%j"),mean, na.rm=TRUE)
		sw<-tapply(data$sw_in,format(time(data),"%j"),mean, na.rm=TRUE)
				
		#replace NAs in TS with averages
		t<-tapply(data$Ta,format(time(data),"%Y"),function(x){ifelse(!is.na(x),x,t)})
		p<-tapply(data$P,format(time(data),"%Y"),function(x){ifelse(!is.na(x),x,p)})
		sw<-tapply(data$sw_in,format(time(data),"%Y"),function(x){ifelse(!is.na(x),x,sw)})
		
		#create the vectors
		t<-do.call(c,t)
		p<-do.call(c,p)
		sw<-do.call(c,sw)
		
		
		# replacemissing
		data$Ta<-na.approx(xts(t,time(data)),na.rm = F)
		data$P<-na.approx(xts(p,time(data)),na.rm = F)
		data$sw_in<-na.approx(xts(sw,time(data)),na.rm = F)
		# total depth in meters
		depth.total<-depth.total*24.5/1000
		result<-list(data,depth.total)	
	}
	
	}
	}
	}
	return(result)
}





