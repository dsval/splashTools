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
	# site_id=395
	# lat<-43.23
	# lon<- -121.81
	
	# grab metadata
	try(meta_sntl<-grabNRCS.elements(site_id = site_id))
	if(is.na(meta_sntl[[1]])){
		result<-NULL
	}else{
	try(upper_yr_sm<-grep("Soil Moisture.*.", meta_sntl[[1]][,1]))
	if(length(upper_yr_sm)==0||is.null(meta_sntl)){
		result<-NULL	
	}else{
		
	# get start day for soil moisture
	
	y1<-as.Date(meta_sntl[[1]][upper_yr_sm,2],format="%Y-%m-%d")
	y1<-format(y1,format="%Y")
	# grab data
	data_sntl<-grabNRCS.data(network=do.call(rbind, strsplit(site_id, ':'))[,1], site_id=do.call(rbind, strsplit(site_id, ':'))[,2], timescale="daily",
		DayBgn = paste0(y1,"-","01","-","01"), DayEnd = '2015-12-31')
		
		
	# precipitation[mm]
	P_raw<-xts(data_sntl$Precipitation.Accumulation..in..Start.of.Day.Values*25.4,as.Date(data_sntl[,1],format="%Y-%m-%d"))
	if(length(P_raw)<=366){
		result<-NULL
	}else{
		
	Ta<-xts((data_sntl$Air.Temperature.Average..degF.-35)*(5/9),time(P_raw))
	Ta[Ta>40]<-NA
	Ta[Ta<=-30]<-NA	
	P<-na.approx(P_raw,na.rm=FALSE)	
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
	# get temperature in celcius
	# pcomp<-xts(data_sntl$Precipitation.Increment..in.*25.4,as.Date(data_sntl[,1],format="%Y-%m-%d"))
	
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
			sm<-rowWeightedMeans(x=sm,w=depth_weigths,na.rm = TRUE)*depth.total*24.5/100
		}
		
		sm[sm>700]<-NA
		sm[sm<0]<-NA
		
		sm<-xts(sm,time(Ta))
		# get snow swe[mm] if swe not available, estimate from snowdepth, assuming 250 Kg/m3 as asnow density
		if(is.null(data_sntl$Snow.Water.Equivalent..in..Start.of.Day.Values)){
			swe<-xts(data_sntl$Snow.Depth..in..Start.of.Day.Values*(250/1000)*24.5,time(Ta))
		}else{
			swe<-xts(data_sntl$Snow.Water.Equivalent..in..Start.of.Day.Values*24.5,time(Ta))
		}
		if(length(swe)==0){
			swe<-Ta*0
		}else{
			swe[swe>2000]<-NA
			swe[swe<0]<-0
		}
		
		
		
		# get mean solar radiation from daymet * daymet provides mean rad per dayligth time, so rad mean= rad_dayligth*dayligthduration/86400
		daymet_data <- daymetr::download_daymet("Daymet",
			lat = lat,
			lon = lon,
			start = y1,
			end = 2015,
			internal = TRUE)
		ind_daym<-as.Date(paste0(daymet_data$data$year,"-",daymet_data$data$yday),format="%Y-%j")
		P_daym<-xts(daymet_data$data$prcp..mm.day.,ind_daym)
		P[is.na(P_raw)]<-P_daym[is.na(P_raw)]
		# get radiation in W/m2
		sw_in<-xts((daymet_data$data$srad..W.m.2.*daymet_data$data$dayl..s.)/86400,ind_daym)
		swe.daymet<-xts(daymet_data$data$swe..kg.m.2.*4,ind_daym)
		swe[is.na(swe)]<-swe.daymet[is.na(swe)]
		data<-merge.xts(P,Ta,sw_in,swe,sm)
		# remove incomplete first year observations, subset from the next complete year	
		if(sum(!is.na(data[,5][y1]))>=365){
			data<-data
		}else{
			y1<-as.character(as.numeric(y1)+1)
			data<-data[paste0(y1,"/")]
		}
		data<-na.approx(data)
		NonNAindex <- which(!is.na(data[,5]))
		firstNonNA <- min(NonNAindex)
		ynonna<-format(time(data[,5][firstNonNA]),format="%Y")
		# remove incomplete first year observations, subset from the next complete year	
		if(sum(!is.na(data[,5][ynonna]))>=365){
			data<-data[paste0(y1,"/")]
		}else{
			ynonna<-as.character(as.numeric(ynonna)+1)
			data<-data[paste0(ynonna,"/")]
		}
		data<-na.approx(data)
		# total depth in meters
		depth.total<-depth.total*24.5/1000
		result<-list(data,depth.total)	
	}
	
	}
	}
	}
	return(result)
}





