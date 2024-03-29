#' cloud2Rad
#' Compute mean solar radiation (W/m2 ) by using the classic AP formula and the elevation effect on the atmosphere's transmittance, as described by Davis et al., (2017)
#' @param  elev_lowres, elevation (masl), Raster layer
#' @param  clds, cloudiness, (fraction), Raster* object with z time dimension
#' @param  ouputdir, directory where to save the files
#' @return  A RasterBrick object with the mean solar radiation (W/m2 ) with the same dimensions as clds
#' @import raster  
#' @keywords splash
#' @export
#' @examples
#' *optional run beginCluster() first, for parallel computing
#' cloud2Rad()
cloud2Rad<-function(clds,elev_lowres,ouputdir=getwd(),inmem=FALSE, ...){
	
	###########################################################################
	# 00. Check if parallel computation is required by the user and if the dimensions of the raster objects match
	###########################################################################
	on.exit(endCluster())
	clcheck<-try(getCluster(), silent=TRUE)
	if(class(clcheck)=="try-error"){
		# If no cluster is initialized, assume only one core will do the calculations, beginCluster(1) saved me the time of coding serial versions of the functions
		beginCluster(1,'SOCK')
		message('Only using one core, use first beginCluster(ncpus) if you need to do the calculations in parallel!!')
		
	}
	###############################################################################################
	# 00. create array for results, get time info
	###############################################################################################
	y<-as.numeric(format(getZ(clds),'%Y'))
	doy<-as.numeric(format(getZ(clds),'%j'))
	ny <-length(doy)
	ind<-getZ(clds)
	# sf
	sf<-calc(clds,fun=function(x){(100-x)/100},filename=paste0(ouputdir,"/","sf.grd"), overwrite=TRUE)
	sf<-setZ(sf,ind)
	# radiation hres
	out<-brick(nrows=nrow(clds), ncols=ncol(clds), crs=crs(clds), nl=ny)
	extent(out)<-extent(clds)
	out<-setZ(out,ind)
	# 1.2 get latitudes
	lat_lr<-elev_lowres*0
	lat.data<-rasterToPoints(elev_lowres)
	lat_lr[!is.na(lat_lr)]<-lat.data[,2]
	rm(lat.data)
	
	# 1.3  calculate slope and aspect
	#terraines<-terrain(elev_lowres, opt=c('slope', 'aspect'), unit='degrees')
	setwd(dirname(rasterTmpFile()))
	writeRaster(elev_lowres,"rawdem.tif",datatype='INT2S',format="GTiff", overwrite=TRUE)
	system("gdaldem slope -s 111120 -compute_edges rawdem.tif slope_deg.tif")
	system("gdaldem aspect -zero_for_flat -compute_edges rawdem.tif aspect_deg.tif")
	terraines<-raster::stack(list(slope='slope_deg.tif',aspect='aspect_deg.tif'))
	
	# 1.3  fillnas slope and aspect
	fillna<-function(ind,x){
		focal(x[[ind]], w = matrix(1,3,3), fun = function(x, i=5) {
			if( is.na(x)[i] ) {
				return(mean(x, na.rm=TRUE ))
			} else {
				return( x[i] )
			}
			}, 
			pad = TRUE, na.rm = FALSE )
	}
	
	# n <- seq(nlayers(terraines)) 
	# terraines<-stack(lapply(X=n, FUN=fillna, x=terraines),quick=TRUE)
	setwd(ouputdir)
	
	###########################################################################
	# 02. Define functions
	###########################################################################
	
	calc_sw_in<-function(sf,lat,elev,slop,asp,y,doy){
		
		###########################################################################
		# 01. Define constants inside functions to avoid exporting one by one to the cluster
		###########################################################################
		kA <- 107           # constant for Rl (Monteith & Unsworth, 1990)
		kalb_sw <- 0.17     # shortwave albedo (Federer, 1968)
		kalb_vis <- 0.03    # visible light albedo (Sellers, 1985)
		kb <- 0.20          # constant for Rl (Linacre, 1968; Kramer, 1957)
		kc <- 0.25          # constant for Rs (Linacre, 1968)
		kd <- 0.50          # constant for Rs (Linacre, 1968)
		kfFEC <- 2.04       # from-flux-to-energy, umol/J (Meek et al., 1984)
		kG <- 9.80665       # gravitational acceleration, m/s^2 (Allen, 1973)
		kGsc <- 1360.8      # solar constant, W/m^2 (Kopp & Lean, 2011)
		kL <- 0.0065        # adiabatic lapse rate, K/m (Cavcar, 2000)
		kMa <- 0.028963     # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
		kMv <- 0.01802      # mol. weight of water vapor, kg/mol (Tsilingiris, 2008)
		kSecInDay <- 86400  # number of seconds in a day
		kPo <- 101325       # standard atmosphere, Pa (Allen, 1973)
		kR <- 8.31447       # universal gas constant, J/mol/K (Moldover et al., 1988)
		kTo <- 288.15       # base temperature, K (Berberan-Santos et al., 1997)
		pir <- pi/180       # pi in radians
		# Paleoclimate variables:
		ke <- 0.01670       # eccentricity of earth's orbit, 2000CE (Berger 1978)
		keps <- 23.44       # obliquity of earth's elliptic, 2000CE (Berger 1978)
		komega <- 283       # lon. of perihelion, degrees, 2000CE (Berger, 1978)
		###########################################################################
		# 02. Define functions
		###########################################################################
		# ************************************************************************
		# Name:     julian_day
		# Inputs:   - double, year (y)
		#           - double, month (m)
		#           - double, day of month (i)
		# Returns:  double, Julian day
		# Features: This function converts a date in the Gregorian calendar
		#           to a Julian day number (i.e., a method of consecutative
		#           numbering of days---does not have anything to do with
		#           the Julian calendar!)
		#           * valid for dates after -4712 January 1 (i.e., jde >= 0)
		# Ref:      Eq. 7.1 J. Meeus (1991), Chapter 7 "Julian Day", Astronomical
		#             Algorithms
		# ************************************************************************
		julian_day <- function(y, m, i) {
			if (m <= 2) {
				y <- y - 1
				m <- m + 12
			}
			a <- floor(y/100)
			b <- 2 - a + floor(a/4)
			
			jde <- floor(365.25*(y + 4716)) + floor(30.6001*(m + 1)) + i + b - 1524.5
			return(jde)
		}
		# ************************************************************************
		# Name:     berger_tls
		# Inputs:   - double, day of year (n)
		#           - double, days in year (N)
		# Returns:  numeric list, true anomaly and true longitude
		# Features: Returns true anomaly and true longitude for a given day.
		# Depends:  - ke ............. eccentricity of earth's orbit, unitless
		#           - komega ......... longitude of perihelion
		#  Ref:     Berger, A. L. (1978), Long term variations of daily insolation
		#             and quaternary climatic changes, J. Atmos. Sci., 35, 2362-2367.
		# ************************************************************************
		berger_tls <- function(n, N) {
			# Variable substitutes:
			xee <- ke^2
			xec <- ke^3
			xse <- sqrt(1 - ke^2)
			
			# Mean longitude for vernal equinox:
			xlam <- (ke/2.0 + xec/8.0)*(1 + xse)*sin(komega*pir) -
			xee/4.0*(0.5 + xse)*sin(2.0*komega*pir) +
			xec/8.0*(1.0/3.0 + xse)*sin(3.0*komega*pir)
			xlam <- 2.0*xlam/pir
			
			# Mean longitude for day of year:
			dlamm <- xlam + (n - 80.0)*(360.0/N)
			
			# Mean anomaly:
			anm <- dlamm - komega
			ranm <- anm*pir
			
			# True anomaly (uncorrected):
			ranv <- ranm + (2.0*ke - xec/4.0)*sin(ranm) +
			5.0/4.0*xee*sin(2.0*ranm) +
			13.0/12.0*xec*sin(3.0*ranm)
			anv <- ranv/pir
			
			# True longitude:
			my_tls <- anv + komega
			# 
			# if (my_tls < 0){
			# 	my_tls <- my_tls + 360
			# } else if (my_tls > 360) {
			# 	my_tls <- my_tls - 360
			# }
			# True anomaly:
			# my_nu <- my_tls - komega
			# if (my_nu < 0){
			# 	my_nu <- my_nu + 360
			# }
			
			my_tls[my_tls < 0]<-my_tls[my_tls < 0] + 360
			my_tls[my_tls > 360]<-my_tls[my_tls > 360] - 360
			my_nu <- my_tls - komega
			my_nu[my_nu < 0]<-my_nu[my_nu < 0]+ 360
			return (c(my_nu, my_tls))
		}
		
		
		# ************************************************************************
		# Name:     dcos
		# Inputs:   double (d), angle in degrees
		# Returns:  double, cosine of angle
		# Features: This function calculates the cosine of an angle (d) given
		#           in degrees.
		# Depends:  pir
		# Ref:      This script is based on the Javascript function written by
		#           C Johnson, Theoretical Physicist, Univ of Chicago
		#           - 'Equation of Time' URL: http://mb-soft.com/public3/equatime.html
		#           - Javascript URL: http://mb-soft.com/believe/txx/astro22.js
		# ************************************************************************
		dcos <- function(d) {
			cos(d*pir)
		}
		
		
		# ************************************************************************
		# Name:     dsin
		# Inputs:   double (d), angle in degrees
		# Returns:  double, sine of angle
		# Features: This function calculates the sine of an angle (d) given
		#           in degrees.
		# Depends:  pir
		# ************************************************************************
		dsin <- function(d) {
			sin(d*pir)
		}
		###########################################################################
		# 02. Define variables
		###########################################################################
		solar <- list()
		n<-doy
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 01. Calculate the number of days in yeark (kN), days
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		kN <- (julian_day(y + 1, 1, 1) - julian_day(y, 1, 1))
		solar$kN <- kN
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 02. Calculate heliocentric longitudes (nu and lambda), degrees
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		my_helio <- berger_tls(n, kN)
		nu <- my_helio[1]
		lam <- my_helio[2]
		solar$nu_deg <- nu
		solar$lambda_deg <- lam
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 03. Calculate distance factor (dr), unitless
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Berger et al. (1993)
		kee <- ke^2
		rho <- (1 - kee)/(1 + ke*dcos(nu))
		dr <- (1/rho)^2
		solar$rho <- rho
		solar$dr <- dr
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 04. Calculate the declination angle (delta), degrees
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Woolf (1968)
		delta <- asin(dsin(lam)*dsin(keps))
		delta <- delta/pir
		solar$delta_deg <- delta
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 05. Calculate variable substitutes (u and v), unitless
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		a<-dsin(delta)*dcos(lat)*dsin(slop)*dcos(asp)-dsin(delta)*dsin(lat)*dcos(slop)
		b<-dcos(delta)*dcos(lat)*dcos(slop)+dcos(delta)*dsin(lat)*dsin(slop)*dcos(asp)
		c<-dcos(delta)*dsin(slop)*dsin(asp)
		d<-b^2+c^2-a^2
		d[d<0]<-0
		sinfirst<-(a*c+b*sqrt(d))/(b^2+c^2)
			
		ru <- -1*a+c*sinfirst
		rv <- b
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 06. Calculate the sunset hour angle (hs), degrees
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Note: u/v equals tan(delta) * tan(lat)
		hs <- acos(-1.0*ru/rv)
		hs <- hs / pir
		hs[ru/rv >= 1.0]<-180 # Polar day (no sunset)
		hs[ru/rv <= -1.0]<-0 # Polar night (no sunrise)
		solar$hs_deg <- hs
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 07. Calculate daily extraterrestrial radiation (ra_d), J/m^2
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# ref: Eq. 1.10.3, Duffy & Beckman (1993)
		ra_d <- (86400/pi)*kGsc*dr*(ru*pir*hs + rv*dsin(hs))
		solar$ra_j.m2 <- ra_d
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 08. Calculate transmittivity (tau), unitless, and sunshine fraction
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# ref:  Eq. 11, Linacre (1968); Eq. 2, Allen (1996)
		# tau_o <- (kc + kd*sf)
		# tau <- tau_o*(1 + (2.67e-5)*elev)
		# solar$tau_o <- tau_o
		# solar$tau <- tau
		# Eq. 2, Allen (1996)
		tau_o = (kc + kd)*(1 + (2.67e-5)*elev)
		#general global AP radiation model Suehrcke, et al., 2013 
		
		tau<-tau_o*(0.1898+(1-0.1898)*sf^0.7410)
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 09. Calculate daily incoming radiation W/m^2
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		rad_in<-ra_d*tau/kSecInDay
		
		return(rad_in)
		
	}
	
	gc()
	
	###############################################################################################
	# 10. set the clusters for rad parallel computing
	###############################################################################################	
	cl <- getCluster()
	
	nodes <- length(cl)
	message('Using cluster with ', nodes, ' nodes')	
	bs <- blockSize(sf, minblocks=nodes*2)
	snow:::clusterExport(cl, c('y','calc_sw_in','doy','bs','sf','elev_lowres','lat_lr','terraines'),envir=environment()) 
	cl
	pb <- pbCreate(bs$n)
	pb <- txtProgressBar(min=1,max = max(bs$n,2), style = 1)
	###############################################################################################
	# 11. create the functions to send to the workers, split the data in chunks
	###############################################################################################	
	clFun_rad <- function(i) {
		sf_lr_block<- getValues(sf, bs$row[i], bs$nrows[i])
		elev_hr_block<-getValues(elev_lowres,bs$row[i], bs$nrows[i])
		lat_hr_block<-getValues(lat_lr,bs$row[i], bs$nrows[i])
		slop_hr_block<-getValues(terraines[[1]],bs$row[i], bs$nrows[i])
		asp_hr_block<-getValues(terraines[[2]],bs$row[i], bs$nrows[i])
		# do calculations
		result<-sapply(1:length(y), function(i) calc_sw_in(sf_lr_block[,i],lat_hr_block,elev_hr_block,slop_hr_block,asp_hr_block,y[i],doy[i]))
		return(result)
	}
			
	###############################################################################################
	# 12. send tasks to the nodes
	###############################################################################################
	for (i in 1:nodes) {
		snow:::sendCall(cl[[i]], clFun_rad, i, tag=i)
	}
	###############################################################################################
	# 13. write to the disk on the go, or save to the ram
	###############################################################################################
	
	if(!inmem){
		out<-writeStart(out,filename=paste0(ouputdir,"/",y[1],"_",y[length(y)],"_",'sw_in',".","nc"),format="CDF",overwrite=TRUE,varname="sw_in", varunit="W/m2",longname='shortwave radiation', xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(y[1]-1,"-",12)), ...)
		
		
	}else {
		matout <- matrix(ncol=nlayers(out), nrow=ncell(out))
		endind<-cumsum(bs$nrows*out@ncols)
		startind<-c(1,endind+1)    
	}
	###############################################################################################
	# 14. receive results from the nodes
	###############################################################################################	
	for (i in 1:bs$n) {
		
		d <- snow:::recvOneData(cl)
		# error?
		if (! d$value$success) {
			stop('error!! check the data...')
		}
		# which block is this?
		b <- d$value$tag
		# cat('received block: ',b,'\n'); flush.console();
		if (!inmem) {
			out <- writeValues(out,d$value$value, bs$row[b])
			
		} else {
			
			matout[startind[b]:endind[b],] <- d$value$value
		}
		
		# need to send more data?
		ni <- nodes + i
		if (ni <= bs$n) {
			snow:::sendCall(cl[[d$node]], clFun_rad, ni, tag=ni)
		}
		setTxtProgressBar(pb,i)
	}
	###############################################################################################
	# 15. close connection with the files, or assign valueas to the raster objects
	###############################################################################################
	
	if (!inmem) {
		out <- writeStop(out)
		
	} else {
		
		out<-setValues(out,matout)
		
	}
	close(pb)
	gc()
	
	return(out)
}
