#' downscaleSolar
#'
#' Downscale solar radiation grid applying corrections for terrain features
#' @param  elev_hres, elevation, a high resolution Raster layer object
#' @param  rad_lowres, solar radiation (w/m2), the low resolution Raster* object with z time dimension
#' @param  ouputdir, directory where to save the files
#' @return  solar radiation (w/m2), Raster* object with the resolution of elev_hres and the same z time dimension as rad_lowres
#' @import raster  
#' @keywords Downscale, solar radiation, terrain
#' @export
#' @examples
#' # *optional run beginCluster() first, for parallel computing
#' downscaleSolar()
downscaleSolar<-function(elev_hres,elev_lowres,rad_lowres,ouputdir=getwd(),inmem=FALSE, ...){
	
	###########################################################################
	# 00. Check if parallel computation is required by the user and if the dimensions of the raster objects match
	###########################################################################
	on.exit(endCluster())
	clcheck<-try(getCluster(), silent=TRUE)
	if(class(clcheck)=="try-error"){
		# If no cluster is initialized, assume only one core will do the calculations, beginCluster(1) saved me the time of coding serial versions of the functions
		beginCluster(1,'SOCK')
		message('Only using one core, use first beginCluster() if you want to run splash in parallel!!')
		
	}
	rasterOptions(tolerance = 0.5)
	compareRaster(rad_lowres, elev_lowres,extent=TRUE, crs=TRUE, res=TRUE, orig=FALSE,rotation=FALSE, values=FALSE, stopiffalse=FALSE, showwarning=TRUE)	
	
	
	###############################################################################################
	# 00. create array for results, get time info
	###############################################################################################
	tmpdir<-dirname(rasterTmpFile())
	setwd(tmpdir)
	y<-as.numeric(format(getZ(rad_lowres),'%Y'))
	doy<-as.numeric(format(getZ(rad_lowres),'%j'))
	ny <-length(doy)
	ind<-getZ(rad_lowres)
	# sf
	sf<-brick(nrows=nrow(elev_lowres), ncols=ncol(elev_lowres), crs=crs(elev_lowres), nl=ny)
	extent(sf)<-extent(elev_lowres)
	sf<-setZ(sf,ind)
	# radiation hres
	out<-brick(nrows=nrow(elev_hres), ncols=ncol(elev_hres), crs=crs(elev_hres), nl=ny)
	extent(out)<-extent(elev_hres)
	out<-setZ(out,ind)
	# 1.2 get latitudes
	lat_lr<-elev_lowres*0
	lat.data<-rasterToPoints(elev_lowres)
	lat_lr[!is.na(lat_lr)]<-lat.data[,2]
	rm(lat.data)
	#### function to get the latitudes from big rasters i.e 1km res global extent
	getlatitude <- function(x, filename, ...) {
		##create array for the results
		out <-raster(x)
		# get the index of the blocks, every block has n rows, bigger the minblocks, smaller the chunk of rows
		bs <- blockSize(x, minblocks=200)
		pb <- pbCreate(bs$n)
		pb <- txtProgressBar(min=1,max = bs$n, style = 1)
		# start writing the outputfile
		out <- writeStart(out, filename, overwrite=TRUE)
		# rsqv<-function(x,y){summary(lm(y~x,na.action=))$r.squared}
		for (i in 1:bs$n) {
			# i=178
			# xmat <- getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
			# ymat <- getValues(y, row=bs$row[i], nrows=bs$nrows[i] )
			xncells<-cellFromRow(x,bs$row[i]:(bs$row[i]+ bs$nrows[i]-1))
			xmat<-getValues(x,bs$row[i], bs$nrows[i])
			xydata<-xyFromCell(x, xncells)
			
			# write the chunk of results, bs$row[i] is putting the results in the correct rows
			out <- writeValues(out, xydata[,2], bs$row[i])
			setTxtProgressBar(pb,i)
		}
		out <- writeStop(out)
		close(pb)
		return(out)
	}
	# 1.2.2 get latitude from the high res dem
	lat_hr<-getlatitude(elev_hres, filename='testlat.grd')
	gc()
	# 1.3  calculate slope and aspect
	rasterOptions(maxmemory=1e9, timer=TRUE, tmptime = 24, chunksize = 1e9,progress='text', overwrite=TRUE,tolerance=0.5,todisk=FALSE)
	terraines<-terrain(elev_hres, opt=c('slope', 'aspect'), unit='degrees')
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
	# cat('fixing nas terrain')
	# terraines<-stack(lapply(X=n, FUN=fillna, x=terraines),quick=TRUE)
	#setwd(ouputdir)
	# rasterOptions(maxmemory=1e8, timer=TRUE, tmptime = 24, chunksize = 1e8, overwrite=TRUE,tolerance=0.5,todisk=FALSE)
	###########################################################################
	# 02. Define functions
	###########################################################################
	calc_sf<-function(sw_in,lat,elev,y,doy){
		
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
		# 05. Calculate variable substitutes (u and v), unitless, assume flat at low spatial resolution
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		ru<- dsin(delta)*dsin(lat)
		rv<- dcos(delta)*dcos(lat)
				
		solar$ru <- ru
		solar$rv <- rv
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
		# 08. Calculate daily incoming radiation J/m^2
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		rad_in<-sw_in*kSecInDay
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 08. Calculate transmittivity (tau), unitless, and sunshine fraction
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# ref:  Eq. 11, Linacre (1968); Eq. 2, Allen (1996)
		# tau_o <- (kc + kd*sf)
		# tau <- tau_o*(1 + (2.67e-5)*elv)
		# solar$tau_o <- tau_o
		# solar$tau <- tau
		# to avoid NA's at polar nigths, tau only defined by elevation, assume clear sky
		tau_o = (kc + kd)*(1 + (2.67e-5)*elev)
		
		tau<-ifelse(rad_in>ra_d | ra_d <= 0,tau_o,rad_in/ra_d)
		#AP sf
		#sf<-((tau/(1 + (2.67e-5)*elev))-kc)/kd
		#general global AP radiation model Suehrcke, et al., 2013 
		#sf = pow(((tau-tau_o*0.1898)/(tau_o*(1-0.1898))),(1/0.7410));
		sf= ((tau-tau_o*0.1898)/(tau_o*(1-0.1898)))^(1/0.7410)
		sf = ifelse(is.na(sf),0,ifelse(sf>1,1,sf))
		return(sf)
		
	}
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
		tau_o <- (kc + kd*sf)
		tau <- tau_o*(1 + (2.67e-5)*elev)
		solar$tau_o <- tau_o
		solar$tau <- tau
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 09. Calculate daily incoming radiation W/m^2
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		rad_in<-ra_d*tau/kSecInDay
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 08. Calculate transmittivity (tau), unitless, and sunshine fraction
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# ref:  Eq. 11, Linacre (1968); Eq. 2, Allen (1996)
		# tau_o <- (kc + kd*sf)
		# tau <- tau_o*(1 + (2.67e-5)*elv)
		# solar$tau_o <- tau_o
		# solar$tau <- tau
		# to avoid NA's at polar nigths, tau only defined by elevation, assume clear sky
		# tau<-ifelse(rad_in>ra_d | ra_d <= 0,(kc + kd)*(1 + (2.67e-5)*elev),rad_in/ra_d)
		# sf<-((tau/(1 + (2.67e-5)*elev))-kc)/kd
		# # kc and kd are not working for the whole world some pixels give these errors
		# sf[sf>1]<-1.0
		# sf[sf<0]<-0.0
		return(rad_in)
		
	}
	###############################################################################################
	# 03. set the clusters for parallel computing sf
	###############################################################################################	
	cl <- getCluster()
	on.exit( returnCluster() )
	nodes <- length(cl)
	message('getting sf with ', nodes, ' nodes')
	bs <- blockSize(rad_lowres, minblocks=nodes)
	parallel:::clusterExport(cl, c('y','calc_sf','doy','bs','rad_lowres','elev_lowres','lat_lr'),envir=environment()) 
	pb <- pbCreate(bs$n)
	pb <- txtProgressBar(min=1,max = bs$n, style = 1)
	###############################################################################################
	# 04. create the functions to send to the workers, split the data in chunks
	###############################################################################################	
	clFun_sf <- function(i) {
		rad_lr_block<- getValues(rad_lowres, bs$row[i], bs$nrows[i])
		elev_lowres_block<-getValues(elev_lowres,bs$row[i], bs$nrows[i])
		lat_lr_block<-getValues(lat_lr,bs$row[i], bs$nrows[i])
		# do calculations
		result<-sapply(1:length(y), function(i) calc_sf(rad_lr_block[,i],lat_lr_block,elev_lowres_block,y[i],doy[i]))
		return(result)
	}
	
	###############################################################################################
	# 05. send tasks to the nodes
	###############################################################################################
	for (i in 1:nodes) {
		parallel:::sendCall(cl[[i]], clFun_sf, i, tag=i)
	}
	###############################################################################################
	# 06. write to the disk on the go, or save to the ram
	###############################################################################################
	
	if(!inmem){
		sf<-writeStart(sf,filename = 'sf_lr.grd')
				
	}else {
		matout <- matrix(ncol=nlayers(sf), nrow=ncell(sf))
		endind<-cumsum(bs$nrows*out@ncols)
		startind<-c(1,endind+1)    
	}
	###############################################################################################
	# 07. receive results from the nodes
	###############################################################################################	
	for (i in 1:bs$n) {
		
		d <- parallel:::recvOneData(cl)
		# error?
		if (! d$value$success) {
			stop('error!! check the data...')
		}
		# which block is this?
		b <- d$value$tag
		# cat('received block: ',b,'\n'); flush.console();
		if (!inmem) {
			sf <- writeValues(sf,d$value$value, bs$row[b])
			
		} else {
			
			matout[startind[b]:endind[b],] <- d$value$value
		}
		
		# need to send more data?
		ni <- nodes + i
		if (ni <= bs$n) {
			parallel:::sendCall(cl[[d$node]], clFun_sf, ni, tag=ni)
		}
		setTxtProgressBar(pb,i)
	}
	###############################################################################################
	# 08. close connection with the files, or assign valueas to the raster objects
	###############################################################################################
	
	if (!inmem) {
		sf <- writeStop(sf)
		
	} else {
		
		sf<-setValues(sf,matout)
		
	}
	close(pb)
	gc()
	# 1.3  fillnas sf
	# n <- seq(nlayers(sf)) 
	# sf<-stack(lapply(X=n, FUN=fillna, x=sf),quick=TRUE)
	
	###############################################################################################
	# 09. project low res raster to hd resolution
	###############################################################################################
	cat('projecting sf')
	#sf_hres<-projectRaster(sf,elev_hres,filename="sf_hres.grd")
	sf_hres<-crop(sf,elev_hres)
	sf_hres<-projectRaster(sf_hres,elev_hres,filename="sf_hres.grd",overwrite=TRUE)
	# fx<-nrow(elev_hres)/nrow(sf_hres)
	# fy<-ncol(elev_hres)/ncol(sf_hres)
	# sf_hres<-disaggregate(sf_hres,c(fy,fx),filename="sf_hres.grd",overwrite=TRUE)
	# gc()
	
	###############################################################################################
	# 10. set the clusters for rad parallel computing
	###############################################################################################	
		
	bs <- blockSize(sf_hres, minblocks=nodes*10)
	message('computing sw with ', nodes, ' nodes')
	parallel:::clusterExport(cl, c('y','calc_sw_in','doy','bs','sf_hres','elev_hres','lat_hr','terraines'),envir=environment()) 
	pb <- pbCreate(bs$n)
	pb <- txtProgressBar(min=1,max = bs$n, style = 1)
	###############################################################################################
	# 11. create the functions to send to the workers, split the data in chunks
	###############################################################################################	
	clFun_rad <- function(i) {
		sf_lr_block<- getValues(sf_hres, bs$row[i], bs$nrows[i])
		elev_hr_block<-getValues(elev_hres,bs$row[i], bs$nrows[i])
		lat_hr_block<-getValues(lat_hr,bs$row[i], bs$nrows[i])
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
		parallel:::sendCall(cl[[i]], clFun_rad, i, tag=i)
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
		
		d <-parallel:::recvOneData(cl)
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
			parallel:::sendCall(cl[[d$node]], clFun_rad, ni, tag=ni)
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
