#' readFluxdata
#'
#' Reads fluxnet data, calculates daily cuantities, including evapotranspiration. VPD, air temperature, and netrad use the daily daytime, albedo and snowmelt are still in experimental stage, the function takes only positive latent heat to compute actual evapotranspiration, the function works for the fluxnet tier 1-2 datasets, ameriflux and europa flux.
#' @param   filename: location of the csv file with subdaily data
#' @param   elev: elevation in m.a.s.l.
#' @import xts
#' @import data.table
#' @keywords fluxnet
#' @return xts dataframe with the following variables:
#' \itemize{
#'         \item \code{aet}: actual evapotranspiration (mm/day)
#'         \item \code{pet}: potential evapotranspiration (mm/day)
#'         \item \code{eeq}: equilibrium evapotranspiration (mm/day)
#'         \item \code{tc}: daily daytime air temperature (C)
#'         \item \code{SW_in}: daily incoming shortwave (W/m^2)
#'         \item \code{pn}: daily precipitation (mm/day)
#'         \item \code{sm}: daily soil moisture (volumetric %)
#'         \item \code{netr}: daily daytime net radiation (MJ/day)
#'         \item \code{VPD}: daily daytime vapour pressure deficit (Pa)
#'         \item \code{VPD_SA}: daily daytime vapour pressure deficit surface-air (Pa)(calculated using Tsurf, measured or inferred from LW)
#'         \item \code{CO2}: daily daytime CO2 concentration (ppm)
#'         \item \code{GPP}: daily GPP (gC/m^2/day)
#'         \item \code{snowmelt}: daily snowmelt (mm/day)---experimental
#'         \item \code{alb}: daily median albedo ---experimental
#'         \item \code{cond}: daily condensation (mm/day) ---experimental, negative LE migth mean refreezing
#' }
#' @export
#' @examples readFluxdata(filename,elev)
readFluxdata<-function(filename,elev=NULL,sensors_d=NULL, SWC='v/v'){
	# testing
	# filename=filenames.fluxnet.subset[36];elev=fluxnet_geo$elev_30m[35]; sensors_d=sensor_depths$`FR-LBr`; SWC='v/v'
	#filename<-filenames.flux.mts[63,1]
	# filename<-"X:/home/WORK/data_input/ameriflux_L4/AMF_L4_NT_US-CPk_HH_2009-2013.csv"
	# elev<-as.numeric(as.character(FLUXNET_2015@data$elv[22]))
	# elev<-918
	
	# fluxnet_data<-read.table(unz(filename.FLX.zip, filename),  header=T, quote="\"", sep=",",na.strings = -9999)
	###############################################################################################
	# 01. Read the data, assign time info, get the time interval, define the time interval in seconds t_conv_f
	###############################################################################################	
	# fluxnet_data<-read.table(filename,  header=T, quote="\"", sep=",",na.strings = -9999)
	fluxnet_data<-data.table::fread(filename,  header=T, quote="\"", sep=",",na.strings = c("NA",'-9999'),colClasses='numeric',integer64="character")
	# fluxnet_data[fluxnet_data==-9999]<-NA
	ind<-strptime(fluxnet_data$TIMESTAMP_START,format="%Y%m%d%H%M",tz="GMT")
	time.freq<-abs(as.numeric(ind[1]-ind[2], units = "hours"))
	t_conv_f<-3600*time.freq
	site<-do.call(rbind,strsplit(basename(filename),'_'))[,2]
	###############################################################################################
	# 02.define the constants
	###############################################################################################
	#
	kA <- 107           # constant for Rl (Monteith & Unsworth, 1990)
	# kalb_sw <- 0.17     # shortwave albedo (Federer, 1968)
	kalb_sw <- 0.17     # shortwave albedo (Federer, 1968)
	kalb_vis <- 0.03    # visible light albedo (Sellers, 1985)
	kb <- 0.20          # constant for Rl (Linacre, 1968; Kramer, 1957)
	kc <- 0.25          # constant for Rs (Linacre, 1968)
	kCw <- 1.05         # supply constant, mm/hr (Federer, 1982)
	kd <- 0.50          # constant for Rs (Linacre, 1968)
	kfFEC <- 2.04       # from-flux-to-energy, umol/J (Meek et al., 1984)
	kfus <- 334000      # latent heat of fusion J/Kg (Monteith & Unsworth, 1990)
	kG <- 9.80665       # gravitational acceleration, m/s^2 (Allen, 1973)
	kGsc <- 1360.8      # solar constant, W/m^2 (Kopp & Lean, 2011)
	kL <- 0.0065        # adiabatic lapse rate, K/m (Cavcar, 2000)
	kMa <- 0.028963     # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
	kMv <- 0.01802      # mol. weight of water vapor, kg/mol (Tsilingiris, 2008)
	kSecInDay <- 86400  # number of seconds in a day
	kp_snow<- 250       # snow density kg/m3 (Bonan et al., 2002)
	kPo <- 101325       # standard atmosphere, Pa (Allen, 1973)
	kR <- 8.31447       # universal gas constant, J/mol/K (Moldover et al., 1988)
	kTo <- 288.15       # base temperature, K (Berberan-Santos et al., 1997)
	kWm <- 150          # soil moisture capacity, mm (Cramer-Prentice, 1988)
	kw <- 0.26          # PET entrainment, (1+kw)*EET (Priestley-Taylor, 1972)
	pir <- pi/180       # pi in radians
	fluidity<-35187037  # fluidity as defined en Hillel (1998), mm^2 (density*gravity)/viscosity
	# Paleoclimate variables:
	ke <- 0.01670       # eccentricity of earth's orbit, 2000CE (Berger 1978)
	keps <- 23.44       # obliquity of earth's elliptic, 2000CE (Berger 1978)
	komega <- 283       # lon. of perihelion, degrees, 2000CE (Berger, 1978)
	ksig <- 5.670374e-8    # Stefan-Boltzmann constant https://physics.nist.gov/cgi-bin/cuu/Value?sigma
	###############################################################################################
	# 03.define the functions
	###############################################################################################
	
	# ************************************************************************
	# Name:     density_h2o
	# Inputs:   - double (tc), air temperature, degrees C
	#           - double (pa), atm pressure, Pa
	# Returns:  double, kg/m^3
	# Features: This function calculates the temperature and pressure
	#           dependent density of pure water
	# * Ref:    Chen, C.T., R.A. Fine, and F.J. Millero (1977), The equation
	#             of state of pure water determined from sound speeds, The
	#             Journal of Chemical Physics 66, 2142;
	#             doi:10.1063/1.434179
	# ************************************************************************
	density_h2o <- function(tc, pa) {
		# Calculate density of water at 1 atm, g/cm^3
		po <- 0.99983952 +
		(6.788260e-5)*tc +
		-(9.08659e-6)*tc*tc +
		(1.022130e-7)*tc*tc*tc +
		-(1.35439e-9)*tc*tc*tc*tc +
		(1.471150e-11)*tc*tc*tc*tc*tc +
		-(1.11663e-13)*tc*tc*tc*tc*tc*tc +
		(5.044070e-16)*tc*tc*tc*tc*tc*tc*tc +
		-(1.00659e-18)*tc*tc*tc*tc*tc*tc*tc*tc
		
		# Calculate the bulk modulus of water at 1 atm, atm
		ko <- 19652.17 +
		148.1830*tc +
		-2.29995*tc*tc +
		0.01281*tc*tc*tc +
		-(4.91564e-5)*tc*tc*tc*tc +
		(1.035530e-7)*tc*tc*tc*tc*tc
		
		# Calculate temperature-dependend coefficients
		ca <- 3.26138 +
		(5.223e-4)*tc +
		(1.324e-4)*tc*tc +
		-(7.655e-7)*tc*tc*tc +
		(8.584e-10)*tc*tc*tc*tc
		
		cb <- (7.2061e-5) +
		-(5.8948e-6)*tc +
		(8.69900e-8)*tc*tc +
		-(1.0100e-9)*tc*tc*tc +
		(4.3220e-12)*tc*tc*tc*tc
		
		# Convert pressure to bar (1 bar = 100000 Pa)
		pbar <- (1e-5)*pa
		
		pw <- (1e3)*po*(ko + ca*pbar + cb*pbar^2)/(ko + ca*pbar + cb*pbar^2 - pbar)
		return(pw)
	}
	
	
	# ************************************************************************
	# Name:     elv2pres
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
	elv2pres <- function(z) {
		kPo*(1 - kL*z/kTo)^(kG*kMa/(kR*kL))
	}
	# ************************************************************************
	# Name:     enthalpy_vap
	# Inputs:   double (tc), air temperature, degrees C
	# Returns:  double, J/kg
	# Features: This function calculates the temperature-dependent enthalpy
	#           of vaporization (latent heat of vaporization)
	# Ref:      Eq. 8, Henderson-Sellers (1984), A new formula for latent heat
	#             of vaporization of water as a function of temperature, Quarterly
	#             Journal of the Royal Meteorological Society, vol. 110, pp. 1186--
	#             1190.
	# ************************************************************************
	enthalpy_vap <- function(tc) {
		1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))^2
	}
	# ************************************************************************
	# Name:     specific_heat
	# Inputs:   double (tc), air temperature, degrees C
	# Returns:  double, specific heat of moist air, J/kg/K
	# Features: This function calculates the spefic heat of moist air
	# Ref:      Tsilingris (2008), Thermophysical and transport properties of
	#           humid air at temperature range between 0 and 100 °C, Energy
	#           Conversion and Management, vol. 49, pp. 1098--1110.
	# ************************************************************************
	specific_heat <- function(tc) {
		
		tc<-ifelse(tc < 0,0,tc)
		cp <- 1.0045714270 +
		(2.050632750e-3)*tc -
		(1.631537093e-4)*tc*tc +
		(6.212300300e-6)*tc*tc*tc -
		(8.830478888e-8)*tc*tc*tc*tc +
		(5.071307038e-10)*tc*tc*tc*tc*tc
		cp <- (1e3)*cp
		
		return(cp)
	}
	
	# ************************************************************************
	# Name:     psychro
	# Inputs:   - double (tc), air temperature, degrees C
	#           - double (pa), atm pressure, Pa
	# Returns:  double, Pa/K
	# Features: This function calculates the temperature and pressure
	#           dependent psychrometric constant
	# Depends:  - enthalpy_vap
	#           - specific_heat
	# Ref:      Allen, R.G., L.S. Pereira, D. Raes, M. Smith (1998),
	#           'Meteorological data,' Crop evapotranspiration - Guidelines
	#           for computing crop water requirements - FAO Irrigation and
	#           drainage paper 56, Food and Agriculture Organization of the
	#           United Nations, Available:
	#           http://www.fao.org/docrep/x0490e/x0490e07.htm
	# ************************************************************************
	psychro <- function(tc, pa) {
		# Calculate the specific heat capacity of water, J/kg/K
		cp <- specific_heat(tc)
		
		# Calculate latent heat of vaporization, J/kg
		lv <- enthalpy_vap(tc)
		
		# Calculate psychrometric constant, Pa/K
		return(cp*kMa*pa/(kMv*lv))
	}
	# ************************************************************************
	# Name:     sat_slope
	# Inputs:   double (tc), degrees C
	# Returns:  double, Pa/K
	# Features: This function calculates the temperature-dependent slope of
	#           the saturation pressure temperature curve using the
	#           methodology presented in the eMast energy.cpp script
	# Ref:      - Eq. 6, Prentice et al. (1993);
	#           - Eq. 13, Allen et al. (1998)
	# ************************************************************************
	sat_slope <- function(tc) {
		(17.269)*(237.3)*(610.78)*exp(17.269*tc/(237.3 + tc))/(237.3 + tc)^2
	}
	
	# ************************************************************************
	# Name:     sat_vap
	# Inputs:   double (tc), degrees C
	# Returns:  double, Pa
	# Features: This function calculates the 
	#           the water vapour pressure at saturation
	#           
	# Ref:      
	#           - Eq. 13, Allen et al. (1998)
	# ************************************************************************
	sat_vap <- function(tc) {
		0.6108*1000*exp((17.27*tc)/(tc+237.3))
	}
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
	###############################################################################################
	# 03.define the daily variables
	###############################################################################################
	
	
	
	# get daily SW in w/m2 
	if(!is.null(fluxnet_data$SW_IN_F)){
		SW_IN_F<-xts(fluxnet_data$SW_IN_F,ind)
	}else{
		if(!is.null(fluxnet_data$SW_IN)){
			SW_IN_F<-xts(fluxnet_data$SW_IN,ind)
		}else{
			SW_IN_F<-xts(rep(NA,length(ind)),ind)
		}
		
	}
	SW_IN_F[SW_IN_F<0]<-NA
	SW_IN_mean<-apply.daily(SW_IN_F,mean,na.rm=TRUE)
	
	# get daily daytime mean temperature
	if(!is.null(fluxnet_data$TA_F)){
		TA_F<-xts(fluxnet_data$TA_F,ind)
	}else{
		if(!is.null(fluxnet_data$TA)){
			TA_F<-xts(fluxnet_data$TA,ind)
		}else{
			TA_F<-xts(rep(NA,length(ind)),ind)
		}
		
	}
	#temperature where LE is negative tcond
	t_cond<-TA_F
	#daytime temperature tc
	TA_F[SW_IN_F<=0]<-NA
	tc<-TA_F
	#mean daytime air temperature TA_F
	TA_F<-apply.daily(TA_F,mean,na.rm=TRUE)
	
	# get daily vpd Pa
		
	if(!is.null(fluxnet_data$VPD_F)){
		vpd<-xts(fluxnet_data$VPD_F*100,ind)
	}else{
		if(!is.null(fluxnet_data$VPD)){
			vpd<-xts(fluxnet_data$VPD*100,ind)
		}else{
			if(!is.null(fluxnet_data$VPD_PI)){
				vpd<-xts(fluxnet_data$VPD_PI*100,ind)
			}else{
				if(!is.null(fluxnet_data$RH)){
					es<-sat_vap(tc)
					vpd<-es*(1-(fluxnet_data$RH/100))
					vpd<-xts(vpd,ind)
				}else{
					vpd<-xts(rep(NA,length(ind)),ind)
				}
			}
		}
		
	}
	
	
	
	
	
	vpd[SW_IN_F<=0]<-NA
	vpd<-apply.daily(vpd,mean,na.rm=TRUE)
	###############################################################################################
	# get daily vpd leaf - air
	###############################################################################################
	if(!is.null(fluxnet_data$T_CANOPY)){
		Tsurf<-xts(fluxnet_data$T_CANOPY,ind)
		Tsurf[SW_IN_F<=0]<-NA
		Tsurf<-apply.daily(Tsurf,mean,na.rm=TRUE)
		#Calculate surfate temperature assuming emissivity 0.97
	}else{
		if(!is.null(fluxnet_data$LW_OUT)){
			if(!is.null(fluxnet_data$LW_IN_F)){
				Tsurf<-(((fluxnet_data$LW_OUT-(1-0.97)*fluxnet_data$LW_IN_F)/(ksig*0.97))^(1/4))-273.15
			}else if(!is.null(fluxnet_data$LW_IN)){
				Tsurf<-(((fluxnet_data$LW_OUT-(1-0.97)*fluxnet_data$LW_IN)/(ksig*0.97))^(1/4))-273.15
			}
			Tsurf[SW_IN_F<=0]<-NA
			# subset only measured data
			if(!is.null(fluxnet_data$LW_IN_F_QC)){
				Tsurf[fluxnet_data$LW_IN_F_QC>1]<-NA
			}
			Tsurf<-xts(Tsurf,ind)
			Tsurf<-apply.daily(Tsurf,mean,na.rm=TRUE)
			
		}else{
			
			Tsurf<-xts(rep(NA,length(vpd)),time(vpd))
		}
	}	
	
	
	# air saturation vapour pressure [Pa]
	es<-sat_vap(TA_F)
	#actual water vapour pressure 
	ea<-es-vpd
	# surface saturation vapour pressure [Pa]
	es_L<-0.6108*1000*exp((17.27*Tsurf)/(Tsurf+237.3))
	#vpd surface air
	vpd_sa<-es_L-ea
	vpd_sa[vpd_sa<0]<-0
	
	
	# get daily co2 ppm
	if(!is.null(fluxnet_data$CO2_F_MDS)){
		co2<-xts(fluxnet_data$CO2_F_MDS,ind)
	}else{
		if(!is.null(fluxnet_data$CO2)){
			co2<-xts(fluxnet_data$CO2,ind)
		}else{
			if(!is.null(fluxnet_data$CO2_1)){
				co2_1<-xts(fluxnet_data$CO2_1,ind)
			}else{
				co2_1<-xts(rep(NA,length(ind)),ind)
			}
			if(!is.null(fluxnet_data$CO2_2)){
				co2_2<-xts(fluxnet_data$CO2_2,ind)
			}else{
				co2_2<-xts(rep(NA,length(ind)),ind)
			}
			co2<-cbind(co2_1,co2_2)
			# assuming 0.5 meter depth
			co2<-rowMeans(co2,na.rm=TRUE)
			co2<-xts(co2,ind)
		}
		
	}
	co2[SW_IN_F<=0]<-NA
	co2<-apply.daily(co2,mean,na.rm=TRUE)
	# get daily ppfd mol/m2
		
	if(!is.null(fluxnet_data$PPFD_IN)){
		ppfd<-xts(fluxnet_data$PPFD_IN*t_conv_f*1e-6,ind)
		ppfd<-apply.daily(ppfd,sum,na.rm=TRUE)
		ppfd[ppfd==0]<-NA
	}else{
		ppfd<-(1e-6)*(kfFEC*(1 - kalb_vis)*(SW_IN_mean*86400))
		
	}
		
	# get daily gpp from umolCO2/m2/s to gC/m2/day
	if(!is.null(fluxnet_data$GPP_NT_VUT_REF)){
		gpp<-xts(fluxnet_data$GPP_NT_VUT_REF*t_conv_f*1e-6,ind)
	}else{
		if(!is.null(fluxnet_data$GPP)){
			gpp<-xts(fluxnet_data$GPP*t_conv_f*1e-6,ind)
		}else if (is.null(fluxnet_data$GPP)){
			gpp<-xts(fluxnet_data$GPP_st_ANN*t_conv_f*1e-6,ind)
		}else{
			gpp<-xts(rep(NA,length(ind)),ind)
		}
		
	}
	
	gpp[gpp<0.0]<-0.0
	gpp<-apply.daily(gpp,sum,na.rm=TRUE)
	gpp<-gpp*12.0107
	gpp[gpp==0.0]<-NA
		
	# Atmospheric pressure, Pa
	if(!is.null(fluxnet_data$PA)){
		patm <-1000*fluxnet_data$PA
		if(!is.null(elev)){
			patm[is.na(patm)]<-elv2pres(elev)
		}
		
		patm<-apply.daily(xts(patm,ind),mean,na.rm=TRUE)
		
	}else if (is.null(fluxnet_data$PA) & !is.null(elev)){
		patm <- elv2pres(elev)
	}else{
		patm<-xts(rep(NA,length(gpp)),time(gpp))
	}
	# 3.1. Calculate water-to-energy conversion (econ) daytime, m^3/J
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Slope of saturation vap press temp curve, Pa/K
	s <- sat_slope(TA_F)
	# Enthalpy of vaporization, J/kg
	lv <- enthalpy_vap(TA_F)
	# Density of water, kg/m^3
	pw <- density_h2o(TA_F, patm)
	# Psychrometric constant, Pa/K
	gam <- psychro(TA_F, patm)
	# get factor to convert energy to water flux J/m^3
	econ <- s/(lv*pw*(s + gam))
	
	
	
	# ET
	
	# get actual evapotranspiration and potential (mm/day) *1000 from m3 to litres, t_conv_f to seconds to half hour
	if(!is.null(fluxnet_data$LE_SSITC_TEST)& is.null(fluxnet_data$LE_F_MDS_QC)){
		LE_QF<-fluxnet_data$LE_SSITC_TEST
	}else if(is.null(fluxnet_data$LE_SSITC_TEST)& !is.null(fluxnet_data$LE_F_MDS_QC)){
		LE_QF<-fluxnet_data$LE_F_MDS_QC
	}else if(is.null(fluxnet_data$LE_F_MDS_QC)){
		LE_QF<-fluxnet_data$LE_fqc
	}
	
	if(!is.null(fluxnet_data$LE_CORR)){
		LE<-xts(fluxnet_data$LE_CORR,ind)
	}else{
		if(!is.null(fluxnet_data$LE_CORR)){
			LE<-xts(fluxnet_data$LE_F_MDS,ind)
		}else{
			if(!is.null(fluxnet_data$LE)){
				LE<-xts(fluxnet_data$LE,ind)
			}else{
				LE<-xts(rep(NA,length(ind)),ind)
			}
		}
		
	}
	# subset only measured data
	LE[LE_QF>1]<-NA
	#aggregate daily negative LE	
	cond_LE<-LE
	t_cond[cond_LE>0]<-NA
	cond_LE[cond_LE>0]<-NA
	cond_LE<-apply.daily(-1*cond_LE*t_conv_f,sum,na.rm=TRUE)
	cond_LE[cond_LE<=0]<-NA
	#mean air temperature when LE is negative
	t_cond<-apply.daily(t_cond,mean,na.rm=TRUE)
	#aggregate daily positive LE
	LE[LE<0]<-NA
	LE<-apply.daily(LE*t_conv_f,sum,na.rm=TRUE)
	LE[LE<=0]<-NA
	# get daily snowmelt
	if(!is.null(fluxnet_data$SW_OUT)){
		alb_day<-apply.daily(xts(fluxnet_data$SW_OUT,ind),sum,na.rm=TRUE)/apply.daily(SW_IN_F,sum,na.rm=T)
		snowmelt<-LE
		snowmelt[alb_day<0.4]<-0
		LE[alb_day>=0.4]<-NA
		aet<-LE*(1/lv)*(1/pw)*1000.0
		snowmelt<-snowmelt*1000/(pw*kfus)
		
		
	}else{
		
		aet<-LE*(1/lv)*(1/pw)*1000.0
		alb_day<-xts(rep(NA,length(aet)),time(aet))
		snowmelt<-xts(rep(NA,length(aet)),time(aet))
	}
	#smooth albedo to monthly
	alb_day<-apply.monthly(alb_day,median)
	alb_day<-na.approx(alb_day,xout=time(aet),na.rm=F)
	# 3.1. Calculate water-to-energy conversion (econ) migthtime, m^3/J
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Slope of saturation vap press temp curve, Pa/K
	s_n <- sat_slope(t_cond)
	# Enthalpy of vaporization, J/kg
	lv_n <- enthalpy_vap(t_cond)
	# Density of water, kg/m^3
	pw_n <- density_h2o(t_cond, patm)
	# daily condensation (mm)
	cond<-cond_LE*(1/lv_n)*(1/pw_n)*1000.0
	
	# get net radiation J/m2
	if(!is.null(fluxnet_data$NETRAD)){
		NETRAD<-xts(fluxnet_data$NETRAD,ind)
		
	}else{
		NETRAD<-xts(rep(NA,length(ind)),ind)
	}
	NETRAD[NETRAD<0]<-NA
	NETRAD<-apply.daily(NETRAD*t_conv_f,sum,na.rm=TRUE)
	NETRAD[NETRAD==0]<-NA
	# get potential evapotrasnpiration
	if(!is.null(fluxnet_data$NETRAD)&!is.null(fluxnet_data$G)){
		AE<-fluxnet_data$NETRAD-fluxnet_data$G
		AE[AE<0]<-NA
		AE<-xts(AE,ind)
	}else if (!is.null(fluxnet_data$NETRAD)& is.null(fluxnet_data$G)){
		
		AE<-fluxnet_data$NETRAD
		AE[AE<0]<-NA
		AE<-xts(AE,ind)
	}else{
		AE<-xts(rep(NA,length(ind)),ind)
	}
	
	AE<-apply.daily(AE*t_conv_f,sum,na.rm=TRUE)
	AE[AE==0]<-NA
	eeq<-AE*econ*1000.0
	pet<-eeq*(1+0.26)
	# !!!!!!!!Error not energy balance clossure, aet>pet, not posible
	pet[aet>pet]<-aet[aet>pet]
	#################################################################################
	######### get soil moisture ######################################################
	#################################################################################
	if(length(grep("^SWC", names(fluxnet_data), value = TRUE))!=0){
		nameswc<-grep("^SWC", names(fluxnet_data), value = TRUE)
		nameswc_qc<-grep("^SWC.*.QC", names(fluxnet_data), value = TRUE)
		nameswc<-nameswc[!(nameswc %in% nameswc_qc)]
		
		sm<-fluxnet_data[,..nameswc]
		if(SWC=='v/v'){
			if(length(nameswc)>1){
				if(is.null(sensors_d)|is.na(sensors_d) | (length(sensors_d)!=length(nameswc))){
					sm<-rowMeans(sm,na.rm=F)
				}else{
					upper_bound=-1*sensors_d
					lowerbound<-c(0,upper_bound[1:(length(upper_bound)-1)])
					depth_inter<-upper_bound-lowerbound
					weight<-depth_inter/sum(depth_inter)
					sm<- apply(sm, 1, weighted.mean, weight)
				}
							
			}
		}else if (SWC=='mm'){
			if(is.null(sensors_d)){
				warning('Impossible to calculate SWC in mm, sensor depth(s) not provided, returning NAs')
				sm<-rep(NA,length(ind))
			}else{
				#m3/m3 * m *1000 = litre/m2,(sm in %, so factor=10)
				upper_bound=-1*sensors_d*10
				lowerbound<-c(0,upper_bound[1:(length(upper_bound)-1)])
				depth_inter<-upper_bound-lowerbound
				sm=t(t(sm) * depth_inter)
				sm<-rowSums(sm,na.rm=F)
			}
			
			
		}
		
				
	}else{
		sm<-rep(NA,length(ind))
	}

	sm<-xts(sm,ind)
	sm<-apply.daily(sm,median,na.rm=TRUE)
	
	

	
	# storage variation
	# for(i in 1:length(sm)){
	# 	j<-i+1
	# 	if(i==1){
	# 		delta_sm<-as.numeric(sm)[j]-as.numeric(sm)[i]
	# 	}else{
	# 		delta_sm<-c(delta_sm,as.numeric(sm)[j]-as.numeric(sm)[i])	
	# 	}
	# }
	# delta_sm<-xts(delta_sm,time(sm))
	
	# delta_sm[delta_sm<0]<-0
	# precipitation
	if(!is.null(fluxnet_data$P_F)){
		P<-xts(fluxnet_data$P_F,ind)
	}else{
		if(!is.null(fluxnet_data$P)){
			P<-xts(fluxnet_data$P,ind)
		}else{
			P<-xts(rep(NA,length(ind)),ind)
		}
	}
	
	P<-apply.daily(P,sum)
	# runoff
	# RO<-P
	# RO[delta_sm>0]<-RO[delta_sm>0]-delta_sm[delta_sm>0]
	# RO[delta_sm<=0]<-0
	# RO[is.na(delta_sm)]<-NA
	
	# bal<-P-aet-RO-delta_sm
	
	
	
	
	rm(fluxnet_data)
	gc()
	result<-merge.xts(aet,pet,eeq,TA_F,SW_IN_mean,P,sm,NETRAD/1e6,vpd,vpd_sa,co2,ppfd,gpp,snowmelt,alb_day,cond,Tsurf)
	index(result)<-as.Date(index(result))
	names(result)<-c("aet","pet","eeq","tc","SW_in","pn","sm","netr","VPD",'VPD_SA',"CO2","PPFD","GPP",'snowmelt','alb','cond','Tsurf')
	
	return(result)
	
}
