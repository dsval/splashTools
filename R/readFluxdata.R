#' readFluxdata
#'
#' Reads fluxnet data, calculates daily cuantities, including evapotranspiration. VPD, air temperature, and netrad use the daily daytime, albedo and snowmelt are still in experimental stage, the function takes only positive latent heat to compute actual evapotrasnpiration, the function works for the fluxnet tier 1-2 datasets, ameriflux and europa flux.
#' @param   filename: location of the csv file with subdaily data
#' @param   elev: elevation in m.a.s.l.
#' @import xts 
#' @keywords splash
#' @return xts dataframe with the following variables:
#' \itemize{
#'         \item \code{aet}: actual evapotranspiration (mm/day)
#'         \item \code{pet}: potential evapotranspiration (mm/day)
#'         \item \code{eeq}: Equilibrium evapotranspiration (mm/day)
#'         \item \code{tc}: daily daytime air temperature (C)
#'         \item \code{SW_in}: daily incoming shortwave (W/m^2)
#'         \item \code{pn}: daily precipitation (mm/day)
#'         \item \code{sm}: daily soil moisture (volumetric %)
#'         \item \code{netr}: daily daytime net radiation (MJ/day)
#'         \item \code{VPD}: daily daytime vapour pressure deficit (Pa)
#'         \item \code{CO2}: daily daytime CO2 concentration (ppm)
#'         \item \code{GPP}: daily GPP (gC/m^2/day)
#'         \item \code{snowmelt}: daily snowmelt (mm/day)---experimental
#'         \item \code{alb}: daily median albedo ---experimental
#' }
#' @export
#' @examples readFluxdata(filename,elev)
readFluxdata<-function(filename,elev){
	# testing
	# filename<-filenames.fluxnet[22]
	# filename.FLX.zip<-filenames.fluxnet[8]
	# elev<-as.numeric(as.character(FLUXNET_2015@data$elv[22]))
	# elev<-as.numeric(FLUXNET_mountain@data$Elev[10])
	
	# fluxnet_data<-read.table(unz(filename.FLX.zip, filename),  header=T, quote="\"", sep=",",na.strings = -9999)
	###############################################################################################
	# 01. Read the data, assign time info, get the time interval, define the time interval in seconds t_conv_f
	###############################################################################################	
	fluxnet_data<-read.table(filename,  header=T, quote="\"", sep=",",na.strings = -9999)
	# fluxnet_data[fluxnet_data==-9999]<-NA
	ind<-strptime(fluxnet_data$TIMESTAMP_START,format="%Y%m%d%H%M",tz="GMT")
	time.freq<-abs(as.numeric(ind[1]-ind[2], units = "hours"))
	t_conv_f<-3600*time.freq
	
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
	tc<-TA_F
	TA_F[SW_IN_F<=0]<-NA
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
	}else{
		ppfd<-(1e-6)*(kfFEC*(1 - kalb_vis)*(SW_IN_mean*86400))
		
	}
		
	# get daily gpp gC/m2
	if(!is.null(fluxnet_data$GPP_NT_VUT_REF)){
		gpp<-xts(fluxnet_data$GPP_NT_VUT_REF*t_conv_f*1e-6,ind)
	}else{
		if(!is.null(fluxnet_data$GPP)){
			gpp<-xts(fluxnet_data$GPP*t_conv_f*1e-6,ind)
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
		patm[is.na(patm)]<-elv2pres(elev)
	}else{
		patm <- elv2pres(elev)
	}
	# 3.1. Calculate water-to-energy conversion (econ), m^3/J
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Slope of saturation vap press temp curve, Pa/K
	s <- sat_slope(tc)
	# Enthalpy of vaporization, J/kg
	lv <- enthalpy_vap(tc)
	# Density of water, kg/m^3
	pw <- density_h2o(tc, patm)
	# Psychrometric constant, Pa/K
	gam <- psychro(tc, patm)
	# get factor to convert energy to water flux J/m^3
	econ <- s/(lv*pw*(s + gam))
	
	
	# ET
	
	# get actual evapotranspiration and potential (mm/day) *1000 from m3 to litres, t_conv_f to seconds to half hour
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
	q_l<-LE
	# get daily snowmelt
	if(!is.null(fluxnet_data$SW_OUT)){
		alb<-fluxnet_data$SW_OUT/SW_IN_F
		alb[fluxnet_data$SW_OUT<0]<-NA
		# alb[SW_IN_F<=0]<-NA
		# alb[fluxnet_data$SW_OUT>SW_IN_F]<-NA
		alb[alb>1]<-1.0
		alb_day<-apply.daily(alb,median,na.rm=TRUE)
		alb_day[is.infinite(alb_day)]<-NA
		snowmelt<-q_l
		snowmelt[snowmelt<0]<-NA
		snowmelt[tc<0]<-NA	
		LE[LE<0]<-NA
		snowmelt<-apply.daily(snowmelt*t_conv_f,sum,na.rm=TRUE)
		# snowmelt[snowmelt==0]<-NA
		snowmelt[alb_day<0.5 | is.na(alb_day)]<-0
		LE<-apply.daily(LE*t_conv_f,sum,na.rm=TRUE)
		LE[!is.na(snowmelt)]<-LE[!is.na(snowmelt)]-snowmelt[!is.na(snowmelt)]
		LE[LE<0]<-0
		aet<-LE*(1/lv)*(1/pw)*1000.0
		snowmelt<-snowmelt*1000/(pw*kfus)
		snowmelt[is.na(alb_day)]<-NA
		
	}else{
		
		LE[LE<0]<-NA
		LE<-apply.daily(LE*t_conv_f,sum,na.rm=TRUE)
		LE[LE<=0]<-NA
		aet<-LE*(1/lv)*(1/pw)*1000.0
		alb_day<-xts(rep(NA,length(aet)),time(aet))
		snowmelt<-xts(rep(NA,length(aet)),time(aet))
	}
	
	
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
	# aet<-xts(aet,ind)
	# aet<-apply.daily(aet,sum)
	
	
	# et_ind<-aet/pet
	# get soil moisture (mm) from volumetric measurements depth 0.5 m each
	if(!is.null(fluxnet_data$SWC_F_MDS_1)){
		sm_1<-fluxnet_data$SWC_F_MDS_1
	}else{
		if(!is.null(fluxnet_data$SWC_1)){
			sm_1<-fluxnet_data$SWC_1
		}else{
			sm_1<-fluxnet_data$P*NA
		}
	}
	if(!is.null(fluxnet_data$SWC_F_MDS_2)){
		sm_2<-fluxnet_data$SWC_F_MDS_2
	}else{
		if(!is.null(fluxnet_data$SWC_2)){
			sm_2<-fluxnet_data$SWC_2
		}else{
			sm_2<-fluxnet_data$P*NA
		}
	}
	if(!is.null(fluxnet_data$SWC_F_MDS_3)){
		sm_3<-fluxnet_data$SWC_F_MDS_3
	}else{
		if(!is.null(fluxnet_data$SWC_3)){
			sm_3<-fluxnet_data$SWC_3
		}else{
			sm_3<-fluxnet_data$P*NA
		}
	}
	sm<-cbind(sm_1,sm_2,sm_3)
	# assuming 0.5 meter depth
	sm<-rowMeans(sm,na.rm=TRUE)
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
	result<-merge.xts(aet,pet,eeq,TA_F,SW_IN_mean,P,sm,NETRAD/1e6,vpd,co2,ppfd,gpp,snowmelt,alb_day)
	names(result)<-c("aet","pet","eeq","tc","SW_in","pn","sm","netr","VPD","CO2","PPFD","GPP",'snowmelt','alb')
	
	return(result)
	
}
