#' unsSWC
#'
#' Estimates water content using an analytca integral of the profile described by Brooks and Corey 1964 and calculates depth to water table if there is any
#' @param  soil_data, vector retrieved from getSoilNRCS() or getSoilSite()
#' @param  uns_depth, depth to where to calculate water content in (m)
#' @param  wn, xts object, total water content of the profile from surface to bedrock (m)
#' @param  ouputdir, directory where to save the files
#' @import xts
#' @keywords splashTools
#' @export
#' @examples
#' unSWC()

unSWC<-function(soil_data,uns_depth,wn){
	# test
	# soil_data=soil_snow[[1]];uns_depth=depth_sm_snow[[1]];wn=sims$`SNTL:1243`[,1]
	# end testing
	UnsWater<-function(psi_m,z_uns,theta_r,theta_s,bub_press,lambda){
		#-----------------------------------------------------------------------
		# Input:    - float, Matric potential (Psi_m), mm
		#           - float, Depth to where to calculate water content (z_uns), m
		#           - float, residual moisture (theta_r), m3 m-3
		#           - float, saturated moisture (theha_s), m3 m-3
		#           - float, bubbling/air entry pressure (bub_press), m3 m-3
		#           - float, slope of the log curve thetha vs Psi_m (lambda), m3 m-3
		# Output:   float, integrated water from surface to depth z_uns, mm
		# Features: Estimates water content using an analytca integral of the profile described by Brooks and Corey 1964
		#          
		# Ref:      Brooks, R.H., Corey, A.T., 1964. Hydraulic properties of porous media. 
		#   		  Hydrology Papers No 17. Colorado State University. doi:10.13031/2013.40684
		#           
		#-----------------------------------------------------------------------
		z_uns<-z_uns*1000
	
		w_uns_z<- theta_r*z_uns+(((psi_m+z_uns)*(theta_r-theta_s)*(bub_press/(psi_m+z_uns))^lambda)/(lambda-1))
		w_uns_0<- theta_r*0+(((psi_m+0)*(theta_r-theta_s)*(bub_press/(psi_m+0))^lambda)/(lambda-1))
		
		w_uns<-w_uns_z-w_uns_0
		return(w_uns)
		
		
	}
	
	# wn<-simdf[,1]
	
	soil_info<-soil_hydro(sand=soil_data[1],clay=soil_data[2],OM=soil_data[3],fgravel =soil_data[4] ,bd = soil_data[5])
	SAT<-soil_info$SAT*(1-soil_data[4]/100)*soil_data[6]*1000
	RES<-soil_info$RES*(1-soil_data[4]/100)*soil_data[6]*1000
	theta_s<-as.numeric(SAT/(soil_data[6]*1000))
	theta_r<-as.numeric(RES/(soil_data[6]*1000))
	z_uns <- uns_depth
	lambda<-as.numeric(1/soil_info$B)
	bub_press<-as.numeric(soil_info$bubbling_p)
	theta_i<-wn/(soil_data[6]*1000)
	theta_i[theta_i>=theta_s]<-theta_s-0.05
	theta_i[theta_i<=theta_r]<-theta_r+0.01
	psi_m = bub_press/((((theta_i-theta_r)/(theta_s-theta_r)))^(1/lambda));
	psi_m<-as.numeric(psi_m)
	uns_wn<-UnsWater(psi_m,z_uns,theta_r,theta_s,bub_press,lambda)
	uns_theta<-uns_wn/(uns_depth*1000)
	uns_theta<-xts(uns_theta,time(wn))
	uns_wn<-xts(uns_wn,time(wn))
	wtd<-(bub_press-psi_m)/1000
	wtd[wtd>soil_data[6]]<-soil_data[6]
	merge.xts(uns_wn,wtd)
	
}
