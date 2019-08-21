#' aggRaster
#'
#' Aggregate raster bricks from daily to monthly in parallel
#' @param  x,ind.months=,ouputdir
#' @import raster  
#' @keywords splash
#' @export
#' @examples
#' splash.grid()	
	
	aggRaster<-function(x,func='mean',inmem=FALSE,varnam='wn',outdir=getwd()){
			###############################################################################################
			# 00. create array for results, fluxes: mm/day, storages (wn, snow): mm
			###############################################################################################
			y<-as.numeric(unique(format(getZ(x),'%Y')))
			ind.months<-seq(getZ(x)[1],getZ(x)[length(getZ(x))], by="month")					
			nm <-length(ind.months)
			# actual soil moisture
			out<-brick(nrows=nrow(x), ncols=ncol(x), crs=crs(x), nl=nm)
			extent(x)<-extent(x)
			out<-setZ(out,ind.months)
			indmonth<-format(getZ(x),'%Y-%m')
			setwd(outdir)
			###############################################################################################
			# 01. set the clusters for parallel computing
			###############################################################################################	
			cl <- getCluster()
			on.exit( returnCluster() )
			nodes <- length(cl)
			bs <- blockSize(x, minblocks=nodes*5)
			parallel::clusterExport(cl, varlist=c('x','func','indmonth','bs'),envir=environment()) 
			pb <- pbCreate(bs$n)
			pb <- txtProgressBar(min=1,max = bs$n, style = 1)
			###############################################################################################
			# 02. create the functions to send to the workers, split the data in chunks
			###############################################################################################	
			clFun <- function(i) {
				
				x_block<-getValues(x,bs$row[i], bs$nrows[i])
				# do calculations
				if(func=='mean'){
					result<-t(sweep(rowsum(t(x_block),indmonth,na.rm = FALSE), 1, table(indmonth), "/") )
				}else{
					result<-t(rowsum(t(x_block),indmonth,na.rm = FALSE))
				}
				
				return(result)
			}
			###############################################################################################
			# 03. send tasks to the nodes
			###############################################################################################
			for (i in 1:nodes) {
				parallel:::sendCall(cl[[i]], clFun, i, tag=i)
			}
			###############################################################################################
			# 04. write to the disk on the go, or save to the ram
			###############################################################################################
			if(varnam=='wn'){
				longname='monthly soil moisture'
			}else if(varnam=='ro'){
				longname='monthly runoff'
			}else if(varnam=='pet'){
				longname='monthly potential evapotranspiration'
			}else if(varnam=='aet'){
				longname='monthly actual evapotranspiration'
			}else if(varnam=='snow'){
				longname='monthly snow water equivalent'
			}else if(varnam=='cond'){
				longname='monthly condensation'
			}else if(varnam=='bflow'){
				longname='monthly baseflow'
			}else if(varnam=='sw_in'){
				longname='monthly solar radiation'
			}
			if(!inmem){
				out<-writeStart(out,filename=paste0(outdir,"/",y[1],"_",y[length(y)],".",varnam,".","nc"),format="CDF",overwrite=TRUE,varname="wn", varunit="mm",
					longname=longname, xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
				
				
			}else {
				matout <- matrix(ncol=nlayers(out), nrow=ncell(out))
				endind<-cumsum(bs$nrows*out@ncols)
				startind<-c(1,endind+1)    
			}
			###############################################################################################
			# 05. receive results from the nodes
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
					out <- writeValues(out,d$value$value, bs$row[b])
					
				} else {
					
					matout[startind[b]:endind[b],] <- d$value$value
				}
				
				# need to send more data?
				ni <- nodes + i
				if (ni <= bs$n) {
					parallel:::sendCall(cl[[d$node]], clFun, ni, tag=ni)
				}
				setTxtProgressBar(pb,i)
			}
			###############################################################################################
			# 06. close connection with the files, or assign valueas to the raster objects
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
