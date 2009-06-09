## given a fixed FDR, compute the power
# assume Poisson(R*C/S),  
#  R: region size; C: average depth-coverage; S: read Size

# type I error
typeI <- c(0.01, 0.001, 0.0001, 0.00001, 0.000001)

                                        
# scale of gamma distribution: need more investigations
gammascale <- 2

# read size
size <- c(25,35,50)

# average coverage
coverage <- c(0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9, 1.0,1.5, 2.0, 3.0, 4.0, 5.0)

# region size -- resolution
region <- c(100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000,7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000, 22000, 24000, 27000, 30000)

# result matrix 
m <- matrix(, nrow =0, ncol=7)
# set the column names
dimnames(m)[[2]] <- c("readSize", "avgDepthCov", "regionSize", "PoissonPower", "GammaPower", "typeIError", "ratio")

# loops
for (s in size) {
   for (e in typeI) {
	for (cov in coverage) {
	  	for (r in region) {  # loop in resolution
    	ratio <- cov / s
    	numavg <- r * cov / s
    
# Poisson model
# cutoff to get type I error 
    pcutoff <- qpois( 1 - e, numavg)
# typeII error 
    ptypeII <- ppois(pcutoff, 1.5*numavg)
    ppower <- round(1 - ptypeII, 4)

# gamma model
    nshape <- numavg / gammascale
    gshape <- 1.5 * nshape
    gcutoff <- qgamma(1 - e, nshape, scale=gammascale)
    gtypeII <- pgamma(gcutoff, shape = gshape, scale=gammascale)
    gpower <- round(1 - gtypeII, 4)
    
    m <- rbind(m, c(s, cov, r, ppower, gpower, e, ratio))

  }
}
}
}
z <- data.frame(m)

write.table(z, "temp.out", quote=FALSE)

