plotPostMT_HDImeans2 = function(paramSampleVec, HDIlow,HDIhi, ylab=NULL ,
                      xlab=NULL , xlim=NULL, main=NULL, credMass=NULL, pltitle=NULL, 
                      showHDI = NULL, colflag = NULL, bfs = NULL, ylim=NULL, bfpos = NULL) {
  if ( is.null(xlab) ) xlab="Parameter"
  if ( is.null(main) ) main=""
  if ( is.null(ylab) ) ylab="Posterior Density"
  if ( is.null(credMass) ) credMass=.95 # 95% HDI as default
  if (is.null(showHDI)) showHDI <- 1
  if (is.null(colflag)) colflag <- 1
  
  if (is.null(xlim)) {
    xlim=range( c( 0 , paramSampleVec))
  }
  varnames=hdiLow=hdiHigh=NULL
  # in case xlim is set to 0, and some value has been given to xrange, center plot symmetrically on zero, using maximal extension in case range is set to 0, and given range otherwise
  if (xlim[1]==0)  {
    maxext <- max(abs(min(paramSampleVec$samples)), abs(max(paramSampleVec$samples)))  #largest extension into positive or negative range
    xlim = c(-maxext, maxext)   #centers plot symmetrically on zero
  }
  
  #HDI = HDIofMCMC( paramSampleVec$samples[paramSampleVec$classify == 2] , credMass )
  HDI <- tapply(paramSampleVec$samples, paramSampleVec$classify, FUN = HDIofMCMC)
  postSummary = matrix( NA , nrow=nlevels(paramSampleVec$classify) , ncol=11 , 
                        dimnames=list( c( 1:nlevels(paramSampleVec$classify)) , 
                                       c("mean","median","mode",
                                         "hdiMass","hdiLow","hdiHigh",
                                         "compVal","pcGTcompVal",
                                         "ROPElow","ROPEhigh","pcInROPE")))             
  postSummary[,"mean"] = aggregate(samples ~ classify, FUN = mean,data = paramSampleVec)[,2]
  postSummary[,"median"] = aggregate(samples ~ classify, FUN = median,data = paramSampleVec)[,2]
  mcmcDensity <- tapply(paramSampleVec$samples, paramSampleVec$classify, FUN = density)
  
  for (i in 1:nlevels(paramSampleVec$classify)) {
    postSummary[i,"mode"] = mcmcDensity[[i]]$x[which.max(mcmcDensity[[i]]$y)]
  }
  
  postSummary[,"hdiMass"]=credMass
  postSummary[,"hdiLow"]=sapply(HDI,function(x) x[1])
  postSummary[,"hdiHigh"]=sapply(HDI,function(x) x[2])
  
  postSummary <- as.data.frame(postSummary)
  
  
  densCurve <- tapply(paramSampleVec$samples, paramSampleVec$classify, FUN = density, adjust = 2)
  myy1 <- (seq(1:length(densCurve)))
  myy2 <- (seq(1:length(densCurve)))
  
  # Display the HDI.
  postSummary$myy1 <- myy1
  postSummary$myy2 <- myy2
  postSummary$varnames <- levels(paramSampleVec$varnames)
  postSummary$varnames <- as.character(postSummary$varnames)
  postSummary$varnames <- factor(postSummary$varnames, levels=unique(postSummary$varnames))
  
  colvals=c("red", "white", "black","grey52", "orange", "lightblue", "lightgreen")[1:nrow(postSummary)]

  # Plot HDIs with means
  if(colflag == 1){
    p1 <- ggplot(postSummary, aes(varnames, mean)) + geom_point(shape = 1, size = 5) + coord_flip() +
      geom_segment(aes(x = seq(1,nrow(postSummary)), xend = seq(1,nrow(postSummary)),
                       y = hdiLow, yend = hdiHigh), colour = "black", size=1.75, lineend = "round")  +
      geom_segment(aes(x = 0.5, xend = nrow(postSummary)+0.5, y = 0, yend = 0), colour = "red", size = 1, lineend = "round")
      #geom_point(aes(x = seq(1,nrow(postSummary)), y=mean, shape = 1), colour = "black", size=5, shape = rep(1,nrow(postSummary)))
  }  else {
    p1 <- ggplot(postSummary, aes(varnames, mean)) + geom_point(shape = 1, size = 5) + coord_flip() +
      geom_segment(aes(x = seq(1,nrow(postSummary)), xend = seq(1,nrow(postSummary)), y = hdiLow, yend = hdiHigh),
                   colour = "black", size=1.75, lineend = "round")# +
      #geom_point(aes(x = seq(1,nrow(postSummary)), y=mean), colour = "black", size=5, shape = rep(1,nrow(postSummary)))
  }
  
  
  myplot<-p1+ggtitle(main) + labs(title = pltitle) + theme_bw() +
    theme(legend.position = c(.80, .85)) + theme(legend.title=element_blank()) +
    theme(legend.text = element_text(size = 12)) +
    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12)) +
    theme(axis.title.y = element_text(size=14),axis.text.y  = element_text(size=12)) +
    theme(plot.title = element_text(size = rel(1.75), hjust = 0)) + theme(legend.position = "none") +
    labs(x=xlab, y=ylab) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          strip.text = element_text(size = 12), axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  # show bfs as numbers
  if (!is.null(bfs)){
    postSummary$bfs <- bfs
    myplot <- myplot + annotate("text", x=seq(1,nrow(postSummary)), y=bfpos, label= bfs, size = 5) + scale_y_continuous(limits = ylim)
  }
  return( myplot )
}