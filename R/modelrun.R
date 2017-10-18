#' Computes Bayes Factors for hierarchical linear models including continuous predictors using the Savage-Dickey density ratio
#' @title Bayes Factors, Posterior Samples, & DIC
#' @name modelrun
#' @param data a \code{data.frame} object with the data to be fitted in the long format.
#' @param dv \code{string} indicating the dependent variable of the model. Has to be normally distributed.
#' @param dat.str a \code{data.frame} object indicating the hierarchical structure in the model with column names "iv" and "type" and an arbitrary number of random variables as the following column names. iv indicates the name of an independent variable as in \code{data}, type its scale of measurement ("cont" for continuous or "cat" for categorical), and the following entries indicate whether a random effect should be modeled for this variable (1) or not (0). Continuous variables have to be entered before categorical variables. The name for the random variable(s) has to be the same as in \code{data}. A categorical variable with n levels is entered as n - 1 simple codes into the model with the first level of the variable being the reference category.
#' @param randvar.ia a \code{list} containing n \code{matrix} objects with n being the number of random variables. In each \code{matrix} the lower triangle can be used to declare the respective two-way interaction as random within the specific random variable. The row- and column- ordering of independent variables is the same as in \code{dat.str}. When not specified, interactions are only modeled as fixed effects by default.
#' @param corstr a \code{list} containing n \code{matrix} objects with n being the number of random variables. In each \code{matrix} the lower triangle can be used to assign correlations between predictors (including the intercept) for each random effect. The first row and column in each \code{matrix} object represents the intercept. The following rows and columns represent the independent variables with the same ordering as in \code{dat.str}. When not specified, no correlations are modeled by default.
#' @param nadapt number of MCMC steps to adapt the sampler (2000 by default).
#' @param nburn number of MCMC steps to burn in the sampler (2000 by default).
#' @param nsteps number of saved MCMC steps in all chains (100'000 by default).
#' @param checkconv indicates that convergence statistics of the main model parameters should be returned in the console and that figures of the chains should be plotted when set to 1 (0 by default).
#' @param mcmc.save.indiv indicates that the chains should be saved in a \code{data.frame} object when set to 1 (0 by default).
#' @param plot.post indicates that the 95 percent highest-density interval of the posterior of the group parameters should be plotted as a figure with the corresponding Bayes Factors when set to 1 (0 by default).
#' @param dic indicates that the deviation information criterion (Spiegelhalter, Best, Carlin, & Linde, 2002) should be computed for a given model when set to 1 (0 by default).
#' @param path defines directory where model is saved as .txt file and model name. Is set to file.path(tempdir(), "model.txt") by default.
#' @details The argument \code{corstr} can be used to model correlations between (a) pairs of predictors and (b) more than two predictors. When both is done within the same random variable, a predictor can only appear in (a) or (b).
#'
#'   \code{modelrun} z-standardizes the dependent variable and the continuous independent variables. To obtain the posteriors in the original scale they have to be retransformed.
#' \subsection{Savage Dickey}{
#' Bayes Factors are computed with the Savage-Dickey density ratio. We use the normal approximation (e.g., Wetzels, Raaijmakers, Jakab, & Wagenmakers, 2009) to estimate the density of the posterior.
#' }
#' @return returns a list with components:
#' \itemize{
#'  \item \code{bf}: a \code{data.frame} object with the Bayes Factor estimates of the group parameters (aka fixed effects).
#'  \item \code{mcmcdf}: a \code{data.frame} object with the saved MCMC chains.
#'  \item \code{dic}: DIC of the fitted model.
#' }
#' @author Thalmann, M., Niklaus, M. Part of this package uses code from John Kruschke.
#' @references
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & van der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64(4), 583.
#'
#' Wetzels, R., Raaijmakers, J. G. W., Jakab, E., & Wagenmakers, E.-J. (2009). How to quantify support for and against the null hypothesis: A flexible WinBUGS implementation of a default Bayesian t test. Psychonomic Bulletin & Review, 16(4), 752-760. https://doi.org/10.3758/PBR.16.4.752
#' @example examples/example.modelrun.R
#'
#' @importFrom stats dnorm sd update aggregate var density median
#' @importFrom reshape melt
#' @importFrom metRology dt.scaled
#' @importFrom rjags coda.samples load.module jags.model
#' @importFrom coda gelman.diag
#' @importFrom methods show
#' @importFrom grid grid.draw
#' @importFrom ggplot2 ggplot aes geom_point geom_segment coord_flip ggtitle labs theme_bw theme annotate scale_y_continuous element_blank element_text element_line rel aes_string
#' @importFrom graphics plot
#' @importFrom grDevices dev.new
#' @importFrom utils globalVariables
#' @export


modelrun <- function(data, dv, dat.str, randvar.ia = NULL, corstr = NULL, nadapt = NULL, nburn = NULL, nsteps = NULL,
                     checkconv = NULL, mcmc.save.indiv = NULL, plot.post = NULL ,dic = NULL,path=NULL){
  if (is.null(nadapt) ) nadapt=2000
  if (is.null(nburn) ) nburn=2000
  if (is.null(nsteps) ) nsteps=100000
  if (is.null(checkconv) ) checkconv=1
  if (is.null(mcmc.save.indiv)) mcmc.save.indiv = 0
  if (is.null(plot.post)) plot.post=0
  if (is.null(dic)) dic=0
  if (is.null(path)) path=file.path(tempdir(), "model.txt")
  ls<-0
  ncont<-sum(dat.str$type=="cont")
  # coding for cat vars when levels(cat)>2
  check <- as.character(dat.str$iv[dat.str$type == "cat"])
  add.vars <- names(dat.str)[3:ncol(dat.str)]
  add.vals <- matrix(dat.str[dat.str$iv %in% check,3:ncol(dat.str)],ncol=1)
  if (is.null(corstr)) corstr=rep(list(matrix(0)),length(add.vars))
  if (any(dat.str$type == "cat")){
    nrcat.act <- sum(dat.str$type == "cat")

    tmp <- lapply(X = data, FUN = levels)
    lvls <- tmp[check]
    codevar <- list()
    # change dat.str df when more than 2 levels for factor
    X.cat <- data.frame(matrix(NA, nrow = nrow(data), ncol = 1))
    for (i in 1:length(lvls)){
      ncodes.i <- length(lvls[[i]])-1
      add <- ncodes.i/(ncodes.i+1)
      codevar[[i]] <- matrix(NA, nrow = nrow(data), ncol = ncodes.i)
      name.i <- vector(length = ncodes.i)
      if(ncodes.i==1){
        for (c in 1:ncodes.i){
          name.i[c] <- paste0(check[i])
          codevar[[i]][,c][data[,check[i]] == levels(data[,check[i]])[c+1]] <- add
          codevar[[i]][,c][is.na(codevar[[i]][,c])] <- -add/ncodes.i
          assign(name.i[c], codevar[[i]][,c])
        }
      }else{
        for (c in 1:ncodes.i){
          name.i[c] <- paste0(check[i],".spl",c)
          codevar[[i]][,c][data[,check[i]] == levels(data[,check[i]])[c+1]] <- add
          codevar[[i]][,c][is.na(codevar[[i]][,c])] <- -add/ncodes.i
          assign(name.i[c], codevar[[i]][,c])
        }
      }
      # update dat.str df
      if(ncodes.i>1){
        dat.str <- dat.str[dat.str$iv!=names(lvls[i]),]
        dat.str.add <- as.data.frame(cbind(name.i, rep("cat", ncodes.i), matrix(unlist(rep(unname(add.vals[i,]), ncodes.i)),nrow = ncodes.i,byrow = T)))
        names(dat.str.add) <- c("iv", "type", c(add.vars))
        dat.str <- rbind(dat.str, dat.str.add)
      }

      change<-names(dat.str[3:ncol(dat.str)])
      if (length(change)>1){
        dat.str[,3:ncol(dat.str)]<-lapply(dat.str[,change],as.numeric)
      } else {dat.str[,3:ncol(dat.str)]<-as.numeric(dat.str[,change])}

      # and update data df
      data <- data[,which(names(data)!=check[i])]
      data <- cbind(data, codevar[[i]])
      names(data)[(ncol(data)-ncodes.i+1):ncol(data)] <- paste0("x",name.i)
      X.cat<-cbind(X.cat,data[(ncol(data)-ncodes.i+1):ncol(data)])
    }
    pre.names<-names(X.cat)[2:ncol(X.cat)]
    X.cat<-data.frame(X.cat[,2:ncol(X.cat)])
    names(X.cat)<-pre.names
    dat.str[,3] <- as.numeric(dat.str[,3])
    # change corstr mat when > 2 levels for factor
    if(any(unlist(lapply(corstr,sum))!=0)){
      ls <- unlist(unname(lapply(lvls, FUN = function(x) length(x))))
      longer <- which(unname(unlist(lapply(lvls, FUN = function(x) length(x)>2))))
      corstr.tmp<-corstr
      counter<-0
      for(i.pre in longer){
        counter<-counter+1
        ind<-1+ncont+i.pre
        for(i in 1:(ls[i.pre]-2)){
          stay1 <- 1+ncont+counter # 1 for intercept
          corstr.tmp<-lapply(corstr.tmp,FUN=function(x)x[1:stay1,1:stay1])
          corstr.tmp<-lapply(corstr.tmp,FUN=function(x)rbind(x, c(x[stay1,1:(stay1-1)],x[stay1,(stay1-1)])))
          corstr.tmp<-lapply(corstr.tmp,FUN=function(x)cbind(x, c(x[1:stay1,stay1],x[stay1,stay1])))
          counter<-counter+1
        }
        if(i.pre<max(longer)){
          remain<-length((ind+1):nrow(corstr.tmp[[1]]))
          tmp<-lapply(corstr,FUN=function(x)matrix(c(x[(ind+1):nrow(x),1:ind],rep(x[(ind+1):nrow(x),ind],remain)),nrow=1))
          corstr.tmp<-mapply(rbind,corstr.tmp,tmp,SIMPLIFY = FALSE)
          corstr.tmp<-lapply(corstr.tmp,FUN=function(x)cbind(x,x[,ncol(x)]))
        }
      }
      corstr<-corstr.tmp
    }

    # same for interaction df
    if(any(unlist(lapply(randvar.ia,sum))!=0)){
      ls <- unlist(unname(lapply(lvls, FUN = function(x) length(x))))
      longer <- which(unname(unlist(lapply(lvls, FUN = function(x) length(x)>2))))
      randvar.ia.tmp<-randvar.ia
      counter<-0
      for(i.pre in longer){
        counter<-counter+1
        ind<-ncont+i.pre
        for(i in 1:(ls[i.pre]-2)){
          stay2 <- ncont+counter # no intercept
          randvar.ia.tmp<-lapply(randvar.ia.tmp,FUN=function(x)x[1:stay2,1:stay2])
          randvar.ia.tmp<-lapply(randvar.ia.tmp,FUN=function(x)rbind(x, c(x[stay2,])))
          randvar.ia.tmp<-lapply(randvar.ia.tmp,FUN=function(x)cbind(x, c(x[,stay2])))
          counter<-counter+1
        }
        if(i.pre<max(longer)){
          remain<-length((ind+1):nrow(randvar.ia.tmp[[1]]))
          tmp<-lapply(randvar.ia,FUN=function(x)matrix(c(x[(ind+1):nrow(x),1:ind],rep(x[(ind+1):nrow(x),counter],remain)),nrow=1))
          randvar.ia.tmp<-mapply(rbind,randvar.ia.tmp,tmp,SIMPLIFY = FALSE)
          randvar.ia.tmp<-lapply(randvar.ia.tmp,FUN=function(x)cbind(x,x[,ncol(x)]))
        }
      }
      randvar.ia<-randvar.ia.tmp
    }
  }
  if (is.null(randvar.ia)) randvar.ia=replicate(length(names(dat.str[3:ncol(dat.str)])),
                                                matrix(0, nrow = nrow(dat.str), ncol = nrow(dat.str)), simplify = F)


  data <- data[!is.na(data[,dv]),]
  randnames <- names(dat.str[3:ncol(dat.str)])
  contnames <- as.character(dat.str$iv[dat.str$type == "cont"])
  catnames <- as.character(dat.str$iv[dat.str$type == "cat"])
  nrcont <- length(contnames)
  nrcat <- length(catnames)

  f.sub<-function(x){
    z=0
    if(x!=0){
      for(i in 1:x){z=z+i}
    }
    return(z)
  }
  nrIA.min<-0
  if(sum(ls)!=0){
    ls.calc <- matrix(ls,nrow=1)
    nrIA.min<-sum(unlist(apply((ls.calc-2),MARGIN = 2,FUN=f.sub)))
  }

  # cont as numeric predictors
  contvars = which(names(data) %in% c(as.character(dat.str$iv[dat.str$type == "cont"]), dv))
  data[,contvars] <- mapply(data[,contvars], FUN = as.character)
  data[,contvars] <- mapply(data[,contvars], FUN = as.numeric)

  # id and cat as factors
  catvars = which(names(data) %in% c(names(dat.str[3:ncol(dat.str)]), as.character(dat.str$iv[dat.str$type == "cat"])))
  if (length(catvars) > 1){
    data[,catvars] <- lapply(data[,catvars], as.factor)
  } else {
    data[,catvars] <- as.factor(data[,catvars])
  }

  for (i in randnames){
    levels(data[,i]) <- 1:length(unique(data[,i]))
  }

  # z-standardize ivs and cont dv
  mDV <- mean(data[,dv])
  sdDV <- sd(data[,dv])
  data$z.DV <- (data[,dv]-mDV)/sdDV
  nxs <- nrow(dat.str)
  X <- data.frame(matrix(NA, nrow = nrow(data), ncol = 1))

  for (i in contnames){
    name.m <- paste("m", i, sep = "")
    name.sd <- paste("sd", i, sep = "")
    name.z <- paste("z", i, sep = "")
    assign(name.m, mean(data[,i]))
    assign(name.sd, sd(data[,i]))
    assign(name.z, (data[,i]-mean(data[,i]))/sd(data[,i]))
    z.val <- as.name(paste("z", i, sep = ""))
    x.name <- paste("x", i, sep = "")
    do.call("<-", list(x.name, z.val))
    X <- do.call("cbind", list(X,  z.val))
    names(X)[which(i == contnames)+1] <- x.name
    data <- do.call("cbind", list(data,z.val))
  }
  if(exists("X.cat")){
    X <- cbind(X,X.cat)
  }
  X <- as.data.frame(X[,2:ncol(X), drop = FALSE])
  Ndata <- length(data$z.DV)

  dataList = list(
    y = data$z.DV ,
    Ndata = nrow(data)
  )

  name.cyc <- vector()
  size.cyc <- vector()
  for (i in randnames){
    name.cyc[i] <- paste("N", i, "s", sep = '')
    size.cyc[i] <- length(levels(data[,i]))
  }
  for (i in 1:length(randnames)){
    dataList[[name.cyc[i]]] <- size.cyc[i]
    dataList[[randnames[i]]] <- as.numeric(data[,randnames[i]])
  }

  for (i in names(X)){
    dataList[[i]] <- X[,which(i == names(X))]
  }

  # write jags model with required parameters
  parameters <- c(modeltext(dat.str, randvar.ia, corstr,path))
  parms.int <- unique(unlist(parameters[["MU"]]))
  parms.int <-parms.int[!is.na(parms.int)]
  parms.int <- parms.int[!grepl("mu0g", x = parms.int)]
  parms.int <- parms.int[!grepl("mu.corr", x = parms.int)]
  params.mon <- unique(c(unlist(parameters[["MU"]]), unlist(parameters[["SIGMA"]])))
  params.mon<-params.mon[!is.na(params.mon)]
  params.mon <- params.mon[!grepl("mu.corr", x = params.mon)]

  if (any(unlist(corstr)==1)){
    params.corr <- unique(unlist(parameters[["RHO"]]))
    params.mu.corr <- unique(unlist(parameters[["mu.corr"]]))
    params.mu.corr <- params.mu.corr[!grepl("mu.corr\\[1\\]", x = params.mu.corr)]

    # prepare dfs for inverse wishart
    wishdf <- parameters[["wishdf"]]
    for (i in 1:length(wishdf)){
      name.wish <- paste0("I", i)
      assign(name.wish, wishdf[i])
      dataList[[name.wish]] <- diag(c(1:wishdf[i]))
    }
  } else {wishdf <- logical(0)}
  #------------------------------------------------------------------------------
  ### RUN THE CHAINS ####
  load.module("dic")
  adaptSteps = nadapt #1000              # Number of steps to "tune" the samplers.
  burnInSteps = nburn # 1000            # Number of steps to "burn-in" the samplers.
  nChains = 3 #3                   # Number of chains to run.
  numSavedSteps= nsteps #15000           # Total number of steps in chains to save.
  thinSteps=1                # Number of steps to "thin" (1=keep every step).
  nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
  add <- NA
  if (mcmc.save.indiv){
    add <- parameters$b.save
  }
  mcmc.save <- c(c(params.mon, "tau"), add)
  if (any(unlist(corstr)==1)){
    mcmc.save <- c(mcmc.save, params.corr, params.mu.corr)
  }
  mcmc.save <- mcmc.save[!is.na(mcmc.save)]
  # Create, initialize, and adapt the model:
  jagsModel = jags.model( path , data=dataList ,
                          n.chains=nChains , n.adapt=adaptSteps )

  # Burn-in:
  cat( "Burning in the MCMC chain...\n" )
  update( jagsModel, n.iter=burnInSteps )
  # The saved MCMC chain:
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples(jagsModel , variable.names=c(mcmc.save,"deviance") ,
                             n.iter=nPerChain , thin=thinSteps)
  mcmcChain <- as.matrix(codaSamples)
  mcmcdf <- as.data.frame(mcmcChain)

  #### convergence diagnostics ####
  if (checkconv == 1){
    # diagnostics look good
    nparms <- length(params.mon)
    cycle <- ceiling(x = nparms/4)
    # cycle through hypermeans and hypersds
    for (i in 1:cycle){
      ind <- c(params.mon[((i-1)*4+1):(i*4)])
      if (any(is.na(ind))){
        ind <- ind[1:(min(which(is.na(ind)))-1)]
      }
      dev.new()
      plot( codaSamples[,ind])
      show( gelman.diag(codaSamples[,ind], multivariate = F))
    }
    # plot correlations for inspection as well
    if (length(wishdf)>0){
      for (dfs in 1:length(wishdf)){
        dummy <- matrix(0, nrow = wishdf[dfs], ncol = wishdf[dfs])
        ind <- which(lower.tri(dummy), arr.ind = TRUE)
        options(warn = -1)
        if (nrow(ind)>4){
          broke <- split(as.data.frame(ind), rep(1:ceiling(nrow(ind)/4),each=4))
          for (i in 1:length(broke)){
            ind <- broke[[i]]
            dev.new()
            plot(codaSamples[,paste0(params.corr[dfs], "[", ind[,1],",",ind[,2],"]")])
            show( gelman.diag(codaSamples[,paste0(params.corr[dfs], "[", ind[,1],",",ind[,2],"]")], multivariate = F))
          }
        } else{
          dev.new()
          plot(codaSamples[,paste0(params.corr[dfs], "[", ind[,1],",",ind[,2],"]")])
          show( gelman.diag(codaSamples[,paste0(params.corr[dfs], "[", ind[,1],",",ind[,2],"]")], multivariate = F))
        }
        options(warn = 0)
      }
    }
    # plot mus out of mu.cor df
    if(exists("params.mu.corr")){
      if(length(params.mu.corr)!=0){
        pl.ind <- match(params.mu.corr, names(mcmcdf))
        cycle <- ceiling(x = length(pl.ind)/4)
        for (i in 1:cycle){
          ind <- c(pl.ind[((i-1)*4+1):(i*4)])
          if (any(is.na(ind))){
            ind <- ind[1:(min(which(is.na(ind)))-1)]
          }
          dev.new()
          plot(codaSamples[,ind])
        }
      }

      # and plot mus not in mu.corr df
      if(any(parameters[["pl.ind"]])){
        ind <- which(parameters[["pl.ind"]]==1)
        show <- c(contnames, catnames)[ind]
        pl.ind <- grep(paste0("\\bmu",show, "\\b",collapse="|"),
                       names(mcmcdf))
        cycle <- ceiling(x = length(pl.ind)/4)
        for (i in 1:cycle){
          ind <- c(pl.ind[((i-1)*4+1):(i*4)])
          if (any(is.na(ind))){
            ind <- ind[1:(min(which(is.na(ind)))-1)]
          }
          plot(codaSamples[,ind])
        }
      }}
    if(length(parameters[["pl.nhcl"]])>0){
      show <- parameters[["pl.nhcl"]]
      cycle <- ceiling(x = length(show)/4)
      for (i in 1:cycle){
        dummy[(cycle*4-3):(cycle*4)] <- NA
        pl <- show[(cycle*4-3):(cycle*4)]
        pl <- pl[1:(min(which(is.na(pl)))-1)]
        plot(codaSamples[,pl])
      }
    }
  }


  #### computation of bfs, and plotting of hdis and bfs ####
  b_post <- as.data.frame(cbind(mcmcdf[,parms.int]))#, mcmcdf$muNsyllG, mcmcdf$muiaG))
  names(b_post) <- parms.int
  # b_post <- b_post[,2:ncol(b_post), drop = FALSE]

  scalecont <- sqrt(2)/4
  scalecat <- 1/2
  bf <- NA
  bf.names<-NA

  #### not from correlation structure ####
  if (any(!grepl("mu.corr", parms.int))){
    # dnorm approach to compute bfs
    counter <- 1
    pl.cont <- contnames[parameters[["pl.ind"]][1:nrcont]]
    if (any(!is.na(pl.cont))) {
      for (i in pl.cont){
        bf <- c(bf,dt.scaled(0,1,0,scalecont)/(dnorm(0, mean(b_post[,paste0("mu",i)]), sd(b_post[,paste0("mu",i)]))+1e-300))
        bf.names<-c(bf.names,i)
        # bf <- c(bf,dt.scaled(0,1,0,scalecont)/(dnorm(0, mean(b_post[,paste0("mu",pl.cont)]), sd(b_post[,paste0("mu",pl.cont)]))+1e-300))
        counter <- counter + 1
      }
    }
    if (length(parameters[["pl.nhclcont"]])>0) {
      for (i in parameters[["pl.nhclcont"]]){
        bf <- c(bf,dt.scaled(0,1,0,scalecont)/(dnorm(0, mean(b_post[,i]), sd(b_post[,i]))+1e-300))
        bf.names<-c(bf.names,i)
        counter <- counter + 1
      }
    }

    pl.cat <- catnames[parameters[["pl.ind"]][(1+nrcont):(nrcat+nrcont)]]
    if (any(!is.na(pl.cat))) {
      for (i in pl.cat){
        bf <- c(bf,dt.scaled(0,1,0,scalecat)/(dnorm(0, mean(b_post[,paste0("mu",i)]), sd(b_post[,paste0("mu",i)]))+1e-300))
        bf.names<-c(bf.names,i)
        # bf <- c(bf,dt.scaled(0,1,0,scalecat)/(dnorm(0, mean(b_post[,paste0("mu",pl.cat)]), sd(b_post[,paste0("mu",pl.cat)]))+1e-300))
        counter <- counter + 1
      }
    }
    if (length((parameters[["pl.nhclcat"]]))) {
      for (i in parameters[["pl.nhclcat"]]){
        bf <- c(bf,dt.scaled(0,1,0,scalecat)/(dnorm(0, mean(b_post[,i]), sd(b_post[,i]))+1e-300))
        bf.names<-c(bf.names,i)
        counter <- counter + 1
      }
    }

    if (any(dat.str$type == "cat")){
      nrIA <- (length(catnames)+length(contnames)-1)*(nrcat.act+length(contnames))/2
      nrIA<-(length(catnames)+length(contnames))*(length(catnames)+length(contnames)-1)/2
      nrIA <- nrIA-nrIA.min
    } else{nrIA <- (length(catnames)+length(contnames)-1)*(length(catnames)+length(contnames))/2}
    if (nrIA > 0) {
      for (k in 1:(nrIA)){
        bf <- c(bf,dt.scaled(0,1,0,scalecat)/(dnorm(0, mean(b_post[,counter]), sd(b_post[,counter])))+1e-300)
        bf.names<-c(bf.names,names(b_post)[counter])
        counter <- counter + 1
      }
    }
    bf <- bf[2:length(bf)]
    bf1 <- prettyNum(bf, digits = 2)
    bf.names<-bf.names[2:length(bf.names)]
    bf.tog1<-data.frame(bf.names,bf=bf1)
    # plotting of 95% hdis with respective bfs
    plot.bs <- melt(b_post,id.vars = NULL)
    names(plot.bs) <- c("classify", "samples")
    plot.bs$classify <- factor(plot.bs$classify)
    plot.bs$varnames <- plot.bs$classify
    plot.bs$varnames <- factor(plot.bs$varnames)
    pl.post <- plotPostMT_HDImeans2(plot.bs, xlab="" , ylab = "Parameter Estimate\n", main="",
                                    ylim = c(min(plot.bs$samples)-.1, max(plot.bs$samples)+.1),
                                    showHDI = 1, colflag = 1, bfs = bf1, bfpos = min(plot.bs$samples))
    if(plot.post==1){
      dev.new()
      grid.draw(pl.post)} # n.b. posteriors are in z-transformed space
    #### from correlation structure ####
  }
  if(exists("params.mu.corr")){
    if (any(grepl("mu.corr", params.mu.corr))){

      bf <- NA
      bf.names<-NA
      comp <- NA
      ind <- (1+1):(1+nrcont)
      ind.plot <- NA
      for (i in ind){
        comp <- append(comp, which(grepl(paste0("\\[",i,"\\]"), params.mu.corr)))
      }
      comp <- comp[2:length(comp)]
      if (!is.na(comp[1])) {
        for (i in ind){
          name <- paste0("mu.corr","\\[",i, "\\]")
          ind.comp <- which(grepl(name,names(mcmcdf)))
          ind.plot <- c(ind.plot, ind.comp)
          bf <- c(bf,dt.scaled(0,1,0,scalecont)/(dnorm(0, mean(mcmcdf[,ind.comp]), sd(mcmcdf[,ind.comp]))+1e-300))
          bf.names<-c(bf.names,parameters[["corrnames"]][i-1])
        }
      }

      comp <- NA
      ind <- (1+nrcont+1):(1+nrcat+nrcont)
      for (i in ind){
        comp <- append(comp, which(grepl(paste0("\\[",i,"\\]"), params.mu.corr)))
      }
      comp <- comp[2:length(comp)]
      if (!is.na(comp[1])) {
        for (i in ind){
          name <- paste0("mu.corr","\\[",i, "\\]")
          ind.comp <- which(grepl(name,names(mcmcdf)))
          ind.plot <- c(ind.plot, ind.comp)
          bf <- c(bf,dt.scaled(0,1,0,scalecat)/(dnorm(0, mean(mcmcdf[,ind.comp]), sd(mcmcdf[,ind.comp]))+1e-300))
          bf.names<-c(bf.names,parameters[["corrnames"]][i-1])
        }
      }


      bf <- bf[2:length(bf)]
      bf2 <- prettyNum(bf, digits = 2)
      bf.names<-bf.names[2:length(bf.names)]
      bf.tog2<-data.frame(bf.names,bf=bf2)
      ind.plot <- ind.plot[2:length(ind.plot)]
      # depending on simulation study, de(-comment) the following line of code
      # bf <- log(bf)
      # plotting of 95% hdis with respective bfs
      plot.bs <- melt(mcmcdf[,ind.plot],id.vars = NULL)
      names(plot.bs) <- c("classify", "samples")
      plot.bs$classify <- factor(plot.bs$classify)
      levels(plot.bs$classify)<-bf.names
      plot.bs$varnames <- plot.bs$classify
      plot.bs$varnames <- factor(plot.bs$varnames)
      pl.post <- plotPostMT_HDImeans2(plot.bs, xlab="" , ylab = "Parameter Estimate\n", main="",
                                      ylim = c(min(plot.bs$samples)-.1, max(plot.bs$samples)+.1),
                                      showHDI = 1, colflag = 1, bfs = bf2, bfpos = min(plot.bs$samples))
      if(plot.post==1){
        dev.new()
        grid.draw(pl.post)} # n.b. posteriors are in z-transformed space
      #
    }
  }
  DIC<-NULL
  if(dic==1){
    meanDev <- mean(mcmcChain[,"deviance"])
    pD <- 0.5*var(mcmcChain[,"deviance"])
    DIC <- meanDev + pD}

  if (exists("bf1") & exists("bf2")){
    return(list(rbind(bf.tog2,bf.tog1), mcmcdf, DIC))
  }
  if (exists("bf1") & !exists("bf2")){
    return(list(bf.tog1, mcmcdf, DIC))
  }
  if (!exists("bf1") & exists("bf2")){
    return(list(bf.tog2, mcmcdf, DIC))
  }
}
