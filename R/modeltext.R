# function to write required jags model as a text file

modeltext = function(dat.str, randvar.ia, corstr,path){

  cont <- as.character(dat.str$iv[dat.str$type == "cont"])
  cat <- as.character(dat.str$iv[dat.str$type == "cat"])
  if (length(randvar.ia) == 1){
    conthcl.pre <- cont[as.logical(dat.str[dat.str$type == "cont",3:ncol(dat.str)])]
    cathcl.pre <- cat[as.logical(dat.str[dat.str$type == "cat",3:ncol(dat.str)])]
    if (length(conthcl.pre)==0){conthcl.pre<-rep(0,length(cont))}
    if (length(cathcl.pre)==0){cathcl.pre<-rep(0,length(cat))}
  } else {
    conthcl.pre <- cont[any(dat.str[dat.str$type == "cont",3:ncol(dat.str)]>0)] # as.logical was rowSums before
    cathcl.pre <- cat[any(dat.str[dat.str$type == "cat",3:ncol(dat.str)]>0)] # as.logical was rowSums before
  }

  randvar <- names(dat.str[3:ncol(dat.str)])

  allnames <- c(cont, cat)
  conthcl <- as.numeric(cont == conthcl.pre)
  cathcl <- as.numeric(cat == cathcl.pre)
  nrhcl <- sum(conthcl, cathcl, unlist(randvar.ia)) + 1 # random slopes plus intercept
  nrcont <- length(cont)
  nrcat <- length(cat)
  nrand <- length(randvar)
  nrIA <- (nrcat+nrcont-1)*(nrcat+nrcont)/2
  # initialize variables for required parameters on different levels of hierarchy
  bcont <- replicate(nrcont, matrix(NA,nrow = nrand, ncol = 3), simplify=F)
  bcat <- replicate(nrcat, matrix(NA,nrow = nrand, ncol = 3), simplify=F)
  # intercepts for random variables
  bI <- replicate(nrand, rep(NA,3), simplify=F)
  bIAs <- replicate(((nrcont+nrcat)^2-(nrcont+nrcat))/2, matrix(NA,nrow = nrand, ncol = 3), simplify=F)
  mucont <- vector("list", 2)
  mucat <- vector("list", 2)
  muIAs <- vector("list", 2)
  precont <- matrix()
  precat <- matrix()
  preIAs <- matrix()
  pl.pre<-numeric()

  taucont <- matrix(NA, nrow = nrcont, ncol = nrand)
  taucat <- matrix(NA, nrow = nrcat, ncol = nrand)
  tauIAs <- matrix(NA, nrow = nrIA, ncol = nrand)
  sigmacont <- matrix(NA, nrow = nrcont, ncol = nrand)
  sigmacat <- matrix(NA, nrow = nrcat, ncol = nrand)
  sigmaIAs <- matrix(NA, nrow = nrIA, ncol = nrand)

  # correlation between random variables: matrices of what order?
  multivar <- rep(list(NULL),nrand)
  bivar <- rep(list(NULL),nrand)
  multicol <- rep(list(1),nrand)
  for (i in 1:nrand){
    threshold <- 1
    test <- which(corstr[[i]]==1, arr.ind = TRUE)
    tab <- table(test)
    tmp <- as.numeric(names(tab)[tab == threshold])
    if (length(tmp)!=0){
      biv.low <- sort(tmp)[1:(length(tmp)/2)]
      biv.up <- sort(tmp)[floor((length(tmp)/2)+1):length(tmp)]
      bivar[[i]] <- matrix(c(biv.low, biv.up), ncol= 2,byrow = TRUE)
    }
    tmp <- as.numeric(names(tab)[tab > threshold])
    if (length(tmp)>=0){
      multicol[[i]] <- max(length(tmp),multicol[[i]])
      if (length(tmp)!=0){
        multivar[[i]] <- matrix(as.numeric(names(tab)[tab > threshold]), nrow = 1, byrow = TRUE)
      } else { multivar[i] <- list(NULL)}
    }
  }

  tau0g<- vector("list",nrand)
  sigma0<-vector("list",nrand)
  mu0g <- vector("list",nrand)
  for (i in 1:length(randvar)){
    if (all(multivar[[i]]!=1)&all(bivar[[i]]!=1)){
      tau0g[[i]] <- paste("tau0g", ".", randvar[i], sep = '')
      sigma0[[i]] <- paste("sigma0", ".", randvar[i], sep = '')
    }
  }

  bimultivar <- mapply(cbind,multivar,bivar,SIMPLIFY=FALSE)
  ind<-lapply(bimultivar,function(x)sum(x==1))
  for(i in 1:length(multivar)){
    if(any(as.logical(unlist(ind)))){
      if (ind[[i]]==0){
        mu0g[[i]] <- "mu.corr[1]"
      }
    } else{
      mu0g[[i]]<-"mu0g"
    }
  }

  # random intercept
  for (i in 1:nrand){
    corrcount <- 0
    if (all(multivar[[i]]!=1) & all(bivar[[i]]!=1)){
      bI[[i]][[1]] <- paste("b", randvar[i], "[", randvar[i],"[i]]", sep = '')
      bI[[i]][[2]] <- paste("b", randvar[i], "[n]", sep = '')
      bI[[i]][[3]] <- paste("b", randvar[i], sep = '')
    } else {
      if (any(multivar[[i]]==1)){
        tmp <- which(multivar[[i]]==1)[1]
        corrcount <- corrcount + 1
        bI[[i]][[1]] <- paste("b", randvar[i], paste(multivar[[i]][tmp,], collapse = ''),
                              "[", randvar[i],"[i],",corrcount, "]" , sep = '')
        bI[[i]][[2]] <- paste("b", randvar[i], paste(multivar[[i]][tmp,], collapse = ''),  "[n", ",", corrcount, "]", sep = '')
        bI[[i]][[3]] <- paste("b", randvar[i], paste(multivar[[i]][tmp,], collapse = ''), sep = '')
      }

      if (any(bivar[[i]] == 1)){
        tmp <- which(bivar[[i]]==1 ,arr.ind = TRUE)[1]
        corrcount <- corrcount + 1
        bI[[i]][[1]] <- paste("b", randvar[i], paste(bivar[[i]][tmp,], collapse = ''),
                              "[", randvar[i],"[i],",corrcount, "]" , sep = '')
        bI[[i]][[2]] <- paste("b", randvar[i], paste(bivar[[i]][tmp,], collapse = ''),  "[n", ",", corrcount, "]", sep = '')
        bI[[i]][[3]] <- paste("b", randvar[i], paste(bivar[[i]][tmp,], collapse = ''), sep = '')
      }
    }
  }

  # where are the intercepts?
  tmp.multi<-unlist(lapply(multivar,FUN=function(x)sum(x==1)))
  tmp.bi<-unlist(lapply(bivar,FUN=function(x)sum(x==1)))

  corrcount <- matrix(c(tmp.multi,tmp.bi),byrow = TRUE,nrow = 2)
  # random slopes and group parameters
  contcount <- matrix(0,nrow = nrcont, ncol = nrand)
  nhclcont <- rep(0,times = nrcont)
  if (nrcont > 0){
    for (ncont in 1:nrcont){
      onlyfix <- 0
      if(conthcl[ncont] == 1){
        for (i in which(dat.str[dat.str$type == "cont",][ncont,3:ncol(dat.str)]==1)){
          if (all(multivar[[i]]!=1+ncont) & all(bivar[[i]]!=1+ncont)){
            bcont[[ncont]][i,1] <- paste("b", cont[ncont], ".", randvar[i],"[", randvar[i],"[i]]",
                                         " * x", cont[ncont], "[i]", sep = "")
            bcont[[ncont]][i,2] <- paste("b", cont[ncont], ".", randvar[i],"[n]", sep = "")
            bcont[[ncont]][i,3] <- paste("b", cont[ncont], ".", randvar[i], sep = "")
            contcount[ncont,i] <- contcount[ncont,i] + 1
          } else {
            if (any(multivar[[i]]==(1+ncont))){
              tmp <- which(multivar[[i]]==1+ncont, arr.ind = TRUE)[1]
              corrcount[1,i] <- corrcount[1,i] + 1
              bcont[[ncont]][i,1] <- paste("b", randvar[i], paste(multivar[[i]][tmp,], collapse = ''),"[",
                                           randvar[i],"[i],",corrcount[1,i], "]", " * x", cont[ncont], "[i]", sep = "")
              bcont[[ncont]][i,2] <- paste("b", randvar[i], paste(multivar[[i]][tmp,], collapse = ''), "[n", ",", corrcount[1,i], "]", sep = "")
              bcont[[ncont]][i,3] <- paste("b", randvar[i], paste(multivar[[i]][tmp,], collapse = ''), sep = "")
            } else if (any(bivar[[i]] == (1+ncont))){
              tmp <- which(bivar[[i]]==(1+ncont) ,arr.ind = TRUE)[1]
              corrcount[2,i] <- corrcount[2,i] + 1
              bcont[[ncont]][i,1] <- paste("b", randvar[i], paste(bivar[[i]][tmp,], collapse = ''),"[",
                                           randvar[i],"[i],",corrcount[2,i], "]", " * x", cont[ncont], "[i]", sep = "")
              bcont[[ncont]][i,2] <- paste("b", randvar[i], paste(bivar[[i]][tmp,], collapse = ''), "[n", ",", corrcount[2,i], "]", sep = "")
              bcont[[ncont]][i,3] <- paste("b", randvar[i], paste(bivar[[i]][tmp,], collapse = ''), sep = "")
            }
          }
        }
      } else {mucont[[1]][[ncont]] <- paste("mu", cont[ncont]," * x", cont[ncont], "[i]", sep = "")
      onlyfix <- 1}
      if (all(unlist(multivar)!=(1+ncont)) & all(unlist(bivar)!=(1+ncont))){
        mucont[[2]][[ncont]] <- paste("mu", cont[ncont], sep = "")
        # plus pre for rescaling as cauchy
        precont[[ncont]] <- paste("pre", cont[ncont], sep = "")
        if (conthcl[ncont]!=1){nhclcont[ncont] <- 1}
        for (i in which(dat.str[dat.str$type == "cont",][ncont,3:ncol(dat.str)]==1)){
          taucont[ncont,i] <- paste("tau", cont[ncont], ".", randvar[i], sep = "")
          # plus sigma for rescaling prec as sd
          sigmacont[ncont,i] <- paste("sigma", cont[ncont], ".", randvar[i], sep = "")
        }}
      # random slopes for cat pred only when required
      if (onlyfix == 0 &  sum(contcount) != sum(dat.str[dat.str$type == "cont",][1:ncont,3:ncol(dat.str)]==1)){
        mucont[[2]][ncont] <- paste0("mu.corr[", (1+ncont),"]")
        for(i in which(contcount == 1, arr.ind = TRUE)[2])
          taucont[ncont,i] <- paste("tau", cont[ncont], ".", randvar[i], sep = "")
        # plus sigma for rescaling prec as sd
        sigmacont[ncont,i] <- paste("sigma", cont[ncont], ".", randvar[i], sep = "")
      }
    }
    pl.pre <- contcount
  }else{bcont <- rm(bcont); mucont <- list(NA, NA); precont <- rm(precont); taucont <- rm(taucont); sigmacont <- rm(sigmacont)}
  catcount <- matrix(0,nrow = nrcat, ncol = nrand)
  nhclcat <- rep(0,times = nrcat)
  if (nrcat > 0){
    for (ncat in 1:nrcat){
      onlyfix <- 0
      # random slopes for continuous predictors
      if(cathcl[ncat] == 1){
        # nrs <- which(dat.str[dat.str$type == "cat",][ncat,3:ncol(dat.str)]==1)
        # nrs <- unique(which(as.matrix(dat.str[dat.str$type == "cat",][,3:ncol(dat.str)]==1),arr.ind = TRUE)[,2])
        # for (i in nrs){
        for (i in which(dat.str[dat.str$type == "cat",][ncat,3:ncol(dat.str)]==1)){
          if (all(multivar[[i]]!=(1+ncont+ncat)) & all(bivar[[i]]!=(1+ncont+ncat))){
            bcat[[ncat]][i,1]  <- paste("b", cat[ncat], ".", randvar[i], "[", randvar[i],"[i]]", " * x", cat[ncat], "[i]", sep = "")
            bcat[[ncat]][i,2] <- paste("b", cat[ncat], ".", randvar[i], "[n]", sep = "")
            bcat[[ncat]][i,3] <- paste("b", cat[ncat], ".", randvar[i], sep = "")
            catcount[ncat, i] <- catcount[ncat, i] + 1
          } else {
            if (any(multivar[[i]]==(1+ncont+ncat))){
              tmp <- which(multivar[[i]]==(1+ncont+ncat), arr.ind = TRUE)[1]
              corrcount[1,i] <- corrcount[1,i] + 1
              bcat[[ncat]][i,1]  <- paste("b", randvar[i],paste(multivar[[i]][tmp,], collapse = ''),
                                          "[", randvar[i],"[i],",corrcount[1,i], "]", " * x", cat[ncat], "[i]", sep = "")
              bcat[[ncat]][i,2] <- paste("b", randvar[i],paste(multivar[[i]][tmp,], collapse = ''), "[n", ",", corrcount[1,i], "]", sep = "")
              bcat[[ncat]][i,3] <- paste("b", randvar[i],paste(multivar[[i]][tmp,], collapse = ''), sep = "")
            } else if (any(bivar[[i]]==(1+ncont+ncat))){
              tmp <- which(bivar[[i]]==(1+ncont+ncat), arr.ind = TRUE)[1]
              corrcount[2,i] <- corrcount[2,i] + 1
              bcat[[ncat]][i,1]  <- paste("b", randvar[i],paste(bivar[[i]][tmp,], collapse = ''),"[",
                                          randvar[i],"[i],",corrcount[2,i], "]", " * x", cat[ncat], "[i]", sep = "")
              bcat[[ncat]][i,2] <- paste("b", randvar[i],paste(bivar[[i]][tmp,], collapse = ''), "[n", ",",
                                         corrcount[2,i], "]", sep = "")
              bcat[[ncat]][i,3] <- paste("b", randvar[i],paste(bivar[[i]][tmp,], collapse = ''), sep = "")
            }
          }
        }
      } else {mucat[[1]][ncat] <- paste("mu", cat[ncat]," * x", cat[ncat], "[i]", sep = "")
      onlyfix <- 1}
      if (all(unlist(multivar)!=(1+nrcont+ncat)) & all(unlist(bivar)!=(1+nrcont+ncat))){
        mucat[[2]][[ncat]] <- paste("mu", cat[ncat], sep = "")
        # plus pre for rescaling as cauchy
        precat[[ncat]] <- paste("pre", cat[ncat], sep = "")
        if (cathcl[ncat]!=1){nhclcat[ncat] <- 1}
        for (i in which(dat.str[dat.str$type == "cat",][ncat,3:ncol(dat.str)]==1)){
          taucat[ncat,i] <- paste("tau", cat[ncat], ".", randvar[i], sep = "")
          # plus sigma for rescaling prec as sd
          sigmacat[ncat,i] <- paste("sigma", cat[ncat], ".", randvar[i], sep = "")
        }
      }
      if (onlyfix == 0 &   sum(catcount) != sum(dat.str[dat.str$type == "cat",][1:ncat,3:ncol(dat.str)]==1)){
        mucat[[2]][ncat] <- paste0("mu.corr[", (1+nrcont+ncat),"]")
        for(i in which(catcount == 1, arr.ind = TRUE)[2])
          taucat[ncat,i] <- paste("tau", cat[ncat], ".", randvar[i], sep = "")
        # plus sigma for rescaling prec as sd
        sigmacat[ncat,i] <- paste("sigma", cat[ncat], ".", randvar[i], sep = "")
      }
    }
  }else{bcat <- rm(bcat); mucat <- list(NA, NA); precat <- rm(precat); taucat <- rm(taucat); sigmacat <- rm(sigmacat) }
  pl.pre <- rbind(pl.pre, catcount)

  pl.ind <- rowSums(pl.pre) == nrand
  # covariance structure if required
  bi.yes <- which(lapply(bivar,length)>0)
  multi.yes <- which(lapply(multivar, length)>0)

  multipart <- replicate(nrand, matrix(NA,nrow = (max(lengths(multivar))), ncol = 1), simplify=F)
  bipart <- replicate(nrand, matrix(NA,nrow = (max(lengths(bivar))/2), ncol = 1), simplify=F)
  # multivar<-lapply(multivar,function(x)as.matrix(x))
  # bivar<-lapply(bivar,function(x)as.matrix(x))
  for (i in multi.yes){
    for (j in 1:nrow(multivar[[i]])){
      multipart[[i]][j,] <- paste("b",randvar[i],paste(unique(multivar[[i]][j,]), collapse = ""),
                                  "[n, 1:",length(multivar[[i]]),
                                  "] ~ dmnorm (mu.corr[c(",paste(unique(multivar[[i]][j,]), collapse = ","),")]",
                                  ", SigmaInv",randvar[i],paste(unique(multivar[[i]][j,]), collapse = ""),
                                  "[1:",length(multivar[[i]]),",1:",length(multivar[[i]]),"])","\n", sep = "")
    }
  }
  for (i in bi.yes){
    for (j in 1:nrow(bivar[[i]])){
      bipart[[i]][j,] <- paste("b",randvar[i],paste(unique(bivar[[i]][j,]), collapse = ""),
                               "[n, 1:",length(bivar[[i]][j,]),
                               "] ~ dmnorm (mu.corr[c(",paste(unique(bivar[[i]][j,]), collapse = ","),")]",
                               ", SigmaInv", randvar[i],paste(as.vector(unique(bivar[[i]][j,])), collapse = ""),
                               "[1:", length(bivar[[i]][j,]),",1:", length(bivar[[i]][j,]),"])","\n", sep = "")
    }
  }
  bipart <- lapply(bipart, function(x) x[!is.na(x)])
  multipart <- lapply(multipart, function(x) x[!is.na(x)])

  mupart.corr <- NA
  pre.corr <- NA
  sigmainv.corr <- NA
  sigma.corr <- NA
  wishdf <- NA
  corr.names<-NA

  for (both in union(unique(unlist(multivar)),unique(unlist(bivar)))){
    if (both == 1){
      mupart.corr <- append(mupart.corr, paste("mu.corr[", both, "] ~ dnorm(0,1)","\n", sep = ""))
    } else {
      mupart.corr <- append(mupart.corr, paste("mu.corr[", both, "] <- pre", both, "* scale",
                                               c("cont","cat")[1+as.numeric(both>(nrcont+1))],"\n", sep = ""))
      pre.corr <- append(pre.corr, paste("pre", both, " ~ dt(0,1,1)","\n", sep = ""))
      corr.names<-c(corr.names,as.character(dat.str$iv[both-1]))
    }
  }
  corr.names<-corr.names[!is.na(corr.names)]
  mupart.corr <- mupart.corr[!is.na(mupart.corr)]
  pre.corr <- pre.corr[!is.na(pre.corr)]
  Icount <- 1
  rho <- NA
  RHO <- NA
  for (i in 1:length(bivar)){
    if (!is.null(bivar[[i]])){
      for (j in 1:nrow(bivar[[i]])){
        long <- length(bivar[[i]][j,])
        sigmainv.corr <- append(sigmainv.corr, paste0("SigmaInv",randvar[i],paste(unique(bivar[[i]][j,]), collapse = ""),
                                                      "[1:",long, ",1:",long, "] ~ dwish(I",
                                                      Icount,"[1:",
                                                      long, ",1:",long, "],", long+1,")","\n"))
        reqsigma <- paste0("Sigma",randvar[i],paste(unique(bivar[[i]][j,]), collapse = ""))
        sigma.corr <- append(sigma.corr, paste0(reqsigma,
                                                "[1:",long, ",1:",long, "] <- inverse(",
                                                "SigmaInv",randvar[i],paste(unique(bivar[[i]][j,]), collapse = ""),
                                                "[1:",long, ",1:",long, "])","\n"))
        rho <- append(rho, paste0("for (i1 in 1:", long, "){\n for(i2 in 1:", long, "){rho",
                                  paste0(randvar[i],paste(unique(bivar[[i]][j,]), collapse = "")),
                                  "[i1,i2] <- ", reqsigma, "[i1,i2]/sqrt(", reqsigma, "[i1,i1] *", reqsigma, "[i2,i2])\n}\n}"))
        RHO <- append(RHO, paste0("rho", paste0(randvar[i],paste(unique(bivar[[i]][j,]), collapse = ""))))
        Icount <- Icount + 1
        wishdf <- c(wishdf, long)
      }
    }
  }
  for (i in 1:length(multivar)){
    if (!is.null(multivar[[i]])){
      for (j in 1:nrow(multivar[[i]])){
        long <- length(multivar[[i]][j,])
        sigmainv.corr <- append(sigmainv.corr, paste0("SigmaInv",randvar[i],paste(unique(multivar[[i]][j,]), collapse = ""),
                                                      "[1:",long, ",1:",long, "] ~ dwish(I", Icount,"[1:",
                                                      long, ",1:",long, "],", long+1,")","\n"))
        reqsigma <- paste0("Sigma",randvar[i],paste(unique(multivar[[i]][j,]), collapse = ""))
        sigma.corr <- append(sigma.corr, paste0(reqsigma,
                                                "[1:",long, ",1:",long, "] <- inverse(",
                                                "SigmaInv",randvar[i],paste(unique(multivar[[i]][j,]), collapse = ""),
                                                "[1:",long, ",1:",long, "])","\n"))
        rho <- append(rho, paste0("for (i1 in 1:", long, "){\n for(i2 in 1:", long, "){rho",
                                  paste0(randvar[i],paste(unique(multivar[[i]][j,]), collapse = "")),
                                  "[i1,i2] <- ", reqsigma, "[i1,i2]/sqrt(", reqsigma, "[i1,i1] *", reqsigma, "[i2,i2])\n}\n}"))
        RHO <- append(RHO, paste0("rho", paste0(randvar[i],paste(unique(multivar[[i]][j,]), collapse = ""))))
        Icount <- Icount + 1
        wishdf <- c(wishdf, long)
      }
    }
  }

  sigmainv.corr <- sigmainv.corr[!is.na(sigmainv.corr)]
  sigma.corr <- sigma.corr[!is.na(sigma.corr)]
  wishdf <- wishdf[!is.na(wishdf)]
  RHO <- RHO[!is.na(RHO)]
  rho <- rho[!is.na(rho)]
  ia.purecont<-vector()
  if (nrcat & nrcont > 0){
    counter <- 0
    for (ncont in 1:(nrcont+nrcat-1)){
      for (ncat in (ncont+1):(nrcont+nrcat)){
        counter <- counter + 1
        # 
        if(sum(startsWith(x = c(allnames[ncont],allnames[ncat]), prefix = substr(allnames[ncont],start = 1, stop = nchar(allnames[ncont])-1)))==2 &
           grepl(pattern=".spl",x=allnames[ncont]) & grepl(pattern=".spl",x=allnames[ncat])){
        }else{
          for (i in 1:nrand){
            if(randvar.ia[[i]][ncat,ncont] == 1){
              bIAs[[counter]][i,1] <- paste("b", allnames[ncont], "x", allnames[ncat], ".", randvar[i], "[", randvar[i],"[i]]",
                                            " * x", allnames[ncont], "[i]"," * x", allnames[ncat], "[i]", sep = "")
              bIAs[[counter]][i,2] <- paste("b", allnames[ncont], "x", allnames[ncat], ".", randvar[i], "[n]", sep = "")
              bIAs[[counter]][i,3] <- paste("b", allnames[ncont], "x", allnames[ncat], ".", randvar[i], sep = "")

              tauIAs[counter,i] <- paste("tau", allnames[ncont], "x", allnames[ncat], ".", randvar[i], sep = "")
              # plus sigma for rescaling prec as sd
              sigmaIAs[counter,i] <- paste("sigma", allnames[ncont], "x", allnames[ncat], ".", randvar[i], sep = "")
            } else {muIAs[[1]][[counter]] <- paste("mu", allnames[ncont],"x", allnames[ncat], " * x",
                                                   allnames[ncont],"[i]", " * x", allnames[ncat],"[i]", sep = "")}

            muIAs[[2]][counter]<- paste("mu", allnames[ncont], "x", allnames[ncat], sep = "")
            # preIAs <- append(preIAs, paste("pre", allnames[ncont], "x", allnames[ncat],  sep = ""))
          }
          preIAs <- append(preIAs, paste("pre", allnames[ncont], "x", allnames[ncat],  sep = ""))
        }
        if(ncont<=nrcont&ncat<=nrcont){
          ia.purecont<-c(ia.purecont,1)
          } else{ia.purecont<-c(ia.purecont,0)}
      }}
    preIAs <- preIAs[2:length(preIAs)]
  } else{muIAs <- rm(muIAs); preIAs <- rm(preIAs); tauIAs <- rm(tauIAs); sigmaIAs <- rm(sigmaIAs)}

  ## assign text for regression formula
  eqparms <- NA
  B <- vector("list", nrand)
  for (i in 1:nrand){
    eqparms <- append(eqparms, bI[[i]][1])
    if (all(multivar[[i]]!=1) & all(bivar[[i]]!=1)){
      B[[i]] <- append(B[[i]], bI[[i]][2])
    }
  }

  # assign random cont effects
  for (i in which(conthcl==1)){
    for (j in which(dat.str[dat.str$type == "cont",][i,3:ncol(dat.str)]==1)){
      eqparms <- append(eqparms, bcont[[i]][j,1])
      if (all(multivar[[j]]!=(1+i)) & all(bivar[[j]]!=(1+i))){
        B[[j]] <- append(B[[j]], bcont[[i]][j,2])
      }}
  }
  # assign fixed cont effects
  for (i in which(conthcl==0)){
    eqparms <- append(eqparms, mucont[[1]][i])
  }

  # same for cat and ia effects
  for (i in which(cathcl==1)){
    for (j in which(dat.str[dat.str$type == "cat",][i,3:ncol(dat.str)]==1)){
      eqparms <- append(eqparms, bcat[[i]][j,1])
      if (all(multivar[[j]]!=(1+sum(conthcl)+i)) & all(bivar[[j]]!=(1+sum(conthcl)+i))){
        B[[j]] <- append(B[[j]], bcat[[i]][j,2])
      }}
  }
  for (i in which(cathcl==0)){
    eqparms <- append(eqparms, mucat[[1]][i])
  }
  # no correlations between random ias
  for (j in 1:nrand){
    for (i in which(randvar.ia[[j]][lower.tri(randvar.ia[[j]])]==1)){
      eqparms <- append(eqparms, bIAs[[i]][j,1])
      B[[j]] <- append(B[[j]], bIAs[[i]][j,2])
    }
  }

  sums <- matrix(0, nrow = nrow(randvar.ia[[1]]), ncol = ncol(randvar.ia[[1]]))
  for (i in 1:nrand){
    sums <- sums + randvar.ia[[i]]
  }
  neitherIA <- which(sums[lower.tri(sums)]==0)
  for (i in neitherIA){
    eqparms <- append(eqparms, muIAs[[1]][i])
  }
  eqparms <- eqparms[!is.na(eqparms)]

  # split up individual effects according to random grouping variable
  TAU <- vector("list", nrand)
  SIGMA <- vector("list", nrand)
  MU <- vector("list", nrand)
  PRE <- vector("list", nrand)

  for (i in 1:nrand){
    TAU[[i]] <- c(tau0g[[i]], taucont[,i], taucat[,i],
                  tauIAs[which(randvar.ia[[i]][lower.tri(randvar.ia[[i]])]==1),i])
    TAU[[i]] <- TAU[[i]][!is.na(TAU[[i]])]
    SIGMA[[i]] <- c(sigma0[i], sigmacont[,i], sigmacat[,i],
                    sigmaIAs[which(randvar.ia[[i]][lower.tri(randvar.ia[[i]])]==1),i])
    SIGMA[[i]] <- SIGMA[[i]][!is.na(SIGMA[[i]])]
  }

  MU<-mu0g

  for (i in which(conthcl==1)){
    for (j in which(dat.str[dat.str$type == "cont",][i,3:ncol(dat.str)]==1)){
      if (all(multivar[[j]]!=(1+i)) & all(bivar[[j]]!=(1+i))){
        MU[[j]] <- append(MU[[j]], mucont[[2]])
      }
    }
  }
  for (i in which(cathcl==1)){
    for (j in which(dat.str[dat.str$type == "cat",][i,3:ncol(dat.str)]==1)){
      if (all(multivar[[j]]!=(1+sum(conthcl)+i)) & all(bivar[[j]]!=(1+sum(conthcl)+i))){
        MU[[j]] <- append(MU[[j]], mucat[[2]])
      }
    }
  }
  for (i in 1:nrand){
    MU[[i]] <- append(MU[[i]], muIAs[[2]][which(randvar.ia[[i]][lower.tri(randvar.ia[[i]])]==1)])
  }
  MU <- lapply(MU, function(x) x[!is.na(x)])
  # part 1 of required model text
  regeq <- matrix()
  for (i in 1:length(eqparms)){
    regeq <- append(regeq, c(eqparms[i], "+","\n"))
  }
  # text modules for jags model text
  likelihood <- paste("y[i] ~ dnorm( mu[i], tau )\n")
  regeq <- paste(regeq[2:(length(regeq)-2)], collapse = '')

  part1 <- paste("for ( i in 1:Ndata ) {\n", likelihood, "mu[i] <- ", regeq,"\n}")

  eff.rand <- vector("list", nrand)
  part2 <- replicate(nrand, 0, simplify=F)


  for (i in 1:length(B)){
    if(!is.null(B[[i]])){
      for (j in 1:length(B[[i]])){
        eff.rand[[i]][j] <- paste(B[[i]][j], "~ dnorm(", MU[[i]][j], ",", TAU[[i]][j],")\n")
      }}
    tmp <- paste(eff.rand[[i]], collapse = ' ')
    part2[[i]] <- paste("for (n in 1:", "N", randvar[i], "s)",
                        sep = '', "{\n",paste(eff.rand[[i]], collapse = ''), collapse = '')
  }
  mucorr <- vector("list", 1)
  for (i in 1:length(MU)){
    mucorr[[i]] <- c(NA, MU[[i]][grepl("corr", MU[[i]])])
    MU[[i]] <- c(NA, MU[[i]][!grepl("corr", MU[[i]])])
  }
  for (i in 1:length(mucont)){
    mucorr <- c(mucorr, mucont[[i]][grepl("corr", mucont[[i]])], mucat[[i]][grepl("corr", mucat[[i]])])
    mucont[[i]] <- c(NA,mucont[[i]][!grepl("corr", mucont[[i]])])
    mucat[[i]] <- c(NA, mucat[[i]][!grepl("corr", mucat[[i]])])
  }

  mucorr <- unique(unlist(lapply(mucorr, function(x) x[!is.na(x)])))
  MU <- lapply(MU, function(x) x[!is.na(x)])

  mucont <- lapply(mucont, function(x) x[!is.na(x)])
  mucat <- lapply(mucat, function(x) x[!is.na(x)])

  PRE <- c(precont, precat, preIAs)
  PRE <- PRE[!is.na(PRE)]
  MU.INT <- c(mu0g, mucont[[2]], mucat[[2]], muIAs[[2]]) #, "mu.corr"

  for (i in 1:length(bipart)){
    if (length(bipart[[i]])>0){
      for (j in 1:length(bipart[[i]])){
        part2[[i]] <- append(part2[[i]], paste(bipart[[i]][j]))
      }}
  }
  for (i in 1:length(multipart)){
    if (length(multipart[[i]])>0){
      for (j in 1:length(multipart[[i]])){
        part2[[i]] <- append(part2[[i]], paste(multipart[[i]][j]))
      }}
    part2[[i]] <- append(part2[[i]], " }")}
  part2 <- paste(append(unlist(part2), "tau ~ dgamma(sg, rg)\n"), collapse = '')

  # now write mu priors
  eff.tauhyp <- NA
  eff.sigmahyp <- NA
  eff.mhyp <- NA
  if (all(unlist(multivar)!=1) & all(unlist(bivar)!=1)){
    eff.mhyp <- paste(MU[[1]][1], "~ dnorm(0,1) \n")
  }
  for (i in 1:nrand){
    if (all(multivar[[i]]!=1) & all(bivar[[i]]!=1)){
      # and tau priors
      eff.tauhyp <- append(eff.tauhyp, paste(TAU[[i]][1], " <- 1/", sigma0[i], "^2\n"))
      eff.sigmahyp <- append(eff.sigmahyp, paste(sigma0[i], " ~ dgamma(1,0.04)\n"))
    }
  }

  if (length(mucont[[2]])>0){
    for (i in 1:length(mucont[[2]])){
      eff.mhyp <- append(eff.mhyp, paste(mucont[[2]][i], " <- ", precont[i], "* scalecont \n"))
    }}

  for (i in 1:nrand){
    for (j in which(dat.str[dat.str$type == "cont",][,2+i]==1)){
      if ((all(multivar[[i]]!=(1+j)) & all(bivar[[i]]!=(1+j)))){
        eff.tauhyp <- append(eff.tauhyp, paste(taucont[j,i], " <- 1/", sigmacont[j,i], "^2\n"))
        eff.sigmahyp <- append(eff.sigmahyp, paste(sigmacont[j,i], " ~ dgamma(1,0.04)\n"))
      }}
  }
  if (length(mucat[[2]])>0){
    for (i in 1:length(mucat[[2]])){
      eff.mhyp <- append(eff.mhyp, paste(mucat[[2]][i], " <- ", precat[i], "* scalecat \n", sep = ""))
    }}

  for (i in 1:nrand){
    for (j in which(dat.str[dat.str$type == "cat",][,2+i]==1)){
      if ((all(multivar[[i]]!=(1+nrcont+j)) & all(bivar[[i]]!=(1+nrcont+j)))){
        eff.tauhyp <- append(eff.tauhyp, paste(taucat[j,i], " <- 1/", sigmacat[j,i], "^2\n", sep = ""))
        eff.sigmahyp <- append(eff.sigmahyp, paste(sigmacat[j,i], " ~ dgamma(1,0.04)\n"))
      }}
  }
  options(warn = -1)
  muIAs[[2]]<-muIAs[[2]][!is.na(muIAs[[2]])]
  options(warn=0)
  if (length(muIAs[[2]])>0){
    for (i in 1:length(muIAs[[2]])){
      if(ia.purecont[i]==0){
        eff.mhyp <- append(eff.mhyp, paste(muIAs[[2]][i], " <- ", preIAs[i], "* scalecat \n", sep = ""))
      } else if(ia.purecont[i]==1){
        eff.mhyp <- append(eff.mhyp, paste(muIAs[[2]][i], " <- ", preIAs[i], "* scalecont \n", sep = ""))
        }
    }}
  for (j in 1:nrand){
    for (k in which(randvar.ia[[j]][lower.tri(randvar.ia[[j]])]==1)){
      eff.tauhyp <- append(eff.tauhyp, paste(tauIAs[k,j], " <- 1/", sigmaIAs[k,j], "^2\n", sep = ""))
      eff.sigmahyp <- append(eff.sigmahyp, paste(sigmaIAs[k,j], " ~ dgamma(1,0.04)\n"))
    }
  }
  if (is.na(eff.tauhyp[1])){
    eff.tauhyp <- eff.tauhyp[2:length(eff.tauhyp)]
  }
  eff.tauhyp <- eff.tauhyp[!is.na(eff.tauhyp)]
  if (is.na(eff.sigmahyp[1])){
    eff.sigmahyp <- eff.sigmahyp[2:length(eff.sigmahyp)]
  }
  eff.sigmahyp <- eff.sigmahyp[!is.na(eff.sigmahyp)]
  if (is.na(eff.mhyp[1])){
    eff.mhyp <- eff.mhyp[2:length(eff.mhyp)]
  }
  eff.mhyp <- eff.mhyp[!is.na(eff.mhyp)]
  prior.pre <- matrix()
  for (i in 1:length(PRE)){
    prior.pre[i] <- paste(PRE[i], "~ dt(0,1,1)\n")
  }
  scalecat <- "scalecat <- 1/2\n"
  scalecont <- "scalecont <- sqrt(2)/4\n"
  part3 <- paste(paste(eff.mhyp, collapse = ' '), paste(eff.tauhyp, collapse = ' '),
                 paste(eff.sigmahyp, collapse = ' '),
                 paste(prior.pre, collapse = ' '), paste(scalecat, collapse = ' '), paste(scalecont, collapse = ' '),
                 paste(mupart.corr, collapse = ' '), paste(pre.corr, collapse = ' '),
                 paste(sigmainv.corr, collapse = ' '), paste(sigma.corr, collapse = ' '),
                 paste("sg <- pow(m,2)/pow(d,2) \nrg <- m/pow(d,2) \nm ~ dgamma(1,1) \nd ~ dgamma(1,1)\n"))
  part4 <- paste(rho, collapse = '')
  # complete model
  modelstring = paste(" model {", part1, part2, part3, part4, "}")
  writeLines(modelstring,con=path)
  bI.save <- matrix(unlist(bI),ncol=3,byrow=TRUE)[,3]
  bcont.save <- unlist(lapply(bcont, function(x) x[,3]))
  bcat.save <- unlist(lapply(bcat, function(x) x[,3]))
  bIAs.save <- unlist(lapply(bIAs, function(x) x[,3]))
  options(warn = -1)
  mcmcsave <- unique(c(bI.save[!is.na(bI.save)], bcont.save[!is.na(bcont.save)],
                       bcat.save[!is.na(bcat.save)], bIAs.save[!is.na(bIAs.save)]))
  options(warn=0)
  parameters <- list(MU = MU.INT, SIGMA = unlist(SIGMA), b.save = mcmcsave)
  parameters[["pl.ind"]] <- pl.ind
  parameters[["pl.nhclcont"]] <- mucont[[2]][nhclcont]
  parameters[["pl.nhclcat"]] <- mucat[[2]][nhclcat]
  if (!is.na(RHO[1])){
    parameters[["RHO"]] <- unlist(RHO)
    parameters[["mu.corr"]] <- mucorr
    parameters[["wishdf"]] <- wishdf
    parameters[["corrnames"]]<-corr.names
  }
  return(parameters)
}
