#All helpers needed to run MOB-ERA-GEE

Xscale <- function(x){
  
  # This function returns X that cor(X) == (t(X)%*%X)
  
  # Below calculation is equal to "norm(, type = "2")" function in R, i.e.,
  # temp <- x
  # for (j in 1:ncol(temp)) { # 0 mean and unit variance (2-norm)
  #	temp[,j] <- temp[,j] / norm(temp[,j,drop=FALSE], type = "2")
  # }
  
  ctz <- scale(x, center=TRUE, scale=FALSE)
  ctz <- as.matrix(ctz)
  covz <- (t(ctz)%*%ctz)/nrow(x)
  Dstz <- sqrt(diag(diag(covz)))
  bz0 <- t(solve(t(Dstz),t(ctz)))
  bz0 <- bz0/sqrt(nrow(x))
  
  return(bz0)
}

linkfun <- function(x, lp, y, weights = NULL, offset = NULL, family = family,
                    mustart = NULL, etastart = NULL, start = NULL)
{
  # Ref: stats::glm() in R
  # In R, family objects: https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/glm.R
  # IWSL with adj DVs: https://www.statistics.ma.tum.de/fileadmin/w00bdb/www/czado/lec2.pdf
  # From linkfun(): adj.DV$w (weights in IWLS) is all 1 when gaussian
  
  x <- as.matrix(x)
  nvars <- ncol(x)
  nobs <- nrow(y)
  if (is.null(weights)) {	weights <- rep.int(1, nobs) }
  if (is.null(offset)) { offset <- rep.int(0, nobs) }
  
  ## Define family:	
  if(is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  
  if(is.function(family)) family <- family()
  
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  ## Get family functions:
  variance <- family$variance
  linkinv  <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv) ) {
    stop("'family' argument seems not to be a valid family object", call. = FALSE)
  }
  mu.eta <- family$mu.eta
  
  unless.null <- function(x, if.null) if(is.null(x)) if.null else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu  <- unless.null(family$validmu,  function(mu) TRUE)
  
  if(is.null(mustart)) {
    ## calculates mustart and may change y and weights and set n (!)
    eval(family$initialize)
  } else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  
  ## Get (initial) eta
  coefold <- NULL
  eta <- lp
  mu <- linkinv(eta)
  mvar <- variance(mu) # the variance as a function of the mean
  
  if (!(validmu(mu) && valideta(eta))) {
    stop("cannot find valid starting values: please specify some", call. = FALSE)
  }
  
  ## Calculate z and w
  good <- weights > 0
  mu.eta.val <- mu.eta(eta)
  if (any(is.na(mu.eta.val[good]))) {stop("NAs in d(mu)/d(eta)")}
  
  # drop observations for which w will be zero
  good <- (weights > 0) & (mu.eta.val != 0)
  z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
  # z: adjusted dependent variable
  w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
  # w: weights
  
  # BETA in GLM: Fisher Scoring
  # beta <- solve(t(X) %*% diag(w) %*% X) %*% (t(X) %*% diag(w) %*% z)
  
  # Outs
  return( list(mu=mu, mvar=mvar, z=z, w=w) )
}

getRi <- function(y, mu, mvar, nt, corstr, Mv = NULL, R = NULL, scale.fix = FALSE, scale.value)
{
  # Reference: Gul Inan, Jianhui Zhou and Lan Wang (2017). PGEE. R package version 1.5.
  # https://CRAN.R-project.org/package=PGEE
  
  N <- nrow(mu) # nobs
  maxclsz <- max(nt)
  sd <- sqrt(mvar)
  res <-(as.vector(y)-mu)/sd 
  
  if(scale.fix==0) {
    phi<-sum(res^2)/(sum(nt))
  } else
    if(scale.fix==1) {
      phi<-scale.value
    }
  
  aindex=cumsum(nt)
  index=c(0,aindex[-length(aindex)])
  
  if (corstr=="independence") 
  {alfa_hat<-0} else 
    if (corstr=="exchangeable") 
    {
      sum1<-0
      sum3<-0
      for ( i in  1:N)          {
        for ( j in  1:nt[i])      {    
          for ( jj in 1:nt[i])      {    
            if  ( j!=jj)              {
              #cat("i",i,"j",j,"jj",jj,"\n")
              sum2<-res[j+index[i]]*res[jj+index[i]]
              #cat("i",i,"j",j,"jj",jj,"\n")
              sum1<-sum1+sum2
              #cat("i",i,"j",j,"jj",jj,"sum2",sum2,"sum1",sum1,"\n")
            }
          }
        }
        sum4<-nt[i]*(nt[i]-1)
        sum3<-sum3+sum4
      } #i
      alfa_hat<-sum1/(sum3*phi)
    } else
      if (corstr=="AR-1") 
      { 
        sum5<-0
        sum6<-0
        for ( i in  1:N)           {
          for ( j in  1:nt[i])       {  
            for ( jj in 1:nt[i])       {  
              if( j>jj && abs(j-jj)==1)  {
                #cat("i",i,"j",j,"jj",jj,"\n")
                sum7<-res[j+index[i]]*res[jj+index[i]]
                sum5<-sum5+sum7           
                #cat("i",i,"j",j,"jj",jj,"sum7",sum7,"sum5", sum5, "\n")
              }
            }
          }
          sum8<-(nt[i]-1)
          sum6<-sum6+sum8
        } #i
        alfa_hat<-sum5/(sum6*phi)
      } else
        if (corstr=="stat_M_dep") 
        {  
          alfa_hat=matrix(0,Mv,1)
          for(m in 1:Mv) {
            sum12<-0
            sum14<-0
            for ( i in  1:N)           {
              for ( j in  1:nt[i])       {  
                for ( jj in 1:nt[i])       {  
                  if( j>jj && abs(j-jj)==m)  {
                    #cat("m",m,"i",i,","j",j,"jj",jj,"\n") 
                    sum11<-res[j+index[i]]*res[jj+index[i]]
                    sum12<-sum12+sum11 
                    #cat("m",m,"i",i,"j",j,"jj",jj,"sum11",sum11,"sum12", sum12, "\n")          
                  } #if
                }
              }
              sum13<-nt[i]-1
              sum14<-sum14+sum13
            } #i
            alfa_hat[m]<-sum12/(sum14*phi) 
          } #m
        }  else
          if (corstr=="non_stat_M_dep") 
          {  
            alfa_hat<-matrix(0,nt[1],nt[1]) #not allowed for unequal number of cluster sizes.
            for( m in 1:Mv)            {
              for ( j in  1:nt[1])       {  
                for ( jj in 1:nt[1])       {  
                  if( j>jj && abs(j-jj)==m)  { 
                    sum16<-0                 
                    for ( i  in 1:N)     {
                      #cat("m",m,"j",j,"jj",jj,"i",i"\n") 
                      sum15<-res[j+index[i]]*res[jj+index[i]]
                      sum16<-sum15+sum16          
                      #cat("j",j,"jj",jj,"i",i,"sum15",sum15,"sum16",sum16,"\n")
                    } #i
                    #cat("j",j,"jj",jj,"sum16",sum16,"\n")
                    alfa_hat[j,jj]<-sum16/(N*phi) 
                  }
                }
              }
            }
          } else
            if (corstr=="unstructured") 
            {  
              alfa_hat<-matrix(0,nt[1],nt[1]) #not allowed for unequal number of cluster sizes.
              for ( j in 1:nt[1])  {  
                for ( jj in 1:nt[1]) {
                  sum20<-0                
                  if (j > jj)          {
                    for ( i  in 1:N )    {
                      #cat("i",i,"j",j,"jj",jj,"\n") 
                      sum21<-res[j+index[i]]*res[jj+index[i]]
                      sum20<-sum21+sum20           
                    } #i
                    #cat("j",j,"jj",jj,"sum20",sum20,"\n")
                    alfa_hat[j,jj]<-sum20/(N*phi) 
                  }
                }
              }
            } else
              if (corstr=="fixed")
              {alfa_hat=NULL
              }
  
  cor1 <- matrix(0, maxclsz, maxclsz)
  
  if (corstr=="independence")                                        
  {cor1<-diag(nt[1])} else
    if (corstr=="exchangeable")                                        
    { for (t1 in 1:nt[i]) {
      for (t2 in 1:nt[i]) {
        if (t1!=t2) 
        {cor1[t1,t2]<-alfa_hat} else 
        {cor1[t1,t2]<-1}
      }
    }
    } else
      if (corstr=="AR-1")                                      
      { for (t1 in 1:nt[i]) {
        for (t2 in 1:nt[i]) {
          cor1[t1,t2]<-alfa_hat^abs(t1-t2)   
        }
      }
      } else
        if (corstr=="stat_M_dep")                                     
        { for (t1 in 1:nt[i]) {
          for (t2 in 1:nt[i]) {
            if (abs(t1-t2)==0)
            {cor1[t1,t2]<-1} else
              for(m in 1:Mv) {
                if (abs(t1-t2)==m)
                {cor1[t1,t2]<-alfa_hat[m]} 
              }
          }
        }
        } else
          if (corstr=="non_stat_M_dep")                                     
          { 
            cor1=alfa_hat+t(alfa_hat)
            diag(cor1)=1
          } else
            if (corstr=="unstructured")                                        
            {cor1=alfa_hat+t(alfa_hat)
            diag(cor1)=1
            } else
              if (corstr=="fixed")
              {cor1=R
              }
  
  return( list(Ehat=cor1, phi=phi) )
  
}

getAVAi <- function(mvar, nt, Rhat, phi)
{	
  # Reference: Gul Inan, Jianhui Zhou and Lan Wang (2017). PGEE. R package version 1.5.
  # https://CRAN.R-project.org/package=PGEE
  
  N <- nrow(mvar) # nobs
  maxclsz <- max(nt)
  aindex <- cumsum(nt)
  index <- c(0,aindex[-length(aindex)])
  
  AVAhat <- array(0, c(maxclsz,maxclsz,N))
  
  for (i in 1:N) {
    
    bigA <- matrix(0,nt[i],nt[i])
    for (j in 1:nt[i]) {
      bigA[j,j] <- mvar[j+index[i]]
    } # for j
    
    ##Inverse of working covariance matrix
    bigV <- sqrt(bigA)%*%Rhat%*%sqrt(bigA)
    invV <- MASS::ginv(bigV)
    #bigV <- phi*invV
    
    AVAhat[1:nt[i],1:nt[i],i] <- invV
  }
  
  return(AVAhat)
}

Errors <- function(x, family){
  
  ## Define family:	
  if(is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  
  if(is.function(family)) family <- family()
  
  ## Get family functions:
  devresid <- family$dev.resid
  linkinv  <- family$linkinv
  
  ## Predicted values
  Pred <- x$X %*% x$W %*% x$B
  di <- devresid(x$Y, linkinv(Pred), wt=1)
  error <- sum(di)
  
  return( list(Error = error) )
  invisible(x)
  # return( list(Dev.resids = di, Error = error) )
}

#-----------------------------
# Calculate CV errors
#-----------------------------
# dev.resids: function giving the deviance residuals as a function of (y, mu, wt)
# family$dev.resids gives the square of the residuals, i.e., (di)^2

CV_Errors <- function(train, test.X, test.Y, family){
  
  ## Define family:	
  if(is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  
  if(is.function(family)) family <- family()
  
  ## Get family functions:
  devresid <- family$dev.resid
  linkinv  <- family$linkinv
  
  ## Predicted values
  Pred <- test.X %*% train$W %*% train$B
  di <- devresid(test.Y, linkinv(Pred), wt=1)
  error <- sum(di)
  
  return( list(CV_Error = error) )
  invisible(x)
  # return( list(Dev.resids = di, Error = error) )
}


MOB.MERA <-  function(y, x, family, corstr, nlv, Pk, lambdaW=0, lambdaB=0, maxit=100, ceps=1e-5, seed = NULL, ...){
  
  if(!is.null(seed)) set.seed(seed) #NEW: setting seed to ensure reproducible resutls in empirical application
  
  # ------------------------------------ #
  # ------ Data inputs & Warnings ------ #
  # ------------------------------------ #
  if(missing(corstr)) corstr="independence"
  if(missing(nlv)) stop("The number of components is not specified")
  if(missing(Pk)) stop("Provide the number of observed predictors (for each component set) in a matrix form")
  if(missing(lambdaW)) lambdaW = 0
  if(missing(lambdaB)) lambdaB = 0
  if(missing(maxit)) maxit = 100		# max no. of iterations
  if(missing(ceps)) ceps = 1e-5		# convergence tolerance
  
  Yfam <- family
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  # No. of repeated measures and samples
  nobs = nrow(y)
  nrepeat = ncol(y)
  n_tot = nobs*nrepeat
  
  # No. repeated measures in each id (used in getcorr() func.)
  nt = as.integer(rep(nrepeat, nobs))
  maxclsz = max(nt)
  
  # Is the form of "family" = formula ?
  ncy = ncol(y)
  if(length(Yfam)!=ncy) stop("The number of dep. variables and family attributes are different!")
  
  # ------------------------------------------------------ #
  # ------ Initialization: Calculate initial values ------ #
  # ------------------------------------------------------ #
  
  # Data Normalization
  X <- Xscale(x)						# zscore(data)/sqrt(ncase-1) (in Matlab)
  colnames(X) <- xnames
  if( any(Yfam == "gaussian") ) {
    y <- Xscale(y)
    colnames(y) <- ynames
  }
  
  # W: Random values from Uniform Dist.
  W00 <- matrix(0, nrow=sum(Pk), ncol=nlv)
  ndset <- ncol(Pk)
  kk = 0
  for (j in 1:ndset){
    k = kk + 1
    kk = kk + Pk[,j]
    if ( Pk[,j]==1 ) { W00[k:kk,j] = 1
    } else { W00[k:kk,j] = 99
    }
  }
  windex <- which(W00 == 99)
  num_windex = length(windex)
  W0 = W00
  W0[windex] <- runif(num_windex)
  
  # F & B: Multivariate linear reg.
  F00 <- X%*%W0
  F0 <- Xscale(F00)	# Normalized (i.e., F0'F0=corr)
  F0_d <- stats::setNames( data.frame(F0), paste( rep("F", ncol(F0)), 1:ncol(F0), sep="" ) )
  #F0_1 <- model.matrix(~., F0_d, contrasts)
  F0_1 <- model.matrix(~., F0_d)
  if(colnames(F0_1)[1]== "(Intercept)") F0 <- F0_1[,-1]
  B0 <- MASS::ginv( t(F0)%*%F0 ) %*% t(F0) %*% y	# (F'*F)\F'*y (in Matlab)
  
  # eta: linear predictor
  lp <- F0%*%B0
  
  # GLM outputs by linkfun()
  # 1) Z: adjusted DVs (column-by-column)
  # 2) mu: linkinv(eta)
  # 3) mvar: variance(mu)
  Z <- mu <- mvar <- matrix( , nrow=nobs, ncol=nrepeat)
  for (r in 1:nrepeat) {		
    adj.DV <- linkfun(X, lp[,r,drop=FALSE], y[,r,drop=FALSE], family = Yfam[r])
    Z[,r] <- adj.DV$z
    mu[,r] <- adj.DV$mu
    mvar[,r] <- adj.DV$mvar
  }
  colnames(Z) <- ynames
  
  # Vi, phi, & Ai: by getRi() & getAVAi
  R_phi_hat <- getRi(y, mu, mvar, nt, corstr)
  AVAi_hat <- getAVAi(mvar, nt, R_phi_hat$Ehat, R_phi_hat$phi)
  
  # ------------------------------------------ #
  # ------ Iterative method starts here ------ #
  # ------------------------------------------ #
  
  # Start values
  oldZ <- Z
  oldAVAi <- AVAi_hat
  oldW <- W0
  oldF <- F0
  oldB <- data.frame(B0)
  
  iter <- 0L
  
  for (iter in 1L:maxit) {
    
    iter <- iter + 1L
    
    if ( iter > maxit ) {
      cat("Message: not converged in", maxit, "iterations (convergence criterion =", ceps, ")", "\n")
      break
    }
    
    #-----------------
    # Step1: Update B
    #----------------- 'oldB' with F0, Qi, lambdaB
    
    vecb <- matrixcalc::vec(t(oldB))
    I_db <- dim(t(oldB))[1] 				# Dimension of I for Qi: Harville, 1997. p.342, (2.11)
    
    fB <- c()								# Obj func for B
    Bterm1 <- Bterm2 <- list()				# For each i, to sum up N matrices later
    
    for ( i in 1L:nrow(oldF) ){
      
      Zi <- oldZ[i,,drop=FALSE]
      Zi <- matrixcalc::vec(Zi)
      
      Fi <- oldF[i,,drop=FALSE]
      Fi <- matrixcalc::vec(Fi)
      
      Qi <- kronecker( t(Fi), diag(I_db) )
      
      Vi_inv <- oldAVAi[1:nt[i],1:nt[i],i]
      
      fB[i] <- t(Zi - Qi%*%vecb) %*% Vi_inv %*% (Zi - Qi%*%vecb)
      Bterm1[[i]] <- t(Qi) %*% Vi_inv %*% Qi
      Bterm2[[i]] <- t(Qi) %*% Vi_inv %*% Zi
      
    }
    
    fB_tot <- sum(fB) + (lambdaB * base::crossprod(vecb))
    Bterm1_tot <- Reduce('+', Bterm1) + (lambdaB * diag( dim(Bterm1[[1]])[1] ))
    Bterm2_tot <- Reduce('+', Bterm2)
    vecB_tot <- MASS::ginv(Bterm1_tot) %*% Bterm2_tot
    
    # Updated estimates: vecB_tot (new "vecb")
    newB <- ks::invvec(vecB_tot, nrow=nrow(oldB), ncol=ncol(oldB), byrow = TRUE)		
    
    #-----------------
    # Step2: Update W
    #----------------- 'oldW' with X, newB, Mi, lambdaW
    
    vecW0 <- matrixcalc::vec(t(oldW))
    
    WXindex <- which(vecW0 != 0)			# Eliminating zeros: Hwang & Takane, 2014. p.72 
    vecW <- vecW0[WXindex,,drop=FALSE]
    
    fW <- c()								# Obj func for W
    Wterm1 <- Wterm2 <- list()				# For each i, to sum up N matrices later
    
    for (i in 1L:nrow(X)){
      
      Zi <- oldZ[i,,drop=FALSE]
      Zi <- matrixcalc::vec(Zi)
      
      WXi <- X[i,,drop=FALSE]
      Mi0 <- kronecker( WXi, t(newB) )
      Mi <- Mi0[,WXindex,drop=FALSE] 		# Eliminating zeros
      
      Vi_inv <- oldAVAi[1:nt[i],1:nt[i],i]
      
      fW[i] <- t(Zi - Mi%*%vecW) %*% Vi_inv %*% (Zi - Mi%*%vecW)
      Wterm1[[i]] <- t(Mi) %*% Vi_inv %*% Mi
      Wterm2[[i]] <- t(Mi) %*% Vi_inv %*% Zi
      
    }
    
    fW_tot <- sum(fW) + (lambdaW * base::crossprod(vecW))
    Wterm1_tot <- Reduce('+', Wterm1) + (lambdaW * diag( dim(Wterm1[[1]])[1] ))
    Wterm2_tot <- Reduce('+', Wterm2)
    vecW_tot <- MASS::ginv(Wterm1_tot) %*% Wterm2_tot
    
    ## Update...
    # 1) W: vecW_tot (new "vecW")
    oldW[which(oldW!=0)] <- vecW_tot
    newW <- oldW
    
    # 2) F, and normalize it
    newF0 <- X%*%newW
    newF <- Xscale(newF0)
    colnames(newF) <- colnames(F0)
    
    # 3) eta: linear predictor
    new_lp <- newF%*%newB
    
    # 4) Z: adj DVs (column-by-column)
    new_Z <- new_mu <- new_mvar <- matrix( , nrow=nobs, ncol=nrepeat)
    for (r in 1:nrepeat) {
      adj.DV <- linkfun(X, new_lp[,r,drop=FALSE], y[,r,drop=FALSE], family = Yfam[r])
      new_Z[,r] <- adj.DV$z
      new_mu[,r] <- adj.DV$mu
      new_mvar[,r] <- adj.DV$mvar
    }
    colnames(new_Z) <- ynames
    
    #-----------------
    # Step3: Update Vi
    #-----------------
    
    new_R_phi <- getRi(y, new_mu, new_mvar, nt, corstr)
    new_AVAi <- getAVAi(new_mvar, nt, new_R_phi$Ehat, new_R_phi$phi)
    
    #-----------------
    # check for convergence
    #-----------------
    
    est_old <- rbind(vecb, vecW)
    est_new <- rbind(vecB_tot, vecW_tot)
    
    if ( sum(abs(est_new - est_old)) < ceps ) {
      
      conv <- TRUE
      cat("Message: converged in", iter, "iterations (convergence criterion =", ceps, ")", "\n")
      break
      
    } else {
      
      oldVi <- new_R_phi$Ehat
      oldZ <- new_Z
      
      oldB <- newB
      oldF <- newF
      
      oldW <- newW	
      
      oldAVAi <- new_AVAi
    }
    
  }
  
  # -------------------- #
  # ------ Errors ------ #
  # -------------------- # Deviance residuals
  
  # Define family:	
  # if(is.character(family)) { family <- get(family, mode = "function", envir = parent.frame()) }
  # if(is.function(family)) family <- family()
  
  # Get family functions:
  # devresid <- family$dev.resid
  # linkinv  <- family$linkinv
  
  # Predicted values
  # Pred <- X %*% newW %*% newB
  # error <- sum( devresid(y, linkinv(Pred), wt=1) )
  
  # --------------------- #
  # ------ Outputs ------ #
  # --------------------- #
  
  colnames(newB) <- colnames(B0)
  rownames(newB) <- colnames(F0)
  colnames(newW) <- colnames(F0)
  rownames(newW) <- colnames(x)
  
  corStr <- corstr
  new_Vi <- new_R_phi$Ehat
  Scale_phi <- new_R_phi$phi
  
#########NEW ADDITIONS FOR MOB ########
  # 1) vec(B)
  vecB <- matrixcalc::vec(t(newB))  # vec(B) in column-major order
  
  # 2) Node objective value (unpenalized loss as long as lambda values are = 0)
  objfun <- fB_tot
  
  # 3) Empirical score contributions (with respect to vec(B))
  estfun <- matrix(NA, nrow=nobs, ncol=length(vecB))
  for (i in 1:nobs) {
    fi <- newF[i,,drop=FALSE]           # 1 x nlv
    ri <- matrix(new_Z[i,,drop=FALSE] - fi %*% newB, ncol=1)  # Q x 1
    Vi_inv <- new_AVAi[1:nt[i],1:nt[i], i] # Q x Q
    # Kronecker product: f_i ⊗ I_Q
    estfun[i,] <- -2 * t(kronecker(fi, diag(nrepeat))) %*% Vi_inv %*% ri
  }
  
  
  output.MOB.MERA <- list(Yfamily=Yfam, X=X, Y=y, Z=Z, # Original values
                          fB_tot=fB_tot, fW_tot=fW_tot, # Error=error,
                          lambdaW = lambdaW, lambdaB = lambdaB,
                          W=newW, B=newB, Scale_phi=Scale_phi, corStr=corStr,
                          Vi=new_Vi, coef = vecB,           # coef vector 
                          objfun = objfun,       # node objective
                          estfun = estfun)         # subject-level score contributions
  output.MOB.MERA
  
}


mob_mera_fit <- function(dat.full, ynames, xnames, Yfam, corstr, nlv, Pk, seed = 123) { #outer function: purpose is to return a funciton that MOB can use at each node
  function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
           estfun = FALSE, object = FALSE, ...) { #inner function: actual node fitting function that mob() calls out at every node
    
    if (is.null(weights)) weights <- rep(1L, nrow(dat.full)) # If no weights provided, assume all rows belong to the node
    
    node_rows <- which(weights > 0L) # Identify which rows are included in the current node
    
    if (length(node_rows) == 0L) stop("No observations in node") # Stop if there are no observations in this node
    
    # Subset the full data set to just the observations in this node and relevant variables
    dat_node <- dat.full[node_rows, c(ynames, xnames), drop = FALSE]
    
    Y_mat <- as.matrix(dat_node[, ynames, drop=FALSE]) #creating outcome and predicor matrices
    X_mat <- as.matrix(dat_node[, xnames, drop=FALSE])
    
    # Fit the custom ERA-GEE model to this node's data
    mob_mera_model <- MOB.MERA(y = Y_mat, x = X_mat, family = Yfam, corstr = corstr, nlv = nlv,
                               Pk = Pk, lambdaW = 0, lambdaB = 0, maxit = 100, ceps = 1e-5, seed = seed)
    
    rval <- list(
      coefficients = as.numeric(mob_mera_model$coef),
      objfun = as.numeric(mob_mera_model$objfun)
    )
    if(estfun) rval$estfun <- mob_mera_model$estfun
    if(object) rval$object <- mob_mera_model
    rval}}
