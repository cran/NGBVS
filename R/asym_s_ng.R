asym_s_ng <- function(y, data,  medstar = 1, numb = 100, burnin = 1, every = 1 ) {

  x <- Rfast::standardise(as.matrix(data))
  mod <- glm( y ~ x, binomial )
  betas <- as.vector( mod$coefficients )
  covs <- vcov(mod)
  solvecovs <- solve( covs, tol = 1e-30 )
 # library(corpcor) ## for the pseudoinverse matrix
  lambdastar=1
  gammasq=1
  #burnin=2000
  #numb=20000
  #every=50
  dm <- dim(x)
  n <- dm[1]   ;   p <- dm[2]
  lambda <- rep( lambdastar, p )
  lambdaaccept <- 0
  lambdacount <- 0
#  alpha=rnorm(1,0,0.1)
#  beta= rnorm(p,0,0.1)
  #medstar=1
  psi <- 2 * lambda * gammasq * 0.01
  numbofits <- burnin + every * numb ## total number of iterations
  holdpsi <- matrix(0, nrow = p, ncol = numb )
  holdbeta <- matrix(0, nrow = p, ncol = numb )
  holdgammasq <- numeric( numb )
  holdlambda <- numeric( numb )
  holdalpha <- numeric( numb )
  lambdasd=0.01
  const <- solvecovs %*% betas
  ##########################################
  for ( i in 1:numbofits ) {
    LAM <- diag( c( 0, 1/psi ) )
    #LAM[LAM>1e+7]=1e+7
    #LAM[LAM<1e-7 & LAM!=0]=1e-7
    covsLAM <- solvecovs + LAM
    #covsLAM[covsLAM <1e-7 & covsLAM!=0]=1e-7
    varstar <- solve( covsLAM, tol = 1e-30 )
    expec <- varstar %*% const
    cholstar <- chol( varstar, pivot = TRUE )
    if ( attr( cholstar, "rank" ) == ncol( cholstar ) ) {
      pivot <- attr( cholstar, "pivot" )
      cholstar <- cholstar[ ,order( pivot ) ]
      randn <- rnorm( p + 1 )
      be <- expec + crossprod(cholstar, randn )     ## expec + t( cholstar ) %*% randn
      alpha <- be[1]
      beta <- be[-1]
	}
    ##########################################
    for ( j in 1:p )  {
      ### main check 1
      if ( beta[j]^2 < 10^(-5) ) {
        check <- 0
        ### subcheck 1 for lambda
        if ( lambda[j] < 0.5 ) {
          while ( check == 0 ) {
            psi[j] <- 1 / rgamma( 1, 0.5 - lambda[j], scale = 2 / beta[j]^2 )
            u <- runif(1)
            check <- as.numeric( u < exp( -0.5 * psi[j] / gammasq ) ) }
        } else {
          while (check == 0) {
            psi[j] <- rgamma( 1, lambda[j] - 0.5, scale = 2 * gammasq)
            u <- runif(1)
            check <- as.numeric( u < exp(-0.5 * beta[j]^2 / psi[j] ) ) } }
      } else psi[j] <- rgig( 1, lambda[j] -0.5, beta[j]^2, 1 / gammasq )
	  } ## end for (j in 1:p)
    #psi[psi>1e+7]=1e+7
    psi[psi < 1e-10 & psi != 0 ] <- 1e-10
    ##########################################
    mupsi <- 2 * lambdastar * gammasq
    newlambdastar <- lambdastar * exp( lambdasd * rnorm( 1 ) )
    newgammasq <- mupsi / ( 2 * newlambdastar )
    newlambda <- rep( newlambdastar, p )
    logaccept<- log( newlambdastar ) - log( lambdastar ) - 142.85 * ( newlambdastar - lambdastar )
    logaccept <- logaccept - p * newlambdastar * log( 2 * newgammasq ) - p * lgamma( newlambdastar )
    logaccept <- logaccept + p * lambdastar * log( 2 * gammasq ) + p * lgamma( lambdastar )
    logaccept <- logaccept + newlambdastar * sum( log( psi ) ) - sum( psi ) / ( 2 * newgammasq )
    logaccept <- logaccept - lambdastar * sum( log( psi ) ) + sum( psi ) / ( 2 * gammasq )
    ##########################################
    accept <- 1
    if ( logaccept < 0 )   accept <- exp( logaccept )
    lambdasd <- lambdasd + ( accept - 0.3 ) / i
    lambdaaccept <- lambdaaccept + accept
    lambdacount <- lambdacount + 1
    u <- runif( 1 )
    if ( u < accept ) {
      lambda <- newlambda
      lambdastar <- newlambdastar
      gammasq <- newgammasq
    }
    ###################
    sha <- 0.5 * sum( psi ) + medstar / ( 2 * lambdastar )
    gammasq <- 1 / rgamma( 1, sum( lambda ) + 2, scale = 1 / sha )
    ########################################
    if ( i > burnin & ( i - burnin ) %% every == 0 ) {
      holdlambda[ ( i - burnin ) / every ] <- lambdastar
      holdgammasq[ ( i - burnin ) / every ] <- gammasq
      holdalpha[ ( i - burnin ) / every ] <- alpha
      holdbeta[ , ( i - burnin ) / every ] <- beta
      holdpsi[ , ( i - burnin ) / every ] <- psi
	}

  } ## iterations are over  \ \  end  for ( i in 1:numbofits )

 # return(holdalpha)
#  return(holdbeta)
#  return(holdpsi)
#  return(holdlambda)
 # return(holdgammasq)

  # write.table(holdalpha,paste("alpha",".txt",sep=""), row.names=F,col.names=F,sep=" ", quote=F)
  #write.table(holdbeta,paste("beta",".txt",sep=""), row.names=F,col.names=F,sep=" ", quote=F)
  #write.table(holdpsi,paste("psi",".txt",sep=""), row.names=F,col.names=F,sep=" ", quote=F)
  #write.table(holdgammasq,paste("gammasq",".txt",sep=""), row.names=F,col.names=F,sep=" ", quote=F)
  #write.table(holdlambda,paste("lambda",".txt",sep=""), row.names=F,col.names=F,sep=" ", quote=F)
  list( alpha = holdalpha, beta = holdbeta, psi = holdpsi, lambda = holdlambda, gammasq = holdgammasq )

}
