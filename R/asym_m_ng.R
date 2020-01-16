asym_m_ng <- function(y, data, FS, medstar = c(0.01,0.0001), numb = 100, burnin = 1, every = 1  ) {

  g1 <- which( FS > 0.5 )    ##### Group1:  FS > 0.5
  g2 <- which( FS == 0.5 )   ##### Group2:  FS = 0.5
  g3 <- which( FS < 0.5 )    ##### Group3:  FS < 0.5
  g4 <- which( is.na( FS ) ) ##### Group4:  FS  = NA

  group <- list( g1, g2, g3, g4 )   ##### Create list for the 4 Groups
  n.g <- length( group )
  len_group <- rep( 0,  n.g )
  for ( i in 1:n.g ) len_group[i] <- length( group[[i]] )     #### number of SNPs in each Group


  x <- Rfast::standardise( as.matrix( data ) )
  mod <- glm( y ~ x, binomial )
  betas <- as.vector( mod$coefficients )
  covs <- vcov( mod )
  solvecovs <- solve( covs, tol = 1e-30 )

  dm <- dim(x)
  n <- dm[1]   ;   p <- dm[2]

  newlambdastar <- lambdastar <- lambdasd <- lambdaaccept <- lambdacount <-  accept <- logaccept <- sha <- gammasq <-  newgammasq <- rep( 1, n.g )

  newlambda <- list(999)
  lambda <- list( rep( lambdastar[1], len_group[1] ), rep( lambdastar[2], len_group[2] ),
                rep( lambdastar[3], len_group[3] ), rep( lambdastar[4], len_group[4] ) )

  LAMBDA <- GAMMASQ <- rep( 1, p )
  w <- 1
  h <- 0
  psi <- rep( 2 * lambda[[1]][1] * gammasq[1] * 0.01, p )

  numbofits <- burnin + every * numb ## total number of iterations
  holdpsi <- matrix( rep( 0, p * numb ), ncol = numb )
  holdbeta <- matrix( rep( 0, p * numb ), ncol = numb )
  holdgammasq <- matrix( rep( 0, n.g * numb ), ncol = numb )
  holdlambda <- matrix( rep( 0, n.g * numb ), ncol = numb )
  holdLAMacsept <- matrix( rep( 0, n.g * numb ), ncol = numb )
  holdalpha <- numeric( numb )
  holdW <- numeric( numb )
  holdH <- numeric( numb )

  lambdasd <- rep( 0.01, n.g )

  const <- solvecovs %*% betas

  #################  Updating beta and alpha #############

  for (i in 1:numbofits ){
    LAM <- diag( c( 0, 1/psi ) )
    #LAM[LAM>1e+7]=1e+7
    #LAM[LAM<1e-7 & LAM!=0]=1e-7
    covsLAM <- solvecovs + LAM
    #covsLAM[covsLAM <1e-7 & covsLAM!=0]=1e-7
    varstar <-  solve( covsLAM )
    expec <- varstar %*% const
    cholstar <- chol( varstar, pivot = TRUE )
    if ( attr( cholstar, "rank" ) == ncol( cholstar ) ) {
      pivot <- attr( cholstar, "pivot" )
      cholstar <- cholstar[ ,order( pivot ) ]
      randn <- rnorm( p + 1 )
      be <- expec + crossprod( cholstar, randn )
      alpha <- be[1]
      beta <- be[-1]
    }
    ################### Updating Psi #######
    for (j in 1:p) {
      ### main check 1
      if ( beta[j]^2 < 10^( -5 ) ) {
        check <- 0
        ### subcheck 1 for lambda
        if ( LAMBDA[j] < 0.5 ) {
          while ( check == 0 ) {
            psi[j] <- 1 / rgamma( 1, 0.5 - LAMBDA[j], scale = 1 / ( 0.5 * beta[j]^2 ) )
            u <- runif( 1 )
            check <- as.numeric( u < exp(- psi[j] / ( 2 * GAMMASQ[j] ) ) ) }
        } else {
          while ( check == 0) {
            psi[j] <- rgamma( 1, LAMBDA[j] - 0.5, scale = 2 * GAMMASQ[j] )
            u <- runif( 1 )
            check <- as.numeric( u < exp( - 0.5 * beta[j]^2 / psi[j] ) ) } }
      } else { psi[j] <- rgig( 1, LAMBDA[j] - 0.5, beta[j]^2, 1 / GAMMASQ[j] ) }
    } ### end for (j in 1:p )

    #psi[psi>1e+4]=1e+4
    psi[ psi < 1e-10 &  psi != 0 ] <- 1e-10
    ########################################## update \lambda  for each Group #######
    for ( g in 1:n.g ) {
      mupsi <- 2 * lambdastar[g] * gammasq[g]
      newlambdastar[g] <- lambdastar[g] * exp( lambdasd[g] *rnorm( 1 ) )
      newgammasq[g] <- mupsi / ( 2 * newlambdastar[g] )
      newlambda[[g]] <- rep( newlambdastar[g], len_group[g] )
      logaccept[g] <- log( newlambdastar[g] ) - log( lambdastar[g] ) - 142.85 * ( newlambdastar[g] - lambdastar[g] )
      logaccept[g] <- logaccept[g] - len_group[g] * newlambdastar[g] * log( 2 * newgammasq[g] ) - len_group[g] * lgamma( newlambdastar[g] )
      logaccept[g] <- logaccept[g] + len_group[g] * lambdastar[g] * log( 2 * gammasq[g] ) + len_group[g] * lgamma( lambdastar[g] )
      logaccept[g] <- logaccept[g] + newlambdastar[g] * sum( log( psi[group[[g]]] ) ) - sum( psi[group[[g]]] ) / ( 2 * newgammasq[g] )
      logaccept[g] <- logaccept[g] - lambdastar[g] * sum( log( psi[group[[g]]] ) ) +sum( psi[group[[g]]] ) / ( 2 * gammasq[g] )
      accept[g] <- 1
      if ( logaccept[g] < 0 ) { accept[g] <- exp( logaccept[g] ) }
      lambdasd[g] <- lambdasd[g] + ( accept[g] - 0.3) / i
      lambdaaccept[g] <- lambdaaccept[g] + accept[g]
      lambdacount[g] <- lambdacount[g] + 1
      u <- runif( 1 )
      if ( u < accept[g] ) {
        lambda[[g]] <- newlambda[[g]]
        lambdastar[g] <- newlambdastar[g]
        gammasq[g] <- newgammasq[g]
      }
      ################### Updating gamma^2 for each Group #######

      if (g == 1) {

        sha[g] <- 0.5 * sum( psi[group[[g]]] ) + medstar[1] / ( 2 * lambdastar[g] )
        gammasq[g] <- 1 / rgamma( 1, sum( unlist( lambda[[g]] ) ) + 2, scale = 1 / sha[g] )

      } else if (g == 2) {
        A <- 0.5 * sum( psi[group[[g]]] ) + medstar[1] / ( 2 * lambdastar[g] )
        B <- 0.5 * sum( psi[group[[g]]] ) + medstar[2] / ( 2 * lambdastar[g] )
        delta1 <- ( w * ( medstar[1])^2 ) / ( ( w * ( medstar[1])^2 ) + ( ( A / B )^( len_group[g] * lambdastar[g] + 2 ) * ( 1 - w ) * ( medstar[2])^2 )   )
        U <- runif( 1 )
        if ( U < delta1 ) {
          sha[g] <- 0.5 * sum( psi[group[[g]]] ) + medstar[1] / ( 2 * lambdastar[g] )
          gammasq[g] <- 1 / rgamma( 1, sum( unlist( lambda[[g]] ) ) + 2, scale = 1 / sha[g] )
        } else {
          sha[g] <- 0.5 * sum( psi[group[[g]]] ) + medstar[2] / ( 2 * lambdastar[g] )
          gammasq[g] <- 1 / rgamma( 1, sum( unlist( lambda[[g]] ) ) + 2, scale = 1 / sha[g] )
        }
      } else if (g == 3) {
        sha[g] <- 0.5 * sum( psi[group[[g]]] ) + medstar[2] / ( 2 * lambdastar[g] )
        gammasq[g] <- 1 / rgamma( 1, sum( unlist( lambda[[g]] ) ) + 2, scale = 1 / sha[g] )

      } else {
        A <- 0.5 * sum( psi[group[[g]]] ) + medstar[1] / ( 2 * lambdastar[g] )
        B <- 0.5 * sum( psi[group[[g]]] ) + medstar[2] / ( 2 * lambdastar[g] )
        delta2 <- ( h * ( medstar[1])^2 ) / ( ( h * (medstar[1])^2 ) + ( ( A / B )^(len_group[g] * lambdastar[g] + 2) * ( 1 - h ) * ( medstar[2])^2 )   )
        U <- runif( 1 )
        if (U < delta2) {
          sha[g] <- 0.5 * sum( psi[group[[g]]] ) + medstar[1] / ( 2 * lambdastar[g] )
          gammasq[g] <- 1 / rgamma( 1, sum( unlist( lambda[[g]] ) ) + 2, scale = 1 / sha[g] )
        } else {
          sha[g] <- 0.5 * sum( psi[group[[g]]] ) + medstar[2] / ( 2 * lambdastar[g] )
          gammasq[g] <- 1 / rgamma( 1, sum( unlist( lambda[[g]] ) ) + 2, scale = 1 / sha[g] )
        }
      }

    }
    LAMBDA <- unlist( lambda )[order( unlist( group ) )]
    GAMMASQ <- c( rep( gammasq[1], len_group[1] ),rep( gammasq[2], len_group[2] ), rep(gammasq[3], len_group[3] ), rep( gammasq[4], len_group[4] ) )
    GAMMASQ <- GAMMASQ[ order( unlist( group ) ) ]

    ################### Updating w for each Group #######

    U <- runif( 1 )
    sha1 <- medstar[1] / ( 2 * lambdastar[2] )
    sha2 <- medstar[2] / ( 2 * lambdastar[2] )
    A <- dgamma( gammasq[2], 2, scale = 1 / ( sha1 ) )
    B <- dgamma( gammasq[2], 2, scale = 1 / ( sha2 ) )
    delta3 <- A / ( A + B )
    if ( U < delta3 ) {
      w <- rbeta( 1, 3, 2 )
    } else {
      w <- rbeta( 1, 2, 3 )
    }

    ################### Updating h for each Group #######

    U <- runif( 1 )
    sha1 <- medstar[1] / ( 2 * lambdastar[4] )
    sha2 <- medstar[2] / ( 2 * lambdastar[4] )
    A <-  dgamma( gammasq[4], 2, scale = 1 / ( sha1 ) )
    B <-  dgamma( gammasq[4], 2, scale = 1 / ( sha2 ) )
    delta4 <- A / ( A + ( 4 * B ) )
    if ( U < delta4 ) {
      h <- rbeta( 1, 2, 4 )
    } else {
      h <- rbeta( 1, 1, 5 )
    }

    if ( i > burnin & ( i - burnin ) %% every == 0 ) {
      holdlambda[ , ( i - burnin ) / every ] <- lambdastar
      holdgammasq[ , ( i - burnin ) / every ] <- gammasq
      holdLAMacsept[ , ( i - burnin ) / every ] <- lambdaaccept
      holdalpha[ ( i - burnin ) / every ] <- alpha
      holdbeta[ , ( i - burnin ) / every ] <- beta
      holdpsi[ , ( i - burnin ) / every ] <- psi
      holdW[ ( i - burnin ) / every ] <- w
      holdH[ ( i - burnin ) / every ] <- h

    }


  } ## iterations are over

  list( alpha = holdalpha, beta = holdbeta, psi = holdpsi, lambda = holdlambda,
        gammasq = holdgammasq, H = holdH, W = holdW)
}

