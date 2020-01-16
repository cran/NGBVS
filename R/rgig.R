rgig <- function(n = 10, lambda = 1, chi = 1, psi = 1) {
    if ( (chi < 0) | (psi < 0) )  stop("Invalid parameters for GIG")
    if ( (chi == 0) & (lambda <= 0) )  stop("Invalid parameters for GIG")
    if ( (psi == 0) & (lambda >= 0) )  stop("Invalid parameters for GIG")

    if ( (chi == 0) & (lambda > 0) ) {
        ## Gamma distribution -> Variance Gamma
        return(rgamma(n, shape = lambda, rate = psi / 2))
    }
    if ( (psi == 0) & (lambda < 0) ) {
        ## Inverse gamma distribution -> Student-t
        return(1 / rgamma(n, shape = - lambda, rate = chi / 2))
    }
    ## IG distribution:
    ## An implementation of the algorithm described in Raible (2000),
    ## copied from the package fBasics.
    if ( lambda == -0.5 )  {
        U <- runif(n)
        V <- rnorm(n)^2
        delta <- sqrt(chi)

        z1 <- function(v, delta, gamma) {
            delta/gamma + v/(2 * gamma^2) - sqrt( v * delta/(gamma^3) + ( v/(2 * gamma^2) )^2 )
        }

        z2 <- function(v, delta, gamma) {
            (delta/gamma)^2/z1(v = v, delta = delta, gamma = gamma)
        }

        pz1 <- function(v, delta, gamma) {
            delta/( delta + gamma * z1(v = v, delta = delta, gamma = gamma) )
        }

        s <- ( 1 - sign(U - pz1(v = V, delta = delta, gamma = sqrt(psi))) ) / 2

        return(z1(v = V, delta = delta, gamma = sqrt(psi)) * s +
               z2(v = V, delta = delta, gamma = sqrt(psi)) * (1 - s))

    } else if (lambda == 1) {
        ## hyp distribution:
        ## An implementation of the algorithm described in Dagpunar (1989),
        ## copied from Rmetrics package 'fBasics'.
        alpha <- sqrt(psi / chi)
        beta <- sqrt(psi * chi)
        m <- sign(beta)

        g <- function(y) {
            return(0.5 * beta * y^3 - y^2 * (0.5 * beta + 2) + y * (-0.5 * beta) + 0.5 * beta)
        }

        upper <- 1
        while ( g(upper) <= 0 )  upper <- 2 * upper

        yM <- uniroot( g, interval = c(0, 1) )$root
        yP <- uniroot( g, interval = c(1, upper) )$root
        a <- (yP - 1) * exp( -0.25 * beta * (yP + 1/yP - 2) )
        b <- (yM - 1) * exp( -0.25 * beta * (yM + 1/yM - 2) )
        ca <- -0.25 * beta * 2
        output <- numeric(n)

        for (i in 1:n) {
            need.value <- TRUE
            while (need.value == TRUE) {
                R1 <- runif(1)
                R2 <- runif(1)
                Y <- 1 + a * R2/R1 + b * (1 - R2)/R1
                if (Y > 0) {
                    if (-log(R1) >= 0.25 * beta * (Y + 1/Y) + ca) {
                       need.value <- FALSE
                    }
                }
            }
            output[i] <- Y
        }
        return(output/alpha)

    } else {
        ## GIG distribution:
        ## An implementation of the algorithm described in Dagpunar (1989),
        ## copied from Rmetrics package 'fBasics'.
        alpha <- sqrt(psi / chi)
        beta <- sqrt(psi * chi)
        m <- ( lambda - 1 + sqrt((lambda - 1)^2 + beta^2) ) / beta

        g <- function(y) {
            0.5 * beta * y^3 - y^2 * (0.5 * beta * m + lambda + 1) +
            y * ((lambda - 1) * m - 0.5 * beta) + 0.5 * beta * m
        }
        upper <- m
        while ( g(upper) <= 0 ) upper <- 2 * upper

        yM <- uniroot( g, interval = c(0, m) )$root
        yP <- uniroot( g, interval = c(m, upper) )$root
        a <- (yP - m) * (yP/m)^( 0.5 * (lambda - 1) ) * exp(-0.25 * beta * (yP + 1/yP - m - 1/m))
        b <- (yM - m) * (yM/m)^(0.5 * (lambda - 1)) * exp(-0.25 * beta * (yM + 1/yM - m - 1/m))
        ca <- -0.25 * beta * (m + 1/m) + 0.5 * (lambda - 1) * log(m)
        output <- numeric(n)

        for (i in 1:n) {
            need.value <- TRUE
            while ( need.value == TRUE ) {
                R1 <- runif(1)
                R2 <- runif(1)
                Y <- m + a * R2/R1 + b * (1 - R2)/R1
                if (Y > 0) {
                    if ( -log(R1) >= -0.5 * (lambda - 1) * log(Y) + 0.25 * beta * (Y + 1/Y) + ca ) {
                        need.value <- FALSE
                    }
                }
            }
            output[i] <- Y
        }
        return(output/alpha)
    }
}
