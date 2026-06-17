
YgVal <- function(cdat, ndat, g) {
  N <- length(cdat) ; ndatcsum <- cumsum(ndat) 
  delhat <- 0 ; tbar <- 0
  for (k in 1:g) {
    sample <- circular(0)
    if (k==1) {low <- 0} else
      if (k > 1) {low <- ndatcsum[k-1]}
    for (j in 1:ndat[k]) { sample[j] <- cdat[j+low] }
    tm1 <- trigonometric.moment(sample, p=1)
    tm2 <- trigonometric.moment(sample, p=2)
    Rbar1 <- tm1$rho; Rbar2 <- tm2$rho ; tbar[k] <- tm1$mu
    delhat[k] <- (1-Rbar2)/(2*Rbar1*Rbar1)
  }
  dhatmax <- max(delhat) ; dhatmin <- min(delhat)
  if (dhatmax/dhatmin <= 4) {
    CP <- 0 ; SP <- 0 ; dhat0 <- 0
    for (k in 1:g) {
      CP <- CP + ndat[k]*cos(tbar[k])
      SP <- SP + ndat[k]*sin(tbar[k])
      dhat0 <- dhat0 + ndat[k]*delhat[k] 
    }
    dhat0 <- dhat0/N
    RP <- sqrt(CP*CP+SP*SP)
    Yg <- 2*(N-RP)/dhat0
    return(Yg) } 
  else if (dhatmax/dhatmin > 4) {
    CM <- 0 ; SM <- 0 ; Yg <- 0
    for (k in 1:g) {
      CM <- CM + (ndat[k]*cos(tbar[k])/delhat[k])
      SM <- SM + (ndat[k]*sin(tbar[k])/delhat[k])
      Yg <- Yg + (ndat[k]/delhat[k]) 
    }
    RM <- sqrt(CM*CM+SM*SM)
    Yg <- 2*(Yg-RM)
    return(Yg) }
}


################################################################################
# Bootstrap version of Watson's test for a common mean direction
################################################################################

YgTestBoot <- function(cdat, ndat, g, indsym, B) {
  N <- length(cdat) ; ndatcsum <- cumsum(ndat) 
  delhat <- 0 ; tbar <- 0 ; centdat <- circular(0)
  for (k in 1:g) {
    sample <- circular(0)  
    if (k==1) {low <- 0} else
      if (k > 1) {low <- ndatcsum[k-1]}
    for (j in 1:ndat[k]) { sample[j] <- cdat[j+low] }
    tm1 <- trigonometric.moment(sample, p=1)
    tm2 <- trigonometric.moment(sample, p=2)
    Rbar1 <- tm1$rho; Rbar2 <- tm2$rho ; tbar[k] <- tm1$mu
    delhat[k] <- (1-Rbar2)/(2*Rbar1*Rbar1)
    if (tbar[k] < 0) {tbar[k] <- tbar[k]+2*pi}
    centsamp <- sample-tbar[k]
    if (indsym == 1) {centsamp <- c(centsamp, -centsamp)}
    centdat <- c(centdat,centsamp)
  }
  centdat <- centdat[-1]
  dhatmax <- max(delhat) ; dhatmin <- min(delhat)
  if (dhatmax/dhatmin <= 4) {
    PorM <- 1 ; CP <- 0 ; SP <- 0 ; dhat0 <- 0
    for (k in 1:g) {
      CP <- CP + ndat[k]*cos(tbar[k])
      SP <- SP + ndat[k]*sin(tbar[k])
      dhat0 <- dhat0 + ndat[k]*delhat[k] 
    }
    dhat0 <- dhat0/N
    RP <- sqrt(CP*CP+SP*SP)
    Yg <- 2*(N-RP)/dhat0
  } else
    if (dhatmax/dhatmin > 4) {
      PorM <- 0 ; CM <- 0 ; SM <- 0 ; Yg <- 0
      for (k in 1:g) {
        CM <- CM + (ndat[k]*cos(tbar[k])/delhat[k])
        SM <- SM + (ndat[k]*sin(tbar[k])/delhat[k])
        Yg <- Yg + (ndat[k]/delhat[k]) 
      }
      RM <- sqrt(CM*CM+SM*SM)
      Yg <- 2*(Yg-RM)
    }
  YgObs <- Yg ; nxtrm <- 1
  
  if (indsym == 0) {
    for (b in 1:B) {
      centsamp <- circular(0) 
      for (k in 1:g) {
        if (k==1) {low <- 0} else
          if (k > 1) {low <- ndatcsum[k-1]}
        for (j in 1:ndat[k]) { centsamp[j] <- centdat[j+low] }
        bootsamp <- sample(centsamp, size=ndat[k], replace=TRUE)
        tm1 <- trigonometric.moment(bootsamp, p=1)
        tm2 <- trigonometric.moment(bootsamp, p=2)
        Rbar1 <- tm1$rho; Rbar2 <- tm2$rho ; tbar[k] <- tm1$mu
        delhat[k] <- (1-Rbar2)/(2*Rbar1*Rbar1)
      }
      if (PorM == 1) {
        CP <- 0 ; SP <- 0 ; dhat0 <- 0
        for (k in 1:g) {
          CP <- CP + ndat[k]*cos(tbar[k])
          SP <- SP + ndat[k]*sin(tbar[k])
          dhat0 <- dhat0 + ndat[k]*delhat[k] 
        }
        dhat0 <- dhat0/N
        RP <- sqrt(CP*CP+SP*SP)
        Yg <- 2*(N-RP)/dhat0
      } else
        if (PorM == 0) {
          CM <- 0 ; SM <- 0 ; Yg <- 0
          for (k in 1:g) {
            CM <- CM + (ndat[k]*cos(tbar[k])/delhat[k])
            SM <- SM + (ndat[k]*sin(tbar[k])/delhat[k])
            Yg <- Yg + (ndat[k]/delhat[k]) 
          }
          RM <- sqrt(CM*CM+SM*SM)
          Yg <- 2*(Yg-RM)
        }
      YgBoot <- Yg
      if (YgBoot >= YgObs) {nxtrm <- nxtrm+1}
    }
    pval <- nxtrm/(B+1)
    return(c(YgObs, pval))
  } else
    
    if (indsym == 1) {
      for (b in 1:B) {
        centsamp <- circular(0) 
        for (k in 1:g) {
          if (k==1) {low <- 0} else
            if (k > 1) {low <- 2*ndatcsum[k-1]}
          for (j in 1:(2*ndat[k])) { centsamp[j] <- centdat[j+low] }
          bootsamp <- sample(centsamp, size=ndat[k], replace=TRUE)
          tm1 <- trigonometric.moment(bootsamp, p=1)
          tm2 <- trigonometric.moment(bootsamp, p=2)
          Rbar1 <- tm1$rho; Rbar2 <- tm2$rho ; tbar[k] <- tm1$mu
          delhat[k] <- (1-Rbar2)/(2*Rbar1*Rbar1)
        }
        if (PorM == 1) {
          CP <- 0 ; SP <- 0 ; dhat0 <- 0
          for (k in 1:g) {
            CP <- CP + ndat[k]*cos(tbar[k])
            SP <- SP + ndat[k]*sin(tbar[k])
            dhat0 <- dhat0 + ndat[k]*delhat[k] 
          }
          dhat0 <- dhat0/N
          RP <- sqrt(CP*CP+SP*SP)
          Yg <- 2*(N-RP)/dhat0
        } else
          if (PorM == 0) {
            CM <- 0 ; SM <- 0 ; Yg <- 0
            for (k in 1:g) {
              CM <- CM + (ndat[k]*cos(tbar[k])/delhat[k])
              SM <- SM + (ndat[k]*sin(tbar[k])/delhat[k])
              Yg <- Yg + (ndat[k]/delhat[k]) 
            }
            RM <- sqrt(CM*CM+SM*SM)
            Yg <- 2*(Yg-RM)
          }
        YgBoot <- Yg
        if (YgBoot >= YgObs) {nxtrm <- nxtrm+1}
      }
      pval <- nxtrm/(B+1)
      return(c(YgObs, pval))
    }
}
