  dyn.load("WHD_R.so")

readParameters <- function(parname) {
    l <- strsplit(gsub(" ","",readLines(parname)),":")
    l <- l[sapply(l,length)==2]
    p <- sapply(l,function(z) as.numeric(z[2]))
   names(p) <- sapply(l,function(z) z[1])
    p <- as.list(p)
    return(p)
}

writeParameters<- function(parfile,par) {
  outfile <- parfile
  l <- length(p)
  cat("Parameter values :\n", file = outfile, append = FALSE)
  for(i in 1:l){
    cat(paste0(names(par[i])," : ",format(p[[i]],scientific=FALSE),"\n"), file = outfile, append = TRUE)
  }
}

orderParameters <- function(par) {    
    par <- with(par,list(alpha,beta,gamma,delta,eta,ksi,phi,psi,omega,u,v,w,l,zcrit))
    names(par) <- c("alpha","beta","gamma","delta","eta","ksi","phi","psi","omega","u","v","w","l","zcrit")
    return(par)
}

F <- function(x,par) {
    n <- length(x)
    par <- orderParameters(par)
    f <- as.double(rep(-1,n))    
    res <- .C("_F",
              as.integer(n),
              as.double(x),
              as.double(f),
              as.double(par))
    
    return(res[[3]])
}

G <- function(y,z,par) {
    n <- length(y)
    if(length(z)!=n) stop("Length of y and z must be the same!")
    par <- orderParameters(par)
    g <- as.double(rep(-1,n))    
    res <- .C("_G",
              as.integer(n),
              as.double(y),
              as.double(z),
              as.double(g),
              as.double(par))
    
    return(res[[4]])
}

dFdx <- function(x,par) {
    n <- length(x)
    par <- orderParameters(par)
    f <- as.double(rep(-1,n))    
    res <- .C("_dFdx",
              as.integer(n),
              as.double(x),
              as.double(f),
              as.double(par))
    return(res[[3]])
}

dGdy <- function(y,z,par) {
    n <- length(y)
    par <- orderParameters(par)
    g <- as.double(rep(-1,n))    
    res <- .C("_dGdy",
              as.integer(n),
              as.double(y),
              as.double(z),
              as.double(g),
              as.double(par))
    
    return(res[[4]])
}

dGdz <- function(y,z,par) {
    n <- length(y)
    par <- orderParameters(par)
    g <- as.double(rep(-1,n))    
    res <- .C("_dGdz",
              as.integer(n),
              as.double(y),
              as.double(z),
              as.double(g),
              as.double(par))
    
    return(res[[4]])
}

H <- function(x,par) {
    n <- length(x)
    par <- orderParameters(par)
    h <- as.double(rep(-1,n))    
    res <- .C("_H",
              as.integer(n),
              as.double(x),
              as.double(h),
              as.double(par))
    
    return(res[[3]])
}

dxyt_dz <- function(x,y,z,par) {
    par <- orderParameters(par)
    var <- c(x,y,0)
    f <- rep(-1,3)
    res <- .C("_dxyt_dz",
              as.double(z),
              as.double(var),
              as.double(f),
              as.double(par))
    return(res[[3]])
}

dxyz_dt <- function(x,y,z,p) {
    p <- orderParameters(p)
    v <- c(x,y,z,0)
    t <- 0
    f <- rep(-1,4)
    res <- .C("_dxyz_dt",
              as.double(t),
              as.double(v),
              as.double(f),
              as.double(p))
    return(res[[3]])
}


dH <- function(x,par) {
    n <- length(x)
    par <- orderParameters(par)
    dh <- as.double(rep(-1,n))    
    res <- .C("_dH",
              as.integer(n),
              as.double(x),
              as.double(dh),
              as.double(par))
    
    return(res[[3]])
}

Pg <- function(x,par) {
    par <- orderParameters(par)
    y <- rep(-1,length(x))
    for(i in 1:length(x)) {
        res <- .C("_P_gamma_eq",
                  as.double(x[i]),
                  as.double(y[i]),
                  as.double(unlist(par)))
        y[i] <- res[[2]]
    }
    return(y)
}

dPg <- function(x,par) {
    par <- orderParameters(par)
    y <- rep(-1,length(x))
    for(i in 1:length(x)) {
        res <- .C("_dP_gamma_eq",
                  as.double(x[i]),
                  as.double(y[i]),
                  as.double(unlist(par)))
        y[i] <- res[[2]]
    }
    return(y)
}

Pa <- function(x,par) {
    par <- orderParameters(par)
    y <- rep(-1,length(x))
    for(i in 1:length(x)) {
        res <- .C("_P_alpha_eq",
                  as.double(x[i]),
                  as.double(y[i]),
                  as.double(unlist(par)))
        y[i] <- res[[2]]
    }
    return(y)
}

dPa <- function(x,par) {
    par <- orderParameters(par)
    y <- rep(-1,length(x))
    for(i in 1:length(x)) {
        res <- .C("_dP_alpha_eq",
                  as.double(x[i]),
                  as.double(y[i]),
                  as.double(unlist(par)))
        y[i] <- res[[2]]
    }
    return(y)
}

gammaClearance <- function(par) {
    par <- orderParameters(par)
    gc <- -1
    res <- .C("_gamma_clearance",
              as.double(gc),
              as.double(unlist(par)))
    return(res[[1]])    
}

gammaSurvival <- function(par) {
    par <- orderParameters(par)    
    gc <- -1
    res <- .C("_gamma_survival",
            as.double(gc),
            as.double(unlist(par)))
  return(res[[1]])        
}

gammaKill <- function(par) {
    par <- orderParameters(par)    
    gk <- -1
    res <- .C("_gamma_kill",
              as.double(gk),
              as.double(unlist(par)))
    return(res[[1]])        
}

alphaKill <- function(par) {
    par <- orderParameters(par)    
    ak <- -1
    res <- .C("_alpha_kill",
              as.double(ak),
              as.double(unlist(par)))
    return(res[[1]])        
}

alphaBlud <- function(target_blud,par,tmax=1000,eps=1e-6,maxiter=50000) {
    par <- orderParameters(par)
    estim <- rep(-1,3)
    test <- -1
    res <- .C("_alpha_blud",
              as.double(estim),
              as.integer(test),
              as.double(target_blud),
              as.double(tmax),
              as.double(eps),
              as.integer(maxiter),
              as.double(par))
    est <- res[[1]]
    attr(est,"code") <- res[[2]]
    return(est)
}

gammaCrit <- function(par,maxg=10) {
    par <- orderParameters(par)
    gcrit <- rep(-1,maxg)
    curv  <- rep(-1,maxg)  
    res <- .C("_gamma_crit",
              as.integer(maxg),
              as.double(gcrit),
              as.integer(curv),
              as.double(unlist(par)))
    v <- unlist(res[[2]][1:res[[1]]])
    return(v[v>0])
}

alphaCrit <- function(par,maxa=10) {
    par <- orderParameters(par)
    acrit <- rep(-1,maxa)
    curv  <- rep(-1,maxa)  
    res <- .C("_alpha_crit",
              as.integer(maxa),
              as.double(acrit),
              as.integer(curv),
              as.double(unlist(par)))
    v <- unlist(res[[2]][1:res[[1]]])
    return(v[v>0])
}

gammaInfection <- function(par,x0=1,dlogy0=0,dlogz0=0,tmax=1000,maxdist=0.001,eps=1e-6,maxiter=100) {
    par <- orderParameters(par)
    estim <- rep(-1,3)
    test <- -1
    res <- .C("_gamma_infection",
              as.double(estim),
              as.integer(test),
              as.double(x0),
              as.double(dlogy0),
              as.double(dlogz0),
              as.double(tmax),
              as.double(maxdist),
              as.double(eps),
              as.integer(maxiter),
              as.double(par));
    estim <- res[[1]]
    attr(estim,"code") <- res[[2]]
    return(estim)
}

alphaInfection <- function(par,x0=1,dlogy0=0,dlogz0=0,tmax=1000,maxdist=0.001,eps=1e-6,maxiter=100) {
    par <- orderParameters(par)
    estim <- rep(-1,3)
    test <- -1
    res <- .C("_alpha_infection",
              as.double(estim),
              as.integer(test),
              as.double(x0),
              as.double(dlogy0),
              as.double(dlogz0),
              as.double(tmax),
              as.double(maxdist),
              as.double(eps),
              as.integer(maxiter),
              as.double(par));
    estim <- res[[1]]
    attr(estim,"code") <- res[[2]]
    return(estim)
}

alphaSPBL <- function(par,alpha_lo,alpha_hi,x0=1,dlogy0=0,dz0=0,tmax=1000,maxdist=0.001,eps=1e-6,maxiter=100) {
    par <- orderParameters(par)
    estim <- rep(-1,3)
    test <- -1
    res <- .C("_alpha_spbl",
              as.double(estim),
              as.integer(test),
              as.double(x0),
              as.double(dlogy0),
              as.double(dz0),
              as.double(alpha_lo),
              as.double(alpha_hi),
              as.double(tmax),
              as.double(eps),
              as.integer(maxiter),
              as.double(par));
    estim <- res[[1]]
    attr(estim,"code") <- res[[2]]
    return(estim)
}

yhomeo <- function(par) {
    par <- orderParameters(par)
    res <- .C("_yhomeo",
              as.double(-1),
              as.integer(-1),
              as.double(unlist(par)))
    if(res[[2]]<0) return(NA)
    return(res[[1]]);
}

equilibrium <- function(par) {
    par <- orderParameters(par)
    x <- y <- z <- re_ev <- im_ev <- im_ev_glob <- as.double(rep(-1,10))
    res <- .C("_equilibrium",
              as.integer(length(x)),
              x,y,z,re_ev,im_ev,im_ev_glob,as.double(par))
    nbeq <- res[[1]]
    if(nbeq<=0) return(NULL)
    eq <- data.frame(x=res[[2]][1:nbeq],y=res[[3]][1:nbeq],z=res[[4]][1:nbeq],
                     re_ev=res[[5]][1:nbeq],im_ev=res[[6]][1:nbeq],im_ev_glob=res[[7]][1:nbeq])
    eq$z[is.nan(eq$z)] <- NA
    return(eq)
}

jacobian <- function(eq,par,logx=FALSE) {
    par <- orderParameters(par)
    dfdy <- rep(-1,9)
    if(logx) {
        res <- .C("_jacobianLogx",
                  as.double(eq),
                  as.double(dfdy),
                  as.double(par))
    } else {
        res <- .C("_jacobian",
                  as.double(eq),
                  as.double(dfdy),
                  as.double(par))
    }
    j <- matrix(res[[2]],ncol=3,byrow=T)
    return(j)
}

dyninf <- function(xinit,yinit,zinit,t,par) {
    par <- orderParameters(par)
    n <- length(t)
    x <- y <- z <- rep(-1,n)
    x[1] <- xinit
    y[1] <- yinit
    z[1] <- zinit
    res <- .C("_dyninf",
              as.integer(n),
              as.double(x),as.double(y),as.double(z),
              as.double(t),
              as.double(unlist(par)))
    l <- 1:res[[1]]
    return(data.frame(t=res[[5]],x=res[[2]][l],y=res[[3]][l],z=res[[4]][l]))
}

tsurv <- function(x,y,z,par,tmax=1000,eps=1e-6,maxiter=1000) {
    par <- orderParameters(par)
    t <- 0;
    res <- .C("_tsurv",
              as.double(x),as.double(y),as.double(z),as.double(t),
              as.double(tmax),as.double(eps),as.integer(maxiter),
              as.double(unlist(par)))
    return(c(res[[1]],res[[2]],res[[3]],res[[4]]));
}

tsurvV2 <- function(x,y,z,par,tmax=1000,eps=1e-6,maxiter=1000) {
    par <- orderParameters(par)
    t <- 0;
    res <- .C("_tsurvV2",
              as.double(x),as.double(y),as.double(z),as.double(t),
              as.double(tmax),as.double(eps),as.integer(maxiter),
              as.double(unlist(par)))
    return(c(res[[1]],res[[2]],res[[3]],res[[4]]));
}
    
survsim <- function(n,par,par_ini,tmax=1000,eps=1e-6,maxiter=1000,logscale=FALSE) {
    par <- orderParameters(par)
    time <- x  <- y <- z <- x0 <- y0 <- z0 <- as.double(rep(-1000,n))
    dead <- as.integer(rep(-1,n))

    
    if(logscale) {
        par_ini <- with(par_ini,c(mu_logx,mu_logy,mu_logz,sigma_logx,sigma_logy,sigma_logz))
    } else {
        par_ini <- with(par_ini,c(mu_x,mu_y,mu_z,sigma_logx,sigma_logy,sigma_logz))
    }    
    logscale <- ifelse(logscale,1,0)  

    
    res <- .C("_survsim",
              as.integer(n),
              as.double(unlist(par)),
              as.double(unlist(par_ini)),
              as.integer(logscale),
              as.double(tmax),
              as.double(eps),
              as.integer(maxiter),
              as.double(x0),as.double(y0),as.double(z0),as.double(x),as.double(y),as.double(z),as.double(time),
              as.integer(dead))
    
    d <- data.frame(x0=res[[8]], y0=res[[9]], z0=res[[10]], x=res[[11]], y=res[[12]], z=res[[13]],time=res[[14]],dead=res[[15]])
    attr(d,"Parameters") = par
    attr(d,"StartingConditions") = par_ini
    return(d)
}

threshold_x <- function(y,z,par,tmax=2000,eps=0.001,max_dist=0.001) {
    par <- orderParameters(par)
    nb <- length(y)
    if(length(z)!=nb) {
        warning("y and z have different length!")
    }
    x <- lo_x <- hi_x <- rep(-1,nb)
    
    res <- .C("_threshold_x",
              as.integer(nb),
              as.double(x),
              as.double(lo_x),
              as.double(hi_x),
              as.double(y),
              as.double(z),
              as.double(tmax),
              as.double(eps),
              as.double(max_dist),
              as.double(par))
    
    return(data.frame(x=res[[2]],lo_x=res[[3]],hi_x=res[[4]],y=res[[5]],z=res[[6]]))
}

threshold_y <- function(x,z,par,tmax=2000,eps=0.001,max_dist=0.001) {
    par <- orderParameters(par)
    nb <- length(x)
    if(length(z)!=nb) {
        warning("x and z have different length!")
    }
    y <- lo_y <- hi_y <- rep(-1,nb)
    
    res <- .C("_threshold_y",
              as.integer(nb),
              as.double(x),
              as.double(y),
              as.double(lo_y),
              as.double(hi_y),
              as.double(z),
              as.double(tmax),
              as.double(eps),
              as.double(max_dist),
              as.double(par))
    
    return(data.frame(x=res[[2]],y=res[[3]],lo_y=res[[4]],hi_y=res[[5]],z=res[[6]]))
}

threshold_z <- function(x,y,par,tmax=2000,eps=0.001,max_dist=0.001) {
    par <- orderParameters(par)
    nb <- length(y)
    if(length(x)!=nb) {
        warning("x and y have different length!")
    }
    z <- lo_z <- hi_z <- rep(-1,nb)
    
    res <- .C("_threshold_z",
              as.integer(nb),
              as.double(x),
              as.double(y),
              as.double(z),
              as.double(lo_z),
              as.double(hi_z),
              as.double(tmax),
              as.double(eps),
              as.double(max_dist),
              as.double(par))
    
    return(data.frame(x=res[[2]],y=res[[3]],z=res[[4]],lo_z=res[[5]],hi_z=res[[6]]))
}

DLEx <- function(x,y,z,par,tmax=1000,eps=1e-6) {
    par <- orderParameters(par)
    n <- length(x)
    if(length(y)!=n || length(z)!=n) {
        y <- rep(y[1],n)
        z <- rep(z[1],n)
    }
    dle <- rep(-1,n)    
    res <- .C("_DLEx",
              as.integer(n),as.double(dle),as.double(x),as.double(y),as.double(z),
              as.double(tmax),as.double(eps),as.double(p))
    return(res[[2]])
}

DLE <- function(x,y,z,par,tmax=1000,eps=1e-6) {
    par <- orderParameters(par)
    n <- length(x)
    if(length(y)!=n || length(z)!=n) {
        y <- rep(y[1],n)
        z <- rep(z[1],n)
    }
    dle <- rep(-1,n)    
    res <- .C("_DLE",
              as.integer(n),as.double(dle),as.double(x),as.double(y),as.double(z),
              as.double(tmax),as.double(eps),as.double(p))
    return(res[[2]])
}

DLEmatrix <- function(x,y,z,par,tmax=1000,eps=1e-6) {
    par <- orderParameters(par)
    L <- rep(0,9)
    res <- .C("_DLEmatrix",
              as.double(L),
              as.double(x),as.double(y),as.double(z),
              as.double(tmax),
              as.double(eps),
              as.double(par))
    L <- matrix(res[[1]],nrow=3,ncol=3,byrow=T)
    return(L)
}

spbl <- function(x0,y0,z0,par,tmax=1000,maxiter=1000,eps=1e-6) {
    par <- orderParameters(par)
    v <- c(x0,y0,z0,0,0,0,0)
    test <- 0
    res <- .C("_spbl",
              as.double(v),
              as.integer(test),
              as.double(tmax),              
              as.double(eps),
              as.integer(maxiter),
              as.double(unlist(par)))
    v <- res[[1]]
    names(v) <- c("xSPBL","ySPBL","zSPBL","tSPBL","xPeak","start","duration")
    attr(v,"code") <- res[[2]]
    return(v)
}

threshold_spbl_x <- function(par,dy=0,dz=0,tmax=1000,eps=1e-6,maxiter=500) {
    par <- orderParameters(par)
    lx <- 0
    res <- .C("_threshold_spbl_x",
              as.double(lx),
              as.double(dy),
              as.double(dz),
              as.double(tmax),
              as.double(eps),
              as.integer(maxiter),
              as.double(par))
   return(res[[1]])           
}

blud <- function(par,tmax=1000,maxiter=1000,eps=1e-6) {
    par <- orderParameters(par)
    v <- rep(-1,4)
    test <- 0
    res <- .C("_blud",
              as.double(v),
              as.integer(test),
              as.double(tmax),              
              as.double(eps),
              as.integer(maxiter),
              as.double(unlist(par)))
    v <- res[[1]]
    attr(v,"code") <- res[[2]]
    return(v)
    if(is.nan(v[3]) || abs(v[3]-p$zcrit)>1e-3) v <- rep(NA,length(v))
    return(v)
}

control <- function(p,x0=0.1,tmax=1000,eps=1e-6,maxiter=1000) {
    p <- orderParameters(p)
    v <- rep(-1,4)
    res <- .C("_dyninf_toControl",
              as.double(v),
              as.double(x0),
              as.double(tmax),
              as.double(eps),
              as.integer(maxiter),              
              as.double(unlist(p)))
    v <- res[[1]]
    dv <- dxyz_dt(v[1],v[2],v[3],p)
    print(dv)
    dv <- sqrt(sum(dv[-4]^2))    
    if(dv < 10*eps) v <- rep(NA,4)
    return(list(v,dv,res[[4]]))
}

maxResponse <- function(par,x0=0.1,tmax=1000,eps=1e-6,maxiter=1000) {
    par <- orderParameters(par)    
    v <- rep(-1,4)
    res <- .C("_dyninf_toMaxResponse",
              as.double(v),
              as.double(x0),
              as.double(tmax),
              as.double(eps),
              as.integer(maxiter),
              as.double(unlist(par)))
    return(c(res[[1]],res[[4]]))
}
