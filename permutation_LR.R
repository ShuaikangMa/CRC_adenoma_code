rm(list=ls())
library(deSolve)
library(parallel)
library(pbapply)
library(deSolve)
library(reshape2)
library(ggplot2)
library(orthopolynom)
library(glmnet)
library(melt)
library(scales)
library(reshape2)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggrepel)
library(mvtnorm)

power_equation1 <- function(x, power_par){power_par[1]*x^power_par[2]}
power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),
                                                   function(c) power_par[c,1]*x^power_par[c,2] )) }
get_legendre_matrix <- function(x,legendre_order){
  legendre_coef <- legendre.polynomials(n = legendre_order, normalized=F)
  legendre_matrix <- as.matrix(as.data.frame(polynomial.values(
    polynomials = legendre_coef, x = scaleX(x, u = -1, v = 1))))
  colnames(legendre_matrix) <- paste0("legendre_",0:legendre_order)
  return(legendre_matrix[,2:(legendre_order+1)])
}

get_legendre_par <- function(y,legendre_order,x){
  legendre_par <- as.numeric(coef(lm(y~get_legendre_matrix(x,legendre_order))))
  return(legendre_par)
}

legendre_fit <- function(par,x){
  legendre_order = length(par)
  fit <- sapply(1:length(par), function(c)
    par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
  legendre_fit <- as.matrix(as.data.frame(polynomial.values(
    polynomials = fit, x = scaleX(x, u = -1, v = 1))))
  x_interpolation <- rowSums(legendre_fit)
  return(x_interpolation)
}

darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

qdODEmod <- function(Time, State, Pars, power_par) {
  nn <- length(Pars)
  n_state <- length(State)
  
  ind_effect <- paste0("alpha", "*", names(State)[1])
  
  if(nn > 1){
    dep_effect <- sapply(2:nn, function(c) paste0("beta", c-1, "*", names(State)[c]))
  } else {
    dep_effect <- character(0)
  }
  dep_effect <- paste0(dep_effect, collapse = "+")
  
  all_effect <- if(nchar(dep_effect) > 0) paste(ind_effect, dep_effect, sep = "+") else ind_effect
  
  expr <- parse(text = all_effect)
  
  with(as.list(c(State, Pars)), {
    dx <- eval(expr)
    if(!is.null(power_par) && n_state > 1){
      dy <- power_par[,1] * power_par[,2] * Time^(power_par[,2]-1)
      if(length(dy) != (n_state-1)) dy <- rep(dy, length.out = n_state-1)
    } else {
      dy <- numeric(0)
    }
    dind <- alpha * x
    ddep_list <- list()
    if(nn > 1){
      for(i in 1:(nn-1)){
        y_name <- if(i <= (n_state-1)) paste0("y", i) else NULL
        beta_name <- paste0("beta", i)
        if(!is.null(y_name) && exists(y_name)) {
          tmp <- paste0(beta_name, "*", y_name)
          expr2 <- parse(text = tmp)
          ddep_list[[i]] <- eval(expr2)
        } else {
          ddep_list[[i]] <- 0
        }
      }
    }
    ret <- c(dx, dy, dind, unlist(ddep_list))
    if(length(ret) < n_state) ret <- c(ret, rep(0, n_state - length(ret)))
    if(length(ret) > n_state) ret <- ret[1:n_state]
    
    return(list(ret))
  })
}

qdODEmod_lgkt <- function(Time, State, Pars, power_par) {
  nn = length(Pars)
  ind_effect = paste0("alpha", "*", names(State)[1])
  dep_effect = sapply(2:nn, function(c) paste0(paste0("beta", c-1), "*", names(State)[c]))
  dep_effect = paste0(dep_effect, collapse = "+")
  all_effect = paste0(ind_effect, "+", dep_effect)
  expr = parse(text = all_effect)
  
  with(as.list(c(State, Pars)), {
    dx = eval(expr)
    dy <- power_par[, 1] * power_par[, 2] * Time^(power_par[, 2] - 1)
    dind = alpha * x
    ddep = sapply(1:(nn-1), function(i) {
      tmp = paste0(paste0("beta", i), "*", paste0("y", i))
      expr2 = parse(text = tmp)
      eval(expr2)
    })
    return(c(dx, dy, dind, ddep))
  })
}

qdODE_ls <- function(pars, data, Time, power_par,bc,lx){
  n = length(pars)
  power_par = as.matrix(power_par)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }

  if(lx ==  "ode"){
    out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                            times = Time, power_par = power_par))
    X = as.numeric(data[1,])
    fit = as.numeric(out[,2])
    ind = as.numeric(out[,(n+2)])
    sse = sum(crossprod(X-fit),sum((ind[ind<0])^2))
  }
  if(lx == "lgkt"){
    rk4 <- function(ode_func, State, Time, Pars, power_par, dt) {
      k1 <- ode_func(Time, State, Pars, power_par)
      k2 <- ode_func(Time + 0.5 * dt, State + 0.5 * dt * k1, Pars, power_par)
      k3 <- ode_func(Time + 0.5 * dt, State + 0.5 * dt * k2, Pars, power_par)
      k4 <- ode_func(Time + dt, State + dt * k3, Pars, power_par)
      
      new_state <- State + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      
      return(new_state)
    }
    
    integrate_rk4 <- function(ode_func, State, Time_range, Pars, power_par, dt) {
      #times <- seq(Time_range[1], Time_range[2], by = dt)
      times = Time
      states <- matrix(NA, nrow = length(times), ncol = length(State))
      colnames(states) <- names(State)
      states[1, ] <- State
      
      current_time <- Time_range[1]
      current_state <- State
      
      for (i in 2:length(times)) {
        current_time <- times[i-1]
        current_state <- rk4(ode_func, current_state, current_time, Pars, power_par, dt)
        states[i, ] <- current_state
      }
      
      return(cbind(times,states))
    }
    Time_range <- c(min(Time), max(Time))
    result <- integrate_rk4(qdODEmod_lgkt, State, Time_range, Pars, power_par, dt=((max(Time)-min(Time))/length(Time))/bc)
    out = as.data.frame(result)
    X = as.numeric(data[1,])
    fit = as.numeric(out[,2])
    ind = as.numeric(out[,(n+2)])
    alpha=1e-1
    beta=1e-1
    sse = sum(crossprod(X-fit),sum((ind[ind<0])^2),alpha*(Pars[1]^2),beta*(sum(Pars[2:n])^2))
  }
  
  
  return(sse)
}

qdODE_fit <- function(pars, data, Time, power_par, LOP_order = 6, new_time = NULL, n_expand = 100,bc,ind_par,lx){
  n = length(pars)
  if (n == 1) {
    Pars <- c(alpha = pars[1])
    State <- c(x = data[1,1], y1 = matrix(numeric(0), nrow = 0, ncol = 1), ind = data[1,1], dep = numeric(0))
  } else if (n == 2) {
    Pars <- c(alpha = pars[1], beta1 = pars[2:n])
    power_par <- t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else {
    Pars <- c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }

  if(lx ==  "ode"){
    out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                            times = Time, power_par = power_par))
    out2 <- data.frame(
      x = out[, "time"],
      y = as.numeric(data[1,]),
      y.fit = out[, "x"],
      ind = if("ind" %in% colnames(out)) out[, "ind"] else NA
    )
    
    dep_cols <- setdiff(colnames(out), c("time", "x", "ind"))
    if(length(dep_cols) > 0){
      out2 <- cbind(out2,dep = out[,(n+3):(ncol(out))])
    } else {
      out2 = out2
    }
  }
  if(lx ==  "lgkt"){
    rk4 <- function(ode_func, State, Time, Pars, power_par, dt) {
      k1 <- ode_func(Time, State, Pars, power_par)
      k2 <- ode_func(Time + 0.5 * dt, State + 0.5 * dt * k1, Pars, power_par)
      k3 <- ode_func(Time + 0.5 * dt, State + 0.5 * dt * k2, Pars, power_par)
      k4 <- ode_func(Time + dt, State + dt * k3, Pars, power_par)
      
      new_state <- State + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      
      return(new_state)
    }
    
    integrate_rk4 <- function(ode_func, State, Time_range, Pars, power_par, dt) {
      times <- seq(Time_range[1], Time_range[2], by = dt)
      states <- matrix(NA, nrow = length(times), ncol = length(State))
      colnames(states) <- names(State)
      states[1, ] <- State
      
      current_time <- Time_range[1]
      current_state <- State
      
      for (i in 2:length(times)) {
        current_time <- times[i-1]
        current_state <- rk4(ode_func, current_state, current_time, Pars, power_par, dt)
        states[i, ] <- current_state
      }
      
      return(cbind(times,states))
    }
    Time_range <- c(min(Time), max(Time))
    result <- integrate_rk4(qdODEmod_lgkt, State, Time_range, Pars, power_par, dt=((max(Time)-min(Time))/length(Time))/bc)
    out = as.data.frame(result)
    out2 = data.frame(x = out[,1], y = power_equation1(seq(Time_range[1], Time_range[2], by = ((max(Time)-min(Time))/length(Time)/bc)),ind_par), y.fit = out[,2],
                      ind = out[,(n+2)], dep = out[,(n+3):(ncol(out))])
  }
  
  if (ncol(out2)>4){
    colnames(out2)[4:ncol(out2)] = c(rownames(data)[1], rownames(data)[2:n])
  }else{
    out2 = out2
  }

  rownames(out2) = NULL
  all_LOP_par = sapply(2:ncol(out2),function(c)get_legendre_par(out2[,c], LOP_order, out2$x))
  
  if (is.null(new_time)) {
    time2 = seq(min(Time), max(Time), length = n_expand)
    out3 = apply(all_LOP_par, 2, legendre_fit, x = time2)
    out3 = cbind(time2, out3)
  } else{
    out3 = apply(all_LOP_par, 2, legendre_fit, x = new_time)
    out3 = cbind(new_time, out3)
  }
  colnames(out3) = colnames(out2)
  
  result = list(fit = out2,
                predict = data.frame(out3),
                LOP_par = all_LOP_par)
  return(result)
}

qdODE_all <- function(result, relationship, i, init_pars = 1, LOP_order = LOP_order, methods = "ls",
                      new_time = NULL, n_expand = 100, maxit = 1e3,bc = bc,lx){
  Time = as.numeric(colnames(result$power_fit))
  variable = c(relationship[[i]]$ind.name, relationship[[i]]$dep.name)
  data = result$power_fit[variable,]
  if (length(variable)<=1) {
    qdODE.est = NA
    result = NA
    return.obj <- append(result, list(ODE.value = NA,
                                      parameters = NA))
  } else{
    power_par = result$power_par[variable,][-1,]
    n = nrow(data)
    pars_int = c(init_pars,relationship[[i]]$coefficient)
    if (methods == "ls") {
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,bc=bc,lx = lx,
                         method = "L-BFGS-B",
                         lower = c(rep(-10,(length(pars_int)))),
                         upper = c(rep(10,(length(pars_int)))),
                         control = list(trace = TRUE, maxit = maxit))
      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time,
                          bc = bc,
                          ind_par = result$power_par[variable,][1,],
                          lx = lx)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    } else{
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,bc=bc,lx = lx ,
                         method = "L-BFGS-B",
                         lower = c(rep(-10,(length(pars_int)))),
                         #lower = c(0, rep(-10,(length(pars_int))-1)),
                         upper = c(rep(10,(length(pars_int)))),
                         control = list(trace = TRUE, maxit = maxit))
      
      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time,
                          bc = bc,
                          ind_par = result$power_par[variable,][1,],
                          lx = lx)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    }
  }
  return(return.obj)
}
get_interaction <- function(data, col,alpha,gamma, scaler,reduction = FALSE ){
  if (nrow(data)==2) {
    return_obj = list(ind.name = rownames(data)[col],
                      dep.name = rownames(data)[-col],
                      coefficient = cor(t(data))[1,2])
    
  } else{
    data <- t(data); name <- colnames(data)
    if(scaler == T){
      data = scale(data)
      # data = log(data)
    }
    y = as.matrix(data[,col])
    x = as.matrix(data[,-col])
    
    n <- ncol(x)
    if (reduction == TRUE) {
      vec <- abs(apply(x, 2, cor, y))
      if (all(is.na(vec))) {
        return_obj = list(ind.name = name[col],
                          dep.name = NA,
                          coefficient = 0)
      } else{
        x = x[,order(vec, decreasing = T)[1:(n/log(n))]]
      }
    }
    
    if ( all(y==0) |  all(y==1) ) {
      return_obj = list(ind.name = name[col],
                        dep.name = NA,
                        coefficient = 0)
    } else{
      ridge_cv <- try(cv.glmnet(x = x, y = y,alpha = 0))
      if ('try-error' %in% class(ridge_cv)) {
        return_obj = list(ind.name = name[col],
                          dep.name = NA,
                          coefficient = 0)
        
      } else{
        ridge_cv <- cv.glmnet(x = x, y = y, type.measure = "mse", nfolds = 10, alpha = 0)
        best_ridge_coef <- abs(as.numeric(coef(ridge_cv, s = ridge_cv$lambda.min))[-1])
        weights = 1 / ((best_ridge_coef + 1e-3)^gamma)
        
        fit <- cv.glmnet(x = x, y = y, alpha = alpha, family = "gaussian", type.measure = "mse",
                         penalty.factor = weights,
                         nfolds = 10, keep = TRUE, thresh=1e-10, maxit=1e6)
        lasso_coef <- coef(fit, s = fit$lambda.min)
        return_obj = list(ind.name = name[col],
                          dep.name = lasso_coef@Dimnames[[1]][lasso_coef@i + 1][-1],
                          coefficient = lasso_coef@x[-1])
        if ( length(return_obj$dep.name)==0 ) {
          print(paste0("No interaction found for ", name[col]))
          tmp = cor(x,y)
          return_obj$dep.name = rownames(tmp)[which.max(abs(tmp))]
          return_obj$coefficient = tmp[which.max(abs(tmp))]*1/3
        }
        
      }
      
    }
    
    
  }
  
  
  return(return_obj)
}

qdODE_parallel <- function(result,alpha,gamma ,scaler,reduction = FALSE, thread = 12, maxit = 1e3,bc=bc,LOP_order,lx,str){
  data = result$power_fit
  relationship = lapply(1:nrow(data),function(c)get_interaction(data, c, alpha,gamma,scaler,reduction = reduction))
  
  cat('Start qdODE test',sep="\n")
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterEvalQ(cl, {require(deSolve)})
  clusterExport(cl, c( "qdODE_ls", "qdODE_fit", "qdODE_all","get_legendre_matrix","power_equation1","qdODEmod",
                       "get_legendre_par","legendre_fit","result","relationship","maxit","power_equation","bc","qdODEmod_lgkt"), envir=environment())
  result = pblapply(1:nrow(data),function(c) qdODE_all(result = result,
                                                       relationship = relationship,
                                                       i = c,
                                                       maxit = maxit,
                                                       bc=bc,
                                                       LOP_order = LOP_order,
                                                       lx = lx
  ), cl = cl)
  stopCluster(cl)
  names(result) = rownames(data)
  names(relationship) = rownames(data)
  return_obj <- list(ode_result = result,
                     relationship = relationship)
  return(return_obj)
}
relationship = lapply(1:nrow(data),function(c)get_interaction(data, c, alpha=1,gamma=1,scaler=T,reduction = F))

compute_strict_LR <- function(i, result, relationship, init_pars=1,
                              LOP_order=6, bc=10, lx="ode", new_time=NULL,
                              n_expand=100, maxit=1000){
  Time <- as.numeric(colnames(result$power_fit))
  variable <- c(relationship[[i]]$ind.name, relationship[[i]]$dep.name)
  data <- result$power_fit[variable, ]
  n <- nrow(data)
  power_par <- result$power_par[variable, , drop = FALSE][-1, , drop = FALSE]
  pars_int <- c(init_pars, relationship[[i]]$coefficient)
  
  H1_fit <- qdODE_all(result, relationship, i, init_pars=init_pars,
                      LOP_order=LOP_order, methods="ls",
                      bc=bc, lx=lx, maxit=maxit)
  pred_H1 <- H1_fit$predict
  y_obs <- pred_H1[,2]
  y_fit_H1 <- pred_H1[,3]

  pars_H0 <- pars_int[1]
  H0_fit <- qdODE_fit(pars=pars_H0, data=data[1,,drop=FALSE], power_par=power_par,
                      Time=Time, bc=bc, ind_par=result$power_par[variable,][1,],
                      lx=lx)
  
  y_fit_H0 <- H0_fit$predict[,3]
  
  sigma2_H1 <- var(y_obs - y_fit_H1)
  sigma2_H0 <- var(y_obs - y_fit_H0)
  
  logL_H1 <- sum(dnorm(y_obs, mean=y_fit_H1, sd=sqrt(sigma2_H1), log=TRUE))
  logL_H0 <- sum(dnorm(y_obs, mean=y_fit_H0, sd=sqrt(sigma2_H0), log=TRUE))
  
  LR <- 2 * (logL_H1 - logL_H0)
  
  return(list(LR=LR, H1_fit=H1_fit, H0_fit=H0_fit))
}
library(parallel)

LR_test_strict <- function(result, relationship, nperm=1000, ncores=12,
                           init_pars=1, LOP_order=6, bc=10, lx="ode",
                           maxit=1000){
  
  nvar <- length(relationship)
  LR_values <- numeric(nvar)

  cl <- makeCluster(ncores)
  clusterEvalQ(cl, {
    library(deSolve)
    library(orthopolynom)
    library(mvtnorm)
  })
  clusterExport(cl, varlist=c("result", "relationship",
                              "compute_strict_LR",
                              "init_pars","LOP_order","bc","lx","maxit",
                              "qdODE_ls","qdODE_fit","qdODE_all","get_legendre_matrix"
                              ,"power_equation1","qdODEmod","get_legendre_par","legendre_fit","power_equation","qdODEmod_lgkt"), 
                
                envir=environment())
  LR_list <- parLapply(cl, 1:nvar, function(i){
    compute_strict_LR(i, result, relationship,
                      init_pars=init_pars, LOP_order=LOP_order,
                      bc=bc, lx=lx, maxit=maxit)$LR
  })
  print("Finished computing LR values.")
  print(LR_list)
  stopCluster(cl)
  
  LR_values <- unlist(LR_list)

  if(nperm > 0){
    p_values <- numeric(nvar)
    
    cl <- makeCluster(ncores)
    clusterEvalQ(cl, {
      library(deSolve)
      library(orthopolynom)
      library(mvtnorm)
      library(pbapply)
    })
    clusterExport(cl, varlist=c("result", "relationship",
                                "compute_strict_LR",
                                "init_pars","LOP_order","bc","lx","maxit",
                                "qdODE_ls","qdODE_fit","qdODE_all","get_legendre_matrix"
                                ,"power_equation1","qdODEmod","get_legendre_par","legendre_fit","power_equation","qdODEmod_lgkt","LR_values","nperm"), 
                  envir=environment())
    

    LR_perm_all <- pblapply(1:nvar, function(i){
      obs_data <- result$power_fit[relationship[[i]]$ind.name, ]
      LR_perm <- numeric(nperm)
      for(p in 1:nperm){
        print(paste("Permutation", p, "for variable", relationship[[i]]$ind.name))
        perm_data <- obs_data[sample(length(obs_data))]
        result_perm <- result
        result_perm$power_fit[relationship[[i]]$ind.name, ] <- perm_data
        LR_perm[p] <- compute_strict_LR(i, result_perm, relationship,
                                        init_pars=init_pars, LOP_order=LOP_order,
                                        bc=bc, lx=lx, maxit=maxit)$LR
      }
      LR_perm
    }, cl = cl)
    # LR_values[i] > quantile(LR_perm, 0.95, names = FALSE)
    
    stopCluster(cl)
    thresold <- numeric(nvar)
    significant <- sapply(1:nvar, function(i){
      LR_values[i] > quantile(LR_perm_all[[i]], 0.95, names = FALSE)
      thresold[i] <- quantile(LR_perm_all[[i]], 0.95, names = FALSE)
    })
    thresold <- sapply(1:nvar, function(i){
      quantile(LR_perm_all[[i]], 0.95, names = FALSE)
    })
    return(list(
      variable = sapply(relationship, function(x) x$ind.name),
      LR = LR_values,
      LR_perm = LR_perm_all,
      significant = significant,
      threshold = thresold
    ))
  } else {
    return(list(
      variable = sapply(relationship, function(x) x$ind.name),
      LR = LR_values,
      LR_perm = NULL,
      significant = rep(NA, nvar)
    ))
  }
}
result = LR_test_strict(result = data, relationship = relationship, nperm=1000, ncores=12,
               init_pars=1, LOP_order=6, bc=10, lx="ode",
               maxit=100)

result$threshold = result$significant
result$significant = result$LR > result$threshold
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

df_perm <- map_dfr(seq_along(result$variable), function(i) {
  data.frame(
    variable = result$variable[i],
    LR = result$LR[i],
    threshold = result$threshold[i],
    significant = ifelse(result$LR_perm[[i]] > result$threshold[i], TRUE, FALSE),
    LR_perm_value = result$LR_perm[[i]]
  )
})

df_obs <- data.frame(
  variable = result$variable,
  LR = result$LR,
  threshold = result$threshold,
  significant = result$significant
)

df_obs$significant <- factor(df_obs$significant, levels = c(TRUE, FALSE))

extract_module_number <- function(x) {
  as.numeric(gsub("M", "", x))
}

df_perm$variable <- factor(df_perm$variable, 
                           levels = unique(df_perm$variable)[order(extract_module_number(unique(df_perm$variable)))])

df_obs$variable <- factor(df_obs$variable, 
                          levels = unique(df_obs$variable)[order(extract_module_number(unique(df_obs$variable)))])

p <- ggplot(df_perm, aes(x = variable, y = LR_perm_value)) +
  geom_boxplot(
    fill = "#f7f7f7", 
    color = "#525252", 
    alpha = 0.8,
    outlier.shape = 21, 
    outlier.fill = "#f03b20", 
    outlier.color = "#525252",
    outlier.size = 1.5,
    outlier.alpha = 0.6,
    width = 0.7
  ) +
  geom_point(
    data = df_obs, 
    aes(x = variable, y = LR, color = significant, fill = significant), 
    shape = 21, 
    size = 3.5,
    stroke = 1
  ) +
  geom_point(
    data = df_obs, 
    aes(x = variable, y = threshold), 
    color = "#2b8cbe", 
    shape = 4, 
    size = 3, 
    stroke = 1.2
  ) +
  scale_color_manual(
    values = c("TRUE" = "#e34a33", "FALSE" = "#2c7fb8"),
    labels = c("TRUE" = "Significant", "FALSE" = "Non-significant"),
    name = "Statistical Significance"
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#fdbb84", "FALSE" = "#7fcdbb"),
    labels = c("TRUE" = "Significant", "FALSE" = "Non-significant"),
    name = "Statistical Significance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "#f0f0f0"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#e0e0e0", fill = NA, linewidth = 0.5)
  ) +
  labs(
    y = "Likelihood Ratio (LR)",
    x = "Module",
    title = "Comparison of Observed and Permutation-based Likelihood Ratios",
    # subtitle = "With 95% significance thresholds",
    caption = "Red points: observed LR values\nBlue crosses: 95% significance thresholds\nBoxplots: permutation distribution"
  )

print(p)
ggsave("/modulLR_permutation_test_plot.png", plot = p, width = 10, height = 6, dpi = 300)
