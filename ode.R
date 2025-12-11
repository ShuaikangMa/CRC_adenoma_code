rm(list = ls())
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
  nn = length(Pars)
  ind_effect = paste0("alpha","*",names(State)[1])
  dep_effect = sapply(2:nn, function(c) paste0(paste0("beta",c-1),"*",names(State)[c]))
  dep_effect = paste0(dep_effect, collapse = "+")
  all_effect = paste0(ind_effect, "+", dep_effect)
  expr = parse(text = all_effect)
  
  with(as.list(c(State, Pars)), {

    dx = eval(expr)
    dy <- power_par[,1]*power_par[,2]*Time^(power_par[,2]-1)
    dind = alpha*x
    for(i in c(1:(nn-1))){
      tmp = paste0(paste0("beta",i),"*",paste0("y",i))
      expr2 = parse(text = tmp)
      assign(paste0("ddep",i),eval(expr2))
    }
    return(list(c(dx, dy, dind, mget(paste0("ddep",1:(nn-1))))))
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
    #Time <- seq(min(Time), max(Time), by = ((max(Time)-min(Time))/length(Time))/bc)
    out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                            times = Time, power_par = power_par))

    X = as.numeric(data[1,])
    fit = as.numeric(out[,2])
    ind = as.numeric(out[,(n+2)])
    sse = sum(crossprod(X-fit),sum((ind[ind<0])^2))
  }
  if(lx == "lgkt"){
    # RK4 implementation
    rk4 <- function(ode_func, State, Time, Pars, power_par, dt) {
      k1 <- ode_func(Time, State, Pars, power_par)
      k2 <- ode_func(Time + 0.5 * dt, State + 0.5 * dt * k1, Pars, power_par)
      k3 <- ode_func(Time + 0.5 * dt, State + 0.5 * dt * k2, Pars, power_par)
      k4 <- ode_func(Time + dt, State + dt * k3, Pars, power_par)
      
      # Update state using RK4 formula
      new_state <- State + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      
      return(new_state)
    }
    
    # Time integration function
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
    out2 = data.frame(x = out[,1], y = as.numeric(data[1,]), y.fit = out[,2],
                      ind = out[,(n+2)], dep = out[,(n+3):(ncol(out))])
  }
  if(lx ==  "lgkt"){
    # RK4 implementation
    rk4 <- function(ode_func, State, Time, Pars, power_par, dt) {
      k1 <- ode_func(Time, State, Pars, power_par)
      k2 <- ode_func(Time + 0.5 * dt, State + 0.5 * dt * k1, Pars, power_par)
      k3 <- ode_func(Time + 0.5 * dt, State + 0.5 * dt * k2, Pars, power_par)
      k4 <- ode_func(Time + dt, State + dt * k3, Pars, power_par)
      
      # Update state using RK4 formula
      new_state <- State + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      
      return(new_state)
    }
    
    # Time integration function
    integrate_rk4 <- function(ode_func, State, Time_range, Pars, power_par, dt) {
      times <- seq(Time_range[1], Time_range[2], by = dt)
      #times = Time
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

  colnames(out2)[4:ncol(out2)] = c(rownames(data)[1], rownames(data)[2:n])
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
qdODE_parallel <- function(result,alpha,gamma ,scaler,reduction = FALSE, thread = 12, maxit = 1e3,bc=bc,LOP_order,lx,str){
  
  # result = k_convert$a;alpha=1;gamma = 1;bc= 5;LOP_order=6;lx = "ode";scaler = T;reduction = F
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

get_interaction <- function(data, col,alpha,gamma, scaler,reduction = FALSE ){
  #col = 1;data = fit_result$power_fit
  # data = all_result$a4$power_fit;col = 6;alpha = 1;gamma = 1;scaler = T;reduction = F

  if (nrow(data)==2) {
    return_obj = list(ind.name = rownames(data)[col],
                      dep.name = rownames(data)[-col],
                      coefficient = cor(t(data))[1,2])
    
  } else{
    data <- t(data); name <- colnames(data)
    if(scaler == T){
      data = scale(data)
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


fun_clu_select <- function(result_fit, result_funclu, i){

  result_funclu$Module.all <- result_funclu$Module.all[order(as.numeric(names(result_funclu$Module.all)))]

  names(result_funclu$Module.all) <- as.character(as.numeric(names(result_funclu$Module.all)) + 1)
  
  cluster.name = rownames(result_funclu$Module.all[[i]])
  original_data = result_fit$original_data[cluster.name,]
  
  times = result_fit$Time
  times_new = seq(min(times),max(times),length = 30)
  
  par.mu = result_fit$power_par[cluster.name,]
  mu.fit = result_fit$power_fit[cluster.name,]
  
  return_obj <- list(original_data = original_data,
                     trans_data = mu.fit,
                     power_par = par.mu,
                     power_fit = mu.fit)
  return(return_obj )
}

qdODEplot_convert <- function(result){

  data = result$predict
  n = ncol(data)
  colnames(data)[4:n] = c(paste0("ind.",colnames(data)[4]),
                          paste0("dep.",colnames(data)[5:n]))
  plot.df = melt(data, id.vars = c("x"))
  name = levels(plot.df[,2])
  ind.name = name[grep("^ind", name)]
  ind.name2 <- sub("^ind\\.", "", ind.name)
  ind.df <- subset(plot.df, plot.df[,2] == ind.name)
  ind.df$type = "ind"
  ind.df$variable = ind.name2
  
  depname = levels(plot.df[,2])[grep("^dep",name )]
  dep.df <- subset(plot.df, plot.df[,2] %in% depname)
  dep.df$type = "dep"

  dep.df$variable <- sub("^dep\\.", "", as.character(dep.df$variable))
  original.df = subset(plot.df, plot.df[,2] == "y")
  original.df$type = "original"
  fit.df = subset(plot.df, plot.df[,2] == "y.fit")
  fit.df$type = "fit"
  plot.df2 = rbind(ind.df, dep.df,fit.df)
  
  name.df = subset(plot.df2, plot.df[,1] == max(plot.df2[,1]))
  name.df = name.df[-nrow(name.df),]
  name.df[,2][name.df[,2] == "y.fit"] = ind.name2
  name.df = name.df[-which(name.df[,4] == "fit"),]
  name.df[,1] = name.df[,1]*1.000
  return_obj = list(plot.df2 = plot.df2,
                    name.df = name.df,
                    ind.name2 = ind.name2)
  return(return_obj)
}

qdODE_plot_base <- function(result,label = 10, show.legend = TRUE,fit_result,point){
  
  result2 = qdODEplot_convert(result)
  plot.df2 = result2$plot.df2

  rownames(fit_result$trans_data) <- gsub(" ", ".", rownames(fit_result$trans_data))
  rownames(fit_result$trans_data) <- gsub(":", ".", rownames(fit_result$trans_data))
  rownames(fit_result$trans_data) <- gsub("\\(", ".", rownames(fit_result$trans_data))
  rownames(fit_result$trans_data) <- gsub("\\)", ".", rownames(fit_result$trans_data))
  scatter.df <- data.frame(x = as.numeric(colnames(fit_result$trans_data)), value = as.numeric(fit_result$trans_data[result2$ind.name2,]))
  name.df = result2$name.df
  ind.name2 = result2$ind.name2
  if(nchar(ind.name2) > 1){
    p <- ggplot(plot.df2, mapping = aes_string(x = "x", y = "value")) +
      geom_line(mapping = aes_string(group = "variable", colour = "type"), size = 1.1,
                show.legend = show.legend) +
      geom_label_repel(
        data = name.df,
        aes(label = variable, colour = type, x = x, y = value),
        show.legend = FALSE,
        direction = "y",          
        fill = NA,                
        box.padding = 0.5,        
        point.padding = 0.3,      
        segment.color = "grey50",
        segment.size = 0.3,       
        min.segment.length = 0.2, 
        label.size = 0.3,
        label.r = unit(0.15, "lines"),
        size = 3,
        nudge_y = 0.1,            
        max.overlaps = 20       
      ) + 
      # ggrepel::geom_text_repel(data = name.df,
      #                          mapping = aes_string(label = "variable", colour = "type", x = "x", y = "value"),
      #                          show.legend = FALSE, direction = "y",
      #                          box.padding = 0.25, segment.color = "grey50",
      #                          size = 2, max.overlaps = Inf) +
      geom_hline(yintercept = 0, size = 0.5, linetype = 'dashed') +
      scale_color_manual(
        name = "Effect Type",
        labels = c("Dep Effect" ,"All Effect", "Ind Effect"),
        values = c("#018A67", "#293890","#BF1D2D")) +
      xlab("Habitat Index") + ylab("Abundance of Individual Modules") +
      ggtitle(result2$ind.name2) + theme_bw() +

      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 0), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(color = "black", fill = NA)  
        # plot.margin = margin(0, 0, 0, 0)
      )
    if (!is.null(point)) {
      p <- p + geom_point(data = scatter.df, 
                          mapping = aes(x = x, y = value), 
                          size = 1, color = "#293890", shape = 16, alpha = 0.4) 
    }
  }
  
  if (is.null(label)) {
    p <- p + scale_x_continuous(limits = c(min(plot.df2$x), max(plot.df2$x) * 1.005))
  } else {
    p.build <- ggplot_build(p)
    xbreaks <- p.build$layout$panel_params[[1]]$x.sec$breaks
    ybreaks <- p.build$layout$panel_params[[1]]$y.sec$breaks
    
    if (length(ybreaks) > 0) {
      ytext <- rep(NA, length(ybreaks))
      for (i in seq_along(ybreaks)) {
        val <- ybreaks[i]
        if (is.na(val)) {
          ytext[i] <- "NA"
        } else if (val == 0) {
          ytext[i] <- "0"
        } else if (val < 0) {
          ytext[i] <- paste0("-", label, "^", abs(val))
        } else {
          ytext[i] <- paste0(label, "^", val)
        }
      }
      ylabel2 <- parse(text = ytext)
    } else {
      ylabel2 <- waiver()
    }
    
    if (length(xbreaks) > 0) {
      xtext <- rep(NA, length(xbreaks))
      for (i in seq_along(xbreaks)) {
        val <- xbreaks[i]
        if (is.na(val)) {
          xtext[i] <- "NA"
        } else if (val == 0) {
          new_val <- min(xbreaks[!is.na(xbreaks) & xbreaks != 0]) - 1
          xbreaks[i] <- new_val
          xtext[i] <- paste0(label, "^", new_val)
        } else if (val < 0) {
          xtext[i] <- paste0("-", label, "^", abs(val))
        } else {
          xtext[i] <- paste0(label, "^", val)
        }
      }
      xlabel2 <- parse(text = xtext)
    } else {
      xlabel2 <- waiver()
    }
    
    if (nchar(ind.name2) %in% c(2, 3)) {
      p <- p + scale_x_continuous(
        breaks = xbreaks,
        labels = xlabel2,
        limits = c(min(plot.df2$x)*0.99, max(plot.df2$x)*1.01)
      ) + scale_y_continuous(
        breaks = ybreaks,
        labels = ylabel2
      ) + theme(axis.text.y = element_text(hjust = 0))
    } else {
      p <- p + scale_x_continuous(
        breaks = xbreaks,
        labels = xlabel2,
        limits = c(min(scatter.df$x)*0.99, max(plot.df2$x)*1.01)
      ) + scale_y_continuous(
        breaks = ybreaks,
        labels = ylabel2
      ) + theme(axis.text.y = element_text(hjust = 0))
    }
  }
  return(p)
}

qdODE_plot_all <- function(result,fit_result,point,label = 10, show.legend = TRUE, nrow = NULL, ncol = NULL,title=NULL,name = NULL){
  
  if(is.null(name)){
    p = lapply(result$ode_result, qdODE_plot_base, label = label, show.legend = show.legend, fit_result = fit_result, point = point)
    ncol = 3
    nrow = 5
    }else{
    p = list()
    for (i in name) {
      print(i)
      p[[i]] <- qdODE_plot_base(result$ode_result[[i]], label = label, show.legend = show.legend, fit_result = fit_result, point = point)
    }
    ncol = 1
    nrow = length(name)
  }
  
  p = lapply(p, "+", xlab(NULL))
  p = lapply(p, "+", ylab(NULL))
  if(nchar(names(p)[1]) == 2 || nchar(names(p)[1]) == 3){
    y_lab <- ggplot() + annotate(geom = "text", size=6,x = 1, y = 1,
                                 label = "Abundance of Individual Modules", angle = 90) +
      coord_cartesian(clip = "off") + theme_void()
  }else{
    y_lab <- ggplot() + annotate(geom = "text", size=6,x = 1, y = 1,
                                 label = "Abundance of Individual Microbe", angle = 90) +
      coord_cartesian(clip = "off") + theme_void()
  }
  
  pp = (y_lab | wrap_plots(p, nrow = nrow, ncol = ncol)) +
    
    plot_annotation(caption = "Habitat Index",
                    #title = title, 
                    theme = theme(plot.caption = element_text(size = 17,hjust=.5),
                                  plot.background = element_rect(fill = "transparent", color = NA),
                                  panel.background = element_rect(fill = "transparent", color = NA))) +
    plot_layout(widths = c(.05, 1),guides = 'collect') + 
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
  return(pp)
}
