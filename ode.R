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

load("a_result_all.RData")
a_fit_result <- result_all
load("h_result_all.RData")
h_fit_result <- result_all

power_equation1 <- function(x, power_par){power_par[1]*x^power_par[2]}
power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),
                                                   function(c) power_par[c,1]*x^power_par[c,2] )) }
get_legendre_matrix <- function(x,legendre_order){
  #生成勒让德多项式系数,其中 n = legendre_order 指定勒让德多项式的阶数，normalized = F 表示不对多项式进行标准化处理
  legendre_coef <- legendre.polynomials(n = legendre_order, normalized=F)
  #polynomial.values 函数用于计算多项式在给定点 x 上的值
  #scaleX(x, u = -1, v = 1) 函数将输入数据 x 线性映射到区间 [-1, 1][−1,1] 上，因为勒让德多项式通常定义在这个区间内
  legendre_matrix <- as.matrix(as.data.frame(polynomial.values(
    polynomials = legendre_coef, x = scaleX(x, u = -1, v = 1))))
  #将列名命名为勒让德阶数
  colnames(legendre_matrix) <- paste0("legendre_",0:legendre_order)
  #只返回除了第一阶外的勒让德值
  return(legendre_matrix[,2:(legendre_order+1)])
}

get_legendre_par <- function(y,legendre_order,x){
  #lm_lx
  #y = out2[,2];legendre_order=6;x=out2$x
  #y为拟合的除时间序列外的y值和变化率，legendre_order为勒让德阶数，x为时间序列
  #将拟合的除时间序列外的y值和变化率与勒让德阶数的值进行线性拟合并提取系数
  #得到截距项与阶数数量相同个数的参数项
  legendre_par <- as.numeric(coef(lm(y~get_legendre_matrix(x,legendre_order))))
  return(legendre_par)
}

legendre_fit <- function(par,x){
  #par = all_LOP_par[,1];x = time2
  #获取线性拟合勒让德的参数数量
  legendre_order = length(par)
  #生成勒让德多项式的线性组合，用于拟合数据
  #勒让德多项式的线性组合，每个多项式的系数由 par 中相应的元素决定
  fit <- sapply(1:length(par), function(c)
    par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
  #生成一个矩阵，每一列表示一个勒让德多项式在输入数据 x 上的值
  legendre_fit <- as.matrix(as.data.frame(polynomial.values(
    polynomials = fit, x = scaleX(x, u = -1, v = 1))))
  #将所有勒让德多项式的值加在一起，得到最终的拟合值或插值值。每一行的和代表了在对应 x 值上的拟合结果
  x_interpolation <- rowSums(legendre_fit)
  #返回最终的拟合值总和
  return(x_interpolation)
}

darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

qdODEmod <- function(Time, State, Pars, power_par) {
  #统计参数1与互作系数的总个数
  nn = length(Pars)
  #设置独立节点名称为"alpha*x"
  ind_effect = paste0("alpha","*",names(State)[1])
  #设置依赖节点的名称为"beta_n*y_n"
  dep_effect = sapply(2:nn, function(c) paste0(paste0("beta",c-1),"*",names(State)[c]))
  #将单个依赖节点名称按+进行组合，例如："beta1*y1+beta2*y2+beta3*y3"
  dep_effect = paste0(dep_effect, collapse = "+")
  #将独立节点名称与总的依赖节点名称进行合并，例如："alpha*x+beta1*y1+beta2*y2+beta3*y3"
  all_effect = paste0(ind_effect, "+", dep_effect)
  #将组合的总体名称解析为表达式
  expr = parse(text = all_effect)
  
  with(as.list(c(State, Pars)), {
    # eval执行parse解析后的表达式或语句
    # 计算独立节点变量拟合值与依赖节点变量拟合值乘以各自系数后的总和值即总体变化率
    dx = eval(expr)
    #计算依赖节点随时间的变化率，该公式由a*t^b对t求导而来
    dy <- power_par[,1]*power_par[,2]*Time^(power_par[,2]-1)
    #计算独立节点变量的变化率
    dind = alpha*x
    #计算每一个依赖节点的独立变化率
    for(i in c(1:(nn-1))){
      tmp = paste0(paste0("beta",i),"*",paste0("y",i))
      expr2 = parse(text = tmp)
      #将每一个计算结果分别赋值给ddep_i
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
  #pars = pars_int;power_par = power_par
  
  #计算设定的参数1与互作系数的总数量
  n = length(pars)
  #将fit的参数ab转换为矩阵形式
  power_par = as.matrix(power_par)
  #判断与节点变量互作变量的个数是否为1
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    #设置参数1为alpha，设置互作系数为beta
    Pars = c(alpha = pars[1], beta = pars[2:n])
    #设置节点变量和 互作变量拟合数据的第一行第一列的数据为x
    #设置除了第一行的第一列的其余数据为y，其格式为纵向
    #设置节点数据为节点变量拟合数据的第一个值
    #设置依赖数据为互作变量个数的零值
    #data为提取的独立变量与依赖变量的幂函数拟合值
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  #求解qdODEmod的微分方程
  #y = State：这是ODE系统的初始状态。State包含每个状态变量的初始值。
  #times = Time：这是时间点的向量，表示在哪些时间点求解ODE系统
  #qdODEmod 接受时间（t）、状态变量（y）、参数（parms）和其他可能的参数（如 power_par），并返回状态变量关于时间的导数（即 dy/dt）
  # out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
  #                         times = Time, power_par = power_par))
  if(lx ==  "ode"){
    #Time <- seq(min(Time), max(Time), by = ((max(Time)-min(Time))/length(Time))/bc)
    out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                            times = Time, power_par = power_par))
    #设置独立节点变量的fit拟合值即真实值为X
    X = as.numeric(data[1,])
    #设置独立节点变量的ode拟合值为fit
    fit = as.numeric(out[,2])
    #设置独立节点变量的ode拟合随时间变化率即独立效应拟合值为ind
    ind = as.numeric(out[,(n+2)])
    #计算误差项，由真实值 X 与模型拟合值 fit 的平方误差和独立效应小于0的平方误差组成
    #在许多实际问题中，独立效应（例如种群密度、资源量等）通常不应为负数。负值在这些情况下是没有物理意义的。因此，模型在求解过程中生成负值是不可取的
    #尽可能减少负值衰减的出现
    
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
    #设置独立节点变量的fit拟合值即真实值为X
    X = as.numeric(data[1,])
    #设置独立节点变量的ode拟合值为fit
    fit = as.numeric(out[,2])
    #设置独立节点变量的ode拟合随时间变化率即独立效应拟合值为ind
    ind = as.numeric(out[,(n+2)])
    #计算误差项，由真实值 X 与模型拟合值 fit 的平方误差和独立效应小于0的平方误差组成
    #在许多实际问题中，独立效应（例如种群密度、资源量等）通常不应为负数。负值在这些情况下是没有物理意义的。因此，模型在求解过程中生成负值是不可取的
    #尽可能减少负值衰减的出现
    # alpha=5e-3
    # beta=1e-3
    alpha=1e-1
    beta=1e-1
    sse = sum(crossprod(X-fit),sum((ind[ind<0])^2),alpha*(Pars[1]^2),beta*(sum(Pars[2:n])^2))
    # sse = sum(crossprod(X-fit),sum((ind[ind<0])^2))
  }
  
  
  return(sse)
}

qdODE_fit <- function(pars, data, Time, power_par, LOP_order = 6, new_time = NULL, n_expand = 100,bc,ind_par,lx){
  #获取独立节点变量与依赖节点变量的总个数
  #pars = qdODE.est$par
  n = length(pars)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    #设置参数1为alpha，设置互作系数为beta
    Pars = c(alpha = pars[1], beta = pars[2:n])
    #设置节点变量和互作变量拟合数据的第一行第一列的数据为x
    #设置除了第一行的第一列的其余数据为y，其格式为纵向
    #设置节点数据为节点变量拟合数据的第一个值
    #设置依赖数据为互作变量个数的零值
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  if(lx ==  "ode"){
    #Time <- seq(min(Time), max(Time), by = ((max(Time)-min(Time))/length(Time))/bc)
    #根据所得参数求解ode
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
  
  # #整理数据，令时间序列为x，独立节点变量的真实值作为y，独立节点变量经ode拟合后的值为y.fit
  #令ode求解后的独立节点变量的变化率值为ind，依赖节点变量的变化率为dep
  
  
  #对整理后的数据进行列命名为独立变量与依赖变量的原始名称
  colnames(out2)[4:ncol(out2)] = c(rownames(data)[1], rownames(data)[2:n])
  #设置行名为空
  rownames(out2) = NULL
  
  #进行勒让德多项式求解
  #最终得到各个拟合值的勒让德求解值
  all_LOP_par = sapply(2:ncol(out2),function(c)get_legendre_par(out2[,c], LOP_order, out2$x))
  
  #判断是否指定新时间名称
  if (is.null(new_time)) {
    #设置新时间序列名称
    time2 = seq(min(Time), max(Time), length = n_expand)
    #将求解到的各阶勒让德系数进行拟合
    out3 = apply(all_LOP_par, 2, legendre_fit, x = time2)
    #将新生成的时间序列与拟合各个阶的总值合并
    out3 = cbind(time2, out3)
  } else{
    out3 = apply(all_LOP_par, 2, legendre_fit, x = new_time)
    out3 = cbind(new_time, out3)
  }
  #将更新后的拟合值的列名重新命名为原始列名
  colnames(out3) = colnames(out2)
  
  # out_plot = out3[,-2]
  # colnames(out_plot)[1] = "time"
  # #  可视化结
  # out_long <- melt(as.data.frame(out_plot), id.vars = "time")
  # 
  # ggplot(out_long, aes(x = time, y = value, color = variable)) +
  #   geom_line() +
  #   labs(title = "ODE Solution", x = "Time", y = "State Variables") +
  #   theme_minimal()
  result = list(fit = out2,
                predict = data.frame(out3),
                LOP_par = all_LOP_par)
  return(result)
}

qdODE_all <- function(result, relationship, i, init_pars = 1, LOP_order = LOP_order, methods = "ls",
                      new_time = NULL, n_expand = 100, maxit = 1e3,bc = bc,lx){
  #i = 1;init_pars = 1; LOP_order = 6; methods = "ls"; new_time = NULL; n_expand = 100; maxit = 1e3;lx = "ode"
  #获取拟合的时间序列
  Time = as.numeric(colnames(result$power_fit))
  #获取节点变量与相应互作变量的名称
  variable = c(relationship[[i]]$ind.name, relationship[[i]]$dep.name)
  #获取节点变量与相应互作变量对应的拟合数据
  data = result$power_fit[variable,]
  #判断是否有互作关系
  if (length(variable)<=1) {
    qdODE.est = NA
    result = NA
    return.obj <- append(result, list(ODE.value = NA,
                                      parameters = NA))
  } else{
    #获取互作变量的ab参数值
    power_par = result$power_par[variable,][-1,]
    #获取节点变量与互作变量的总数
    n = nrow(data)
    #获取互作变量与节点变量之间的相关系数
    pars_int = c(init_pars,relationship[[i]]$coefficient)
    #判断拟合方式并用optim进行拟合
    if (methods == "ls") {
      #利用optim对qdODE_ls的sse进行拟合使其最小化，并获得优化后的相关系数参数
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,bc=bc,lx = lx,
                         method = "L-BFGS-B",
                         lower = c(rep(-10,(length(pars_int)))),
                         upper = c(rep(10,(length(pars_int)))),
                         control = list(trace = TRUE, maxit = maxit))
      #利用勒让德多项式进行求解拟合数值
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
  #将结果名修改为原始物种名
  names(result) = rownames(data)
  #将关系拟合结果名修改为原始物种名
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
                          #dep.name 变量包含了被选中的自变量的名称
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

load("k_convert.RData")

a1<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=1)
h1<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=1)
a2<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=2)
h2<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=2)
a3<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=3)
h3<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=3)
a4<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=4)
h4<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=4)
a5<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=5)
h5<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=5)
a6<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=6)
h6<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=6)
a7<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=7)
h7<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=7)
a8<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=8)
h8<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=8)
a9<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=9)
h9<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=9)
a10<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=10)
h10<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=10)
a11<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=11)
h11<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=11)
a12<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=12)
h12<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=12)
a13<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=13)
h13<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=13)
a14<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=14)
h14<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=14)
a15<-fun_clu_select(result_fit = a_fit_result,result_funclu = k_convert$a,i=15)
h15<-fun_clu_select(result_fit = h_fit_result,result_funclu = k_convert$b,i=15)


qdODE_parallel_result_a<-qdODE_parallel(k_convert$a, reduction = F, scaler = T,
                                       thread = 12, maxit = 1e3,alpha = 1,gamma = 1,bc= 5,LOP_order=6,lx = "ode")
qdODE_parallel_result_h<-qdODE_parallel(k_convert$b, reduction = F, scaler = T,
                                        thread = 12, maxit = 1e3,alpha = 1,gamma = 1,bc= 5,LOP_order=6,lx = "ode")



ode.M1_a = qdODE_parallel(a1,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M1_h = qdODE_parallel(h1,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)

ode.M2_a = qdODE_parallel(a2,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M2_h = qdODE_parallel(h2,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)

ode.M3_a = qdODE_parallel(a3,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M3_h = qdODE_parallel(h3,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)


ode.M4_a = qdODE_parallel(a4,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M4_h = qdODE_parallel(h4,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)


ode.M5_a = qdODE_parallel(a5,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M5_h = qdODE_parallel(h5,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)

ode.M6_a = qdODE_parallel(a6,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M6_h = qdODE_parallel(h6,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)


ode.M7_a = qdODE_parallel(a7,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M7_h = qdODE_parallel(h7,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)


ode.M8_a = qdODE_parallel(a8,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M8_h = qdODE_parallel(h8,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)


ode.M9_a = qdODE_parallel(a9,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M9_h = qdODE_parallel(h9,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)

ode.M10_a = qdODE_parallel(a10,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M10_h = qdODE_parallel(h10,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)

ode.M11_a = qdODE_parallel(a11,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M11_h = qdODE_parallel(h11,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)

ode.M12_a = qdODE_parallel(a12,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M12_h = qdODE_parallel(h12,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)

ode.M13_a = qdODE_parallel(a13,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = F,reduction = F)
ode.M13_h = qdODE_parallel(h13,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = F,reduction = F)

ode.M14_a = qdODE_parallel(a14,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = F,reduction = F)
ode.M14_h = qdODE_parallel(h14,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = F,reduction = F)

ode.M15_a = qdODE_parallel(a15,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)
ode.M15_h = qdODE_parallel(h15,alpha=1,gamma = 1,bc= 5,LOP_order=6,lx = "ode",scaler = T,reduction = F)

qdODEplot_convert <- function(result){
  #提取经过Legendre拟合后的预测数据
  data = result$predict
  #计算预测值的总列数
  n = ncol(data)
  #对独立变量与依赖变量的变化率的列名进行重命名
  colnames(data)[4:n] = c(paste0("ind.",colnames(data)[4]),
                          paste0("dep.",colnames(data)[5:n]))
  #对数据以x值作为依据进行展开为三列，分别是x值，变量名以及变量值
  plot.df = melt(data, id.vars = c("x"))
  #取变量名列的各重复变量名的唯一值，即除了x以外的所有变量名
  name = levels(plot.df[,2])
  #提取名称中含有ind的作为独立变量名
  ind.name = name[grep("^ind", name)]
  #删除独立变量名前加的ind. 只保留单纯的独立变量
  # ind.name2 = strsplit(ind.name,split = "\\.")[[1]][2]
  ind.name2 <- sub("^ind\\.", "", ind.name)
  #提取数据值转为一列后的数据值中独立变量对应的数据
  ind.df <- subset(plot.df, plot.df[,2] == ind.name)
  #添加独立变量数据框中一列类型为ind
  ind.df$type = "ind"
  #将变量名一列的名称修改为独立变量的真实名称
  ind.df$variable = ind.name2
  
  #提取数据框中的依赖变量名
  depname = levels(plot.df[,2])[grep("^dep",name )]
  #提取与依赖变量名对应的相关数据
  dep.df <- subset(plot.df, plot.df[,2] %in% depname)
  #为依赖变量数据框添加类别为dep
  dep.df$type = "dep"
  #将依赖变量数据框的名称修改为原依赖变量名
  # dep.df$variable = sapply(strsplit(as.character(dep.df$variable),"\\."),"[",2)
  dep.df$variable <- sub("^dep\\.", "", as.character(dep.df$variable))
  #提取经过Legendre拟合后的独立变量的真实值作为初始值
  original.df = subset(plot.df, plot.df[,2] == "y")
  #将该值定义为初始值
  original.df$type = "original"
  #提取经过Legendre拟合后的幂函数拟合值作为fit值
  fit.df = subset(plot.df, plot.df[,2] == "y.fit")
  #将该数据库类型定义为fit
  fit.df$type = "fit"
  #将Legendre拟合后的独立变量变化率依赖变量变化率以及幂函数拟合值合并
  plot.df2 = rbind(ind.df, dep.df,fit.df)
  
  #提取在时间序列最大值时，各变量的名称数据以及类别
  name.df = subset(plot.df2, plot.df[,1] == max(plot.df2[,1]))
  #删除最后一行空值
  name.df = name.df[-nrow(name.df),]
  #将类别为fit的名称修改为独立变量名
  name.df[,2][name.df[,2] == "y.fit"] = ind.name2
  #删除类别为fit的行
  name.df = name.df[-which(name.df[,4] == "fit"),]
  #？？？？第一列最高时间值乘以1.002
  name.df[,1] = name.df[,1]*1.000
  #返回整理后包含独立变量和依赖变量变化率以及幂函数拟合值的数据框
  #最大值时间点时的独立变量与依赖变量的名称数值与类别
  #独立变量名称
  return_obj = list(plot.df2 = plot.df2,
                    name.df = name.df,
                    ind.name2 = ind.name2)
  return(return_obj)
}

qdODE_plot_base <- function(result,label = 10, show.legend = TRUE,fit_result,point){
  #label = 10; show.legend = TRUE; nrow = NULL; ncol = NULL;point = T;result = qdODE_parallel_result$ode_result$M1
  # label = 10; show.legend = TRUE; nrow = NULL; ncol = NULL;point = T;
  # result = result$ode_result$Lacticaseibacillus
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
    
    # ---- y轴label ----
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
    
    # ---- x轴label ----
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
    
    # ---- 最终添加 ----
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
  #result = qdODE_parallel_result_a;label = 10; show.legend = TRUE; nrow = NULL; ncol = NULL;fit_result = a_fit_result;point = T
  #result = ode.M12;label = 10; show.legend = TRUE; nrow = NULL; ncol = NULL;fit_result = a_fit_result;point = T
  # name = C("M1","M2","M3","M4","M5")
  # name = NULL
  # p = lapply(result$ode_result, qdODE_plot_base, label = label, show.legend = show.legend,fit_result = fit_result,point = point)
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

qdODE_a = qdODE_plot_all(qdODE_parallel_result_a,fit_result = a_fit_result,point = T,name = c("M11","M12","M13","M14","M15"))
qdODE_h = qdODE_plot_all(qdODE_parallel_result_h,fit_result = h_fit_result,point = T,name = c("M11","M12","M13","M14","M15"))
ggsave("a_model.png", plot = qdODE_a, width = 4, height = 10, dpi = 300
       ,limitsize = FALSE,bg = "transparent")
ggsave("h_model.png", plot = qdODE_h, width = 4, height = 10, dpi = 300
       ,limitsize = FALSE,bg = "transparent")


m1 = qdODE_plot_all(ode.M1_a,a_fit_result,point = T)
ggsave("a_m1_plot.png", plot = m1, width = 30, height = 18, dpi = 300,limitsize = FALSE)
m1 = qdODE_plot_all(ode.M1_h,h_fit_result,point = T)
ggsave("h_m1_plot.png", plot = m1, width = 30, height = 18, dpi = 300,limitsize = FALSE)

m2 = qdODE_plot_all(ode.M2_a,a_fit_result,point = T)
ggsave("a_M2_plot.png", plot = m2, width = 20, height = 12, dpi = 300,limitsize = FALSE)
m2 = qdODE_plot_all(ode.M2_h,h_fit_result,point = T)
ggsave("h_M2_plot.png", plot = m2, width = 20, height = 12, dpi = 300,limitsize = FALSE)


m3 = qdODE_plot_all(ode.M3_a,a_fit_result,point = T,name = c("Phocaeicola"))
ggsave("a_M3_allotuPhocaeicola.png", plot = m3, width = 6, height = 4, dpi = 300,limitsize = FALSE)
m3 = qdODE_plot_all(ode.M3_h,h_fit_result,point = T,name = c("Phocaeicola"))
ggsave("h_M3_allotuPhocaeicola.png", plot = m3, width = 6, height = 4, dpi = 300,limitsize = FALSE)

m4 = qdODE_plot_all(ode.M4_a,a_fit_result,point = T,name = c("Methylorubrum"))
ggsave("a_M4_allMethylorubrum.png", plot = m4, width = 6, height = 4, dpi = 300,limitsize = FALSE)
m4 = qdODE_plot_all(ode.M4_h,h_fit_result,point = T,name = c("Methylorubrum"))
ggsave("h_M4_allMethylorubrum.png", plot = m4, width = 6, height = 4, dpi = 300,limitsize = FALSE)

m12 = qdODE_plot_all(ode.M12_a,a_fit_result,point = T,name = c("Lacticaseibacillus"))
ggsave("a_M12_Lacticaseibacillus.png", plot = m12, width = 6, height = 4, dpi = 300,limitsize = FALSE)
m12 = qdODE_plot_all(ode.M12_h,h_fit_result,point = T,name = c("Lacticaseibacillus"))
ggsave("h_M12_Lacticaseibacillus.png", plot = m12, width = 6, height = 4, dpi = 300,limitsize = FALSE)


biqdODE_plot_base <- function(result1, result2, label = 10, show.legend = FALSE, remove.label = FALSE){
  resulta = qdODEplot_convert(result1)
  resultb = qdODEplot_convert(result2)
  plot.df1 = resulta$plot.df2
  name.df1 = resulta$name.df
  ind.name2 = resulta$ind.name2
  name.df1$x = name.df1$x*0.99
  
  plot.df2 = resultb$plot.df2
  name.df2 = resultb$name.df
  name.df2$x = name.df2$x*0.99
  
  lower1 = min(plot.df1[,3])
  upper1 = max(plot.df1[,3])
  lower2 = min(plot.df2[,3])
  upper2 = max(plot.df2[,3])
  y_min = round(min(lower1, lower2),1)-0.05
  y_max = round(max(upper1, upper2),1)+0.05
  
  p1 = ggplot(plot.df1, mapping = aes_string(x = "x", y = "value")) +
    geom_line(mapping = aes_string(group = "variable", colour = "type"), size = 1.1,
              show.legend = show.legend) +
    geom_text(name.df1, mapping = aes_string(label = "variable", colour = "type",
                                             x = "x", y = "value"), size = 2,show.legend = FALSE) +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    scale_color_manual(
      name = "Effect Type",
      labels = c("Dep Effect", "All Effect", "Ind Effect"),
      values = c("green", "blue", "red")) + theme_bw() +
    xlab("Habitat Index") + ylab("Niche Index") +
    theme(axis.text.y = element_text(hjust = 0)) + xlab(NULL)+
    ggtitle(ind.name2) + theme(plot.title = element_text(hjust = 1))
  
  
  s1 = p1 + scale_y_continuous(limits = c(y_min, y_max))
  
  p2 = ggplot(plot.df2, mapping = aes_string(x = "x", y = "value")) +
    geom_line(mapping = aes_string(group = "variable", colour = "type"), size = 1.1,
              show.legend = show.legend) +
    geom_text(name.df2, mapping = aes_string(label = "variable", colour = "type",
                                             x = "x", y = "value"), size = 2,show.legend = FALSE) +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    scale_color_manual(
      name = "Effect Type",
      labels = c("Dep Effect", "All Effect", "Ind Effect"),
      values = c("green", "blue", "red")) + theme_bw() +
    theme(axis.text.y = element_text(hjust = 0)) + xlab(NULL)
  
  s2 = p2 + scale_y_continuous(limits = c(y_min, y_max))
  
  if (is.null(label)) {
    p1 = p1 + scale_x_continuous(limits = c(min(plot.df1$x), max(plot.df1$x)*1.005))+
      theme(axis.title.x = element_text(vjust=1))
    p2 = p2 + scale_x_continuous(limits = c(min(plot.df2$x), max(plot.df2$x)*1.005))
  } else {
    xlabel1 = ggplot_build(s1)$layout$panel_params[[1]]$x.sec$breaks
    xlabel2 = ggplot_build(s2)$layout$panel_params[[1]]$x.sec$breaks
    
    ylabel = ggplot_build(s1)$layout$panel_params[[1]]$y.sec$breaks
    
    xlabel1.2 = parse(text= paste(label,"^", xlabel1, sep="") )
    xlabel2.2 = parse(text= paste(label,"^", xlabel2, sep="") )
    if (0 %in% ylabel) {
      pos0 = which(ylabel==0)
      text=paste(label,"^", ylabel, sep="")
      text[pos0] = "0"
      if (any(na.omit(ylabel)<0)) {
        pos1 = which(ylabel<0)
        text[pos1] = paste0("-",label,"^",abs(ylabel[pos1]))
        ylabel2 = parse(text=text)
      } else{ylabel2 = parse(text=text)}
    } else{
      ylabel2 = parse(text= paste(label,"^", ylabel, sep="") )
    }
    p1 = p1 + scale_x_continuous(labels = xlabel1.2, limits = c(min(xlabel1)-0.05, max(xlabel1)+0.05)) +
      scale_y_continuous(limits = c(y_min,y_max),labels = ylabel2) +
      theme(axis.text.y = element_text(hjust = 0))+
      theme(plot.margin = unit(c(0,0,0,0),"lines"))+
      theme(axis.title.x = element_text(vjust=1))
    
    p2 = p2 + scale_x_continuous(labels = xlabel2.2, limits = c(min(xlabel2)-0.05, max(xlabel2)+0.05)) +
      scale_y_continuous(limits = c(y_min,y_max),labels = ylabel2) +
      theme(axis.text.y = element_text(hjust = 0)) +
      theme(axis.text.y = element_blank(), axis.ticks.length.y = unit(-0.1,"cm")) +
      ylab(NULL)+
      theme(plot.margin = unit(c(0,0,0,0),"lines"))
    
  }
  
  if (remove.label == TRUE) {
    p1 = p1 + theme(axis.title=element_blank())
    p2 = p2 + theme(axis.title=element_blank())
  }
  
  pp = p1+p2
  return(pp)
}

biqdODE_plot_all <- function(result1, result2, label = 10, show.legend = FALSE,
                             remove.label = TRUE, nrow = NULL, ncol = NULL){
  n = length(result1$ode_result)
  p = lapply(1:n, function(c)
    biqdODE_plot_base(result1$ode_result[[c]], result2$ode_result[[c]],
                      label = label, show.legend = show.legend, remove.label = remove.label))
  
  p = lapply(p, "+", xlab(NULL))
  p = lapply(p, "+", ylab(NULL))
  
  y_lab <- ggplot() + annotate(geom = "text", size=7,x = 1, y = 1,
                               label = "Niche Index", angle = 90) +
    coord_cartesian(clip = "off") + theme_void()
  pp = (y_lab | wrap_plots(p, nrow = nrow, ncol = ncol)) +
    plot_annotation(caption = "Habitat Index",
                    theme = theme(plot.caption = element_text(size = 20,hjust=.5))) +
    plot_layout(widths = c(.03, 1),guides = 'collect')
  
  return(pp)
}
p<- biqdODE_plot_all(ode.M12_a, ode.M12_h, label = 10, show.legend = TRUE)
ggsave("C:/Users/Administrator/Desktop/biqdODE_M12_plot.png", plot = p, width = 50, height = 15, dpi = 300,limitsize = FALSE,bg = "transparent")
