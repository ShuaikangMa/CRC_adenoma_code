rm(list = ls())
library(ggplot2)
library(reshape2)
library(parallel)
library(mvtnorm)
library(pbapply)
library(deSolve)
library(orthopolynom)
library(glmnet)

data <- read.csv("rawdata.csv")
rownames(data) <- data[,1]
data <- data[,-1]
h_data <- data[,c(grep("^H",colnames(data)))]
a_data <- data[,c(grep("^A",colnames(data)))]


filter_by_mean_relative_abundance <- function(abundance_matrix, threshold = 0.1) {
  if (!is.matrix(abundance_matrix) && !is.data.frame(abundance_matrix)) {
    stop("Input must be a matrix or data.frame.")
  }
  
  rel_abund_mat <- sweep(abundance_matrix, 2, colSums(abundance_matrix), FUN = "/")
  
  mean_rel_abund <- rowMeans(rel_abund_mat)
  
  keep_microbes <- mean_rel_abund >= threshold

  filtered_matrix <- abundance_matrix[keep_microbes, , drop = FALSE]
  
  return(filtered_matrix)
}

h_data = filter_by_mean_relative_abundance(h_data, threshold = 0.00001)
a_data = filter_by_mean_relative_abundance(a_data, threshold = 0.00001)

data_cleaning <- function(data){
  x = round(ncol(data))
  zero_counts = rowSums(data == 0)
  keep_no = which(zero_counts <= x *0.7)
  data2 = data[keep_no,]
  return(data2)
}
h_data <- data_cleaning(h_data)
a_data <- data_cleaning(a_data)

get_common_otus <- function(otu_data1, otu_data2) {
  common_otus <- intersect(rownames(otu_data1), rownames(otu_data2))
  otu_data1_filtered <- otu_data1[common_otus, , drop = FALSE]
  otu_data2_filtered <- otu_data2[common_otus, , drop = FALSE]
  return(list(otu_data1 = otu_data1_filtered, otu_data2 = otu_data2_filtered))
}
h_a_result <- get_common_otus(otu_data1 = h_data, otu_data2 = a_data)
h_result <- h_a_result$otu_data1
a_result <- h_a_result$otu_data2
# write.csv(a_result, file = "a.csv")
# write.csv(h_result, file = "h.csv")

get_nozero_otus_result <- function(otu_data_bec,trans = log10){

  zero_counts <- rowSums(otu_data_bec == 0)
  otu_data = otu_data_bec
  site = order(colSums(otu_data))
  data_filtered = otu_data[,site]

  x = trans(colSums(data_filtered)+1)
  X <- data_filtered
  X = trans(X+1)
  
  zero_counts1 <- rowSums(X == 0)
  data_n <- X[which(zero_counts1<=50), ]
  error_function <- function(x, y, pars){
    y_fit = pars[1]*x^pars[2]
    error = sum((y-y_fit)^2)
    return(error)
  }
  
  power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),
                                                     function(c) power_par[c,1]*x^power_par[c,2] ) )}
  power_equation_base <- function(x, y){
    
    x <- as.numeric(x)
    y <- as.numeric(y)
    
    lmFit <- lm( log(y+1)  ~ log(x+1))
    coefs <- coef(lmFit)
    a <- exp(coefs[1])
    b <- coefs[2]
    pars <- optim(c(a = a, b = b), error_function, x = x, y = y, method = "L-BFGS-B")$par
    model <- try(nls(y~a*x^b,start = list(a = pars[1], b = pars[2]),
                     control = nls.control(maxiter = 10e4,tol = 1e-6, minFactor = 1e-200,warnOnly = TRUE)))
    if( 'try-error' %in% class(model)) {
      result = NULL
    }
    else{
      result = model
    }
    return(result)
  }
  power_equation_all <- function(x,y, maxit=1e4){
    result <- power_equation_base(x,y)
    iter <- 1
    while( is.null(result) && iter <= maxit) {
      iter <- iter + 1 
      try(result <- power_equation_base(x,y))
    }
    return(result)
  }
  power_equation_fit <- function(data_n,x, n=100, trans = log, thread = 12) {
    trans_data = data_n
    colnames(trans_data) = x
    
    core.number <- thread
    cl <- makeCluster(getOption("cl.cores", core.number))
    clusterExport(cl, c("power_equation_all", "power_equation_base", "trans_data", "x","error_function"), envir = environment())
    all_model = parLapply(cl = cl, 1:nrow(data_n), function(c) power_equation_all(x, trans_data[c,]))
    stopCluster(cl)
    
    
    names(all_model) = rownames(data_n)
    no = which(sapply(all_model, length)>=1)
    all_model2 = all_model[no]
    data2 = data_filtered[no,]
    trans_data2 = trans_data[no,]
    
    new_x = seq(min(x), max(x), length = n)
    power_par = t(vapply(all_model2, coef, FUN.VALUE = numeric(2), USE.NAMES = TRUE))
    power_fit = t(vapply(all_model2, predict, newdata = data.frame(x=new_x),
                         FUN.VALUE = numeric(n), USE.NAMES = TRUE))
    
    colnames(power_fit) = new_x
    result = list(original_data = data2, trans_data = trans_data2,
                  power_par = power_par, power_fit = power_fit,
                  Time = x)
    return(result)
  }
  r1 = power_equation_fit(data_n,x, thread = 12)
  
  return(r1)
}
# result1 <- get_nozero_otus_result(otu_data_bec = h_result)
result1 <- get_nozero_otus_result(otu_data_bec = a_result)
params = result1$power_par

get_zero_otus_result <- function(otu_data_bec,trans = log10){
  #otu_data_bec = a_result;trans = log10
  otu_data = otu_data_bec
  data_filtered = otu_data[,order(colSums(otu_data))]
  site = order(colSums(data_filtered))
  data_filtered = data_filtered[,site]
  
  
  x = trans(colSums(data_filtered)+1)
  X <- data_filtered
  X = trans(X+1)
  
  zero_counts1 <- rowSums(X == 0)
  data_n <- X[which((zero_counts1<=length(colnames(otu_data)))&zero_counts1>50),]
  colnames(data_n) <- x
  
  rownum = length(row.names(data_n))
  params1 <- data.frame(matrix(nrow = rownum, ncol = 5))
  model = list()
  for (i in 1:rownum) {
    i = 1
    print(i)
    y=data_n[i,]
    y0=y[which(y==0)]
    y1=y[-which(y==0)]
    
    
    x0=x[which(y==0)]
    x1=x[-which(y==0)]
    
    beta1 = c(0.1,0.1)
    beta2 = c(0.1,0.1)
    sigma = c(0.1)
    
    JointLogLik <-  function(par,x0,x1,y0,y1){
      beta1 = par[1:2]
      beta2 = par[3:4]
      sigma = par[5]
      
      p0 = exp(beta1[1]*x0^beta1[2])/(1+exp(beta1[1]*x0^beta1[2]))
      p1 = 1/(1+exp(beta1[1]*x1^beta1[2]))
      
      LogLik1 = sum(log(p0))
      LogLik2 = sum(log(p1)) + 
        sum(dnorm(y1-beta2[1]*x1^beta2[2], mean = 0, sd = abs(sigma), log = T))
      LogLik = LogLik1 + LogLik2
      return(-LogLik)
    }
    par = c(beta1,beta2,sigma)
    par_hat = optim(par,JointLogLik,x0=as.numeric(x0),x1=as.numeric(x1),y0=as.numeric(y0),y1=as.numeric(y1), 
                    method = "Nelder-Mead",
                    # method = "BFGS",
                    control = list(maxit = 2e4, 
                                   trace = T, 
                                   parscale = rep(1e-2,5),
                                   #factr = 1e-300,
                                   pgtol = 1e-300))
    if ('try-error' %in% class(par_hat) ){
      par_hat = optim(par,JointLogLik,x0=as.numeric(x0),x1=as.numeric(x1),y0=as.numeric(y0),y1=as.numeric(y1), 
                      # method = "Nelder-Mead",
                      method = "L-BFGS",
                      control = list(maxit = 8e3, 
                                     trace = T, 
                                     parscale = rep(1e-2,5),
                                     #factr = 1e-300,
                                     pgtol = 1e-300))
    }
    if (  par_hat$convergence != 0 ){
      while(par_hat$convergence != 0){
        par = par_hat$par
        par_hat = optim(par,JointLogLik,x0=as.numeric(x0),x1=as.numeric(x1),y0=as.numeric(y0),y1=as.numeric(y1), 
                        method = "Nelder-Mead", 
                        #method = "BFGS",
                        control = list(maxit = 2e3, 
                                       trace = F, 
                                       parscale = rep(1e-2,5),
                                       #factr = 1e-300,
                                       pgtol = 1e-300))
      }
    }
    # par_hat$par[c(3)] = par_hat$par[c(3)]*3/4
    model[[i]] <-  par_hat
    params1[i,] <- par_hat$par

    rownames(params1)[i] <- rownames(data_n)[i]
  }
  for (i in 1:length(model)) {
    if (model[[i]]$convergence != 0){
      print(i)
      y=data_n[i,]
      y0=y[which(y==0)]
      y1=y[-which(y==0)]
      
      x0=x[which(y==0)]
      x1=x[-which(y==0)]
      while(model[[i]]$convergence != 0){
        par = model[[i]]$par
        par_hat = optim(par,JointLogLik,x0=as.numeric(x0),x1=as.numeric(x1),y0=as.numeric(y0),y1=as.numeric(y1),
                        method = "Nelder-Mead",
                        #method = "BFGS",
                        control = list(maxit = 2e3,
                                       trace = F,
                                       parscale = rep(1e-2,5),
                                       #factr = 1e-300,
                                       pgtol = 1e-300))
        model[[i]] = par_hat
      }
      model[[i]] = par_hat
      params1[i,] = par_hat$par
      # result2 <- list(power_par = params1,model = model)
    }
    
  }
  return(list(power_par = params1, model = model))
}

# result2 <- get_zero_otus_result(otu_data_bec = h_result)
result2 <- get_zero_otus_result(otu_data_bec = a_result)
params1 = result2$power_par

get_all_otus_result <- function(otu_data_bec,params,params1,type,trans = log10){
  otu_data = otu_data_bec
  data_filtered = otu_data[,order(colSums(otu_data))]
  
  x = trans(colSums(data_filtered)+1)
  X <- data_filtered
  
  zero_counts1 <- rowSums(X == 0)
  data_n <- X[which(zero_counts1==0), ]
  data1 <- data_filtered
  data1 <- data1[c(which(zero_counts1<=50),which((zero_counts1<=length(colnames(otu_data)))&zero_counts1>50)),]
  X <- X[c(which(zero_counts1<=50),which((zero_counts1<=length(colnames(otu_data)))&zero_counts1>50)),]
  data1 <- trans(data1+1)
  
  power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),
                                                     function(c) power_par[c,1]*x^power_par[c,2] ) )}
  
  data_fit <- data.frame(matrix(nrow = length(row.names(data_filtered)), ncol = length(colnames(data_filtered))))
  data_fit[1:nrow(params),] <- power_equation(x,params)
  data_fit[(nrow(params)+1):(nrow(params)+nrow(params1)),] <- power_equation(x,params1[,3:4])
  power_par <- data.frame(matrix(nrow = nrow(data_filtered), ncol = 2))
  power_par[1:nrow(params),] <- params
  power_par[(nrow(params)+1):(nrow(params)+nrow(params1)),] <- params1[,3:4]
  
  colnames(power_par) <- c("a","b")
  colnames(data1) <- as.numeric(x)
  colnames(data_fit) <- as.numeric(x)
  data_fit <- apply(data_fit, 2, as.numeric)
  rownames(data_fit) <- rownames(data1)
  power_par <- apply(power_par, 2, as.numeric)
  rownames(power_par) <- rownames(data1)
  
  fit_result = list(original_data = X, trans_data = data1,
                    power_par = power_par, power_fit = data_fit,
                    Time =x)
}
result_all <- get_all_otus_result(otu_data_bec = a_result,params = params,params1 = params1)
# result_all <- get_all_otus_result(otu_data_bec = a_result,params = params,params1 = params1)

power_equation_plot <- function(result, label = 10,title,n=10,params,params1){
  #result = result_all;label = 10;title = "A fit result";n = 10;params = params;params1 = params1
  data1 = result[[2]]
  data2 = result[[4]]
  # no = sample(1:nrow(data1),n)
  # no = rownames(a4$original_data)
  # no = c("Carjivirus")
  df_original =  reshape2::melt(as.matrix(data1))
  df_fit = reshape2::melt(as.matrix(data2))
  
  label = 10
  p <- ggplot() +
    geom_point(df_original, mapping = aes_string(x = "Var2", y = "value"),colour = "#015493",
               show.legend = F, alpha = 0.5, shape = 1) +
    geom_line(df_fit, mapping = aes_string(x = "Var2", y = "value"),colour = "#015493", size = 1.25, show.legend = F)  +
    facet_wrap(~Var1) +
    xlab("Habitat Index") + ylab("Niche Index") + theme(axis.title=element_text(size=18)) +
    theme_bw() +
    geom_text(data = df_fit, aes(label = Var1, x = ((min(Var2)+max(Var2))/2), y = max(df_original$value) * 0.9), show.legend = FALSE, check_overlap = TRUE, size = 4)+
    theme_bw() +
    theme(axis.title=element_text(size=15),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10,hjust = 0),
          panel.spacing = unit(0.0, "lines"),
          plot.margin = unit(c(1,1,1,1), "lines"),
          strip.background = element_blank(),
          plot.background = element_blank(),
          strip.text = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    ggtitle(title)
  
  
  if (is.null(label)) {
    p = p
  } else {
    xlabel = ggplot_build(p)$layout$panel_params[[1]]$x.sec$breaks
    ylabel = ggplot_build(p)$layout$panel_params[[1]]$y.sec$breaks
    xlabel2 = parse(text= paste(label,"^", xlabel, sep="") )
    if (is.na(ylabel[1]) | ylabel[1] == 0) {
      ylabel2 = parse(text=c(0,paste(label,"^", ylabel[2:length(ylabel)], sep="")))
    } else{
      ylabel2 = parse(text= paste(label,"^", ylabel, sep="") )
    }
    p = p + scale_x_continuous(labels = xlabel2) + scale_y_continuous(labels = ylabel2)
  }
  return(p)
}

# p <- power_equation_plot(result_all, label = 10, title = "H fit result", params = params,params1 = params1)
p <- power_equation_plot(result_all, label = 10, title = "A fit result", params = params,params1 = params1)

# ggsave("./H_fit_result1.png", plot = p, width = 30, height = 20)
ggsave("./A_fit_result1.png", plot = p, width = 30, height = 20)


data_match <- function(result1, result2){
  matchname = intersect(rownames(result1$original_data),rownames(result2$original_data))

  new_result1 = list(original_data = result1$original_data[matchname,],
                     trans_data = result1$trans_data[matchname,],
                     power_par = result1$power_par[matchname,],
                     power_fit = result1$power_fit[matchname,],
                     Time = result1$Time)

  new_result2 = list(original_data = result2$original_data[matchname,],
                     trans_data = result2$trans_data[matchname,],
                     power_par = result2$power_par[matchname,],
                     power_fit = result2$power_fit[matchname,],
                     Time = result2$Time)
  result = list(dataset1 = new_result1, dataset2 = new_result2)
  return(result)
}
match_data = data_match(result1 = a_fit_result, result2 = h_fit_result)

bipower_equation_plot <- function(result, label = 10, n = 5, show.legend = TRUE,
                                  color1 = "#f57170", color2 = "#10ddc2", 
                                  title = "", seed = NULL, point_size = 1.5,
                                  line_size = 1.0, facet_cols = min(n, 5),
                                  plot_ratio = c(3, 2), panel_height = 1) {
  #result = match_data;label = 10;n = 5;show.legend = TRUE;color1 = "#f57170"; color2 = "#10ddc2"
  #title = ""; seed = NULL;point_size = 1.5;line_size = 1.0;facet_cols = min(n, 5);plot_ratio = c(3, 2);panel_height = 1
  # Set seed for reproducibility if provided
  if (!missing(seed)) {
    set.seed(seed)
  }
  
  # Load required libraries
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  require(gtable)
  require(grid)
  
  # Prepare data
  adenoma_original <- result$dataset1[[2]]
  adenoma_fitted <- result$dataset1[[4]]
  sample_indices <- sample(1:nrow(adenoma_original), n)

  
  # Process adenoma data
  df_adenoma_original <- melt(as.matrix(adenoma_original[sample_indices, ]))
  df_adenoma_original$type <- "Adenoma group"
  df_adenoma_fitted <- melt(as.matrix(adenoma_fitted[sample_indices, ]))
  df_adenoma_fitted$type <- "Adenoma group"
  df_adenoma_residuals <- df_adenoma_original
  df_adenoma_residuals$value <- df_adenoma_original$value - df_adenoma_fitted$value
  
  # Process control data
  control_original <- result$dataset2[[2]]
  control_fitted <- result$dataset2[[4]]
  df_control_original <- melt(as.matrix(control_original[sample_indices, ]))
  df_control_original$type <- "Control group"
  df_control_fitted <- melt(as.matrix(control_fitted[sample_indices, ]))
  df_control_fitted$type <- "Control group"
  df_control_residuals <- df_control_original
  df_control_residuals$value <- df_control_original$value - df_control_fitted$value
  
  # Combine datasets
  df_fit <- rbind(df_adenoma_fitted, df_control_fitted)
  df_original <- rbind(df_adenoma_original, df_control_original)
  df_residuals <- rbind(df_adenoma_residuals, df_control_residuals)
  
  # Convert sample IDs to character
  df_fit$Var1 <- as.character(df_fit$Var1)
  df_original$Var1 <- as.character(df_original$Var1)
  df_residuals$Var1 <- as.character(df_residuals$Var1)
  
  # Custom theme for consistent styling
  transparent_theme <- theme_minimal(base_size = 12) +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      strip.background = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, color = "black"),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 10, face = "bold"),
      panel.spacing = unit(1, "lines"),
      plot.margin = unit(c(5.5, 5.5, 2, 5.5), "points")  # Standard margins
    )
  
  # Create main plot
  main_plot <- ggplot() +
    geom_point(
      data = df_original,
      aes(x = Var2, y = value, color = type),
      alpha = 0.4, shape = 16, size = point_size,
      show.legend = show.legend
    ) +
    geom_line(
      data = df_fit,
      aes(x = Var2, y = value, color = type),
      size = line_size,
      show.legend = show.legend
    ) +
    facet_wrap(~ Var1, ncol = facet_cols) +
    labs(
      x = "",
      y = expression(paste("Niche Index")),
      color = "Study Group",
      linetype = "Study Group",
      title = title
    ) +
    scale_color_manual(
      values = c("Adenoma group" = color1, "Control group" = color2),
      labels = c("Adenoma", "Control")
    ) +
    scale_linetype_manual(
      values = c("Adenoma group" = "solid", "Control group" = "dashed"),
      labels = c("Adenoma", "Control")
    ) +
    transparent_theme +
    theme(
      legend.position = "right",
      panel.spacing.x = unit(12, "points")  # Consistent horizontal spacing
    )
  
  # Create residuals plot
  residuals_plot <- ggplot(df_residuals) +
    geom_point(
      aes(x = Var2, y = value, color = type),
      alpha = 0.4, shape = 16, size = point_size,
      show.legend = FALSE
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray30") +
    facet_wrap(~ Var1, ncol = facet_cols) +
    labs(
      x = expression(paste("Habitat Index (", log[10], ")")),
      y = "Residuals"
    ) +
    scale_color_manual(
      values = c("Adenoma group" = color1, "Control group" = color2)
    ) +
    transparent_theme +
    theme(
      strip.text = element_blank(),  # Remove facet labels for residuals
      panel.spacing.x = unit(12, "points")  # Match main plot spacing
    )
  
  # Format axis labels with exponential notation if label is provided
  if (!is.null(label)) {
    format_labels <- function(breaks) {
      if (label == 10) {
        parse(text = paste("10^", breaks, sep = ""))
      } else {
        parse(text = paste(label, "^", breaks, sep = ""))
      }
    }
    
    main_plot <- main_plot + 
      scale_x_continuous(labels = format_labels) +
      scale_y_continuous(labels = format_labels)
    
    residuals_plot <- residuals_plot +
      scale_x_continuous(labels = format_labels)
  }
  
  # Convert both plots to grobs
  g1 <- ggplotGrob(main_plot)
  g2 <- ggplotGrob(residuals_plot)
  
  # Get the max width from both plots
  max_width <- grid::unit.pmax(g1$widths, g2$widths)
  
  # Set equal widths for both plots
  g1$widths <- as.list(max_width)
  g2$widths <- as.list(max_width)
  
  # Set equal heights for panels
  panel_heights <- unit(rep(panel_height, n), "null")
  g1$heights[g1$layout$t[grep("panel", g1$layout$name)]] <- panel_heights
  g2$heights[g2$layout$t[grep("panel", g2$layout$name)]] <- panel_heights
  
  # Combine plots with precise alignment
  combined <- arrangeGrob(
    g1, g2,
    ncol = 1,
    heights = plot_ratio
  )
  
  # Return the combined plot
  return(combined)
}

plot_result <- bipower_equation_plot(
  match_data, 
  n = 4,
  facet_cols = 4,         # Number of columns
  panel_height = 1,     # Relative panel height
  plot_ratio = c(3, 2)    # Main to residual height ratio
)

ggsave("bipower_equation_plot.png", plot = plot_result, width = 12, height = 5,limitsize = FALSE)

