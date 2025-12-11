rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

normalization <- function(x, z = 0.2){(x-min(x))/(max(x)-min(x))+z}

network_conversion <- function(result){
  n = ncol(result$fit)
  effect.mean = apply(result$fit,2,mean)[4:n]
  effect.predict.mean = apply(result$predict,2,mean)[4:n]
  effect.total = colSums(result$fit)[4:n]
  
  temp = matrix(NA,nrow = n-3, ncol=3)
  colnames(temp) = c("From", "To", "Effect")
  temp[,1] = colnames(result$fit)[4:n]
  temp[,2] = colnames(result$fit)[4]
  temp[,3] = effect.predict.mean
  if (nrow(temp)==2) {
    temp = t(data.frame(temp[-1,]))
  } else{
    temp = data.frame(temp[-1,])
  }
  temp2 = matrix(NA,nrow = n-3, ncol=3)
  colnames(temp2) = c("From", "To", "Effect")
  temp2[,1] = colnames(result$fit)[4:n]
  temp2[,2] = colnames(result$fit)[4]
  temp2[,3] = effect.total
  if (nrow(temp2)==2) {
    temp2 = t(data.frame(temp2[-1,]))
  } else{
    temp2 = data.frame(temp2[-1,])
  }

  output <- list(ind.name = colnames(result$fit)[4],
                 dep.name = colnames(result$fit)[5:n],
                 ODE.par = result$ODE.value,
                 ind.par = result$LOP_par[,3],
                 dep.par = result$LOP_par[,4:(n-1)],
                 effect.mean = effect.predict.mean,
                 effect.total = effect.total,
                 effect.all = result$fit,
                 edge = temp,
                 edge.total = temp2,
                 ind.effect = effect.predict.mean[1])
  return(output)
}

network_maxeffect <- function(result){
  module = result[[1]]
  after <- data.frame(do.call(rbind, lapply(module, "[[", "edge.total")))
  after$Effect = as.numeric(after$Effect)
  Module.maxeffect = max(abs(after$Effect))
  
  if (length(result) == 2) {
    after <- data.frame(do.call(rbind, lapply(result[[2]], "[[", "edge")))
    after$Effect = as.numeric(after$Effect)
    res.maxeffect = max(abs(after$Effect))
  } else{
    after = lapply(2:length(result),function(c) data.frame(do.call(rbind, lapply(result[[c]], "[[", "edge"))))
    after = do.call(rbind,after)
    after$Effect = as.numeric(after$Effect)
    res.maxeffect = max(abs(after$Effect))
  }
  maxeffect = max(c(Module.maxeffect,res.maxeffect))
  return(maxeffect)
}

network_plot <- function(result, title = NULL, maxeffect = NULL, type = NULL,data_name,plot_type,point_size,title_size,select_type,plot_size,change_name= NULL,change_x= NULL,change_y= NULL){
  
  extra <- sapply(result,"[[", "ind.effect")
  
  if (is.null(type)) {
    after <- data.frame(do.call(rbind, lapply(result, "[[", "edge")))
  } else{
    after <- data.frame(do.call(rbind, lapply(result, "[[", "edge.total")))
    
  }

  after$Effect = as.numeric(after$Effect)
  rownames(after) = NULL

  after$edge.colour = NA
  for (i in 1:nrow(after)) {
    if(after$Effect[i]>=0){
      after$edge.colour[i] = "#FE433C"
    } else{
      after$edge.colour[i] = "#0095EF"
    }
  }

  nodes <- data.frame(unique(after[,2]),unique(after[,2]),extra)
  colnames(nodes) <- c("id","name","ind_effect")
  
  nodes$influence <- aggregate(Effect ~ To, data = after, sum)[,2]
  nodes$node.colour = NA


  for (i in 1:nrow(nodes)) {
    if(nodes$influence[i]>=0){
      nodes$node.colour[i] = "#F39B7FFF"
    } else{
      nodes$node.colour[i] = "#8491B4FF"
    }
  }
  data = after
  degree_data <- data %>%
    group_by(Node = From) %>%
    summarise(
      Out_positive = sum(Effect > 0),
      Out_negative = sum(Effect < 0)
    ) %>%
    full_join(
      data %>% 
        group_by(Node = To) %>% 
        summarise(
          In_positive = sum(Effect > 0),  # 正Effect的入度
          In_negative = sum(Effect < 0)   # 负Effect的入度
        ),
      by = "Node"
    ) %>%
    mutate(
      Out_positive = replace_na(Out_positive, 0),
      Out_negative = replace_na(Out_negative, 0),
      In_positive = replace_na(In_positive, 0),
      In_negative = replace_na(In_negative, 0),
      Total = Out_positive + Out_negative + In_positive + In_negative
    ) %>%
    arrange(desc(Total)) %>% 
    mutate(Node = fct_inorder(Node))

  plot_data <- degree_data %>%
    select(Node, Out_positive, Out_negative, In_positive, In_negative) %>%
    mutate(
      Out_positive = -Out_positive, 
      Out_negative = -Out_negative  
    ) %>%
    pivot_longer(
      cols = c("Out_positive", "Out_negative", "In_positive", "In_negative"),
      names_to = c("DegreeType", "EffectType"),
      names_sep = "_",
      values_to = "Count"
    ) %>%
    mutate(
      DegreeType = factor(DegreeType, levels = c("In", "Out")),
      EffectType = factor(EffectType, levels = c("positive", "negative")),
      Label = ifelse(Count == 0, "", as.character(abs(Count))),
      LabelPos = ifelse(Count > 0, 
                        Count - 0.3 * sign(Count), 
                        Count + 0.3 * sign(Count))
    )
  
  y_breaks <- unique(c(seq(min(plot_data$Count)-24, 0, by = 1),
                       seq(0, max(plot_data$Count)+1, by = 1)))
  
  p<-ggplot(plot_data, aes(x = Node, y = Count, fill = interaction(DegreeType, EffectType))) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3, position = "stack") +
    geom_hline(yintercept = 0, color = "gray30", linewidth = 0.5) +
    
    scale_fill_manual(
      values = c(
        "Out.positive" = "#FE433C",  
        "Out.negative" = "#0095EF",  
        "In.positive" = "#FE433C",   
        "In.negative" = "#0095EF"   
      ),
      labels = c(
        "Out.positive" = "Out-degree (positive)",
        "Out.negative" = "Out-degree (negative)",
        "In.positive" = "In-degree (positive)",
        "In.negative" = "In-degree (negative)"
      ),
      name = "Degree Type"
    ) +
    
    scale_y_continuous(
      breaks = y_breaks,
      labels = abs,
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    
    labs(
      x = "Module",
      y = "Degree Count"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      text = element_text(family = "Arial"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15, 
                                 margin = ggplot2::margin(t = 5)),
      axis.title = element_text(size = 15, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, 
                                   margin = ggplot2::margin(b = 15)),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      # plot.background = element_rect(fill = "white", color = NA)
      panel.background = element_blank(),        
      plot.background = element_blank(),     
      panel.grid = element_blank(),          
      panel.border = element_blank()     
    ) +
    
    theme(panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5))
  save_path = "./"
  save_all = paste0(save_path, title,".png")
  ggsave(save_all, p, bg = "transparent", width = 16, height = 9, dpi = 300)
  
  if (is.null(maxeffect)) {
    after[,3] <- normalization(abs(after[,3]))
    nodes[,3:4] <- normalization(abs(nodes[,3:4]))
  } else{
    after[,3] <- (abs(after[,3]))/maxeffect*1.5+0.1
    nodes[,3:4] <- (abs(nodes[,3:4]))/maxeffect*1.5+0.1
  }
  
  net <- graph_from_data_frame( d=after,vertices = nodes,directed = T )
  if (plot_type == "fr"){
    l <- layout_with_fr(net)
  }
  if (plot_type == "kk"){
    l <- layout_with_kk(net)
  }
  if (plot_type == "circle"){
    l <- layout_in_circle(net)
  }
  if (plot_type == "randomly"){
    l <- layout_randomly(net)
  }
  if (plot_type == "grid"){
    l <- layout_on_grid(net)
  }
  if (plot_type == "drl"){
    l <- layout_with_drl(net)
  }
  row.names(l) = names(result)
  if(!is.null(change_name)){
    if(length(change_name) == 1){
      l[which(row.names(l)==change_name),] = c(l[which(row.names(l)==change_name),][1]+change_x,l[which(row.names(l)==change_name),][2]+change_y)
    }
    if(length(change_name) == 2){
      l[which(row.names(l)==change_name[1]),] = c(l[which(row.names(l)==change_name[1]),][1]+change_x,l[which(row.names(l)==change_name[1]),][2]+change_y)
      l[which(row.names(l)==change_name[2]),] = c(l[which(row.names(l)==change_name[2]),][1]-change_x,l[which(row.names(l)==change_name[2]),][2]+change_y)
    }
  }

  save_path="./"
  plot_name=paste(save_path,data_name,"network_plot_", plot_type,select_type ,".png", sep = "")
  CairoPNG(plot_name, width=plot_size, height=plot_size, bg="white", res=300)
  plot.igraph(net,
              # edge.arrow.size = arrow_scaled,
              edge.arrow.size = E(net)$Effect*3, 
              edge.arrow.width = E(net)$Effect*2, 
              edge.arrow.mode = 2,       
              edge.arrow.curved = 0.1,  
              vertex.label=V(net)$name,
              vertex.label.color="black",
              vertex.shape="circle",
              vertex.label.cex=V(net)$ind_effect*title_size,
              vertex.size=V(net)$ind_effect*6+point_size,
              edge.curved=0.05,
              edge.color=E(net)$edge.colour,
              edge.frame.color=E(net)$edge.colour,
              edge.width=E(net)$Effect*5,
              vertex.color=V(net)$node.colour,
              vertex.frame.color='black',  
              layout=l,
              main=title,
              margin=c(-.05,-.05,-.05,-.05))
  dev.off()
  
}
