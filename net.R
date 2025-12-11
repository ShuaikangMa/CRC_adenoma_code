rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

normalization <- function(x, z = 0.2){(x-min(x))/(max(x)-min(x))+z}


network_conversion <- function(result){
  #result = qdODE_parallel_result$ode_result$Anaerosalibacter
  n = ncol(result$fit)
  #对ode求解后的拟合值按列进行求平均值，提取最后的独立变量与依赖变量的变化率值的均值
  effect.mean = apply(result$fit,2,mean)[4:n]
  #对勒让德多项式拟合后的拟合值按列求平均值，提取最后的独立变量与依赖变量的变化率值的均值
  effect.predict.mean = apply(result$predict,2,mean)[4:n]
  #对ode求解后的拟合值按列进行求和，提取最后的独立变量与依赖变量的变化率值的总和值
  effect.total = colSums(result$fit)[4:n]
  
  #创建勒让德拟合互作列表，行数为独立变量与依赖变量的个数总和，列数为上述三个均值
  temp = matrix(NA,nrow = n-3, ncol=3)
  #设置互作的个体类别
  colnames(temp) = c("From", "To", "Effect")
  #设置第一列为互作列分别为独立变量与各个依赖变量
  temp[,1] = colnames(result$fit)[4:n]
  #设置第二列为独立变量名称
  temp[,2] = colnames(result$fit)[4]
  #设置第三列为勒让德拟合后独立变量与依赖变量均值
  temp[,3] = effect.predict.mean
  if (nrow(temp)==2) {
    temp = t(data.frame(temp[-1,]))
  } else{
    #删除第一行独立变量与独立变量之间互作
    temp = data.frame(temp[-1,])
  }
  #ode求解拟合互作列表，行数为独立变量与依赖变量的个数总和，列数为上述三个均值
  temp2 = matrix(NA,nrow = n-3, ncol=3)
  #设置互作的个体类别
  colnames(temp2) = c("From", "To", "Effect")
  #设置第一列为互作列分别为独立变量与各个依赖变量
  temp2[,1] = colnames(result$fit)[4:n]
  #设置第二列为独立变量名称
  temp2[,2] = colnames(result$fit)[4]
  #设置第三列为ode拟合后独立变量与依赖变量均值
  temp2[,3] = effect.total
  if (nrow(temp2)==2) {
    temp2 = t(data.frame(temp2[-1,]))
  } else{
    #删除第一行独立变量与独立变量之间互作
    temp2 = data.frame(temp2[-1,])
  }
  #              设置独立变量名称
  output <- list(ind.name = colnames(result$fit)[4],
                 #设置依赖变量名称
                 dep.name = colnames(result$fit)[5:n],
                 #设置optim(qdODEls)后最小的误差值为ODE.par
                 ODE.par = result$ODE.value,
                 #独立变量各阶勒让德线性拟合系数作为独立变量参数
                 ind.par = result$LOP_par[,3],
                 #各依赖变量各阶勒让德线性拟合系数作为依赖变量参数
                 dep.par = result$LOP_par[,4:(n-1)],
                 #将勒让德拟合后各变化率的均值作为效应均值
                 effect.mean = effect.predict.mean,
                 #将ode求解后的各变化率的均值作为效应总值
                 effect.total = effect.total,
                 #将所有的ode求解后的拟合值作为总体效应值
                 effect.all = result$fit,
                 #勒让德拟合后的变化率均值作为边
                 edge = temp,
                 #ode求解拟合后的变化率均值为边总和
                 edge.total = temp2,
                 #将勒让德拟合后的独立效应的均值做为独立效应值
                 ind.effect = effect.predict.mean[1])
  return(output)
}

net_all_a = lapply(qdODE_parallel_result_a$ode_result, network_conversion)
net_all_h = lapply(qdODE_parallel_result_h$ode_result, network_conversion)

net_m2_a = lapply(ode.M2_a$ode_result, network_conversion)
net_m2_h = lapply(ode.M2_h$ode_result, network_conversion)

net_m3_a = lapply(ode.M3_a$ode_result, network_conversion)
net_m3_h = lapply(ode.M3_h$ode_result, network_conversion)

net_m4_a = lapply(ode.M4_a$ode_result, network_conversion)
net_m4_h = lapply(ode.M4_h$ode_result, network_conversion)


net_m5_a = lapply(ode.M5_a$ode_result, network_conversion)
net_m5_h = lapply(ode.M5_h$ode_result, network_conversion)

net_m6_a = lapply(ode.M6_a$ode_result, network_conversion)
net_m6_h = lapply(ode.M6_h$ode_result, network_conversion)


net_m8_a = lapply(ode.M8_a$ode_result, network_conversion)
net_m8_h = lapply(ode.M8_h$ode_result, network_conversion)

net_m9_a = lapply(ode.M9_a$ode_result, network_conversion)
net_m9_h = lapply(ode.M9_h$ode_result, network_conversion)

net_m10_a = lapply(ode.M10_a$ode_result, network_conversion)
net_m10_h = lapply(ode.M10_h$ode_result, network_conversion)


net_m12_h = lapply(ode.M12_h$ode_result, network_conversion)
net_m12_a = lapply(ode.M12_a$ode_result, network_conversion)
 
net_m13_h = lapply(ode.M13_h$ode_result, network_conversion)
net_m13_a = lapply(ode.M13_a$ode_result, network_conversion)


net_m15_h = lapply(ode.M15_h$ode_result, network_conversion)
net_m15_a = lapply(ode.M15_a$ode_result, network_conversion)

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


library(igraph)
library(Cairo)
network_plot <- function(result, title = NULL, maxeffect = NULL, type = NULL,data_name,plot_type,point_size,title_size,select_type,plot_size,change_name= NULL,change_x= NULL,change_y= NULL){
  # result = net_m12_a;type = NULL;maxeffect = NULL
  #result = net_all_a;type = NULL;maxeffect = NULL
  #ind effect control node size
  #提取节点的信息
  extra <- sapply(result,"[[", "ind.effect")
  
  #提取边的信息，查看有无指定类型
  if (is.null(type)) {
    #勒让德拟合后的独立变量与依赖变量之间变化率均值作为边
    after <- data.frame(do.call(rbind, lapply(result, "[[", "edge")))
  } else{
    #ode拟合后的独立变量与依赖变量之间变化率均值作为边
    after <- data.frame(do.call(rbind, lapply(result, "[[", "edge.total")))
    
  }
  #将其数值部分转变为数值类型
  after$Effect = as.numeric(after$Effect)
  #将行名修改为空白
  rownames(after) = NULL

  #定义边的颜色，红色为促进，蓝色为抑制
  after$edge.colour = NA
  for (i in 1:nrow(after)) {
    if(after$Effect[i]>=0){
      after$edge.colour[i] = "#FE433C"
    } else{
      after$edge.colour[i] = "#0095EF"
    }
  }

  #创建关于节点信息的数据框
  #after[,2]提取边数据中的目标节点并获取唯一值分别设置为id名和名字
  #将各个物种的效应值作为第三列，该操作将依赖变量与独立变量的名称相对应故可以使用独立变量的值
  nodes <- data.frame(unique(after[,2]),unique(after[,2]),extra)
  colnames(nodes) <- c("id","name","ind_effect")
  
  #将 after 数据框按照 To 列中的每个唯一值进行分组，然后对每组的 Effect 值进行求和，并取其数值列
  #可以得到每个节点对其他节点的总的影响程度
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
  # 首先计算每个节点的出度和入度，并标记Effect的正负
  degree_data <- data %>%
    # 计算出度（From节点）并标记Effect正负
    group_by(Node = From) %>%
    summarise(
      Out_positive = sum(Effect > 0),  # 正Effect的出度
      Out_negative = sum(Effect < 0)   # 负Effect的出度
    ) %>%
    full_join(
      # 计算入度（To节点）并标记Effect正负
      data %>% 
        group_by(Node = To) %>% 
        summarise(
          In_positive = sum(Effect > 0),  # 正Effect的入度
          In_negative = sum(Effect < 0)   # 负Effect的入度
        ),
      by = "Node"
    ) %>%
    # 处理NA值
    mutate(
      Out_positive = replace_na(Out_positive, 0),
      Out_negative = replace_na(Out_negative, 0),
      In_positive = replace_na(In_positive, 0),
      In_negative = replace_na(In_negative, 0),
      Total = Out_positive + Out_negative + In_positive + In_negative
    ) %>%
    arrange(desc(Total)) %>%  # 按总度数降序排列
    mutate(Node = fct_inorder(Node))  # 锁定排序
  
  # 转换为绘图格式
  plot_data <- degree_data %>%
    select(Node, Out_positive, Out_negative, In_positive, In_negative) %>%
    mutate(
      Out_positive = -Out_positive,  # 出度转为负值
      Out_negative = -Out_negative   # 出度转为负值
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
      Label = ifelse(Count == 0, "", as.character(abs(Count))),  # 0值不显示标签
      LabelPos = ifelse(Count > 0, 
                        Count - 0.3 * sign(Count), 
                        Count + 0.3 * sign(Count))
    )
  
  # 生成y轴刻度hall-5
  y_breaks <- unique(c(seq(min(plot_data$Count)-24, 0, by = 1),
                       seq(0, max(plot_data$Count)+1, by = 1)))
  
  # 绘制高质量柱状图
  p<-ggplot(plot_data, aes(x = Node, y = Count, fill = interaction(DegreeType, EffectType))) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3, position = "stack") +  # 使用stack确保正负分开显示
    geom_hline(yintercept = 0, color = "gray30", linewidth = 0.5) +
    
    # 颜色设置（正红负蓝）
    scale_fill_manual(
      values = c(
        "Out.positive" = "#FE433C",  # 正Effect出度
        "Out.negative" = "#0095EF",  # 负Effect出度
        "In.positive" = "#FE433C",   # 正Effect入度
        "In.negative" = "#0095EF"    # 负Effect入度
      ),
      labels = c(
        "Out.positive" = "Out-degree (positive)",
        "Out.negative" = "Out-degree (negative)",
        "In.positive" = "In-degree (positive)",
        "In.negative" = "In-degree (negative)"
      ),
      name = "Degree Type"
    ) +
    
    # y轴设置
    scale_y_continuous(
      breaks = y_breaks,
      labels = abs,
      expand = expansion(mult = c(0.1, 0.1))  # 调整y轴范围
    ) +
    
    # 坐标轴和标题
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
      panel.background = element_blank(),        # 移除面板背景
      plot.background = element_blank(),        # 移除绘图区域背景
      panel.grid = element_blank(),             # 移除所有网格线
      panel.border = element_blank()            # 移除边框
    ) +
    
    # 添加边框
    theme(panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5))
  save_path = "C:/Users/Administrator/Desktop/"
  save_all = paste0(save_path, title,".png")
  ggsave(save_all, p, bg = "transparent", width = 16, height = 9, dpi = 300)
  
  # after_df = cbind(after, nodes[match(after$To, nodes$name), c("influence", "node.colour")])
  # write.csv(after, file = "C:/Users/Administrator/Desktop/m5_a_after.csv", row.names = F)
  # write.csv(nodes, file = "C:/Users/Administrator/Desktop/m5_a_nodes.csv", row.names = F)
  if (is.null(maxeffect)) {
    after[,3] <- normalization(abs(after[,3]))
    nodes[,3:4] <- normalization(abs(nodes[,3:4]))
  } else{
    after[,3] <- (abs(after[,3]))/maxeffect*1.5+0.1
    nodes[,3:4] <- (abs(nodes[,3:4]))/maxeffect*1.5+0.1
  }
  
  #将数据进行图像化转变
  net <- graph_from_data_frame( d=after,vertices = nodes,directed = T )
  #将图像化转换后的数据进行画图布局调整
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

  save_path="C:/Users/Administrator/Desktop/"
  plot_name=paste(save_path,data_name,"network_plot_", plot_type,select_type ,".png", sep = "")
  CairoPNG(plot_name, width=plot_size, height=plot_size, bg="white", res=300)
  plot.igraph(net,
              # edge.arrow.size = arrow_scaled,
              edge.arrow.size = E(net)$Effect*3,       # 箭头头部大小
              edge.arrow.width = E(net)$Effect*2,       # 箭头宽度 
              edge.arrow.mode = 2,          # 箭头方向 (0=无,1=后端,2=前端,3=双向)
              edge.arrow.curved = 0.1,      # 箭头弯曲度
              #节点标签使用节点名称
              vertex.label=V(net)$name,
              #节点标签颜色为黑色
              vertex.label.color="black",
              #节点形状为圆形
              vertex.shape="circle",
              #节点标签大小与 ind_effect 成比例
              vertex.label.cex=V(net)$ind_effect*title_size,
              #节点大小与 ind_effect 成比例
              vertex.size=V(net)$ind_effect*6+point_size,
              #边的弯曲程度为 0.05
              edge.curved=0.05,
              #边颜色根据 edge.colour 设置
              edge.color=E(net)$edge.colour,
              #边框颜色根据 edge.colour 设置
              edge.frame.color=E(net)$edge.colour,
              #边宽度与 Effect 成比例
              edge.width=E(net)$Effect*5,
              #节点颜色根据 node.colour 设置
              vertex.color=V(net)$node.colour,
              # 设置为NA去除节点周围的圆圈边
              vertex.frame.color='black',  
              #使用随机布局
              layout=l,
              #设置图形标题
              main=title,
              #设置图形边距
              margin=c(-.05,-.05,-.05,-.05))
  
  dev.off()
  
}
#fr kk circle randomly grid
network_plot(net_all_a, title = "Adenoma Module Network",data_name = "A",title_size = 4
             ,plot_type = "circle",point_size = 20,select_type = "penaltyall",plot_size = 5000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_all_h, title = "Healthy Module Network",data_name = "H",title_size = 4
             ,plot_type = "circle",point_size = 20,select_type = "penaltyall",plot_size = 5000
             ,change_name = NULL  ,change_x = 1,change_y = -1)

network_plot(net_m1, title = "M1 Network",data_name = "",title_size = 2
             ,plot_type = "circle",point_size = 7,select_type = "penaltym1",plot_size = 8000
             ,change_name = NULL  ,change_x = 1,change_y = -1)


network_plot(net_m2_a, title = "Adenoma M2 Network",data_name = "A",title_size = 3
             ,plot_type = "circle",point_size = 5,select_type = "penaltym2",plot_size = 8000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_m2_h, title = "Healthy M2 Network",data_name = "H",title_size = 3
             ,plot_type = "circle",point_size = 5,select_type = "penaltym2",plot_size = 8000
             ,change_name = NULL  ,change_x = 1,change_y = -1)




network_plot(net_m2, title = "M2 Network",data_name = "",title_size = 3
             ,plot_type = "circle",point_size = 1,select_type = "penaltym2",plot_size = 22000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_m3, title = "M3 Network",data_name = "",title_size = 2
             ,plot_type = "circle",point_size = 2,select_type = "penaltym3",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)


network_plot(net_m3_a, title = "Adenoma M3 Network",data_name = "A",title_size = 2
             ,plot_type = "grid",point_size = 14,select_type = "penaltym3",plot_size = 8000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_m3_h, title = "Healthy M3 Network",data_name = "H",title_size = 2
             ,plot_type = "grid",point_size = 14,select_type = "penaltym3",plot_size = 8000
             ,change_name = NULL  ,change_x = 1,change_y = -1)



network_plot(net_m4_a, title = "Adenoma M4 Network",data_name = "A1",title_size = 2
             ,plot_type = "circle",point_size = 14,select_type = "penaltym4",plot_size = 8000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_m4_h, title = "Healthy M4 Network",data_name = "H1",title_size = 2
             ,plot_type = "circle",point_size = 14,select_type = "penaltym4",plot_size = 8000
             ,change_name = NULL  ,change_x = 1,change_y = -1)


network_plot(net_m5_a, title = "Adenoma M5 Network",data_name = "A",title_size = 2
             ,plot_type = "circle",point_size = 5,select_type = "penaltym5",plot_size = 8000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_m5_h, title = "Healthy M5 Network",data_name = "H",title_size = 2
             ,plot_type = "circle",point_size = 5,select_type = "penaltym5",plot_size = 8000
             ,change_name = NULL  ,change_x = 1,change_y = -1)

network_plot(net_m6_a, title = "Adenoma M6 Network",data_name = "A",title_size = 2
             ,plot_type = "circle",point_size = 5,select_type = "penaltym6",plot_size = 8000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_m6_h, title = "Healthy M6 Network",data_name = "H",title_size = 2
             ,plot_type = "circle",point_size = 5,select_type = "penaltym6",plot_size = 8000
             ,change_name = NULL  ,change_x = 1,change_y = -1)


network_plot(net_m7, title = "M7 Network",data_name = "",title_size = 2
             ,plot_type = "circle",point_size = 14,select_type = "penaltym7",plot_size = 5000
             ,change_name = NULL  ,change_x = 1,change_y = -1)

network_plot(net_m8_a, title = "Adenoma M8 Network",data_name = "A",title_size = 2
             ,plot_type = "circle",point_size = 9,select_type = "penaltym8",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_m8_h, title = "Healthy M8 Network",data_name = "H",title_size = 2
             ,plot_type = "circle",point_size = 9,select_type = "penaltym8",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)

network_plot(net_m9_a, title = "Adenoma M9 Network",data_name = "A",title_size = 2
             ,plot_type = "circle",point_size = 9,select_type = "penaltym9",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_m9_h, title = "Healthy M9 Network",data_name = "H",title_size = 2
             ,plot_type = "circle",point_size = 9,select_type = "penaltym9",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)

network_plot(net_m10_a, title = "Adenoma M10 Network",data_name = "A",title_size = 2
             ,plot_type = "circle",point_size = 9,select_type = "penaltym10",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_m10_h, title = "Healthy M10 Network",data_name = "H",title_size = 2
             ,plot_type = "circle",point_size = 9,select_type = "penaltym10",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)


network_plot(net_m12_h, title = "Healthy M12 Network",data_name = "H",title_size = 2
             ,plot_type = "grid",point_size = 9,select_type = "penaltym12",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_m12_a, title = "Adenoma M12 Network",data_name = "A",title_size = 2
             ,plot_type = "grid",point_size = 9,select_type = "penaltym12",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)


network_plot(net_m13_h, title = "Healthy M12 Network",data_name = "H",title_size = 2
             ,plot_type = "grid",point_size = 9,select_type = "penaltym12",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_m13_a, title = "Adenoma M12 Network",data_name = "A",title_size = 2
             ,plot_type = "grid",point_size = 9,select_type = "penaltym12",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)

network_plot(net_m15_h, title = "Healthy M15 Network",data_name = "H",title_size = 2
             ,plot_type = "grid",point_size = 9,select_type = "penaltym12",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
network_plot(net_m15_a, title = "Adenoma M15 Network",data_name = "A",title_size = 2
             ,plot_type = "grid",point_size = 9,select_type = "penaltym12",plot_size = 10000
             ,change_name = NULL  ,change_x = 1,change_y = -1)
