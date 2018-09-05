
###################################
### ggplot tile heatmap wrapper ###
###################################

tile_heatmap_wrapper<-function(input_matrix, 
  output_file, 
  width=10, 
  height=4, 
  units="in",
  colour_clip=4, 
  cluster='both', 
  xlab='x', 
  ylab='y', 
  xtick_labels=NULL, 
  ytick_labels=NULL,
  colour_type='continuous', 
  colour_low='blue', 
  colour_high='red', 
  colour_mid='white', 
  colour_limits=NULL, 
  mono=F, 
  na_colour="grey50",
  xaxis_angle=330, 
  xaxis_hjust=0, 
  xaxis_vjust=NULL, 
  xaxis_size=5, 
  omit_xtext=F, 
  omit_xticks=F,
  yaxis_angle=NULL, 
  yaxis_hjust=NULL, 
  yaxis_vjust=NULL, 
  yaxis_size=NULL, 
  omit_ytext=F, 
  omit_yticks=F, 
  plot_title='', 
  input_matrix_text=NULL, 
  text_size=0.25, 
  highlight_regions=NULL, 
  x_breaks=waiver(), 
  y_breaks=waiver(), 
  plot = T){

  ### variables 
  # input_matrix: matrix of heatmap values (required)
  # output_file: plot output file path
  # width: plot width in "units"
  # height: plot height in "units"
  # units: plot size units ("in", "cm", or "mm")
  # colour_clip: maximum absolute value of colour scale
  # cluster: heirarchically cluster ("none", "row", "column", "none")
  # xlab: x-axis label
  # ylab: y-axis label
  # xtick_labels: display labels for x ticks
  # ytick_labels: display labels for y ticks
  # colour_type: colour scale type ("continuous", "categorical")
  # colour_low: colour scale lower limit colour -- passed to scale_colour_gradient2
  # colour_high: colour scale upper limit colour -- passed to scale_colour_gradient2 
  # colour_mid: colour scale zero colour -- passed to scale_colour_gradient2
  # colour_limits: upper and lower value limits of colour scale -- passed to scale_colour_gradient2
  # mono: use monotype font (True, False)
  # na_colour: colour to use for NA values
  # xaxis_angle: rotation angle for x tick labels -- passed to element_text
  # xaxis_hjust: horizontal justification of x tick labels (in [0, 1]) -- passed to element_text
  # xaxis_vjust: vertical justification of x tick labels (in [0, 1]) -- passed to element_text
  # xaxis_size: text size of x tick labels (in pts) -- passed to element_text
  # omit_xtext: omit x tick labels (True, False)
  # omit_xticks: omit x ticks (True, False)
  # yaxis_angle: rotation angle for y tick labels -- passed to element_text
  # yaxis_hjust: horizontal justification of y tick labels (in [0, 1]) -- passed to element_text
  # yaxis_vjust: vertical justification of y tick labels (in [0, 1]) -- passed to element_text
  # yaxis_size: text size of y tick labels (in pts) -- passed to element_text
  # omit_ytext: omit y tick labels (True, False)
  # omit_yticks: omit y ticks (True, False)
  # plot_title: main title for plot
  # input_matrix_text: matrix of heatmap text
  # text_size: size of heatmap text
  # highlight_regions: list of highlighted regions of form: list("red" = list("region1" = c(_min_, _max_), "region2" = c(_min_, _max_), ...), "blue" = list("region3" = c(_min_, _max_), ...), ...)
  # x_breaks: x-axis breaks (for displaing xtick_labels)
  # y_breaks: y-axis breaks (for displaing ytick_labels)
  # plot: whether to plot the heatmap (True, False)

  require(ggplot2)

  order_row<-rev(1:dim(input_matrix)[1])
  order_col<-1:dim(input_matrix)[2]
  if(cluster %in% c('both', 'row')){
    d <- dist(input_matrix, method = "euclidean") # distance matrix
    order_row <- hclust(d, method="ward")$order
  }
  if(cluster %in% c('both', 'column')){
    d <- dist(t(input_matrix), method = "euclidean") # distance matrix
    order_col<- hclust(d, method="ward")$order    
  }
  plot_df<-melt(input_matrix[order_row,order_col])
  colnames(plot_df)<-c('y', 'x', 'value')
  plot_df$label<-""
  if(!is.null(input_matrix_text)){
    plot_df_text<-melt(input_matrix_text[order_row,order_col])
    colnames(plot_df_text)<-c('y', 'x', 'label')
    plot_df$label<-plot_df_text$label
  }
  if(colour_type=='continuous' & colour_clip){
    plot_df$value[plot_df$value>colour_clip]<-colour_clip
    plot_df$value[plot_df$value<(-colour_clip)]<-(-colour_clip)    
  }
  p <- ggplot(plot_df, aes(x, y)) + geom_tile(aes(fill = value)) + geom_text(aes(label = label), size=text_size) +
    # theme_bw() + 
    theme(axis.text.x=list(element_text(angle = xaxis_angle, hjust = xaxis_hjust, vjust = xaxis_vjust, size = xaxis_size, family=c('', 'mono')[as.numeric(mono)+1]), element_blank())[[as.numeric(omit_xtext)+1]],
          axis.text.y=list(element_text(angle = yaxis_angle, hjust = yaxis_hjust, vjust = yaxis_vjust, size = yaxis_size, family=c('', 'mono')[as.numeric(mono)+1]), element_blank())[[as.numeric(omit_ytext)+1]],
          axis.ticks.x=list(element_line(), element_blank())[[as.numeric(omit_xticks)+1]],
          axis.ticks.y=list(element_line(), element_blank())[[as.numeric(omit_yticks)+1]]) + 
    xlab(xlab) + ylab(ylab) + labs(title = plot_title)
  if(!is.null(highlight_regions)){
    for(i in names(highlight_regions)){
      for(j in names(highlight_regions[[i]])){
        p <- p + geom_rect(data = NULL, mapping = aes_now(xmin=highlight_regions[[i]][[j]][1]-0.5, xmax=highlight_regions[[i]][[j]][2]+0.5, ymin=highlight_regions[[i]][[j]][1]-0.5, ymax=highlight_regions[[i]][[j]][2]+0.5), fill = NA, colour = i)
      }
    }
  }
  #xtick labels specified
  if(!is.null(xtick_labels)){
    if(is.numeric(plot_df$x)){
      #X is numeric
      p <- p + scale_x_continuous(breaks=x_breaks, labels=xtick_labels)
    }else{
      #X is discrete
      p <- p + scale_x_discrete(breaks=x_breaks, labels=xtick_labels)
    }
  }else{
    if(is.numeric(plot_df$x)){
      p <- p + scale_x_continuous(breaks=x_breaks)
    }
  }
  #ytick labels specified
  if(!is.null(ytick_labels)){
    if(is.numeric(plot_df$y)){
      #Y is numeric
      p <- p + scale_y_continuous(breaks=y_breaks, labels = ytick_labels)
    }else{
      #Y is discrete
      p <- p + scale_y_discrete(breaks=y_breaks, labels=ytick_labels)
    }
  }else{
    if(is.numeric(plot_df$y)){
      p <- p + scale_y_continuous(breaks=y_breaks)
    }
  }
  if(colour_type=='continuous'){
    p <- p + scale_fill_gradient2(low = colour_low, high = colour_high, mid = colour_mid, midpoint = 0, limits=colour_limits, na.value=na_colour)    
  }
  if(colour_type=='categorical'){
    p <- p + scale_fill_brewer(palette='Set1')    
  }
  if(plot){
    ggsave(file=output_file, width=width, units=units, height=height)
  }
  return(p)
}