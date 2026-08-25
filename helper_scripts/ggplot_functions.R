library(ggplot2)

#' Customize the theme of a ggplot object
#'
#' @param plot A ggplot object
#' @param base_theme The base theme to use for the plot. Default is theme_bw().
#' @param title Title for the plot
#' @param plot_title_size Size of the plot title in pts. Default is 24.
#' @param xlab The x-axis title
#' @param ylab The y-axis title
#' @param axis_title_size Size of the axes titles in pts. Default is 20.
#' @param axis_text_size Size of the axes text in pts. Default is 18.
#' @param x_labels Optional labels to use for the x-axis
#' @param format_x_labels Logical value indicating whether to format the X axis labels by replacing "_" with a space and using title case. Default is FALSE.
#' @param x_labels_angle Number of degrees to rotate x-axis labels. Default is 0. 
#' @param show_legend A logical value indicating whether to show the legend or not. Default is TRUE. 
#' @param fill_title The title for fill in the legend
#' @param color_title The title for color in the legend
#' @param colors An optional vector of colors to use for the color mapping. 
#' @param fill_colors An optional vector of colors to use for the fill. 
#' @param fill_labels An optional vector of labels for the fill legend
#' @param legend_title_size Size of the legend titles in pts. Default is 20.
#' @param legend_text_size Size of the legend text in pts. Default is 18.
#' @param legend_key_size Size of the legend key in cm. Default is 1.
#' @param format_fill_labels Logical value indicating whether to format the legend labels by replacing "_" with a space and using title case. Default is FALSE.
#' @param scale_x A scale to use for the x-axis. If none provided, one is chosen by default.
#' @param scale_y A scale to use for the y-axis. If none provided, one is chosen by default. 
#' @param facet A variable to use for faceting. Default is not to facet. 
#' @param facet_nrow Number of rows to use for faceting. 
#' @param facet_ncol Number of rows to use for faceting. 
#' @param facet_scales scale for facet_wrap(). Default is "free". 
#' @param strip_text_size Size of the facet labels. Default is 20.
#' @param facet_labels Labels to use for panels when facetting
#' @return A ggplot object
#' @export
customize_ggplot_theme = function(plot, base_theme = theme_bw(), title = NULL, plot_title_size = 24,
  xlab = NULL, ylab = NULL, axis_title_size = 20, axis_text_size = 18, x_labels = NULL, format_x_labels = F, x_labels_angle = 0,  
  show_legend = T, fill_title = NULL, color_title = NULL, fill_colors = NULL, colors = NULL, fill_labels = waiver(), 
  legend_title_size = 20, legend_text_size = 18, legend_key_size = 1, format_fill_labels = F, 
  scale_x = NULL, scale_y = NULL, 
  facet = NULL, facet_nrow = NULL, facet_ncol = NULL, facet_scales = "free", strip_text_size = 20, facet_labels = NULL){
  
  plot = plot +
    base_theme +
    theme(plot.title = element_text(hjust = 0.5, size = plot_title_size), 
	    axis.title = element_text(size = axis_title_size), axis.text = element_text(size = axis_text_size), 
  	  legend.title = element_text(size = legend_title_size), legend.text = element_text(size = legend_text_size), 
      legend.key.size = unit(legend_key_size, "cm"), 
      strip.text = element_text(size = strip_text_size)) +
    labs(title = title, x = xlab, y = ylab, fill = fill_title,  color = color_title)
  
  # Remove legend if specified
  if(!show_legend){
    plot = plot + theme(legend.position = "None")
  }
  
  # Rotate x-axis lables if specified
  if(x_labels_angle != 0){
    plot = plot + theme(axis.text.x = element_text(angle = x_labels_angle, hjust = 1))
  }
  
  # Use x_labels if provided
  if(!is.null(x_labels)){
    plot = plot + scale_x_discrete(labels = x_labels)
  }
  
  # Format x axis labels if specified
  if(format_x_labels){
    xlabels = ggplot_build(plot)$layout$panel_params[[1]]$x$get_labels()
    new_xlabels = stringr::str_to_title(gsub("_", " ", xlabels))
    plot = plot + scale_x_discrete(labels = new_xlabels)
  }
  
  # Find the type of geom the plot uses
  plot_geom = class(plot$layers[[1]]$geom)
  
  # Change fill colours and labels if specified
  if(!is.null(fill_colors)){
    plot = plot + scale_fill_manual(values = fill_colors, labels = fill_labels)
  } else {
    plot = plot + scale_fill_discrete(labels = fill_labels)
  }
  
  # Change colors if specified
  if(!is.null(colors)){
    plot = plot + scale_color_manual(values = colors)
  }
  
  # If scale_x and scale_y not provided. they are inferred
  if(is.null(scale_x)){
    if("GeomBar" %in% plot_geom & !"GeomCol" %in% plot_geom){
      scale_x = scale_x_continuous(expand = c(0,0), labels = scales::comma) 
    }
  }
  
  if(is.null(scale_y)){
    if("GeomBar" %in% plot_geom){
      scale_y = scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::comma)
    }
  }
  
  plot = plot + scale_x + scale_y
  
  if(!is.null(facet)){
    if(!is.null(facet_labels)){
      original_facet_labels = sort(unique(plot$data[[facet]]))
      labeller = as_labeller(setNames(facet_labels, original_facet_labels))
    } else {
      labeller = "label_value"
    }
    plot = plot + facet_wrap(facet, nrow = facet_nrow, ncol = facet_ncol, labeller = labeller, scales = facet_scales)
  }
  
  return(plot)
  
}

#' Encode p-values with symbols
#'
#' @param p_values A vector of p-values
#' @param cutpoints A vector of cutpoints for significance. Default is c(0, 0.001, 0.01, 0.05, 1)
#' @param symbol
#' @return A vector of symbols representing p-values
#' @export
sig_sym = function(p_values, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbol = "*"){
  
  symbols = sapply(3:0, function(x) paste(rep(symbol, x), collapse = ""))
  
  significance_symbols = as.character(symnum(p_values, corr = FALSE, na = FALSE, cutpoints = cutpoints, symbols = symbols))
  return(significance_symbols)
}
