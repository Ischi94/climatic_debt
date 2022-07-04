# general configurations for plots

#set ggplot output
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 10))

ggplot2::theme_update(text = element_text(colour = "grey20", size = 10), 
                      legend.text = element_text(colour = "grey20", size = 10),
                      legend.title = element_text(colour = "grey20", size = 10),
                      axis.text = element_text(colour = "grey20", size = 10),
                      axis.ticks = element_line(colour = "grey50"),
                      strip.text = element_text(colour = "grey20", size = 10),
                      panel.grid.minor = element_blank(),  
                      panel.grid.major = element_blank(),
                      plot.title = element_text(colour = "grey20", size = 10))

# define output sizes
image_width <- 183
image_height <- 100
image_units <- "mm"

# define pallets

# define common color
colour_lavender = "#9678B6"
colour_brown = "#9C7673"
colour_green = "#105157"
colour_grey = "#D3CFC2"
colour_coral = "coral2"
