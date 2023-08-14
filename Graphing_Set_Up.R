library(ggplot2)
#set up the theme for all the graphs
theme_set(theme_bw()) 
theme_update(axis.text.x = element_text(size = 12),
             axis.text.y = element_text(size = 12),
             axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14),
             legend.title = element_text(size = 12),
             legend.text = element_text(size = 10),
             legend.spacing = unit(0.25, "cm"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             panel.spacing = unit(0, "lines"),
             legend.key = element_rect(fill = "white"),
             panel.spacing.y = unit(-0.25, "lines"),
             panel.border = element_rect(colour = "black", 
                                         fill = NA, size = 1),
             strip.text.x = element_text(size = 18, colour = "black", 
                                         face = "bold"))


colfunc <- colorRampPalette(rev(c("red","yellow", "grey40", "springgreen","royalblue")))
cols.d  <- data.frame(
  
  vals = seq(-1, 7.5, length = 50) , colo = colfunc(50)
)



elv <- cov %>% dplyr::select(elev, elev1)
