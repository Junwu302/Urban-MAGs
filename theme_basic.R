library(ggplot2)
theme_0 = theme_bw() + theme(axis.title = element_text(size=18, color = "black"),
                                      axis.text = element_text(size=16, color = "black"),
                                      plot.margin = unit(c(.5, 1, .5, .5), "cm"),
                                      legend.position = "none")

theme_1 = theme_bw() + theme(axis.title = element_text(size=18, color = "black"),
                             axis.text.x = element_text(angle = 45,size=14, hjust = 1,
                                                        vjust = 1, color = "black"),
                             axis.text.y = element_text(size=16, color = "black"),
                             plot.margin = unit(c(.5, 1, .5, .5), "cm"),
                             legend.position = "none")

theme_2 = theme_bw() + theme(axis.title = element_text(size=18, color = "black"),
                             axis.text = element_text(size=16, color = "black"),
                             plot.margin = unit(c(.5, 1, .5, .5), "cm"),
                             legend.text = element_text(size=16),
                             legend.title = element_text(size=18))

theme_3 = theme_bw() + theme(axis.title = element_text(size=18, color = "black"),
                             axis.text.x = element_text(angle = 45,size=14, hjust = 1,
                                                        vjust = 1, color = "black"),
                             axis.text.y = element_text(size=16, color = "black"),
                             plot.margin = unit(c(.5, 1, .5, .5), "cm"),
                             legend.text = element_text(size=16),
                             legend.title = element_text(size=18))
