library(tidyverse)


xs <- seq(0, 10, 0.1)

weakS2_func <- function(x){
    ((x - 5)/10)^2 + 6
}
weakS1_func <- function(x){
    -(((x - 5)/10)^2) + 4
}
figure1_df <- data.frame('x' = xs,
                         'S1' = weakS1_func(xs),
                         'S2' = weakS2_func(xs)) %>%
    pivot_longer(-x) %>%
    mutate(Type="Weak")


strongS2_func <- function(x){
    ((x - 5)/4)^2 + 5.5 
}
strongS1_func <- function(x){
    -(((x - 5)/4)^2) + 4.5
}

figure2_df <- data.frame('x' = xs,
                         'S1' = strongS1_func(xs),
                         'S2' = strongS2_func(xs)) %>%
    pivot_longer(-x) %>%
    mutate(Type="Strong")

crossingS2_func <- function(x){
    -0.4*x + 7.1
}
crossingS1_func <- function(x){
    0.4*x + 2.9
}

figure3_df <- data.frame('x' = xs,
                         'S1' = crossingS1_func(xs),
                         'S2' = crossingS2_func(xs)) %>%
    pivot_longer(-x) %>%
    mutate(Type="Crossing")


figureTotal_df <- bind_rows(figure1_df, figure2_df, figure3_df) %>%
    mutate(Type = factor(Type, levels = c("Weak", "Strong", "Crossing")))

ggplot(data = figureTotal_df, aes(x = x, y = value, color = name)) +
    geom_line(size = 1.5) +
    scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0,10)) +
    labs(y = "Energy", x = "Time") +
    facet_wrap(~Type) +
    scale_color_manual(breaks = c("S1", "S2"), values = c("blue", "red")) +
    theme_bw() +
    theme(axis.text=element_blank(),
          strip.text = element_text(size=10),
          axis.title=element_text(size=15),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text.align = 0,
          legend.position = "top") 
ggsave("~/potentialparadox.github.io/Thesis/Images/naCrossings.png", width = 7.5, height = 4)
    
    
