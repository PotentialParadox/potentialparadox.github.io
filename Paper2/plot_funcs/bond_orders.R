library(tidyverse)
library(png)
library(grid)
library(ggimage)

setwd('~/Documents/paper2/bond_orders')

read_bondorder <- function(filename){
    read_csv(filename) %>%
        mutate(Bond=paste(Atom1, Atom2, sep='-')) %>%
        select(Trajectory, Timefs, Bond, BondOrder) %>%
        group_by(Timefs, Bond) %>%
        summarise(MeanBondOrder = mean(BondOrder)) %>%
        mutate(Bond = case_when(
            Bond == "14-15" ~ "d6",
            Bond == "15-16" ~ "d5",
            Bond == "16-17" ~ "d4",
            Bond == "6-7" ~ "d1",
            Bond == "7-8" ~ "d2",
            Bond == "8-9" ~ "d3",
            TRUE ~ "NA"
        ))
}

vacuum <- read_bondorder('bond-order-vacuum.csv') %>%
    mutate(NSolvent = "Vacuum")
ch3oh_5s <- read_bondorder('bond-order-ch3oh-5s.csv') %>%
    mutate(NSolvent = "5")
ch3oh_10s <- read_bondorder('bond-order-ch3oh-10s.csv') %>%
    mutate(NSolvent = "10")


bondorders = bind_rows(
    vacuum,
    ch3oh_5s,
    ch3oh_10s
)

img_df = data.frame("t" = c(500), "y" = c(1.2), "Bond" = c("d2"), "Image" = "ppvno2.png", "NSolvent" = c("Vacuum"))
ggplot() +
    geom_line(data = bondorders, aes(x = Timefs, y = MeanBondOrder, color = NSolvent), size=1.5) +
    geom_image(data=img_df, aes(x=t, y = y, image = Image), size = 1) +
    facet_wrap(~ Bond) +
    labs(x = "Time (fs)", y = "Bond Order", 
         color = "Number QM Solvents") +
    theme_bw() +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20),
          strip.text.x=element_text(size=20),
          legend.title = element_text(size=20),
          legend.text = element_text(size = 18),
          legend.text.align = 0,
          legend.position = c(0.85, 0.9)
    )
ggsave("~/potentialparadox.github.io/Paper2/Images/bond_order/solvent_comparison.png", width = 13, height = 10)
