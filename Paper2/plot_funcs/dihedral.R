library(tidyverse)
setwd('~/Documents/paper2/dihedrals/')

collect_dihedrals <- function (runtime, infile, state, location, solvent) {
  data <- read_csv(infile, col_names = FALSE)
  nsteps <- ncol(data)
  time_labels <- (0:(nsteps-1))*runtime/(nsteps-1)
  names(data) <- as.character(time_labels)
  data %>% rowid_to_column("Trajectory") %>%
    gather(key = "Time", value = "Angle", -Trajectory) %>%
    group_by(Time) %>%
    summarise(MeanAngle = mean(Angle)) %>%
    mutate(State = state) %>%
    mutate(Location = location) %>%
    mutate(Time = as.numeric(Time)) %>%
    mutate(SolventID = solvent) %>%
    mutate(Solvent = if_else(location == '', solvent, sprintf("%s-%s", solvent, location))) %>%
    mutate(Description = if_else(location == '', state, sprintf("%s-%s", state, location)))
}


#################################
# PPV3-NO2 SM Solvent Comparison
#################################

# Near states
# dihs = getDihedrals(args.n_trajs, suffix, [[7, 8, 9, 10],
                                          ## [5, 6, 7, 8]])
vacuum_near <- collect_dihedrals(runtime = 1,
                                infile = 'dihedral-vacuum-near.csv',
                                state = "Sm",
                                solvent = "Vacuum",
                                location = "Near")
ch3oh_near <- collect_dihedrals(runtime = 1,
                                infile = 'dihedral-ch3oh-near.csv',
                                state = "Sm",
                                solvent = "0",
                                location = "Near")
ch3oh_5s_near <- collect_dihedrals(runtime = 1,
                                   infile = 'dihedral-ch3oh-5s-near.csv',
                                   state = "Sm",
                                   solvent = "5",
                                   location = "Near")
ch3oh_10s_near <- collect_dihedrals(runtime = 1,
                                    infile = 'dihedral-ch3oh-10s-near.csv',
                                    state = "Sm",
                                    solvent = "10",
                                    location = "Near")
dihedral_near <- bind_rows(
    vacuum_near,
    ch3oh_near,
    ch3oh_5s_near,
    ch3oh_10s_near
)

# Far states
        ## dihs = getDihedrals(args.n_trajs, suffix, [[18, 17, 16, 15],
        ##                                            [16, 15, 14, 13]])
vacuum_far <- collect_dihedrals(runtime = 1,
                                infile = 'dihedral-vacuum-far.csv',
                                state = "Sm",
                                solvent = "Vacuum",
                                location = "Far")
ch3oh_far <- collect_dihedrals(runtime = 1,
                                infile = 'dihedral-ch3oh-far.csv',
                                state = "Sm",
                                solvent = "0",
                                location = "Far")
ch3oh_5s_far <- collect_dihedrals(runtime = 1,
                                   infile = 'dihedral-ch3oh-5s-far.csv',
                                   state = "Sm",
                                   solvent = "5",
                                   location = "Far")
ch3oh_10s_far <- collect_dihedrals(runtime = 1,
                                    infile = 'dihedral-ch3oh-10s-far.csv',
                                    state = "Sm",
                                    solvent = "10",
                                    location = "Far")
dihedral_far <- bind_rows(
    vacuum_far,
    ch3oh_far,
    ch3oh_5s_far,
    ch3oh_10s_far
)

dihedral_data <- bind_rows(
    dihedral_near,
    dihedral_far
) %>%
    mutate(Location = factor(Location, levels = c("Near", "Far"))) %>%
    mutate(SolventID = factor(SolventID, levels = c("Vacuum", "0", "5", "10")))

dihedral_data %>%
    ggplot(aes(x = Time, y = MeanAngle, color = SolventID)) +
    geom_line(size = 1.5) +
    theme_bw() +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    labs(x = "Time (ps)",
         y = expression(paste('Angle ', theta, ' (Degrees)')),
         color = "Number QM Solvents") +
    theme(axis.text=element_text(size=15),
          strip.text = element_text(size=12),
          axis.title=element_text(size=20),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.text.align = 0,
          legend.position = "top") +
    facet_wrap(~ Location)
ggsave("~/potentialparadox.github.io/Paper2/Images/dihedral/solvent_comparison.png", width = 7.5, height = 5)
