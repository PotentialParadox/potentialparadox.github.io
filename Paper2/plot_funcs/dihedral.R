library(tidyverse)
setwd('~/Documents/paper2/dihedrals/')

collect_dihedrals <- function (runtime, infile, description) {
  data <- read_csv(infile, col_names = FALSE)
  nsteps <- ncol(data)
  time_labels <- (0:(nsteps-1))*runtime/(nsteps)
  names(data) <- as.character(time_labels)
  data %>% rowid_to_column("Trajectory") %>% 
    gather(key = "Time", value = "Angle", -Trajectory) %>% 
    group_by(Time) %>% 
    summarise(MeanAngle = mean(Angle)) %>% 
    mutate(Description = description) %>% 
    mutate(Time = as.numeric(Time))
}


####################
# PPV3-NO2
####################

# Sm - Vacuum from ~/backup4/TestRuns/Paper2/ppv-no2-validation/vacuum-bla
# Near angles are [[7,8,9,10],[16,15,14,13]
# Far angles are [[18,17,16,15], [16,15,14,13]]
# Only 260 samples
sm_ppv3_no2_vac_near <- collect_dihedrals(runtime = 1,
                                         infile = 'dihedral-ppv3-no2-sm-1.csv',
                                         description = "Sm-Near")

sm_ppv3_no2_vac_far <- collect_dihedrals(runtime = 1,
                                         infile = 'dihedral-ppv3-no2-sm-2.csv',
                                         description = "Sm-Far")

dihedral_data <- bind_rows(
  sm_ppv3_no2_vac_near,
  sm_ppv3_no2_vac_far
) %>% mutate(Description = as.factor(Description))

dihedral_data %>% ggplot(aes(x = Time, y = MeanAngle, color = Description)) +
  geom_point() +
  geom_smooth(method = 'loess', formula = y ~ x) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)
)
        
