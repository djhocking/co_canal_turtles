##### Multi-Session SCR Model Creation JAGS ####
##

## Also dubbed stratified population model

#### Adding in Variation in N per site ####
library(dplyr)
library(tidyr)
library(rgdal)
library(sp)

testing <- TRUE
run_date <- Sys.Date()

# number of possible individuals per site
M <- 7000 
if(testing) {
  M <- 700
}

min_cap_rate <- 0.015

site_num_spatial <- as.matrix(cbind(c(2,4,6,7,1,9,8,3,5,10,11,12), 
c("A","C","D","E","F","G","J","K","L","M","N","O")))

Sites <- read.csv(file = "Data/trapids_sites.csv", header = TRUE, stringsAsFactors = FALSE)
Sites$site_num <- as.integer(site_num_spatial[ , 1])
coords <- read.csv(file = "Data/coords.csv", stringsAsFactors = FALSE)
EDF <- read.csv(file = "Data/EDF.csv", stringsAsFactors = FALSE)
n_traps_site <- read.csv(file = "Data/Max_Traps_Site.csv", stringsAsFactors = FALSE) # number of traps per site
n_traps_site$site_num <- as.integer(site_num_spatial[ , 1])
n_traps_site <- n_traps_site[order(n_traps_site$site_num), ]
n_traps <- n_traps_site$max_traps
# K <- max(EDF$day)
n_days <- max(EDF$day)

########## Process Trap Locations ###########
trap_locs_degrees <- coords
trap_locs_degrees$trap <- 1:nrow(trap_locs_degrees)
trap_num <- trap_locs_degrees$trap


# convert to utm to have distance in meters
coords_dd = SpatialPoints(coords[ , c("lon", "lat")], proj4string=CRS("+proj=longlat"))
coords_utm <- spTransform(coords_dd, CRS("+init=epsg:26917"))

trap_locs <- coords_utm
trap_locs <- as.data.frame(trap_locs)
trap_locs <- cbind(trap_num, trap_locs)
colnames(trap_locs) = c("trap_id", "easting", "northing")
# trap_locs as single vector with distance between

## Creating trap location vector per site using coordinates (sp package required)
# order of sites - (2,4,6,7,1,9,8,3,5,10,11,12)
trap_dist_2 <- spDistsN1(coords_utm[1:8, ], coords_utm[1, ])
trap_dist_4 <- spDistsN1(coords_utm[9:18, ], coords_utm[9, ])
trap_dist_6 <- spDistsN1(coords_utm[19:26, ], coords_utm[19, ])
trap_dist_7 <- spDistsN1(coords_utm[27:40, ], coords_utm[27, ])
trap_dist_1 <- spDistsN1(coords_utm[41:47, ], coords_utm[41, ])
trap_dist_9 <- spDistsN1(coords_utm[48:54, ], coords_utm[48, ])
trap_dist_8 <- spDistsN1(coords_utm[61:70, ], coords_utm[61, ])
trap_dist_3 <- spDistsN1(coords_utm[71:80, ], coords_utm[71, ])
trap_dist_5 <- spDistsN1(coords_utm[81:90, ], coords_utm[81, ])
trap_dist_10 <- spDistsN1(coords_utm[91:102, ], coords_utm[91, ])
trap_dist_11 <- spDistsN1(coords_utm[103:112, ], coords_utm[103, ])
trap_dist_12 <- spDistsN1(coords_utm[113:122, ], coords_utm[113, ])

trap_dist_list <- list(trap_dist_1, trap_dist_2, trap_dist_3, trap_dist_4,
                       trap_dist_5, trap_dist_6, trap_dist_7, trap_dist_8,
                       trap_dist_9, trap_dist_10, trap_dist_11, trap_dist_12)

trap_locs <- matrix(NA, 12, max(n_traps_site$max_traps))
for (i in 1:12) {
  trap_locs[i, 1:n_traps_site$max_traps[i]] <- trap_dist_list[[i]] / 100
}

xlim <- matrix(NA, 12, 2)
for(i in 1:12){
  xlim[i, 1:2] <- c(min(trap_dist_list[[i]]) - 250, max(trap_dist_list[[i]]) + 250) / 100 # need to have buffer on each side without being negative. 
}

####### EDF FILE ########

# only a small number of CPIC were not sexed (~1%). They technically could be used but it would be really difficult with NIMBLE and the two other zero tricks and site looks that are already being applied. For this analysis, we will exclude those individuals from the data. Same with NA of which there are only two (SODO).

# filter(EDF, is.na(sex))
# filter(EDF, sex == "U")

# only do analysis for Males and Females. Skip for unrecorded and juvenile because insufficient data for calculating separate home range sizes and such for juveniles.
EDF <- EDF %>%
  filter(sex != "U",
         !is.na(sex))

# Take out sites H and I and get just CPIC data
EDF_CPIC <- EDF %>%
  filter(site != "H" & site != "I" & species == "CPIC")

old <- c("A", "C", "D", "E", "F", "G", "J", "K", "L", "M", "N", "O")
new <- c(2,4,6,7,1,9,8,3,5,10,11,12)

EDF_CPIC$site <- as.integer(as.character(factor(EDF_CPIC$site, old, new)))


## subtract 6 from trap ids > = 61 (Sites H and I)
# EDF_CPIC$trap_id_edited <- ifelse(EDF_CPIC$trap_id >= 61, EDF_CPIC$trap_id - 6, EDF_CPIC$trap_id - 0)

##### Want EM ARRAY with ijk with index for site ########
EM_CPIC <- EDF_CPIC %>%
  group_by(site, ind, trap, day, sex) %>%
  select(site, ind, trap, day, sex) %>%
  mutate(count = 1) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  group_by(site) %>%
  mutate(id_site = as.integer(as.factor(ind))) %>%
  ungroup() %>%
  mutate(site_num = as.integer(as.factor(site)))
         # ,
         # site_id_ck = paste0(site, "_", ind))

old <- 1:12
new <- c(2,4,6,7,1,9,8,3,5,10,11,12)

EM_CPIC$site_num <- as.integer(as.character(factor(EM_CPIC$site_num, old, new)))


# get number of unique individuals caught per site, with spatially relevant site IDs
# inds <- EDF_CPIC %>%
#   group_by(site) %>%
#   select(site, ind) %>%
#   distinct() %>%
#   summarise(individuals = n())

# y[i,j,k,l] individual x trap x occassion x site
# to be the same size for an array will need to augment at least up to that while building array. Maybe could use nimbleList(). Not sure it's worth it since augmenting anyway and not sure how to use within BUGS code.

# Find number of animals per site
# filter by site
# filter by day
# individuals x trap - expand combos, max number of individuals per site
# spread
# fill into array (loop by site and day)

n_ind_site <- EM_CPIC %>%
  group_by(site, site_num) %>%
  select(site, site_num, id_site) %>%
  distinct() %>%
  summarise(n = max(id_site))

# assume capture rate of min_cap_rate (~0.03) to get max individuals to augment (M[g])
df_M <- n_ind_site %>%
  mutate(M = n / min_cap_rate + 10)
df_M

n_ind_site[order(n_ind_site$site_num), ] # orderin by the old site number and not the new site

n_sites <- length(unique(n_ind_site$site))
n_days <- 4

EM_CPIC_expanded <- EM_CPIC %>%
  expand(nesting(site, ind), trap, day) %>%
  left_join(EM_CPIC) %>%
  select(-sex) %>%
  ungroup() %>%
  group_by(site) %>%
  mutate(id_site = as.integer(as.factor(ind))) %>%
  ungroup() %>%
  mutate(site_num = as.integer(as.factor(site))) %>%
  select(-site) %>%
  left_join(n_traps_site) %>%
  mutate(count = if_else(is.na(count) & trap <= max_traps, 0, count)) # %>%

# expected sizes for each individual at each site x trap x day
n_ind_site$n * 14 * 4
sum(n_ind_site$n * 14 * 4)
sum(n_ind_site$n * 14 * 4) == nrow(EM_CPIC_expanded)

summary(EM_CPIC_expanded)


em_cpic_wide <- EM_CPIC_expanded %>%
  arrange(trap, site_num, id_site, day) %>%
  pivot_wider(names_from = trap, values_from = count, names_prefix = "trap_")

# expected sizes for each individual at each site x day
sum(n_ind_site$n * 4) == nrow(em_cpic_wide)
summary(em_cpic_wide)

# Separate individual characteristics to use if wanted
ind_covs <- EDF %>%
  group_by(site, ind, species, sex) %>%
  select(site, ind, sex, species, carapace, mass, trap) %>%
  summarise(carapace = mean(carapace),
            mass = mean(mass),
            n_measurements = n(),
            mean_trap = mean(trap))


ind_covs_cpic <- ind_covs %>%
  filter(species == "CPIC",
         !(site %in% c("H", "I")))

# Don't need to augment the data if using the zeros trick. They don't contribute to the likelihood - but to be in an array they need to be the same size so expand (augment) to the largest number of turtles caught at any site (to avoid using nesting)
EM_array <- array(NA, dim = c(max(n_ind_site$n), max(n_traps), n_days, n_sites))
for(l in 1:n_sites) {
  for(k in 1:n_days) {
    tmp <- em_cpic_wide %>%
    filter(site_num == l,
           day == k) %>%
      dplyr::select(starts_with("trap"))
    EM_array[1:nrow(tmp), , k, l] <- as.matrix(tmp)
  }
}

prod(dim(EM_array)) # size of array for the likelihood
length(EM_array[!is.na(EM_array)]) # Number of datapoints contributing to the likelihood

# get starting values 1 if indiviudal caught on any day at any trap for each site
# this needs to be the size of the augment not just the size of the EM_array captures
psi_st <- 0.2
Z_st <- matrix(rbinom(max(M), 1, psi_st), n_sites, max(M))
for(l in 1:n_sites) {
  Z_st[l, 1:n_ind_site$n[l]] <- 1 ##?
}
psi_sex_st <- 0.5

# Start values for s (activity centers) of augments (from random uniform constrained by state space size)

s_st_obs <- EM_CPIC_expanded %>%
  left_join(ind_covs_cpic) %>%
  filter(count != 0) %>%
  group_by(site_num, id_site) %>%
  select(site_num, id_site, trap) %>%
  summarise(loc = ((mean(trap) - 1) * 25) / 100 ) 

# make augmented starting location dataframe - will have NA for augments
s_st <- data.frame(expand.grid(site_num = 1:12, id_site = 1:max(M))) %>%
  left_join(s_st_obs) %>%
  as_tibble()

# randomly assign locations within the limits of the site to augmented individuals
s_st <- s_st %>%
  left_join(n_traps_site) %>%
  mutate(loc = ifelse(is.na(loc), runif(1, -1.5, ((((max_traps - 1) * 25) + 150) / 100)), loc))

summary(s_st)

# Select useful rows and spread to site x individual for analysis
s_st <- s_st %>%
  select(site_num, id_site, loc) %>%
  pivot_wider(names_from = id_site, values_from = loc) %>%
  ungroup() %>%
  select(-site_num) %>%
  as.matrix()

##### Sex vector divided by site

sex_array <- array(NA, dim = c(M, n_sites))

sex <- EM_CPIC %>%
  group_by(site_num, id_site, sex) %>%
  select(site_num, id_site, sex) %>%
  summarise_all(mean) %>%
  right_join(expand.grid(id_site = 1:max(M), site_num = 1:n_sites)) %>%
  mutate(sex2 = ifelse(sex == "M", 0, NA_integer_),
         sex2 = ifelse(sex == "F", 1, sex2)) %>%
  select(-sex) %>%
  rename(sex = sex2) 

summary(sex)
unique(sex$sex)

sex <- sex %>%
  pivot_wider(names_from = id_site, values_from = sex) %>%
  ungroup() %>%
  select(-site_num) %>%
  as.matrix()

#### Behavior Matrix ######

# individual x day x site 0 until caught, 1's after caught, don't necessarily need for augmented individuals

EM_CPIC_expanded
em_cpic_wide

recaps <- EM_CPIC_expanded %>%
  group_by(site, ind, id_site, site_num, day) %>%
  summarise(count = sum(count, na.rm = TRUE)) %>%
  group_by(site, ind, id_site, site_num) %>%
  mutate(caps = cumsum(count),
         recap = caps - count)

recaps <- recaps %>%
  ungroup() %>%
  select(site_num, id_site, day, recap) %>%
  mutate(recap = if_else(recap > 1, 1, recap)) %>%
  pivot_wider(names_from = day, values_from = recap, names_prefix = "day_")

recaptured <- array(0, dim = c(max(n_ind_site$n), n_days, n_sites))
for(l in 1:n_sites) {
  tmp <- recaps %>%
    filter(site_num == l) %>%
    select(starts_with("day")) %>%
    as.matrix()
  recaptured[1:nrow(tmp), 1:ncol(tmp), l] <- tmp
}

##### Augmented Zeros #####

augs <- matrix(0, n_sites, max(M))

#### Create matrix changing site IDs to be spatially relevant ####
# Matrix with current site IDs (which are sequential temporally) and spatially sequential site IDs


########## SAVE ALL OBJECTS NEEDED FOR MODEL ##########

if(!dir.exists("Data/Derived")) dir.create("Data/Derived", recursive = TRUE)

if(testing) {
  save(recaptured, Z_st, s_st, trap_locs, augs, sex, psi_st, psi_sex_st, EM_array, n_days, n_sites, n_traps_site, n_ind_site, M, xlim, run_date, df_M, file = paste0("Data/Derived/all_site_testing_", run_date, ".RData"))
} else {
  save(recaptured, Z_st, s_st, trap_locs, augs, sex, psi_st, psi_sex_st, EM_array, n_days, n_sites, n_traps_site, n_ind_site, M, xlim, df_M, file = "Data/Derived/all_site.RData") # other objects needed?
}

rm(list = ls())

