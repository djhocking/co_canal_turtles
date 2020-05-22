##### Multi-Session SCR Model Creation JAGS ####
##

## Also dubbed stratified population model

#### Adding in Variation in N per site ####
library(dplyr)
library(tidyr)
library(rgdal)
library(sp)

testing <- FALSE
run_date <- Sys.Date()
Species <- c("CPIC", "PRUB", "CSER", "SODO")[4] # change index depending on species

# number of possible individuals per site
M <- switch(Species,
            CPIC = 7000,
            PRUB = 1000,
            CSER = 500,
            SODO = 1000
)

if(testing) {
  M <- 100
}

min_cap_rate <- 0.01

##### Functions #####

nowt <- function(x = NULL) x

#####################

##### Load Data #####
site_index <- read.csv(file = "Data/site_index.csv", header = TRUE, stringsAsFactors = FALSE)

sites <- read.csv(file = "Data/trapids_sites.csv", header = TRUE, stringsAsFactors = FALSE)
sites <- left_join(sites, site_index, by = "site")
coords <- read.csv(file = "Data/coords.csv", stringsAsFactors = FALSE)
EDF <- read.csv(file = "Data/EDF.csv", stringsAsFactors = FALSE)
EDF <- left_join(EDF, site_index, by = "site")
n_traps_site <- read.csv(file = "Data/Max_Traps_Site.csv", stringsAsFactors = FALSE) # number of traps per site
n_traps_site <- left_join(n_traps_site, site_index, by = "site")
n_traps_site <- n_traps_site[order(n_traps_site$site_num), ]
n_traps <- n_traps_site$max_traps
# K <- max(EDF$day)
n_days <- max(EDF$day)

##### Make Combos of site-trap-days #####
# this will be used later to expand capture histories

combo_site_trap <- tidyr::expand_grid(site_num = n_traps_site$site_num, 
                                      trap = 1:max(n_traps)) %>%
  dplyr::left_join(n_traps_site) %>%
  dplyr::filter(trap <= max_traps)

combo_site_trap_day <- tidyr::expand_grid(site_num = n_traps_site$site_num, 
                                          trap = 1:max(n_traps),
                                          day = 1:n_days) %>%
  dplyr::left_join(n_traps_site) %>%
  dplyr::filter(trap <= max_traps)

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
  xlim[i, 1:2] <- c(min(trap_dist_list[[i]]) - 1000, max(trap_dist_list[[i]]) + 1000) / 100 # need to have buffer on each side without being negative. 
}

####### EDF FILE ########

# only a small number of CPIC were not sexed (~1%). They technically could be used but it would be really difficult with NIMBLE and the two other zero tricks and site looks that are already being applied. For this analysis, we will exclude those individuals from the data. Same with NA of which there are only two (SODO).

# filter(EDF, is.na(sex))
# filter(EDF, sex == "U")

# only do analysis for Males and Females. Skip for unrecorded and juvenile because insufficient data for calculating separate home range sizes and such for juveniles.
EDF <- EDF %>%
  filter(sex != "U",
         !is.na(sex))

#Take out sites H and I and get data only for 1 species
EDF_Sp <- EDF %>%
  filter(site != "H" & site != "I" & species == Species,
         ind != 70 | species != "CSER")
# EDF_Sp

##### Want EM ARRAY with ijk with index for site ########


EM <- EDF_Sp %>%
  group_by(site_num, ind, trap, day, sex) %>%
  select(site_num, ind, trap, day, sex) %>%
  mutate(count = 1) %>%
  summarise(count = sum(count)) %>%  ##?
  ungroup() %>%
  group_by(site_num) %>%
  mutate(id_site = as.integer(as.factor(ind))) %>%
  ungroup() %>%
  left_join(n_traps_site) %>% # changed to left join because missing will be handled with grid_expand
  select(-c(site, max_traps)) %>%
  # mutate(count = ifelse(is.na(count), 0, count)) %>%
  # mutate(day = ifelse(is.na(day), 1, day)) %>%
  # mutate(id_site = ifelse(is.na(id_site), 1, id_site)) %>% # is this needed???
  # mutate(trap = ifelse(is.na(trap), 13, trap)) %>%
  # mutate(trap = ifelse(trap == 13 & site_num == 5, 14, trap)) %>%
  # mutate(trap = ifelse(trap == 13 & site_num == 6, 11, trap)) ## change
  nowt()


# get number of unique individuals caught per site, with spatially relevant site IDs

# y[i,j,k,l] individual x trap x occassion x site
# to be the same size for an array will need to augment at least up to that while building array. Maybe could use nimbleList(). Not sure it's worth it since augmenting anyway and not sure how to use within BUGS code.

# Find number of animals per site
# filter by site
# filter by day
# individuals x trap - expand combos, max number of individuals per site
# spread
# fill into array (loop by site and day)

n_sites <- length(sites$site_num)
site_num_all <- as.data.frame(sites$site_num)
colnames(site_num_all) <- "site_num"

n_ind_site <- EM %>%
  group_by(site_num) %>%
  select(site_num, id_site) %>%
  distinct() %>%
  summarise(n = max(id_site)) %>%
  full_join(site_num_all) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>% # unnecessary now...
  arrange(site_num) %>%
  nowt()

# n_sites <- length(unique(n_ind_site$site_num))
n_sites <- 12
n_days <- 4

# assume capture rate of min_cap_rate (~0.03) to get max individuals to augment (M[g])
df_M <- n_ind_site %>%
  mutate(M = if_else(n > 0, n / min_cap_rate + 10, 500))
df_M

# 
# EM_expanded <- EM %>%
#   expand(nesting(site_num, ind), trap, day) %>%
EM_expanded <- combo_site_trap_day %>%
  left_join(EM) %>%
  # select(-id_site) %>%
  # mutate(ind = if)
  ungroup() %>%
  # unnest() %>%
  expand(nesting(site_num, ind, id_site), trap, day) %>%
  left_join(EM) %>%
  select(-sex) %>%
  ungroup() %>%
  # unnest() %>%
  # group_by(site_num) %>%
  # mutate(id_site = as.integer(as.factor(ind))) %>%
  # ungroup() %>%
  # mutate(site_num = as.integer(as.factor(site))) %>%
  left_join(n_traps_site) %>%
  select(-site) %>%
  mutate(count = if_else(is.na(count) & trap <= max_traps & !is.na(ind), 0, count)) %>%
  arrange(site_num, ind, trap, day) %>%
  nowt()

# expected sizes for each individual at each site x trap x day
n_ind_site$n * 14 * 4
sum(n_ind_site$n * 14 * 4)
sum(n_ind_site$n * 14 * 4) == nrow(EM_expanded) # FALSE FOR Species not at each site

summary(EM_expanded)

dplyr::filter(EM_expanded, site_num == 1 & day == 1)

em_wide <- EM_expanded %>%
  ungroup() %>%
  # mutate(ind = ifelse(is.na(ind), "A", ind),
  #        id_site = ifelse(is.na(id_site), -1, id_site)) %>%
  arrange(trap, site_num, id_site, day) %>%
  pivot_wider(names_from = trap, values_from = count, names_prefix = "trap_") # id_cols = c(site_num, day, ind), 

# dplyr::filter(em_wide, site_num == 1 & day == 1)
# expected sizes for each individual at each site x day
sum(n_ind_site$n * 4) == nrow(em_wide)
summary(em_wide)

# Separate individual characteristics to use if wanted
ind_covs <- EDF %>%
  group_by(site_num, site, ind, species, sex) %>%
  select(site_num, site, ind, sex, species, carapace, mass, trap) %>%
  summarise(carapace = mean(carapace),
            mass = mean(mass),
            n_measurements = n(),
            mean_trap = mean(trap)) %>%
  arrange(site_num, species, ind) %>%
  nowt()


ind_covs_sp <- ind_covs %>%
  filter(species == Species,
         !(site %in% c("H", "I")))

# Don't need to augment the data if using the zeros trick. They don't contribute to the likelihood - but to be in an array they need to be the same size so expand (augment) to the largest number of turtles caught at any site (to avoid using nesting)
EM_array <- array(NA, dim = c(max(n_ind_site$n), max(n_traps), n_days, n_sites))
for(l in 1:n_sites) {
  for(k in 1:n_days) {
    tmp <- em_wide %>%
      filter(site_num == l,
             day == k,
             !is.na(ind)) %>%
      dplyr::select(starts_with("trap"))
    if(nrow(tmp) == 0) next
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

s_st_obs <- EM_expanded %>%
  left_join(ind_covs_sp) %>%
  filter(count != 0) %>%
  group_by(site_num, id_site) %>%
  select(site_num, id_site, trap) %>%
  arrange(site_num, id_site, trap) %>%
  summarise(loc = ((mean(trap) - 1) * 25) / 100 )

# make augmented starting location dataframe - will have NA for augments
s_st <- data.frame(expand.grid(site_num = 1:12, id_site = 1:max(M))) %>%
  left_join(s_st_obs) %>%
  as_tibble() %>%
  arrange(site_num, id_site)

# randomly assign locations within the limits of the site to augmented individuals
s_st <- s_st %>%
  left_join(n_traps_site) %>%
  ungroup() %>%
  as.data.frame()

n_missing <- s_st %>%
  filter(is.na(loc)) %>%
  nrow()

# s_st <- s_st %>%
#   mutate(loc = if_else(is.na(loc), runif(n_missing, -1.5, ((((min(n_traps_site$max_traps) - 1) * 25) + 150) / 100)), loc)) 

s_st[is.na(s_st)] <- runif(n_missing, -1.5, ((((min(n_traps_site$max_traps) - 1) * 25) + 150) / 100))

summary(s_st)
hist(s_st$loc)

# Select useful columns and spread to site x individual for analysis
s_st <- s_st %>%
  select(site_num, id_site, loc) %>%
  pivot_wider(names_from = id_site, values_from = loc) %>%
  ungroup() %>%
  select(-site_num) %>%
  as.matrix()

##### Sex vector divided by site

sex_array <- array(NA, dim = c(M, n_sites))

sex <- EM %>%
  group_by(site_num, id_site, sex) %>%
  select(site_num, id_site, sex) %>%
  summarise_all(mean) %>%
  right_join(expand.grid(id_site = 1:max(M), site_num = 1:n_sites)) %>%
  mutate(sex2 = ifelse(sex == "M", 0, NA_integer_),
         sex2 = ifelse(sex == "F", 1, sex2)) %>%
  arrange(site_num, id_site) %>%
  select(-sex) %>%
  rename(sex = sex2) 

summary(sex)
unique(sex$sex)

sex <- sex %>%
  ungroup() %>%
  as.data.frame() %>%
  pivot_wider(names_from = id_site, values_from = sex) %>%
  ungroup() %>%
  select(-site_num) %>%
  as.matrix()

#### Behavior Matrix ######

# individual x day x site 0 until caught, 1's after caught, don't necessarily need for augmented individuals

EM_expanded
em_wide

recaps <- EM_expanded %>%
  group_by(site_num, ind, id_site, day) %>%
  summarise(count = sum(count, na.rm = TRUE)) %>%
  group_by(site_num, ind, id_site) %>%
  mutate(caps = cumsum(count),
         recap = caps - count)

recaps <- recaps %>%
  ungroup() %>%
  select(site_num, id_site, day, recap) %>%
  mutate(recap = if_else(recap > 1, 1, recap)) %>%
  dplyr::filter(!is.na(id_site)) %>%
  arrange(site_num, id_site, day) %>%
  pivot_wider(names_from = day, values_from = recap, names_prefix = "day_")

recaptured <- array(0, dim = c(max(n_ind_site$n), n_days, n_sites))
for(l in 1:n_sites) {
  tmp <- recaps %>%
    filter(site_num == l) %>%
    select(starts_with("day")) %>%
    as.matrix()
  if(nrow(tmp) == 0) next
  recaptured[1:nrow(tmp), 1:ncol(tmp), l] <- tmp
}

##### Augmented Zeros #####

augs <- matrix(0, n_sites, max(M))

###### Number of individuals caught per day across sites #####

caps_day <- EM_expanded %>%
  ungroup() %>%
  group_by(day) %>%
  summarise(caps = sum(count, na.rm = TRUE))

#### Create matrix changing site IDs to be spatially relevant ####
# Matrix with current site IDs (which are sequential temporally) and spatially sequential site IDs



########## SAVE ALL OBJECTS NEEDED FOR MODEL ##########

if(!dir.exists(paste0("Data/Derived/", Species, "/"))) dir.create(paste0("Data/Derived/", Species, "/"), recursive = TRUE)

if(testing) {
  save(Species, recaptured, Z_st, s_st, trap_locs, augs, sex, psi_st, psi_sex_st, EM_array, n_days, n_sites, n_traps_site, n_ind_site, M, xlim, run_date, df_M, caps_day, file = paste0("Data/Derived/", Species, "_all_site_testing_", run_date, ".RData"))
} else {
  save(Species, recaptured, Z_st, s_st, trap_locs, augs, sex, psi_st, psi_sex_st, EM_array, n_days, n_sites, n_traps_site, n_ind_site, M, xlim, run_date, df_M, caps_day, file = paste0("Data/Derived/", Species, "_all_site_reg_", run_date, ".RData"))
}

rm(list = ls())
