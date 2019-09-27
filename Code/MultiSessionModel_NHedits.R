##### Multi-Session SCR Model Creation JAGS ####
##

## Also dubbed stratified population model

#### Adding in Variation in N per site ####
library(dplyr)
library(tidyr)
library(rjags)
library(parallel)
library(zoo)
library(rgdal)
library(sp)
library(utils)
# library(reshape)
# library(plyr)

testing <- FALSE
run_date <- Sys.Date()

n_traps <- 14

# number of possible individuals per site

M <- 5000 

if(testing) {
  M <- 700
}

Sites <- read.csv(file = "Data/trapids_sites.csv", header = TRUE)

max_trap_csv <- read.csv(file = "Data/Max_Traps_Site.csv")
max_trap <- max_trap_csv$max_traps

coords <- read.csv(file = "Data/coords.csv")
str(coords)
summary(coords)

trap_locs_degrees <- coords
trap_locs_degrees$trap <- 1:nrow(trap_locs_degrees)
trap_num <- trap_locs_degrees$trap ## Change?? 122 right now


# convert to utm to have distance in meters
coords_dd = SpatialPoints(coords[ , c("lon", "lat")], proj4string=CRS("+proj=longlat"))
coords_utm <- spTransform(coords_dd, CRS("+init=epsg:26917"))

trap_locs <- coords_utm
trap_locs <- as.data.frame(trap_locs)
trap_locs <- cbind(trap_num, trap_locs)
colnames(trap_locs) = c("trap_id", "easting", "northing")
# trap_locs as single vector with distance between

## Creating trap location vector per site using coordinates (sp package required)
trap_dist_A <- spDistsN1(coords_utm[1:8, ], coords_utm[1, ])
trap_dist_C <- spDistsN1(coords_utm[9:18, ], coords_utm[9, ])
trap_dist_D <- spDistsN1(coords_utm[19:26, ], coords_utm[19, ])
trap_dist_E <- spDistsN1(coords_utm[27:40, ], coords_utm[27, ])
trap_dist_F <- spDistsN1(coords_utm[41:47, ], coords_utm[41, ])
trap_dist_G <- spDistsN1(coords_utm[48:54, ], coords_utm[48, ])
trap_dist_J <- spDistsN1(coords_utm[61:70, ], coords_utm[61, ])
trap_dist_K <- spDistsN1(coords_utm[71:80, ], coords_utm[71, ])
trap_dist_L <- spDistsN1(coords_utm[81:90, ], coords_utm[81, ])
trap_dist_M <- spDistsN1(coords_utm[91:102, ], coords_utm[91, ])
trap_dist_N <- spDistsN1(coords_utm[103:112, ], coords_utm[103, ])
trap_dist_O <- spDistsN1(coords_utm[113:122, ], coords_utm[113, ])

trap_dist_list <- list(trap_dist_A, trap_dist_C, trap_dist_D, trap_dist_E,
                       trap_dist_F, trap_dist_G, trap_dist_J, trap_dist_K,
                       trap_dist_L, trap_dist_M, trap_dist_N, trap_dist_O)
# trap_locs <- trap_dist_list

trap_locs <- matrix(NA, 12, max(max_trap))
for (i in 1:12) {
trap_locs[i, 1:max_trap[i]] <- trap_dist_list[[i]] / 100
}



xlim <- matrix(NA, 12, 2)
for(i in 1:12){
  xlim[i, 1:2] <- c(min(trap_dist_list[[i]]) - 150, max(trap_dist_list[[i]]) + 150) / 100 # need to have buffer on each side without being negative. Just added 50 to the end for testing but will have to think through
}


####### EDF FILE ########

EDF <- read.csv(file = "Data/EDF.csv", stringsAsFactors = FALSE)
head(EDF)
summary(EDF)

#Take out sites H and I

EDF_CPIC <- EDF %>%
  filter(site != "H" & site != "I" & species == "CPIC")
EDF_CPIC

## subtract 6 from trap ids > = 61 (Sites H and I)

#Add a new column for integer session values (session = site)

EDF_CPIC$trap_id_edited <- ifelse(EDF_CPIC$trap_id >= 61, EDF_CPIC$trap_id - 6, EDF_CPIC$trap_id - 0)
EDF_CPIC$site_num <- as.integer(as.factor(EDF_CPIC$site))

summary(EDF_CPIC)


####### WORK ON THIS SECTION NEXT 3_11_19 ########

# traplocsA <- traplocsA / 100
# matrixA <- matrixA / 100 # scale for computational purposes
# 
n_traps_site <- read.csv(file = "Data/Max_Traps_Site.csv") # number of traps per site
n_traps <- n_traps_site$max_traps

########

#N <- nrow(EDF_CPIC[which(EDF_CPIC$recap == "N"), ])

N_persite <- list()
N <- list()
for(i in 1:12) {
  N_persite[[i]] <- nrow(EDF_CPIC[which(EDF_CPIC$recap == "N" & EDF_CPIC$site_num == i), ])
  N[[i]] <- N_persite[[i]]
}
n_ind <- N
n_ind_total <- nrow(EDF_CPIC[which(EDF_CPIC$recap == "N"), ])
do.call(sum, N)
K <- max(EDF_CPIC$day) # trap nights per session
# buffer <- 1 # check literature to make sure doesn't need to be larger

##n_ind <- length(unique(EDF_CPIC$ind)) ## Needs to match up with N? 3 off...
## should't use unique as does not count those that were b/w sites; why counting? Codes are different...

traps_per_site <- read.csv(file = "Data/trapids_sites.csv")


# Make encounter histories with number of times each individual is captured in each trap
##### Want EM ARRAY with ijk with index for site ########

########

str(EDF)
EDF_CPIC

EM_CPIC <- EDF_CPIC %>%
  group_by(site_num, ind, trap_id_edited, day) %>%
  select(site_num, ind, trap_id_edited, day) %>%
  mutate(count = 1) %>%
  summarise_all(sum) %>%
  #spread(trap_id_edited, count, fill = 0) %>%
  ungroup() %>%
  mutate(id = as.integer(as.factor(ind)))

EM_CPIC$site_trap <- ave(EM_CPIC$trap_id_edited, EM_CPIC$site_num, FUN = function(x) as.numeric(factor(x)))


site_trap_combos <- expand.grid(site_num = 1:12, site_trap = 1:14) %>%
  arrange(site_num, site_trap) # ... need to add in extra traps per site, need to go through and label trap ids with 14 trap "gap" per site


foo <- site_trap_combos  %>%
   left_join(EM_CPIC)

foo <- foo %>%
left_join(select(max_trap_csv, -site)) %>%
  # mutate(count = ifelse(site_trap > max_trap, NA_integer_, count)) %>%
  mutate(count = ifelse(site_trap <= .$max_traps & is.na(count), 0, count))

#add in augments

foo_spread <- foo %>%
spread(site_trap, count, fill = 0)

full_df <- tidyr::expand(foo_spread, id, day, site_num)
full_df <- na.omit(full_df)

EM <- left_join(full_df, foo_spread)


EM <- as.data.frame(EM, stringsAsFactors = FALSE)
EM <- na.omit(EM)
EM <- as.data.frame(EM, stringsAsFactors = FALSE)

############

#n_traps <- 14
J <- n_traps
num_sites <- max(EDF_CPIC$site_num)
G <- num_sites

###########

###############
#19_4_19

EM_array <- array(NA, dim = c(M, (max(max_trap) + 1), K, G))

# make 4D array: individual (M) x trap (n_traps) x day (K) x site (G)
# split by day
for(k in 1:K) {
  for(g in 1:G) {
  foo <- EM[(which(EM[]$day == k & EM$site_num == g)), ]
  foo_less <- select(foo, -c(site_num, ind, day, trap_id_edited, max_traps))
  df_aug <- as.data.frame(matrix(0, nrow = (M - nrow(foo_less)), ncol = (max(max_trap) + 1)), stringsAsFactors = FALSE)
  # df_aug$site_num <- g
  # df_aug <- df_aug[ , c(ncol(df_aug), 1:(ncol(df_aug)-1))]
  colnames(df_aug) <- colnames(foo_less)
  foo_augment <- bind_rows(foo_less, df_aug)
  foo_augment$id <- ifelse(foo_augment$id == 0, NA, foo_augment$id)
  EM_array[ , , k, g] <- as.matrix(foo_augment)
  }
}

foob <- EM_array[ , -1, , ]

target <- c(1:M)

  for (k in 1:K) {
    for (g in 1:G) {
    food <- EM_array[ , , k, g]
    foob_arranged <- food[match(target, food[ ,1]), ]
    foob_arranged[is.na(foob_arranged)] <- 0
    foob_arranged[ , 1] <- ifelse(foob_arranged[ , 1] == 0, NA, foob_arranged[ ,1])
    foob_arranged <- foob_arranged[ , -1]
    #foob_arranged <- select()
    #foob_arranged <- food[order(food[ , 1]), ]
    foob[ , , k, g] <- as.matrix(foob_arranged)
    }
  }

EM_array <- foob ## WHY 8's in first column??

## Need to divide M by 4, thus change ids?
#### 500 per day per site? or 500 per site (500/4 per day?), also it does not keep the same id for the same indiiduals between days... so recap calc will not work?

## EM_array_2 --> collapsed day so could make sst vector with all unique individuals caught per site; made sure each individual was only counted (recaps not counted)

# get starting values 1 if indiviudal caught on any day at any trap for each site


z <- matrix(NA, G, M)
bar <- matrix(NA, K, M)
for(g in 1:G) {
  for(k in 1:K) {
    foog <- as.matrix(EM_array[ , , k, g])
    bar[k, ] <- apply(foog, 1, max, na.rm = TRUE)
  }
  z[g, ] <- apply(bar, 2, max, na.rm = TRUE)
}

#Start values for s (activity centers) of augments (from random uniform constrained by state space size)
#X <- traplocsA
# Now populated by starting positions uniformally placed within state space
# For every individual that is not 0, change starting x point to mean of traps associated with encounters for that individual; leaves 0's there from the augmented population and also puts in activity center for augmented individuals that were randomly given an encounter history (caught at least 1 time)

sum_caps <- apply(EM_array, c(1,2,4), sum)  ## c(1,2): dimensions to apply function to; 1 = rows, 2 = columns; collapsed day in this instance
## this is wrong, it doesn't distinguish rows per day as different individuals thus the max number of individuals caught one day becomes the total number of caught inviduals here... Need to sum each individual per for each site


traplocsE <- as.matrix(trap_locs[4, ])
row.names(traplocsE) <- NULL
colnames(traplocsE) <- NULL

unif_array <- array(NA, dim = c(M, G))

for (g in 1:G){
x <- runif(M, min = xlim[g, ][1], max = xlim[g, ][2])
unif_array[ , g] <- as.matrix(x)
}


sst <- matrix(NA, M, G)

for(g in 1:G){
    sst[ , g] <- (sum_caps[ , , g] %*% traplocsE) / (ifelse(rowSums(sum_caps[ , , g]) > 0, rowSums(sum_caps[ , , g]), 1))
}

for(i in 1:M){
  for(g in 1:G){
sst[i, g] <- ifelse(sst[i , g] == 0, unif_array[i, g], sst[i, g])
  }
}


##### Sex vector divided by site

sex_array <- array(NA, dim = c(M, G))

# make 4D array: individual (M) x trap (n_traps) x day (K) x site (G)
# split by day
EM_CPIC_sex <- EDF_CPIC %>%
  group_by(site_num, ind, trap_id_edited, day, sex) %>%
  select(site_num, ind, trap_id_edited, day, sex) %>%
  mutate(count = 1) %>%
  summarise_all(sum) %>%
  #spread(trap_id_edited, count, fill = 0) %>%
  ungroup() %>%
  mutate(id = as.integer(as.factor(ind)))


######
EM_array2 <- matrix(NA_integer_, M, G)
target <- c(1:M)

# make 4D array: individual (M) x trap (n_traps) x day (K) x site (G)
# split by day
for(k in 1:K) {
  for(g in 1:G) {
    foo <- EM_CPIC_sex[(which(EM_CPIC$site_num == g)), ]
    foo_less <- select(foo, -c(site_num, ind, day, trap_id_edited, count))
    foo_less <- foo_less %>%
      distinct %>%
      mutate(sex = ifelse(sex == "U", NA, sex),
             sex = ifelse(sex == "M", 0, sex),
             sex = ifelse(sex == "F", 1, sex)) %>%
      mutate(sex = as.integer(sex))
    df_aug <- as.data.frame(matrix(NA_integer_, nrow = (M - nrow(foo_less)), ncol = 2), stringsAsFactors = FALSE)
    # df_aug$site_num <- g
    # df_aug <- df_aug[ , c(ncol(df_aug), 1:(ncol(df_aug)-1))]
    colnames(df_aug) <- colnames(foo_less)
    foo_augment <- bind_rows(foo_less, df_aug)
    target <- as.data.frame(1:M)
    colnames(target) <- "id"
    foo_augment <- left_join(target, foo_augment, by = "id")
    foo_augment <- select(foo_augment, -id)
    #foob_arranged <- foo_augment[match(target, foo_augment[ , 2]), ]
    EM_array2[ , g] <- as.matrix(foo_augment)
  }
}

Sex <- t(EM_array2)


#### Behavior Matrix ######

BM <- EDF_CPIC %>%
  group_by(site_num, ind, day, recap) %>%
  select(site_num, ind, day, recap) %>%
  ungroup() %>%
  mutate(id = as.integer(as.factor(ind)))

str(BM)
BM

BM <- as.data.frame(BM, stringsAsFactors = FALSE)

BM$behav <- ifelse(BM$recap == "N", 1, 0)

BM_less <- select(BM, -recap)

B_array <- array(NA, dim = c(M, 2, K, G))

for (i in 1:M){
  for(k in 1:K){
    for(g in 1:G){
      foo <- BM_less[(which(BM_less[]$day == k & BM_less$site_num == g)), ]
      foo_less <- select(foo, -c(site_num, ind, day))
      df_aug <- as.data.frame(matrix(0, nrow = (M - nrow(foo_less)), ncol = 2), stringsAsFactors = FALSE)
      colnames(df_aug) <- colnames(foo_less)
      foo_augment <- bind_rows(foo_less, df_aug)
      foo_augment$id <- ifelse(foo_augment$id == 0, NA, foo_augment$id)
      B_array[ , , k, g] <- as.matrix(foo_augment)
    }
  }
}

foob <- array(NA, dim = c(M, K, G))

target <- c(1:M)

for (k in 1:K) {
  for (g in 1:G) {
    food <- B_array[ , , k, g]
    foob_arranged <- food[match(target, food[ ,1]), ]
    foob_arranged <- foob_arranged[ , -1]
    foob_arranged[is.na(foob_arranged)] <- 0
    #foob_arranged[ , 1] <- ifelse(foob_arranged[ , 1] == 0, NA, foob_arranged[ ,1])
    #foob_arranged <- select()
    #foob_arranged <- food[order(food[ , 1]), ]
    foob[ , k, g] <- as.matrix(foob_arranged)
  }
}

# make all after first capture = 1
for(g in 1:G) {
  for(m in 1:M) {
    for(k in 2:K) {
      foob[m, k, g] <- ifelse(foob[m, k-1, g] >= 1, 2, foob[m, k, g])
    }
  }
}

foob <- ifelse(foob == 1, 0, foob)
C <- ifelse(foob == 2, 1, foob)

n_sites <- G


max_ind_sp <- max(EM_array[ , , , ], na.rm = TRUE)

########## SAVE ALL OBJECTS NEEDED FOR MODEL ##########

if(!dir.exists("Data/Derived")) dir.create("Data/Derived", recursive = TRUE)

if(testing) {
  save(z, sst, n_sites, EM_array, Sex, trap_locs, K, M, xlim, max_trap, C, G, run_date, max_ind_sp, file = paste0("Data/Derived/all_site_testing_", run_date, ".RData"))
} else {
  save(z, sst, n_sites, EM_array, Sex, trap_locs, K, M, xlim, max_trap, C, G, run_date, max_ind_sp, file = "Data/Derived/all_site.RData") # other objects needed?
}

