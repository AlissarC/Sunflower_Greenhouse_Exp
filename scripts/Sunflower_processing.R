
###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(LeafArea)
library(rJava)

ij.path <- "../ImageJ/"
imagepath.disk <- "../input/leaves_scan/leaves_disk/" # 10 disks per plant at three time points 
imagepath.leaves <- "../input/leaves_scan/leaves/" ## 3-8 leaves per plant at the end of the experiment 

leaves_disk <- run.ij(path.imagej = ij.path,
                      set.directory = imagepath.disk,
                      distance.pixel = 118,
                      known.distance = 1, low.size = 0.01)
view(leaves_disk)

leaves_area <- run.ij(path.imagej = ij.path,
                      set.directory = imagepath.leaves,
                      distance.pixel = 118,
                      known.distance = 1, low.size = 0.01)

view(leaves_area)


leaves_disk <- separate(leaves_disk, col = sample, into=c("ID_trt", "n_trt", "Drought_trt", "plant_id", "time", "type"), sep= "_", remove=FALSE)
view(leaves_disk) # Verification 
names(leaves_disk)

leaves_area <- separate(leaves_area, col = sample, into=c("ID_trt", "n_trt", "Drought_trt", "plant_id", "leaf"), sep= "_", remove=FALSE)
view(leaves_area) # Verification 
names(leaves_area)

leaves_disk$ID_trt = NULL # remove the ID_trt, it is not useful
leaves_disk$sample = NULL # remove the samples, it is not useful anymore
leaves_disk$plant_id<- tolower(leaves_disk$plant_id) # make capital P small p
leaves_disk$type= NULL # remove the type, it is not useful anymore
names(leaves_disk)[names(leaves_disk) == "total.leaf.area"] <- "leaves_disk_area"  # rename total.leaf.area"


leaves_area$ID_trt = NULL # remove the ID_trt, it is not useful
leaves_area$sample = NULL # remove the samples, it is not useful anymore
leaves_area$plant_id <- tolower(leaves_area$plant_id) # make capital P small p
leaves_area$type= NULL # remove the type, it is not useful anymore

aggregate(leaf ~ plant_id + n_trt + Drought_trt,
          data = leaves_area,
          FUN = function(x) sum(!is.na(x)))
names(leaves_area)

####################### calculate the leaf area for disks used to measure nitrogen and C isotopes####################
####################################################################################################################

#Join files containing the number and biomass of disks: nb_disk_plant_time ###################### 
## In this file the biomass is for the dried disks (leaves_disk_mg_noChl) used for N and C isotopes measurements
#we scanned 8-10 disks for each plant, some were used to chlorophyll measurements and some for SLA

nb_disk_plant_time <-  read.csv("../input/nb_disk_plant_time.csv")
names(nb_disk_plant_time)

## add p to plant_id
nb_disk_plant_time$plant_id <- paste0("p", nb_disk_plant_time$plant_id)

leaves_disk_area= left_join(nb_disk_plant_time,leaves_disk, by = c("n_trt", "Drought_trt", "plant_id", "time" )) 
names(leaves_disk_area)
view(leaves_disk_area)

#### Calculate SLA ############

### disk area for disks dried and used for SLA (the disk_area in this df is for all disks including the disks used of chlorophyll) 

leaves_disk_area$leaf_area_for_SLA = (leaves_disk_area$nb_disk_SLA*leaves_disk_area$leaves_disk_area)/leaves_disk_area$nb_disks_scanned
hist(leaves_disk_area$leaf_area_for_SLA)
## no disks scanned at t1, then no SLA
## area in cm2, dry weight of disks in mg
##*10 to convert from mg/cm2 to g/m2, calculation of SLA and LMA
leaves_disk_area$LMA = (leaves_disk_area$leaves_disk_mg_noChl/leaves_disk_area$leaf_area_for_SLA)*10 
leaves_disk_area$SLA = 1/leaves_disk_area$LMA #m2/g

hist(leaves_disk_area$LMA)
hist(leaves_disk_area$SLA)
names(leaves_disk_area)

## rearrange columns
leaves_disk_area <- leaves_disk_area[, c(
  "time", "n_trt", "Drought_trt", "plant_id",
  setdiff(names(leaves_disk_area), c("time", "n_trt", "Drought_trt", "plant_id"))
)]
names(leaves_disk_area)


### load Isotopes, leaf nitrogen, and A-Ci data ################
A_Ci_N_isotopes <-  read.csv("../input/A_Ci_N_isotopes.csv") 
names(A_Ci_N_isotopes)
nrow(A_Ci_N_isotopes) #192
table(A_Ci_N_isotopes$time)

A_Ci_N_isotopes$plant_id <- paste0("p",A_Ci_N_isotopes$plant_id)

################ Big Delta calculus #####################################
delta_air = -8
A_Ci_N_isotopes$big_D13 = (delta_air - A_Ci_N_isotopes$deltaC13)/(1+A_Ci_N_isotopes$deltaC13*0.001)
hist(A_Ci_N_isotopes$big_D13)

######################## Chi calculus ###############################
a = 4.4
b = 28
A_Ci_N_isotopes$chi= (A_Ci_N_isotopes$big_D13 -a)/(b-a)
hist(A_Ci_N_isotopes$chi)

names(A_Ci_N_isotopes)
view(A_Ci_N_isotopes)

names(leaves_disk_area)
view(leaves_disk_area)
################ Merge with SLA/LMA data ################## 

Iso_leaf_traits = left_join(leaves_disk_area, A_Ci_N_isotopes, by = c("time", "n_trt", "Drought_trt", "plant_id"))
                                                           
names(Iso_leaf_traits)
view(Iso_leaf_traits)

######## calculate Nmass gN/gleaf
Iso_leaf_traits$Nmass= Iso_leaf_traits$percentage_N/100
hist(Iso_leaf_traits$Nmass)

######## calculate Narea gN/m2
Iso_leaf_traits$Narea= Iso_leaf_traits$Nmass * Iso_leaf_traits$LMA
hist(Iso_leaf_traits$Narea)

############################# call functions to calculate beta and chi #################################################
install.packages ("R.utils") # to read source function
library (R.utils)
sourceDirectory("../scripts/functions", modifiedOnly = FALSE, verbose = TRUE)
################ calculate atmospheric pressure (Pa) from elevation (m)#####################
Iso_leaf_traits$z = 976
calc_patm = function(z) {
  
  kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
  kTo = 298.15   # base temperature, K (Prentice, unpublished)
  kL = 0.0065    # temperature lapse rate, K/m (Allen, 1973)
  kG = 9.80665   # gravitational acceleration, m/s**2 (Allen, 1973)
  kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
  kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
  
  patm = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))
  
  patm
}

Iso_leaf_traits$patm<- calc_patm(Iso_leaf_traits$z)

################ calculate  nstar (unitless relative viscosity of h2o at temperature relative to 25°C) ####################
calc_nstar = function(temp, z){ # temp in °C and z in m
  
  patm = calc_patm(z)
  
  # viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa) 
  ns      = calc_viscosity_h2o( temp, z )  # Pa s 
  ns25    = calc_viscosity_h2o( 25, z )  # Pa s 
  nstar = ns / ns25                       # (unitless)
  
  nstar
  
}
# Apparently it call automatically calc_viscositt_h2o 
Iso_leaf_traits$nstar<- calc_nstar(temp = Iso_leaf_traits$tmp, z = Iso_leaf_traits$z)

########################## calculate gammastar (Pa) ###############################################
calc_gammastar_pa = function(temp, z) {
  
  patm = calc_patm(z)
  rat = calc_patm(z) / calc_patm(0)
  gammastar25 = 4.332 * rat  # Pa
  Hgm=37830 # J mol-1
  R = 8.314        # J K-1 mol-1
  O2 = 2.09476e5 # ppm
  O2_0 = O2 * 1e-6 * calc_patm(0)
  O2_z = O2 * 1e-6 * calc_patm(z)
  
  temp_k = 273.15+ temp
  
  gStar_pa = gammastar25*exp((Hgm/R)*(1/298.15-1/temp_k))
  
  gStar_pa
  
}
Iso_leaf_traits$gammastar_pa <-calc_gammastar_pa(temp = Iso_leaf_traits$tmp, z = Iso_leaf_traits$z)

######## calcualte the Michaelis-Menton coefficient (Pa) for Rubisco from temperature##############################
calc_km_pa = function(temp, z) {
  
  patm = calc_patm(z) 
  rat = patm / calc_patm(0)
  
  R = 8.314        
  O2 = 2.09476e5      
  Kc25 = 41.03 * rat 
  Ko25 = 28210 * rat 
  Hkc = 79430  
  Hko = 36380 
  
  temp_k = 273.15 + temp
  
  Kc_pa =Kc25 * exp(Hkc * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  Ko_pa =Ko25* exp(Hko * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  
  O2_pa = O2 * (1e-6) * patm 
  
  Km_pa = Kc_pa * (1 + O2_pa/Ko_pa)
  
  Km_pa 
  
}

Iso_leaf_traits$Km<-calc_km_pa(temp = Iso_leaf_traits$tmp, z = Iso_leaf_traits$z)

names(Iso_leaf_traits)
##################### Beta calculus ######################################
Iso_leaf_traits$D <- Iso_leaf_traits$VPD*1000
Iso_leaf_traits$ca <- 420
Iso_leaf_traits$beta_num <- 1.6 * Iso_leaf_traits$nstar * Iso_leaf_traits$D * ((Iso_leaf_traits$chi-Iso_leaf_traits$gammastar_pa/Iso_leaf_traits$ca) ^ 2)
Iso_leaf_traits$beta_denom <- ((1 - Iso_leaf_traits$chi)^2) * (Iso_leaf_traits$Km + Iso_leaf_traits$gammastar_pa)

Iso_leaf_traits$betaleaf <- Iso_leaf_traits$beta_num / Iso_leaf_traits$beta_denom

summary(Iso_leaf_traits$beta_num)
hist(Iso_leaf_traits$beta_num)
hist(Iso_leaf_traits$betaleaf)
plot(Iso_leaf_traits$beta_num, Iso_leaf_traits$chi)
plot(Iso_leaf_traits$beta_denom, Iso_leaf_traits$chi)
plot(Iso_leaf_traits$chi, Iso_leaf_traits$betaleaf)
plot(Iso_leaf_traits$chi, log(Iso_leaf_traits$betaleaf))

#################### Load transpiration and drought data #################
transpiration <- read.csv("../input/transpiration_all.csv")
names(transpiration)
transpiration$plant_id <- paste0("p", transpiration$plant_id)

# Extract the values for each individual plant at the end of the experiment of the cumulative Transpiration, Nitrogen and DS daily
transpiration_final <- transpiration %>%
  group_by(n_trt, Drought_trt, plant_id) %>%
  filter(DAD == max(DAD, na.rm = TRUE)) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  select(n_trt,Drought_trt,plant_id,DAD,
         transpiration_cumul_L,Cumulative_N_mg, Sum_DS_daily)
transpiration_final
names(transpiration_final)
transpiration_final <- transpiration_final %>%
  rename(DAD_final = DAD)
view(transpiration_final)

################### calculate the total leaf area per plant ################

########### open the file containing dimensions of leaves scanned at t4 ###########
leaves_WL <- read.csv("../input/Leaves_scan_dimensions.csv")
leaves_WL
leaves_WL <- separate(leaves_WL, col = sample, into=c("ID_trt", "n_trt", "Drought_trt", "plant_id", "leaf"), sep= "_", remove=FALSE)
view(leaves_WL) # Verification 
leaves_WL$ID_trt = NULL # remove the ID_trt, it is not useful
leaves_WL$sample = NULL # remove the samples, it is not useful anymore
leaves_WL$plant_id<- tolower(leaves_WL$plant_id) # make capital P small p

################# Join with the file leaf area file 
names(leaves_area)
leaves_WL_area = left_join(leaves_WL,leaves_area, by = c("n_trt", "Drought_trt", "plant_id", "leaf")) 
names(leaves_WL_area)
view(leaves_WL_area) 

leaves_WL_area$WxL = leaves_WL_area$Length * leaves_WL_area$Width  ## calculate length * width

plot(leaves_WL_area$WxL, leaves_WL_area$total.leaf.area) ## linear relationship

model <- lm(total.leaf.area ~ WxL, data = leaves_WL_area)

summary(model)
summary(model)$r.squared ## 0.98
coef(model)

#Equation : leaf_area = 0.692*WxL + 3.1065

### equation deduced from the relationship between leaf scans and the product of their width and length : R2=0.98

### now loaf the file containing the length and width of each leaf/plant/date
morpho <- read.csv("../input/morpho.csv")
morpho
########### calculate length X width for each leaf/plant/time
morpho$WxL = morpho$length * morpho$width
hist(morpho$WxL)

## calculate the leaf area for each time using the deduced equation
morpho$leaf_area = 0.692*morpho$WxL + 3.1065 
hist(morpho$leaf_area)

names(morpho)
morpho$plant_id
#################  calculate the plant surface area by summing all leaves/plant/time 
morpho_all_sum_plant <- morpho %>%  
  group_by (n_trt,Drought_trt,plant_id, time) %>% 
  summarize(plant_surface_area = sum (leaf_area, na.rm = TRUE))
morpho_all_sum_plant

view(morpho_all_sum_plant)
names(morpho_all_sum_plant)


############################ Load biomass data ######################
biomass <- read.csv("../input/Biomass_final.csv")
names (biomass)
biomass$plant_id <- paste0("p", biomass$plant_id)

#################### Join biomass to plant surface area ###################
biomass_surface = left_join(biomass,morpho_all_sum_plant, by = c("n_trt", "Drought_trt", "plant_id",
                                                                   "time")) 
names(biomass_surface)

### join to transpiration #############
biomass_surface_transpiration = left_join(biomass_surface,transpiration_final, by = c("n_trt", "Drought_trt", "plant_id")) 
names(biomass_surface_transpiration)
view(biomass_surface_transpiration)
##### Join to iso_N_A_Ci data ####################
names(Iso_leaf_traits)
names(biomass_surface_transpiration)
sunflower_all_data = left_join(Iso_leaf_traits, biomass_surface_transpiration, by = c("time","n_trt", "Drought_trt", "plant_id"))
sunflower_all_data
names(sunflower_all_data)
view(sunflower_all_data)

#### calculate sum biomass ###########

sunflower_all_data$Abg_all_dry = sunflower_all_data$dry_leaves_petioles+sunflower_all_data$dry_stem+sunflower_all_data$dry_flowers+sunflower_all_data$dry_senesced_leaves
hist(sunflower_all_data$Abg_all_dry) # Dry Aboveground biomass in g

sunflower_all_data$Abg_all_fresh = sunflower_all_data$fresh_leaves_petioles+sunflower_all_data$fresh_stem+sunflower_all_data$fresh_flowers
hist(sunflower_all_data$Abg_all_fresh) # fresh Aboveground biomass g 

sunflower_all_data$biomass_all= sunflower_all_data$Abg_all_dry+sunflower_all_data$dry_roots
hist(sunflower_all_data$biomass_all) # Total dry  biomass 

########### calculation of the cost of transpiration Ecost gdry roots/liter water ###################
sunflower_all_data$Ecost = sunflower_all_data$dry_roots/sunflower_all_data$transpiration
hist(sunflower_all_data$Ecost)

########### calculation of Ecost considering water stored in the plant gdry roots/liter water  ###################
sunflower_all_data$water_stored = sunflower_all_data$Abg_all_fresh-sunflower_all_data$Abg_all_dry
hist(sunflower_all_data$water_stored)

sunflower_all_data$Ecost_bis = sunflower_all_data$dry_roots/(sunflower_all_data$transpiration + sunflower_all_data$water_stored)
hist(sunflower_all_data$Ecost_bis)

plot(sunflower_all_data$Ecost_bis, sunflower_all_data$Ecost)

################### calculate Ncost g Roots / gN leaves ######################## 

sunflower_all_data$N_total_leaves = sunflower_all_data$Nmass * sunflower_all_data$dry_leaves_petioles
hist(sunflower_all_data$N_total_leaves)
hist(sunflower_all_data$Nmass)


sunflower_all_data$Ncost = sunflower_all_data$dry_roots/sunflower_all_data$N_total_leaves
hist(sunflower_all_data$Ncost)

sunflower_all_data$betaplant = sunflower_all_data$Ncost/sunflower_all_data$Ecost
hist(sunflower_all_data$betaplant)

plot(log(sunflower_all_data$betaplant), log(sunflower_all_data$betaleaf))

########## calculate the standardized cumulative drought severity DSI  ########
sunflower_all_data$DSI= sunflower_all_data$Sum_DS_daily/sunflower_all_data$DAD
hist(sunflower_all_data$DSI)
write.csv(sunflower_all_data, "../output/sunflower_all.csv", row.names = FALSE)
names(sunflower_all_data)

#### clean the final output file  
sunflower_all_data$date.x = NULL 
sunflower_all_data$leaves_disk_mg_noChl = NULL 
sunflower_all_data$nb_disks_scanned = NULL 
sunflower_all_data$nb_disk_SLA = NULL 
sunflower_all_data$leaves_disk_area = NULL 
sunflower_all_data$leaf_area_for_SLA = NULL 
sunflower_all_data$nb_disk_SLA = NULL 
sunflower_all_data$aci_RMSE = NULL 
sunflower_all_data$tmp = NULL
sunflower_all_data$VPD = NULL
sunflower_all_data$z = NULL
sunflower_all_data$patm = NULL
sunflower_all_data$nstar = NULL
sunflower_all_data$gammastar_pa = NULL
sunflower_all_data$Km = NULL
sunflower_all_data$D = NULL
sunflower_all_data$ca = NULL
sunflower_all_data$beta_num = NULL
sunflower_all_data$beta_denom = NULL

names(sunflower_all_data)

write.csv(sunflower_all_data, "../output/sunflower_clean.csv", row.names = FALSE)



