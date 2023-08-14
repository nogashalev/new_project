#######################################################################
# Using random effect as an underreporting factor to scale incidence  #
# Author: Noga Shalev                                                 #
# Date: 2023_8_9                                                      #
# Contents: Part A: SETUP                                             #
#           Part B: MODEL                                             #
#           Part C: PLOTS                                             #
#######################################################################

# A. Setup

#clear environment
rm(list = ls())

#load libraries
library("data.table")
library("lme4")
library("ggplot2")
library("dplyr")

#load functions
source("/ihme/cc_resources/libraries/current/r/get_population.R")
source("/ihme/cc_resources/libraries/current/r/get_location_metadata.R")
source("/ihme/cc_resources/libraries/current/r/get_model_results.R")
source("/ihme/cc_resources/libraries/current/r/get_outputs.R")

#load objects
release_id <- 9
acause    <- "whooping"
cause_id  <- 339
me_id <- 1424
years <-1940:2023

#load data

#location file
locations <- get_location_metadata(release_id = 9, location_set_id=22)[, .(location_name, region_name, super_region_name, location_id, ihme_loc_id, level)]

#load regression file from pertussis 01 script, most recent best version
regress <- fread("/snfs1/WORK/12_bundle/whooping/00_documentation/models/02.12.21/model_inputs/incidence_regression_input.csv")

#load predictions file for same date
predictions <- fread("/snfs1/WORK/12_bundle/whooping/00_documentation/models/02.12.21/01_case_predictions_from_model.csv")

#load model results
results <- get_outputs("cause", cause_id = cause_id,
                       release_id = release_id,
                       compare_version_id = 7904,
                       measure_id = c(6), #incidence
                       metric = 1, #count
                       year_id = years,
                       location_id = unique(locations$location_id),
                       sex_id = 3,
                       age_group_id = 22)


# B. Model
#original model call
me_model <- lmer(ln_inc ~ ln_unvacc + (1 | ihme_loc_id), data=regress)
summary(me_model)
residuals <- resid(me_model)

plot(residuals ~ fitted(me_model), main="Residuals vs. Fitted", cex = 0.1, pch = 20)
abline(h=0, col="red")

qqnorm(residuals, cex = 0.1, pch = 20)
qqline(residuals, col="red")

#plotting against the predictor resulted in mistmatched length of residuals versus predictor... 
# plot(residuals ~ ln_unvacc, data=regress, main="Residuals vs. Predictor", cex = 0.1, pch = 20)
# abline(h=0, col="red")

# Make dataframe of REs
ranefs <- data.frame(ranef(me_model))
setDT(ranefs)

ranefs <- ranefs[, c("term", "grpvar") := NULL]
setnames(ranefs, "grp", "ihme_loc_id")

# Merge random effects with locations
ranefs <- locations[ranefs, on = .(ihme_loc_id = ihme_loc_id)]

# plot ranefs by location
ranefs %>% ggplot(mapping = aes(x = reorder(super_region_name, condval), y = condval)) +
  geom_boxplot()+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none")+
  labs(title = "Random effect by super-region", 
       y = "Random effect", 
       x = "")

#create subset of high income locations
high_inc <- ranefs[super_region_name == "High-income", ]
col_5 <- viridis::viridis(5, direction = -1)
high_inc %>% ggplot(mapping = aes(reorder(x = location_name, condval), y = condval)) +
  geom_point(size = 1, aes(color = region_name))+
  scale_color_manual(values = col_5)+
  coord_flip()+
  theme_bw()+
  labs(title = "Random effects for high-income locations", 
       y = "Random effect", 
       x = "")

# Calculate mean of random effect for high income region 
mean_high_inc_rf <- mean(ranefs$condval[ranefs$super_region_name=="High-income"])

# merge ranefs with regress file
merged <- ranefs[regress, on = .(ihme_loc_id = ihme_loc_id, location_id = location_id)]
# there are 241 rows of missing RE estimates these are in 8 locations:
# Antigua and Barbuda, Saint Kitts and Nevis, United States Virgin Islands, Cook Islands, 
# Nauru, Niue, Palau, Tuvalu   
#remove NA rows for now
merged <- merged[!is.na(condval), ]

#also remove year id rows prior to 1980
merged <- merged[year_id >= 1980, ]

#Adjust incidence using formula: incidence_adjusted/incidence_unadjusted = exp(mean_RE_high - REc)
merged[, adjusted_incidence := (exp(mean_high_inc_rf - condval))*incidence]

# Plot adjusted vs unadjusted incidence
merged %>% ggplot(mapping = aes(x = incidence, y = adjusted_incidence))+
  geom_point(size = 0.3)+
  theme_bw()+
  facet_wrap(~super_region_name)+
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.2)+
  labs(title = "Adjusted vs unadjusted pertussis incidence, using random effect scalar", 
       x = "Unadjusted", 
       y = "Adjusted")

# for high income region substitute unadjusted incidence as incidence, ie no correction
merged[, final_incidence := ifelse(super_region_name=="High-income", incidence, adjusted_incidence)]

merged %>% ggplot(mapping = aes(x = incidence, y = final_incidence))+
  geom_point(size = 0.3)+
  theme_bw()+
  facet_wrap(~super_region_name)+
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.2)+
  labs(title = "Adjusted vs unadjusted pertussis incidence, using random effect scalar", 
       x = "Unadjusted", 
       y = "Adjusted")

# merged$year_id <- as.factor(merged$year_id)
# col_40 <- viridis::viridis(40, direction = -1)
# 
# by_year <- merged %>% ggplot(mapping = aes(x = incidence, y = final_incidence))+
#   geom_point(size = 0.3, aes(color = year_id))+
#   scale_color_manual(values = col_40)+
#   theme_bw()+
#   facet_wrap(~super_region_name)+
#   geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.2)+
#   theme(legend.position = "none")

by_location <- merged %>% ggplot(mapping = aes(x = incidence, y = final_incidence))+
  geom_point(size = 0.3, aes(color = location_name))+
  # scale_color_manual(values = col_40)+
  theme_bw()+
  facet_wrap(~super_region_name)+
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.2)+
  theme(legend.position = "none")

plotly::ggplotly(by_year)
plotly::ggplotly(by_location)

#because of the large spread accross high income REs create new version using only DR locations
stars <- fread("/snfs1/WORK/03_cod/01_database/02_programs_nocode/source_metadata/smp/gbd2020_best/stars_by_iso3_time_window.csv")
setDT(stars)

#keep star rating used in COD models
full_time <- stars[time_window == "full_time_series", ][, c("ihme_loc_id", "location_level", "location_name") := NULL]
full_time <- locations[full_time, on = .(location_id = location_id)]

full_time <- full_time[stars >=4 & level == 3, ] #these are all in high income locations
star_locs <- unique(full_time$location_name) #29! [nb there were 35 in high income, so 7 not stars!]

#identify which locations are not modelled
high_inc_locs <- unique(merged$location_name[merged$super_region_name=="High-income"])
high_not_star <- setdiff(high_inc_locs, star_locs)
# [1] "Brunei Darussalam" "Republic of Korea" "Andorra"          
# [4] "Cyprus"            "Monaco"            "San Marino"

# create a flag for starred location in the merged folder
merged[, star := factor(ifelse(location_name %in% star_locs, 1, 0))]

by_star <- merged[super_region_name == "High-income"] %>% ggplot(mapping = aes(reorder(x = location_name, condval), y = condval)) +
  geom_point(size = 1, aes(color = star))+
  # scale_color_manual(values = col_5)+
  coord_flip()+
  theme_bw()

plotly::ggplotly(by_star)

#create a new RE mean using starred locations only!
star_mean <- mean(merged$condval[merged$super_region_name == "High-income" & merged$star == 1])

#adjust incidence by star RE mean
merged[, adjusted_incidence_2 := (exp(star_mean - condval)) * incidence]
#unadjust starred regions
merged[, final_incidence_2 := ifelse(star == 1, incidence, adjusted_incidence_2)]

by_loc_2 <- merged %>% ggplot(mapping = aes(x = incidence, y = final_incidence_2))+
  geom_point(size = 0.3, aes(color = location_name))+
  # scale_color_manual(values = col_40)+
  theme_bw()+
  facet_wrap(~super_region_name)+
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.2)+
  theme(legend.position = "none")

plotly::ggplotly(by_loc_2)

#look at values of adjustment factor
merged[, `:=` (re_adjust_1 = exp(mean_high_inc_rf - condval), re_adjust_2 = exp(star_mean - condval))]
hist(merged$re_adjust_1, breaks = 100)
hist(merged$re_adjust_2, breaks = 100)

#let's examine where there are outliers 2SD beyond the mean
outlier_re <- merged[re_adjust_2 > mean(re_adjust_2) + 2* sd(re_adjust_2)] 
unique(outlier_re$location_name)# just Egypt!

#Create different RE scalars by decade
merged <- merged[year_id >= 1980, ]
merged[, decade := floor(year_id / 10) * 10]
by_decade <- split(merged, merged$decade)

#run lme model for each decade: write a for loop here!!!
me_model_80 <- lmer(ln_inc ~ ln_unvacc + (1 | ihme_loc_id), data=by_decade[[1]])
me_model_90 <- lmer(ln_inc ~ ln_unvacc + (1 | ihme_loc_id), data=by_decade[[2]])
me_model_00 <- lmer(ln_inc ~ ln_unvacc + (1 | ihme_loc_id), data=by_decade[[3]])
me_model_10 <- lmer(ln_inc ~ ln_unvacc + (1 | ihme_loc_id), data=by_decade[[4]])

time_models <- c(me_model_80,me_model_90, me_model_00, me_model_10)

# ranefs_ <- data.frame()
# for (i in time_models) {
#   ranefs_[i] <- data.frame(unlist(ranef(i)))
# }

# for (i in 1:4) {
#   ranefs_[i] <- ranef(unlist(time_models[i]))
# }

ranefs_80 <- data.frame(ranef(me_model_80))
ranefs_90 <- data.frame(ranef(me_model_90))
ranefs_00 <- data.frame(ranef(me_model_00))
ranefs_10 <- data.frame(ranef(me_model_10))

ranefs_80$decade <- "1980"
ranefs_90$decade <- "1990"
ranefs_00$decade <- "2000"
ranefs_10$decade <- "2010"

time_ranefs <- rbind(ranefs_80, ranefs_90, ranefs_00, ranefs_10)
setDT(time_ranefs)
time_ranefs <- time_ranefs[, .(grp, condval, decade)]

time_ranefs <- locations[time_ranefs, on = .(ihme_loc_id = grp)]

# compare ranefs by decade

time_ranefs %>% ggplot(mapping = aes(x = decade, y = condval))+
  geom_boxplot()+
  facet_wrap(~super_region_name)
#create mean starred ranefs by decade

time_ranefs_star <- time_ranefs[location_name %in% star_locs,][, .(mean_ranef = mean(condval)), by = decade]
merged$decade <- as.character(merged$decade)
merged <- time_ranefs_star[merged, on = .(decade = decade)]

merged[, re_adjust_time := exp(mean_ranef- condval)]
merged[, time_incidence := re_adjust_time * incidence]
merged[, time_incidence := ifelse(location_name %in% star_locs, incidence, time_incidence)]

merged %>% ggplot(mapping = aes(x = incidence, y = time_incidence))+
  geom_point(size = 0.3)+
  theme_bw()+
  facet_wrap(~super_region_name)+
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.2)+
  labs(title = "Adjusted vs unadjusted pertussis incidence, using time-varying random effect scalar", 
       x = "Unadjusted", 
       y = "Adjusted")

merged %>% ggplot(mapping = aes(x = final_incidence_2, y = time_incidence))+
  geom_point(size = 0.3)+
  theme_bw()+
  facet_wrap(~super_region_name)+
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.2)+
  labs(title = "Two adjustments of pertussis incidence, using random effect scalar", 
       x = "No time adjustment", 
       y = "Time adjustment")

#merge with modelled results to compare
results <- results[!is.na(val), .(location_id, year_id, val)]

results <- results[location_id %in% merged$location_id, ] #there are only 9 modelled years in the results... 

merged_results <- merged[results, on = .(location_id = location_id, year_id = year_id)]

merged_results <- merged_results[!is.na(location_name), ]

output_comp <- merged_results %>% ggplot(mapping = aes(x = val, y = time_incidence))+
  geom_point(aes(color = location_name))+
  theme_bw()+
  theme(legend.position = "none")
plotly::ggplotly(output_comp)

incidence_comp <- merged_results %>% ggplot(mapping = aes(x = val, y = incidence))+
  geom_point(aes(color = location_name))+
  theme_bw()+
  theme(legend.position = "none")

incidence_comp
plotly::ggplotly(output_comp)
#OK that comparison is flawed... but confirmed that india is modelled as millions of cases! 

#Merge with model predictions, should do a rowise mean but just chose the first draw column for now
predictions <- predictions[, 1:3]
setnames(predictions, "case_draw_0", "predicted_incidence") 

predictions <- predictions[location_id %in% merged$location_id]

merged_predictions <- predictions[merged, on = .(location_id = location_id, year_id = year_id)]

predictions_comp <- merged_predictions %>% ggplot(mapping = aes(x = predicted_incidence, y = time_incidence))+
  geom_point(size = 0.5, aes(color = location_name))+
  theme_bw()+
  theme(legend.position = "none")+
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.2)
# geom_point(y = merged$incidence)

predictions_comp
plotly::ggplotly(predictions_comp)

#use get model results: but these are incidence rates!
# pertussis <- get_model_results('epi', me_id, release_id = release_id)

#mimic prediction using the formula applied to draws without the draws using intercept beta2 a beta 3! 
# draws[, (draw_cols) := lapply(draw_nums, function(x) { exp( betas[1, x] +
# ( betas[2, x] * ln_unvacc ) +
# ( betas[3, x] ) ) * ( pop / 100000 ) } )]

merged_predictions[is.na(pop), pop := ukpop]

intercept <- as.numeric(fixef(me_model)["(Intercept)"])
beta1 <- as.numeric(fixef(me_model)["ln_unvacc"])

merged_predictions[, formula_incidence := exp(intercept + beta1*ln_unvacc + max(ranefs$condval)) * pop/100000]

plot(merged_predictions$formula_incidence, merged_predictions$predicted_incidence, cex = 0.5, pch = 20)
plot(merged_predictions$formula_incidence, merged_predictions$time_incidence, cex = 0.5, pch = 20)

#prepare for STGPR
merged_predictions[, `:=` (nid = 83133, measure = "continuous", age_group_id = 22, sex_id = 3, is_outlier = 0)]
merged_predictions[, `:=` (measure_id = 6, sample_size = pop, variance = 0.001, me_name = "Pertussis incidence")]

#need to create separate files now because for STGPR the target variable has to be labelled val
time_re_target <- copy(merged_predictions)
time_re_target <- setnames(time_re_target, "time_incidence", "val")

predict_by_formula <- copy(merged_predictions)
predict_by_formula <- setnames(predict_by_formula, "formula_incidence", "val")

# need to limit to desired columns taking this from the tetanus PAB model!
columns_keep <- colnames(tt)
#subset columns to keep
time_re_target <- time_re_target[, ..columns_keep, with = F]
predict_by_formula <- predict_by_formula[, ..columns_keep, with = F]

fwrite(time_re_target, "/snfs1/WORK/12_bundle/whooping/27633/pertussis_model_with_time_re.csv")
fwrite(predict_by_formula, "/snfs1/WORK/12_bundle/whooping/27633/pertussis_model_with_formula_re.csv")

#!!! NB path to data has no quotes! Didnt work until I a) put in pop as sample size, b) limited the file to the above columns. 
# I did all of this at the same time so it's unclear to me what fixed it! 

#vetting shows increase in numbers in high income lcoations in mid 1990s
merged_predictions$year_id <- as.factor(merged_predictions$year_id)

merged_predictions[super_region_name == "High-income"] %>% ggplot(mapping = aes(x = year_id, y = incidence))+
  geom_point(aes(color = location_name))+
  facet_wrap(~super_region_name)

##launch STGPR
source("/ihme/code/st_gpr/central/src/stgpr/api/public.R")
status <- get_model_status(206588)


#practice STGPR prior to launch requested new MeID 8_8_2023
stub <- fread("/mnt/team/cc/pub/stgpr_walkthrough/stgpr_stub_data.csv")
config <- fread("/mnt/team/cc/pub/stgpr_walkthrough/stgpr_stub_config.csv")
file.edit("/mnt/team/cc/pub/stgpr_walkthrough/stgpr_example_workflow_script.R")












