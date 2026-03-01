#########################################
# packages necessary
R.version.string
#########################################
library(ggplot2)
library(car)
library(emmeans)
library(nlme)
library(piecewiseSEM)
library(multcomp)
library(r2glmm)
library(viridis)
############## load data #################################################################### 
stat_all <- read.csv("../output/sunflower.csv")
names(stat_all)

stat_t4 = subset(stat_all, time =="t4")
table(stat_t4$time)
table(stat_t4$Cumul_N)
betaleaf_lm_t4 <- lm(log(betaleaf)~ Cumul_N*Sdt_DD,
                       data = stat_t4)
Anova(betaleaf_lm_t4) 

stat_t3 = subset(stat_all, time =="t3")
table(stat_t3$time)
table(stat_t3$Cumul_N)
betaleaf_lm_t3 <- lm(log(betaleaf)~ Cumul_N*Sdt_DD,
                     data = stat_t3)
Anova(betaleaf_lm_t3) 

stat_t2 = subset(stat_all, time =="t2")
table(stat_t2$time)
table(stat_t2$Cumul_N)
betaleaf_lm_t2 <- lm(log(betaleaf)~ Cumul_N*Sdt_DD,
                     data = stat_t2)
Anova(betaleaf_lm_t2) 

#############################  Analysis  #############################################
names(stat_t4)
plot(stat_t4$chi, stat_t4$betaleaf)

######################## Beta leaf ######################################################################################### 
hist(stat_t4$betaleaf)
hist(log(stat_t4$betaleaf))
nrow(stat_t4)
betaleaf_lm_t4 <- lm(log(betaleaf)~ Cumul_N*Sdt_DD,
                    data = stat_t4)
Anova(betaleaf_lm_t4) ## drought and N effects 

plot(betaleaf_lm_t4)
shapiro.test(residuals(betaleaf_lm_t4)) # 0.64
outlierTest(betaleaf_lm_t4)  ### No outliers
AIC(betaleaf_lm_t4) #14.99
vif(betaleaf_lm_t4) # no colinerarity
plot(resid(betaleaf_lm_t4) ~ fitted(betaleaf_lm_t4)) ## good
r.squaredGLMM(betaleaf_lm_t4) ## marginal : 0.27, conditional: 0.27
summary(betaleaf_lm_t4) 
qqline(residuals(betaleaf_lm_t4))
densityPlot(residuals(betaleaf_lm_t4))

residuals <- resid(betaleaf_lm_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad


############################ set this for all figures in this script #####################################
cumul_N_levels <- c(540, 1620)
# prediction grid
Sdt_DD_range_t4 <- seq(from = min(stat_t4$Sdt_DD, na.rm = TRUE), to = max(stat_t4$Sdt_DD, na.rm = TRUE),
  length.out = 400)
#######################################################################################################

newdat_t4_betaleaf <- expand.grid(Sdt_DD  = Sdt_DD_range_t4, Cumul_N = cumul_N_levels)

# Predictions + 95% confidence interval 
pred_log_betaleaf <- predict(betaleaf_lm_t4, newdata = newdat_t4_betaleaf, interval = "confidence")
newdat_t4_betaleaf$fit <- exp(pred_log_betaleaf[, "fit"])
newdat_t4_betaleaf$lwr <- exp(pred_log_betaleaf[, "lwr"])
newdat_t4_betaleaf$upr <- exp(pred_log_betaleaf[, "upr"])
names(newdat_t4_betaleaf)

#################################################################################################
slopes_betaleaf<-  emtrends(betaleaf_lm_t4, ~ Cumul_N, var = "Sdt_DD",
         at = list(Cumul_N = c(540, 1620)))
summary(slopes_betaleaf, infer = c(TRUE, TRUE))
################## Figures : betaleaf ################################################

# colors
pal_line <- viridis(100, option = "viridis", direction = -1)
col_540  <- pal_line[10]
col_1620 <- pal_line[90]

betaleaf_plot_t4 <- ggplot(stat_t4, aes(x = Sdt_DD, y = betaleaf)) +
  
  geom_point(aes(shape = factor(drought), color = Cumul_N),
             size = 2.4, alpha = 0.60) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_viridis_c(option = "viridis", direction = -1, name = "N level", end = 0.95) +
  
  geom_ribbon(data = subset(newdat_t4_betaleaf, Cumul_N == 540),aes(x = Sdt_DD, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,fill = col_540, alpha = 0.10, colour = NA) +
  
  geom_ribbon(data = subset(newdat_t4_betaleaf, Cumul_N == 1620),aes(x = Sdt_DD, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,fill = col_1620, alpha = 0.10, colour = NA) +
  
  # fit lines
  geom_line(data = subset(newdat_t4_betaleaf, Cumul_N == 540),aes(x = Sdt_DD, y = fit),
    inherit.aes = FALSE,color = col_540, linewidth = 2.4, alpha = 0.95, linetype = "dashed") +
  
  geom_line(data = subset(newdat_t4_betaleaf, Cumul_N == 1620),aes(x = Sdt_DD, y = fit),
    inherit.aes = FALSE,color = col_1620, linewidth = 2.4, alpha = 0.95) +
  
  coord_cartesian(xlim = range(Sdt_DD_range_t4)) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    axis.title.y = element_text(size = 18, colour = "black"),
    axis.title.x = element_text(size = 18, colour = "black"),
    axis.text    = element_text(size = 14, colour = "black"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(linewidth = 0.8, colour = "black")) +
  
  ylab(expression(italic('β')['leaf'])) +
  xlab(expression(italic('DSI')))

betaleaf_plot_t4

############################################ Cost Nitrogen whole plant ##################################
hist(stat_t4$Ncost)
hist(log(stat_t4$Ncost))

Cost_nitrogen_lm_t4 <- lm(Ncost ~ Cumul_N*Sdt_DD,
                            data = stat_t4)
Anova(Cost_nitrogen_lm_t4)# N and interaction effects 

shapiro.test(residuals(Cost_nitrogen_lm_t4)) # 0.22
outlierTest(Cost_nitrogen_lm_t4)  ### 
AIC(Cost_nitrogen_lm_t4) # 517.788
vif(Cost_nitrogen_lm_t4) # no colinerarity
plot(resid(Cost_nitrogen_lm_t4) ~ fitted(Cost_nitrogen_lm_t4)) ## good
r.squaredGLMM(Cost_nitrogen_lm_t4) ## marginal : 0.238, conditional: 0.24
summary(Cost_nitrogen_lm_t4) 
qqnorm(residuals(Cost_nitrogen_lm_t4))
qqline(residuals(Cost_nitrogen_lm_t4))
densityPlot(residuals(Cost_nitrogen_lm_t4))

residuals <- resid(Cost_nitrogen_lm_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

newdat_t4_Ncost <- expand.grid(Sdt_DD  = Sdt_DD_range_t4, Cumul_N = cumul_N_levels)

# Predictions + 95% confidence interval 
pred_Ncost <- predict(Cost_nitrogen_lm_t4,
                      newdata = newdat_t4_Ncost,
                      interval = "confidence")

newdat_t4_Ncost$fit <- pred_Ncost[, "fit"]
newdat_t4_Ncost$lwr <- pred_Ncost[, "lwr"]
newdat_t4_Ncost$upr <- pred_Ncost[, "upr"]

#################################################################################################
slopes_Ncost<-  emtrends(Cost_nitrogen_lm_t4, ~ Cumul_N, var = "Sdt_DD",
                             at = list(Cumul_N = c(540, 1620)))
summary(slopes_Ncost, infer = c(TRUE, TRUE))
################## Figures : Ncost ################################################

# colors
pal_line <- viridis(100, option = "viridis", direction = -1)
col_540  <- pal_line[10]
col_1620 <- pal_line[90]

Ncost_plot_t4 <- ggplot(stat_t4, aes(x = Sdt_DD, y = Ncost)) +
  
  geom_point(aes(shape = factor(drought), color = Cumul_N),
             size = 2.4, alpha = 0.60) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_viridis_c(option = "viridis", direction = -1, name = "N level", end = 0.95) +
  
  geom_ribbon(data = subset(newdat_t4_Ncost, Cumul_N == 540),aes(x = Sdt_DD, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,fill = col_540, alpha = 0.10, colour = NA) +
  
  geom_ribbon(data = subset(newdat_t4_Ncost, Cumul_N == 1620),aes(x = Sdt_DD, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,fill = col_1620, alpha = 0.10, colour = NA) +
  
  # fit lines
  geom_line(data = subset(newdat_t4_Ncost, Cumul_N == 540),aes(x = Sdt_DD, y = fit),
            inherit.aes = FALSE,color = col_540, linewidth = 2.4, alpha = 0.95) +
  
  geom_line(data = subset(newdat_t4_Ncost, Cumul_N == 1620),aes(x = Sdt_DD, y = fit),
            inherit.aes = FALSE,color = col_1620, linewidth = 2.4, alpha = 0.95, linetype = "dashed") +
  
  coord_cartesian(xlim = range(Sdt_DD_range_t4)) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    axis.title.y = element_text(size = 18, colour = "black"),
    axis.title.x = element_text(size = 18, colour = "black"),
    axis.text    = element_text(size = 14, colour = "black"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(linewidth = 0.8, colour = "black")) +
  
  ylab(expression(italic('b')['plant'])) +
  xlab(expression(italic('DSI')))

Ncost_plot_t4


######################## Ecost ############################################
hist(stat_t4$Ecost)
hist(log(stat_t4$Ecost))

Cost_trans_roots_lm_t4 <- lm(Ecost~ Cumul_N*Sdt_DD,
                               data = stat_t4)

Anova(Cost_trans_roots_lm_t4) # Nothing

summary(Cost_trans_roots_lm_t4)
shapiro.test(residuals(Cost_trans_roots_lm_t4)) 
outlierTest(Cost_trans_roots_lm_t4)  ### p = 0.11
AIC(Cost_trans_roots_lm_t4) # 94.13
vif(Cost_trans_roots_lm_t4) # no colinerarity
plot(resid(Cost_trans_roots_lm_t4) ~ fitted(Cost_trans_roots_lm_t4)) ## good
r.squaredGLMM(Cost_trans_roots_lm_t4) ## marginal : 0.05, conditional: 0.05
summary(Cost_trans_roots_lm_t4) 
qqnorm(residuals(Cost_trans_roots_lm_t4))
qqline(residuals(Cost_trans_roots_lm_t4))
densityPlot(residuals(Cost_trans_roots_lm_t4))

residuals <- resid(Cost_trans_roots_lm_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad


newdat_t4_Ecost <- expand.grid(Sdt_DD  = Sdt_DD_range_t4, Cumul_N = cumul_N_levels)

# Predictions + 95% confidence interval 
pred_Ecost <- predict(Cost_trans_roots_lm_t4,
                      newdata = newdat_t4_Ecost,
                      interval = "confidence")

newdat_t4_Ecost$fit <- pred_Ecost[, "fit"]
newdat_t4_Ecost$lwr <- pred_Ecost[, "lwr"]
newdat_t4_Ecost$upr <- pred_Ecost[, "upr"]

#################################################################################################

################## Figures : Ecost ################################################

# colors
pal_line <- viridis(100, option = "viridis", direction = -1)
col_540  <- pal_line[10]
col_1620 <- pal_line[90]

Ecost_plot_t4 <- ggplot(stat_t4, aes(x = Sdt_DD, y = Ecost)) +
  
  geom_point(aes(shape = factor(drought), color = Cumul_N),
             size = 2.4, alpha = 0.60) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_viridis_c(option = "viridis", direction = -1, name = "N level", end = 0.95) +
  
  geom_ribbon(data = subset(newdat_t4_Ecost, Cumul_N == 540),aes(x = Sdt_DD, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,fill = col_540, alpha = 0.10, colour = NA) +
  
  geom_ribbon(data = subset(newdat_t4_Ecost, Cumul_N == 1620),aes(x = Sdt_DD, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,fill = col_1620, alpha = 0.10, colour = NA) +
  
  # fit lines
  geom_line(data = subset(newdat_t4_Ecost, Cumul_N == 540),aes(x = Sdt_DD, y = fit),
            inherit.aes = FALSE,color = col_540, linewidth = 2.4, alpha = 0.95, linetype = "dashed") +
  
  geom_line(data = subset(newdat_t4_Ecost, Cumul_N == 1620),aes(x = Sdt_DD, y = fit),
            inherit.aes = FALSE,color = col_1620, linewidth = 2.4, alpha = 0.95, linetype = "dashed") +
  
  coord_cartesian(xlim = range(Sdt_DD_range_t4)) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    axis.title.y = element_text(size = 18, colour = "black"),
    axis.title.x = element_text(size = 18, colour = "black"),
    axis.text    = element_text(size = 14, colour = "black"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(linewidth = 0.8, colour = "black")) +
  
  ylab(expression(italic('a')['plant'])) +
  xlab(expression(italic('DSI')))

Ecost_plot_t4

######################################################################################################################
############################## betaplant ####################################################
names(stat_t4)
hist(log(stat_t4$betaplant))
betaplant_lm_t4 <- lm(log(betaplant)~ Cumul_N*Sdt_DD,
                    data = stat_t4)
Anova(betaplant_lm_t4) # N effects and interaction. plus mariginal drought effects 
summary(betaplant_lm_t4) 
shapiro.test(residuals(betaplant_lm_t4)) 
outlierTest(betaplant_lm_t4)  
AIC(betaplant_lm_t4) 
vif(betaplant_lm_t4) # no colinerarity
plot(resid(betaplant_lm_t4) ~ fitted(betaplant_lm_t4)) ## good
r.squaredGLMM(betaplant_lm_t4) 
summary(betaplant_lm_t4) 
qqnorm(residuals(betaplant_lm_t4))
qqline(residuals(betaplant_lm_t4))
densityPlot(residuals(betaplant_lm_t4))

residuals <- resid(betaplant_lm_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad
summary(betaplant_lm_t4)
range(stat_t4$betaplant)

newdat_t4_betaplant<- expand.grid(Sdt_DD  = Sdt_DD_range_t4, Cumul_N = cumul_N_levels)

# Predictions + 95% confidence interval 
pred_log_betaplant <- predict(betaplant_lm_t4, newdata = newdat_t4_betaplant, interval = "confidence")
newdat_t4_betaplant$fit <- exp(pred_log_betaplant[, "fit"])
newdat_t4_betaplant$lwr <- exp(pred_log_betaplant[, "lwr"])
newdat_t4_betaplant$upr <- exp(pred_log_betaplant[, "upr"])
names(newdat_t4_betaplant)

#################################################################################################
#################################################################################################
slopes_betaplant<-  emtrends(betaplant_lm_t4, ~ Cumul_N, var = "Sdt_DD",
                            at = list(Cumul_N = c(540, 1620)))
summary(slopes_betaplant, infer = c(TRUE, TRUE))
################## Figures : betaplant ################################################

# colors
pal_line <- viridis(100, option = "viridis", direction = -1)
col_540  <- pal_line[10]
col_1620 <- pal_line[90]

betaplant_plot_t4 <- ggplot(stat_t4, aes(x = Sdt_DD, y = betaplant)) +
  
  geom_point(aes(shape = factor(drought), color = Cumul_N),
             size = 2.4, alpha = 0.60) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_viridis_c(option = "viridis", direction = -1, name = "N level", end = 0.95) +
  
  geom_ribbon(data = subset(newdat_t4_betaplant, Cumul_N == 540),aes(x = Sdt_DD, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,fill = col_540, alpha = 0.10, colour = NA) +
  
  geom_ribbon(data = subset(newdat_t4_betaplant, Cumul_N == 1620),aes(x = Sdt_DD, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,fill = col_1620, alpha = 0.10, colour = NA) +
  
  # fit lines
  geom_line(data = subset(newdat_t4_betaplant, Cumul_N == 540),aes(x = Sdt_DD, y = fit),
            inherit.aes = FALSE,color = col_540, linewidth = 2.4, alpha = 0.95) +
  
  geom_line(data = subset(newdat_t4_betaplant, Cumul_N == 1620),aes(x = Sdt_DD, y = fit),
            inherit.aes = FALSE,color = col_1620, linewidth = 2.4, alpha = 0.95, linetype = "dashed") +
  
  coord_cartesian(xlim = range(Sdt_DD_range_t4)) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    axis.title.y = element_text(size = 18, colour = "black"),
    axis.title.x = element_text(size = 18, colour = "black"),
    axis.text    = element_text(size = 14, colour = "black"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(linewidth = 0.8, colour = "black")) +
  
  ylab(expression(italic('β')['plant'])) +
  xlab(expression(italic('DSI')))

betaplant_plot_t4

####################### Structural equation model ################################
plot(stat_t4$betaleaf, stat_t4$vcmax_tleaf)

sunflower_beta <- psem(
  lm1 <- lm(Ncost~ Sdt_DD + Cumul_N , data = stat_t4),
  lm2 <- lm(Ecost ~Sdt_DD, data = stat_t4),
  lm3 <- lm(betaplant~Ncost + Ecost, data = stat_t4),
  lm4 <- lm(betaleaf~ betaplant+ Sdt_DD + Cumul_N , data = stat_t4),
  lm5 <- lm(vcmax_tleaf~betaleaf + betaplant, data = stat_t4),
  lm6 <- lm(Narea~vcmax_tleaf, data = stat_t4))

summary(sunflower_beta)
plot(sunflower_beta)



