
######################## Export tables #####################################
library(emmeans)
library(multcomp)
library(dplyr)
library(officer)
library(flextable)


N_levels <- c(540, 1620)
var_name <- "DSI"   

betaplant_lm_t4       <- lm(log(betaplant) ~ Cumul_N * DSI, data = stat_t4)
Cost_nitrogen_lm_t4   <- lm(Ncost ~ Cumul_N * DSI, data = stat_t4)
Ecost_lm_t4           <- lm(Ecost ~ Cumul_N * DSI, data = stat_t4)     
betaleaf_lm_t4        <- lm(log(betaleaf) ~ Cumul_N * DSI, data = stat_t4)

extract_anova_fp <- function(model, response_label) {
  tab <- as.data.frame(car::Anova(model, type = 2))
  
  out <- tab %>%
    dplyr::mutate(
      Response  = response_label,
      Predictor = rownames(tab),
      F_value   = round(`F value`, 3),
      P_value   = round(`Pr(>F)`, 4),
      Signif    = dplyr::case_when(
        `Pr(>F)` < 0.001 ~ "***",
        `Pr(>F)` < 0.01  ~ "**",
        `Pr(>F)` < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::select(Response, Predictor, F_value, P_value, Signif) %>%
    dplyr::filter(Predictor != "(Intercept)")
  
  out
}

anova_all <- bind_rows(
  extract_anova_fp(betaleaf_lm_t4,      "betaleaf (log)"),
  extract_anova_fp(Cost_nitrogen_lm_t4, "Ncost"),
  extract_anova_fp(Ecost_lm_t4,         "Ecost"),
  extract_anova_fp(betaplant_lm_t4,     "betaplant (log)")
)

doc <- read_docx()
doc <- body_add_par(doc, "Type II ANOVA (F and P values)", style = "heading 1")
doc <- body_add_flextable(doc, flextable(anova_all))
print(doc, target = "../tables/ANOVA_F_P_tables.docx")


#extract slope tables + comparison 
make_slope_outputs <- function(model, response_label, N_levels, var_name = "DSI") {
  
  tr <- emtrends(
    model,
    ~ Cumul_N,
    var = var_name,
    at = list(Cumul_N = N_levels)
  )
  
  slopes_tbl <- as.data.frame(summary(tr, infer = TRUE)) %>%
    mutate(
      Response  = response_label,
      Slope_SE  = sprintf("%.3f \u00B1 %.3f", !!sym(paste0(var_name, ".trend")), SE),
      P_value   = round(p.value, 4)
    ) %>%
    dplyr::select(Response, Cumul_N, Slope_SE, P_value)
  
  cld_tbl <- as.data.frame(cld(tr, Letters = letters)) %>%
    mutate(
      Response = response_label,
      Slope_SE = sprintf("%.3f \u00B1 %.3f%s",
                         !!sym(paste0(var_name, ".trend")), SE, .group)
    ) %>%
    dplyr::select(Response, Cumul_N, Slope_SE, .group)

  comp_tbl <- as.data.frame(pairs(tr)) %>%
    mutate(
      Response = response_label,
      P_value  = round(p.value, 4)
    ) %>%
    dplyr::select(Response, contrast, estimate, SE, df, t.ratio, P_value)
  
  list(
    slopes_tbl = slopes_tbl,
    cld_tbl    = cld_tbl,
    comp_tbl   = comp_tbl
  )
}

out_betaleaf   <- make_slope_outputs(betaleaf_lm_t4,      "betaleaf (log)",   N_levels, var_name)
out_Ncost      <- make_slope_outputs(Cost_nitrogen_lm_t4, "Ncost",            N_levels, var_name)
out_Ecost      <- make_slope_outputs(Ecost_lm_t4,         "Ecost",            N_levels, var_name)
out_betaplant  <- make_slope_outputs(betaplant_lm_t4,     "betaplant (log)",  N_levels, var_name)

#Combine tables
slopes_all <- bind_rows(
  out_betaleaf$slopes_tbl,
  out_Ncost$slopes_tbl,
  out_Ecost$slopes_tbl,
  out_betaplant$slopes_tbl
)

cld_all <- bind_rows(
  out_betaleaf$cld_tbl,
  out_Ncost$cld_tbl,
  out_Ecost$cld_tbl,
  out_betaplant$cld_tbl
)

comparisons_all <- bind_rows(
  out_betaleaf$comp_tbl,
  out_Ncost$comp_tbl,
  out_Ecost$comp_tbl,
  out_betaplant$comp_tbl
)

print(slopes_all)
print(cld_all)
print(comparisons_all)

#Export to Word 
doc <- read_docx()

doc <- body_add_par(doc, "Slopes at representative N levels (test slope != 0)", style = "heading 1")
doc <- body_add_flextable(doc, flextable(slopes_all))

doc <- body_add_par(doc, "Slopes with grouping letters (slope differences between N levels)", style = "heading 1")
doc <- body_add_flextable(doc, flextable(cld_all))

doc <- body_add_par(doc, "Pairwise comparison of slopes (540 vs 1620)", style = "heading 1")
doc <- body_add_flextable(doc, flextable(comparisons_all))

print(doc, target = "../tables/Slope_tables_time4.docx")