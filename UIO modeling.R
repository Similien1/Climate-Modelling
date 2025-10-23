
# ============================================================================
# DIARRHEA MODELING PIPELINE - CLEANED & ORGANIZED
# ============================================================================

# 1. SETUP -------------------------------------------------------------------
# (Optional) setwd() and package installation
# setwd("D:/your/path")

# 2. LOAD PACKAGES ----------------------------------------------------------
#if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, sf, tidyr, ggplot2, lubridate, ggtext, INLA, spdep,
  readxl, rio, here, stringr, INLAOutputs, viridis, patchwork,
  rnaturalearthdata, plotly, gganimate, leaflet, rnaturalearth, usdm,
  RColorBrewer, gifski, forecast, shiny, car, corrplot
)

# ------------------ Additional cleaning and improvements --------------------
# The full cleaned and organized script will be inserted here in the next step.
# ============================================================================
# DIARRHEA MODELING PIPELINE — CLEANED & STRUCTURED
# ============================================================================

# 1. SETUP -------------------------------------------------------------------
# (Optional) setwd("D:/your/path")


# 3. LOAD & CLEAN SPATIAL DATA ----------------------------------------------
Rwa_data <- st_read("shapefile/rwa_adm3_2006_NISR_WGS1984_20181002.shp") %>%
  mutate(
    Sec_ID = row_number(),
    ADM3_EN = case_when(
      ADM3_EN == "Mageregere" ~ "Mageragere",
      ADM3_EN == "Shyrongi"   ~ "Shyorongi",
      ADM3_EN == "Ririma"     ~ "Rilima",
      TRUE ~ ADM3_EN
    )
  ) %>%
  st_make_valid()

# Build adjacency graph for INLA
nb_obj <- poly2nb(Rwa_data, queen = FALSE)
nb2INLA("Adj_Map.graph", nb_obj)
Rwa_adj <- file.path(getwd(), "Adj_Map.graph")

# Plot Rwanda sectors
ggplot(Rwa_data) +
  geom_sf(color = "blue", fill = "white") +
  coord_sf() + theme_bw()

# 4. LOAD CLINICAL & MALARIA DATA -------------------------------------------
U5_diarrhea <- read_excel("U_5_diarrhea_climate.xlsx")
Malaria_raw <- read_excel("malaria_modelling_Oslo.xlsx")

U5_diarrhea <- U5_diarrhea %>%
  mutate(
    Sector = str_to_sentence(str_replace(Sector, "\\s*\\([^\\)]+\\)", "")),
    District = str_to_sentence(District)
  )

# 5. CLEAN U5 & MALARIA DATA ------------------------------------------------
U5_diarrhea_clean <- U5_diarrhea %>%
  group_by(District, Year) %>%
  mutate(across(c(Population, U_5_Population, U_5_Ratio, Diarrhoea_cases),
                ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x))) %>%
  ungroup() %>%
  rename(district = District, sector = Sector) %>%
  mutate(across(c(district, sector), ~ str_to_lower(str_squish(.x))),
         Year = as.integer(Year), Month = as.integer(Month))

Malaria_sub <- Malaria_raw %>%
  select(District, Sector, Year, Month, `Malaria Cases`) %>%
  mutate(
    district = str_to_lower(str_squish(District)),
    sector   = str_to_lower(str_squish(Sector)),
    Malaria_Cases = `Malaria Cases`
  ) %>%
  select(district, sector, Year, Month, Malaria_Cases)

# 6. MERGE DISEASE DATASETS -------------------------------------------------
Merged <- left_join(U5_diarrhea_clean, Malaria_sub,
                    by = c("district", "sector", "Year", "Month")) %>%
  replace_na(list(Malaria_Cases = 0, Diarrhoea_cases = 0)) %>%
  mutate(
    district = str_to_sentence(district),
    sector   = str_to_sentence(sector)
  )

# Join with spatial polygons
Geo_data <- merge(Rwa_data, Merged,
                  by.x = c("ADM3_EN", "ADM2_EN"),
                  by.y = c("sector", "district"),
                  all.x = FALSE)

# 7. PREPARE FOR MODELLING --------------------------------------------------
dat <- Geo_data %>%
  arrange(Sec_ID, Year, Month) %>%
  mutate(
    ID              = Sec_ID,
    Malaria_Cases   = replace_na(Malaria_Cases, 0),
    Diarrhoea_cases = replace_na(Diarrhoea_cases, 0),
    ID.time         = rep(1:(10*12), times = n_distinct(Sec_ID)),
    ID.space.time   = row_number()
  ) %>%
  rename(
    Prec              = `Precipitation (ERA5-Land)`,
    Min_temperature   = `Min_ temperature`,
    Relative_humidity = `Relative humidity (ERA5-Land)`
  )

# 7.1 LAGGED COVARIATES ------------------------------------------------------
dat <- dat %>%
  arrange(Sec_ID, Year, Month) %>%
  group_by(Sec_ID) %>%
  mutate(
    Max_temp_lag1 = lag(Max_temperature, 1),
    Max_temp_lag2 = lag(Max_temperature, 2),
    Prec_lag1     = lag(Prec, 1),
    Prec_lag2     = lag(Prec, 2),
    Min_temp_lag1 = lag(Min_temperature, 1),
    Min_temp_lag2 = lag(Min_temperature, 2),
    RH_lag1       = lag(Relative_humidity, 1),
    RH_lag2       = lag(Relative_humidity, 2),
    Air_temp_lag1 = lag(Air_Temperature, 1),
    Air_temp_lag2 = lag(Air_Temperature, 2)
  ) %>%
  ungroup()

# 8. MULTICOLLINEARITY CHECK ------------------------------------------------
df_clim <- dat %>% st_drop_geometry() %>%
  select(Max_temperature, Prec, Min_temperature, Relative_humidity, Air_Temperature) %>%
  na.omit()

par(mar = c(0, 0, 2, 0))
corrplot(cor(df_clim), method = "color", addCoef.col = "black", tl.cex = 0.8)
title("Correlation Matrix of Climate Covariates")

# VIF diagnostics
df_num <- as.data.frame(lapply(df_clim, as.numeric))
vif_res <- usdm::vif(df_num)
print(vif_res)

# 2) Priors (PC priors commonly used for BYM2)
hyper_bym2 <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)),  # P(1/sd > 1) = 0.01
  phi  = list(prior = "pc",      param = c(0.5, 2/3))  # P(phi < 0.5) = 2/3
)


formula_bym <- as.integer(Diarrhoea_cases) ~ 1 +
  f(ID, model = "bym", graph = Rwa_adj, scale.model = TRUE, constr = TRUE) +
  f(ID.time, model = "rw1", constr = TRUE, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(ID.space.time, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

formula_bym_climate <- update(formula_bym,
                              . ~ . + scale(Max_temperature) + scale(Prec) +
                                scale(Min_temperature) + scale(Relative_humidity)
)

formula_bym_climate_lags <- update(formula_bym,
                                   . ~ . +
                                     scale(Max_temp_lag1) + scale(Max_temp_lag2) +
                                     scale(Prec_lag1) + scale(Prec_lag2) +
                                     scale(Min_temp_lag1) + scale(Min_temp_lag2) +
                                     scale(RH_lag1) + scale(RH_lag2)
)

formula_bym2_climate <- as.integer(Diarrhoea_cases) ~ 1 +
  scale(Max_temperature) + scale(Prec) +
  scale(Min_temperature) + scale(Relative_humidity) +
  f(ID, model = "bym2", graph = Rwa_adj, scale.model = TRUE, constr = TRUE,
    hyper = hyper_bym2) +
  f(ID.time, model = "rw1", constr = TRUE, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(ID.space.time, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

# 10. FIT MODELS -------------------------------------------------------------
fit_inla_model <- function(formula, data, offset_var) {
  inla(
    formula = formula,
    family  = "poisson",
    data    = data,
    offset  = log(offset_var / 1000),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    control.predictor = list(compute = TRUE)
  )
}

spatial_model           <- fit_inla_model(formula_bym, dat, dat$U_5_Population)
spatial_model_clima     <- fit_inla_model(formula_bym_climate, dat, dat$U_5_Population)
spatial_model_lags      <- fit_inla_model(formula_bym_climate_lags, dat, dat$U_5_Population)
spatial_model_bym2      <- fit_inla_model(formula_bym2_climate, dat, dat$U_5_Population)

# 10.1 Summaries
summary(spatial_model)
FixedEffects(spatial_model)
######&&&&&&&&&&&&&&&&&&
summary(spatial_model_clima)
FixedEffects(spatial_model_clima)
FixedEffects(spatial_model_lags)
FixedEffects(spatial_model_bym2) 


# 11. EVALUATE MODELS --------------------------------------------------------
get_model_metrics <- function(model, observed) {
  pred <- model$summary.fitted.values$mean
  waic <- model$waic$waic
  dic  <- model$dic$dic
  mlcpo <- -mean(log(model$cpo$cpo), na.rm = TRUE)
  rmse <- sqrt(mean((observed - pred)^2))
  mae  <- mean(abs(observed - pred))
  tibble(WAIC = waic, DIC = dic, MeanLogCPO = mlcpo, RMSE = rmse, MAE = mae)
}

models <- list(
  NoClimate    = spatial_model,
  Climate      = spatial_model_clima,
  Climate_Lags = spatial_model_lags,
  Climate_BYM2 = spatial_model_bym2
)

metrics_df <- bind_rows(lapply(models, get_model_metrics, dat$Diarrhoea_cases), .id = "Model")
print(metrics_df)




# 13. TEMPORAL EFFECTS ------------------------------------------------------
n_time <- length(spatial_model_clima$marginals.random$ID.time)

RR_temp <- sapply(spatial_model_clima$marginals.random$ID.time, function(m)
  inla.emarginal(function(x) exp(x), m))

CI_temp <- t(sapply(spatial_model_clima$marginals.random$ID.time, function(m)
  c(
    lo = inla.qmarginal(0.025, inla.tmarginal(exp, m)),
    hi = inla.qmarginal(0.975, inla.tmarginal(exp, m))
  )
))

temporal_df <- tibble(
  Time = seq_len(n_time),
  RR   = RR_temp,
  Low  = CI_temp[,"lo"],
  High = CI_temp[,"hi"],
  Date = as.Date("2015-01-01") %m+% months(seq_len(n_time) - 1)
)

# Plot temporal trend
ggplot(temporal_df, aes(x = Date, y = RR)) +
  geom_line(color = "darkgreen") +
  geom_ribbon(aes(ymin = Low, ymax = High), alpha = 0.2) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  labs(title = "Temporal Trend of Rate Ratios", x = "Date", y = "Rate Ratio") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Now build the dataframe
spatial_df <- tibble(
  Sec_ID = unique(dat$Sec_ID),
  RR     = RR_spatial,
  PP     = PP_spatial
) %>%
  mutate(
    RR_cat = cut(RR, breaks = c(0, 0.5, 1, 1.2, 2, 3, 6, 18), include.lowest = TRUE),
    PP_cat = cut(PP, breaks = c(0, 0.5, 0.95, 1), include.lowest = TRUE)
  )

# 14. SPATIAL MAPS (RR & PP) ------------------------------------------------
map_rr <- left_join(Rwa_data, spatial_df, by = "Sec_ID")
# Get number of spatial units
n_sectors <- length(unique(dat$Sec_ID))

# Extract only first n marginals (structured effects)
RR_spatial <- sapply(spatial_model_clima$marginals.random$ID[1:n_sectors],
                     function(m) inla.emarginal(function(x) exp(x), m))

PP_spatial <- sapply(spatial_model_clima$marginals.random$ID[1:n_sectors],
                     function(m) 1 - inla.pmarginal(0, m))



# RR Map
RR<-ggplot(map_rr) +
  geom_sf(aes(fill = RR_cat), color = NA) +
  scale_fill_brewer(palette = "PuOr", na.value = "grey90") +
  labs(title = "Spatial Posterior RR", fill = "RR Category") +
  theme_void() + theme(legend.position = "bottom")

# PP Map
PP<- ggplot(map_rr) +
  geom_sf(aes(fill = PP_cat), color = NA) +
  scale_fill_viridis_d(option = "plasma", direction = -1, na.value = "grey90") +
  labs(title = "Spatial Posterior Probability", fill = "PP Category") +
  theme_void() + theme(legend.position = "bottom")

RR|PP
# 15.1 Monthly mean climate and pop
climate_means <- dat %>%
  st_drop_geometry() %>%
  group_by(Sec_ID, Month) %>%
  summarise(
    Max_temperature   = mean(Max_temperature, na.rm = TRUE),
    Prec              = mean(Prec, na.rm = TRUE),
    Min_temperature   = mean(Min_temperature, na.rm = TRUE),
    Relative_humidity = mean(Relative_humidity, na.rm = TRUE),
    U_5_Population    = last(U_5_Population),
    .groups = "drop"
  )

# 15.2 Define next 6 months
last_obs_date <- as.Date(paste0(max(dat$Year), "-", sprintf("%02d", max(dat$Month)), "-01"))
future_months <- seq(from = last_obs_date %m+% months(1), by = "month", length.out = 6)

# 15.3 Create future data frame
future_df <- expand.grid(Sec_ID = unique(dat$Sec_ID), Date = future_months) %>%
  mutate(
    Year  = year(Date),
    Month = month(Date)
  ) %>%
  left_join(climate_means, by = c("Sec_ID", "Month")) %>%
  mutate(
    ADM3_EN = Rwa_data$ADM3_EN[match(Sec_ID, Rwa_data$Sec_ID)],
    ADM2_EN = Rwa_data$ADM2_EN[match(Sec_ID, Rwa_data$Sec_ID)],
    Diarrhoea_cases = NA,
    ID              = Sec_ID,
    ID.time         = max(dat$ID.time) + row_number(),
    ID.space.time   = max(dat$ID.space.time) + row_number(),
    offset          = log(U_5_Population / 1000)
  )

# 15.4 Combine past + future
dat_extended <- bind_rows(dat, future_df)

# 15.5 Forecast using fitted model
model_forecast <- inla(
  formula         = formula_bym2_climate,
  family          = "poisson",
  data            = dat_extended,
  offset          = dat_extended$offset,
  control.compute = list(dic = FALSE, waic = FALSE),
  control.predictor = list(compute = TRUE)
)

# 15.6 Extract forecast results
is_future <- is.na(dat_extended$Diarrhoea_cases)
pred_vals <- model_forecast$summary.fitted.values[is_future, c("mean", "0.025quant", "0.975quant")]

future_results <- future_df %>%
  bind_cols(pred_vals) %>%
  rename(
    Predicted_Mean  = mean,
    Lower_95_CI     = `0.025quant`,
    Upper_95_CI     = `0.975quant`
  )

# 15.7 Faceted maps of forecast (All Rwanda)
rwanda_grid <- Rwa_data %>%
  left_join(future_results, by = "Sec_ID") %>%
  mutate(MonthYear = format(Date, "%b %Y"))

ggplot(rwanda_grid) +
  geom_sf(aes(fill = Predicted_Mean), color = NA) +
  facet_wrap(~ MonthYear, ncol = 3) +
  scale_fill_viridis_c(option = "plasma", name = "Predicted Mean") +
  labs(title = "Forecasted Diarrhoea Cases Across Rwanda", subtitle = "Next 6 Months") +
  theme_void() +
  theme(
    strip.text    = element_text(size = 10, face = "bold"),
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.title  = element_text(size = 11),
    legend.text   = element_text(size = 10)
  )

# 1. Aggregate total cases per month across all sectors
country_level_df <- dat %>%
  st_drop_geometry() %>%
  mutate(Predicted = spatial_model_clima$summary.fitted.values$mean) %>%
  mutate(Date = as.Date(paste0(Year, "-", sprintf("%02d", Month), "-01"))) %>%
  group_by(Date) %>%
  summarise(
    Observed = sum(Diarrhoea_cases, na.rm = TRUE),
    Predicted = sum(Predicted, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Reshape for ggplot
plot_total_df <- country_level_df %>%
  pivot_longer(cols = c(Observed, Predicted),
               names_to = "Type", values_to = "Cases")

# 3. Plot
ggplot(plot_total_df, aes(x = Date, y = Cases, color = Type)) +
  geom_line(size = 1) +
  geom_point(size = 2.5) +
  geom_text(
    data = filter(plot_total_df, Type == "Observed"),
    aes(label = round(Cases, 0)),
    vjust = -0.7, size = 3, color = "black", show.legend = FALSE
  ) +
  geom_text(
    data = filter(plot_total_df, Type == "Predicted"),
    aes(label = round(Cases, 0)),
    vjust = 1.2, size = 3, color = "blue", show.legend = FALSE
  ) +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "blue")) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(
    title = "National Diarrhea Cases: Observed vs Predicted",
    subtitle = "Monthly Total (All Sectors)",
    x = "Month", y = "Number of Cases", color = ""
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

########################################################


library(dplyr)
library(scales)

ggplot(plot_total_df, aes(x = Date, y = Cases, color = Type)) +
  geom_line(size = 1) +
  geom_point(size = 2.5) +
  geom_text(
    data = dplyr::filter(plot_total_df, Type == "Observed"),
    aes(label = paste0("Obs=", scales::comma(round(Cases, 0)))),
    vjust = -0.7, size = 3, color = "black", show.legend = FALSE
  ) +
  geom_text(
    data = dplyr::filter(plot_total_df, Type == "Predicted"),
    aes(label = paste0("Pred=", scales::comma(round(Cases, 0)))),
    vjust = 1.2, size = 3, color = "blue", show.legend = FALSE
  ) +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "blue")) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(title = "National Diarrhea Cases: Observed vs Predicted",
       subtitle = "Monthly Total (All Sectors)",
       x = "Month", y = "Number of Cases", color = "") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


####################################################
# 1. Add predicted mean + 95% CI to dataset
dat_ci <- dat %>%
  st_drop_geometry() %>%
  mutate(
    Predicted_Mean = spatial_model_clima$summary.fitted.values$mean,
    Predicted_Low  = spatial_model_clima$summary.fitted.values$`0.025quant`,
    Predicted_High = spatial_model_clima$summary.fitted.values$`0.975quant`,
    Date           = as.Date(paste0(Year, "-", sprintf("%02d", Month), "-01"))
  )

# 2. Aggregate to country-level totals per month
country_level_ci <- dat_ci %>%
  group_by(Date) %>%
  summarise(
    Observed       = sum(Diarrhoea_cases, na.rm = TRUE),
    Predicted      = sum(Predicted_Mean, na.rm = TRUE),
    Predicted_Low  = sum(Predicted_Low, na.rm = TRUE),
    Predicted_High = sum(Predicted_High, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(year(Date) >= 2024 & year(Date) <= 2025)  # Restrict to 2024–2025

# 3. Reshape for plotting
plot_total_df <- country_level_ci %>%
  pivot_longer(cols = c(Observed, Predicted),
               names_to = "Type", values_to = "Cases")

# 4. Plot with confidence ribbon
ggplot() +
  # Confidence interval for predicted
  geom_ribbon(
    data = country_level_ci,
    aes(x = Date, ymin = Predicted_Low, ymax = Predicted_High),
    fill = "blue", alpha = 0.2
  ) +
  # Lines for predicted and observed
  geom_line(data = plot_total_df, aes(x = Date, y = Cases, color = Type), size = 1) +
  geom_point(data = plot_total_df, aes(x = Date, y = Cases, color = Type), size = 2.5) +
  # Labels
  geom_text(
    data = filter(plot_total_df, Type == "Observed"),
    aes(x = Date, y = Cases, label = round(Cases, 0)),
    vjust = -0.7, size = 3, color = "black", show.legend = FALSE
  ) +
  geom_text(
    data = filter(plot_total_df, Type == "Predicted"),
    aes(x = Date, y = Cases, label = round(Cases, 0)),
    vjust = 1.5, size = 3, color = "blue", show.legend = FALSE
  ) +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "blue")) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(
    title = "National Diarrhea Cases in Rwanda (2024–2025)",
    subtitle = "Observed vs Predicted with 95% Confidence Intervals",
    x = "Month", y = "Number of Cases", color = ""
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
##### to solve the iss of uge cI

