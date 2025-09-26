# ===============================================================
# Local Trend Analysis of Amazonian Fish in Reservoirs
# ===============================================================

# DATA PREPARATION
# ---------------------------------------------------------------
library(tidyverse)
library(parameters)
library(mgcv)         
library(ggeffects)    
library(gratia)       
library(emmeans)      
library(performance)  
library(patchwork)    
library(MASS)         
library(DHARMa)      
library(jtools)
library(MuMIn)
library(readxl)
library(openxlsx)
library(effects)
library(gt)

options(scipen = 99)

# Resolve function conflicts
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)

# Minimal figure theme
fig_theme <- theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  strip.background = element_blank(),
  strip.text.y = element_text(),
  legend.background = element_blank(),
  legend.key = element_blank()
)

# Load data
setwd("C:/Aymar Backup/aymar/Aymar/Pós/Doutorado")
summary_data <- read.xlsx("parana_study_case.xlsx")

analysis_data <- summary_data %>%
  mutate(
    reservoir = as.factor(reservoir),
    station   = as.factor(station),
    year      = as.numeric(year)
  )

# ---------------------------------------------------------------
# 2. GLOBAL TREND ANALYSIS
# ---------------------------------------------------------------

# Global Poisson GLM (linear trend)
m_global_linear <- glm(
  Amazon ~ year + reservoir,
  family = poisson(link = "log"), na.action = na.omit,
  data = analysis_data
)

# Global GAM (smooth trend)
m_global_gam <- gam(
  Amazon ~ s(year) + reservoir,
  family = poisson(link = "log"),
  method = "ML", na.action = na.omit,
  data = analysis_data
)

# Model comparison
model_comp <- compare_performance(
  Linear = m_global_linear,
  GAM = m_global_gam,
  metrics = c("AICc", "R2")
)
model_comp
m_global_gam <- gam(
  Amazon ~ s(year) + reservoir,
  family = poisson(link = "log"),
  na.action = na.omit,
  data = analysis_data
)
# Check model diagnostics
check_overdispersion(m_global_gam)
summary(m_global_gam)

# Derivatives (rates of change)
dev_global <- gratia::derivatives(m_global_gam, term = "s(year)")
print(dev_global, n = Inf)

# Predictions for visualization
global_pred <- ggpredict(m_global_gam, terms = "year [all]") %>%
  as_tibble() %>%
  rename(year = x) %>%
  as.data.frame()

# ---------------------------------------------------------------
# 3. LOCAL ANALYSES (by reservoir)
# ---------------------------------------------------------------
# Split dataset
reservoir_data <- split(analysis_data, analysis_data$reservoir)

#Três Irmãos Reservoir
tir_reservoir <- reservoir_data[["Três Irmãos"]]

m_tir_linear <- glm(
  Amazon ~ year,
  family = poisson(link = "log"), na.action = na.omit,
  data = tir_reservoir
)

m_tir_gam <- gam(
  Amazon ~ s(year),
  family = poisson(link = "log"),
  method = "ML", na.action = na.omit,
  data = tir_reservoir
)

model_comp_tir <- compare_performance(
  Linear = m_tir_linear,
  GAM = m_tir_gam,
  metrics = c("AICc", "R2")
)
model_comp_tir
m_tir_gam <- gam(
  Amazon ~ s(year),
  family = poisson(link = "log"),
  na.action = na.omit,
  data = tir_reservoir)

tir_pred <- ggpredict(m_tir_gam, terms = "year [all]") %>%
  as_tibble() %>%
  rename(year = x) %>%
  as.data.frame()

#### LOCAL Jupia
jup_reservoir <- subset(analysis_data, reservoir == "Jupiá")

m_jup_linear <- glm(Amazon ~ year ,
                    family = poisson(link = "log"), na.action = na.omit,
                    data = jup_reservoir)
m_jup_gam <- gam(Amazon ~ s(year),
                 family = poisson(link = "log"),
                 method = "ML", na.action = na.omit,
                 data = jup_reservoir)

model_comp_jup <- compare_performance(
  Linear = m_jup_linear,
  GAM = m_jup_gam,
  metrics = c("AICc", "R2")
)
model_comp_jup
AICc(m_tir_linear)
AICc(m_tir_gam)

m_jup_gam <- gam(Amazon ~ s(year),
                 family = poisson(link = "log"),
                 na.action = na.omit,
                 data = jup_reservoir)
dev_jup<- gratia::derivatives(m_jup_gam, select = "s(year)")
print(dev_jup,n=Inf)
jup_pred <- ggpredict(m_jup_gam, terms = "year [all]") |> 
  as.data.frame() |> 
  dplyr::rename(year = x)


###LOCAL porto
porto_reservoir <- subset(analysis_data, reservoir == "Porto primavera")

m_ppr_linear <- glm(Amazon ~ year ,
                    family = poisson(link = "log"), na.action = na.omit,
                    data = porto_reservoir)

m_ppr_gam <- gam(Amazon ~ s(year),
                 family = poisson(link = "log"),
                 method = "ML", na.action = na.omit,
                 data = porto_reservoir)

model_comp_ppr <- compare_performance(
  Linear = m_ppr_linear,
  GAM = m_ppr_gam,
  metrics = c("AICc", "R2"))
model_comp_ppr
AICc(m_ppr_linear)
AICc(m_ppr_gam)

porto_pred = ggpredict(m_ppr_linear, terms = "year [all]") |> 
  as.data.frame() |> 
  dplyr::rename(year = x)


###LOCAL Ilha
isa_reservoir <- subset(analysis_data, reservoir == "I. Solteira")

m_isa_linear <- glm(Amazon ~ year ,
                    family = poisson(link = "log"), na.action = na.omit,
                    data = isa_reservoir)

m_isa_gam <- gam(Amazon ~ s(year),
                 family = poisson(link = "log"),
                 method = "ML", na.action = na.omit,
                 data = isa_reservoir)

model_comp_isa <- compare_performance(
  Linear = m_isa_linear,
  GAM = m_isa_gam,
  metrics = c("AICc", "R2"))
(model_comp_isa)
AICc(m_isa_linear)
AICc(m_isa_gam)

isa_pred <- ggpredict(m_isa_linear, terms = "year [all]") |> 
  as.data.frame() |> 
  rename(year = x)

# ---------------------------------------------------------------
# 4. SUMMARIZE RESULTS
# ---------------------------------------------------------------

results_summary <- data.frame(
  Reservoir = c("Global", "Três Irmãos", "Jupiá", "Porto Primavera", "I. Solteira"),
  AICc_Linear = c(model_comp$AICc[1], model_comp_tir$AICc[1], model_comp_jup$AICc[1],
                  model_comp_ppr$AICc[1], model_comp_isa$AICc[1]),
  AICc_GAM = c(model_comp$AICc[2], model_comp_tir$AICc[2], model_comp_jup$AICc[2],
               model_comp_ppr$AICc[2], model_comp_isa$AICc[2]),
  ΔAICc = sapply(list(model_comp, model_comp_tir, model_comp_jup, model_comp_ppr, model_comp_isa),
                 function(x) x$AICc[2] - x$AICc[1]),
  Best_Model = c(
    ifelse(model_comp$AICc[2] < model_comp$AICc[1] - 2, "GAM", "Linear"),
    ifelse(model_comp_tir$AICc[2] < model_comp_tir$AICc[1] - 2, "GAM", "Linear"),
    ifelse(model_comp_jup$AICc[2] < model_comp_jup$AICc[1] - 2, "GAM", "Linear"),
    ifelse(model_comp_ppr$AICc[2] < model_comp_ppr$AICc[1] - 2, "GAM", "Linear"),
    ifelse(model_comp_isa$AICc[2] < model_comp_isa$AICc[1] - 2, "GAM", "Linear")
  )
)

metricas_modelos <- bind_rows(
  model_comp %>% mutate(Reservoir = "Global"),
  model_comp_tir %>% mutate(Reservoir = "Três Irmãos"),
  model_comp_jup %>% mutate(Reservoir = "Jupiá"),
  model_comp_ppr %>% mutate(Reservoir = "Porto Primavera"),
  model_comp_isa %>% mutate(Reservoir = "I. Solteira")
) %>% 
  select(Reservoir, R2, R2_Nagelkerke) %>% 
  group_by(Reservoir) %>% 
  mutate(Model = c("Linear", "GAM")) %>% 
  tidyr::pivot_wider(names_from = Model, values_from = c(R2, R2_Nagelkerke))%>%
  select(where(~ !all(is.na(.))))
results_summary <- results_summary %>%
  left_join(metricas_modelos, by = "Reservoir") %>%
  dplyr::mutate(across(where(is.numeric), round, 2))

knitr::kable(results_summary, digits = 2)

# ---------------------------------------------------------------
# 5. BINOMIAL MODELS (Proportion of Amazonian Fish)
# ---------------------------------------------------------------

m_prop_linear <- glm(
  cbind(prop_amazon * total, (1 - prop_amazon) * total) ~ year + reservoir,
  family = binomial(link = "logit"), data = summary_data
)

m_prop_gam <- gam(
  cbind(prop_amazon * total, (1 - prop_amazon) * total) ~ s(year) + reservoir,
  data = summary_data,
  method = "ML", family = binomial("logit")
)

model_comp_prop <- compare_performance(
  Linear = m_prop_linear,
  GAM = m_prop_gam,
  metrics = c("AICc", "R2")
)
model_comp_prop
# Predictions
global_pred_prop <- ggpredict(m_prop_linear, terms = "year [all]") %>%
  as_tibble() %>%
  rename(year = x) %>%
  as.data.frame()

### LOCAL
reservoir_data <- split(summary_data, summary_data$reservoir)
names(reservoir_data)
# "Três Irmãos"
tir_reservoir <- reservoir_data[["Três Irmãos"]]
tir_prop_linear <- glm(cbind(prop_amazon * total, (1 - prop_amazon) * total) ~ year ,
                       family = binomial(link = "logit"), data = tir_reservoir)

tir_prop_gam <- gam(cbind(prop_amazon * total, (1 - prop_amazon) * total) ~ s(year), 
                    data = tir_reservoir,
                    method = "ML", family = binomial("logit"))
model_comp_tir_prop <- compare_performance(
  Linear = tir_prop_linear,
  GAM = tir_prop_gam,
  metrics = c("AICc", "R2")
)
model_comp_tir_prop

tir_pred_prop <- ggpredict(tir_prop_linear, terms = "year [all]") |> 
  as_tibble() |> 
  rename(year = x) %>% as.data.frame(.)

#####JUP
jup_reservoir <- reservoir_data[["Jupiá"]]
jup_prop_linear <- glm(cbind(prop_amazon * total, (1 - prop_amazon) * total) ~ year ,
                       family = binomial(link = "logit"), data = jup_reservoir)

jup_prop_gam <- mgcv::gam(cbind(prop_amazon * total, (1 - prop_amazon) * total) ~ s(year), 
                          data = jup_reservoir,
                          method = "ML", family = binomial("logit"))

model_comp_jup_prop <- performance::compare_performance(
  Linear = jup_prop_linear,
  GAM = jup_prop_gam,
  metrics = c("AICc", "R2"))
model_comp_jup_prop

jp_pred_prop <- ggpredict(jup_prop_linear, terms = "year [all]") |> 
  as_tibble() |> 
  rename(year = x) %>% as.data.frame(.)

#### ISA 
isa_reservoir <- reservoir_data[["I. Solteira"]]
isa_prop_linear <- glm(cbind(prop_amazon * total, (1 - prop_amazon) * total) ~ year ,
                       family = binomial(link = "logit"), data = isa_reservoir)

isa_prop_gam <- gam(cbind(prop_amazon * total, (1 - prop_amazon) * total) ~ s(year), 
                    data = isa_reservoir,
                    family = binomial("logit"))

model_comp_isa_prop <- compare_performance(
  Linear = isa_prop_linear,
  GAM = isa_prop_gam,
  metrics = c("AICc", "R2"))
model_comp_isa_prop
isa_pred_prop <- ggpredict(isa_prop_linear, terms = "year [all]") |> 
  as_tibble() |> 
  rename(year = x) %>% as.data.frame(.)

### PORTO PRIMAVERA
ppr_reservoir <- reservoir_data[["Porto primavera"]]
ppr_prop_linear <- glm(cbind(prop_amazon * total, (1 - prop_amazon) * total) ~ year ,
                       family = binomial(link = "logit"), data = ppr_reservoir)

ppr_prop_gam <- gam(cbind(prop_amazon * total, (1 - prop_amazon) * total) ~ s(year), 
                    data = ppr_reservoir,
                    method = "REML", family = binomial("logit"))

model_comp_ppr_prop <- compare_performance(
  Linear = ppr_prop_linear,
  GAM = ppr_prop_gam,
  metrics = c("AICc", "R2"))
model_comp_ppr_prop
ppr_pred_prop<- ggpredict(ppr_prop_linear, terms = "year [all]") |> 
  as_tibble() |> 
  rename(year = x) %>% as.data.frame(.)

# Combine results for binomial models
results_summary_prop <- data.frame(
  Reservoir = c("Global", "Três Irmãos", "Jupiá", "Porto Primavera", "I. Solteira"),
  AICc_Linear = c(model_comp_prop$AICc[1], model_comp_tir_prop$AICc[1], 
                  model_comp_jup_prop$AICc[1], model_comp_ppr_prop$AICc[1],
                  model_comp_isa_prop$AICc[1]),
  AICc_GAM = c(model_comp_prop$AICc[2], model_comp_tir_prop$AICc[2],
               model_comp_jup_prop$AICc[2], model_comp_ppr_prop$AICc[2],
               model_comp_isa_prop$AICc[2]),
  ΔAICc = sapply(list(model_comp_prop, model_comp_tir_prop, model_comp_jup_prop,
                      model_comp_ppr_prop, model_comp_isa_prop),
                 function(x) x$AICc[2] - x$AICc[1]),
  Best_Model = c(
    ifelse(model_comp_prop$AICc[2] < model_comp_prop$AICc[1] - 2, "GAM", "Linear"),
    ifelse(model_comp_tir_prop$AICc[2] < model_comp_tir_prop$AICc[1] - 2, "GAM", "Linear"),
    ifelse(model_comp_jup_prop$AICc[2] < model_comp_jup_prop$AICc[1] - 2, "GAM", "Linear"),
    ifelse(model_comp_ppr_prop$AICc[2] < model_comp_ppr_prop$AICc[1] - 2, "GAM", "Linear"),
    ifelse(model_comp_isa_prop$AICc[2] < model_comp_isa_prop$AICc[1] - 2, "GAM", "Linear")
  )
)
results_summary_prop


