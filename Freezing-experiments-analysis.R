library(magrittr)
library(readxl)
library(dplyr)
library(ggplot2)
library(lme4)

# Read the data:
Names <- c('SampleID', 'Species', 'Age', 'Continent', 'Region', 'Country',
          'Plate', 'Individual', 'Survival')
fr <- read_xlsx('Supplementary File S3.xlsx',
                sheet = 'Freezing survival data',
                range = cell_rows(6:553),
                col_names = Names)
fr %>% filter(!grepl('PERMAFROST', Species)) %>% nrow # 404 individuals

# Find the species with zero survival rate:
Zsp <- fr %>% group_by(Species) %>% summarize(Surv = sum(Survival)) %>%
  filter(Surv == 0) %$% Species %>% as.character

# Leave only the species with a positive survival rate:
fr %<>% filter(!Species %in% Zsp)
fr %>% filter(!grepl('PERMAFROST', Species)) %>% nrow # 330 individuals

# Median survival rate:
fr %>%
  group_by(Age, Plate) %>%
  summarize(Surv = sum(Survival)/n()*100) %>% 
  group_by(Age) %>% 
  summarize(
    Median = median(Surv),
    MinSurv = min(Surv),
    MaxSurv = max(Surv)
    )

  
# Plotting --------------------------------------------------------------------
fr %>%
  mutate(Species = sub('SP1', 'Hprim14', Species)) %>% 
  group_by(Region, Species, Plate) %>%
  summarize(Surv = sum(Survival)/n()*100) %>% 
  ggplot(aes(Species, Surv, fill = Region)) + geom_boxplot() +
  scale_fill_manual(values = c(
    'arctic' = 'deepskyblue',
    'temporal' = 'beige')) +
  labs(y = 'Survival rate, %', x = NULL) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   color = c(1,1,1,1,1,1,2,1)))


# Modeling --------------------------------------------------------------------
## Model controls:
modContr <- glmerControl(
  check.conv.singular = .makeCC(action = 'message', tol = 1e-5),
  check.conv.grad = .makeCC('warning', tol = 1e-2))

## Intercept-only model, Plate random:
fr.glmer00 <- fr %>%
  glmer(Survival ~ 1 + (1|Plate), data = .,
        family = 'binomial', control = modContr)
summary(fr.glmer00)

## Species fixed, Plate random:
fr.glmer01 <- fr %>%
  glmer(Survival ~ Species + (1|Plate), data = .,
        family = 'binomial', control = modContr)
summary(fr.glmer01)
anova(fr.glmer01, fr.glmer00)

## Region fixed, Plate random:
fr.glmer02 <- fr %>%
  glmer(Survival ~ Region + (1|Plate), data = .,
        family = 'binomial', control = modContr)
summary(fr.glmer02)
anova(fr.glmer02, fr.glmer00)

## Intercept-only; Species and Plate random; Plate nested in Species:
fr.glmer03 <- fr %>%
  glmer(Survival ~ 1 + (1|Species/Plate), data = .,
        family = 'binomial', control = modContr)
summary(fr.glmer03)
anova(fr.glmer03, fr.glmer00)

## Region fixed; Species and Plate random; Plate nested in Species:
fr.glmer04 <- fr %>%
  glmer(Survival ~ Region + (1|Species/Plate), data = .,
        family = 'binomial', control = modContr)
summary(fr.glmer04)
anova(fr.glmer04, fr.glmer03)


# Pairwise comparisons --------------------------------------------------------

# Diff. in prop. between the contemp. and ancient Arctic species:
fr %>% filter(Region == 'arctic') %>%
  group_by(Species) %>% 
  summarise(N = n(), Surv = sum(Survival)) %$% 
  prop.test(Surv, N, alternative = 'less') # A. glauca has smaller median
# Xsq = 0.63587; p = 0.2126

# Diff. in prop. between the ancient species and the sister contemp. species,
# A. vaga isolate Hprim14, GenBank KU860644.1:
fr %>% filter(grepl('Unamur|AL3_15', SampleID)) %>%
  group_by(Species) %>% 
  summarise(N = n(), Surv = sum(Survival)) %$% 
  prop.test(Surv, N, alternative = 'less') # Hprim14 has smaller median
# Xsq = 2.8572; p = 0.04548
