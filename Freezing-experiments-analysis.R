# This script assumes all the data in the working directory
library(magrittr)
library(dplyr)
library(lme4)
# Also requires ggplot2, ape, and phytools

# Read the data:
fr <- read.table('Freezing-experiments-data-with-metadata.csv',
                 header = T, sep = ',', skip = 14)[1:4]
fr.meta <- read.table('Freezing-experiments-data-with-metadata.csv',
                      header = T, sep = ',', skip = 1, nrows = 11)
fr %>% filter(grepl('AL3_15', SampleID)) %>% nrow # 144 individuals
fr %>% filter(!grepl('AL3_15', SampleID)) %>% nrow # 402 individuals

# Find the species with zero survival rate:
Zsp <- fr %>% group_by(SampleID) %>% summarize(Surv = sum(Survival)) %>%
  filter(Surv == 0) %$% SampleID %>% as.character

# Leave only the species with a positive survival rate:
fr %<>% filter(!SampleID %in% Zsp)
fr %>% filter(!grepl('AL3_15', SampleID)) %>% nrow # 318 individuals

# Get Region and Species from the table metadata
fr %<>%
  mutate(Region = fr.meta[match(SampleID, fr.meta$SampleID),'Region']) %>% 
  mutate(Species = fr.meta[match(SampleID, fr.meta$SampleID),'Species'])
  
# Median survival rate:
fr %>%
  group_by(SampleID, Plate) %>%
  summarize(Surv = sum(Survival)/n()*100) %>% 
  summarize(
    Median = median(Surv),
    MinSurv = min(Surv),
    MaxSurv = max(Surv)
    )

  
# Plotting --------------------------------------------------------------------
fr %>%
  group_by(Region, Species, Plate) %>%
  summarize(Surv = sum(Survival)/n()*100) %>% 
  ggplot2::ggplot(aes(Species, Surv, fill = Region)) + geom_boxplot() +
  scale_fill_manual(values = c(
    'Arctic' = 'deepskyblue',
    'Temporal' = 'beige')) +
  labs(y = 'Survival rate, %', x = NULL) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12,
                                   color = c(1,1,1,1,1,1,2,1)))


# Modeling --------------------------------------------------------------------
## Testing for the presence of the phylogenetic signal:
phy <- ape::read.nexus('Adineta-phylogeny.nxs')
fr_ag <- fr %>% group_by(SampleID) %>%
  summarise(Surv = sum(Survival)/n()) %>%
  as.data.frame
fr_ag <- fr_ag[fr_ag$SampleID %in% phy$tip.label, ] # not all guys are on the tree
fr_ag_v <- setNames(fr_ag$Surv, fr_ag$SampleID)
phytools::phylosig(phy, fr_ag_v, method = 'K', test = T)
phytools::phylosig(phy, fr_ag_v, method = 'lambda', test = T)
## No phylogenetic signal detected, so we use mixed-effects models
## assuming uncorrelated errors

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
fr %>% filter(Region == 'Arctic') %>%
  group_by(SampleID) %>% 
  summarise(N = n(), Surv = sum(Survival)) %$% 
  prop.test(Surv, N)
# Xsq = 0.63587; p = 0.4252

# Diff. in prop. between the ancient species and the sister contemp. species,
# A. vaga isolate Hprim14, GenBank KU860644.1:
fr %>% filter(grepl('Unamur|AL3_15', SampleID)) %>%
  group_by(SampleID) %>% 
  summarise(N = n(), Surv = sum(Survival)) %$% 
  prop.test(Surv, N)
# Xsq = 2.8572; p = 0.09096
