library(Bchron)
ages <- BchronCalibrate(
  ages = c(26260, 20158, 18228),
  ageSds = c(210, 99, 78),
  calCurves = c('intcal13', 'intcal13', 'intcal13')
)
summary(ages, prob = 99)
sapply(ages, function(x) round(mean(x$ageGrid)))
sapply(ages, function(x) round(sd(x$ageGrid)))
