library(rbacon)
library(dbplyr)

#functions to compute weighted mean and standard deviation (from SDMTools)
wt.mean <- function(x, wt) {
  s = which(is.finite(x*wt)); wt = wt[s]; x = x[s] #remove NA info
  sum(wt * x)/sum(wt)
}

wt.sd <- function(x, wt) {
  s = which(is.finite(x + wt)); wt = wt[s]; x = x[s] #remove NA info
  xbar = wt.mean(x,wt) #get the weighted mean
  sqrt(sum(wt *(x-xbar)^2)*(sum(wt)/(sum(wt)^2-sum(wt^2))))
}


#params
thick<-10
d.by<-1
core.name="tram20.v3"
acc.mean=20
acc.shape=1.5
mem.strength=10
mem.mean=0.5
ssize=8000
# mem.strength * (1 - mem.mean) #should be smaller than 1

x <- Bacon(
  coredir = "/home/graciela/Nextcloud/saco.csic.es/1 graciela/projects/Plan nacional/PYCACHU/PYC-TRAM20/CHRONO/preliminary/bacon model/Bacon_runs",
  core.name,
  acc.mean=acc.mean,
  acc.shape=acc.shape,
  mem.mean=mem.mean,
  mem.strength=mem.strength,
  thick=thick,
  d.by=d.by,
  ask=FALSE,
  suggest=FALSE,
  plot.pdf=TRUE,
  depths.file=FALSE,
  normal=FALSE,
  rotate.axes=TRUE,
  ssize=ssize
  )

#extracting calibration object from info
calibration <- info$calib$probs

#every object in calibration is a matrix. The first column is the age, and the second column is the probability
calibration[[1]]

#getting ages and probs from the matrix
ages <- calibration[[1]][, 1]
probs <- calibration[[1]][, 2]

#calibration curve
plot(ages, probs, type = "l")

#weighted mean of first date
wt.mean(x = ages, wt = probs)

#weighted standard deviation of first date
wt.sd(x = ages, wt = probs)

#the important part

#weighted mean for every date
calibration.mean <- lapply(
  calibration,
  FUN = function(x) wt.mean(x = x[, 1], wt = x[, 2])
  ) %>%
  unlist()

#weighted sd for every date
calibration.sd <- lapply(
  calibration,
  FUN = function(x) wt.sd(x = x[, 1], wt = x[, 2])
  ) %>%
  unlist()

#to the dates data frame
df <- read.table(
  "/home/graciela/Nextcloud/saco.csic.es/1 graciela/projects/Plan nacional/PYCACHU/PYC-TRAM20/CHRONO/preliminary/bacon model/Bacon_runs/tram20.v3/tram20.v3.csv",
  header = TRUE,
  sep = ","
  ) %>%
  dplyr::mutate(
    calibration.mean = calibration.mean,
    calibration.sd = calibration.sd
  )
calibrated<-write.csv(df,"/home/graciela/Nextcloud/saco.csic.es/1 graciela/projects/Plan nacional/PYCACHU/PYC-TRAM20/CHRONO/preliminary/bacon model/Bacon_runs/tram20.v3/", row.names = FALSE)
save(calibrated.csv, file=calibrated.csv.RData)


