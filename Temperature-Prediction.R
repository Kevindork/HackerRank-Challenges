############################
# Temperature-Prediction.R #
############################

# A simple algorithm for predicting missing monthly average high and low 
# temperatures, based on kernel smoothing deviations from the norm and an
# ad-hoc bias correction

### Parameters
# Number of Observations
n <- as.numeric(readLines("Temperature-Prediction-Input.txt", 1))
# One sided kernel bandwidth, in years
b <- 2
# Minimum number of nonempty entries required in each kernel
m <- 12

### Reading and Cleaning Data
sol <- read.table("Temperature-Prediction-Output.txt", header=F)[[1]]
infile <- read.table("Temperature-Prediction-Input.txt", skip=1, header=T, 
                     na.strings=paste("Missing_", 1:1500, sep=""))
infile$month <- match(infile$month, month.name)
infile$year <- infile$yyyy + (infile$month - 1) / 12
i.tmax.na <- is.na(infile$tmax)
i.tmin.na <- is.na(infile$tmin)

### Overall average temperatures by month
tmax.months <- aggregate(infile$tmax[!i.tmax.na], 
                         by=list(infile$month[!i.tmax.na]), FUN=mean)
tmin.months <- aggregate(infile$tmin[!i.tmin.na], 
                         by=list(infile$month[!i.tmin.na]), FUN=mean)

### Deviations from the monthly averages
infile$tmax.std <- infile$tmax - tmax.months[infile$month, 2]
infile$tmax.std <- infile$tmax.std / sd(infile$tmax.std)
infile$tmin.std <- infile$tmin - tmin.months[infile$month, 2]
infile$tmin.std <- infile$tmin.std / sd(infile$tmin.std)

### Kernel smoother
EstTemp <- function(i, year, i.missing, std, month, min.kernsize=m, bw=b) {
  t.kern <- year[i] + c(-bw,bw)
  i.kern <- which(year > t.kern[1] & year < t.kern[2])
  na.kern <- i.missing[i.kern] 
  while(sum(!na.kern) < min.kernsize) {
    bw <- bw + 1/4
    t.kern <- year[i] + c(-bw, bw)
    i.kern <- which(year > t.kern[1] & year < t.kern[2])
    na.kern <- i.missing[i.kern]  
  }
  w <- dnorm(year[i.kern], year[i], bw / 2) * !na.kern
  w <- w / sum(w)
  return(c(dev=sum(w * std[i.kern], na.rm=T), 
           avg=month[infile$month[i], 2]))
}

### 1st pass prediction (Kernel smoothing)
i.nonmissing <- !i.tmax.na & !i.tmin.na
n.nonmissing <- sum(i.nonmissing)
pred.all <- data.frame(array(dim=c(n, 8)))
names(pred.all) <- c("pred.max", "err.max", "dev.max", "avg.max", 
                     "pred.min", "err.min", "dev.min", "avg.min")
for(i in 1:n) {
  pred.all[i, 3:4] <- EstTemp(i, infile$year, i.tmax.na, infile$tmax.std, 
                              tmax.months)
  pred.all[i, 7:8] <- EstTemp(i, infile$year, i.tmin.na, infile$tmin.std, 
                              tmin.months)
}
pred.all[, 1] <- rowSums(pred.all[, 3:4])
pred.all[, 5] <- rowSums(pred.all[, 7:8])
pred.all[, 2] <- pred.all[, 1] - infile$tmax
pred.all[, 6] <- pred.all[, 5] - infile$tmin

### 2nd pass prediction (Bias correction)
lm.ymax.xmin <- lm(err.max ~ err.min + I(err.min^2), data=pred.all)
lm.ymin.xmax <- lm(err.min ~ err.max + I(err.max^2), data=pred.all)

bias.max <- predict(lm.ymax.xmin, new=pred.all, se.fit=F)
bias.min <- predict(lm.ymin.xmax, new=pred.all, se.fit=F)
pred.max <- pred.all[which(i.tmax.na), 1] - bias.max[which(i.tmax.na)]
pred.min <- pred.all[which(i.tmin.na), 5] - bias.min[which(i.tmin.na)]

pred.max <- data.frame(temp=pred.max, index=which(i.tmax.na))
pred.min <- data.frame(temp=pred.min, index=which(i.tmin.na))
pred <- rbind(pred.max, pred.min)
pred <- pred[order(pred$index),]

# Training Set Score
# 1-mean(abs(pred$temp - sol))/5

write(pred$temp, stdout(), sep="\n")
