######################################################################
## estimation of the truncated skew normal distribution
######################################################################

nd <- rnorm(1000)
fnorm(c(0, 10),  nd)

optim(c(0, 10), fnorm, x = nd)

snd <- nd * pnorm(-0.2 * nd)
hist(snd, breaks = 100)

tsnd <- snd[snd < 0.5]
hist(tsnd, breaks = 100)

ftnorm <- function(beta, x){
    -sum(dtnorm(x, beta[1], beta[2], upper = 0.5, log = TRUE) *
         ptnorm(beta[3] * x, beta[1], beta[2], upper = 0.5, log = TRUE))
}

est <- optim(c(0, 10, 0.5), fnorm, x = tsnd)

hist(tnd, breaks = 100, freq = FALSE)
curve(dtnorm(x, est$par[1], est$par[2], upper = 0.5) *
      pnorm(x * est$par[3] * x, est$par[1], est$par[2]),
      add = TRUE, col = "red")



library(copula)
library(nortest)

## Generate the gumbel copula
myCop.gumbel <- archmCopula(family = "gumbel", 0.1, dim = 2, param = 2)

## Specify the marginal of the joint distribution, here the marginal
## is set to N(0, 25) and N(10, 1). Feel free to change them
myMvd <- mvdc(copula = myCop.gumbel,
              margins = c("norm", "norm"),
              paramMargins = list(list(mean = 0, sd = 20),
                             list(mean = 10, sd = 1)))

## Generate samples fromthe joint distribution, a small sample of 100
## has taken to minimise the effect of CLT.
set.seed(123)
myJoint <- rmvdc(myMvd, 10000)

## test the marginal of the joint distribution, Anderson-Darling is
## one of the most stringent test for normality compare to many well
## know test such as KS, JB, shapiro ....
ad.test(myJoint[, 1])
ad.test(myJoint[, 2])

## This plot shows that the joint is clearly non-normal and has the
## typical gumbel copula shape
plot(hexbin(myJoint[, 1], myJoint[, 2]))

## Now take the difference and test, and we see no evidence of
## normality
normdiff <- myJoint[, 2] - myJoint[, 1]
ad.test(normdiff)
hist(normdiff, freq = FALSE)
curve(dnorm(x, mean(normdiff), sd(normdiff)), col = "red", add = TRUE)


rnorm(1000, 10, 1)


## Toy example for gumbel copula with log-normal distribution
## (Taken from the documentation of copula::fitMvdc)

## Specify the copula
gumbel.cop <- gumbelCopula(3, dim=2)
myMvd <- mvdc(gumbel.cop, c("lnorm","lnorm"), list(list(meanlog = 1.17),
                                                   list(meanlog = 0.76)))
## Generate some random sample to test
x <- rmvdc(myMvd, 1000)

## Fit the random sample
fit <- fitMvdc(x, myMvd, c(1,1,2))
fit

##  Jun Yan (2007). Enjoy the Joy of Copulas: With a Package copula.
## Journal of Statistical Software, 21(4), 1-21. URL
##  http://www.jstatsoft.org/v21/i04/.

##  Ivan Kojadinovic, Jun Yan (2010). Modeling Multivariate Distributions
##  with Continuous Margins Using the copula R Package. Journal of
##  Statistical Software, 34(9), 1-20. URL
##  http://www.jstatsoft.org/v34/i09/.







######################################################################
## Ross's random walk codes
######################################################################

rw2d1 <- function(n) {
    xpos = ypos = numeric(n)
    xdir = c(TRUE, FALSE)
    pm1 = c(1, -1)
    for(i in 2:n)
        if (sample(xdir, 1)) {
            xpos[i] = xpos[i - 1] + sample(pm1, 1)
            yplos = ypos[i - 1]
        }
    else {
        xpos[i] = xpos[i - 1]
        ypos[i] = ypos[i - 1] + sample(pm1, 1)
    }
    list(x = xpos, y = ypos)
}

system.time(rw2d1(100000))

rw2d2 <- function(n) {
    steps = sample(c(-1, 1), n - 1, replace = TRUE)
    xdir = sample(c(TRUE, FALSE), n - 1, replace = TRUE)
    xpos = c(0, cumsum(ifelse(xdir, steps, 0)))
    ypos = c(0, cumsum(ifelse(xdir, 0, steps)))
    list(x = xpos, y = ypos)
}

system.time(rw2d2(100000))

rw2d3 <- function(n) {
    xsteps = c(-1, 1, 0, 0)
    ysteps = c( 0, 0, -1, 1)
    dir = sample(1:4, n - 1, replace = TRUE)
    xpos = c(0, cumsum(xsteps[dir]))
    ypos = c(0, cumsum(ysteps[dir]))
    list(x = xpos, y = ypos)
}

system.time(rw2d3(100000))


Rprof()
for(i in 1:100)
    pos = rw2d2(100000)
Rprof(NULL)




######################################################################
## Populate MoE table
######################################################################

dates <- as.numeric(unlist(read.csv(file = "dates.csv", header = FALSE)))
member <- as.numeric(unlist(read.csv(file = "members.csv", header = FALSE)))


pop <- data.frame(member = rep(member, each = length(dates)),
                  dates = rep(dates, length(member)))


######################################################################
## check whether scaling part of the data changes coefficient of
## others
######################################################################


x <- rnorm(100, 5, 10)
y <- 20 + 3 * x + rnorm(100, 0, 20)
test.df <- data.frame(x = x, y = y)

plot(test.df)

formula <- as.formula("y ~ x")

summary(lm(formula, data = test.df))

formula <- as.formula("y ~ poly(x, 1)")


######################################################################
## Detemnine breaks in plot
######################################################################

td <- data.frame(a = rnorm(1000, 1000, 100), b = rnorm(1000, 1000, 100))

labPlot <- function(mydata, n.breaks, mylab){
    require(labeling)
    mybreaks = heckbert(min(mydata$a), max(mydata$b), n.breaks)
    bt = length(mybreaks)
    qplot(x = a, y = b, data = mydata) +
        scale_y_continuous(limits = range(mybreaks),
                           breaks = mybreaks,
                           labels = c(mybreaks[-bt], mylab))
}

labPlot(td, 5, "Cm")

labPlot(td, 3, "Cm")
labPlot(td, 5, "Cm")
labPlot(td, 10, "Cm")


######################################################################
## More new investigation from Matthieu
######################################################################

library(ggplot2)
DF <- data.frame(x=runif(100),
                 fac = factor(sample(c("A","B"), size=100, replace=TRUE)),
                 longFac = factor(rep(c("A long long long name",
                                        "Another long name",
                                        "Me too I am small, but less",
                                        "I am long also", "I'm short"),
                 each = 20)))

qplot(x=longFac, y = x, data = DF, geom = "bar", stat = "identity",
      fill = fac, position="fill")

ggplot(data = DF, aes(x = longFac)) + geom_bar(aes(x * 100))


qplot(x = fac, y = x, data = DF, geom = "bar", stat = "identity",
      fill = longFac, position = "fill") +
    opts(legend.position = "top",legend.direction = "horizontal") +
    scale_color_discrete(breaks = unique(DF$longFac),
                         labels = c("a", "b", "c", "d", "e"))


######################################################################
## Investigate the speed of different subset
######################################################################

library(compiler)
library(rbenchmark)

test.df <- data.frame(names = rep(LETTERS, 1000),
                      value = rep(1:26, 1000))



fa <- function(df){
    df[which(df$value == max(df$value)), ]
}
fa2 <- cmpfun(fa)

fb <- function(df){
    df[which.max(df$value), ]
}
fb2 <- cmpfun(fb)

fc <- function(df){
    subset(df, subset = value == max(df$value, na.rm = T))
}
fc2 <- cmpfun(fc)

benchmark(fa(test.df), fa2(test.df), fb(test.df), fb2(test.df),
          fc(test.df), fc2(test.df), order = "relative",
          replications = 10000)



######################################################################
## Geom Bar and stats identity
######################################################################

library(ggplot2)
ggplot(diamonds, aes(color, fill=cut)) + geom_bar(stats = "identity") +
    coord_flip()

ggplot(diamonds, aes(color, fill=cut)) + geom_bar() + coord_flip()

qplot(cut, meanprice, geom="bar")
dev.new()
qplot(cut, meanprice, geom="bar", stat = "identity")


######################################################################
## Avoiding loop with multiple objects for Filipo
######################################################################

## Create fake data for each year
for(i in 1:10){
    assign(paste("a_", i, sep = ""), data.frame(a = LETTERS, b = 1:26/i))
}

## Try to minus a2 ~ a9 from a1, replace this with merge_recurse
ot <- a_2
for(i in 3:10){
    ot <- join(ot, get(paste("a_", i, sep = "")), by = "a", type = "full")
}

new.ot <- data.frame(a = ot["a"], b.sum = rowSums(ot[, -1]))

full.ot <- join(a_1, new.ot, by = "a", type = "full")

(final.ot <- data.frame(a = full.ot["a"], b = full.ot["b"] - full.ot["b.sum"]))

library(reshape)
merge_recurse(dfs = get(paste("a_", 1:10, sep = "")), by = "a")


merge_recurse(dfs = get(paste("a_", 1:10, sep = "")))



merge_recurse(list(a_2, a_3, a_4), by = "a")


merge(merge(merge(a_1, a_2, by = "a"), a_3, by = "a"), a_4, by = "a")


######################################################################
## Use existing data.frame
######################################################################

## Create fake data
faostat.df <- data.frame(cbind(1:10, matrix(rnorm(1000), nc = 100)))
colnames(faostat.df) <- c("faostat",
                          paste(rep(LETTERS[1:10], each = 10),
                                rep(2000:2010, length.out = 10),
                                sep = "-"))

## Extract the A component
afaostat.df <- faostat.df[c("faostat", colnames(faostat.df)
                              [grep("A", colnames(faostat.df))])]
aclass.df <- melt(afaostat.df, id = "faostat")
aclass.df$variable <- substr(aclass.df$variable, 3, 6)
colnames(aclass.df) <- c("faostat", "year", "FD")

## Extract the non-A component
nonafaostat.df <- faostat.df[colnames(faostat.df)
                             [-grep("A", colnames(faostat.df))]]
nonaclass.df <- melt(nonafaostat.df, id = "faostat")
colnames(nonaclass.df) <- c("faostat", "class", "value")

## sum up the variables by year
nonasum.df <- aggregate(nonaclass.df$value,
                        list(nonaclass.df$faostat,
                             substr(nonaclass.df$class, 3, 6)),
                        sum)
colnames(nonasum.df) <- c("faostat", "year", "sum")

## merge back with the A-class and the take the difference
final.df <- merge(aclass.df, nonasum.df, by = c("faostat", "year"))
final.df$diff <- final.df$FD - final.df$sum

######################################################################
## Testing NPMLE
######################################################################

mix.efron = normmix(mu=c(-10.9, -7.0, -4.9, -1.8, -1.1, 0.0, 2.4, 6.1),
                    pi=c(1.5, 1.3, 5.6, 12.3, 13.6, 60.8, 2.7, 2.2))

set.seed(1)
x = rnormmix(n=1000, mix=mix.efron)
cnm.normmix(x, tol=1e-5)


######################################################################
## Checking expression
######################################################################

y <- parse(text = "x^2")

plot(1)
text(0.9, 0.9, y)

a <- "perc of CO[2]"

a <- "CO[2]"
expr <- eval(parse(text = paste("quote(", a, ")", sep = "")))
exprtxt <- "Perc of "
text<-"This is the text of the legend:"
year<- "The year is 2000"

myLab <- substitute(paste(text, exprtxt, expr, year, sep = " "),
                    list(text = text, exprtxt = exprtxt,
                         expr = expr, year = year))

plot(1)
text(0.9, 0.9, labels = myLab)

text(1.1, 1.1, labels =
     bquote(paste(.(text), .(expr), .(year), sep = " ")))


## Break character text into expression and txt
split <- unlist(strsplit("Perc of CO[2] produced", " "))

splitIndex <- which(attr(regexpr("[[:punct:]]", split),
                         "match.length") == 1)

## Break the text into 3 chunks
## The expression chunk
expr <- eval(parse(text = paste("quote(", split[splitIndex], ")",
                   sep = " ")))

## The text before expression
if(splitIndex != 1){
    startTxt <- paste(paste(split[c(1:(splitIndex - 1))],
                            collapse = " "), " ", sep = "")
}else{
    startTxt <- ""
}
## The text after expression
if(splitIndex != length(split)){
    endTxt <- paste(" ", paste(split[c((splitIndex + 1):length(split))],
                               collapse = " "), sep = "")
}else{
    endTxt <- ""
}

plot(1, type = "n")
text(1, 1, labels = bquote(paste(.(startTxt), .(expr), .(endTxt),
               sep = "")))


######################################################################
## Investigate ylim issue
######################################################################


DF <- data.frame(x = rnorm(100, mean = 5),
                 y = sample(letters[1:2], size=100, replace = TRUE))

qplot(x=y, y=x, data=DF, geom="bar", stat="identity")

qplot(x = y, y = x, data = DF, geom = "bar", stat = "identity",
      ylim = c(100, 300))


######################################################################
## To test to use if or swith
######################################################################

istat <- function(char){
    if(char == "a") {print("a")}
    else if(char == "b") {print("b")}
    else if(char == "c") {print("c")}
}

sstat <- function(char){
    switch(char,
          a = print("a"),
          b = print("b"),
          c = print("c")
      )
}

library(rbenchmark)
mychar <- "c"

benchmark(istat(mychar), sstat(mychar), replications = 1e5,
          order = "relative")

######################################################################
## Is wide or long format more efficient in terms of storing data
######################################################################

## Clearly shows that the data should be stored as long format.
l.mat <- matrix(rnorm(8e7), nc = 8)
object.size(l.mat)
l.df <- data.frame(l.mat)
object.size(l.df)

w.mat <- matrix(rnorm(1e7), nc = 10000)
object.size(w.mat)
w.df <- data.frame(w.mat)
object.size(w.df)



######################################################################
## Neural network example
######################################################################

###
### prepare data
###
library(mlbench)
data(BostonHousing)

# inspect the range which is 1-50
summary(BostonHousing$medv)


##
## model linear regression
##

lm.fit <- lm(medv ~ ., data=BostonHousing)

lm.predict <- predict(lm.fit)

# mean squared error: 21.89483
mean((lm.predict - BostonHousing$medv)^2)

plot(BostonHousing$medv, lm.predict,
    main="Linear regression predictions vs actual",
    xlab="Actual")


##
## model neural network
##
require(nnet)

# scale inputs: divide by 50 to get 0-1 range
nnet.fit <- nnet(medv/50 ~ ., data=BostonHousing, size=2)

# multiply 50 to restore original scale
nnet.predict <- predict(nnet.fit)*50

# mean squared error: 16.40581
mean((nnet.predict - BostonHousing$medv)^2)

plot(BostonHousing$medv, nnet.predict,
    main="Neural network predictions vs actual",
    xlab="Actual")




######################################################################
## Test the caret package and random forest
######################################################################

library(caret)
library(randomForest)
## TODO: Check what foreach is


######################################################################
## Passing arguements
######################################################################

## Don't work
x.lm <- function(formula, data, ...)
{
  eval(substitute(lm(f, data), list(f = formula)))
}

## Work
x.lm <- function(formula, data, ...)
{
    Call <- match.call(expand.dots = TRUE)
    Call[[1]] <- as.name("lm")
    print(names(Call))
    Call$formula <- as.formula(terms(formula))
    eval(Call)
}
x.lm(formula = formula("y1 ~ x1"), data = anscombe, model = FALSE)

######################################################################
## Parallel computin
######################################################################
library(multicore)
require(foreach) #loads package foreach
require(doMC) #loads both doMC and multicore
search()   #make sure all 3 packages are loaded
multicore:::detectCores(all.tests=TRUE)  #to see how many cores are available


######################################################################
## Comparison between matrices
######################################################################

test.mat1 <- matrix(1:20, nc = 5)

test.mat2 <- rbind(test.mat1[sample(1:5, 2), ], matrix(101:120, nc = 5))

compMat <- function(mat1, mat2){
    nr1 <- nrow(mat1)
    nr2 <- nrow(mat2)
    mat2[duplicated(rbind(mat1, mat2))[(nr1 + 1):(nr1 + nr2)], ]
}

compMat(test.mat1, test.mat2)




######################################################################
## POSIXct time stamp
######################################################################

x <- 1472562988 + 1:10; tz <- rep("EST",10)

# Case 1: Works as documented
ct <- as.POSIXct(x, tz=tz[1], origin="1960-01-01")

# Case 2: Fails
ct <- as.POSIXct(x, tz=tz, origin="1960-01-01")

##If case 2 worked, it'd be a little easier to process paired (time,
##time zone) vectors from different time zones.


######################################################################
## SQL script for SYB
######################################################################


#SELECT * FROM warehouse.dbo.data
#  WHERE DomainCode='QC'
#  AND ItemCode=1754
#  AND ElementCode=5510



setwd("c:/Documents and Settings/Kao/Desktop/Dropbox/AdaMat/SYB")


library(RJDBC)
drv <- JDBC("com.microsoft.sqlserver.jdbc.SQLServerDriver",
            "sqljdbc4.jar")
conn <- dbConnect(drv, "jdbc:sqlserver://FAOSTAT-PROD\\Production",
                  "Warehouse", "w@reh0use")

P3.FEED.FAO.ESS.GPCPIN.CRPS.df <-
    dbGetQuery(conn, "SELECT * FROM warehouse.dbo.data WHERE DomainCode = 'QI' AND ItemCode =  2041 AND ElementCode = 434 AND AreaCode = 5000")

P3.FEED.FAO.ESS.GPCPIN.LSTK.df <-
    dbGetQuery(conn, "SELECT * FROM warehouse.dbo.data WHERE DomainCode = 'QI' AND ItemCode =  2044 AND ElementCode = 434 AND AreaCode = 5000")

P3.FEED.FAO.ESS.GPCPIN.FOOD.df <-
    dbGetQuery(conn, "SELECT * FROM warehouse.dbo.data WHERE DomainCode = 'QI' AND ItemCode =  2054 AND ElementCode = 434 AND AreaCode = 5000")

P3.FEED.FAO.ESS.GPCPIN.NFOOD.df <-
    dbGetQuery(conn, "SELECT * FROM warehouse.dbo.data WHERE DomainCode = 'QI' AND ItemCode =  2057 AND ElementCode = 434 AND AreaCode = 5000")


P4.ENV.FAO.CC.CE.QP.df <-
    dbGetQuery(conn, "SELECT * FROM warehouse.dbo.data WHERE DomainCode = 'QC'AND ItemCode = 1717 AND ElementCode = 5510")


## Manipulate the data
CRPS <- P3.FEED.FAO.ESS.GPCPIN.CRPS.df[c("AreaCode", "Year", "Value")]
LSTK <- P3.FEED.FAO.ESS.GPCPIN.LSTK.df[c("AreaCode", "Year", "Value")]
FOOD <- P3.FEED.FAO.ESS.GPCPIN.FOOD.df[c("AreaCode", "Year", "Value")]
NFOOD <- P3.FEED.FAO.ESS.GPCPIN.NFOOD.df[c("AreaCode", "Year", "Value")]

colnames(CRPS)[3] <- "CRPS"
colnames(LSTK)[3] <- "LSTK"
colnames(FOOD)[3] <- "FOOD"
colnames(NFOOD)[3] <- "NFOOD"

GPCPIN.df <- Reduce("merge", list(CRPS, LSTK, FOOD, NFOOD))

write.csv(GPCPIN.df,
          "c:/Documents and Settings/Kao/Desktop/Dropbox/Adamat/SYB/Data/EXT.DATA/P3.FEED.FAO.ESS.GPCPIN.AGG.csv")

cp <- P4.ENV.FAO.CC.CE.QP.df[c("AreaCode", "Year", "Value")]

write.csv(cp, "c:/Documents and Settings/Kao/Desktop/Dropbox/Adamat/SYB/Data/EXT.DATA/P4.ENV.FAO.CC.CE.QP.csv")

write.csv(CRPS, "c:/Documents and Settings/Kao/Desktop/Dropbox/Adamat/SYB/Data/EXT.DATA/P3.FEED.FAO.ESS.GPCPIN.CRPS.csv")




    P3.FEED.FAO.ESS.FD.IXv.df <-
        dbGetQuery(conn, "SELECT * FROM warehouse.dbo.data WHERE DomainCode = 'TP' AND ItemCode =  1982 AND ElementCode = 5922")

P3.FEED.FAO.ESS.FD.IXv.df <-
    P3.FEED.FAO.ESS.FD.IXv.df[c("AreaCode", "Year", "Value")]
write.csv(P3.FEED.FAO.ESS.FD.IXv.df, "c:/Documents and Settings/Kao/Desktop/Dropbox/Adamat/SYB/Data/EXT.DATA/P3.FEED.FAO.ESS.FD.IXv.csv")



######################################################################
## parallel computing
######################################################################


library(parallel)
cl <- makeCluster(3)
parLapply(cl, 1:3, sqrt)
stopCluster(cl)


########################################################################
## Skewness in S&P 500
########################################################################

library(zoo)
library(tseries)
source("http://www.portfolioprobe.com/R/blog/kurtskew.R")
spx.skew250 <- rollapply(spxret, 250, pp.skew, align='right')

########################################################################
## Parallel R Loops for Windows and Linux
########################################################################

library(foreach)
library(doMC)
registerDoMC(4)

## No parallel
system.time(foreach(i=1:10) %do% {i})
## Paralleled
system.time(foreach(i=1:10) %dopar% {i})

foreach(i=letters[1:26], .combine = c) %dopar% {i}


## Output as column matrix
col.mat <- foreach(i=1:10, .combine = cbind) %do% { i }

## Output as row matrix
row.mat <- foreach(i=1:10, .combine = rbind) %do% { i }

## Output as duplicated matrix
rc.mat <- foreach(i=1:10, .combine = rbind) %do% { c(i,i) }


## Test fibonacci sequence
n <- 5000
np <- double(n)
np[1] = 0
np[2] = 1
## The original loop
system.time(
for(i in 3:n){
  np[i] = np[i - 1] + np[i - 2]
})
## The foreach loop
system.time(foreach(i = 3:n, .combine = rbind) %do% {
  np[i] = np[i - 1] + np[i - 2]}
)
## The paralleled foreach loop
system.time(foreach(i = 3:n, .combine = rbind) %dopar% {
  np[i] = np[i - 1] + np[i - 2]}
)

n = 1e3
system.time({
    test = foreach(i = 1:n, .combine = c) %dopar% {
        tmp = rnorm(100)
        tmp2 = mean(tmp)
    }
})

system.time({
    test = for(i in 1:n){
        tmp = rnorm(100)
        tmp2 = mean(tmp)
    }
})





#########################################################################
## Bagging and basic ensemble model
########################################################################

## Generate random numbers
set.seed(10)
y <- c(1:1000)
x1 <- c(1:1000) * runif(1000,min = 0,max = 2)
x2 <- c(1:1000) * runif(1000,min = 0,max = 2)
x3 <- c(1:1000) * runif(1000,min = 0,max = 2)

## break the data into train and test set
all_data <- data.frame(y, x1, x2, x3)
positions <- sample(nrow(all_data), size = floor((nrow(all_data)/4)*3))
training <- all_data[positions,]
testing <- all_data[-positions,]

## Build the model
lm_fit <- lm(y ~ x1 + x2 + x3,data = training)
predictions <- predict(lm_fit, newdata = testing)
error <- sqrt((sum(testing$y - predictions)^2)/nrow(testing))

## Build the model with bagging
library(foreach)
length_divisor <- 4
iterations <- 1000
system.time(
foreach(m = 1:iterations, .combine=cbind) %do% {
  training_positions <-
    sample(nrow(training), size = floor((nrow(training)/length_divisor)))
  train_pos <- 1:nrow(training) %in% training_positions
  lm_fit <- lm(y ~ x1 + x2 + x3, data = training[train_pos,])
  predict(lm_fit, newdata = testing)
}
predictions <- rowMeans(predictions)
error <- sqrt((sum(testing$y - predictions)^2)/nrow(testing))
)

## Function for bagging
bagging<-function(training, testing, length_divisor = 4,
                  iterations = 1000){
  predictions <- foreach(m = 1:iterations, .combine = cbind) %do% {
    training_positions <- sample(nrow(training),
                                 size = floor((nrow(training)/length_divisor)))
    train_pos <- 1:nrow(training) %in% training_positions
    lm_fit <- lm(y ~ x1 + x2 + x3, data = training[train_pos,])
    predict(lm_fit,newdata=testing)
  }
  rowMeans(predictions)
}


library(randomForest)
rf_fit <- randomForest(y ~ x1 + x2 + x3, data = training, ntree = 500)
predictions <- predict(rf_fit, newdata = testing)
error <- sqrt((sum((testing$y-predictions)^2))/nrow(testing))
error


length_divisor <- 6
iterations <- 5000
predictions <- foreach(m=1:iterations,.combine=cbind) %do% {
  training_positions <-
    sample(nrow(training), size = floor((nrow(training)/length_divisor)))
  train_pos <- 1:nrow(training) %in% training_positions
  lm_fit <- lm(y ~ x1 + x2 + x3, data = training[train_pos,])
  predict(lm_fit, newdata = testing)
}
lm_predictions <- rowMeans(predictions)

library(randomForest)
rf_fit <- randomForest(y ~ x1 + x2 + x3,data = training, ntree = 500)
rf_predictions <- predict(rf_fit, newdata = testing)
predictions <- (lm_predictions+rf_predictions)/2
error <- sqrt((sum((testing$y-predictions)^2))/nrow(testing))
error

predictions<-(lm_predictions+rf_predictions*9)/10
error<-sqrt((sum((testing$y-predictions)^2))/nrow(testing))
error

library(e1071)
svm_fit <- svm(y ~ x1 + x2 + x3, data = training)
svm_predictions <- predict(svm_fit, newdata = testing)
error<-sqrt((sum((testing$y-svm_predictions)^2))/nrow(testing))
error

length_divisor <- 6
iterations <- 5000
predictions <- foreach(m=1:iterations,.combine=cbind) %do% {
  training_positions <-
    sample(nrow(training), size = floor((nrow(training)/length_divisor)))
  train_pos<-1:nrow(training) %in% training_positions
  svm_fit <- svm(y ~ x1 + x2 + x3, data = training[train_pos,])
  predict(svm_fit,newdata=testing)
}
svm2_predictions <- rowMeans(predictions)
error<-sqrt((sum((testing$y-svm2_predictions)^2))/nrow(testing))
error

predictions <- (svm_predictions+rf_predictions)/2
error <- sqrt((sum((testing$y-predictions)^2))/nrow(testing))
error

predictions <- (svm_predictions*2+rf_predictions)/3
error <- sqrt((sum((testing$y-predictions)^2))/nrow(testing))
error


########################################################################
## Facebook API access
########################################################################

## Obtain your token from Facebook graphic API
## (https://developers.facebook.com/tools/explorer?method=GET&path=543951338)
myFBToken <- "AAACEdEose0cBAIrCJDsrKrAdIrqbJmOz3BHGMuCtvvswzrPbhfX78tdc0XqmDvZBDcAoXoTK3eSPQrSCxWLWUkI8yBWsZAkh2tIVvwZBAZDZD"

## Function from Romain Francois to scrap photo url using the facebook
## graphic API

require(RCurl)
require(rjson)

facebook <-  function( path = "me", access_token = myFBToken, options){
  if( !missing(options) ){
    options <- sprintf( "?%s",
                       paste( names(options), "=", unlist(options),
                             collapse = "&", sep = "" ) )
  } else {
    options <- ""
  }
  data <- getURL(sprintf("https://graph.facebook.com/%s%s&access_token=%s",
                         path, options, access_token ))
    fromJSON(data)
}

## Download the photos
dir.create( "FBphotos" )
photos <- facebook( "me/photos/")
sapply( photos$data, function(x){
    url <- x$source
    download.file( url, file.path( "FBphotos", basename(url) ) )
})



########################################################################
## Testing SQL query
########################################################################

test7 <- data.frame(g = letters[c(1, 1, 2, 2, 3, 3)], v = c(1, 2, 1, 3, 4, 1))
check <- sqldf("SELECT * FROM test7 GROUP BY g HAVING v = max(v)")

########################################################################
## Dirk's ls()
########################################################################

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}
# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

######################################################################
## Function equivalent to trim in SQL, removing space in front or
## training of a character string
######################################################################

trim <- function(x){
  gsub("^[[:space:]]+|[[:space:]]+$", "", x)
}

########################################################################
## Write a function to forecast a repetitive sequence
########################################################################

seq.ts <- c(rep(c(1, 9), 3), 5)

seqForecast <- function(ts, n.ahead){
  T = length(ts)
  for(i in 3:T){
    if(ts[i] != ts[i - 2]){
      ts[i] = ts[i - 2]
      warning("Data point ", i, " is not a typical data point, please check")
    }
  }
  pred = rep(ts[(T - 1):(T - 2)], n.ahead/2)
  pred
}

########################################################################
## Example of resistant regression
########################################################################

## The idea is to trim off the tail of observations, so the estimate
## obtain are not affected by these potential outliers


library(quantmod)
end <- format(Sys.Date(),"%Y-%m-%d")
start <- format(as.Date("2003-09-18"),"%Y-%m-%d")
dat0 = as.matrix(getSymbols("AAPL", src="yahoo", from = start,
  to = end, auto.assign = FALSE))
dat1 = as.matrix(getSymbols("spy", src = "yahoo", from = start,
  to = end, auto.assign = FALSE))
n = NROW(dat0)
ret = (dat0[2:n,6]/dat0[1:(n-1),6] - 1)
ret_spy = (dat1[2:n,6]/dat1[1:(n-1),6] - 1)

##  Resistant Regression function: mltsreg = modified least trimmed
##  regression
mltsreg <- function(x = ret, y = ret_spy, k){
  trimmedx = x[x < quantile(x,k) & x > quantile(x,1 - k)]
  trimmedy = y[x < quantile(x,k) & x > quantile(x,1 - k)]
  MtrimmedReg = lm(trimmedx ~ trimmedy)
  MtrimmedReg
}

## plot:
plot(ret ~ ret_spy,
     main = "AAPL Over SPY - Resistant Regression, Fitted Values",
     xlab = "Market Returns", ylab = " Returns")
s = seq(.8, .95, 0.05)
for.legend = NULL
col1 = c(1:length(s))
for(i in 1:length(s)){
  abline(mltsreg(k = s[i])$coef[1:2], col = col1[i], lwd = 1.5)
  for.legend[i] = c(paste(signif(1 - s[i], 2), "Percent Trimmed, ",
              "Beta =" ,paste(signif(mltsreg(k = s[i])$coef[2],2) )))
}
abline(lm(ret~ret_spy)$coef[1:2], col = 'gold', lwd = 3)
legend("bottomright", c(for.legend,'OLS, Beta = 0.97'), bty = "n",
       col = c(col1,'gold'), lty = c(1), lwd = c(rep(1.5,5,),3))


########################################################################
## Determining machine episilon (Ross Ihaka)
########################################################################

## R uses IEEE 754 Arithmetic standard
macheps =
  function(){
    {
      eps = 1
      while(1 + eps/2 != 1)
        eps = eps/2
      eps
    }
  }
macheps()

## Cancellation error, the subtraction of nearly equal number is a major
## source of inaccuracy in numerical error.


########################################################################
## Write a function for calculating least squares growth rate
########################################################################

lsgr <- function(x){
  ## Test function for generating the least squares growth rate
  ##
  ## Args:
  ##   x: The time series for the growth rate to be calculated
  ##
  ## Returns:
  ##   The least squares growth rate of the time series
  ##
  ## NOTES:
  ##   Missing values are omitted.
  ##
  ## Error handling:
  if(sum(is.na(x)) > 0.5 * length(x))
    stop("Over 50% of the data are missing")

  t <- 1:length(x)
  lsgr <- (exp(coef(lm(log(x) ~ t))[2]) - 1) * 100
  lsgr
}

x <- c(abs(rnorm(30)))
misInd <- sample(c(1:30), 5)
x[misInd] <- NA

lsgr(x)


x <- c(1243020281.54778, 1351538833.70195, 1252430786.64695,
       1373814378.00948, 1233133288.98527, 1439421907.48375,
       1544739823.87934, 1603322832.32072, 1665776322.98803,
       1716453871.17722, 1743506520.28913, 1656767701.30891,
       1684932752.23116, 1777604053.60388, 1875372276.55209,
       1924131955.74244, 2072290116.33461, 2192482943.08202,
       2381036476.18707, 2561995248.37729)
myGr <- lsgr(x)

myVal <- x[1] * (1 + myGr/100)^(1:20)
pdf("lsgr.pdf")
plot(x, type = "b")
lines(myVal, col = "red")
points(myVal, col = "red")
dev.off()

agMean <- function(x){
  T <- length(x)
  start <- 1
  end <- x[T]/x[1]
  a <-(start + end)/T
  g <- sqrt(start * end)^(1/T)
  while(a != g){
    a = 0.5 * (a + g)
    g = sqrt(a * g)
  }
  a
}

agMean(x)




########################################################################
## Some map projection plots
########################################################################

## Load required packages
require(akima)
require(maps)
require(mapproj)

## Loading functions from "www.menugget.blogspot.com"
source("menugget.R")

## create xyz data
set.seed(1)
mecca <- data.frame(lon=39.826138, lat=21.422508)
x <- runif(10000, min=-180, max=180)
y <- runif(10000, min=-90, max=90)
z <- earth.dist(x, y, mecca$lon, mecca$lat)
xyz_data <- cbind(x,y,z)

###interpolate data
INTERP <- interp(x, y, z, xo=seq(-180, 180, length = 100), yo=seq(-90, 90, length = 100))

###basic plot
image(INTERP, col=rainbow(100))
map("world",add=TRUE, col=1)

###projection plots
#make polygons
polys<-vector("list", length(as.vector(INTERP$z)))
for(i in 1:length(polys)){
 lonx <- pos2coord(pos=i, dim.mat=dim(INTERP$z))[1]
 laty <- pos2coord(pos=i, dim.mat=dim(INTERP$z))[2]
 ifelse(laty < length(INTERP$y), neigh_y<-c(laty+1,lonx), neigh_y<-c(laty-1,lonx))
 ifelse(lonx < length(INTERP$x), neigh_x<-c(laty,lonx+1), neigh_x<-c(laty,lonx-1))
 dist_y <- earth.dist(INTERP$x[lonx], INTERP$y[laty], INTERP$x[neigh_y[2]], INTERP$y[neigh_y[1]])
 dist_x <- earth.dist(INTERP$x[lonx], INTERP$y[laty], INTERP$x[neigh_x[2]], INTERP$y[neigh_x[1]])
 s1 = new.lon.lat(INTERP$x[lonx], INTERP$y[laty], 180, dist_y/2)
 s3 = new.lon.lat(INTERP$x[lonx], INTERP$y[laty], 0, dist_y/2)
 s2 = new.lon.lat(INTERP$x[lonx], INTERP$y[laty], 270, dist_x/2)
 s4 = new.lon.lat(INTERP$x[lonx], INTERP$y[laty], 90, dist_x/2)
 polys[[i]] = cbind(c(s2[1], s2[1], s4[1], s4[1]), c(s1[2], s3[2], s3[2], s1[2]))
}

#values to colors
pal <- color.palette(rev(c("purple","black","red","yellow")), n.steps.between=c(2,5,20))
COL <- val2col(INTERP$z, col=pal(100))

#plot parameters
WIDTHS=c(4,4,1)
HEIGHTS=c(4)
OMI=c(0.1,0.1,0.1,0.1)

#define plot device and layout
png("xyz_2_projections.png", width=sum(WIDTHS)+sum(OMI[c(2,4)]), height=sum(HEIGHTS)+sum(OMI[c(1,3)]), units="in",res=400)
#x11(width=sum(WIDTHS)+sum(OMI[c(2,4)]), height=sum(HEIGHTS)+sum(OMI[c(1,3)]))
layout(matrix(c(1:3),nrow=1,ncol=3), width=WIDTHS, height=HEIGHTS, respect=TRUE)
layout.show(3)

#1st plot - perspective projection
XLIM=c(-180,180)
YLIM=c(-90,90)
PROJ="perspective"
ORIENT=c(mecca$lat,mecca$lon,0)
PAR=2
par(mai=c(0.1,0.1,0.5,0.1))
map("world", projection=PROJ, orientation=ORIENT, par=PAR, xlim=XLIM, ylim=YLIM, col=NA)
for(i in 1:length(polys)){
 polygon(mapproject(x=polys[[i]][,1], y=polys[[i]][,2]), col=COL[i], border=COL[i], lwd=0.1)
}
map("world", add=TRUE, projection="", xlim=XLIM, ylim=YLIM, col="grey")
map.grid(c(-180, 180, -80, 80), nx=12, ny=6, labels=FALSE, col="grey")
mtext("Perspective projection", side=3, line=1)

#2nd plot - stereographic projection
XLIM=c(-180,180)
YLIM=c(-90,90)
PROJ="stereographic"
ORIENT=c(mecca$lat,mecca$lon,0)
PAR=NULL
par(mai=c(0.1,0.1,0.5,0.1))
map("world", projection=PROJ, orientation=ORIENT, par=PAR, xlim=XLIM, ylim=YLIM, col=NA)
for(i in 1:length(polys)){
 polygon(mapproject(x=polys[[i]][,1], y=polys[[i]][,2]), col=COL[i], border=COL[i], lwd=0.1)
}
map("world", add=TRUE, projection="", xlim=XLIM, ylim=YLIM, col="grey")
map.grid(c(-180, 180, -80, 80), nx=12, ny=6, labels=FALSE, col="grey")
mtext("Stereographic projection", side=3, line=1)

#3rd plot - color scale
par(mai=c(0.1,0.1,0.5,0.7))
image.scale(INTERP$z, col=pal(100), horiz=FALSE, xaxt="n", yaxt="n", xlab="", ylab="")
box()
axis(4)
mtext("Distance to Mecca [km]", side=4, line=3)

#close device
dev.off()


########################################################################
## Fit distribution and generate the mixture distribution.
########################################################################

library(MASS)

## Generate three sample of log-normal distribution
d1 <- rlnorm(100, 1, 0.5)
d2 <- rlnorm(100, 2, 0.5)
d3 <- rlnorm(100, 3, 0.5)

## Fit the distribution
coef1 <- coef(fitdistr(d1, densfun = "log-normal"))
coef2 <- coef(fitdistr(d2, densfun = "log-normal"))
coef3 <- coef(fitdistr(d3, densfun = "log-normal"))

# Initialise parameters, change the n if you want a different sample size.
n <- 10000
mysamp <- double(n)

## Determine which sample to draw from
u <- runif(n)

## Draw the sample
mix <- table(as.numeric(cut(u, c(0, 0.4, 0.8, 1))))
mysamp <- c(rlnorm(mix[1], meanlog = coef1["meanlog"], sdlog = coef1["sdlog"]),
            rlnorm(mix[2], meanlog = coef2["meanlog"], sdlog = coef2["sdlog"]),
            rlnorm(mix[3], meanlog = coef3["meanlog"], sdlog = coef3["sdlog"]))
hist(mysamp, breaks = 500)



########################################################################
## Bootstrap fit a bivariate density
########################################################################

rv <- rmvst(1000, dim = 2)
colnames(rv) <- c("x", "y")
plot(rv)

## Parametric: I need more thought on this, you will need to integrate
##             the density to find the highest associated x value
##             corresponding to the 5th percentile of y.
mvFit(rv)


## Non-parametric

x <- rnorm(1000)
y <- 2 + 2*x + rnorm(1000)
rv <- cbind(x, y)
plot(rv)

biDen <- function(x, y, alpha, n.sim){
    raw <- cbind(x, y)
    est <- matrix(0, nc = 2, nr = n.sim)
    for(i in 1:n.sim){
        ## Take a random sample from the data and sort
        ind <- sample(1:NROW(raw), NROW(raw) * 0.3, replace = TRUE)
        samp <- raw[ind, ]
        samp <- data.frame(samp[order(samp[, 2]), ])

        ## Find the value of x associated with 5th percentile of y
        ind2 <- NROW(samp) * alpha
        est[i, 1] <-samp[ind2, 1]
        rst <- lm(samp[, 2] ~samp[, 1], data = samp)
        est[i, 2] <- predict(rst,
                             newdata = data.frame(quantile(samp[, 1], 0.05)))
    }
    ## Sort the estimate
    ## The mean is the bootstrapped estimate of the value of x most
    ## likely to correspond to the 5th percentile of y and the
    ## confidence interval of it

    ##
    est
    ## list(mean = mean(est), lowerci = est[n.sim * 0.05],
    ##      upperci = est[n.sim * 0.95])
}

alpha = 0.05
x5pct <- biDen(rv[, 1], rv[, 2], alpha = alpha, n.sim = 1000)
plot(rv)
abline(h = quantile(rv[, 2], alpha), col = "red", lty = 2)
abline(v = mean(x5pct[1, ]), col = "red", lty = 2)
abline(v = mean(x5pct[2, ]), col = "red", lty = 2)

abline(v = x5pct$lowerci, col = "orange", lty = 2)
abline(v = x5pct$upperci, col = "orange", lty = 2)
abline(v = quantile(rv[, 1], alpha), col = "blue", lty = 2)



########################################################################
## Sending email with R
########################################################################

library("sendmailR")


from <- sprintf("<Project1@%s>", Sys.info()[4])
to <- "<kobegoya@hotmail.com>"

from <- "<mkao006@gmail.com>"
subject <- "Test Email From R"

#create list with text of body as first element
#second list element is R object to attach using the mime_part() function
body <- list("Email sent from R. Dataframe attached.",mime_part(data.frame(x=rnorm(1000),y=rnorm(1000)),name="output"))

sendmail(from, to, subject, body, control=list(smtpServer="ASPMX.L.GOOGLE.COM"))


########################################################################
## testing some theory
########################################################################



## Testing the variance of a mixture
n1 = 10000
n2 = 10000
N = n1 + n2
s1 = 3
s2 = 5
a <- rnorm(n1, 0, s1)
b <- rnorm(n2, 0, s2)
test <- c(a, b)
sd(test)
sqrt((n1/N)^2 * s1^2 + (n2/N)^2 * s2^2)
sqrt((n1/N) * s1^2 + (n2/N) * s2^2)

## Looks like not squaring the weight is the correct formula...

## Testing the skewed version

## Testing the variance of a mixture
n1 = 10000
n2 = 10000
N = n1 + n2
s1 = 0.2
s2 = 0.8
xi = 3
a <- rsnorm(n1, 0, s1, xi = xi)
b <- rsnorm(n2, 0, s2, xi = xi)
test <- c(a, b)
hist(test, breaks = 100, freq = FALSE)
sd(test)
sqrt((n1/N) * s1^2 + (n2/N) * s2^2)


curve(.dsnorm(x, xi = 2), -3, 3)
curve(.dsnorm2(x, xi = 2), add = TRUE, col = "red", lty = 2)

disc <- function(pr, pt){
    list(pr = pr, pt = pt)
}

dsnorm2 <- function(x, mean = 0,  sd = 1, xi = 1.5){
    .dsnorm2(x, xi = xi)
}

.dsnorm2 <- function (x, xi){
    m1 = 2/sqrt(2 * pi)
    mu = m1 * (xi - 1/xi)
    sigma = sqrt((1 - m1^2) * (xi^2 + 1/xi^2) + 2 * m1^2 - 1)
    z = x * sigma + mu
    Xi = xi^sign(z)
    g = 2/(xi + 1/xi)
    Density = g * dnorm(x = z/Xi)
    Density * sigma
}

.rsnorm2 <- function (n, xi){
    weight = xi/(xi + 1/xi)
    z = runif(n, -weight, 1 - weight)
    Xi = xi^sign(z)
    Random = -abs(rnorm(n))/Xi * sign(z)
    Random
}

.rsnorm <- function (n, xi)
{
    weight = xi/(xi + 1/xi)
    z = runif(n, -weight, 1 - weight)
    Xi = xi^sign(z)
    Random = -abs(rnorm(n))/Xi * sign(z)
    m1 = 2/sqrt(2 * pi)
    mu = m1 * (xi - 1/xi)
    sigma = sqrt((1 - m1^2) * (xi^2 + 1/xi^2) + 2 * m1^2 - 1)
    print(sigma)
    Random = (Random - mu)/sigma
    Random
}

test <- .rsnorm(10000, xi = 3)
test2 <- .rsnorm2(10000, xi = 3)
par(mfrow = c(1, 2))
hist(test, breaks = 100)
hist(test2, breaks = 100)

dsnorm2 <- function(x, xi){
    m1 = 2/sqrt(2 * pi)
    mu = m1 * (xi - 1/xi)
    sigma = sqrt((1 - m1^2) * (xi^2 + 1/xi^2) + 2 * m1^2 - 1)
    z = x * sigma + mu
    Xi = xi * Heaviside(z) + 1/xi * Heaviside(-z)
    g = 2/(xi + 1/xi)
    Density = g * dnorm(x = z/Xi)
    Density * sigma
}

curve(.dsnorm(x, xi = 3), -3, 3)
curve(dsnorm2(x, xi = 3), add = TRUE, col = "red")

ll <- function(x, theta, sigma, xi){
    mu = theta * sigma * sqrt(2/pi) * (xi - 1/xi)
    sd = sqrt(theta^2 * sigma^2 * ((1 - 2/pi) * (xi^2 + 1/xi^2)+ 4/pi - 1))
    z = x * sd + mu
    log(2) + log(sd) + log(xi + 1/xi) - 0.5 * log(2 * pi) - log(theta) - log(sigma) -
            (z * (xi * Heaviside(z) + 1/xi * Heaviside(-z)))^2/(2 * theta^2 * sigma^2)
}

inc <- 1e-5
sigma = 2
theta = 1
xi = 2
#curve(ll(x, theta = 1, sigma = 1, xi = 1), -3, 3)
curve((ll(x, theta = 1 + inc, sigma = 1, xi = 1) -
      ll(x, theta = 1 - inc, sigma = 1, xi = 1))/(2 * inc), -100, 100,
      ylim = c(-100, 100))
abline(h = 0, col = "red")
abline(h = , col = "red")
2 * sigma * ((1 - 2/pi) * (xi^2 + 1/xi^2) + 4/pi - 1) - 1/theta
2 * theta  * ((1 - 2/pi) * (xi^2 + 1/xi^2) + 4/pi - 1) - 1/sigma

curve(dsnorm(x, 0, 1.1, 1.5) + dsnorm(x, 0, 0.2, 1.5), -3, 3)

curve(((2/(xi + 1/xi)) * (dnorm(xi * x) * Heaviside(-x) + dnorm(x/xi) * Heaviside(x))), -3, 3)
curve(dsnorm(x, sd = 1.35, xi), add = TRUE, col = "red")



dmnorm <- function(x, varmix = disc(1, 1)){
  n = length(x)
  n.mix = length(varmix$pt)
  d = 1/sqrt(2 * pi * matrix(rep(varmix$pt^2, each = n), nc = n.mix)) *
    exp(-0.5 * x^2/matrix(rep(varmix$pt^2, each = n), nc = n.mix))
  d %*% varmix$pr
}

dmsnorm <- function(x, varmix = disc(1, 1), xi = 1, demean = TRUE){
    z = x * (xi * Heaviside(-x) + 1/xi * Heaviside(x))
    (2/(xi + 1/xi)) * dmnorm(z, varmix = varmix)
}

dmnorm2 <- function(x, varmix = disc(1, 1), xi = 1){
  n = length(x)
  n.mix = length(varmix$pt)
  d = 1/sqrt(2 * pi * matrix(rep(varmix$pt^2, each = n), nc = n.mix)) *
      exp(-0.5 * (x - c((sqrt(2 * varmix$pt^2/pi) %*% varmix$pr) * (xi - 1/xi)))^2/
          matrix(rep(varmix$pt^2, each = n), nc = n.mix))
  d %*% varmix$pr
}


dmsnorm2 <- function(x, varmix = disc(1, 1), xi = 1, demean = TRUE){
    z = x * (xi * Heaviside(-x) + 1/xi * Heaviside(x))
    (2/(xi + 1/xi)) * dmnorm2(z, varmix = varmix)
}

rmnorm <- function(n, varmix = disc(1, 1)){
    w <- runif(n)
    breaks = cumsum(varmix$pr)
    m <- double(n)
    for(i in 1:n){
        m[i] = sum(w[i] > breaks) + 1
    }
    rnorm(n, sd = varmix$pt[m])
}


rmsnorm <-function(n, varmix = disc(1, 1), xi = 1){
    weight = xi/(xi + 1/xi)
    z = runif(n, -weight, 1 - weight)
    Xi = xi^sign(z)
    Random = -abs(rmnorm(n, varmix = varmix))/Xi * sign(z)
    m1 = c(sqrt(2 * varmix$pt^2/pi) %*% varmix$pr)
    m2 = sum(varmix$pr * varmix$pt^2)
    mu = m1 * (xi - 1/xi)
    sigma = sqrt((m2 - m1^2) * (xi^2 + 1/xi^2) + 2 * m1^2 - m2)
#    Rand = (Random - mu)/sigma
    Random
}

xi = 2
n = 100000
test <- rmnorm(n, varmix = disc(pt = c(1, 2), pr = c(0.5, 0.5)))
hist(test, breaks = 300, freq = FALSE)
curve(0.5 * dnorm(x) + 0.5 * dnorm(x, sd = 2), add = TRUE, col = "red")

mix = disc(pt = c(1, 2), pr = c(0.5, 0.5))
test2 <- rmsnorm(n, varmix = mix, xi = xi)
hist(test2, breaks = 200, freq = FALSE)
curve(dmsnorm(x, varmix = mix, xi = xi), add = TRUE, col = "red")
curve(dmsnorm2(x, varmix = mix, xi = xi), add = TRUE, col = "red")


m1 = c(sqrt(2 * mix$pt^2/pi) %*% mix$pr)
m2 = sum(mix$pr * mix$pt^2)
mu = m1 * (xi - 1/xi)
sigma = sqrt((m2 - m1^2) * (xi^2 + 1/xi^2) + 2 * m1^2 - m2)

test3 = (test2 - mu)/sigma
mean(test3)
sd(test3)
hist(test3, breaks = 200, freq = FALSE)

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
pmnorm <- function(q, varmix = disc(1, 1)){
    z = outer(q, sqrt(2 * varmix$pt^2), FUN = "/")
    d = 0.5 * (1 + erf(z))
    d %*% varmix$pr
}

curve(pmnorm(x, varmix = disc(1, 1)), -5, 5)
curve(pnorm(x), add = TRUE, col = "red", lty = 3)
curve(pmnorm(x, varmix = disc(pt = c(0.5, 2), pr = c(0.5, 0.5))), add = TRUE, col = "blue")

pmsnorm <- function (q, varmix = disc(1, 1), xi = 1){
    m1 = 2/sqrt(2 * pi)
    mu = m1 * (xi - 1/xi)
    sigma = sqrt((1 - m1^2) * (xi^2 + 1/xi^2) + 2 * m1^2 - 1)
    z = q * sigma + mu
    Xi = xi^sign(z)
    g = 2/(xi + 1/xi)
    Probability = Heaviside(z) - sign(z) * g * Xi * pmnorm(q = -abs(z)/Xi, varmix = varmix)
    Probability
}


pmsnorm2 <- function (q, varmix = disc(1, 1), xi = 1){
    ## m1 = 2/sqrt(2 * pi)
    ## mu = m1 * (xi - 1/xi)
    ## sigma = sqrt((1 - m1^2) * (xi^2 + 1/xi^2) + 2 * m1^2 - 1)
    ## z = q * sigma + mu
    Xi = xi^sign(q)
    g = 2/(xi + 1/xi)
    Probability = Heaviside(q) - sign(q) * g * Xi * pmnorm(q = -abs(q)/Xi, varmix = varmix)
    Probability
}

curve(pmnorm(x, varmix = disc(1, 1)), -5, 5)
curve(pmsnorm(x, varmix = disc(1, 1)), add = TRUE, col = "red", lty = 3)
curve(pmsnorm(x, varmix = disc(1, 1), xi = 2), add = TRUE)
curve(psnorm(x, xi = 2), add = TRUE, col = "green", lty = 3)
curve(pmsnorm2(x, varmix = disc(1, 1), xi = 2), add = TRUE, n = 10001)
curve(pmsnorm2(x, varmix = disc(1, 1), xi = 2) - .2, add = TRUE, n = 10001)
curve(dmsnorm(x, xi = 2), -5, 5)
abline(v = 0, col = "red")
integrate(function(x) dmsnorm(x, xi = 2), lower = -Inf, upper = -.5)
pmsnorm2(-0.5, xi = 2)
abline(h = 0.2, col = "red")

ierf <- function(x) {
    qnorm((1 + x) /2) / sqrt(2)
}

qmnorm <- function(p, varmix = disc(1, 1), xi = 1){
    outer(ierf(2 * p - 1), sqrt(2 * varmix$pt^2), FUN = "*") %*% varmix$pr
}



curve(qnorm(x), 0, 1)
curve(qmnorm(x), add = TRUE, col = "red", lty = 3)
curve(qmnorm(x, varmix = disc(pt = c(0.5, 2), pr = c(0.5, 0.5))), add = TRUE, col = "red", lty = 3)
curve(qmnorm(x, varmix = disc(pt = c(0.5, 2), pr = c(0.5, 0.5)), xi =2), add = TRUE, col = "red", lty = 3)

qnorm2 <- function(p, mean = 0, sd = 1){
    prob = double(length(p))
    for(i in 1:length(p)){
        prob[i] = uniroot(function(x) pnorm(x, mean = mean, sd = sd) - p[i], interval = c(-10, 10))$root
    }
    prob
}

curve(qnorm2(x), 1e-5, 1 - 1e-5)
curve(qnorm(x), add = TRUE, col = "red")


library(raster)
waldo = brick(system.file("external/DepartmentStore.grd", package = "raster"))




curve(dlogis(x), -5, 5)
curve(dnorm(x), add = TRUE, col = "red")



########################################################################
## Looking at NA's
########################################################################

# store the one-but-least C-integer. The L in the end forces the number
# to be "integer", not "numeric"
x <- -2147483647L
typeof(x)

# adding an integer works fine, since we move further into the range:
typeof(x+1L)

# substracting an integer gives a warning telling us that the result is out-of-range:
typeof(x-1L)


# substracting a non-integer 1 ("numeric") yields a non-integer:
typeof(x-1)

typeof(x+1L)

########################################################################
## ggmap and the wine data
########################################################################

require(ggmap)
require(mapproj)

## source data from: http://www.napawineproject.com/project-notes/index.asp
## Excel file exported as CSV with header rows removed

wine <- read.csv("Napa-Winery-List62012.csv", stringsAsFactors=FALSE,
                 header = TRUE)
nwines <- nrow(wine)

# wineries that accept tastings, either by appointment or walk-ins
tasting <- wine$App=="No" | wine$App=="Yes"
# sum(tasting)
# wine$Address[tasting]

# wineries with valid addresses (P.O. Boxes not counted as valid)
good.address <- wine$Address != ""
good.address[grep("Box", wine$Address)] <- F
# wine$Address[!good.address]

# wineries you can visit to taste
visit <- good.address & tasting
# sum(visit) # 354 total

## filter data frame to visitable wineries
visit.wine <- wine[visit,]

## paste address fields into a single CA address for geocoding
addresses <- with(visit.wine,
 paste(Address, City, "CA", sep=", ")
)
## some non-ASCII characters in source data will mess up geocoding
addresses <- iconv(addresses, to="ASCII",sub="")

## convert addresses to lat/long
# locs <- geocode(addresses)
## I got a "there is no connection 3" error when trying to do all at
## once, so breaking up into batches

loc1 <- geocode(addresses[1:100])
loc2 <- geocode(addresses[101:200])
loc3 <- geocode(addresses[201:300])
loc4 <- geocode(addresses[301:length(addresses)])
locs <- rbind(loc1, loc2, loc3, loc4)

visit.wine$lon <- locs$lon
visit.wine$lat <- locs$lat

## Couldn't get the watercolor maps to work
SHmap <- qmap(c(lon=map.center$lon, lat=map.center$lat), source="stamen",
              maptype="watercolor")

map.center <- geocode("Saint Helena, CA")
SHmap <- qmap(c(lon=map.center$lon, lat=map.center$lat), source="stamen",
              zoom=8)
SHmap + geom_point(
  aes(x=lon, y=lat, colour=App), data=visit.wine[1:10, ]) +
  scale_colour_manual(values=c("dark blue","orange"))+
  labs(colour="Appointment Required")

library(fGarch)
test.df <- data.frame(lon = map.center$lon + rged(1000, mean = 0, sd = 0.1,
                        nu = 1),
                      lat = map.center$lat + rged(1000, mean = 0, sd = 0.1,
                        nu = 1))
test.df$dist = sqrt((test.df$lon - map.center$lon)^2 +
  (test.df$lat - map.center$lat)^2)

test <- qmap(c(lon=map.center$lon, lat = map.center$lat),
             source = "stamen", zoom = 1)

test + geom_point(aes(x = lon, y = lat, colour = dist), data = test.df)
## test <- qmap(c(lon=map.center$lon, lat = map.center$lat),
##              source = "stamen", zoom = 7)

pdf(file = "randPointsMap.pdf")
print(test + geom_point(aes(x = lon, y = lat, colour = dist), data = test.df))
graphics.off()
system("evince randPointsMap.pdf&")




########################################################################
## Rgaphviz
########################################################################

library(rgraphviz)

set.seed(123)
V <- letters[1:10]
M <- 1:4
g1 <- randomGraph(V, M, 0.2)
plot(g1)
plot(g1, "neato")

g2 <- addEdge("a", "g", g1, 1)
plot(g2)


construct <- new("graphNEL", nodes = c("Construction", "Growth",
                               "Least Squares", "Geometric",
                               "Absolute Changes", "Index"),
                 edgemode = "directed")

## Construct second level
construct2 <- addEdge("Construction", "Growth", construct, 1)
construct2 <- addEdge("Construction", "Index", construct2, 20)

## Third level
construct3 <- addEdge("Growth", "Geometric", construct2, 1)
construct3 <- addEdge("Growth", "Least Squares", construct3, 1)
construct3 <- addEdge("Growth", "Absolute Changes", construct3, 20)
plot(construct3)

defAttrs <- getDefaultAttrs()



########################################################################
## Improve y-axis problem when a the scale are very different
########################################################################

test.df <- data.frame(type = rep(c("Im", "Ex"), each = 60),
                      region = rep(rep(month.abb[1:6], each = 10), 2),
                      category = rep(letters[1:10], 12),
                      trade = rnorm(120), stringsAsFactors = FALSE)


ggplot(test.df, aes(x = region, y = trade)) +
  geom_bar(aes(fill = category)) + facet_wrap(~type, nc = 1, scale = "free_y")



########################################################################
## Title: Play around with Lucia's data
## Time : 24-07-2012
########################################################################

## Install the packages online
install.packages("foreign")
install.packages("MASS")
install.packages("GB2")


## Load this library to read file from other statistical package
library(foreign)
## This library is used for fitdistr() function for fitting distribution
library(MASS)
## This library allows one to fit GB2 distribution
library(GB2)

## Read the data
lucia.df <- read.dta("lucia.dta")

## Assign the data matrix to a single vector, since this is the only
## variable of interest
myVar <- lucia.df$dec_ae

## Histogram of the variable
hist(myVar, breaks = 100, freq = FALSE)

## Fit the log-normal distribution then super-impose on the data
ln.fit <- fitdistr(myVar, "lognormal")
curve(dlnorm(x, meanlog = ln.fit$estimate[1], sdlog = ln.fit$estimate[2]),
      add = TRUE, col = "blue")

## Fit the log-normal distribution with weights
wllLN <- function(beta, data, weights){
  ll = -sum(weights * dlnorm(data, meanlog = beta[1],
    sdlog = beta[2], log = TRUE))
  ll
}

wln.fit <- nlminb(start = c(8, 0.5), wllLN, weights = lucia.df$iwght,
                  data = myVar)
curve(dlnorm(x, meanlog = wln.fit$par[1], sdlog = wln.fit$par[2]),
      add = TRUE, col = "blue")


## By the look of the plot, the gb2 doesn't seem to fit very well.
bg2.fit <- mlfit.gb2(z = myVar, w = lucia.df$iwght)
curve(dgb2(x, shape1 = bg2.fit[[2]]$par[1], scale = bg2.fit[[2]]$par[2],
           shape2 = bg2.fit[[2]]$par[3], shape3 = bg2.fit[[2]]$par[4]),
      add = TRUE, col = "red")



## We will take a closer look at the log-normal fit
## qqplot
qqplot(sort(myVar), qlnorm(1:length(myVar)/length(myVar),
                           meanlog = ln.fit$estimate[1],
                           sdlog = ln.fit$estimate[2]),
       xlab = "Sample Quantile",
       ylab = "Theoretical Quantile")
abline(a = 0, b = 1, col = "red")
## Comments (Michael): The qqplot appears to be fine, the tail of the
##                     data however does look to be heavier than what
##                     the log-normal would predict


## qqplot of the gb2 distribution
qqplot(sort(myVar), qgb2(1:length(myVar)/length(myVar),
                         shape1 = bg2.fit[[2]]$par[1],
                         scale = bg2.fit[[2]]$par[2],
                         shape2 = bg2.fit[[2]]$par[3],
                         shape3 = bg2.fit[[2]]$par[4]),
       xlab = "Sample Quantile",
       ylab = "Theoretical Quantile")
abline(a = 0, b = 1, col = "red")


## Plot the empirical cumulative distribution function and theoretical
## distribution function
plot(ecdf(myVar))
curve(plnorm(x, meanlog = ln.fit$estimate[1], sdlog = ln.fit$estimate[2]),
      add = TRUE, col = "red")
## Comments (Michael): The comparison between the empirical CDF and the
##                     fitted CDF seems to be good, however from the
##                     above plot we know the problem lies in both tails
##                     of the data.


## Kolmogorov Smirnov test
ks.test(unique(myVar), function(x) plnorm(x, meanlog = ln.fit$estimate[1],
                                          sdlog = ln.fit$estimate[2]))
## Comments (Michael): The Kolmogorov Smirnov test is rejected, however
##                     depending on the purpose of your
##                     analysis/thesis. The test may reject the
##                     log-normal distribution but from the plot it
##                     still seems a good distritbuion.


########################################################################
## Title : Play around how to aggregate countries
## Time  : 22-07-2012
########################################################################


## NOTE (Michael): The aggregation should not be based on geographical
##                 location, rather it is the representation of the
##                 sovereign state. Should represent at the time of the
##                 state. (e.g. China should not have Hong Kong prior to
##                 1998 and Macao to 1999)
##
## NOTE (Michael): No need to aggregate Serbia and Montenegro, since
##                 they are recognised as independent sovereign state.
##
## NOTE (Michael): Kosovo is not yet recognised by the UN, with divided
##                 decision.
##
## NOTE (Michael): Taiwan and Kosovo with limited recognition and
##                 disbuted territory will be discarded.
##
## NOTE (Michael): Probably dont need a mean, since only sum and
##                 weighted mean are required for aggregating disputed
##                 area.
##
## NOTE (Michael): The function can not be vectorized, because the
##                 country case differ with each variable.
##
## NOTE (Michael): Can not discard the sub-region after the computation
##                 immediately, since it might be used later for
##                 constructing other variables. Should do it at the end
##                 of member aggregation.
##
## NOTE (Michael): Need to aggregat the sub-regions but at the same time
##                 keep the sub-regions for weights. Then the
##                 construction and the second level of aggregation is
##                 computed with mutually exclusive set. (Update: The
##                 final data would not be mutually exclusive due to
##                 regional aggregation having FAOST_CODE)

## Just for China, Hong Knog, Macau and Taiwan

cty.df <- read.csv(file = "Countries.csv", header = TRUE, na.string = "",
                   stringsAsFactors = FALSE)

test <- getWDI("AG.LND.AGRI.ZS")

## Use the rural density (%) and total population for example
test.df <- getWDItoSYB(c("SP.URB.TOTL.IN.ZS", "SP.POP.TOTL"),
                       c("SP.URB.TOTL.IN.ZS", "SP.POP.TOTL"), date = 2010)
work.df <- test.df$data[-258, ]
work.df <- rbind(work.df, c(351, NA, NA, NA))
work.df <- rbind(work.df, c(357, NA, NA, NA))


## The aggregation function, data is the data.frame, variable is the
## variable of interest, weight is the weighting variable.
aggCHMT <- function(data, variable, weight = NULL,
                    aggMethod = c("sum", "weightMean")){
  aggMethod = match.arg(aggMethod)
  if(aggMethod == "sum"){
    if(is.na(data[data$FAOST_CODE == 351, variable])){
      if(!is.na(data[data$FAOST_CODE == 357, variable])){
        agg <- sum(data[data$FAOST_CODE %in% c(357, 96, 128),
                        variable], na.rm = TRUE)
      } else if(!is.na(data[data$FAOST_CODE == 41, variable])){
        agg <- sum(data[data$FAOST_CODE %in% c(41, 96, 128, 214),
                        variable], na.rm = TRUE)
      } else {
        warning("not enough data to compute China, NA is assigned")
        agg <- NA
      }
      data[data$FAOST_CODE == 351, variable] <- agg
    }
  }
  if(aggMethod == "weightMean"){
    if(is.null(weight)) stop("Weighting variable not specified")
    if(is.na(data[data$FAOST_CODE == 351, variable])){
      if(!is.na(data[data$FAOST_CODE == 357, variable])){
        unweight <- sum(data[data$FAOST_CODE %in% c(357, 96, 128),
                             variable] *
                        data[data$FAOST_CODE %in% c(357, 96, 128),
                             weight], na.rm = TRUE)
        CHMTweight <- sum(data[data$FAOST_CODE %in% c(357, 96, 128),
                               weight], na.rm = TRUE)
        agg = unweight/CHMTweight
      } else if(!is.na(data[data$FAOST_CODE == 41, variable])){
        unweight <- sum(data[data$FAOST_CODE %in% c(41, 96, 128, 214),
                             variable] *
                        data[data$FAOST_CODE %in% c(41, 96, 128, 214),
                             weight], na.rm = TRUE)
        CHMTweight <- sum(data[data$FAOST_CODE %in% c(41, 96, 128, 214),
                               weight], na.rm = TRUE)
        agg = unweight/CHMTweight
        ## data[data$FAOST_CODE %in% c(41, 96, 128, 214), variable] <- NA
      } else {
        warning("not enough data to compute China, NA is assigned")
        agg <- NA
      }
      data[data$FAOST_CODE == 351, variable] <- agg
    }
  }
  data
}



## Construct the population first
work2.df <- aggCHMT(work.df, "SP.POP.TOTL")
## Then construct the average of the rural density
work3.df <- aggCHMT(work2.df, "SP.URB.TOTL.IN.ZS", weight = "SP.POP.TOTL",
                    agg = "weightMean")


## Aggregating Kosovo (275) and Serbia (272)
aggKS <- function(data, variable, weight, agg = c("sum", "weightMean")){
  agg = match.arg(agg)
  if(!is.na(data[data$FAOST_CODE == 275, variable])){
    if(agg == "sum"){
      tmp <- sum(data.frame(data[data$FAOST_CODE %in%
                                 c(272, 275), variable, drop = FALSE]))
    } else if(agg == "weightMean"){
      if(missing(weight) || !is.character(weight) || !is.numeric(weight))
        stop("Weighting variable not specified")
      wt <- ifelse(is.character(weight),
                   data[data$FAOST_CODE %in% c(272, 275), weight],
                   weight)
      unweight <- sum(data.frame(data[data$FAOST_CODE %in% c(272, 275),
                                      variable] * wt))
      KSweight <- sum(data[data$FAOST_CODE %in% c(272, 275), weight,
                           drop = FALSE])
      tmp <- unweight/KSweight
    }
  }
  data[data$FAOST_CODE == 272, variable] <- tmp
  data
}


work4.df <- aggKS(work3.df, "SP.URB.TOTL.IN.ZS", weight = "SP.POP.TOTL",
                  agg = "weightMean")


## Aggregating Tansania (215) and Zanzibar (3698)
aggTS <- function(data, variable, weight = NULL,
                  agg = c("sum", "weightMean")){
  agg = match.arg(agg)
  if(!is.na(data[data$FAOST_CODE == 3698, variable])){
    if(agg == "sum"){
      tmp <- colSums(data.frame(data[data$FAOST_CODE %in% c(215, 3698),
                                     variable]))
      data[data$FAOST_CODE %in% c(215, 3698), variable] <- NA
    } else if(agg == "weightMean"){
      if(is.null(weight)) stop("Weighting variable not specified")
      unweight <- colSums(data.frame(data[data$FAOST_CODE %in% c(215, 3698),
                                          variable] *
                                     data[data$FAOST_CODE %in% c(215, 3698),
                                          weight]))
      TSweight <- colSums(data[data$FAOST_CODE %in% c(215, 3698), weight])
      tmp <- unweight/TSweight
      data[data$FAOST_CODE %in% c(215, 3698), variable] <- NA
    }
    data[data$FAOST_CODE == 215, variable] <- tmp
  }
  data
}


## Aggregating West bank (245) Gaza strip (76), only do it if Palestine
## (299) is absent
aggGW <- function(data, variable, weight = NULL,
                  agg = c("sum", "weightMean")){
  if(is.na(data[data$FAOST_CODE == 299, variable])){
    agg = match.arg(agg)
    if(!is.na(data[data$FAOST_CODE == 76, variable])){
      if(agg == "sum"){
        tmp <- colSums(data.frame(data[data$FAOST_CODE %in% c(245, 76),
                                       variable]))
        data[data$FAOST_CODE %in% c(245, 76), variable] <- NA
      } else if(agg == "weightMean"){
        if(is.null(weight)) stop("Weighting variable not specified")
        unweight <- colSums(data.frame(data[data$FAOST_CODE %in% c(245, 76),
                                            variable] *
                                       data[data$FAOST_CODE %in% c(245, 76),
                                            weight]))
        TSweight <- colSums(data[data$FAOST_CODE %in% c(245, 76), weight])
        tmp <- unweight/TSweight
        data[data$FAOST_CODE %in% c(245, 76), weight] <- NA
      }
      data[data$FAOST_CODE == 299, variable] <- tmp
    }
  }
  data
}

testFunc <- function(){
  data(countryProfile, envir = environment())
  head(countryProfile)
}

## This is the first level aggregation, should probably do the whole thing.
## (1) Aggregate sub-region
## (2) Assign NA to aggregated region
## (3) Discard non-member countries.
aggCountry <- function(data, variable, weight = NULL,
                      agg = c("sum", "weightMean")){
  agg = match.arg(agg)
  ## Aggregate all the sub regions first
  n.var = length(variable)
  if(!all(n.var == c(length(weight), length(agg))))
    stop("Length of inputs are not the same")
  ## Probably should avoid the use of data
  aggC <- data
  for(i in 1:n.var){
    aggC <- aggCHMT(aggC, variable = variable[i], weight = weight[i],
                    agg = agg[i])
    aggC <- aggKS(aggC, variable = variable[i], weight = weight[i],
                  agg = agg[i])
    aggC <- aggTS(aggC, variable = variable[i], weight = weight[i],
                  agg = agg[i])
    aggC <- aggGW(aggC, variable = variable[i], weight = weight[i],
                  agg = agg[i])
  }
  ## Think of a better way to do this, might be able to use comment or
  ## attributes
  aggC[c(41, 96, 128, 214, 275, 299, 357, 3968), variable] <- NA
  data(countryProfile, envir = environment())
  countryStatus <- subset(countryProfile, select = c("FAOST_CODE", "STATUS"))
  final.df <- merge(countryStatus, aggC, all.x = TRUE, by = "FAOST_CODE")
  member <- subset(final.df, subset = STATUS == "Member State")
  nonmember <- subset(final.df, subset = STATUS != "Member State")
  allNA <- apply(nonMember[, variable], 2, function(x) all(is.na(x)))
  if(!all(allNA))
    warning(paste(variable[allNA], " is not fully constructed.\n"))
  member
}


## This is just a preliminary function
aggRegion <- function(data, inFaoCode, outFaoCode,  var, weightVar = NULL, year,
                      aggMethod = c("sum", "mean", "weighted.mean")){
  if(!outFaoCode %in% unique(data$FAOST_CODE)){
    nr = NROW(data)
    data[nr + 1, "FAOST_CODE"] <- outFaoCode
    data[nr + 1, "Year"] <- year
  }
  data[data$FAOST_CODE == outFaoCode & data$Year == year, var] <-
    apply(data[data$FAOST_CODE %in% inFaoCode & data$Year == year,
               var, drop = FALSE], 2, aggMethod, w = data$weightVar)
  data
}

aggRegion(test.df, inFaoCode = 1:5, outFaoCode = 5000,
          var = "SP.POP.TOTL", year = 1995, aggMethod = "sum")

region.df <- data.frame(reg = rep(c("Eastern Africa", "Middle Africa",
                          "Northern Africa", "Southern Africa",
                          "Western Africa" ,"Caribbean" ,"Central America",
                          "South America", "Northern America", "Central Asia",
                          "Eastern Asia", "Southern Asia",
                          "South-Eastern Asia", "Western Asia",
                          "Eastern Europe", "Northern Europe",
                          "Southern Europe", "Western Europe",
                          "Australia and New Zealand",
                          "Melanesia Micronesia", "Polynesia"), 3),
                        subReg = 1:(21 * 3),
                        runif(63, 0, 1))



aggRegion <- function(data, variable, weight = NULL, regionCode,
                  agg = c("sum", "mean", "weightMean")){

}




########################################################################
## Title: Test the use of units
## Date: 02-08-2012
########################################################################

library(Hmisc)

test.df <- data.frame(x = ts(1:12, start = c(2000, 1), frequency = 12),
                      y = ts(1:12, start = c(2001, 1), frequency = 12))
units(test.df$x) = "cm"
units(test.df$y) = "m"
comment(test.df) <- "this is a test data set"

describe(test.df)
contents(test.df)

library(IRanges)
test2.df <- DataFrame(x = 1:10, y = letters[1:10])
metadata(test2.df) <- list(units=list(x = "cm", y = "m"))

str(subset(test2.df, select = x))

########################################################################
## Title: Cache results using environment to speed up computation
## Date: 02-08-2012
########################################################################

mkFact <- function (){
  .env <- environment()
  function(n) {
    nm <- paste(".", n, sep = "/")
    if (n <= 0)
      return(1)
    if (exists(nm, envir = .env))
      return(get(nm, envir = .env))
    ans <- n * Recall(n - 1)
    assign(nm, ans, envir = .env)
    ans
  }
}
Fact <- mkFact()

Fact(4)
objects(env = environment(Fact), all = TRUE)
Fact(6)
objects(env = environment(Fact), all = TRUE)



########################################################################
## Title: Regression with triangular distribution
## Date: 05-08-2012
########################################################################

dtriangle <- function(x, l = -1, m = 0, u = 1){
  h = 2/(u - l)
  for(i in length(x)){
    if(x < l | x > u){
      d = 0
    } else if(x <= m){
      d = h/((m - l) * (x - l))
    } else {
      d = -h/((u - m) * (x - m))
    }
  }
  d
}


########################################################################
## Title: GoogleVis test
## Date: 16-08-2012
########################################################################

## Plot world wide earth quakes of the last 30 days with magnitude >= 4.0
library(XML)
library(googleVis)

## Get earthquake data of the last 30 days
eq = readHTMLTable(readLines("http://www.iris.edu/seismon/last30.html"),
  colClasses=c("factor", rep("numeric", 4), "factor"), which=2)

## Create a lat:long location variable
eq$loc = paste(eq$LAT, eq$LON, sep=":")

## Create a geo chart
G <- gvisGeoChart(eq, "loc", "DEPTH km", "MAG",
                   options = list(displayMode = "Markers",
                   colorAxis = "{colors:['purple', 'red', 'orange', 'grey']}",
                   backgroundColor = "lightblue"), chartid = "test")

## Display geo chart in local web browser
fn <- plot(G)

## The location of the html
fn

## Subset the map to only Iran
IR <- gvisGeoChart(eq, "loc", "DEPTH km", "MAG",
                  options=list(displayMode="Markers", region="IR",
                  colorAxis="{colors:['purple', 'red', 'orange', 'grey']}",
                  backgroundColor="lightblue"), chartid="Iran")
plot(IR)



########################################################################
## Title: Test how tcl/tk works in producing a GUI interface
## Date: 16-08-2012
## NOTE: The tcltk somehow ends the R session, not sure why.
########################################################################

require(tcltk)
require(cairoDevice)
mydialog <- function(){

  xvar <- tclVar("")
  yvar <- tclVar("")
  zvar <- tclVar("")

  tt <- tktoplevel()
  tkwm.title(tt,"MYTEST")
  x.entry <- tkentry(tt, textvariable=xvar)
  y.entry <- tkentry(tt, textvariable=yvar)
  z.entry <- tkentry(tt, textvariable=zvar)

  reset <- function() {
    tclvalue(xvar)<-""
    tclvalue(yvar)<-""
    tclvalue(zvar)<-""
  }

  reset.but <- tkbutton(tt, text="Reset", command=reset)

  submit <- function() {
    x <- as.numeric(tclvalue(xvar))
    y <- as.numeric(tclvalue(yvar))
    z <- as.numeric(tclvalue(zvar))
    tkmessageBox(message=paste("x + y + z = ", x+y+z, ""))
  }
  submit.but <- tkbutton(tt, text="submit", command=submit)

  quit.but <- tkbutton(tt, text = "Close Session",
                       command = function() {
                         q(save = "no")
                         tkdestroy(tt)
                       }
                       )

  tkgrid(tklabel(tt,text="Put your variables.."),columnspan=3, pady = 10)
  tkgrid(tklabel(tt,text="x variable"), x.entry, pady= 10, padx= 10)
  tkgrid(tklabel(tt,text="y variable"), y.entry, pady= 10, padx= 10)
  tkgrid(tklabel(tt,text="z variable"), z.entry, pady= 10, padx= 10)
  tkgrid(submit.but, reset.but, quit.but, pady= 10, padx= 10)

}

test = mydialog()


## A different attempt

## Get the function from here
## source("https://raw.github.com/gist/3342915/89d8942a9ae3a82138507f359dd5fcd13e04b939/varEntryDialog.r")

## This doesn't crash, so can compare and investigate with the previous function
vals <- varEntryDialog(vars=c('Variable1', 'Variable2'))


########################################################################
## Title: Parallel process
## Date: 20-08-2012
########################################################################

library(parallel)

detectCores()
nonpar.test <- function(text.var, gc.rate=10){
    ntv <- length(text.var)
    require(parallel)
    pos <-  function(i) {
        paste(sapply(strsplit(tolower(i), " "), nchar), collapse=" | ")
    }
    lapply(seq_len(ntv), function(i) {
            x <- pos(text.var[i])
            if (i%%gc.rate==0) gc()
            return(x)
        }
    )
}

nonpar.test(rep("I wish I ran in parallel.", 20))

par.test <- function(text.var, gc.rate=10){
    ntv <- length(text.var)
    require(parallel)
    pos <-  function(i) {
        paste(sapply(strsplit(tolower(i), " "), nchar), collapse=" | ")
    }
#======================================
    cl <- makeCluster(mc <- getOption("cl.cores", 4))
    clusterExport(cl=cl, varlist=c("text.var", "ntv", "gc.rate", "pos"),
        envir=environment())
    parLapply(cl, seq_len(ntv), function(i) {
#======================================
            x <- pos(text.var[i])
            if (i%%gc.rate==0) gc()
            return(x)
        }
    )
}

par.test(rep("I wish I ran in parallel.", 20))

system.time(pos(rajSPLIT$dialogue, parallel=T))

system.time(pos(rajSPLIT$dialogue, progress.bar =F))



########################################################################
## Title: Improve SYB plots
## Date: 24-08-2012
########################################################################

test.df = data.frame(a = c(runif(40), rnorm(40), rexp(40),
                       rchisq(40, df = 1), rbinom(40, 1, 0.5)),
  b = rep(c("Uniform", "Normal", "Exponential", "Chi-squared", "Binomial"),
    each = 40))

## violin plot
ggplot(data = test.df, aes(x = b, y = a)) + geom_violin() + geom_jitter() +
  xlab(NULL) + ylab(NULL)

## Box plot
ggplot(data = test.df, aes(x = b, y = a)) + geom_boxplot() + geom_jitter() +
  xlab(NULL) + ylab(NULL)


test2.df = data.frame(a = rnorm(3000), b = rnorm(3000))

ggplot(data = test2.df, aes(x = a, y = b)) + geom_point(alpha = 1/3) +
  xlab(NULL) + ylab(NULL)


## Mosaic plots/tree map
library(portfolio)
mosaic.df = read.csv("http://datasets.flowingdata.com/post-data.txt")

map.market(id=mosaic.df$id, area=mosaic.df$views, group=mosaic.df$category,
           color=mosaic.df$comments, main="Test Mosaic")




library(FAOSYB)
rur = getWDI(indicator = "SP.RUR.TOTL")
urb = getWDI(indicator = "SP.URB.TOTL")

rur.df = na.omit(subset(rur, Country %in% c("High income", "Middle income",
  "Low income")))[, -2]


urb.df = na.omit(subset(urb, Country %in% c("High income", "Middle income",
  "Low income")))[, -2]

all.df = merge(rur.df, urb.df)
mall.df = melt(all.df, c("Country", "Year"))
mall.df$variable = as.character(mall.df$variable)
mall.df$Country = as.character(mall.df$Country)

## Doesn't seem to work very well since the scales are way too different
map.market(id = mall.df$Country, group = mall.df$Year,
           area = mall.df$value, color = mall.df$value)


########################################################################
## Title: Spatial data analysis
## Date: 25-08-2012
########################################################################

library(sp)
getClass("Spatial")
getClass("CRS")

m = matrix(c(0, 0, 1, 1), nc = 2, dimnames = list(NULL, c("min", "max")))
crs = CRS(projargs = as.character(NA))

S = Spatial(bbox = m, proj4string = crs)


## Read the location of the CRAN depositories
CRAN.df = read.table("http://www.asdar-book.org/datasets/CRAN051001a.txt",
  header = TRUE)
CRAN.mat = cbind(CRAN_df$long, CRAN_df$lat)
row.names(CRAN.mat) = 1:NROW(CRAN.mat)

## Set the projectiong
allCRS = CRS("+proj=longlat +ellps=WGS84")
CRAN.sp = SpatialPoints(CRAN.mat, proj4string = llCRS)

## Plot the data (no boundaries)
plot(CRAN.sp)


## Accessor function
bbox(CRAN.sp)
proj4string(CRAN.sp)
coordinates(CRAN.sp)

########################################################################
## Title: Scrapping member countries
## Date: 29-08-2012
########################################################################

library(XML)

## Need to find settings to read the table properly
url = "http://www.fao.org/legal/home/fao-members/en/"
test = readHTMLTable(url, stringsAsFactors = FALSE)

members = test[[1]]$V1



########################################################################
## Title: Test the efficiency of matrix and data.frame efficiency
## Date: 31-08-2012
########################################################################

test.df = data.frame(matrix(rnorm(1000000), nc = 100))

f1 = function(df){
  df * df
}

f2 = function(df){
  tmp = matrix(unlist(test.df), nc = NCOL(df))
  tmp2 = tmp * tmp
  data.frame(tmp2)
}

## It seems converting data.frame to matrix takes a long time already so
## better to just use multiplicationusing data frame

########################################################################
## Title: Test the transparence and overlap of ggplot
## Date: 31-08-2012
########################################################################

library(plyr)
library(reshape)

n = 500
set.seed(587)
x <- rnorm(n, 0.8, 1.2)
e <- rnorm(n, 0, 2) * (abs(x)^1.5 + 0.5)
y <- 8 * x - x^3 + e
plot(x, y)
test = loess(y ~ x)
testFit = arrange(data.frame(x = test$x, f = test$fitted), x)
lines(testFit)

test.df = data.frame(x =x , y = y)
plot(x, y, pch = 19, cex = 0.6)
bn = 100
sn = 1000
base = arrange(subset(test.df, select = x, drop = FALSE), x)
for(i in 1:bn){
  tmp = arrange(test.df[sample(1:nrow(test.df), sn,
           replace = TRUE), ], x)
  smth = with(tmp, loess(y ~ x))
  lines(smth$x, smth$fitted, col = rgb(0, 0, 0.5, alpha = 0.3),
        lwd = 3)
  lines(smth$x, smth$fitted, col = rgb(0, 0, 0.8, alpha = 0.3),
        lwd = 2)
  lines(smth$x, smth$fitted, col = rgb(1, 1, 1, alpha = 0.3),
        lwd = 1)
}

## use the bivariate density approach
library(ks)
data(faithful)
H <- Hpi(x=faithful)
fhat <- kde(x=faithful, H=H)
plot(fhat, display="filled.contour2")
points(faithful, cex=0.5, pch=16)

test = expand.grid(x= seq(0, 1, length = 100), y = seq(0, 1, length = 100))
test$dist = sqrt((test$x + 0.5)^2 + (test$y - 0.5)^6)
## test$dist = 1/NROW(test)
test$dist = test$dist/sum(test$dist)



test.mat = test[sample(1:NROW(test), size = 100, prob = test$dist,
  replace = TRUE), c("x", "y")]
H = Hpi(x = test.mat)
fhat = kde(x = test.mat, H = H)
plot(fhat, display = "image", col = rev(heat.colors(100)))
plot(fhat, display="filled.contour2", cont=seq(1, 100,by=1))
points(test.mat, pch = 19)



########################################################################
## Title: Test the speed of the ddply using parallel + idata.frame
## Date: 05-09-2012
########################################################################

library(plyr)
library(reshape)

## test idata.frame
test.df = data.frame(a = rep(letters[1:10], 100000), v1 = rnorm(1000000),
  v2 = rnorm(1000000), v3 = rnorm(1000000))

test.df = data.frame(a = rep(letters[1:10], 100000),
  matrix(rnorm(1000000 * 7), nc = 7))
colnames(test.df) = c("name", paste("v", 1:7, sep = ""))

itest.df = idata.frame(test.df)

system.time(
  agg1.df <- ddply(test.df, .variables = .(name), .fun = function(x) sum(x$v1))
  )

## idata.frame seems to increase the speed by a margin
system.time(
  iagg1.df <- ddply(itest.df, .variables = .(name), .fun = function(x) sum(x$v1))
  )

## Parallel does not increase the efficiency when the execution time is small
system.time(
  pagg1.df <- ddply(test.df, .variables = .(name), .fun = function(x) sum(x$v1),
                    .parallel = TRUE)
  )


system.time(
  piagg1.df <- ddply(test.df, .variables = .(name), .fun = function(x) sum(x$v1),
                     .parallel = TRUE)
  )


test4.df = ddply(test.df, .(name),
  .fun = function(x) colSums(x[, c("v5", "v6"), drop = FALSE] *
    x[, c("v7", "v7"), drop = FALSE])/colSums(x[, c("v7", "v7"), drop = FALSE]))

## immutable data frame doesn't seem to allow you to do subsetting etc.
test4.df = ddply(itest.df, .(name),
  .fun = function(x) colSums(x[, c("v5", "v6"), drop = FALSE] *
    x[, c("v7", "v7"), drop = FALSE])/colSums(x[, c("v7", "v7"), drop = FALSE]))


########################################################################
## Title: Work to produce fluorescent looking color
## Date: 06-09-2012
########################################################################

n = 1000
x = rnorm(n)
y = 2 * x + rnorm(n, sd = abs(x/2))

plot(x, y, cex = 0.5, pch = 19)
test.df = data.frame(x = x, y = y)
nb = 100
for(i in 1:nb){
  newSamp = test.df[sample(1:n, 20, replace = TRUE), ]
  bt.lm = lm(y ~ x, data = newSamp)
  ## abline(a = coef(bt.lm)[1], b = coef(bt.lm)[1],
  ##        col = rgb(0, 0, 0.5, alpha = 0.3), lwd = 3)
  abline(coef = coef(bt.lm),
         col = rgb(0, 0, 0.5, alpha = 0.1), lwd = 3)
}

plot(1:100, col = rgb(0, 0, 1:100/100, alpha = 0.8), cex = 3, pch = 19)


varLines = function(x, y, lwd){
  new.n = 300 * length(x)
  d = cbind(x, y)
  ds = spline(d, n = new.n)
  print(str(ds))
  myCol = rgb(0, 0, 1:new.n/new.n, alpha = 0.5)
  k = sqrt(lwd) + 1
  myLwd = c(abs(seq(-sqrt(k) + 1, sqrt(k) - 1, length = new.n)) + 1)^2
  for(i in 1:(new.n - 1)){
    lines(c(ds$x[i], ds$x[i + 1]), c(ds$y[i], ds$y[i + 1]), col = myCol[i],
          lwd = myLwd[i])
  }
}

x = c(2, 3, 4, 5, 4)
y = c(2, 3, 5, 3, 2)
plot(x, y, pch = 19)
varLines(x, y, lwd = 50)

## Should check out polygon for efficiency

########################################################################
## Title: Test data.table
## Date: 07-09-2012
########################################################################

library(data.table)

## Create data.tables
test.df = data.table(x = c("a", "a", "b", "b", "b"), y = rnorm(5))
CARS = data.table(cars)

## This shows us what data.table is available in the ls() and also their
## size.
tables()

## Now add key to the table
setkey(test.df, x)

## Now we can also subset the rows
test.df["c", ]

library(FAOSYB)

test.df = data.frame(FAOST_CODE = rep(c(1:200000), each = 50),
  Year = rep(1961:2010, 200000), value1 = rnorm(200000 * 50),
  value2 = rnorm(200000 * 50),  stringsAsFactors = FALSE)
test2.df = test.df
colnames(test2.df)[3:4] = c("value3", "value4")

test.dt = data.table(test.df)
test2.dt = data.table(test2.df)
setkey(test.dt, FAOST_CODE, Year)

## Subset by the key
test.dt[J(5, 2010), ]

## Fail
gc(reset = TRUE)
system.time(
  test.dt[, transform(test.dt, b = sum(value1)), by = Year]
)
gc()

## Success
gc(reset = TRUE)
a = merge(test.df,test2.df)
gc()

########################################################################
## Title: Test whether it is better to merge large data.frame or
##        multiple small data.frame.
## Date: 07-09-2012
########################################################################

dfNames = paste("df", 1:1000, sep = "")

for(i in 1:length(dfNames)){
  assign(dfNames[i], {data.frame(FAOST_CODE = rep(1:200, each = 5),
                                 Year = rep(1996:2010, 200),
                                 Value = rnorm(1000))
                    }
         )
}

p1.df = Reduce(function(x, y) merge(x, get(y), by = c("FAOST_CODE", "Year")),
  as.list(paste("df", 2:500, sep = "")), init = df1)
p2.df = Reduce(function(x, y) merge(x, get(y), by = c("FAOST_CODE", "Year")),
  as.list(paste("df", 502:1000, sep = "")), init = df501)

colnames(p1.df) = c("FAOST_CODE", "Year", paste("V", 1:500, sep = ""))
colnames(p2.df) = c("FAOST_CODE", "Year", paste("V", 501:1000, sep = ""))

## Investigate time
system.time({
  bigMerge = merge(p1.df, p2.df, by = c("FAOST_CODE", "Year"))
})

system.time({
  base = df1
  for(i in 1:length(dfNames)){
  base = merge(base, get(dfNames[i]), by = c("FAOST_CODE", "Year"))
}
})


########################################################################
## Title: Simulation of missing value
## Date: 10-09-2012
########################################################################

library(WDI)

## Start with Total population
## Download data
pop = WDI(indicator = "SP.POP.TOTL", start = 2010, end = 2010)
pop = na.omit(pop)

## Order data
pop = pop[order(pop$SP.POP.TOTL, decreasing = TRUE), ]

## Subset  country by 30,000 rule
pop = pop[pop$SP.POP.TOTL >= 30000, ]

## Remove regions
pop = pop[-grep("income|devel|OECD|[E|e]urope|[A|a]sia|[A|a]Africa|[P|p]acific|[N|n]meriaUN|World|Heav|Euro", pop$country), ]

## Extract top country assume they are never missing
top = 50
mainCountry = pop[1:top, ]
otherCountry = pop[(top + 1):NROW(pop), ]

## Start the simulation based on maximum 50% missing
n.bt = 5000
missingSum = double(n.bt)
missingProp = double(n.bt)
for(i in 1:n.bt){
  missingProp[i] = runif(1, min = 0, max = 0.5)
  missingSum[i] = sum(otherCountry[sample(1:NROW(otherCountry),
              NROW(otherCountry) * (1 - missingProp[i]), replace = FALSE),
              "SP.POP.TOTL"])
}

simTotal = missingSum + sum(mainCountry[, "SP.POP.TOTL"])
hist(simTotal, breaks = 200)
abline(v = sum(pop[, "SP.POP.TOTL"]), col = "red")
summary(simTotal)

## Start the simulation based on maximum 20% missing
n.bt = 5000
missingSum = double(n.bt)
missingProp = double(n.bt)
for(i in 1:n.bt){
  missingProp[i] = runif(1, min = 0, max = 0.2)
  missingSum[i] = sum(otherCountry[sample(1:NROW(otherCountry),
              NROW(otherCountry) * (1 - missingProp[i]), replace = FALSE),
              "SP.POP.TOTL"])
}

simTotal = missingSum + sum(mainCountry[, "SP.POP.TOTL"])
hist(simTotal, breaks = 200)
abline(v = sum(pop[, "SP.POP.TOTL"]), col = "red")
summary(simTotal)


## Now investigate urban population
## Download data
upop = WDI(indicator = "SP.RUR.TOTL.ZS", start = 2010, end = 2010)
upop = na.omit(upop)

## Merge with population
apop = merge(upop, pop)

## Order data
apop = apop[order(apop$SP.POP.TOTL, decreasing = TRUE),]

## Compute true value
trueMean = with(apop, weighted.mean(x = SP.RUR.TOTL.ZS, w = SP.POP.TOTL))

## Extract top country assume they are never missing
top = 30
mainCountry = apop[1:top, ]
otherCountry = apop[(top + 1):NROW(apop), ]


## Start the simulation based on maximum 50% missing
n.bt = 5000
missingMean = double(n.bt)
missingProp = double(n.bt)
for(i in 1:n.bt){
  missingProp[i] = runif(1, min = 0, max = 0.5)
  ind = sample(1:NROW(otherCountry), NROW(otherCountry) * (1 - missingProp[i]),
    replace = FALSE)
  tmp.df = rbind(otherCountry[ind, ], mainCountry)
  missingMean[i] = weighted.mean(x = tmp.df[, "SP.RUR.TOTL.ZS"],
              w = tmp.df[, "SP.POP.TOTL"], na.rm = TRUE)
}


hist(missingMean, breaks = 200)
abline(v = trueMean, col = "red")
summary(missingMean)


## Start the simulation based on maximum 20% missing
n.bt = 5000
missingMean = double(n.bt)
missingProp = double(n.bt)
for(i in 1:n.bt){
  missingProp[i] = runif(1, min = 0, max = 0.2)
  ind = sample(1:NROW(otherCountry), NROW(otherCountry) * (1 - missingProp[i]),
    replace = FALSE)
  tmp.df = rbind(otherCountry[ind, ], mainCountry)
  missingMean[i] = weighted.mean(x = tmp.df[, "SP.RUR.TOTL.ZS"],
              w = tmp.df[, "SP.POP.TOTL"], na.rm = TRUE)
}

hist(missingMean, breaks = 200)
abline(v = trueMean, col = "red")
summary(missingMean)



## Now investigate urban population
## Download data
ugdp = WDI(indicator = "NY.GDP.PCAP.CD", start = 2010, end = 2010)
ugdp = na.omit(ugdp)

## Merge with gdpulation
agdp = merge(ugdp, pop)

## Order data
agdp = agdp[order(agdp$SP.POP.TOTL, decreasing = TRUE),]

## Compute true value
trueMean = with(agdp, weighted.mean(x = NY.GDP.PCAP.CD, w = SP.POP.TOTL))

## Extract top country assume they are never missing
agdp = agdp[order(agdp$NY.GDP.PCAP.CD, decreasing = TRUE),]
mainCountry = agdp[agdp$SP.POP.TOTL >= 5000000 | agdp$NY.GDP.PCAP.CD >= 3000, ]
otherCountry = agdp[!(agdp$SP.POP.TOTL >= 5000000 |
  agdp$NY.GDP.PCAP.CD >= 3000), ]


## Start the simulation based on maximum 50% missing
n.bt = 5000
missingMean = double(n.bt)
missingProp = double(n.bt)
for(i in 1:n.bt){
  missingProp[i] = runif(1, min = 0, max = 0.5)
  ind = sample(1:NROW(otherCountry), NROW(otherCountry) * (1 - missingProp[i]),
    replace = FALSE)
  tmp.df = rbind(otherCountry[ind, ], mainCountry)
  missingMean[i] = weighted.mean(x = tmp.df[, "NY.GDP.PCAP.CD"],
              w = tmp.df[, "SP.POP.TOTL"], na.rm = TRUE)
}


hist(missingMean, breaks = 200)
abline(v = trueMean, col = "red")
summary(missingMean)



## Start the simulation based on maximum 50% missing
n.bt = 5000
missingMean = double(n.bt)
missingProp = double(n.bt)
for(i in 1:n.bt){
  missingProp[i] = runif(1, min = 0, max = 0.2)
  ind = sample(1:NROW(otherCountry), NROW(otherCountry) * (1 - missingProp[i]),
    replace = FALSE)
  tmp.df = rbind(otherCountry[ind, ], mainCountry)
  missingMean[i] = weighted.mean(x = tmp.df[, "NY.GDP.PCAP.CD"],
              w = tmp.df[, "SP.POP.TOTL"], na.rm = TRUE)
}


hist(missingMean, breaks = 200)
abline(v = trueMean, col = "red")
summary(missingMean)

## Investigate why the population is so different

########################################################################
## Title: investigate the speed of mean and .Internal(mean)
## Date: 12-09-2012
########################################################################

tf1 = function(x) mean(x)
tf2 = function(x) .Internal(mean(x))

library(rbenchmark)

benchmark(
  tf1 = tf1(1:10000),
  tf2 = tf2(1:10000),
  replications = 10^(4:6),
  order = c("replications", "elapsed")
  )

## It is slightly faster, but only for very large N

########################################################################
## Title: Replacement by Mean or by Median
## Date: 13-09-2012
########################################################################

library(FAOSYB)

## Download data
test = getWDI()

## Remove regions
test = test[-grep("income|devel|OECD|[E|e]urope|[A|a]sia|[A|a]Africa|[P|p]acific|[N|n]meriaUN|World|Heav|Euro|Not", test$Country), ]
pop = subset(test, select = "SP.POP.TOTL", subset = Year == 2010, drop = TRUE)

## obtain real sum
trueSum = sum(pop)

## Assume that the top 50 countries are reported
pop = sort(pop, decreasing = TRUE)
mainPop = sum(pop[1:50])
otherPop = pop[51:length(pop)]

## Simulation with replace by mean
n.sim = 50000
simSum = double(n.sim)
for(i in 1:length(simSum)){
  missProp = runif(1, 0.7, 1)
  tmp = otherPop
  ind = sample(1:length(otherPop), missProp * length(otherPop), replace = FALSE)
  impValue = mean(tmp[ind])
  tmp[-ind] = impValue
  simSum[i] = sum(tmp)
}
totalSimPopMean = mainPop + simSum

## Plot simulation
hist(totalSimPopMean, breaks = 300)
abline(v = trueSum, col = "red")
abline(v = mean(totalSimPopMean), col = "blue")

## Simulation with replacement by median
n.sim = 50000
simSum = double(n.sim)
for(i in 1:length(simSum)){
  missProp = runif(1, 0.7, 1)
  tmp = otherPop
  ind = sample(1:length(otherPop), missProp * length(otherPop), replace = FALSE)
  impValue = median(tmp[ind])
  tmp[-ind] = impValue
  simSum[i] = sum(tmp)
}

totalSimPopMed = mainPop + simSum
hist(totalSimPopMed, breaks = 300)
abline(v = trueSum, col = "red")
abline(v = mean(totalSimPopMed), col = "blue")

## Distribution of the true population
hist(pop, breaks = 100)
abline(v = mean(pop), col = "red")
abline(v = median(pop), col = "blue")

## The simulation shows that replacement by mean is asymptotic unbiased
## while replacement by median is biased significantly. However, this
## only holds if the data is symmetric and MCAR.




pop = subset(country.df, select = c(FAOST_CODE, P1.DEM.UN.WPP.POP.TOT),
  subset = Year == 2010, drop = TRUE)

dt2010.df = subset(country.df, subset = Year == 2010)
dt2010.df = dt2010.df[, -c(1:2)]

dt2010.df[which(dt2010.df == Inf, arr.ind = TRUE)[,1],
          which(dt2010.df == Inf, arr.ind = TRUE)[,2]] = NA

dt2010.df = dt2010.df[,
  !sapply(dt2010.df, function(x) sum(is.na(x)) == length(x))]

pdf(file = "dataHist.pdf")
for(i in 1:NCOL(dt2010.df)){
  hist(dt2010.df[, i], breaks = 100, main = colnames(dt2010.df)[i])
  abline(v = mean(dt2010.df[, i], na.rm = TRUE), col = "red")
  legend("topleft", legend = round(sum(na.omit(dt2010.df[, i]) <
                      mean(dt2010.df[, i], na.rm = TRUE))/
         length(na.omit(dt2010.df[, i])), 2))
}
graphics.off()


## Read in data from old SYB
p1d.df = read.csv(file = "~/Dropbox/AdaMat/SYB/Output/csv/P1_DEM.csv",
  header = TRUE)
p1m.df = read.csv(file = "~/Dropbox/AdaMat/SYB/Output/csv/P1_MAC.csv",
  header = TRUE)
p1r.df = read.csv(file = "~/Dropbox/AdaMat/SYB/Output/csv/P1_RES.csv",
  header = TRUE)
p2.df = read.csv(file = "~/Dropbox/AdaMat/SYB/Output/csv/P2.csv", header = TRUE)
p3.df = read.csv(file = "~/Dropbox/AdaMat/SYB/Output/csv/P3.csv", header = TRUE)
p4.df = read.csv(file = "~/Dropbox/AdaMat/SYB/Output/csv/P4.csv", header = TRUE)

all.df = Reduce(merge, list(p1d.df, p1m.df, p1r.df, p2.df, p3.df, p4.df))

pdf(file = "dataHist.pdf")
par(mfrow = c(2, 2))
for(i in 12:NCOL(all.df)){
  if(mode(all.df[, i]) == "numeric"){
    hist(all.df[, i], breaks = 100, main = colnames(all.df)[i])
    abline(v = mean(all.df[, i], na.rm = TRUE), col = "red")
    legend("topleft", legend = round(sum(na.omit(all.df[, i]) <
                        mean(all.df[, i], na.rm = TRUE))/
                        length(na.omit(all.df[, i])), 2))
  }
}
graphics.off()
system("evince dataHist.pdf&")

## It seems there is still large amount of highly skewed data.


## Read the name file, lets lets just work with the first 100 variables
library(WDI)
WDI = read.csv(file = "http://dl.dropbox.com/u/18161931/WorldBankIndicators.csv",
  stringsAsFactors = FALSE)
WDI = WDI[-c(1:10), ]

## Download and merge the data. Some vairables are not collected in 2010
## and thus they are discarded
WDI.df = WDI(indicator = WDI$WDI_NAME[1], start = 2000, end = 2010)
for(i in 2:NROW(WDI)){
  print(i)
  tmp = WDI(indicator = WDI$WDI_NAME[i], start = 2000, end = 2010)
  if(!inherits(tmp, "try-error") &
     (sum(is.na(tmp[, WDI$WDI_NAME[i]])) != NROW(tmp)))
    print(str(tmp))
    WDI.df = merge(WDI.df, tmp, by = c("iso2c", "country", "year"),
      all = TRUE)
}

## Produce histogram to examine the suitability of mean imputation
pdf(file = "dataDist.pdf")
for(i in 3:NCOL(WDI.df)){
  hist(WDI.df[, i], breaks = 100, main = colnames(WDI.df)[i], xlab = NULL)
  abline(v = mean(WDI.df[, i], na.rm = TRUE), col = "red")
  pctBelowMean = round(100 * sum(na.omit(WDI.df[, i]) <
    mean(WDI.df[, i], na.rm = TRUE))/length(na.omit(WDI.df[, i])), 2)
  legend("topright", legend = paste(pctBelowMean,
                       "% of data are below the mean", sep = ""))
}
graphics.off()

sWDI.df = WDI.df[-grep("income|devel|OECD|[E|e]urope|[A|a]sia|[A|a]Africa|[P|p]acific|[N|n]meriaUN|World|Heav|Euro", WDI.df$country), ]

sparsity.df = cbind(sWDI.df[, 1:3], sparsity = apply(sWDI.df[, -c(1:3)], 1,
                              function(x) sum(is.na(x))/length(x)))

asp.df = arrange(sparsity.df, year, sparsity)

ggplot(asp.df, aes(x = year, y = sparsity)) +
  geom_line(aes(col = country))

spar.df = with(asp.df, aggregate(sparsity, list(country), mean))
colnames(spar.df) = c("country", "sparsity")
spar.df = arrange(spar.df, sparsity)

########################################################################
## Title: Text mining and world cloud
## Date: 14-09-2012
########################################################################

#install.packages(c("wordcloud","tm"),repos="http://cran.r-project.org")
library(wordcloud)
library(tm)

## Word cloud 1
wordcloud("May our children and our children's children to a
thousand generations, continue to enjoy the benefits conferred
upon us by a united country, and have cause yet to rejoice under
those glorious institutions bequeathed us by Washington and his
compeers.", colors = brewer.pal(6, "Dark2"), random.order = FALSE)

## Word cloud 2
data(SOTU)
SOTU <- tm_map(SOTU, function(x) removeWords(tolower(x), stopwords()))
wordcloud(SOTU, colors = brewer.pal(6, "Dark2"), random.order = FALSE)

## A big improvement of the version 2.2 wordcloud is to avoid text overlay
states <- c('Alabama', 'Alaska', 'Arizona', 'Arkansas',
	'California', 'Colorado', 'Connecticut', 'Delaware',
	'Florida', 'Georgia', 'Hawaii', 'Idaho', 'Illinois',
	'Indiana', 'Iowa', 'Kansas', 'Kentucky', 'Louisiana',
	'Maine', 'Maryland', 'Massachusetts', 'Michigan', 'Minnesota',
	'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada',
	'New Hampshire', 'New Jersey', 'New Mexico', 'New York',
	'North Carolina', 'North Dakota', 'Ohio', 'Oklahoma', 'Oregon',
	'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota',
	'Tennessee', 'Texas', 'Utah', 'Vermont', 'Virginia', 'Washington',
	'West Virginia', 'Wisconsin', 'Wyoming')

loc <- rmvnorm(50,c(0,0),matrix(c(1,.7,.7,1),ncol=2))

plot(loc[,1],loc[,2],type="n")
text(loc[,1],loc[,2],states)

textplot(loc[,1],loc[,2],states)

plot(loc[,1],loc[,2],type="n")
nc <- wordlayout(loc[,1],loc[,2],states,cex=50:1/20)
text(nc[,1] + .5*nc[,3],nc[,2]+.5*nc[,4],states,cex=50:1/20)


########################################################################
## Title: Singular value decomposition
## Date: 14-09-2012
########################################################################

## SVD to find a generalized inverse of a non-full-rank matrix
library(MASS)

a <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1), 9, 4)

a.svd <- svd(a)
a.svd$d

ds <- diag(1/a.svd$d[1:3])
u <- a.svd$u
v <- a.svd$v
us <- as.matrix(u[, 1:3])
vs <- as.matrix(v[, 1:3])

## A = v %*% d %*% u (generalized inverse)
(a.ginv <- vs %*% ds %*% t(us))

## Same as a (ABA = A)
round(a %*% a.ginv %*% a)

# using the function ginv defined in MASS
ginv(a)


## Image processing
library(ReadImages)
x = read.jpeg("pansy.jpg")

plot(x, useRaster = TRUE)

r <- imagematrix(x, type = "grey")

plot(r, useRaster = TRUE)


r.svd <- svd(r)
d <- diag(r.svd$d)
dim(d)

u <- r.svd$u
v <- r.svd$v
plot(1:length(r.svd$d), r.svd$d)

# first approximation
u1 <- as.matrix(u[-1, 1])
v1 <- as.matrix(v[-1, 1])
d1 <- as.matrix(d[1, 1])
l1 <- u1 %*% d1 %*% t(v1)
l1g <- imagematrix(l1, type = "grey")
plot(l1g, useRaster = TRUE)

# more approximation
depth <- 5
us <- as.matrix(u[, 1:depth])
vs <- as.matrix(v[, 1:depth])
ds <- as.matrix(d[1:depth, 1:depth])
ls <- us %*% ds %*% t(vs)
lsg <- imagematrix(ls, type = "grey")

plot(lsg, useRaster = TRUE)


# Last approximation
depth <- 100
ufs <- as.matrix(u[, 1:depth])
vfs <- as.matrix(v[, 1:depth])
dfs <- as.matrix(d[1:depth, 1:depth])
lfs <- ufs %*% dfs %*% t(vfs)
lfsg <- imagematrix(lfs, type = "grey")

plot(lfsg, useRaster = TRUE)
sum(c(object.size(ufs), object.size(vfs), object.size(dfs)))
dev.new()
plot(imagematrix(u %*% d %*%t(v), type = "grey"))
sum(c(object.size(u), object.size(v), object.size(d)))

## Although the plots only look minorly different through naked eye, but
## the size of the object is only 1/5th


## Test on my own profile pic
test = read.jpeg(file = "profilePic.jpg")

plot(test, useRaster = TRUE)

svdApprox = function(x, prop = 0.5){
  svd.x = svd(x)
  l = prop * length(svd.x$d)
  d = diag(svd.x$d)
  print(l)
  v = svd.x$v[, 1:l]
  u = svd.x$u[, 1:l]
  d = d[1:l, 1:l]
  approx = u %*% d %*% t(v)
  imagematrix(approx)
}

test2 = svdApprox(test[, , 1], prop = 0.01)
test3 = svdApprox(test[, , 2], prop = 0.01)
test4 = svdApprox(test[, , 3], prop = 0.01)

test5 = abind(test2, test3, test4, along = 3)

## An approximation of the plot
plot(imagematrix(test5), useRaster = TRUE)


########################################################################
## Title: Simulation of preferential attachment
## Date: 15-09-2012
########################################################################


library(animation)
saveHTML({
  N = 1000
  for(i in 1:100){
    m = 100 * i
    rlt = rep(1, N)
    for(i in 1:m){
      selection = sample(x = 1:N, size = 1, prob = rlt)
      rlt[selection] = rlt[selection] + 1
    }
    rlt = rlt - 1
        par(mfrow = c(1, 2))
    bin = hist(rlt, breaks =  m/10,
      main = paste("Histogram of preferential network\n with ", m,
        " individuals choosing\nbetween 1000 clubs", sep = ""),
      xlab = "Number of Individuals", ylab = "Frequency")
    plot(log(bin$breaks[-1]), log(bin$counts),
         main = "log-log histogram", xlab = "", ylab = "")
    lm.df = data.frame(lbin = log(bin$breaks[-1]), lcounts = log(bin$counts))
    lm.df = lm.df[lm.df$lcounts != -Inf & lm.df$lcounts != 0 , ]
    fit = lm(lcounts ~lbin, data = lm.df)
    abline(coef = coef(fit))
  }
})

########################################################################
## Title: Plotting network
##        (http://www.ats.ucla.edu/stat/r/faq/snplot.htm)
## Date: 09-15-2012
########################################################################

library(igraph)
x <-read.table("http://www.ats.ucla.edu/stat/r/faq/mat25.txt", header=FALSE)

network = as.matrix(x)
g1 = graph.adjacency(network)



test[sample(1:10000, 300, replace = FALSE)] = 1
gtest = graph.adjacency(test)
plot(gtest)

## Lets grow a erdos-renyi network

n.vertex = 50
n.edge = 100
system("rm /tmp/RtmpVNNohv/images/*")
test = matrix(0, nc = n.vertex, nr = n.vertex)
edge = sample(1:vertex^2, n.edge, replace = FALSE)
saveHTML({
  for(i in 1:n.edge){
    if(i == 1)
      plot.new()
    test[edge[i]] = 1
    gtest = graph.adjacency(test, mode = "undirected")
    plot(gtest)
  }
})

## Compute betweeness
b1 = betweenness(g1, directed = FALSE)
b1

# compute closeness
c1 = closeness(g1, mode="out")
c1

# compute degree
d1 = degree(g1, mode="out")

plot(g1)


test = as.matrix(expand.grid(1:100, 1:100))
test = test[-sample(1:10000, 2000, replace = FALSE), ]
test = test[test[, 1] != test[, 2], ]


test.gh = graph(test, directed = FALSE)
plot(test.gh)


xlist<-read.graph("http://www.ats.ucla.edu/stat/r/faq/elist1.txt",
                  format="edgelist")
xlist

########################################################################
## Title: Heat map of NBA player profiles
## Date: 20-09-2012
########################################################################

nba = read.csv("http://datasets.flowingdata.com/ppg2008.csv")
nba$Name <- with(nba, reorder(Name, PTS))

library(ggplot2)
library(reshape)
nba.m <- melt(nba)
nba.m <- ddply(nba.m, .(variable), transform,
               rescale = scale(value))

p <- ggplot(nba.m, aes(variable, Name)) + geom_tile(aes(fill = rescale),
     colour = "white") + scale_fill_gradient(low = "white",
     high = "steelblue")

base_size <- 9
p + theme_grey(base_size = base_size) + labs(x = "",
    y = "") + scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + opts(legend.position = "none",
    axis.ticks = theme_blank(), axis.text.x = theme_text(size = base_size *
        0.8, angle = 330, hjust = 0, colour = "grey50"))


########################################################################
## Title: Boxplot VS barplot
## Date: 22-09-2012
########################################################################

library(FAOSYB)
pop = getWDI()

test = aggRegion(pop,
  relationDF = FAOcountryProfile[, c("ISO2_WB_CODE", "MOTHER_M49_CODE")],
  aggVar = "SP.POP.TOTL")


test.df = data.frame(reg = rep(letters[1:5], 200), norm = rnorm(100),
  exp = rexp(1000, rate = 0.3))

agg.df = ddply(.data = test.df, .variables = .(reg),
  .fun = function(x) colMeans(x[, c("norm", "exp")]))

magg.df = melt(agg.df, "reg")
colnames(magg.df) = c("reg", "type", "value")


ggplot(magg.df, aes(x = reg, y = value)) + geom_bar() +
  facet_wrap(~type, scales = "free_y")

dev.new()
mtest.df = melt(test.df, "reg")
colnames(mtest.df) = c("reg", "type", "value")
ggplot(mtest.df, aes(x = reg, y = value)) + geom_boxplot() + geom_jitter() +
  facet_wrap(~type, scales = "free_y")


test = rexp(10000, rate = 1/10)

########################################################################
## Title: Simple imputation by mean
## Date: 24-09-2012
########################################################################

library(plyr)
library(reshape)

test.df = data.frame(rtCode = c(1, 1, 2, 2, 2),
  cmdCode = c(1001, 1001, 10011000, 10011000, 10011000),
  trade = rexp(5))
test.df[c(2, 5), "trade"] = NA
test.df

meanImputation = function(data, group, var){
  ## Compute the mean
  mean.df = ddply(.data = data, .variable = group,
    .fun = function(x) mean(x[, var], na.rm = TRUE))
  colnames(mean.df)[3] = "Mean"

  ## merge back
  full.df = merge(data, mean.df)

  ## Replace missing value with group mean
  full.df[is.na(full.df[, var]), var] = full.df[is.na(full.df[, var]), "Mean"]
  full.df
}

meanImputation(test.df, group = c("rtCode", "cmdCode"), var = "trade")


########################################################################
## Title: Linear mix model
## Date: 25-09-2012
########################################################################

library(MASS)
data(oats)
names(oats) = c('block', 'variety', 'nitrogen', 'yield')
oats$mainplot = oats$variety
oats$subplot = oats$nitrogen

library(nlme)
m1.nlme = lme(yield ~ variety*nitrogen,
                      random = ~ 1|block/mainplot,
                      data = oats)

summary(m1.nlme)


########################################################################
## Title: Partial matching string
## Date: 29-09-2012
########################################################################

## look up matches of one dataframe in another dataframe.
## the strings to be matched are comprised of 1 or more words
## and seperated by white space.
## method: match strings that have the highest fraction of words that match up

d1 <- read.csv("http://s.telegraph.co.uk/graphics/conrad/PercentageUsingTheNet.csv",
               header = T, sep = ",", encoding = "UTF-8")
d2 <- read.csv("http://www.iso.org/iso/country_names_and_code_elements_txt",
               header = T, sep = ";", encoding = "UTF-8")

## strings to be compared d1$ECONOMY and d2$Country.Name
mystr.1 <- as.character(d1$ECONOMY)
mystr.2 <- as.character(d2$Country.Name)
mystr.3 <- as.character(d2$ISO.3166.1.alpha.2.code)

## remove punctuation and multiple spaces
mystr.1 <- tolower(gsub("[^[:alnum:][:space:]]", "", mystr.1))
mystr.1 <- gsub("\\s+", " ", mystr.1)
mystr.2 <- tolower(gsub("[^[:alnum:][:space:]]", "", mystr.2))
mystr.2 <- gsub("\\s+", " ", mystr.2)

mystr.1.spl <- strsplit(mystr.1, " ")
(mystr.1 <- unlist(lapply(mystr.1.spl, function(x) paste(x, collapse = " "))))

mystr.2.spl <- strsplit(mystr.2, " ")
(mystr.2 <- unlist(lapply(mystr.2.spl, function(x) paste(x, collapse = " "))))

## function that finds matching words in string (words seperated by single space!)
n.wordshared <- function(x, y) {
    sum(!is.na(match(unlist(strsplit(x, " ")),
                     unlist(strsplit(y, " ")))
         )
        )
    }

## example
n.wordshared(x = "hello world", y = "good bye world")

## function that calculates fraction of shared words
fr.wordshared <- function(x, y) {
                     n.wordshared(x, y) / (length(unique(unlist(strsplit(x, " "))))
                                           + length(unique(unlist(strsplit(y, " ")))))
                          }

## example
fr.wordshared(x = "hello world", y = "good bye world")

mydf <- data.frame(str1 = mystr.1, mymatch = "", match.iso = "",
                   stringsAsFactors = F)

## now look up every element of string 1 in string 2
## and if there are matching words assign match to dataframe
for (i in 1:nrow(mydf)) {
   xx <- sapply(mystr.2, fr.wordshared, y = mystr.1[i])
   if (sum(xx) == 0) {
     mydf$mymatch[i] <- NA
     mydf$match.iso[i] <- NA
     } else {
     mydf$mymatch[i] <- paste(names(which(xx == max(xx))), collapse = "; ")
     mydf$match.iso[i] <- paste(mystr.3[as.integer(which(xx == max(xx)))], collapse = "; ")
   }
}

## see result
print(mydf)

## these are the multiple matches
(aa <- mydf[grep(";", mydf$mymatch), ])

## these were not matched
(bb <- mydf[is.na(mydf$mymatch), ])

########################################################################
## Title: Power law and log-log distribution
## Date: 03-10-2012
########################################################################

x = 2^(0:100)
y1 = x^-2
y2 = x^-2.1

plot(x, 1y)
plot(log(x), log(y1))
lines(log(x), log(y2))


########################################################################
## Title: Mix models
## Date: 11-10-2012
########################################################################

library(lme4)
library(MEMSS)

print(dotplot(reorder(Rail, travel) ~ travel, Rail, xlab = "Travel time (ms)",
              ylab = "Rail"))
Rm1ML <- lmer(travel ~ 1 + (1|Rail), Rail, REML = FALSE, verbose = TRUE)
Rm1ML

## This is the model matrix Z
Rm1ML@Zt

## The deviances
Rm1ML@deviance



## multiple random effect model
print(sm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))

## Extract the variance and covariance matrix
show(vc <- VarCorr(sm1))


## Model with two nested grouping factor
print(Om1 <- lmer(yield ~ nitro + Variety + (1|Block/Variety), Oats),
      corr = FALSE)

print(Om1a <- lmer(yield ~ nitro + (1|Block/Variety), Oats), corr = FALSE)

print(sm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy),
      corr = FALSE)

print(Om2 <- lmer(yield ~ nitro + (1|Variety:Block) + (nitro|Block), Oats),
      corr = FALSE)



########################################################################
## Title: Test multivaritate state space model
## Date: 15-10-2012
########################################################################

library(FAOSYB)
library(MARSS)
library(reshape2)
pop = getWDItoSYB()$data
## colnames(pop)[4] = "pop"
cpop = cast(data = pop, formula = Country + ISO2_WB_CODE~ Year)

pop.mat = as.matrix(cpop[, -c(1:3)])

## Takes way too long to fit.
pop.ss = MARSS(pop.mat)

B1 = matrix(list("b", 0, 0, "b"), 2, 2)
U1 = matrix(0, 2, 1)
Q1 = matrix(c("q11", "q12", "q12", "q22"), 2, 2)
Z1 = matrix(c(1, 0, 1, 0, 1, 0), 3, 2)
A1 = matrix(list("a1", 0, 0), 3, 1)
R1 = matrix(list("r1", 0, 0, 0, "r2", 0, 0, 0, "r2"), 3, 3)
pi1 = matrix(0, 2, 1)
V1 = diag(1, 2)

model.list = list(B = B1, U = U1, Q = Q1, Z = Z1,
  A = A1, R = R1, x0 = pi1, V0 = V1)

y = matrix(c(1, 2, NA, NA, 3.2, 8, 2, 5, 2, NA, 5.1, 5, 1, NA, 2, 2.2, NA, 7),
  nr = 3, byrow = TRUE)


y.fit = MARSS(y = y, model = model.list)



data(graywhales)
years = graywhales[, 1]
loggraywhales = log(graywhales[,2])
kem = MARSS(loggraywhales)



u <- rnorm(25)
myMod <- dlmModReg(u, dV = 14.5)

buildFun <- function(x) {
  m <- dlmModPoly(1, dV = exp(x[1]))
  m$JW <- matrix(1)
  m$X <- matrix(exp(x[2]), nc = 1, nr = length(Nile))
  j <- which(time(Nile) == 1899)
  m$X[j,1] <- m$X[j,1] * (1 + exp(x[3]))
  return(m)
}
fit <- dlmMLE(Nile, parm = c(0,0,0), build = buildFun)

dlmNileJump <- buildFun(fit$par)

nileJumpFilt <- dlmFilter(Nile, dlmNileJump)
nileJumpSmooth = dlmSmooth(Nile, dlmNileJump)
plot(Nile, type = 'o', col = "seagreen")
lines(dropFirst(nileJumpFilt$m), type = 'o', pch = 20, col = "brown")
lines(nileJumpSmooth$s, type = 'o', pch = 20, col = "purple")


popFun = function(x){
  dlmModARMA(ar = x[1], sigma2 = x[2], m0 = x[3])
}

testData = na.omit(pop[pop$Country == "China", "SP.POP.0014.TO.ZS",
  drop = TRUE])
plot(testData, type = "o", col = "seagreen")


splitFirstObs = function(x){
  x[is.na(x)] = 0
  cx = cumsum(x)
  cx == 0
}

t = 1:length(testData)
for(i in 1:100){
  tmp = testData
  missInd = sample(1:length(tmp), 20, replace = FALSE)
  tmp[missInd] = NA
  ava.t = !splitFirstObs(tmp)
  tmp = tmp[ava.t]
  testEst = dlmMLE(y = tmp,
    parm = c(0, var(testData),
      ifelse(is.na(tmp[1]), mean(tmp, na.rm = TRUE), tmp[1])),
    build = popFun)
  testFit = popFun(testEst$par)
  testFilter = dlmFilter(tmp, testFit)
  lines(t[ava.t], testFilter$f, type = "o", pch = 20, col = "brown")
}


########################################################################
## Title: Test Matrix creation speed
## Date: 23-10-2012
########################################################################

doInstall <- TRUE  # Change to FALSE if you don't want packages installed.
toInstall <- c("rbenchmark", "ggplot2")
if(doInstall){install.packages(toInstall, repos = "http://cran.r-project.org")}
lapply(toInstall, library, character.only = TRUE)

simpleLoop <- function(vec = 1:10, nc = 5){  # Expected to be slow
  out <- c()
  for(ii in 1:nc){
    out <- cbind(out, vec)
    }
  return(out)
  }

preAllocated <- function(vec = 1:10, nc = 5){  # Expected to be faster
  out <- matrix(NA, ncol = nc, nrow = length(vec))
  for(ii in 1:nc){
    out[, ii] <- vec
    }
  return(out)
  }

addSweep <- function(vec = 1:10, nc = 5){  # "Sweeping" the vector through
  out <- matrix(0, ncol = nc, nrow = length(vec))
  sweep(out, 1, vec, "+")
  }

repVector <- function(vec = 1:10, nc = 5){  # Expected to be fastest
  replicate(nc, vec)
  }

outerProduct <- function(vec = 1:10, nc = 5){  # Just figured this one out,
  outer(vec, rep(1, nc))                       # matrix algebra!
  }

matCreate = function(vec, nc = 5){
    matrix(vec, nc = nc)
}

simpleLoop(1:15, 15)  # Examples
repVector(1:15, 15)
outerProduct(1:15, 15)  # etc...

setVector <- 1:100  # These parameters pass onto benchmark test
setColumns <- 100

# Benchmark testing
benchmarkResults <- within(benchmark(
    outerProduct = outerProduct(setVector, setColumns),
    repVector = repVector(setVector, setColumns),
    addSweep = addSweep(setVector, setColumns),
    preAllocated = preAllocated(setVector, setColumns),
    simpleLoop = simpleLoop(setVector, setColumns),
    matCreate = matCreate(setVector, setColumns),
    replications = 10 ^ (3:5),
    columns = c('test', 'replications', 'elapsed', 'relative'),
    order = c('relative', 'elapsed', 'replications')),
                           { average = elapsed / replications })

benchmarkResults


########################################################################
## Title: Spatial data analysis
## Date: 29-10-2012
########################################################################

## Read the data
cran.df = read.table("http://www.asdar-book.org/datasets/CRAN051001a.txt",
    header = TRUE)
cran.mat = with(cran.df, cbind(long, lat))
row.names(cran.mat) = 1:nrow(cran.mat)

## Create the spatial object
llCRS = CRS("+proj=longlat +ellps=WGS84")
cran.sp = SpatialPoints(cran.mat, proj4string = llCRS)

## Create spatial data frame
cran.spdf = SpatialPointsDataFrame(cran.mat, cran.df, proj4string = llCRS,
    match.ID = TRUE)



########################################################################
## Title: plotGoogleMaps
## Date: 9-11-2012
########################################################################

library(plotGoogleMaps)

data(meuse)
coordinates(meuse) = ~x+y
proj4string(meuse) = CRS("+init=epsg:28992")
m = bubbleGoogleMaps(meuse, zcol = "zinc", filename = "myMap.htm")
m2 = bubbleGoogleMaps(meuse, zcol = "lead", filename = "myMap.htm")
m = bubbleGoogleMaps(meuse, zcol = "cadmium", filename = "myMapCadmium.htm",
    layerName = "Bubble plot - meuse", colPalette = terrain.colors(5),
    strokeColor="")


ell = data.frame(E = c(7456263,7456489,7456305),
                 N = c(4954146,4952978,4952695),
                 A = c(2.96,4.55,7.10),
                 B = c(2.35,2.11,2.29),
                 teta = c(28.35242,41.04491,38.47216))

coordinates(ell) = ~E+N
proj4string(ell) = CRS("+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m")

m = ellipseGoogleMaps(ell, filename = "Ellipse.htm", mapTypeId = "ROADMAP")


library(RgoogleMaps)
require(RColorBrewer)
plotclr <-brewer.pal(8,"YlOrRd")
plotclrAlpha = AddAlpha(plotclr,0.5)
plot(1:16, col = c(plotclr, plotclrAlpha), pch = 19)


data(NYleukemia)
population =  NYleukemia$data$population
cases =  NYleukemia$data$cases
mapNY =  GetMap(center=c(lat=42.67456, lon=-76.00365),
    destfile = "NYstate.png", maptype = "mobile", zoom=9)
ColorMap(100*cases/population, mapNY, NYleukemia$spatial.polygon,
         add = FALSE,alpha = 0.35, log = TRUE, location = "topleft")


test = GetMap(center = c(lat = 0, lon = 0),
    zoom = 1)
ColorMap(100*cases/population, test, NYleukemia$spatial.polygon,
         add = FALSE,alpha = 0.35, log = TRUE, location = "topleft")

library(ggmap)



########################################################################
## Title: More on data.table
## Date: 2012-12-03
########################################################################

library(data.table)

## Example data from ChainLadder package
## Data is also stored in a Google Fusion table
url <- paste0("http://www.google.com/fusiontables/api/query?",
              "sql=SELECT+*+FROM+1SL7c4TwyI1YxuQELc0R3PjsYC3TwhP3o7k_NZzc")
myData <- read.csv(url, stringsAsFactors=TRUE)
head(myData)

myData <- data.table(myData)

## Add a column with cumulative values. No copying required as I use
## the ':=' operator, which works by reference
myData[order(dev), cvalue:=cumsum(value), by=list(origin, lob)]

## Plot for each line of business, grouped by the
## origin year against development year
require(lattice)
IncrPlot <- xyplot(value/1e3 ~ dev | lob, groups=origin,
                   data=myData, main="Incremental developments",
                   scales="free", t="l", as.table=TRUE, layout=c(3,4))

CumPlot <- xyplot(cvalue/1e3 ~ dev | lob, groups=origin,
                  data=myData, main="Cumulative developments",
                  scales="free", t="l", as.table=TRUE, layout=c(3,4))

print(IncrPlot, position=c(0, 0, 0.5, 1), more=TRUE)
print(CumPlot, position=c(0.5, 0, 1, 1))

## Get latest development period for each origin year and line of business
latestData <- myData[, .SD[max(dev)] , by=list(origin, lob) ]
head(latestData)



## Compute new columns
head(myData[, test2 := cumsum(value), by = list(lob, origin)], 30)
head(myData[, test3 := min(na.omit(value)), by = list(origin)], 30)

head(myData[, list(origin)])

myData[, list(mean = mean(value), sd = sd(value)), by = list(origin)]

## Juset use subset to select the columns within a function
myData[0]





########################################################################
## Title: Analysis on Ethanol price
## Date: 2012-10-12
########################################################################

## Load the data
price = c(256.4539509, 253.0103323, 238.4680218, 221.3614295,
204.8028305, 265.7775962, 311.9164897, 308.6385168, 302.0694839,
318.3226139, 352.268386, 360.9401209, 375.4628175, 378.4099161,
383.9545443, 392.784195, 397.5325892, 428.1213504, 438.918173,
427.636953, 457.0196865, 399.4136503, 395.5303745, 341.6413987,
313.0843292, 314.4814121, 285.3876799, 233.3696616, 119.9630217,
141.7500589, 155.5623494, 196.3637513, 229.4625896, 251.289821,
213.8123262, 167.3298679, 197.1651085, 192.2812134, 196.2751106,
158.7049077, 212.5895108, 294.058646, 312.2727467, 306.9955505,
300.3393169, 246.584474, 202.4623712, 209.2861045, 253.7511122,
231.6354806, 212.4739085, 238.4443347, 274.9498517, 314.1621453,
334.5082248, 308.4537737, 306.0770432, 293.0944382, 288.5887663,
260.41877, 220.4400935, 224.4706271, 205.5273835, 208.1694923,
182.0554726, 266.0416519, 309.1569521, 318.6518091, 296.2521305,
235.7184512, 249.0379145, 280.3955744, 266.9097046, 257.4328707,
224.1437189, 217.7147482, 212.2407043, 302.7581227, 309.7844237,
278.8109404, 272.1748998, 247.3041921, 231.2236569, 226.7766045,
239.8048683, 223.2149923, 224.8639848, 234.0798574, 256.3646129,
240.5365164, 221.0765932)

## Assuming that the latest data refering to October 2012
price.ts = ts(price, frequency = 12, start = c(2005, 4))

## There seems to be a collapse in price in 2007 and some seasonality
plot(price.ts)

## use STL with a periodic frame of 12 to decompose the time series,
## the trend also supports the collpase of the margin in 2007.
price.stl = stl(price.ts, s.window = 12)
plot(price.stl)

## It looks like that the error can be bi-modal, however the important
## part is that the market can take very big hit sa shown by the large
## negative return.
hist(price.stl$time.series[, "remainder"], breaks = 30)


########################################################################
## Title: TFX: An R Interface to the TrueFX Web API
## Date: 2012-12-11
## Source:http://rpubs.com/gsee/TFX
########################################################################

install.packages("TFX", repos = "http://cran.stat.auckland.ac.nz")

library(TFX)
QueryTrueFX()

## Note: if you don't see millisecond resolution, set
##       options(digits.secs=3) or higher.


########################################################################
## Title: Scraping functions from github
## Date: 2012-12-11
########################################################################

## Function for sourcing github codes
source_github = function(url, ...) {
    require(RCurl)
    sapply(c(url, ...), function(u) {
        eval(parse(text = getURL(u, followlocation = TRUE,
                   cainfo = system.file("CurlSSL", "cacert.pem",
                   package = "RCurl"))), envir = .GlobalEnv)
    })
}

## URL for the raw file
url = "https://raw.github.com/MatthieuStigler/tsDyn/predRoll/tsDyn/R/"

## Functions within the folder


funcs = c("BBCtest.R", "KapShinTest.R", "TVAR.sim.R",
          "TVAR_LRtest.R", "TVARestim.R", "TVECM.HSTest.R", "TVECM.R",
          "TVECM.sim.R", "TVECMSeoTest.R", "VECM_symbolic.R", "aar.R",
          "accuracy.R", "acf.custom.R", "autopairs.R", "autotriples.R",
          "autotriples.rgl.R", "delta.R", "far.R", "guiSystem.R",
          "isLinear.R", "lineVar.R", "linear.R", "llar.R", "lstar.R",
          "misc.R", "miscSETAR.R", "nlVar-methods.R", "nlar-methods.R",
          "nnet.R", "predict.R", "predict_rolling.R", "rank.select.R",
          "rank.test.R", "resVar.R", "selectLSTAR.R", "selectSETAR.R",
          "setar.R", "setar.sim.R", "setarTest.R", "star.R", "tgarch.R",
          "tsDynGUI.R", "vec2var.tsDyn.R", "zzz.R")

## Load all the code
source_github(url = paste0(url, funcs))



########################################################################
## Title: Example of state-space model
## Date: 2012-12-21
########################################################################

## Load library dlm and load the manual.
library(dlm)
vignette("dlm")

## Load the UK gas consumption data
data(UKgas)

## Plot the data
plot(UKgas)

## Log transform the data to try to alleviate the increase in variance.
logUKgas = log(UKgas)

## Still show sign of heteroscedasticity.
plot(logUKgas)

## NOTE: Using differencing and power transformation to satisfy the
## stationarity assumption is generally not a good approach.

## Use non-stationary time series model such as state-space
## model. (seaonal smoothing window is set to 4 because it is
## quarterly data)
UKgas.stl = stl(UKgas, s.window = 4)
plot(UKgas.stl)

## NOTE: The model does not assume stationarity and is capable of
## modelling non-stationary time series. In addition, it handles
## missing value naturally.




########################################################################
## Title: Random forests and conditional trees
## Date: 2013-01-02
########################################################################

library(randomForest)
library(party)
library(VIM)

## Load the titanic data
url = "http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/titanic3.csv"
titanic.df = read.csv(url, header = TRUE, stringsAsFactors = FALSE)
titanic.df$sex = as.factor(titanic.df$sex)


## Does not show much, because these are bi-lateral relationships.
plot(titanic.df[, c("pclass", "survived", "age", "sibsp", "parch", "fare")])


## Impute missing value with hot-deck imputation.
n.miss = sum(is.na(titanic.df$age))
set.seed(587)
impTitanic.df = hotdeck(data = titanic.df, variable = "age",
    ord_var = c("parch", "sibsp"), domain_var = "pclass")
impTitanic.df = hotdeck(data = impTitanic.df, variable = "fare",
    domain_var = "pclass")


## obtain training and validation set
n = NROW(impTitanic.df)
ind = sample(1:n, n * 0.6, replace = FALSE)
titTrain.df = impTitanic.df[ind, ]
titValid.df = impTitanic.df[-ind, ]

## Glm
tit.glm = glm(survived ~ pclass + sex:age + sibsp + parch + fare,
              family = binomial(logit), data = titTrain.df)
titValid.df$glmpred = as.numeric(predict(tit.glm, titValid.df,
    type = "response") >= 0.5)
with(titValid.df, table(survived, glmpred))
sum(titValid.df$survive != titValid.df$glmpred)/NROW(titValid.df)


## Random Forest
tit.rf = randomForest(formula = as.factor(survived) ~ pclass + sex + age +
                      sibsp + parch + fare, data = titTrain.df, ntree = 5000,
                      importance = TRUE, na.action = na.omit)


titValid.df$rfpred = predict(tit.rf, titValid.df)
with(titValid.df, table(survived, rfpred))
sum(titValid.df$survive != titValid.df$rfpred)/NROW(titValid.df)

## Conditional tree
tit.ct = ctree(survived ~ pclass + sex + age + sibsp + parch + fare,
               data = titTrain.df)
## plot(tit.ct)

titValid.df$ctpred = as.numeric(predict(tit.ct, titValid.df) >= 0.5)
with(titValid.df, table(survived, ctpred))
sum(titValid.df$survive != titValid.df$ctpred)/NROW(titValid.df)





## Lets bootstrap the miss-classification rate
n.iter = 100
n = NROW(impTitanic.df)
glmError = double(n.iter)
rfError = double(n.iter)
ctError = double(n.iter)

for(i in 1:n.iter){
    print(i)
    ind = sample(1:n, n * 0.6, replace = FALSE)
    titTrain.df = impTitanic.df[ind, ]
    titValid.df = impTitanic.df[-ind, ]

    ## Glm
    tit.glm = glm(survived ~ pclass + sex:age + sibsp + parch + fare,
    family = binomial(logit), data = titTrain.df)
    titValid.df$glmpred = as.numeric(predict(tit.glm, titValid.df,
    type = "response") >= 0.5)
    glmError[i] = sum(titValid.df$survive != titValid.df$glmpred)/
        NROW(titValid.df)


    ## Random Forest
    tit.rf = randomForest(formula = as.factor(survived) ~ pclass + sex + age +
    sibsp + parch + fare, data = titTrain.df, ntree = 5000,
    importance = TRUE, na.action = na.omit)
    titValid.df$rfpred = predict(tit.rf, titValid.df)
    rfError[i] = sum(titValid.df$survive != titValid.df$rfpred)/NROW(titValid.df)

    ## Conditional tree
    tit.ct = ctree(survived ~ pclass + sex + age + sibsp + parch + fare,
    data = titTrain.df)
    titValid.df$ctpred = as.numeric(predict(tit.ct, titValid.df) >= 0.5)
    ctError[i] = sum(titValid.df$survive != titValid.df$ctpred)/NROW(titValid.df)
}

allError.df = data.frame(glmError = glmError, rfError = rfError,
    ctError = ctError)
apply(allError.df, 2, mean)


########################################################################
## Title: Test the Berkeley earth package
## Date: 2013-01-07
########################################################################

library(BerkeleyEarth)

test = downloadBerkeley(BerkeleyUrls[1, , drop = FALSE])


########################################################################
## Title: Intro to Rcpp
## Date: 2013-01-10
########################################################################

library(Rcpp)

cppFunction('
  int add(int x, int y, int z) {
    int sum = x + y + z;
    return sum;
  }'
)

add # like a regular R function, printing displays info about the function
add(1, 2, 3)



########################################################################
## Title: Social Network Analysis Labs in R and SoNIA (Stanford)
##        (http://sna.stanford.edu/rlabs.php)
## Date: 2012-01-25
########################################################################

##########################################################################
# You may cite these labs as follows: McFarland, Daniel, Solomon Messing,
# Mike Nowak, and Sean Westwood. 2010. "Social Network Analysis
# Labs in R." Stanford University.
##########################################################################


##########################################################################
# LAB 1 - Introductory Lab
# The point of this lab is to introduce students to the packages of
# SNA and Igraph, to cover some basic R commands, to load and manage
# data, to generate graph visualizations, and to export the data for
# use elsewhere.
##########################################################################

###
# 0. R BASICS
###


# To install all packages need for Social Network Analysis
# Labs in R, uncomment and run the following code:

#source("http://sna.stanford.edu/setup.R")

# You only need to run this once per computer!

# To load the packages from , you need to call the "library"
# command. Note that you need to do this each session; packages
# don't load automatically by default (though you can set this
# as a preference if you'd like).

# For this lab, we will use the "igraph" package.
# A manual is available at
# http://cran.r-project.org/web/packages/igraph/igraph.pdf.
library(igraph)

# Sometimes, different packages overlap in functionality and
# cause unexpected behavior when both are loaded simultaneously.
# If you ever want to remove an existing library, use the
# "detach" command:
#
# detach(package:igraph)

# IMPORTANT NOTE: Unlike in most languages, R objects are numbered
# from 1 instead of 0, so if you want the first element in a
# vector, you would reference it by vector_name[1]. HOWEVER,
# igraph objects are numbered starting from 0. This can lead to
# lots of confusion, since it's not always obvious at first which
# objects are native to R and which belong to igraph.


###
# 1. LOADING DATA
###

# The <- operator sets a variable equal to something. In this case,
# we will set a number of basic R data structures, called "data
# frames," to hold the contents of the files we will open.
#
# read.table() is the most common R command for loading data from
# files in which values are in tabular format. The function loads
# the table into a data frame object, which is the basic data type
# for most operations in R. By default, R assumes that the table
# has no header and is delimited by any white space; these
# settings are fine for our purposes here.
#
# One handy aspect of R is that you can read in data from a URL
# directly by referencing the URL in the read.table() function,
# as follows:
advice_data_frame <- read.table('http://sna.stanford.edu/sna_R_labs/data/Krack-High-Tec-edgelist-Advice.txt')
friendship_data_frame <- read.table('http://sna.stanford.edu/sna_R_labs/data/Krack-High-Tec-edgelist-Friendship.txt')
reports_to_data_frame <- read.table('http://sna.stanford.edu/sna_R_labs/data/Krack-High-Tec-edgelist-ReportsTo.txt')

# If the files you want to work with are on your local machine,
# the easiest way to access them is to first set your working
# directory via the setwd() command, and then reference the
# files by name:
#
# setwd('path/to/your_directory')
# your_data_frame <- read.table('your_file_name')

# Note that when you set a variable equal to something, if all
# goes well R will not provide any feedback. To see the data we
# just loaded, it's necessary to call the variables directly.
advice_data_frame

# Since this is a bit long, we can see just the top six rows via
# head()...
head(friendship_data_frame)

# ... or the bottom six rows via tail().
tail(reports_to_data_frame)

# To view your data in a spreadsheet-like window, use the command 'fix()'.
fix(reports_to_data_frame)

# The attribute data for this lab is in a comma-separated-value
# (CSV) file. read.csv() loads a CSV file into a data frame
# object. In this case, we do have a header row, so we set
# header=T, which tells R that the first row of data contains
# column names.
attributes <- read.csv('http://sna.stanford.edu/sna_R_labs/data/Krack-High-Tec-Attributes.csv', header=T)
attributes

# Other commands may be used to load data from files in different
# formats. read.delim() is a general function for loading any
# delimited text file. The default is tab-delimited, but this can
# be overridden by setting the "sep" parameter. For example:
#
#     f <- read.delim("tab_delimited_file.txt")
#     f <- read.delim("colon_delimited_file.txt", sep=':')
#
# The 'foreign' package will allow you to read a few other
# custom data types, such as SPSS files via read.spss() and
# STATA files via read.dta().

# When data files are part of an R package you can read them as
# follows:
#
# data(kracknets, package = "NetData")

# In the future, we will load data this way. However, it is useful
# to get a sense of how things often must be done in R.


###
# 2.2. LOADING GRAPHS
###

# For convenience, we can assign column names to our newly
# imported data frames. c() is a common generic R function that
# combines its arguments into a single vector.
colnames(advice_data_frame) <- c('ego', 'alter', 'advice_tie')
head(advice_data_frame)

colnames(friendship_data_frame) <- c('ego', 'alter', 'friendship_tie')
head(friendship_data_frame)

colnames(reports_to_data_frame) <- c('ego', 'alter', 'reports_to_tie')
head(reports_to_data_frame)

# Take a look at each data frame using the 'fix()" function. Note that you'll
# need to close each fix window before R will evaluate the next line of code.
fix(advice_data_frame)
fix(friendship_data_frame)
fix(reports_to_data_frame)

# Before we merge these data, we need to make sure 'ego' and 'alter' are the
# same across data sets. We can compare each row using the == syntax.
# The command below should return TRUE for every row if all ego rows
# are the same for advice and friendship:
advice_data_frame$ego == friendship_data_frame$ego

# That's a lot of output to sort through. Instead, we can just have R return
# which row entries are not equal using the syntax below:
which(advice_data_frame$ego != friendship_data_frame$ego)

# Repeat for other variables
which(advice_data_frame$alter != friendship_data_frame$alter)
which(reports_to_data_frame$alter != friendship_data_frame$alter)
which(reports_to_data_frame$ego != friendship_data_frame$ego)

# Now that we've verified they are all the same, we can combine them into
# a single data frame.
krack_full_data_frame <- cbind(advice_data_frame,
    friendship_data_frame$friendship_tie,
    reports_to_data_frame$reports_to_tie)
head(krack_full_data_frame)

# Notice that the last two variable names are now
# "reports_to_data_frame$reports_to_tie"
# and "friendship_data_frame$friendship_tie".
# That's a little long. We can rename them
# as follows:

names(krack_full_data_frame)[4:5] <- c("friendship_tie",
    "reports_to_tie")
head(krack_full_data_frame)

# Another way to build the data frame is to use R's
# data.frame syntax from the start:
krack_full_data_frame <- data.frame(ego = advice_data_frame[,1],
    alter = advice_data_frame[,2],
    advice_tie = advice_data_frame[,3],
    friendship_tie = friendship_data_frame[,3],
    reports_to_tie = reports_to_data_frame[,3])
head(krack_full_data_frame)


# Now let's move on to some data processing.

# Reduce to non-zero edges so that the edge list only contains
# actual ties of some type.
krack_full_nonzero_edges <- subset(krack_full_data_frame,
    (advice_tie > 0 | friendship_tie > 0 | reports_to_tie > 0))
head(krack_full_nonzero_edges)

# Now we can import our data into a "graph" object using igraph's
# graph.data.frame() function. Coercing the data into a graph
# object is what allows us to perform network-analysis techniques.
krack_full <- graph.data.frame(krack_full_nonzero_edges)
summary(krack_full)

# By default, graph.data.frame() treats the first two columns of
# a data frame as an edge list and any remaining columns as
# edge attributes. Thus, the 232 edges appearing in the summary()
# output refer to the 232 pairs of vertices that are joined by
# *any type* of tie. The tie types themselves are listed as edge
# attributes.

# To get a vector of edges for a specific type of tie, use the
# get.edge.attribute() function.
get.edge.attribute(krack_full, 'advice_tie')
get.edge.attribute(krack_full, 'friendship_tie')
get.edge.attribute(krack_full, 'reports_to_tie')

# If you would like to symmetrize the network, making all
# asymmetric ties symmetric, use the as.undirected()
# function:
krack_full_symmetrized <- as.undirected(krack_full, mode='collapse')
summary(krack_full_symmetrized)



###
# 3. ADDING VERTEX ATTRIBUTES TO A GRAPH OBJECT
###

# One way to add the attributes to your graph object is to iterate
# through each attribute and each vertex. This means that we will
# add one attribute at a time to each vertex in the network.
#
# V(krack_full) returns a list of the IDs of each vertex in the
# graph. names(attributes) returns a list of the column names in
# the attributes table. The double-for loop tells R to repeat the
# code between the brackets once for each attribute and once for
# each vertex.
for (i in V(krack_full)) {
    for (j in names(attributes)) {
        krack_full <- set.vertex.attribute(krack_full,
                                           j,
                                           index = i,
                                           attributes[i + 1, j])
    }
}

# A shorter way is to just read in attribute names when you
# create the graph object:

# First create a vector of vertex labels, in this case 1:n
attributes = cbind(1:length(attributes[,1]), attributes)

krack_full <- graph.data.frame(d = krack_full_nonzero_edges,
                               vertices = attributes)

# Note that we now have 'AGE,' 'TENURE,' 'LEVEL,' and 'DEPT'
# listed alongside 'name' as vertex attributes.
summary(krack_full)

# We can see a list of the values for a given attribute for all of
# the actors in the network.
get.vertex.attribute(krack_full, 'AGE')
get.vertex.attribute(krack_full, 'TENURE')
get.vertex.attribute(krack_full, 'LEVEL')
get.vertex.attribute(krack_full, 'DEPT')


###
# 4. VISUALIZE THE NETWORKS
###

# We can use R's general-purpose plot() method to generate custom
# visualizations of the network.

# R only lets us look at one plot at a time.  To make our work easier
# we will save our plots as PDF files.  To jus create a plot execute
# the code between the PDF function and "dev.off()".

# In order to save PDF files we must tell R where to put them.  We do
# this with the setwd() command.  You must put the full path to the
# folder where you will output the files here.

# In OS X you can get this information by selecting the folder, right
# clicking and selecting "Get Info."  The path is listed under "Where."

# In Windows you can get this information by selecting the folder, right
# clicking and selecting "Properties."  The path information is listed
# "location".

# example: setwd("/Users/seanwestwood/Desktop/lab_1")
setwd("")

# First, let's plot the network with all possible ties.
pdf("1.1_Krackhardt_Full.pdf")
plot(krack_full)
dev.off()

# This is a bit of a jumble, so let's look at the networks for
# single edge types.

# advice only
krack_advice_only <- delete.edges(krack_full,
    E(krack_full)[get.edge.attribute(krack_full,
    name = "advice_tie") == 0])
summary(krack_advice_only)
pdf("1.2_Krackhardt_Advice.pdf")
plot(krack_advice_only)
dev.off()

# friendship only
krack_friendship_only <- delete.edges(krack_full,
    E(krack_full)[get.edge.attribute(krack_full,
    name = "friendship_tie") == 0])
summary(krack_friendship_only)
pdf("1.3_Krackhardt_Friendship.pdf")
plot(krack_friendship_only)
dev.off()

# reports-to only
krack_reports_to_only <- delete.edges(krack_full,
    E(krack_full)[get.edge.attribute(krack_full,
    name = "reports_to_tie") == 0])
summary(krack_reports_to_only)
pdf("1.4_Krackhardt_Reports.pdf")
plot(krack_reports_to_only)
dev.off()

# Still kind of messy, so let's clean things up a bit. For
# simplicity, we'll focus on reports_to ties for now.

# First, we can optimize the layout by applying the layout
# algorithm to the specific set of ties we care about. Here
# we'll use Fruchterman-Rheingold; other options are
# described in the igraph help page for "layout," which
# can be accessed by entering ?layout.

reports_to_layout <- layout.fruchterman.reingold(krack_reports_to_only)
pdf("1.5_Krackhardt_Reports_Fruchterman_Reingold.pdf")
plot(krack_reports_to_only,
     layout=reports_to_layout)
dev.off()

# Now let's color-code vertices by department and clean up the
# plot by removing vertex labels and shrinking the arrow size.
dept_vertex_colors = get.vertex.attribute(krack_full,"DEPT")
colors = c('Black', 'Red', 'Blue', 'Yellow', 'Green')
dept_vertex_colors[dept_vertex_colors == 0] = colors[1]
dept_vertex_colors[dept_vertex_colors == 1] = colors[2]
dept_vertex_colors[dept_vertex_colors == 2] = colors[3]
dept_vertex_colors[dept_vertex_colors == 3] = colors[4]
dept_vertex_colors[dept_vertex_colors == 4] = colors[5]

pdf("1.6_Krackhardt_Reports_Color.pdf")
plot(krack_reports_to_only,
    layout=reports_to_layout,
    vertex.color=dept_vertex_colors,
    vertex.label=NA,
    edge.arrow.size=.5)
dev.off()
# Now let's set the vertex size by tenure.
tenure_vertex_sizes = get.vertex.attribute(krack_full,"TENURE")

pdf("1.7_Krackhardt_Reports_Vertex_Size.pdf")
plot(krack_reports_to_only,
     layout=reports_to_layout,
     vertex.color=dept_vertex_colors,
     vertex.label=NA,
     edge.arrow.size=.5,
     vertex.size=tenure_vertex_sizes)
dev.off()

# Now let's incorporate additional tie types. We'll use the
# layout generated by the reports-to ties but overlay the
# advice and friendship ties in red and blue.

tie_type_colors = c(rgb(1,0,0,.5), rgb(0,0,1,.5), rgb(0,0,0,.5))
E(krack_full)$color[ E(krack_full)$advice_tie==1 ] = tie_type_colors[1]
E(krack_full)$color[ E(krack_full)$friendship_tie==1 ] = tie_type_colors[2]
E(krack_full)$color[ E(krack_full)$reports_to_tie==1 ] = tie_type_colors[3]
E(krack_full)$arrow.size=.5
V(krack_full)$color = dept_vertex_colors
V(krack_full)$frame = dept_vertex_colors

pdf("1.8_Krackhardt_Overlayed_Ties.pdf")
plot(krack_full,
     layout=reports_to_layout,
     vertex.color=dept_vertex_colors,
     vertex.label=NA,
     edge.arrow.size=.5,
     vertex.size=tenure_vertex_sizes)


# Add a legend. Note that the plot window must be open for this to
# work.
legend(1,
       1.25,
       legend = c('Advice',
                  'Friendship',
                  'Reports To'),
       col = tie_type_colors,
       lty=1,
       cex = .7)
dev.off()

# Another option for visualizing different network ties relative
# to one another is to overlay the edges from one tie type on the
# structure generated by another tie type. Here we can use the
# reports-to layout but show the friendship ties:

pdf("1.9_Krackhardt_Overlayed_Structure.pdf")
plot(krack_friendship_only,
     layout=reports_to_layout,
     vertex.color=dept_vertex_colors,
     vertex.label=NA,
     edge.arrow.size=.5,
     vertex.size=tenure_vertex_sizes,
     main='Krackhardt High-Tech Managers')
dev.off()


###
# 5. EXPORT THE NETWORK
###

# The write.graph() function exports a graph object in various
# formats readable by other programs. There is no explicit
# option for a UCINET data type, but you can export the graph
# as a Pajek object by setting the 'format' parameter to 'pajek.'
# Note that the file will appear in whichever directory is set
# as the default in R's preferences, unless you previously
# changed this via setwd().
write.graph(krack_full, file='krack_full.dl', format="pajek")

# For a more general file type (e.g., importable to Excel),
# use the "edgelist" format. Note that neither of these will
# write the attributes; only the ties are maintained.
write.graph(krack_full, file='krack_full.txt', format="edgelist")




####################################################################
# LAB 2: Methodological beginnings - Density, Reciprocity, Triads, #
# Transitivity, and heterogeneity. Node and network statistics.    #
####################################################################


# NOTE: if you have trouble because some packages are not installed,
# see lab 1 for instructions on how to install all necessary packages.
# Also see Lab 1 for prior functions.


##############################################################
#
# Lab 2
#
# The purpose of this lab is to acquire basic cohesion
# metrics of density, reciprocity, reach, path distance,
# and transitivity. In addition, we'll develop triadic
# analyses and a measure of ego-network heterogenity.
#
##############################################################



###
# 1. SET UP SESSION
###
install.packages("NetData")

library(igraph)
library(NetData)


###
# 2. LOAD DATA
###

# We would ordinarily need to follow the same proceedure we did for the Krackhardt data
# as we did in lab 1; see that lab for detail.

data(kracknets, package = "NetData")

# Reduce to non-zero edges and build a graph object
krack_full_nonzero_edges <- subset(krack_full_data_frame, (advice_tie > 0 | friendship_tie > 0 | reports_to_tie > 0))
head(krack_full_nonzero_edges)

krack_full <- graph.data.frame(krack_full_nonzero_edges)
summary(krack_full)

# Set vertex attributes
for (i in V(krack_full)) {
    for (j in names(attributes)) {
        krack_full <- set.vertex.attribute(krack_full, j, index=i, attributes[i+1,j])
    }
}
summary(krack_full)

# Create sub-graphs based on edge attributes
krack_advice <- delete.edges(krack_full, E(krack_full)[get.edge.attribute(krack_full,name = "advice_tie")==0])
summary(krack_advice)

krack_friendship <- delete.edges(krack_full, E(krack_full)[get.edge.attribute(krack_full,name = "friendship_tie")==0])
summary(krack_friendship)

krack_reports_to <- delete.edges(krack_full, E(krack_full)[get.edge.attribute(krack_full,name = "reports_to_tie")==0])
summary(krack_reports_to)


###
# 3. NODE-LEVEL STATISTICS
###

# Compute the indegree and outdegree for each node, first in the
# full graph (accounting for all tie types) and then in each
# tie-specific sub-graph.
deg_full_in <- degree(krack_full, mode="in")
deg_full_out <- degree(krack_full, mode="out")
deg_full_in
deg_full_out

deg_advice_in <- degree(krack_advice, mode="in")
deg_advice_out <- degree(krack_advice, mode="out")
deg_advice_in
deg_advice_out

deg_friendship_in <- degree(krack_friendship, mode="in")
deg_friendship_out <- degree(krack_friendship, mode="out")
deg_friendship_in
deg_friendship_out

deg_reports_to_in <- degree(krack_reports_to, mode="in")
deg_reports_to_out <- degree(krack_reports_to, mode="out")
deg_reports_to_in
deg_reports_to_out

# Reachability can only be computed on one vertex at a time. To
# get graph-wide statistics, change the value of "vertex"
# manually or write a for loop. (Remember that, unlike R objects,
# igraph objects are numbered from 0.)

reachability <- function(g, m) {
    reach_mat = matrix(nrow = vcount(g),
                       ncol = vcount(g))
    for (i in 1:vcount(g)) {
        reach_mat[i,] = 0
        this_node_reach <- subcomponent(g, (i - 1), mode = m)

        for (j in 1:(length(this_node_reach))) {
            alter = this_node_reach[j] + 1
            reach_mat[i, alter] = 1
        }
    }
    return(reach_mat)
}

reach_full_in <- reachability(krack_full, 'in')
reach_full_out <- reachability(krack_full, 'out')
reach_full_in
reach_full_out

reach_advice_in <- reachability(krack_advice, 'in')
reach_advice_out <- reachability(krack_advice, 'out')
reach_advice_in
reach_advice_out

reach_friendship_in <- reachability(krack_friendship, 'in')
reach_friendship_out <- reachability(krack_friendship, 'out')
reach_friendship_in
reach_friendship_out

reach_reports_to_in <- reachability(krack_reports_to, 'in')
reach_reports_to_out <- reachability(krack_reports_to, 'out')
reach_reports_to_in
reach_reports_to_out


# Often we want to know path distances between individuals in a network.
# This is often done by calculating geodesics, or shortest paths between
# each ij pair. One can symmetrize the data to do this (see lab 1), or
# calculate it for outward and inward ties separately. Averaging geodesics
# for the entire network provides an average distance or sort of cohesiveness
# score. Dichotomizing distances reveals reach, and an average of reach for
# a network reveals what percent of a network is connected in some way.

# Compute shortest paths between each pair of nodes.
sp_full_in <- shortest.paths(krack_full, mode='in')
sp_full_out <- shortest.paths(krack_full, mode='out')
sp_full_in
sp_full_out

sp_advice_in <- shortest.paths(krack_advice, mode='in')
sp_advice_out <- shortest.paths(krack_advice, mode='out')
sp_advice_in
sp_advice_out

sp_friendship_in <- shortest.paths(krack_friendship, mode='in')
sp_friendship_out <- shortest.paths(krack_friendship, mode='out')
sp_friendship_in
sp_friendship_out

sp_reports_to_in <- shortest.paths(krack_reports_to, mode='in')
sp_reports_to_out <- shortest.paths(krack_reports_to, mode='out')
sp_reports_to_in
sp_reports_to_out


# Assemble node-level stats into single data frame for export as CSV.

# First, we have to compute average values by node for reachability and
# shortest path. (We don't have to do this for degree because it is
# already expressed as a node-level value.)
reach_full_in_vec <- vector()
reach_full_out_vec <- vector()
reach_advice_in_vec <- vector()
reach_advice_out_vec <- vector()
reach_friendship_in_vec <- vector()
reach_friendship_out_vec <- vector()
reach_reports_to_in_vec <- vector()
reach_reports_to_out_vec <- vector()

sp_full_in_vec <- vector()
sp_full_out_vec <- vector()
sp_advice_in_vec <- vector()
sp_advice_out_vec <- vector()
sp_friendship_in_vec <- vector()
sp_friendship_out_vec <- vector()
sp_reports_to_in_vec <- vector()
sp_reports_to_out_vec <- vector()

for (i in 1:vcount(krack_full)) {
    reach_full_in_vec[i] <- mean(reach_full_in[i,])
    reach_full_out_vec[i] <- mean(reach_full_out[i,])
    reach_advice_in_vec[i] <- mean(reach_advice_in[i,])
    reach_advice_out_vec[i] <- mean(reach_advice_out[i,])
    reach_friendship_in_vec[i] <- mean(reach_friendship_in[i,])
    reach_friendship_out_vec[i] <- mean(reach_friendship_out[i,])
    reach_reports_to_in_vec[i] <- mean(reach_reports_to_in[i,])
    reach_reports_to_out_vec[i] <- mean(reach_reports_to_out[i,])

    sp_full_in_vec[i] <- mean(sp_full_in[i,])
    sp_full_out_vec[i] <- mean(sp_full_out[i,])
    sp_advice_in_vec[i] <- mean(sp_advice_in[i,])
    sp_advice_out_vec[i] <- mean(sp_advice_out[i,])
    sp_friendship_in_vec[i] <- mean(sp_friendship_in[i,])
    sp_friendship_out_vec[i] <- mean(sp_friendship_out[i,])
    sp_reports_to_in_vec[i] <- mean(sp_reports_to_in[i,])
    sp_reports_to_out_vec[i] <- mean(sp_reports_to_out[i,])
}

# Next, we assemble all of the vectors of node-levelvalues into a
# single data frame, which we can export as a CSV to our working
# directory.
node_stats_df <- cbind(deg_full_in,
                       deg_full_out,
                       deg_advice_in,
                       deg_advice_out,
                       deg_friendship_in,
                       deg_friendship_out,
                       deg_reports_to_in,
                       deg_reports_to_out,

                       reach_full_in_vec,
                       reach_full_out_vec,
                       reach_advice_in_vec,
                       reach_advice_out_vec,
                       reach_friendship_in_vec,
                       reach_friendship_out_vec,
                       reach_reports_to_in_vec,
                       reach_reports_to_out_vec,

                       sp_full_in_vec,
                       sp_full_out_vec,
                       sp_advice_in_vec,
                       sp_advice_out_vec,
                       sp_friendship_in_vec,
                       sp_friendship_out_vec,
                       sp_reports_to_in_vec,
                       sp_reports_to_out_vec)

write.csv(node_stats_df, 'krack_node_stats.csv')

# Question #1 - What do these statistics tell us about
# each network and its individuals in general?

###
# 3. NETWORK-LEVEL STATISTICS
###

# Many initial analyses of networks begin with distances and reach,
# and then move towards global summary statistics of the network.
#
# As a reminder, entering a question mark followed by a function
# name (e.g., ?graph.density) pulls up the help file for that function.
# This can be helpful to understand how, exactly, stats are calculated.

# Degree
mean(deg_full_in)
sd(deg_full_in)
mean(deg_full_out)
sd(deg_full_out)

mean(deg_advice_in)
sd(deg_advice_in)
mean(deg_advice_out)
sd(deg_advice_out)

mean(deg_friendship_in)
sd(deg_friendship_in)
mean(deg_friendship_out)
sd(deg_friendship_out)

mean(deg_reports_to_in)
sd(deg_reports_to_in)
mean(deg_reports_to_out)
sd(deg_reports_to_out)


# Shortest paths
# ***Why do in and out come up with the same results?
# In and out shortest paths are simply transposes of one another;
# thus, when we compute statistics across the whole network they have to be the same.

mean(sp_full_in[which(sp_full_in != Inf)])
sd(sp_full_in[which(sp_full_in != Inf)])
mean(sp_full_out[which(sp_full_out != Inf)])
sd(sp_full_out[which(sp_full_out != Inf)])

mean(sp_advice_in[which(sp_advice_in != Inf)])
sd(sp_advice_in[which(sp_advice_in != Inf)])
mean(sp_advice_out[which(sp_advice_out != Inf)])
sd(sp_advice_out[which(sp_advice_out != Inf)])

mean(sp_friendship_in[which(sp_friendship_in != Inf)])
sd(sp_friendship_in[which(sp_friendship_in != Inf)])
mean(sp_friendship_out[which(sp_friendship_out != Inf)])
sd(sp_friendship_out[which(sp_friendship_out != Inf)])

mean(sp_reports_to_in[which(sp_reports_to_in != Inf)])
sd(sp_reports_to_in[which(sp_reports_to_in != Inf)])
mean(sp_reports_to_out[which(sp_reports_to_out != Inf)])
sd(sp_reports_to_out[which(sp_reports_to_out != Inf)])

# Reachability
mean(reach_full_in[which(reach_full_in != Inf)])
sd(reach_full_in[which(reach_full_in != Inf)])
mean(reach_full_out[which(reach_full_out != Inf)])
sd(reach_full_out[which(reach_full_out != Inf)])

mean(reach_advice_in[which(reach_advice_in != Inf)])
sd(reach_advice_in[which(reach_advice_in != Inf)])
mean(reach_advice_out[which(reach_advice_out != Inf)])
sd(reach_advice_out[which(reach_advice_out != Inf)])

mean(reach_friendship_in[which(reach_friendship_in != Inf)])
sd(reach_friendship_in[which(reach_friendship_in != Inf)])
mean(reach_friendship_out[which(reach_friendship_out != Inf)])
sd(reach_friendship_out[which(reach_friendship_out != Inf)])

mean(reach_reports_to_in[which(reach_reports_to_in != Inf)])
sd(reach_reports_to_in[which(reach_reports_to_in != Inf)])
mean(reach_reports_to_out[which(reach_reports_to_out != Inf)])
sd(reach_reports_to_out[which(reach_reports_to_out != Inf)])

# Density
graph.density(krack_full)
graph.density(krack_advice)
graph.density(krack_friendship)
graph.density(krack_reports_to)

# Reciprocity
reciprocity(krack_full)
reciprocity(krack_advice)
reciprocity(krack_friendship)
reciprocity(krack_reports_to)

# Transitivity
transitivity(krack_full)
transitivity(krack_advice)
transitivity(krack_friendship)
transitivity(krack_reports_to)

# Triad census. Here we'll first build a vector of labels for
# the different triad types. Then we'll combine this vector
# with the triad censuses for the different networks, which
# we'll export as a CSV.

census_labels = c('003',
                  '012',
                  '102',
                  '021D',
                  '021U',
                  '021C',
                  '111D',
                  '111U',
                  '030T',
                  '030C',
                  '201',
                  '120D',
                  '120U',
                  '120C',
                  '210',
                  '300')
tc_full <- triad.census(krack_full)
tc_advice <- triad.census(krack_advice)
tc_friendship <- triad.census(krack_friendship)
tc_reports_to <- triad.census(krack_reports_to)

triad_df <- data.frame(census_labels,
                       tc_full,
                       tc_advice,
                       tc_friendship,
                       tc_reports_to)
triad_df

# To export any of these vectors to a CSV for use in another program, simply
# use the write.csv() command:
write.csv(triad_df, 'krack_triads.csv')

# Question #2 - (a) How do the three networks differ on network statictics?
# (b) What does the triad census tell us? Can you calculate the likelihood of
# any triad's occurrence? (c) See the back of Wasserman and Faust and its section
# on triads. Calculate the degree of clustering and hierarchy in Excel.
# What do we learn from that?




###
# 4. HETEROGENEITY
###

# Miller and McPherson write about processes of homophily and
# here we take a brief look at one version of this issue.
# In particular, we look at the extent to which each actor's
# "associates" (friend, advisor, boos) are heterogenous or not.

# We'll use a statistic called the IQV, or Index of Qualitative
# Variation. This is just an implementation of Blau's Index of
# Heterogeneity (known to economists as the Herfindahl-Hirschman
# index), normalized so that perfect heterogeneity (i.e., equal
# distribution across categories) equals 1.

# NOTE that this code only works with categorical variables that
# have been numerically coded to integer values that ascend
# sequentially from 0; you may have to recode your data to get this
# to work properly.
# We are interested in many of the attributes of nodes.  To save
# time and to make our lives better we are going to create a function
# that will provide an IQV statistic for any network and for
# any categorical variable.  A function is a simple way to
# create code that is both reusable and easier to edit.

# Functions have names and receive arguments.  For example,
# anytime you call table() you are calling the table function.
# We could write code to duplicate the table function for each
# of our variables, but it is faster to write a single general tool
# that will provide frequencies for any variable. If I have
# a dataframe with the variable gender and I want to see the
# split of males and females I would pass the argument
# "dataframe$gender" to the table function. We follow a
# similar model here. Understanding each step is less important
# than understanding the usefulness and power of functions.

get_iqvs <- function(graph, attribute) {

#we have now defined a function, get_iqvs, that will take the
# graph "graph" and find the iqv statistic for the categorical
# variable "attribute." Within this function whenever we use the
#variables graph or attribute they correspond to the graph and
# variable we passed (provided) to the function

mat <- get.adjacency(graph)

# To make this function work on a wide variety of variables we
# find out how many coded levels (unique responses) exist for
# the attribute variable programatically

    attr_levels = get.vertex.attribute(graph,
                                       attribute,
                                       V(graph))

    num_levels = length(unique(attr_levels))
    iqvs = rep(0, nrow(mat))

# Now that we know how many levels exist we want to loop
# (go through) each actor in the network. Loops iterate through
# each value in a range.  Here we are looking through each ego
# in the range of egos starting at the first and ending at the
# last.  The function nrow provides the number of rows in an
# object and the ":" opperand specifies the range.  Between
# the curly braces of the for loop ego will represent exactly
# one value between 1 and the number of rows in the graph
# object, iterating by one during each execution of the loop.

    for (ego in 1:nrow(mat)) {

        # initialize actor-specific variables
        alter_attr_counts = rep(0, num_levels)
        num_alters_this_ego = 0
        sq_fraction_sum = 0

# For each ego we want to check each tied alter for the same
# level on the variable attribute as the ego.

        for (alter in 1:ncol(mat)) {

            # only examine alters that are actually tied to ego
            if (mat[ego, alter] == 1) {

                num_alters_this_ego = num_alters_this_ego + 1

                # get the alter's level on the attribute
                alter_attr = get.vertex.attribute(graph,
                    attribute, (alter - 1))

                # increment the count of alters with this level
                # of the attribute by 1
                alter_attr_counts[alter_attr + 1] =
                    alter_attr_counts[alter_attr + 1] + 1
            }
        }

        # now that we're done looping through all of the alters,
        # get the squared fraction for each level of the attribute
        # out of the total number of attributes
        for (i in 1:num_levels) {
            attr_fraction = alter_attr_counts[i] /
                num_alters_this_ego
            sq_fraction_sum = sq_fraction_sum + attr_fraction ^ 2
        }

        # now we can compute the ego's blau index...
        blau_index = 1 - sq_fraction_sum

        # and the ego's IQV, which is just a normalized blau index
        iqvs[ego] = blau_index / (1 - (1 / num_levels))
    }

# The final part of a function returns the calculated value.
#  So if we called get_iqvs(testgraph, gender) return would
# provide the iqvs for gender in the test graph.  If we are also
# intersted in race we could simply change the function call
# to get_iqvs(testgraph, race).  No need to write all this
# code again for different variables.

    return(iqvs)
}



# For this data set, we'll look at homophily across departments,
# which is already coded 0-4, so no recoding is needed.

advice_iqvs <- get_iqvs(krack_advice, 'DEPT')
advice_iqvs

friendship_iqvs <- get_iqvs(krack_friendship, 'DEPT')
friendship_iqvs

reports_to_iqvs <- get_iqvs(krack_reports_to, 'DEPT')
reports_to_iqvs

# Question #3 - What does the herfindahl index reveal about
# attribute sorting in networks? What does it mean for each network?


#####
# Extra-credit: What might be a better way to test the occurrence
# of homophily or segregation in a network? How might we code that in R?
#####

#####
# Tau statistic (code by Sam Pimentel)
#####


#R code for generating random graphs:
#requires packages ergm, intergraph

#set up weighting vectors for clustering and hierarchy
clust.mask <- rep(0,16)
clust.mask[c(1,3,16)] <- 1
hier.mask <- rep(1,16)
hier.mask[c(6:8,10:11)]  <- 0

#compute triad count and triad proportion for a given weighting vector
mask.stat <- function(my.graph, my.mask){
    n.nodes <- vcount(my.graph)
    n.edges <- ecount(my.graph)
    #set probability of edge formation in random graph to proportion of possible edges present in original
    p.edge <- n.edges/(n.nodes*(n.nodes +1)/2)
    r.graph <- as.network.numeric(n.nodes, density = p.edge)
    r.igraph <- as.igraph(r.graph)
    tc.graph <- triad.census(r.igraph)
    clust <- sum(tc.graph*my.mask)
    clust.norm <- clust/sum(tc.graph)
    return(c(clust,clust.norm))
}

#build 100 random graphs and compute their clustering and hierarchy measurements to create an empirical null distribution
emp.distro <- function(this.graph){
  clust <- matrix(rep(0,200), nrow=2)
  hier <- matrix(rep(0,200),nrow=2)
  for(i in c(1:100)){
     clust[,i] <- mask.stat(this.graph, clust.mask)
     hier[,i] <- mask.stat(this.graph, hier.mask)
  }
  my.mat <- rbind(clust, hier)
  rownames(my.mat) <- c("clust.ct", "clust.norm", "hier.ct", "hier.ct.norm")
  return(my.mat)
}

#fix randomization if desired so results are replicable
#set.seed(3123)
#compute empirical distributions for each network
hc_advice <- emp.distro(krack_advice)
hc_friend <- emp.distro(krack_friendship)
hc_report <- emp.distro(krack_reports_to)

#find empirical p-value
get.p <- function(val, distro)
{
    distro.n <- sort(distro)
    distro.n <- distro.n - median(distro.n)
    val.n <- val - median(distro.n)
    p.val <- sum(abs(distro.n) > abs(val.n))/100
    return(p.val)
}
get.p(198, hc_full[1,])
get.p(194, hc_advice[1,])
get.p(525, hc_friend[1,])
get.p(1003, hc_report[1,])
get.p(979, hc_full[3,])
get.p(1047, hc_advice[3,])
get.p(1135, hc_friend[3,])
get.p(1314, hc_report[3,])

#generate  95% empirical confidence intervals for triad counts

#clustering
c(sort(hc_advice[1,])[5], sort(hc_advice[1,])[95])
c(sort(hc_friend[1,])[5], sort(hc_friend[1,])[95])
c(sort(hc_report[1,])[5], sort(hc_report[1,])[95])

#hierarchy
c(sort(hc_advice[3,])[5], sort(hc_advice[3,])[95])
c(sort(hc_friend[3,])[5], sort(hc_friend[3,])[95])
c(sort(hc_report[3,])[5], sort(hc_report[3,])[95])



########################################################################
## Title: SPARQL with R in less than 5 min
## Date: 2013-02-02
## Source: http://www.programmingr.com/content/sparql-with-r/
########################################################################

library(SPARQL) # SPARQL querying package
library(ggplot2)

# Step 1 - Set up preliminaries and define query
# Define the data.gov endpoint
endpoint <- "http://services.data.gov/sparql"

# create query statement
query <-
"PREFIX  dgp1187: <http://data-gov.tw.rpi.edu/vocab/p/1187/>
SELECT ?ye ?fi ?ac
WHERE {
?s dgp1187:year ?ye .
?s dgp1187:fires ?fi .
?s dgp1187:acres ?ac .
}"

# Step 2 - Use SPARQL package to submit query and save results to a data frame
qd <- SPARQL(endpoint,query)
df <- qd$results

# Step 3 - Prep for graphing

# Numbers are usually returned as characters, so convert to numeric and create a
# variable for "average acres burned per fire"
str(df)
df <- as.data.frame(apply(df, 2, as.numeric))
str(df)

df$avgperfire <- df$ac/df$fi

# Step 4 - Plot some data
ggplot(df, aes(x=ye, y=avgperfire, group=1)) +
geom_point() +
stat_smooth() +
scale_x_continuous(breaks=seq(1960, 2008, 5)) +
xlab("Year") +
ylab("Average acres burned per fire")

ggplot(df, aes(x=ye, y=fi, group=1)) +
geom_point() +
stat_smooth() +
scale_x_continuous(breaks=seq(1960, 2008, 5)) +
xlab("Year") +
ylab("Number of fires")

ggplot(df, aes(x=ye, y=ac, group=1)) +
geom_point() +
stat_smooth() +
scale_x_continuous(breaks=seq(1960, 2008, 5)) +
xlab("Year") +
ylab("Acres burned")

# In less than 5 mins we have written code to download just
# the data we need and have an interesting result to explore!


########################################################################
## Title: Arc diagram
## Date: 2013-02-05
## Source: http://www.r-bloggers.com/arc-diagrams-in-r-les-miserables/
########################################################################

## load devtools
library(devtools)

## install arcdiagram
install_github('arcdiagram', username='gastonstat')

## load arcdiagram
library(arcdiagram)


## location of 'gml' file
mis_file = "lesmiserables.txt"

## read 'gml' file
mis_graph = read.graph(mis_file, format="gml")

## get edgelist
edgelist = get.edgelist(mis_graph)

## get vertex labels
vlabels = get.vertex.attribute(mis_graph, "label")

## get vertex groups
vgroups = get.vertex.attribute(mis_graph, "group")

## get vertex fill color
vfill = get.vertex.attribute(mis_graph, "fill")

## get vertex border color
vborders = get.vertex.attribute(mis_graph, "border")

## get vertex degree
degrees = degree(mis_graph)

## get edges value
values = get.edge.attribute(mis_graph, "value")

## load reshape
library(reshape)

## data frame with vgroups, degree, vlabels and ind
x = data.frame(vgroups, degrees, vlabels, ind=1:vcount(mis_graph))

## arranging by vgroups and degrees
y = arrange(x, desc(vgroups), desc(degrees))

## get ordering 'ind'
new_ord = y$ind

## plot arc diagram
arcplot(edgelist, ordering=new_ord, labels=vlabels, cex.labels=0.8,
        show.nodes=TRUE, col.nodes=vborders, bg.nodes=vfill,
        cex.nodes = log(degrees)+0.5, pch.nodes=21,
        lwd.nodes = 2, line=-0.5,
        col.arcs = hsv(0, 0, 0.2, 0.25), lwd.arcs = 1.5 * values)





########################################################################
## Title: Check the data of ASTI
## Date: 2013-02-05
########################################################################

library(FAOSTAT)

tmp = read.csv("asti_agr_spending.csv", header = TRUE, encoding = "UTF-8")

tmp2 = fillCountryCode(country = "country", data = tmp)

## Check which China they are using
tmp2[tmp2$country = "China", "FAOST_CODE"] = 41
tmp2[tmp2$country = "Congo, Republic of", "FAOST_CODE"] = 46
tmp2[tmp2$country = "Czech Rep", "FAOST_CODE"] = 167
tmp2[tmp2$country = "Gambia, The", "FAOST_CODE"] = 75
tmp2[tmp2$country = "Korea, Rep", "FAOST_CODE"] = 117
tmp2[tmp2$country = "Slovak Rep", "FAOST_CODE"] = 199
tmp2[tmp2$country = "Sudan", "FAOST_CODE"] = 206
tmp2[tmp2$country = "United States", "FAOST_CODE"] = 231
tmp2[tmp2$country = "Vietnam", "FAOST_CODE"] = 237




########################################################################
## Title: Relearning box plot and labelling outliers
## Date: 2013-02-05
########################################################################

if(!is.element("FAOSTAT", .packages()))
    install.packages("FAOSTAT")

library(FAOSTAT)

## Download data on Cassava production
cp.lst = getFAOtoSYB(name = "cassava_production", domainCode = "QC",
    itemCode = 125, elementCode = 5510)

## Take the country level data, take only data for 2011 and remove NA's
cp.df = cp.lst$entity[!is.na(cp.lst$entity$cassava_production) &
                      cp.lst$entity$Year == 2011, ]

## Merge with the country profile to obtain the name
ccp.df = merge(cp.df, FAOcountryProfile[, c("FAOST_CODE", "ABBR_FAO_NAME")],
    all.x = TRUE)

## Merge with the regional pofile to obtain the UNSD M49 macro region
## composition.
rcp.df = merge(ccp.df, FAOregionProfile[, c("FAOST_CODE", "UNSD_MACRO_REG")],
    all.x = TRUE)

## Compute the quantile
qrcp.df = ddply(.data = rcp.df, .variables = .(UNSD_MACRO_REG), transform,
    lQntl = quantile(cassava_production, probs = 0.25, na.rm = TRUE),
    uQntl = quantile(cassava_production, probs = 0.75, na.rm = TRUE))

## Compute the lower and upper bound
brcp.df = ddply(.data = qrcp.df, .variables = .(UNSD_MACRO_REG), transform,
    lBound = lQntl - 1.5 * (uQntl - lQntl),
    uBound = uQntl + 1.5 * (uQntl - lQntl))

## Remove the country names if it is within the bounds
with(brcp.df, {
    brcp.df[cassava_production <= uBound &
            cassava_production >= lBound, "ABBR_FAO_NAME"] <<- NA
})


## Plot the data
set.seed(587)
ggplot(data = brcp.df, aes(x = UNSD_MACRO_REG, y = cassava_production)) +
    geom_boxplot(outlier.colour = NA) +
    geom_text(aes(label = ABBR_FAO_NAME), size = 2,
              position = position_jitter(width = 0.1)) +
    labs(x = NULL, y = NULL, title = "Production of Cassava by region")


########################################################################
## Title: Checking the relationship between mortality and population
##        growth.
## Date: 2013-02-08
########################################################################

library(FAOSTAT)

check = getWDItoSYB(indicator = c("SP.POP.TOTL", "SP.DYN.TFRT.IN",
                    "SP.DYN.IMRT.IN", "SP.DYN.LE00.IN"))
check2 = translateCountryCode(check$entity, from = "ISO2_WB_CODE",
    to = "FAOST_CODE")

check3 = merge(check2, grConstruct(data = check2, origVar = "SP.POP.TOTL",
    type = "geo", n = 10))
check3 = merge(check3, FAOregionProfile[, c("FAOST_CODE", "UNSD_MACRO_REG")],
    all.x = TRUE)

ggplot(data = check3[check3$Year == 2010, ],
       aes(x = SP.POP.TOTL.GR1, y = SP.DYN.IMRT.IN)) +
    geom_point(aes(col = UNSD_MACRO_REG))

ggplot(data = check3[check3$Year == 2010, ],
       aes(x = SP.POP.TOTL.GR1, y = SP.DYN.TFRT.IN)) +
    geom_point(aes(col = UNSD_MACRO_REG))

ggplot(data = check3,
       aes(x = SP.POP.TOTL.GR1, y = SP.DYN.LE00.IN)) +
    geom_point(aes(col = UNSD_MACRO_REG))

ggplot(data = check3[check3$Country == "Italy", ],
       aes(x = Year, y = SP.DYN.TFRT.IN)) +
    geom_point(size = 1.5) + geom_line()

ggplot(data = check3[check3$Country == "Italy", ],
       aes(x = SP.POP.TOTL.GR1, y = SP.DYN.TFRT.IN)) +
    geom_path(aes(col = Year))

ggplot(data = check3[check3$Country == "South Africa", ],
       aes(x = SP.POP.TOTL.GR10, y = SP.DYN.TFRT.IN)) +
    geom_path(aes(col = Year))




########################################################################
## Title: How are the rest of the Oceania differ from Australia and
##        New Zealand in terms of agriculture?
## Date: 2013-02-11
########################################################################

library(FAOSTAT)

## Download the data
ocQuery = list(name = c("agriLand", "arableLand", "permtCrop", "permtPast"),
               domainCode = c("RL", "RL", "RL", "RL"),
               elementCode = c(5110, 5110, 5110, 5110),
               itemCode = c(6610, 6621, 6650, 6655))

raw.lst = getFAOtoSYB(query = ocQuery)

## use country level data
raw.df = subset(raw.lst$entity)

## merge to obatain country name
raw2.df = merge(raw.df, FAOcountryProfile[, c("FAOST_CODE", "ABBR_FAO_NAME")])

## merge with regional profile
region.df = merge(raw2.df, FAOregionProfile[, c("FAOST_CODE", "UNSD_MACRO_REG")],
    all.x = TRUE)


## subset to get Oceania
oceania.df = subset(region.df, subset = UNSD_MACRO_REG == "Oceania")

## Substitute NA with zeros
oceania.df[is.na(oceania.df)] = 0


## Plot the data

## We can see that Australia dwarfs every other country in the region.
ggplot(data = oceania.df, aes(x = Year, y = agriLand)) +
    geom_line() + facet_wrap(~ABBR_FAO_NAME, scales = "fixed")

ggplot(data = oceania.df, aes(x = Year, y = agriLand)) +
    geom_line() + facet_wrap(~ABBR_FAO_NAME, scales = "free_y")

ggplot(data = oceania.df, aes(x = Year, y = arableLand)) +
    geom_line() + facet_wrap(~ABBR_FAO_NAME, scales = "free_y")

ggplot(data = oceania.df, aes(x = Year, y = permtCrop)) +
    geom_line() + facet_wrap(~ABBR_FAO_NAME, scales = "free_y")

ggplot(data = oceania.df, aes(x = Year, y = permtPast)) +
    geom_line() + facet_wrap(~ABBR_FAO_NAME, scales = "free_y")


check = subset(region.df, ABBR_FAO_NAME == "Italy")
check2 = subset(region.df, ABBR_FAO_NAME == "New Zealand")


## Examine agricultural land
ggplot(data = region.df, aes(x = Year, y = agriLand)) +
    geom_line(aes(col = ABBR_FAO_NAME)) + guides(color = guide_legend(nrow = 15)) +
    theme(legend.position = "top")


mdiff = function(x, ...){
    c(NA, diff(x, ...))
}

region2.df = ddply(region.df, .(FAOST_CODE), transform,
    dagriLand = mdiff(agriLand))

class.df = subset(region.df, subset = !is.na(agriLand) & !is.na(arableLand) &
    !is.na(permtCrop) & !is.na(permtPast) & Year == 2009)

## class.df = ddply(class.df, transform, parablLand = arablLand/agriLand, ppermtCrop = permtCrop/agriLand, ppermtPast = permtPast/agriLand)

class.df$parableLand = with(class.df, arableLand/agriLand)
class.df$ppermtCrop = with(class.df, permtCrop/agriLand)
class.df$ppermtPast = with(class.df, permtPast/agriLand)

clusterData.df = scale(class.df[, c("arableLand", "permtCrop", "permtPast",
    "parableLand", "ppermtCrop", "ppermtPast")])
wss <- (nrow(clusterData.df)-1)*sum(apply(clusterData.df,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(clusterData.df,
  	 centers=i)$withinss)
plot(1:15, wss, type = "b")


class.df$cluster = kmeans(clusterData.df, centers = 6)$cluster
with(class.df, aggregate(arableLand, FUN = mean, by = list(cluster)))
with(class.df, aggregate(permtCrop, FUN = mean, by = list(cluster)))
with(class.df, aggregate(permtPast, FUN = mean, by = list(cluster)))

with(class.df, aggregate(parableLand, FUN = mean, by = list(cluster)))
with(class.df, aggregate(ppermtCrop, FUN = mean, by = list(cluster)))
with(class.df, aggregate(ppermtPast, FUN = mean, by = list(cluster)))
table(class.df$clust)



## cluster 1 are megat agriculture
class.df[class.df$cluster == 1, "ABBR_FAO_NAME"]
table(class.df[class.df$cluster == 1, "UNSD_MACRO_REG"])

## cluster 2 are tiny countries which are more crop oriented.
class.df[class.df$cluster == 2, "ABBR_FAO_NAME"]
table(class.df[class.df$cluster == 2, "UNSD_MACRO_REG"])

## cluseter 3 are relative balanced
class.df[class.df$cluster == 3, "ABBR_FAO_NAME"]
table(class.df[class.df$cluster == 3, "UNSD_MACRO_REG"])

## cluster 4 is mainly pasture, which potentialy means they are farming country?
class.df[class.df$cluster == 4, "ABBR_FAO_NAME"]
table(class.df[class.df$cluster == 4, "UNSD_MACRO_REG"])

## cluster 5 is mainly arableLand
class.df[class.df$cluster == 5, "ABBR_FAO_NAME"]
table(class.df[class.df$cluster == 5, "UNSD_MACRO_REG"])

## cluster 6 is might be the rest, nothing particular interesting with them.
class.df[class.df$cluster == 6, "ABBR_FAO_NAME"]
table(class.df[class.df$cluster == 6, "UNSD_MACRO_REG"])



########################################################################
## Title: Use data.table with ggplot2
## Date: 2013-02-16
########################################################################

library(data.table)
library(ggplot2)

DF = data.frame(x = c("b","b","b","a","a"), v = rnorm(5))
DT = data.table(DF)
DT2 = DT

setkey(DT, x)
DT["b", ]
DT["b", mult = "first"]
DT["b", mult = "last"]

ggplot(data = DT, aes(x = x, y = v)) + geom_point()


## Examples from the data.tabl documentation
grpsize = ceiling(1e7/26^2)

dft1 = system.time(
    DF <- data.frame(
        x = rep(LETTERS,each = 26 * grpsize),
        y = rep(letters,each = grpsize),
        v = runif(grpsize * 26^2),
        stringsAsFactors = FALSE)
)

dft2 = system.time(
    ans1 <- DF[DF$x == "R" & DF$y == "h", ]
)


DT = data.table(DF)
setkey(DT,x,y)
dtt1 = system.time(
    ans2 <- DT[J("R","h")]
)


## This shows how fast the computation compares between all the
## aggregate functions. The data table is so much faster.
system.time(
    DT[, j = sum(v), by = x]
    )

system.time(
    aggregate(DF$v, by = list(DF$x), FUN = sum)
    )

system.time(
    ddply(.data = DF, .variables = .(x), .fun = function(x) sum(x$v))
    )


system.time(
    tapply(DT$v, DT$x, sum)
)



########################################################################
## Title: Bigcor: Large correlation matrices in R
## Date: 2013-02-27
## Source: http://www.r-bloggers.com/bigcor-large-correlation-matrices-in-r/
########################################################################

MAT <- matrix(rnorm(40000 * 28), nrow = 28)
cor(MAT)

corMAT <- ff(vmode = "single", dim = c(40000, 40000))

bigcor <- function(x, nblocks = 10, verbose = TRUE, ...)
{
    library(ff, quietly = TRUE)
    NCOL <- ncol(x)

    ## test if ncol(x) %% nblocks gives remainder 0
    if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")

    ## preallocate square matrix of dimension
    ## ncol(x) in 'ff' single format
    corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))

    ## split column numbers into 'nblocks' groups
    SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))

    ## create all unique combinations of blocks
    COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
    COMBS <- t(apply(COMBS, 1, sort))
    COMBS <- unique(COMBS)

    ## iterate through each block combination, calculate correlation matrix
    ## between blocks and store them in the preallocated matrix on both
    ## symmetric sides of the diagonal
    for (i in 1:nrow(COMBS)) {
        COMB <- COMBS[i, ]
        G1 <- SPLIT[[COMB[1]]]
        G2 <- SPLIT[[COMB[2]]]
        if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
        flush.console()
        COR <- cor(MAT[, G1], MAT[, G2], ...)
        corMAT[G1, G2] <- COR
        corMAT[G2, G1] <- t(COR)
        COR <- NULL
    }

    gc()
    return(corMAT)
}

bigcor(MAT)

########################################################################
## Title:
## Date:
########################################################################



install.pandoc <- function() {
    ## Credit for the source:
    ## http://stackoverflow.com/questions/15071957/is-it-possible-to-install-pandoc-on-windows-using-an-r-command

    if(require(XML)) {
        page = readLines('http://code.google.com/p/pandoc/downloads/list',
                          warn = FALSE)
        pagetree = htmlTreeParse(page, error = function(...){},
            useInternalNodes = TRUE, encoding = "UTF-8")
        url = xpathSApply(pagetree, '//tr[2]//td[1]//a ', xmlAttrs)[1]
        url = paste('http', url, sep = ':')
    } else {
        page = readLines('http://code.google.com/p/pandoc/downloads/list',
            warn = FALSE)
        pat = "//pandoc.googlecode.com/files/pandoc-[0-9.]+-setup.exe"
        target.line = grep(pat, page, value = TRUE)
        m = regexpr(pat, target.line)
        url = regmatches(target.line, m)
        ## (The http still needs to be prepended.
        url = paste('http', url, sep = ':')
    }

    pandoc_exe = tempfile(fileext = ".exe")
    download.file(url, pandoc_exe, mode = "wb")
    system(pandoc_exe)
    unlink(pandoc_exe)
}

## simply run this command to use the function:
install.pandoc()


########################################################################
## Title: Geometric growth rate vs least squares growth rate
## Date: 2013-02-28
########################################################################

library(FAOSTAT)

## Data with exponential growth rate
t = 1:10
y = exp(1.05 * t)


pop.df = getWDItoSYB(indicator = "SP.POP.TOTL")

pop2.df = translateCountryCode(pop.df$entity, from = "ISO2_WB_CODE",
    to = "FAOST_CODE")

test = grConstruct(data = pop2.df, origVar = "SP.POP.TOTL", type = "ls",
    newVarName = "SP.POP.TOTL.LSGR10", n = 10)
pop3.df = merge(pop2.df, test)
test = grConstruct(data = pop2.df, origVar = "SP.POP.TOTL", type = "geo",
    newVarName = "SP.POP.TOTL.GEOGR10", n = 10)
pop4.df = merge(pop3.df, test)

pop4.df$diff = abs(pop4.df$SP.POP.TOTL.LSGR10 - pop4.df$SP.POP.TOTL.GEOGR10)
check.df = pop4.df[which(pop4.df$diff >= 1), ]
## There doesn't appear to be significant differences if there is an
## exponential trend.




########################################################################
## Title: Extract data for Adam
## Date: 2013-04-11
########################################################################


setwd("/home/mk/Dropbox/SYBproject/GlobalBooks/SYB/2013")
## UN_CODE 156 should not be matched with Taiwan.
## Load the data manualy

fsi.df = test[test$Area == "Territory" & test$FAOST_CODE != 214,
  c("UN_CODE", "Year", "DT.OUT.UTST.POP.SH", "DO.OUT.ACPU.POP.NO")]
tfsi.df = translateCountryCode(fsi.df, from = "UN_CODE", to = "ISO3_WB_CODE")
stfsi.df = subset(tfsi.df, !is.na(ISO3_WB_CODE))
finalFSI.df = merge(stfsi.df,
  FAOcountryProfile[, c("UN_CODE", "OFFICIAL_FAO_NAME")], all.x = TRUE)



other.df = test2[test2$Area == "Territory",
  c("UN_CODE", "Year", "OFFICIAL_FAO_NAME",
"RL.AREA.ARBL.HA.SHP", "SI.POV.GAPS", "SI.POV.GAP2", "SI.POV.NAGP",
"SI.POV.RUGP", "SI.POV.URGP", "SI.POV.DDAY", "SI.POV.2DAY",
"SI.POV.NAHC", "SI.POV.RUHC", "SI.POV.URHC")]
tother.df = translateCountryCode(other.df, from = "UN_CODE", to = "ISO3_WB_CODE")
stother.df = subset(tother.df, !is.na(ISO3_WB_CODE))


final.df = merge(finalFSI.df, stother.df, all = TRUE)
final.df$OFFICIAL_FAO_NAME = NULL
final.df = merge(final.df, FAOcountryProfile[, c("ISO3_WB_CODE", "LAB_NAME")],
  all.x = TRUE)
final.df = final.df[final.df$Year %in% 2000:2012, ]

FAOcountryProfile[FAOcountryProfile$UN_CODE %in%
                  unique(finalFSI.df$UN_CODE[!(finalFSI.df$UN_CODE %in%
                                               other.df$UN_CODE)]), ]

write.csv(final.df, "~/Desktop/adamRequest.csv", row.names = FALSE, na = "")


########################################################################
## Title: Piece wise regression spline
########################################################################

N <- 40 # number of sampled points
K <- 5  # number of knots

piece.formula <- function(var.name, knots) {
  formula.sign <- rep(" - ", length(knots))
  formula.sign[knots < 0] <- " + "
  paste(var.name, "+",
        paste("I(pmax(", var.name, formula.sign, abs(knots), ", 0))",
              collapse = " + ", sep=""))
}

f <- function(x) {
    2 * sin(6 * x)
}

set.seed(1)
x <- seq(-1, 1, len = N)
y <- f(x) + rnorm(length(x))

knots <- seq(min(x), max(x), len = K + 2)[-c(1, K + 2)]
model <- lm(formula(paste("y ~", piece.formula("x", knots))))

plot(x, y)
abline(v = knots)

par(mar = c(4, 4, 1, 1))
plot(x, y)
lines(x, f(x))
new.x <- seq(min(x), max(x) ,len = 10000)
points(new.x, predict(model, newdata = data.frame(x = new.x)),
      col = "red", pch = ".")
points(knots, predict(model, newdata = data.frame(x = knots)),
       col = "red", pch = 18)






########################################################################
## Title: Testing lme model on wheat production
## Date: 2013-04-26
########################################################################


library(FAOSTAT)
library(nlme)

## Function to compute changes from year to year
diffv = function(x){
    T = length(x)
    if(sum(!is.na(x)) >= 2){
      tmp = c(x[2:T]/x[1:(T - 1)])
    } else {
      tmp = rep(NA, length(x) - 1)
    }
    tmp - 1
}

mySpline = function(x){
    tmp = try(smooth.spline(x))
    if(!inherits(tmp, "try-error")){
        final = fitted(tmp)
    } else {
        final = rep(NA, length(x))
    }
    as.numeric(final)
}

myItem = 15
myItemName =
    unique(FAOmetaTable$itemTable[FAOmetaTable$itemTable$itemCode ==
                                  myItem, "itemName"])
tmp.df = getFAOtoSYB(name = "Yield", domainCode = "QC",
    itemCode = myItem, elementCode = 5419)$entity
## colnames(tmp.df) = c("FAOST_CODE", "Year", "Area", "Yield", "Production",
##         "Seed")
## mtmp.df = melt(tmp.df, id.var = c("FAOST_CODE", "Year"))
mmtmp.df = merge(tmp.df, FAOregionProfile[, c("FAOST_CODE",
    "UNSD_MACRO_REG_CODE", "UNSD_SUB_REG_CODE")], all.x = TRUE)
mmtmp.df = merge(mmtmp.df, unique(na.omit(FAOregionProfile[,
    c("UNSD_MACRO_REG_CODE", "UNSD_MACRO_REG")])), all.x = TRUE)
mmtmp.df = merge(mmtmp.df, unique(na.omit(FAOregionProfile[,
    c("UNSD_SUB_REG_CODE", "UNSD_SUB_REG")])), all.x = TRUE)
mmtmp.df$UNSD_MACRO_REG_CODE= NULL
mmtmp.df$UNSD_SUB_REG_CODE= NULL
why.dt = data.table(mmtmp.df)

## why.dt = tmp.dt[variable == "Yield", ]
why.dt[is.na(UNSD_SUB_REG), UNSD_SUB_REG := "Unknown"]
why.dt[is.na(UNSD_SUB_REG), UNSD_MACRO_REG := "Unknown"]
why.dt[, yieldCh := c(NA, diffv(Yield)), by = "FAOST_CODE"]
why.dt[, groupYieldCh := mean(yieldCh, na.rm = TRUE),
       by = c("UNSD_SUB_REG", "Year")]
why.dt[is.na(groupYieldCh), groupYieldCh := 0]
why.dt[, averageYield := mean(Yield, na.rm = TRUE),
       by = c("UNSD_SUB_REG", "Year")]
why.dt[, medianYield := median(Yield, na.rm = TRUE),
       by = c("UNSD_SUB_REG", "Year")]
why.dt[, nCountry := length(unique(FAOST_CODE)),
       by = c("UNSD_SUB_REG", "Year")]
why.dt[, splineYield :=
       ## fitted(smooth.spline(x = cbind(averageYield, medianYield))),
       fitted(smooth.spline(x = Year, y = averageYield)),
       by = "UNSD_SUB_REG"]


why.fit = lme(Yield ~ Year * UNSD_SUB_REG + averageYield,
    random = ~1 + Year|FAOST_CODE, data = why.dt, na.action = na.omit)
why.dt[!is.na(Yield), fittedGroup := fitted(why.fit, level = 0)]
why.dt[!is.na(Yield), fittedSeries := fitted(why.fit, level = 1)]


ggplot(why.dt, aes(x = Year, y = Yield)) +
    ## geom_line(aes(x = Year, y = fittedGroup), col = "black", lwd = 1.1) +
    ## geom_line(aes(x = Year, y = averageYield), col = "red", lwd = 1.1) +
    ## geom_line(aes(x = Year, y = medianYield), col = "green", lwd = 1.1) +
    geom_line(aes(col = factor(FAOST_CODE))) +
    geom_line(aes(x = Year, y = splineYield),
              col = "steelblue",
              lwd = 1.5, alpha = 0.7) +
    scale_color_manual(values = c(
                       rep("black", length(unique(why.dt$FAOST_CODE))))) +
    facet_wrap(~UNSD_SUB_REG, ncol = 4) +
    theme(legend.position = "none") +
    labs(y = paste0("Yield of ", myItemName, " (", myItem, ")",
         " by sub-region"), x = NULL)

ggplot(why.dt,
       aes(x = Year, y = Yield - fittedSeries)) +
    geom_line(aes(col = factor(FAOST_CODE))) +
    facet_wrap(~UNSD_SUB_REG, ncol = 4) +
    theme(legend.position = "none")

ggplot(why.dt, aes(x = Year, y = groupYieldCh)) +
    geom_line() +
    facet_wrap(~UNSD_SUB_REG, ncol = 4) +
    theme(legend.position = "none")



ggplot(why.dt, aes(x = Year, y = averageYield)) +
    geom_line() +
    facet_wrap(~UNSD_SUB_REG, ncol = 4) +
    theme(legend.position = "none")


########################################################################
## Title: Example for FAOSTAT analysis section (Part 1 - growth)
## Date: 2013-05-03
########################################################################

library(FAOSTAT)

wp.lst = getFAOtoSYB(name = "wheatProd", domainCode = "QC",
  elementCode = 5510, itemCode = 15)
wp.df = wp.lst$entity

usa.df = wp.df[wp.df$FAOST_CODE == 231, ]

expgr = function(x, n = 1){
  T = length(x)
  if(sum(is.na(x)) == T){
    expgr = rep(NA, T)
    warning("All values are NA")
  } else {
    firstObs = ifelse(any(is.na(x)), min(which(!is.na(x))), 1)
    if(n > T - firstObs - 1){
      expgr = rep(NA, T)
      warning("Time series does not have sufficient values")
    } else {
      if(sum(is.na(x[firstObs:T])) > 0.5 * (T - firstObs + 1)){
        expgr = rep(NA, T)
        warning("Over 50% of the data are missing")
      } else {
        expgr = double(T)
        expgr[1:(firstObs + n - 1)] = NA
        expgr[(firstObs + n):T] = (log(x[(firstObs + n):T]/
                                       x[firstObs:(T - n)])/n) * 100
      }
    }
  }
  expgr
}

compareGrowth = function(Data, year, value, n){
  l = length(Data[, value])
  Data[, "ls"] = lsgr(Data[, value], n = n)
  Data[, "geo"] = geogr(Data[, value], n = n)
  Data[, "exp"] = expgr(Data[, value], n = n)
  par(mfrow = c(2, 1))
  maximum = max(abs(Data[, c("ls", "geo", "exp")]), na.rm = TRUE)
  plot(Data[, year], Data[, value], type = "b", xlab = "", ylab = value)
  myFormulae = paste0("log(", value, ") ~ ", year)
  lsFit = exp(fitted(lm(myFormulae, data = Data)))
  lines(Data[, year], lsFit, col = "red")
  lines(Data[c(1, l), year], Data[c(1, l), value], col = "steelblue")
  expR = log(Data[l, value]/Data[1, value])/l
  lines(Data[, year], Data[1, value] * exp(expR * 1:l), col = "orange")
  legend("topleft", legend = c("Least squares", "Geometric", "Exponential"),
         col = c("red", "steelblue", "orange"), lwd = 1, lty = 1, bty = "n")
  plot(Data[, year], Data[, value], type = "n",
       ylim = c(-maximum, maximum),
       xlab = "", ylab = "")
  lines(Data[, year], Data[, "ls"], col = "red")
  lines(Data[, year], Data[, "geo"], col = "steelblue")
  lines(Data[, year], Data[, "exp"], col = "orange")
  abline(h = 0, lty = 2)
}

compareGrowth(usa.df, "Year", "wheatProd", 10)

wheatExVal.lst = getFAOtoSYB(name = "wheatExVal", domainCode = "TP",
  itemCode = 15, elementCode = 5922)

wev.df = wheatExVal.lst$entity
usa.df = wev.df[wev.df$FAOST_CODE == 231, ]
compareGrowth(usa.df, "Year", "wheatExVal", 20)

## Function to compute the growth rates
computeGr = function(series, n.period){
  require(FAOSTAT)
  n = length(series)
  cbind(lsgr = lsgr(series, n = n.period), geogr = geogr(series, n = n.period))
}

## Function to compute growth rate fit
computeGrFit = function(series){
  observed = range(which(!is.na(series)))
  t = observed[1]:observed[2]
  ls_fit = exp(fitted(lm(log(series) ~ t)))
  geo_fit = c(rep(NA, observed[1] - 1),
    series[observed[1]] * cumprod(rep((series[observed[2]]/series[observed[1]])^
    (1/(observed[2] - observed[1])), length(t))),
    rep(NA, length(series) - observed[2]))
  cbind(lsfit = ls_fit, geofit = geo_fit)
}

wheat.ts = usa.df$wheatProd

gr = computeGr(wheat.ts, n.period = 10)
fit = computeGrFit(wheat.ts)

plot(wheat.ts, type = "b")
lines(fit[, 1], col = "red")
lines(fit[, 2], col = "blue")

plot(gr[, 1], type = "l", col = "red")
lines(gr[, 2], col = "blue")


########################################################################
## Title: LME chapter 1
## Date: 2013-05-06
## Link: http://lme4.r-forge.r-project.org/book/Ch1.pdf
########################################################################

library(lme4)
str(Dyestuff)

## Fit the first model
fm1 = lmer(Yield ~ 1 + (1|Batch), data = Dyestuff)
print(fm1)

fm1 = lmer(Yield ~ 1 + (1|Batch), data = Dyestuff, REML = FALSE)
print(fm1)

## Fit the second model
(fm2 <- lmer(Yield ~ 1 + (1|Batch), Dyestuff2))
(fm2ML <- update(fm2, REML = FALSE))

## Fit first model with MLE and results printed
fm1ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE, verbose = TRUE)
env(fm1ML)$Lambda

## Fix and random effect coefficient
fixef(fm1ML)
ranef(fm1ML)

pr1 <- profile(fm1ML)




########################################################################
## Title: Example for FAOSTAT analysis section (Part 2 - Trend)
## Date: 2013-05-03
########################################################################

library(FAOSTAT)
library(forecast)

wp.lst = getFAOtoSYB(name = "wheatProd", domainCode = "QC",
  elementCode = 5510, itemCode = 15)
wp.df = wp.lst$entity
wheat.ts = ts(usa.df$wheatProd, freq = 1, start = c(1961, 1))
usa.df = wp.df[wp.df$FAOST_CODE == 231, ]
par(mfrow = c(2, 1))

## Smoothing/filter
wheat.los = with(usa.df, lowess(Year, wheatProd))
with(usa.df, plot(Year, wheatProd, type = "b"))
lines(wheat.los, col = "red")
plot(usa.df$wheatProd - wheat.los$y, type = "b")

wheat.ma = ma(usa.df$wheatProd, order = 5)
with(usa.df, plot(Year, wheatProd, type = "b"))
lines(usa.df$Year, wheat.ma, col = "red")
plot(usa.df$Year, usa.df$wheatProd - wheat.ma, type = "b")

wheat.spl = with(usa.df, spline(Year, wheatProd, n = length(wheatProd)))
with(usa.df, plot(Year, wheatProd, type = "b"))
with(wheat.spl, lines(x = x, y = y, col = "red"))
plot(usa.df$Year, usa.df$wheatProd - wheat.spl$y, type = "b")

myKernel = kernel("daniell", c(1, 2))
wheat.ker = kernapply(wheat.ts, myKernel)
plot(wheat.ts, type = "b")
lines(wheat.ker, col = "red")
plot(wheat.ts - wheat.ker, type = "b")

## Model based
wheat.lm = tslm(wheat.ts ~ trend)
with(usa.df, plot(Year, wheatProd, type = "b"))
lines(usa.df$Year, fitted(wheat.lm), col = "red")
plot(usa.df$Year, residuals(wheat.lm), type = "b")

wheat.arima = auto.arima(usa.df$wheatProd)
with(usa.df, plot(Year, wheatProd, type = "b"))
lines(usa.df$Year, fitted(wheat.arima), col = "red")
plot(usa.df$Year, residuals(wheat.arima), type = "b")

wheat.sts = StructTS(usa.df$wheatProd, type = "trend")
with(usa.df, plot(Year, wheatProd, type = "b"))
lines(usa.df$Year, fitted(wheat.sts)[, 1], col = "red")
plot(usa.df$Year, residuals(wheat.sts), type = "b")

wheat.net = nnetar(usa.df$wheatProd)
with(usa.df, plot(Year, wheatProd, type = "b"))
lines(usa.df$Year, fitted(wheat.net), col = "red")
plot(usa.df$Year, residuals(wheat.net), type = "b")





########################################################################
## Title: Simulate 3d data for plot
## Date: 2013-05-17
########################################################################

library(scatterplot3d)

smallFarm.df = data.frame(x = rep("small", 100),
  y = runif(100, 0, 30),
  z = runif(100, 0, 0.3),
  production = rnorm(100, 1, 0.3))
medFarm.df = data.frame(type = rep("Medium", 100),
  size = runif(100, 100, 350),
  commercialization = runif(100, 0.2, 0.5),
  production = rnorm(100, 2, 0.5))
largeFarm.df = data.frame(type = rep("Large", 100),
  size = runif(100, 500, 1000),
  commercialization = runif(100, 0.8, 1),
  production = rnorm(100, 3, 1))
unknownFarm.df = data.frame(type = rep("Unknown", 20),
  size = runif(20, 0, 1000),
  commercialization = runif(20, 0, 1),
  production = rnorm(20, 2, 2))

n = 30
hobbyFarm.df = data.frame(x = runif(n, 0, 0.05),
  y = rep(0, n), z = rep(0, n))
famFarm.df = data.frame(x = runif(n, 0.05, 0.2),
  y = rep(0, n), z = runif(n, 0, 0.2))
medFarm.df = data.frame(x = runif(n, 0.2, 0.4),
  y = runif(n, 0, 0.8), z = runif(n, 0.2, 0.5))
comFarm.df = data.frame(x = runif(n, 0.08, 10),
  y = runif(n, 0.9, 1), z = runif(n, 0.5, 1))

finalFarm.df = rbind(hobbyFarm.df, famFarm.df, medFarm.df, comFarm.df)
finalFarm.df$type = factor(rep(c("hobby", "family", "medium", "commercial"),
  each = n), levels = c("hobby", "family", "medium", "commercial"))

finalFarm.df$color = with(finalFarm.df,
  rgb(rev(as.numeric(type))/length(unique(type)), 0, 0, alpha = 0.75))

pdf(file = "farm3dClassificationTest.pdf", width = 10)
farm3d = with(finalFarm.df,
     scatterplot3d(x, y, z, color = color, angle = 40, type = "h", pch = 19,
                   box = FALSE, xlab = "Scale of operation",
                   ylab = "Extent of hired Farm labour",
                   zlab = "Market participation",
                   main = 'Clustering of Crop Farms Managed by Households in "X" Agro-climatic Zone'))
legend(farm3d$xyz.convert(8, 0.2, 0.3), pch = 19, col = unique(finalFarm.df$color),
       legend = unique(finalFarm.df$type), bty = "n")
graphics.off()
system("evince farm3dClassificationTest.pdf&")


########################################################################
## Title: missForest package
## Date: 2013-05-31
########################################################################

library(missForest)

data(iris)
iris.mis = prodNA(iris, noNA = 0.1)
summary(iris.mis)
iris.imp = missForest(iris.mis)

iris.imp = missForest(iris.mis, verbose = TRUE)



########################################################################
## Title: geometric vs least square growth fit
## Date: 2013-06-04
########################################################################

## real = Exponential growth
n = 100
r = 0.02
b0 = 100
y = b0 + 10 * exp(r * 1:n) + rnorm(n, sd = 20)
T = 1:n
plot(T, y, type = "l", ylim = c(0, max(y)))
egr = (y[n]/y[1])^(1/n)
lines(T, y[1] * egr^T, col = "red")
y.fit = lm(log(y) ~ T)
lines(T, exp(fitted(y.fit)), col = "steelblue")
lgr = (y[n] - y[1])/n
lines(T, y[1] +  lgr * T, col = "orange")
sgr = (y[n] - y[1])/n^2
lines(T, y[1] +  sgr * T^2, col = "brown")
cgr = (y[n] - y[1])/n^3
lines(T, y[1] +  cgr * T^3, col = "grey")
lines(b0 + 10 * exp(r * 1:n), col = "green")
c(r, lgr, sgr, cgr, egr, exp(coef(y.fit)[2]) - 1)


## Real = Trend
n = 100
r = 0.02
b0 = 100
T = 1:n
y = b0 + 10 * T + rnorm(n, sd = 20)
plot(T, y, type = "l", ylim = c(0, max(y)))
egr = (y[n]/y[1])^(1/n)
lines(T, y[1] * egr^T, col = "red")
y.fit = lm(log(y) ~ T)
lines(T, exp(fitted(y.fit)), col = "steelblue")
lgr = (y[n] - y[1])/n
lines(T, y[1] +  lgr * T, col = "orange")
sgr = (y[n] - y[1])/n^2
lines(T, y[1] +  sgr * T^2, col = "brown")
cgr = (y[n] - y[1])/n^3
lines(T, y[1] +  cgr * T^3, col = "grey")
lines(b0 + 10 * T, col = "green")
c(r, lgr, sgr, cgr, egr, exp(coef(y.fit)[2]) - 1)


## Testing Benford's distribution for exponential growth
n = 5000
r = 0.02
b0 = 0
y = b0 + 10 * exp(r * 1:n) + rnorm(n, sd = 20)
ychar = as.character(log(y))
hist(as.numeric(substr(ychar, 1, 1)), freq = FALSE, breaks = 20)
points(1:9 + 0.75, log(1 + 1/(1:9), base = 10), pch = 19, col = "red")


test = 1:n * cumprod(sample(c(1/2, 2), n, replace = TRUE))
plot.ts(test)

cdist = as.integer(substr(as.character(test), 1, 1))
hist(cdist[cdist != 0], breaks = 20, freq = FALSE)
points(1:9, log(1 + 1/(1:9), base = 10), pch = 19, col = "red")





########################################################################
## Title: Check out how many people downloaded my package
## Date: 2013-06-14
########################################################################


## (1) Download all files
## ---------------------------------------------------------------------

# Here's an easy way to get all the URLs in R
start <- as.Date('2012-10-01')
today <- as.Date('2013-06-10')

all_days <- seq(start, today, by = 'day')

year <- as.POSIXlt(all_days)$year + 1900
urls <- paste0('http://cran-logs.rstudio.com/', year, '/', all_days, '.csv.gz')

# only download the files you don't have:
missing_days <- setdiff(as.character(all_days),
                        tools::file_path_sans_ext(dir("CRANlogs"), TRUE))

dir.create("CRANlogs")

for (i in 42:length(missing_days)) {
  print(paste0(i, "/", length(missing_days)))
  download.file(urls[i], paste0('CRANlogs/', missing_days[i], '.csv.gz'))
}


## (2) Load all file into a single file
## ---------------------------------------------------------------------


file_list <- list.files("CRANlogs", full.names=TRUE)

logs <- list()
for (file in file_list) {
	print(paste("Reading", file, "..."))
	logs[[file]] <- read.table(file, header = TRUE, sep = ",", quote = "\"",
         dec = ".", fill = TRUE, comment.char = "", as.is=TRUE)
}

# rbind together all files
library(data.table)
dat <- rbindlist(logs)

# add some keys and define variable types
dat[, date:=as.Date(date)]
dat[, package:=factor(package)]
dat[, country:=factor(country)]
dat[, weekday:=weekdays(date)]
dat[, week:=strftime(as.POSIXlt(date),format="%Y-%W")]

setkey(dat, package, date, week, country)

save(dat, file="CRANlogs/CRANlogs.RData")

# for later analyses: load the saved data.table
# load("CRANlogs/CRANlogs.RData")


## (3) Overall Anlaysis
## ---------------------------------------------------------------------

library(ggplot2)
library(plyr)

str(dat)

# Overall downloads of packages
d1 <- dat[, length(week), by=package]
d1 <- d1[order(V1), ]
d1[package=="TripleR", ]
d1[package=="psych", ]

# plot 1: Compare downloads of selected packages on a weekly basis
agg1 <- dat[J(c("TripleR", "RSA")),
            length(unique(ip_id)), by=c("week", "package")]

ggplot(agg1, aes(x=week, y=V1, color=package, group=package)) +
  geom_line() + ylab("Downloads") + theme_bw() +
  theme(axis.text.x  = element_text(angle=90, size=8, vjust=0.5))


agg1 <- dat[J(c("psych", "TripleR", "RSA")),
            length(unique(ip_id)), by=c("week", "package")]

ggplot(agg1, aes(x=week, y=V1, color=package, group=package)) +
  geom_line() + ylab("Downloads") + theme_bw() +
  theme(axis.text.x  = element_text(angle=90, size=8, vjust=0.5))




########################################################################
## Title: Visualise correlation matrix for FSI
## Date: 2013-06-22
########################################################################

library(plyr)
library(reshape2)
library(ggplot2)
library(FAOSTAT)
library(car)

tmp.df = read.csv(file = "fsi_correlation_matrix.csv", stringsAsFactors = FALSE,
  check.names = FALSE)
fsiCor.df = subset(tmp.df, Variable != "" & Tests == "Pearson Correlation")
fsiCor.df$Tests = NULL
fsiCor.df[, -1] = lapply(fsiCor.df[, -1],
           FUN = function(x) as.numeric(gsub("\\*", "", x)))

tmp = data.matrix(fsiCor.df[, -1])
tmp[upper.tri(tmp, diag = TRUE)] = NA
fsiCor.df[, -1] = tmp

mfsiCor.df = melt(fsiCor.df, id.var = "Variable")
colnames(mfsiCor.df) = c("Var1", "Var2", "rho")
mfsiCor.df$rho = round(mfsiCor.df$rho * 100)
mfsiCor.df$Var2 = factor(mfsiCor.df$Var2, levels = unique(mfsiCor.df$Var1))
mfsiCor.df$Var1 = factor(mfsiCor.df$Var1, levels = rev(unique(mfsiCor.df$Var1)))
mfsiCor.df$label = ifelse(abs(mfsiCor.df$rho) >= 75,
  as.character(mfsiCor.df$rho), "")
mfsiCor.df = subset(mfsiCor.df, !(Var1 %in% c("Dietary energy supply", "Average Dietary Energy Supply Adequacy", "Share of dietary energy supply derived from cereals, roots and tubers", "Average protein supply", "Average supply of protein of animal origin", "Poverty headcount ratio at 1.25$ a day", "Depth of the food deficit", "Prevalence of food inadequacy", "Average Value of Food Production")) & !(Var2 %in% c("Dietary energy supply", "Average Dietary Energy Supply Adequacy", "Share of dietary energy supply derived from cereals, roots and tubers", "Average protein supply", "Average supply of protein of animal origin", "Poverty headcount ratio at 1.25$ a day", "Depth of the food deficit", "Prevalence of food inadequacy", "Cereal import dependency ratio")), droplevels = TRUE)

## mfsiCor.df = subset(mfsiCor.df, Var2 != "Poverty headcount ratio at 1.25$ a day" &
##   Var1 != "Average Dietary Energy Supply Adequacy", droplevels = TRUE)
mfsiCor.df$Var1 = droplevels(mfsiCor.df$Var1)
mfsiCor.df$Var2 = droplevels(mfsiCor.df$Var2)

pdf(file = "fsiCorGreen.pdf", width = 10, height = 10)
print(ggplot(data = mfsiCor.df, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill = rho), size = 10) +
      scale_fill_gradient2(na.value = "white", mid = "grey90",
                           low = muted("red", c = 300),
                           high = "darkgreen", breaks = c(-100, 0, 100),
                           limits = c(-100, 100)) +
      geom_text(aes(label = label), size = 2.5, col = "white") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
              legend.position = "top",
            panel.background = element_rect(fill = "white"),
            axis.ticks.y = element_line(colour = "white")) +
      labs(x = NULL, y = NULL, fill = "Strength of Correlation"))
graphics.off()
system("evince fsiCorGreen.pdf&")


postscript(file = "fsiCorGreen.eps", width = 10, height = 10,
           paper = "special", horizontal = FALSE)
print(ggplot(data = mfsiCor.df, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill = rho), size = 10) +
      scale_fill_gradient2(na.value = "white", mid = "grey90",
                           low = muted("red", c = 300),
                           high = "darkgreen", breaks = c(-100, 0, 100),
                           limits = c(-100, 100)) +
      geom_text(aes(label = label), size = 2.5, col = "white") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
              legend.position = "top",
            panel.background = element_rect(fill = "white"),
            axis.ticks.y = element_line(colour = "white")) +
      labs(x = NULL, y = NULL, fill = "Strength of Correlation"))
graphics.off()
system("evince fsiCorGreen.eps&")

pdf(file = "fsiCorBlue.pdf", width = 10, height = 10)
print(ggplot(data = mfsiCor.df, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill = rho), size = 10) +
      scale_fill_gradient2(na.value = "white", mid = "grey90",
                           low = muted("red", c = 300),
                           breaks = c(-100, 0, 100),
                           limits = c(-100, 100)) +
      geom_text(aes(label = label), size = 2.5, col = "white") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
              legend.position = "top",
            panel.background = element_rect(fill = "white"),
            axis.ticks.y = element_line(colour = "white")) +
      labs(x = NULL, y = NULL, fill = "Strength of Correlation"))
graphics.off()
system("evince fsiCorBlue.pdf&")

postscript(file = "fsiCorBlue.eps", width = 10, height = 10,
           paper = "special", horizontal = FALSE)
print(ggplot(data = mfsiCor.df, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill = rho), size = 10) +
      scale_fill_gradient2(na.value = "white", mid = "grey90",
                           low = muted("red", c = 300),
                           breaks = c(-100, 0, 100),
                           limits = c(-100, 100)) +
      geom_text(aes(label = label), size = 2.5, col = "white") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
              legend.position = "top",
            panel.background = element_rect(fill = "white"),
            axis.ticks.y = element_line(colour = "white")) +
      labs(x = NULL, y = NULL, fill = "Strength of Correlation"))
graphics.off()
system("evince fsiCorBlue.eps&")



########################################################################
## Title: Testing the twitter API
## Date: 2013-07-13
########################################################################

library(twitteR)
library(tm)
library(wordcloud)
library(Snowball)
library(SnowballC)

reqURL <- "https://api.twitter.com/oauth/request_token"
accessURL <- "http://api.twitter.com/oauth/access_token"
authURL <- "http://api.twitter.com/oauth/authorize"
consumerKey <- "cl3jHTIy3Ysm9roLliOw"
consumerSecret <- "b33KilDqKbqvYxHE9OA4oBPkgDVUGwwOOKLfLY0uvoo"
twitCred <- OAuthFactory$new(consumerKey=consumerKey,
                             consumerSecret=consumerSecret,
                             requestURL=reqURL,
                             accessURL=accessURL,
                             authURL=authURL)
download.file(url="http://curl.haxx.se/ca/cacert.pem",
              destfile="cacert.pem")
twitCred$handshake(cainfo = system.file("CurlSSL", "cacert.pem",
                     package = "RCurl"))

registerTwitterOAuth(twitCred)

faoSearch = searchTwitteR("#FAO", n = 1000, cainfo = "cacert.pem",
    lang = "en")


getTimeStamp = function(tweets){
    tweets$created
}

test = lapply(faoSearch, FUN = getTimeStamp)



sofiSearch = searchTwitteR("#SOFI2013", n = 1000, cainfo = "cacert.pem",
    lang = "en")

## sofiSearch = searchTwitteR("#SOFI2013", n = 1000, cainfo = "cacert.pem")

tweetTime = as.POSIXct(sapply(sofiSearch,
    FUN = function(x) as.POSIXct(getTimeStamp(x))), tz = "GMT",
    origin = "1970-01-01")

## Plot the number of tweets against time
tweetTime.df = data.frame(tweetTime = tweetTime, length(tweetTime):1)
plot(tweetTime.df, type = "l", ylab = "Number of tweets about SOFI 2013")
plot(tweetTime.df[tweetTime.df$tweetTime >=
                  as.POSIXct("2013-10-01 00:00:00 GMT", tz = "GMT"), ],
     type = "l", ylab = "Number of tweets about SOFI 2013")


sofiSearch.df = do.call("rbind", lapply(sofiSearch, as.data.frame))
myCorpus <- Corpus(VectorSource(sofiSearch.df$text))
myCorpus <- tm_map(myCorpus, tolower)
myCorpus <- tm_map(myCorpus, removePunctuation)
removeURL <- function(x) gsub("http[[:alnum:]]*", "", x)
myCorpus <- tm_map(myCorpus, removeURL)
myStopwords <- c(stopwords("english"), "available", "via")
myStopwords <- setdiff(myStopwords, c("r", "big"))
finalCorpus <- tm_map(myCorpus, removeWords, myStopwords)

myCorpusCopy <- finalCorpus
# stem words
myCorpus <- tm_map(myCorpus, stemDocument)
# inspect documents (tweets) numbered 11 to 15
# inspect(myCorpus[11:15])
# The code below is used for to make text fit for paper width

for (i in 11:15) {
    cat(paste("[[", i, "]] ", sep=""))
    writeLines(strwrap(myCorpus[[i]], width=73))
}

myCorpus <- tm_map(myCorpus, stemCompletion, dictionary=myCorpusCopy)
inspect(myCorpus[11:15])

myTdm <- TermDocumentMatrix(myCorpus, control=list(wordLengths=c(1,Inf)))

findFreqTerms(myTdm, lowfreq=10)


termFrequency <- rowSums(as.matrix(myTdm))
termFrequency <- subset(termFrequency, termFrequency>=10)
library(ggplot2)
library(reshape2)
library(plyr)
freq.df = data.frame(word = names(termFrequency), freq = termFrequency)
freq.df = arrange(freq.df, freq)
freq.df$word = factor(freq.df$word, levels = freq.df$word)
ggplot(data = freq.df[freq.df$freq >= 30, ], aes(x = word, y = freq)) +
    geom_bar() + coord_flip()

findAssocs(myTdm, "hunger", 0.25)


m <- as.matrix(myTdm)
wordFreq <- sort(rowSums(m), decreasing=TRUE)
# word cloud
set.seed(375) # to make it reproducible
grayLevels <- gray( (wordFreq+10) / (max(wordFreq)+10) )
wordcloud(words=names(wordFreq), freq=wordFreq, min.freq=3,
          random.order=F, colors=grayLevels)



########################################################################
## Title: Playing around with canonical correlation
## Date: 2013-09-16
########################################################################

require(ggplot2)
require(GGally)
require(CCA)

mm <- read.csv("http://www.ats.ucla.edu/stat/data/mmreg.csv")
colnames(mm) <- c("Control", "Concept", "Motivation", "Read", "Write",
                  "Math",  "Science", "Sex")
summary(mm)
xtabs(~Sex, data = mm)

psych <- mm[, 1:3]
acad <- mm[, 4:8]

ggpairs(psych)
ggpairs(acad)


cc1 <- cc(psych, acad)


########################################################################
## Title: Generate the famous correlation plot
## Link: http://en.wikipedia.org/wiki/File:Correlation_examples2.svg
## Date: 2013-09-20
########################################################################

library(mvtnorm)
library(RSVGTipsDevice)

MyPlot <- function(xy, xlim = c(-4, 4), ylim = c(-4, 4), eps = 1e-15) {
   title = round(cor(xy[,1], xy[,2]), 1)
   if (sd(xy[,2]) < eps) title = "" # corr. coeff. is undefined
   plot(xy, main = title, xlab = "", ylab = "",
        col = "darkblue", pch = 16, cex = 0.2,
        xaxt = "n", yaxt = "n", bty = "n",
        xlim = xlim, ylim = ylim)
}

MvNormal <- function(n = 1000, cor = 0.8) {
   for (i in cor) {
      sd = matrix(c(1, i, i, 1), ncol = 2)
      x = rmvnorm(n, c(0, 0), sd)
      MyPlot(x)
   }
}

rotation <- function(t, X) return(X %*% matrix(c(cos(t), sin(t), -sin(t), cos(t)), ncol = 2))

RotNormal <- function(n = 1000, t = pi/2) {
   sd = matrix(c(1, 1, 1, 1), ncol = 2)
   x = rmvnorm(n, c(0, 0), sd)
   for (i in t)
      MyPlot(rotation(i, x))
}

Others <- function(n = 1000) {
   x = runif(n, -1, 1)
   y = 4 * (x^2 - 1/2)^2 + runif(n, -1, 1)/3
   MyPlot(cbind(x,y), xlim = c(-1, 1), ylim = c(-1/3, 1+1/3))

   y = runif(n, -1, 1)
   xy = rotation(-pi/8, cbind(x,y))
   lim = sqrt(2+sqrt(2)) / sqrt(2)
   MyPlot(xy, xlim = c(-lim, lim), ylim = c(-lim, lim))

   xy = rotation(-pi/8, xy)
   MyPlot(xy, xlim = c(-sqrt(2), sqrt(2)), ylim = c(-sqrt(2), sqrt(2)))

   y = 2*x^2 + runif(n, -1, 1)
   MyPlot(cbind(x,y), xlim = c(-1, 1), ylim = c(-1, 3))

   y = (x^2 + runif(n, 0, 1/2)) * sample(seq(-1, 1, 2), n, replace = TRUE)
   MyPlot(cbind(x,y), xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))

   y = cos(x*pi) + rnorm(n, 0, 1/8)
   x = sin(x*pi) + rnorm(n, 0, 1/8)
   MyPlot(cbind(x,y), xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))

   xy1 = rmvnorm(n/4, c( 3,  3))
   xy2 = rmvnorm(n/4, c(-3,  3))
   xy3 = rmvnorm(n/4, c(-3, -3))
   xy4 = rmvnorm(n/4, c( 3, -3))
   MyPlot(rbind(xy1, xy2, xy3, xy4), xlim = c(-3-4, 3+4), ylim = c(-3-4, 3+4))
}

output <- function() {
   devSVGTips(width = 7, height = 3.2) # remove first and last line for no svg exporting
   par(mfrow = c(3, 7), oma = c(0,0,0,0), mar=c(2,2,2,0))
   MvNormal(800, c(1.0, 0.8, 0.4, 0.0, -0.4, -0.8, -1.0));
   RotNormal(200, c(0, pi/12, pi/6, pi/4, pi/2-pi/6, pi/2-pi/12, pi/2));
   Others(800)
   dev.off() # remove first and last line for no svg exporting
}

## Generate the plot
output()


########################################################################
## Title: Parallel line graph
## Date: 2013-09-23
########################################################################


pdf(file = "parallelCord.pdf")
test.df = data.frame(ProjectBudget = rnorm(100, 100, 20))
test.df$OO = ifelse(test.df$ProjectBudget >= 100,
    test.df$ProjectBudget * 3, test.df$ProjectBudget)
test.df$SO = test.df$OO + rnorm(100, 0, 40)
parcoord(test.df)
abline(h = 0.5, col= "red")
graphics.off()




########################################################################
## Title: Playing around with EMD and SSA
## Date: 2013-09-23
## Link: http://www.r-bloggers.com/wheres-the-magic-emd-and-ssa-in-r/
########################################################################

data(gas, package = "forecast")
library(EMD)
library(hht)
ee <- EEMD (c(gas), time (gas), 250, 100, 6, "trials")
eec <- EEMDCompile ("trials", 100, 6)



########################################################################
## Title: Test between merge_recurse and Reduce merge
## Date: 2013-10-01
########################################################################


## Create artificial data
artificial.lst = list()
for(i in 1:20){
    artificial.lst[[i]] = tips
}

merge_reduce = function(dfs){
    tmp = Reduce(f = function(x, y) merge(x, y),  x = dfs[-1],
           init = dfs[[1]])
    tmp
}

## Look like recurse is faster by 3.5 times
benchmark(merge_recurse(artificial.lst),
          merge_reduce(artificial.lst), order = "relative",
          replications = 1)


########################################################################
## Title: package and book for lme4
## Date: 2013-10-08
## Source: http://lme4.r-forge.r-project.org/lMMwR/lrgprt.pdf
########################################################################

library(lme4)
(fm01 = lmer(Yield ~ 1 + (1|Batch), Dyestuff))
(fm01ML <- lmer(Yield ~ 1 + (1|Batch), Dyestuff, REML=FALSE))
pr01 <- profile(fm01ML)
xyplot(pr01, aspect = 1.3, absVal = TRUE, conf = 0.95)
confint(pr01)


(fm02 <- lmer(Yield ~ 1 + (1|Batch), Dyestuff2))
(fm02ML <- update(fm02, REML=FALSE))
summary(fm02a <- lm(Yield ~ 1, Dyestuff2))

fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML=FALSE, verbose=1)


splom(pr01)



## Longitudinal study
str(sleepstudy)
print(xyplot(Reaction ~ Days | Subject, sleepstudy, aspect = "xy",
             layout = c(6,3), type = c("g", "p", "r"),
             index.cond = function(x,y) coef(lm(y ~ x))[1],
             xlab = "Days of sleep deprivation",
             ylab = "Average reaction time (ms)"))

(fm06 = lmer(Reaction ~ 1 + Days + (1 + Days|Subject), sleepstudy,
              REML=FALSE))




########################################################################
## Title: Rotate matrix by 90 degree
## Date: 2013-10-19
########################################################################


test = matrix(1:10, nc = 5)
t(test[, 5:1])



########################################################################
## Title: Splines
## Date: 2013-10-23
########################################################################

library(splines)
data(UKgas)
x = c(log(UKgas))
x.bs = bs(x, df = 7)
plot(x, type = "l")
for(i in 1:NCOL(x.bs)){
    lines(x.bs[, i], col = rgb(i/NCOL(x.bs), 0, 0, alpha = 0.5))
}
lines(colSums(x.bs))

plot(x, type = "l")
lines(fitted(lm(x ~ bs(x, df = 1))), col = "red")

myDim = 20
x.bs = bs(x, df = myDim)
par(mfrow = c(myDim, 1), mar = rep(0, 4))
for(i in 1:NCOL(x.bs)){
    plot(x.bs[, i], col = rgb(i/NCOL(x.bs), 0, 0, alpha = 0.5), type = "l")
}



########################################################################
## Title; Illustration of the SImpson's paradox
## Date:2013-10-24
########################################################################

x = rnorm(10)
y = 2 + 3 * x
y1 = y - 5
y2 = y1 - 10
y.df = data.frame(x, y, y1)

pdf("simpson_paradox.pdf")
plot(y ~ x, ylim = range(y, y1, y2), col = "black", pch = 19,
     xlab = "% of children under five years of age who are stunt",
     ylab = "Prevalence of undernourishment")
points(y1 ~ x, col = "blue", pch = 19)
points(y2 ~ x, col = "orange", pch = 19)
legend("topleft", legend = c("2000", "2005", "2010"), pch = 19,
       bty = "n", col = c("black", "blue", "orange"))
graphics.off()



########################################################################
## Title: Generate binomial distribution with MCMC
## Date: 2013-10-28
########################################################################

dbnorm = function(x, d, log = TRUE) {
    if (log)
        (-t(x) %*% d %*% x)/2
    else exp(-t(x) %*% d %*% x/2)
}

bnorm.mh = function(x0 = c(-5, 5), rho = 0.8,
    n = 500, FUN = dbnorm) {
    d = solve(cbind(c(1, rho), c(rho, 1)))
    chain = matrix(nrow = n + 1, ncol = 2)
    chain[1, ] = x = x0
    for (i in 2:(n + 1)) {
        z = x + runif(2, -1, 1)
        p = min(1, exp(FUN(z, d, log = TRUE) -
            FUN(x, d, log = TRUE)))
        if (runif(1) <= p)
            chain[i, ] = x = z
        else chain[i, ] = x
    }
    chain
}

plot.bnorm = function(...) {
    m = bnorm.mh(...)
    plot(m[, 1], m[, 2], type = "l", xlim = c(-5, 5), ylim = c(-5, 5),
         xlab = "x", ylab = "y", col = "blue")
    points(m[1, 1], m[1, 2], pch = 15, col = "red")
}

plot.bnorm(x0 = c(-5, 5), n = 1000)

plot(test, pch = 19, col = rgb(0, 0, 0, alpha = 0.01))
points(0, 0, col = "blue")
points(colMeans(test)[1], colMeans(test)[2], col = "red")


########################################################################
## Title: Robust Regression
## Date: 2013-10-29
########################################################################

library(car)
mod.ls = lm(prestige ~ income + education, data = Duncan)

mod.ls.2 = update(mod.ls, subset = - c(6, 16))
summary(mod.ls.2)

library(MASS)
mod.huber = rlm(prestige ~ income + education, data = Duncan)
summary(mod.huber)

plot(mod.huber$w, ylab = "Huber Weight")
smallweights = which(mod.huber$w < 0.8)
showLabels(1:45, mod.huber$w, rownames(Duncan), id.method = smallweights,
           id.cex = 0.6)

mod.bisq = rlm(prestige ~ income + education, data = Duncan,
    method = "MM")
summary(mod.bisq)

plot(mod.bisw$w, ylab = "Bisquare weight")
showLabels(1:45, mod.bisq$w, rownames(Duncan),
           id.method = which(mod.bisq$w < 0.8), id.cex = 0.6)

mod.lst = ltsreg(prestige ~ income + education, data = Duncan,
    method = "lts")
print(mod.lst)


data(coleman)
set.seed(0)
## Default for a very long time:
summary(m1 <- lmrob(Y ~ ., data=coleman))
summary(lm1 <- lm(Y ~., data = coleman))

plot(m1$rweights, ylab = "MM weights")
showLabels(1:20, m1$rweights, rownames(coleman),
           id.method = which(m1$rweights < 0.8), id.cex = 0.6)
plot(m1)

test = coleman
test$outlier = rep("black", NROW(coleman))
test[which(m1$rweights < 0.8), "outlier"] = "red"

# Scatterplot Matrices from the car Package
library(car)
scatterplotMatrix(~salaryP + fatherWc + sstatus + teacherSc +
                   motherLev + Y|outlier, data=test,
                   main="coleman data with outlier")


########################################################################
## Title: Data manipulation for pseticide data
## Date: 2013-11-05
########################################################################

library(reshape2)
library(lattice)

## Read the data
cameroonPesticide.df = read.csv(file = "cameroon_pesticide.csv",
    stringsAsFactors = FALSE)

## Normalize the data and extract the year and type from the label
normalized.df = melt(cameroonPesticide.df, id.var = "Item.description")
normalized.df$Year =
    as.numeric(gsub("[^0-9]", "", normalized.df$variable))
normalized.df$type = substring(normalized.df$variable,
    sapply(gregexpr("\\.", normalized.df$variable),
           function(x) x[[1]]) + 1)

## Plot the data
xyplot(value ~ Year|Item.description,
       data = normalized.df, type = c("g", "l"))

individualItem.df = normalized.df[!normalized.df$Item.description %in%
    c("Insecticides", "Herbicides", "Fungicides&Bactericides"), ]
aggregateItem.df = normalized.df[normalized.df$Item.description %in%
    c("Insecticides", "Herbicides", "Fungicides&Bactericides"), ]


xyplot(value ~ Year|Item.description,
       data = individualItem.df[individualItem.df$Item.description !=
           "Mineral Oils", ],
       type = c("g", "l"))

xyplot(value ~ Year|Item.description, data = aggregateItem.df,
       type = c("g", "l"))


########################################################################
## Title: A gentle introduction to Rcpp
## Date: 2013-11-06
########################################################################

## Faithful data density estimation
xx = faithful$eruptions
fit = density(xx)
plot(fit)

fit1 = density(xx)
fit2 = replicate(10000, {
    x = sample(xx, replace = TRUE)
    density(x, from = min(fit1$x), to = max(fit1$x))$y
})
fit3 = apply(fit2, 1, quantile, c(0.025, 0.975))
plot(fit1, ylim = range(fit3))
polygon(c(fit1$x, rev(fit1$x)),
        c(fit3[1, ], rev(fit3[2, ])),
        col = "grey", border = FALSE)
lines(fit1)



## Fibonacci numbers
fibR = function(n){
    if(n == 0) return(0)
    if(n == 1) return(1)
    return(fibR(n - 1) + fibR(n - 2))
}

## extern "C" SEXP fibWrapper(SEXP xs){
##     int x = Rcpp::as<int>(xs);
##     int fib = fibonacci(x);
##     return(Rcpp::wrap(fib));
## }

incltxt = "
int fibonacci(const int x){
    if(x == 0) return(0);
    if(x == 1) return(1);
    return fibonacci(x - 1) + fibonacci(x - 2);
}"

## Use inline
library(inline)
fibRcpp = cxxfunction(signature(xs = "int"),
    plugin = "Rcpp",
    incl = incltxt,
    body = "
int x = Rcpp::as<int>(xs);
return Rcpp::wrap(fibonacci(x));
")

fibR(10)
fibRcpp(10)



benchmark(fibR(20), fibRcpp(20), order = "relative", replications = 300)


########################################################################
## Title: Rcpp test
## Date: 2013-11-07
########################################################################

library(inline)

signDiffR = function(x){
    sign(diff(x))
}

incltxt = "
int fibonacci(const int x){
    if(x == 0) return(0);
    if(x == 1) return(1);
    return fibonacci(x - 1) + fibonacci(x - 2);
}"

## Use inline

fibRcpp = cxxfunction(signature(xs = "int"),
    plugin = "Rcpp",
    incl = incltxt,
    body = "
int x = Rcpp::as<int>(xs);
return Rcpp::wrap(fibonacci(x));
")



########################################################################
## Title: Iterative imputation
## Date: 2013-11-012
########################################################################

b0 = 2
b1 = 3
x = rnorm(100)
y = b0 + b1 * x + rnorm(100)
plot(y ~ x)

miss.prop = 0.2
x.miss = sample(1:length(x), size = length(x) * miss.prop)
y.miss = sample(1:length(y), size = length(x) * miss.prop)
x.new = x
x.new[x.miss] = NA
y.new = y
y.new[y.miss] = NA
plot(y.new ~ x.new)

foo = function(x, y){
    tmp = data.table(x = x, y = y)
    for(i in 1:100){
        y.fit = lm(y ~ x, data = tmp)
        tmp[!is.na(x) & is.na(y),
            y := predict(y.fit, newdata = tmp[!is.na(x) & is.na(y), ])]
        x.fit = lm(x ~ y, data = tmp)
        tmp[!is.na(y) & is.na(x),
            x := predict(x.fit, newdata = tmp[!is.na(y) & is.na(x), ])]
    }
    list(x = tmp$x, y = tmp$y)
}

test = foo(x = x.new, y = y.new)
with(test, plot(y ~ x))
points(y ~ x, col = "red", pch = 16)

plot(x ~ test$x)
plot(y ~ test$y)

mean(abs(x[is.na(x.new)] - test$x[is.na(x.new)]), na.rm = TRUE)
mean(abs(y[is.na(y.new)] - test$y[is.na(y.new)]), na.rm = TRUE)

with(test, lm(y ~ x))

########################################################################
## Title: Cyclone Haiyan with RGoogleVis
## Date: 2013-11-13
## Source: http://lamages.blogspot.it/2013/11/googlevis-047-with-rstudio-integration.html
########################################################################

library(XML)
url <- "http://www.gdacs.org/Cyclones/report.aspx?eventid=41058&episodeid=28&eventtype=TC"
dat <- readHTMLTable(readLines(url), which=5)
dat$latlon <- dat[,8]
levels(dat$latlon) <- sapply(
  strsplit(levels(dat[,8]), ",\n        "),
       function(x) paste(x[2], x[1], sep=":")
)
dat$Category <- factor(dat$Category, levels=levels(dat$Category)[c(6,7,1:5)],
                       ordered=TRUE)
dat$cat <- as.numeric(dat$Category)
dat$Gust_kmh <- dat[,6]
levels(dat$Gust_kmh) <- sapply(strsplit(levels(dat[,6]), "km"),
                               function(x) gsub(" ", "",x[1]))
dat$Gust_kmh <- as.numeric(as.character(dat$Gust_kmh))

library(googleVis)
M <- gvisGeoChart(dat, "latlon", sizevar="cat",
                  colorvar="Gust_kmh",
                  options=list(region='035',
                               backgroundColor="lightblue",
                               datalessRegionColor="grey"))
plot(M)



########################################################################
## Title: test rWBclimate
## Date: 2013-11-14
## source: https://github.com/ropensci/rWBclimatea
########################################################################



library(rWBclimate)
usa.dat <- get_model_temp(locator = "USA", type = "mavg",
                          start = 2080, end = 2100)
usa.dat.bcc <- usa.dat[usa.dat$gcm == "bccr_bcm2_0", ]
usa.dat.had <- usa.dat[usa.dat$gcm == "ukmo_hadcm3", ]
## Add a unique ID to each for easier plotting
usa.dat.bcc$ID <- paste(usa.dat.bcc$scenario, usa.dat.bcc$gcm, sep = "-")
usa.dat.had$ID <- paste(usa.dat.had$scenario, usa.dat.had$gcm, sep = "-")
plot.df <- rbind(usa.dat.bcc, usa.dat.had)

ggplot(plot.df, aes(x = as.factor(month), y = data, group = ID,
                    colour = gcm,  linetype = scenario)) +
    geom_point() + geom_path() +
    ylab("Average temperature in degrees C \n between 2080 and 2100") +
    xlab("Month") + theme_bw()


allISO3country = na.omit(FAOcountryProfile[, "ISO3_WB_CODE"])
omitList = c("ALA", "ADO", "ATA", "BES", "BVT")
allISO3country = allISO3country[!allISO3country %in% omitList]
test = get_historical_precip(locator = allISO3country,
    time_scale = "month")


test = try(get_historical_precip(locator = "BVT",
    time_scale = "month"))



########################################################################
## Title: Compare numerical integration and analytical computation
## Date: 2013-11-14
########################################################################

library(sn)
library(rbenchmark)
analytical = function(x){
    psn(x)
}

numerical = function(x){
    integrate(function(y) dsn(y), -Inf, x)
}

benchmark(analytical(100), numerical(100), order = "relative",
          replications = 10000)


integrate(function(y) dsn(y), -Inf, 50)


########################################################################
## Title: Rewriting the ADECUP function
## Date: 2013-11-14
########################################################################

library(sn)
dec = 2755.3605000
mder = 1729.7813050
cv = 0.2774188
sk = 0.8536067
cpPar = c(mean = dec, sd = cv * dec, gamma = sk)

myDsn = function(x, cpPar){
    dpPar = cp.to.dp(cpPar)
    2/dpPar["scale"] * dnorm((x - dpPar["location"])/dpPar["scale"]) *
         pnorm(dpPar["shape"] * ((x - dpPar["location"])/dpPar["scale"]))
}

dsn(dec, dp = cp.to.dp(cpPar))
myDsn(dec, cpPar = cpPar)

integrate(function(y) y * myDsn(x = y, par = cpPar), -Inf, 0)

meansn = function(cpPar){
    dpPar = cp.to.dp(cpPar)
    dpPar["location"] + dpPar["scale"] * sqrt(2/pi) *
        (dpPar["shape"]/sqrt((1 + dpPar["shape"]^2)))
}

foo = function(cpPar){
    dpPar = cp.to.dp(cpPar)
    a = meansn(cpPar)
    b = (dpPar["scale"]/sqrt(2 * pi) *
         (1 + (dpPar["shape"]/sqrt(1 + dpPar["shape"]^2))) +
         (dpPar["location"] * (pi/2 - atan(dpPar["location"])))/(2 * pi))
    print(a)
    print(b)
    a - b
}
foo(cpPar = cpPar)

integrate(function(y) y * dsn(y, dp = cp.to.dp(cpPar)), 0, mder)

(foo(cpPar = cpPar) + integrate(function(y) y * dsn(y, dp = cp.to.dp(cpPar)), 0, mder)[["value"]])/psn(mder, dp = cp.to.dp(cpPar))

## real answer
integrate(function(y) y * dsn(y, dp = cp.to.dp(cpPar)), -Inf, mder)[["value"]]/integrate(function(y) dsn(y, dp = cp.to.dp(cpPar)), -Inf, mder)[["value"]]



testPar = c(mean = 0, sd = 1, gamma = 0)
foo = function(cpPar){
    dpPar = cp.to.dp(cpPar)
    a = meansn(cpPar = cpPar)
    b = 1/(sqrt(2 * pi)) *
        (1 + (dpPar["shape"]/(1 + dpPar["shape"]^2)))
    ## print(a)
    ## print(b)
    a + b
}

integrate(function(y) y * dsn(y, shape = 2), -Inf, 0)
integrate(function(y) y * dsn(y, shape = 2), 0, Inf)
foo(cpPar = dp.to.cp(c(mean = 0, sd = 1, skewness = 2)))


adecup = function(mder, cpPar){
    dpPar = cp.to.dp(cpPar)
    k = (mder - dpPar["location"])/dpPar["scale"]
    print(k)
    k2 = foo(cpPar = cpPar)
    if(k > 0){
        tmp = (k2 + integrate(function(y)
                              y * dsn(y, shape = dpPar["shape"]),
                              0, k)[["value"]])/
                                  psn(k, shape = dpPar["shape"])
    } else {
        tmp = (k2 - integrate(function(y)
                              y * dsn(y, shape = dpPar["shape"]),
                              k, 0)[["value"]])/
                                  psn(k, shape = dpPar["shape"])
    }
    tmp
}
ADECUP(dec = 2755.3605000, mder = 1729.7813050,
       cv = 0.2774188, sk = 0.8536067)
adecup(mder = mder, cpPar = cpPar)



########################################################################
## Title: Testing foreach package
## Date: 2013-11-15
########################################################################

library(foreach)
x = foreach(i = 1:3) %do% sqrt(i)


test.df = data.frame(a = rnorm(10), b = rnorm(10))

test = foreach(i = 1:NROW(test.df)) %do% {
    with(test.df, a + b)
}


## Generate 1000 random variables

test = foreach(i = 1:10, .combine = "c") %do% {
    rnorm(100)
}

test = function(){
    tmp = foreach(i = 1:10, .combine = "c") %do% {
        rnorm(100)
    }
}


library(rbenchmark)
benchmark(rnorm(100), test, order = "relative",
          replications = 10000)



########################################################################
## Title: Relationship of FSI variables.
## Date: 2013-11-20
########################################################################

library(network)
library(reshape2)
test = read.csv(file = "FSIconstruction13.csv")
test2 = test[, c("STS_ID", "STS_ID_CONSTR1", "STS_ID_CONSTR2",
    "STS_ID_WEIGHT")]
test3 = melt(test2, id.var = "STS_ID")
test.net = network(unique(test3[test3$value != "", c(3, 1)]))

pdf(file = "~/Desktop/FSIvariableNetwork.pdf", width = 20, height = 20)
plot(test.net, displaylabels = TRUE,vertex.col = "skyblue",
     edge.col = rgb(0, 0, 0, alpha = 0.5), arrowhead.cex = 0.5)
graphics.off()
system("evince ~/Desktop/FSIvariableNetwork.pdf&")

########################################################################
## Title: Missing proportion of area sown in the statistical working
##        system
## Date: 2013-11-25
########################################################################
areaSown.df = read.csv(file = "swsAreaSown.csv", header = TRUE,
    stringsAsFactors = FALSE)

areaSown.df$validData = with(areaSown.df, ifelse(Symb %in% c(" ", "*"),
    Num, NA))

test = ddply(.data = areaSown.df, .variables = .(Item.Code, Item.Name),
      .fun = function(x) (sum(is.na(x$validData))/length(x$validData)) *
    100)
colnames(test)[3] = "missing proportion of area sown"
test2 = ddply(.data = areaSown.df, .variables = .(Item.Code, Item.Name),
      .fun = function(x) length(unique(x$Area.Code)))
colnames(test2)[3] = "Number of producing country"
test3 = merge(test, test2)
test3 = test3[order(test3[, 3]), ]
write.csv(test3, file = "~/Desktop/missingProfileAreaSown.csv",
          row.names = FALSE)


########################################################################
## Title: Test of amelia 2
## Date: 2013-11-29
## Source: http://cran.r-project.org/web/packages/Amelia/vignettes/amelia.pdf
########################################################################

library(Amelia)
data(freetrade)
summary(freetrade)
a.out = amelia(freetrade, m = 5, ts = "year", cs = "country")
hist(a.out$imputations[[3]]$tariff, col="grey", border="white")
a.out.more <- amelia(freetrade, m = 10, ts = "year", cs = "country", p2s=0)



test = swsToImputationDataTable("~/github/sws_imputation/case_study/sws_data/pepperSUA.csv", denormalizer = "Element.Code")

final.dt = test[, list(areaName, Year, areaNum, productionNum,
    yieldNum, unsdMacroReg, unsdSubReg)]

ameliaTest = amelia(final.dt, m = 5, ts = "Year", cs = "areaName",
    p2s = 2, idvars = c("unsdMacroReg", "unsdSubReg"), polytime = 2,
    intercs = TRUE)

xyplot(productionNum ~ Year|areaName, data = ameliaTest$imputations$imp1,
       type = c("g", "l"))

xyplot(productionNum ~ Year|areaName, data = test,
       type = c("g", "l"))



########################################################################
## Title: Use filter to compute weighted rolling average
## Date: 2014-01-04
########################################################################

x = 1:5
w = seq(0.1, 0.5, by = 0.1)
f = rep(1/3, 3)

filter(x * 1/3, filter = w)


k = 3
c(rep(NA, k), diff(cumsum(1:5)/k, lag = k))




########################################################################
## Title: Simple ensembling prediction
## Date: 2014-01-05
########################################################################


data(sunspot.year)
plot(sunspot.year)
sunspot.yearTsp = tsp(sunspot.year)
lines(loess(sunspot.year ~
            seq(sunspot.yearTsp[1], sunspot.yearTsp[2],
                by = 1/sunspot.yearTsp[3]),
            data = sunspot.year, span = 0.1)$fitted ~
      seq(sunspot.yearTsp[1], sunspot.yearTsp[2],
          by = 1/sunspot.yearTsp[3]),
      col = "red", lty = 2)
lines(loess(sunspot.year ~
            seq(sunspot.yearTsp[1], sunspot.yearTsp[2],
                by = 1/sunspot.yearTsp[3]),
            data = sunspot.year, span = 0.3)$fitted ~
      seq(sunspot.yearTsp[1], sunspot.yearTsp[2],
          by = 1/sunspot.yearTsp[3]),
      col = "green", lty = 2)
lines(loess(sunspot.year ~
            seq(sunspot.yearTsp[1], sunspot.yearTsp[2],
                by = 1/sunspot.yearTsp[3]),
            data = sunspot.year, span = 0.5)$fitted ~
      seq(sunspot.yearTsp[1], sunspot.yearTsp[2],
          by = 1/sunspot.yearTsp[3]),
      col = "blue", lty = 2)




########################################################################
## Title: Gibbs sampler for bivariate normal distribution from 783
## Date: 2014-01-16
########################################################################


bnorm.gibbs = function(x0 = c(-5, 5), rho = 0.8, n = 500) {
    s = sqrt(1 - rho^2)
    chain = matrix(nrow = n + 1, ncol = 2)
    chain[1, ] = x0
    for (i in 2:(n + 1)) {
        chain[i, 1] = rnorm(1, rho * chain[i - 1, 2], s)
        chain[i, 2] = rnorm(1, rho * chain[i, 1], s)
    }
    chain
}

plot.bchain = function(chain) {
    plot(chain[, 1], chain[, 2], type = "l",
         xlim = c(-5, 5), ylim = c(-5, 5),
         xlab = "x", ylab = "y", col = "blue")
    points(chain[1, 1], chain[1, 2], pch = 15,
           col = "red")
}


plot.bchain(bnorm.gibbs(x0 = c(-5, 5)))


########################################################################
## Title: Playing around with some beer data
## Date: 2014-01-16
########################################################################

test = read.csv("http://www.craftbeeranalytics.com/uploads/3/3/8/9/3389428/ratebeer_beerjobber.txt", sep = "\t", stringsAsFactors = FALSE)

by(data = test$abv, INDICES = test$style,
   FUN = function(x) mean(x, na.rm = TRUE))



histogram(~abv|style, data = test, breaks = 50)
histogram(~score.overall|style, data = test, breaks = 50)
histogram(~score.by.style|style, data = test, breaks = 50)


########################################################################
## Title: Test scripts for under water color adjustment
## Date: 2014-01-20
########################################################################

library(jpeg)
x = readJPEG("colorAdjustmentSample.jpg")
plot.new()
plot.window(xlim = c(0, 100), ylim = c(0, 100))
rasterImage(x, 0, 0, 100, 100)

## Lets try adjusting the white balance first

colorHist = function(image, breaks = 300){
    opar = par()
    par(mfrow = c(3, 1), mar = rep(0, 4))
    hist(image[, , 1], axes = FALSE, col = "red", xlab = "", main = "",
         breaks = breaks, xlim = c(0, 1))
    hist(image[, , 2], axes = FALSE, col = "green", xlab = "", main = "",
         breaks = breaks, xlim = c(0, 1))
    hist(image[, , 3], axes = FALSE, col = "blue", xlab = "", main = "",
         breaks = breaks, xlim = c(0, 1))
    par(opar)
}
colorHist(x, breaks = 70)


colorMatrix = function(image){
    opar = par()
    par(mfrow = c(3, 3), mar = rep(0, 4))
    hist(image[, , 1], axes = FALSE, breaks = 300, xlim = c(0, 1),
         col = "red", main = "")
    image(t(matrix(rev(x[, , 1] - x[, , 2]),  nc = NCOL(x[, , 1]))),
          useRaster = TRUE, axes = FALSE)
    image(t(matrix(rev(x[, , 1] - x[, , 3]),  nc = NCOL(x[, , 1]))),
          useRaster = TRUE, axes = FALSE)
    image(t(matrix(rev(x[, , 1]),  nc = NCOL(x[, , 1]))),
          useRaster = TRUE, col = rgb(0:100/100, 0, 0), axes = FALSE)
    hist(image[, , 2], axes = FALSE, breaks = 300, xlim = c(0, 1),
         col = "green", main = "")
    image(t(matrix(rev(x[, , 2] - x[, , 3]),  nc = NCOL(x[, , 1]))),
          useRaster = TRUE, axes = FALSE)
    image(t(matrix(rev(x[, , 2]),  nc = NCOL(x[, , 1]))),
          useRaster = TRUE, col = rgb(0, 0:100/100, 0), axes = FALSE)
    image(t(matrix(rev(x[, , 3]),  nc = NCOL(x[, , 1]))),
          useRaster = TRUE, col = rgb(0, 0, 0:100/100), axes = FALSE)
    ## plot(1:5, type = "n", axes = FALSE)
    ## plot(1:5, type = "n", axes = FALSE)
    hist(image[, , 3], axes = FALSE, breaks = 300, xlim = c(0, 1),
         col = "blue", main = "")
    par(opar)
}
colorMatrix(x)


tmp = x
tmp[, , 1] = (tmp[, , 1])/max(tmp[, , 1])
tmp[, , 2] = tmp[, , 2]/max(tmp[, , 2])
tmp[, , 3] = tmp[, , 3]/max(tmp[, , 3])
colorMatrix(tmp)
colorHist(tmp, breaks = 70)

plot.new()
plot.window(xlim = c(0, 110), ylim = c(0, 100))
rasterImage(x, 0, 0, 50, 100)
rasterImage(tmp, 60, 0, 110, 100)



tmp = x
tmp[, , 1] = (tmp[, , 1] - min(tmp[, , 1]))/max(tmp[, , 1])
## tmp[, , 2] =
##     (tmp[, , 2] - min(tmp[, , 2]))/(max(tmp[, , 2]/max(tmp[, , 1])))
## tmp[, , 2] = ifelse(tmp[, , 2] >= 0.1, tmp[, , 2] - 0.1, 0)
## tmp[, , 3] = (tmp[, , 3] - min(tmp[, , 3]))/max(tmp[, , 3])
## tmp[, , 3] = ifelse(tmp[, , 3] >= 0.1, tmp[, , 3] - 0.1, 0)
## tmp[, , 3] =
##     (tmp[, , 3] - min(tmp[, , 3]))/(max(tmp[, , 3]/max(tmp[, , 1])))
tmp[, , 2] = log(tmp[, , 2] + 1)
tmp[, , 3] = log(tmp[, , 3] + 1)
## tmp = ifelse(tmp <= 0.9, tmp + 0.1, 1)
colorHist(tmp, breaks = 70)
## tmp = array(tmp,
##     matrix(0, nrow = NROW(tmp[, , 1]), nc = NCOL(tmp[, , 1])),
##     dim = c(2880, 3840 ,4))
plot.new()
plot.window(xlim = c(0, 110), ylim = c(0, 100))
rasterImage(x, 0, 0, 50, 100)
rasterImage(tmp, 60, 0, 110, 100)
## rasterImage(test, 60, 0, 110, 100)



tmp = x
tmp[, , 1] = exp(tmp[, , 1]) - min(exp(tmp[, , 1]))
colorHist(tmp, breaks = 70)
## tmp = array(tmp,
##     matrix(0, nrow = NROW(tmp[, , 1]), nc = NCOL(tmp[, , 1])),
##     dim = c(2880, 3840 ,4))

plot.new()
plot.window(xlim = c(0, 110), ylim = c(0, 100))
rasterImage(x, 0, 0, 50, 100)
rasterImage(tmp, 60, 0, 110, 100)
## rasterImage(test, 60, 0, 110, 100)



test = tmp
## test[, , 3] = 0
test[, , 2] = ifelse(test[, , 2] > 0.4 & test[, , 2] < 0.6,
        test[, , 2] + 0.05, test[, , 2])
test[, , 1] = ifelse(test[, , 1] > 0.4 & test[, , 1] < 0.6,
        test[, , 1] + 0.2, test[, , 1])
plot.new()
plot.window(xlim = c(0, 110), ylim = c(0, 100))
rasterImage(tmp, 0, 0, 50, 100)
rasterImage(test, 60, 0, 110, 100)


plot.new()
plot.window(xlim = c(0, 110), ylim = c(0, 100))
rasterImage(x, 0, 0, 50, 100)
rasterImage(test, 60, 0, 110, 100)


########################################################################
## Title: Testing turboEM
## Date: 2014-01-29
########################################################################

library(turboEM)
poissmix.dat = data.frame(death = 0:9,
    freq = c(162,267,271,185,111,61,27,8,3,1))
y = poissmix.dat$freq


fixptfn = function(p, y) {
    pnew = rep(NA,3)
    i = 0:(length(y)-1)
    denom = p[1]*exp(-p[2])*p[2]^i + (1 - p[1])*exp(-p[3])*p[3]^i
    zi = p[1]*exp(-p[2])*p[2]^i / denom
    pnew[1] = sum(y*zi)/sum(y)
    pnew[2] = sum(y*i*zi)/sum(y*zi)
    pnew[3] = sum(y*i*(1-zi))/sum(y*(1-zi))
    p = pnew
    ## print(pnew)
    return(pnew)
}


objfn = function(p, y) {
    i = 0:(length(y)-1)
    loglik = y*log(p[1]*exp(-p[2])*p[2]^i/exp(lgamma(i+1)) +
        (1 - p[1])*exp(-p[3])*p[3]^i/exp(lgamma(i+1)))
    return ( -sum(loglik) )
}


res = turboem(par=c(0.5, 1, 3), fixptfn=fixptfn, objfn=objfn,
               method=c("em", "squarem", "pem"), y=y)
options(digits=13)
res





lmImpute = function(formula, data, n.iter = 100, tol = 1e-10){
    response = as.character(formula[[2]])
    missIndex = which(is.na(data[, response]))
    ## Initialize missing values
    data[, response] = na.approx(data[, response], na.rm = FALSE)
    init = lm(formula = formula, data = data)
    data[is.na(data[, response]), response] =
        predict(init, data[is.na(data[, response]), ])
    ## plot(1:52, data[, response])
    ll = -Inf
    for(i in 1:n.iter){
        print(i)
        fit = lm(formula = formula, data = data)
        print(coef(fit))
        ## lines(1:52, fitted(fit), col = "blue")
        l1 = logLik(fit)
        if((l1 - ll) < tol){
            break
        } else {
            data[missIndex, response] = fitted(fit)[missIndex]
            ll = l1
        }
        plot(1:NROW(data), data[, response])
        lines(fitted(fit))
    }
    list(impute = data[, response], fitted = fitted(fit))
}

## This is wrong
turboLmImpute = function(formula, data, n.iter = 100, tol = 1e-3){
    response = as.character(formula[[2]])
    missIndex = which(is.na(data[, response]))
    ## Initialize missing values
    data[, response] = na.approx(data[, response], na.rm = FALSE)
    init = lm(formula = formula, data = data)
    data[is.na(data[, response]), response] =
        predict(init, data[is.na(data[, response]), ])
    ## plot(1:52, data[, response])
    ll = -Inf
    fixptfn = function(p, y)
        data[, response] = p
        fitted(lm(formula = formula, data = data))

    objfn = function(p, y) sum(abs(p - y))
    res = turboem(par = data[missIndex, response], fixptfn = fixptfn,
        objfn = objfn, method = c("em", "squarem", "pem"),
        y = data[missIndex, response])
    res
}

test2 = turboLmImpute(formula = wheatYield_231 ~ ., data = tmp)


## Download data
library(FAOSTAT)
library(plyr)
library(zoo)
library(nnls)
library(swsImputation)

wheatYield.df = getFAOtoSYB(name = "wheatYield", domainCode = "QC",
    itemCode = 15, elementCode = 5419)$entity

## Select data which
cwheatYield.df = dcast(wheatYield.df, Year ~ FAOST_CODE,
    value.var = "wheatYield")
colnames(cwheatYield.df)[-1] =
    paste("wheatYield", colnames(cwheatYield.df)[-1], sep = "_")

## For simplicity we only work with countries that has full data
fwheatYield.df =
    cwheatYield.df[, which(colSums(is.na(cwheatYield.df)) == 0)]

tmp = fwheatYield.df[, unique(c("wheatYield_231",
    sample(colnames(fwheatYield.df), 10)))]

missProp = 0.8
missing = sample(1:NROW(fwheatYield.df),
    size = round(missProp * NROW(tmp)))
tmp[missing, "wheatYield_231"] = NA

system.time({
    test1 = lmImpute(formula = wheatYield_231 ~ ., data = tmp)
    })

system.time({
    test2 = turboLmImpute(formula = wheatYield_231 ~ ., data = tmp)
    })

tmp = fwheatYield.df[, unique(c("wheatYield_231",
    sample(colnames(fwheatYield.df), 10)))]
missProp = 0.7
missing = sample(1:NROW(fwheatYield.df),
    size = round(missProp * NROW(tmp)))
tmp[missing, "wheatYield_231"] = NA
test1 = lmImpute(formula = wheatYield_231 ~ ., data = tmp)
noEM.fit = lm(wheatYield_231 ~., data = tmp)
noEM.pred = predict(noEM.fit, newdata = tmp)
with(tmp, plot(1:NROW(tmp), wheatYield_231, pch = 19, cex = 2))
lines(1:52, test1$fitted, col = "steelblue", pch = 19)
lines(1:52, noEM.pred, col = "red", pch = 19)


########################################################################
## Title: package mi
## Date: 2014-02-13
########################################################################

library(mi)
data(CHAIN)

## Analyze the missing pattern
missing.pattern.plot(CHAIN, gray.scale = TRUE)
mp.plot(CHAIN, y.order = TRUE, x.order = TRUE, gray.scale = TRUE)

## Create imputation information matrix
info = mi.info(CHAIN)

## Preprocess data
IMP = mi.preprocess(CHAIN, info = info)
attr(IMP, "mi.info.preprocessed")

## Get the formula
info$imp.formula

## carry out the imputation
IMP = mi(CHAIN, info = info, check.coef.convergence = TRUE, n.iter = 50)

## Check the convergence, doesn't seem to work
plot(IMP@bugs)
plot(IMP@mcmc[, , 1])

## pooling estimates
fit = lm.mi(h39b.W1 ~ age.W1 + c28.W1 + pcs.W1 + mcs37.W1 +
             + b05.W1 + haartadhere.W1, IMP)
display(fit)



########################################################################
## Title: Forecast test
## Date: 2014-03-01
########################################################################

library(forecast)




########################################################################
## Title: Estimating growth rate in the presence of missing values
## Date: 2014-03-04
########################################################################

library(mi)

## Load the data, the production domain of Romania
romaniaMaize.df = structure(list(
    year = 1961:2011,
    production = c(5739600, 4932400, 6022700, 6691700,
        5877000, 8021600, 6857879, 7105287, 7675808, 6535512,
        7850300, 9816700, 7397200, 7439600, 9240661,
        11583200, 10113559, 9671200, 11768100, 10563300,
        8321500, 10550000, 10296000, 9891000, 11903200,
        10900900, 7526900, 7182200, 6761800, 6809604,
        10497338, 6828270, 7987450, 9343000, 9923132,
        9607944, 12686700, 8623370, 10934800, 4897600,
        9119200, 8399800, 9576985, 14541564, 10388499,
        8984729, 3853920, 7849080, 7973258, 9042032,
        11717591),
    areaHarvested = c(3428400, 3106766, 3379300, 3319100,
        3305800, 3287700, 3221075, 3344013, 3293068,
        3084036, 3131400, 3196500, 2956800, 2963000,
        3304691, 3377600, 3317671, 3178798, 3311291,
        3287560, 3077900, 2764387, 2935120, 3090530,
        3090100, 2858100, 2787100, 2579400, 2733400, 2466735,
        2574999, 3335920, 3065682, 2983400, 3109236, 3277041,
        3037700, 3128900, 3013400, 3049400, 2974000, 2761223,
        3119104, 3196130, 2609110, 2512944, 2263080, 2432210,
        2333501, 2094249, 2587102),
    seed = c(124000, 135000, 133000, 132000, 139000, 137000,
        140000, 140000, 135000, 135000, 138000, 135000,
        133000, 140000, 120000, 107000, 134000, 139000,
        133000, 135000, 124000, 144000, 148000, 145000,
        135000, 130000, 122000, 130000, 139000, 145000,
        180000, 180000, 142500, 125600, 130900, 90587,
        124753, 77000, 76000, 69907, 70692, 72785, 79065,
        67000, 64600, 63829, 65123, 62851, 59985, 53880,
        64744)),
    .Names = c("year", "production", "areaHarvested", "seed"),
    row.names = 6110:6160, class = "data.frame")

## Simulate missing value in production
n = NROW(romaniaMaize.df)
pctMiss = 0.5
myseed = 587
set.seed(myseed)
missIndex = sample(n, n * pctMiss)
romaniaMaize.df$productionSim = romaniaMaize.df$production
romaniaMaize.df[missIndex, "productionSim"] = NA

jpeg(file = "~/tmp/mi_illustration_1.jpg", width = 720)
with(romaniaMaize.df, plot(year, production, type = "b",
                           col = ifelse(is.na(productionSim),
                               "red", "black"), pch = 19,
                           ylim = c(0, max(production))))
legend("bottomleft", legend = c("observed", "simulated missing"),
       col = c("black", "red"), pch = 19, bty = "n")
graphics.off()

## Create missing information matrix
info = mi.info(romaniaMaize.df)
info = mi.info.update.include(info,
    list(production = FALSE, year = TRUE))

## I assume that the production has a trend and is determined by the
## area harvested and seed utilized.
info = mi.info.update.imp.formula(info,
    list(productionSim = "productionSim ~ year + areaHarvested + seed"))


## multiple imputation via mi
romaniaMaize.mi = mi(romaniaMaize.df, info = info, n.imp = 50,
    n.iter = 100, seed = myseed, run.past.convergence = TRUE)

## Compute the single imputation
romaniaMaize.df$imputed = romaniaMaize.df$productionSim
romaniaMaize.df[is.na(romaniaMaize.df$imputed), "imputed"] =
    rowMeans(sapply(romaniaMaize.mi@imp,
                    FUN = function(chain) chain$productionSim@random))

## Plot and compare the imputations
jpeg(file = "~/tmp/mi_illustration_2.jpg", width = 720)
with(romaniaMaize.df, plot(year, production, type = "n",
                          ylim = c(0, max(production, imputed))))
sapply(mi.completed(romaniaMaize.mi), FUN = function(chain)
       lines(romaniaMaize.df$year, chain$productionSim,
             col = rgb(0.5, 0.5, 0.5, alpha = 0.3)))
with(romaniaMaize.df, lines(year, production,
                          ylim = c(0, max(production, imputed))))
with(romaniaMaize.df, points(year, production, pch = 19,
                          ylim = c(0, max(production, imputed)),
                          col = is.na(productionSim) + 1))
with(romaniaMaize.df, lines(year, imputed, col = "steelblue", lwd = 2))
legend("bottomleft", legend = c("observed", "simulated missing",
                         "single imputation", "multiple imputation"),
       col = c("black", "red", "steelblue", "grey"), pch = 19, bty = "n")
graphics.off()

## Fit the regression on single/multiple imputed data set
single.fit = lm(log(romaniaMaize.df$imputed) ~ year,
    data = romaniaMaize.df)
multiple.fit = lm.mi(log(productionSim) ~ year, romaniaMaize.mi)

## Compare the growth rate
(exp(coef(single.fit)[2]) - 1) * 100
((exp(coef(multiple.fit)[2]) - 1) * 100)

## Compare the P-value
2 * pt(single.fit$coef[2]/sqrt(diag(vcov(single.fit)))[2],
       df = single.fit$df.residual, lower.tail = FALSE)
2 * pt(multiple.fit@mi.pooled$coef[2]/multiple.fit@mi.pooled$se[2],
       df = single.fit$df.residual, lower.tail = FALSE)




########################################################################
## Title: model simulation on USA maize
## Date: 2014-03-04
########################################################################

library(forecast)

usaMaize = c(91388000, 91604000, 102093008, 88504000, 104216928,
105861408, 123458304, 113022816, 119055936, 105471120, 143420656,
141733312, 144041760, 119420320, 148361072, 159751184, 165234544,
184613008, 201383008, 168647008, 206222000, 209180000, 106030000,
194880000, 225453008, 208943008, 181142000, 125194000, 191319008,
201532000, 189866496, 240719008, 160984992, 255292992, 187968992,
234527008, 233867008, 247882000, 239548992, 251852210, 241375035,
227765357, 256227304, 299873563, 282260662, 267501056, 331175072,
307142010, 332548610, 316164930, 313948610)


for(i in 1:100){
    start = sample(length(usaMaize) - 1, 1)
    remain = length(usaMaize) - start
    end = start + sample(remain, 1)
    tmp = usaMaize[start:end]
    fit = try(auto.arima(tmp))
    if(!inherits(fit, "try-error"))
        print(coef(fit))
}



########################################################################
## Title: EM-algorithm for estimating AR(1)
## Date: 2014-03-11
########################################################################

emar = function(x, n.iter = 100, tol = 1e-5, ...){
    require(forecast)
    missInd = which(is.na(x))
    ## Provide initial value for z (missing value)
    x.initimp = na.locf(na.locf(na.approx(x, na.rm = FALSE),
        na.rm = FALSE), na.rm = FALSE, fromLast = TRUE)
    old.ll = -Inf
    for(i in 1:n.iter){
        cat("Iteration:", i, "\n")
        ## Update the parameter (M-step)
        x.fit = arima(x = x.initimp, ..., method = "ML")
        fit.ll = logLik(x.fit)
        if((fit.ll - old.ll) > tol){
            ## Update the z (E-step) - fitted value of arima is
            ## one-step ahead forecast.
            x.initimp[missInd] = fitted(x.fit)[missInd]
            old.ll = fit.ll
        } else {
            break
        }
    }
    list(x = x, imputed = x.initimp, ll = fit.ll, model = x.fit)
}

n = 100
pctMiss = 0.5
raw = arima.sim(model = list(ar = 0.8), n = n, sd = 2)
simMiss = sample(n, n * pctMiss)
sim = raw
sim[simMiss] = NA
sim.imp = emar(sim, order = c(1, 0, 0), include.mean = FALSE)
plot(raw, type = "l")
points(raw, col = is.na(sim) + 1, pch = 19)
points(sim.imp$imputed, col = ifelse(is.na(sim), "steelblue", "black"),
       pch = 19)
sim.imp$model
arima(raw, order = c(1, 0, 0), include.mean = FALSE)
arima(sim, order = c(1, 0, 0), include.mean = FALSE, method = "ML")


########################################################################
## Title: Adding transparency to color
## Date: 2014-03-12
## Source:
########################################################################

library(MASS)
library(RColorBrewer)
# Source the colorRampAlpha file

addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

# colorRampPaletteAlpha()
colorRampPaletteAlpha <- function(colors, n=32, interpolate='linear') {
  # Create the color ramp normally
  cr <- colorRampPalette(colors, interpolate=interpolate)(n)
  # Find the alpha channel
  a <- col2rgb(colors, alpha=T)[4,]
  # Interpolate
	if (interpolate=='linear') {
		l <- approx(a, n=n)
	} else {
		l <- spline(a, n=n)
	}
	l$y[l$y > 255] <- 255 # Clamp if spline is > 255
	cr <- addalpha(cr, l$y/255.0)
	return(cr)
}

# ----------
# scalars:
col1 <- "red"
col2 <- rgb(1,0,0)
addalpha(col2, 0.8)
addalpha(col2,0.8)

# scalar alpha with vector of colors:
col3 <- c("red", "green", "blue", "yellow")
addalpha(col3, 0.8)
plot(rnorm(1000), col=addalpha(brewer.pal(11,'RdYlGn'), 0.5), pch=16)

# alpha and colors vector:
alpha <- seq.int(0, 1, length.out=4)
addalpha(col3, alpha)

# Simple example
x <- seq.int(0, 2*pi, length=1000)
y <- sin(x)
plot(x, y, col=addalpha(rep("red", 1000), abs(sin(y))))

# with RColorBrewer
x <- seq.int(0, 1, length.out=100)
z <- outer(x,x)
c1 <- colorRampPalette(brewer.pal(11, 'Spectral'))(100)
c2 <- addalpha(c1,x)
par(mfrow=c(1,2))
image(x,x,z,col=c1)
image(x,x,z,col=c2)

# colorRampPaletteAlpha()
# Create normally distributed data
x <- rnorm(1000)
y <- rnorm(1000)
k <- kde2d(x,y,n=250)

# Sample colors with alpha channel
col1 <- addalpha("red", 0.5)
col2 <-"green"
col3 <-addalpha("blue", 0.2)
cols <- c(col1,col2,col3)

# colorRampPalette ditches the alpha channel
# colorRampPaletteAlpha does not
cr1 <- colorRampPalette(cols)(32)
cr2 <- colorRampPaletteAlpha(cols, 32)

par(mfrow=c(1,2))
plot(x, y, pch=16, cex=0.3)
image(k$x,k$y,k$z,col=cr1, add=T)
plot(x, y, pch=16, cex=0.3)
image(k$x,k$y,k$z,col=cr2, add=T)

# Linear vs. spline interpolation
cr1 <- colorRampPaletteAlpha(cols, 32, interpolate='linear') # default
cr2 <- colorRampPaletteAlpha(cols, 32, interpolate='spline')
plot(x, y, pch=16, cex=0.3)
image(k$x,k$y,k$z,col=cr1, add=T)
plot(x, y, pch=16, cex=0.3)
image(k$x,k$y,k$z,col=cr2, add=T)


########################################################################
## Title: Bayesian Search Theory
## Date: 2014-03-17
## Source: http://jacobsimmering.com/2014/03/14/BayesianSearching.html
########################################################################

library(ggplot2)
d <- data.frame(x = rep(seq(-30, 30), each = 61),
                y = rep(seq(-30, 30), times = 61))
d$PrP <- dnorm(d$x, 0, sqrt(100)) * dnorm(d$y, 0, sqrt(100))
ggplot(d, aes(x = x, y = y, z = PrP)) +
    geom_point(aes(alpha = PrP)) + stat_contour()


detectionPower <- function(x, y, dx = 10, dy = 5,
                           p0 = 0.975, d = 0.925) {
    x2 <- x - dx
    y2 <- y - dx
    r <- sqrt(x2^2 + y2^2)
    power <- p0 * d^r
}
d$PrD <- detectionPower(d$x, d$y)
ggplot(d, aes(x = x, y = y, z = PrD)) +
    geom_point(aes(alpha = PrD)) + stat_contour(binwidth = 0.1)

d$valueOfSearch <- d$PrP * d$PrD
ggplot(d, aes(x = x, y = y, z = d$valueOfSearch)) +
    geom_point(aes(alpha = d$valueOfSearch)) +
    stat_contour()

nd <- data.frame(x = rep(d$x, 3), y = rep(d$y, 3),
                 value = c(d$PrP, d$PrD, d$valueOfSearch),
                 metric = gl(3, nrow(d),
                     labels = c("PDF of Object", "Detection Prob",
                         "Value of Search")))
ggplot(nd, aes(x = x, y = y, z = value)) + stat_contour() +
    facet_grid(. ~ metric) +
    scale_x_continuous(limits = c(-30, 30)) +
    scale_y_continuous(limits = c(-30,  30))

bayesUpdate <- function(searched, p0, pD) {
    (p0 * (1 - searched * pD))/(1 - p0 + p0 * (1 - pD))
}
d$searched <- rank(-1 * d$valueOfSearch) <= 100
d$newSearchValue <- bayesUpdate(d$searched, d$PrP, d$PrD)
nd <- data.frame(x = rep(d$x, 2), y = rep(d$y, 2),
                 valueOfSearch = c(d$valueOfSearch, d$newSearchValue),
                 searched = rep(d$searched, 2),
                 search = rep(c("Before Any Searching",
                     "First Wave"), each = nrow(d)))
nd$searched[nd$searched == FALSE] <- NA
ggplot(nd, aes(x = x, y = y, z = valueOfSearch)) + stat_contour() +
    facet_grid(. ~  search) +
    geom_point(aes(color = searched, alpha = valueOfSearch))



searchCount <- rep(0, nrow(d))
probInSearchArea <- numeric(1000)
probFindingInGrid <- numeric(1000)
p0 <- d$PrP * d$PrD
pD <- d$PrD
for (i in 1:100) {
    searchLocations <- rank(-1 * p0) <= 100
    searchCount <- searchCount + searchLocations
    probInSearchArea[i] <- sum(p0[searchLocations])
    probFindingInGrid[i] <- sum(p0)
    p0 <- bayesUpdate(searchLocations, p0, pD)
    nSearches <- data.frame(x = d$x, y = d$y, count = searchCount)
    print(ggplot(nSearches, aes(x = x, y = y, z = count)) +
          stat_contour() +
          geom_point(aes(alpha = count)))

}

nSearches <- data.frame(x = d$x, y = d$y, count = searchCount)
ggplot(nSearches, aes(x = x, y = y, z = count)) + stat_contour() +
    geom_point(aes(alpha = count))

############################################################
## Title: Using Genetic Algorithms in Quantitative Trading
## Date: 2014-03-19
## Source: https://gist.github.com/thertrader/9509411
############################################################
library(PerformanceAnalytics)
library(rgenoud)
library(quantmod)
library(TTR)


###############
outputPath <- getwd()
theInstrument <- "SPY"

data <- getSymbols(Symbols = theInstrument,
                   src = "yahoo",
                   from = "2000-01-01",
                   auto.assign = FALSE)

colnames(data) <- c("open","high","low","close","volume","adj.")


###############
fitnessFunction <- function(xx=c(1,1,1,1)){
  print(xx)
  rtn <- ROC(data[,"close"],n=1)
  rsi <- RSI(data[,"close"],n=xx[1],maType="SMA")
  smas <- SMA(data[,"close"],n=xx[3])
  smal <- SMA(data[,"close"],n=xx[4])
  aa <- cbind(data[,"close"],rtn,rsi,smas,smal)
  colnames(aa) <- c("close","rtn","rsi","smas","smal")

  isData <- aa[index(aa) < "2011-01-01"]

  posBuySignal <- which(isData[,"rsi"] <= (1 - xx[2]) & isData[,"smas"] > isData[,"smal"]) + 1
  if (length(posBuySignal) == 0)
    posBuySignal <- NULL
  posSellSignal <- which(isData[,"rsi"] > xx[2] & isData[,"smas"] < isData[,"smal"]) + 1
  if (length(posSellSignal) == 0)
    posSellSignal <- NULL
  allSignals <- c(posBuySignal,posSellSignal)
  allSignals <- allSignals[which(allSignals <= nrow(isData))]

  if (!is.null(allSignals) && length(allSignals) >= 50)
    theStat <- SharpeRatio.annualized(isData[sort(allSignals),"rtn"])
  if (is.null(allSignals) | length(allSignals) < 50)
    theStat <- 0

  return(theStat)
}


###############
tradingStatistics <- function(isOrOos = TRUE, xx = c(1,1,1,1)){
  print(xx)
  rtn <- ROC(data[,"close"],n=1)
  rsi <- RSI(data[,"close"],n=xx[1],maType="SMA")
  smas <- SMA(data[,"close"],n=xx[3])
  smal <- SMA(data[,"close"],n=xx[4])
  aa <- cbind(data[,"close"],rtn,rsi,smas,smal)
  colnames(aa) <- c("close","rtn","rsi","smas","smal")

  if (isOrOos == TRUE)
    sampleData <- aa[index(aa) < "2011-01-01"]
  if (isOrOos == FALSE)
    sampleData <- aa[index(aa) >= "2011-01-01"]

  posBuySignal <- which(sampleData[,"rsi"] <= (1 - xx[2]) & sampleData[,"smas"] > sampleData[,"smal"]) + 1
  if (length(posBuySignal) == 0)
    posBuySignal <- NULL
  posSellSignal <- which(sampleData[,"rsi"] > xx[2] & sampleData[,"smas"] < sampleData[,"smal"]) + 1
  if (length(posSellSignal) == 0)
    posSellSignal <- NULL
  allSignals <- c(posBuySignal,posSellSignal)
  allSignals <- allSignals[which(allSignals <= nrow(sampleData))]

  totalRtn <- sum(sampleData[sort(allSignals),"rtn"])
  numberOfTrades <- length(sampleData[sort(allSignals),"rtn"])
  hitRatio <- length(which(sampleData[sort(allSignals),"rtn"] > 0))/numberOfTrades

  return(list(totalRtn=totalRtn,numberOfTrades=numberOfTrades,hitRatio=hitRatio))
}


###########
optimum <- genoud(fitnessFunction,
                  nvars = 4,
                  max = TRUE,
                  pop.size = 30,
                  max.generations = 50,
                  wait.generations = 15,
                  hard.generation.limit = TRUE,
                  starting.values = c(5,70,30,100),
                  MemoryMatrix = TRUE,
                  Domains = matrix(c(5,50,10,50,
                                     50,90,50,200),
                                   nrow=4,ncol=2),
                  default.domains = 4,
                  solution.tolerance = 0.00001,
                  gr = NULL,
                  boundary.enforcement = 2,
                  lexical = FALSE,
                  gradient.check = FALSE,
                  data.type.int = TRUE,
                  hessian = FALSE,
                  unif.seed = 812821,
                  int.seed = 53058,
                  print.level = 2,
                  share.type = 0,
                  instance.number = 0,
                  output.path = paste(outputPath,theInstrument,"_Summary.txt",sep=""),
                  output.append = FALSE,
                  project.path = paste(outputPath,theInstrument,"_Details.txt",sep=""),
                  P1=2, P2=8, P3=8, P4=6, P5=6, P6=6, P7=8, P8=6, P9=0,
                  P9mix = NULL,
                  BFGSburnin = 0,
                  BFGSfn = NULL,
                  BFGShelp = NULL,
                  cluster = FALSE,
                  balance = FALSE,
                  debug = FALSE,
                  control = list())

solution <- optimum$par



########################################################################
## Title: Example from Jason Pushon
## Date: 2014-03-25
########################################################################


# Generate data
set.seed(1234)
n = 10e3

male.prob = 0.10
male.avgclaim = 800
female.prob = 0.075
female.avgclaim = 600

# male probability of having a claim = 0.05 & average claim = $1000
male.claims = rbinom(n, 1, male.prob)*rgamma(n, 3, 1/male.avgclaim)

# female probability of having a claim = 0.01 & average claim = $1000
female.claims = rbinom(n, 1, female.prob)*rgamma(n, 3, 1/female.avgclaim)

#   premiums =   (prob of claim x average claim)
male.prem = rgamma(n, 5, 1/(male.prob*male.avgclaim))
male.prem[male.prem<(male.prob*male.avgclaim)] <-
    (male.prob*male.avgclaim)

female.prem = rgamma(n, 5, 1/(female.prob*female.avgclaim))
female.prem[female.prem<(female.prob*female.avgclaim)] <-
    (female.prob*female.avgclaim)

male.data = cbind(gender='male',claims=male.claims, premium=male.prem)
female.data = cbind(gender='female', claims=female.claims,
    premium=female.prem)

data = as.data.frame(rbind(male.data, female.data))
data$claims <- as.numeric(as.character(data$claims))
data$premium <- as.numeric(as.character(data$premium))
data$lossratio <- data$claims/data$premium
mean(data$lossratio)
sum(data$claims)/sum(data$premium)



# Average loss ratio for males (policy weighted)
mean(data[data$gender =='male',]$lossratio)

## NOTE(Michael): weighted average loss ratio for male
with(data[data$gender == "male", ],
     weighted.mean(x = lossratio, w = premium))

# Aggregated loss ratio for males
sum(data[data$gender =='male',]$claims)/
    sum(data[data$gender =='male',]$premium)

# Average loss ratio for females (policy weighted)
mean(data[data$gender =='female',]$lossratio)

## NOTE(Michael): weighted average loss ratio for female
with(data[data$gender == "female", ],
     weighted.mean(x = lossratio, w = premium))

# Aggregated loss ratio for females
sum(data[data$gender =='female',]$claims)/
    sum(data[data$gender =='female',]$premium)

par(mfrow=c(2,2))
hist(data$premium)
hist(data[data$claims > 0,]$claims)
hist(data$lossratio)
hist(data[data$claims > 0,]$lossratio)

xyplot(claims ~ premium|gender, data = data)

data$claimed = as.numeric(data$claims != 0)
histogram(~lossratio|gender, data = data)

xyplot(math ~ socst|female + ses, data = ml,
       panel = function(x, y){
           panel.xyplot(x, y)
           panel.lmline(x, y)
       }
   )



########################################################################
## Title: Testing turboEM part 2
## Date: 2014-03-26
########################################################################

library(turboEM)
poissmix.dat = data.frame(death = 0:9,
    freq = c(162,267,271,185,111,61,27,8,3,1))
y = poissmix.dat$freq


## Function to update parameter given data
fixptfn = function(p, y) {
    pnew = rep(NA,3)
    i = 0:(length(y)-1)
    denom = p[1]*exp(-p[2])*p[2]^i + (1 - p[1])*exp(-p[3])*p[3]^i
    zi = p[1]*exp(-p[2])*p[2]^i / denom
    pnew[1] = sum(y*zi)/sum(y)
    pnew[2] = sum(y*i*zi)/sum(y*zi)
    pnew[3] = sum(y*i*(1-zi))/sum(y*(1-zi))
    p = pnew
    print(pnew)
    return(pnew)
}

## Function to compute likelihood
objfn = function(p, y) {
    i = 0:(length(y)-1)
    loglik = y*log(p[1]*exp(-p[2])*p[2]^i/exp(lgamma(i+1)) +
        (1 - p[1])*exp(-p[3])*p[3]^i/exp(lgamma(i+1)))
    return ( -sum(loglik) )
}


res = turboem(par=c(0.8, 1, 3), fixptfn=fixptfn, objfn=objfn,
               method=c("em", "squarem", "pem"), y=y)
options(digits=13)
res

res1 = turboem(par=c(0.5, 1, 3), fixptfn=fixptfn, objfn=objfn,
               method=c("em", "squarem", "pem"), y = y,
    control.run = list(keep.objfval = TRUE))
plot(res1)



llar1 = function(x, rho, sigma2){
    T = length(x)
    -(T - 1)/2 * log(2 * pi) - (T - 1)/2 * log(sigma2) -
        sum((x[-1] - rho * x[-T])^2)/(2 * sigma2)
}

y = arima.sim(model = list(ar = 0.8), n = 100)
y.fit = arima(y, order = c(1, 0, 0), include.mean = FALSE)
logLik(y.fit)
llar1(x = y, rho = coef(y.fit)["ar1"], sigma2 = y.fit$sigma2)

y = arima.sim(model = list(ar = 0.8), n = 100)
y.fit = arima(y, order = c(1, 0, 0), include.mean = FALSE)
predict(y.fit, n.ahead = 10)$pred
forecast(y.fit, h = 10)$mean
plot(forecast(y.fit, h = 10))



########################################################################
## Title: Demonstration for imputation presentation
## Date: 2014-03-26
########################################################################


## Example for multiple imputation
## ---------------------------------------------------------------------

emar = function(x, n.iter = 100, tol = 1e-5, ...){
    require(forecast)
    missInd = which(is.na(x))
    ## Provide initial value for z (missing value)
    x.initimp = na.locf(na.locf(na.approx(x, na.rm = FALSE),
        na.rm = FALSE), na.rm = FALSE, fromLast = TRUE)
    old.ll = -Inf
    for(i in 1:n.iter){
        cat("Iteration:", i, "\n")
        ## Update the parameter (M-step)
        x.fit = arima(x = x.initimp, ...)
        fit.ll = logLik(x.fit)
        if((fit.ll - old.ll) > tol){
            ## Update the z (E-step)
            x.initimp[missInd] = fitted(x.fit)[missInd]
            old.ll = fit.ll
        } else {
            break
        }
    }
    list(x = x, imputed = x.initimp, ll = fit.ll, model = x.fit)
}

## Simulate an ar(1) process
n = 100
pctMiss = 0.5
raw = arima.sim(model = list(ar = 0.8), n = n)

## Create missing value
simMiss = sample(n, n * pctMiss)
sim = raw
sim[simMiss] = NA

## Impute the time serie and estimate the ar(1) model
sim.imp = emar(sim, order = c(1, 0, 0), include.mean = FALSE)

## Plot the result and return the model
par(mfrow = c(3, 1))
plot(raw, type = "l", main = "raw data")
plot(sim, type = "b", main = "simulated")
plot(sim.imp$imputed, type = "l", main = "imputed",
ylab = "imputed")
sim.imp$model


## Example for multiple imputation
## ---------------------------------------------------------------------

library(mi)

## Load the data, the production domain of Romania
romaniaMaize.df = structure(list(
    year = 1961:2011,
    production = c(5739600, 4932400, 6022700, 6691700,
        5877000, 8021600, 6857879, 7105287, 7675808, 6535512,
        7850300, 9816700, 7397200, 7439600, 9240661,
        11583200, 10113559, 9671200, 11768100, 10563300,
        8321500, 10550000, 10296000, 9891000, 11903200,
        10900900, 7526900, 7182200, 6761800, 6809604,
        10497338, 6828270, 7987450, 9343000, 9923132,
        9607944, 12686700, 8623370, 10934800, 4897600,
        9119200, 8399800, 9576985, 14541564, 10388499,
        8984729, 3853920, 7849080, 7973258, 9042032,
        11717591),
    areaHarvested = c(3428400, 3106766, 3379300, 3319100,
        3305800, 3287700, 3221075, 3344013, 3293068,
        3084036, 3131400, 3196500, 2956800, 2963000,
        3304691, 3377600, 3317671, 3178798, 3311291,
        3287560, 3077900, 2764387, 2935120, 3090530,
        3090100, 2858100, 2787100, 2579400, 2733400, 2466735,
        2574999, 3335920, 3065682, 2983400, 3109236, 3277041,
        3037700, 3128900, 3013400, 3049400, 2974000, 2761223,
        3119104, 3196130, 2609110, 2512944, 2263080, 2432210,
        2333501, 2094249, 2587102),
    seed = c(124000, 135000, 133000, 132000, 139000, 137000,
        140000, 140000, 135000, 135000, 138000, 135000,
        133000, 140000, 120000, 107000, 134000, 139000,
        133000, 135000, 124000, 144000, 148000, 145000,
        135000, 130000, 122000, 130000, 139000, 145000,
        180000, 180000, 142500, 125600, 130900, 90587,
        124753, 77000, 76000, 69907, 70692, 72785, 79065,
        67000, 64600, 63829, 65123, 62851, 59985, 53880,
        64744)),
    .Names = c("year", "production", "areaHarvested", "seed"),
    row.names = 6110:6160, class = "data.frame")

## Simulate missing value in production
n = NROW(romaniaMaize.df)
pctMiss = 0.5
myseed = 587
set.seed(myseed)
missIndex = sample(n, n * pctMiss)
romaniaMaize.df$productionSim = romaniaMaize.df$production
romaniaMaize.df[missIndex, "productionSim"] = NA

with(romaniaMaize.df, plot(year, production, type = "b",
                           col = ifelse(is.na(productionSim),
                               "red", "black"), pch = 19,
                           ylim = c(0, max(production))))
legend("bottomleft", legend = c("observed", "simulated missing"),
       col = c("black", "red"), pch = 19, bty = "n")

## Create missing information matrix
info = mi.info(romaniaMaize.df)
info = mi.info.update.include(info,
    list(production = FALSE, year = TRUE))

## I assume that the production has a trend and is determined by the
## area harvested and seed utilized.
info = mi.info.update.imp.formula(info,
    list(productionSim = "productionSim ~ year + areaHarvested + seed"))


## multiple imputation via mi
romaniaMaize.mi = mi(romaniaMaize.df, info = info, n.imp = 50,
    n.iter = 100, seed = myseed, run.past.convergence = TRUE)

## Compute the single imputation
romaniaMaize.df$imputed = romaniaMaize.df$productionSim
romaniaMaize.df[is.na(romaniaMaize.df$imputed), "imputed"] =
    rowMeans(sapply(romaniaMaize.mi@imp,
                    FUN = function(chain) chain$productionSim@random))

## Plot and compare the imputations
with(romaniaMaize.df, plot(year, production, type = "n",
                          ylim = c(0, max(production, imputed))))
sapply(mi.completed(romaniaMaize.mi), FUN = function(chain)
       lines(romaniaMaize.df$year, chain$productionSim,
             col = rgb(0.5, 0.5, 0.5, alpha = 0.3)))
with(romaniaMaize.df, lines(year, production,
                          ylim = c(0, max(production, imputed))))
with(romaniaMaize.df, points(year, production, pch = 19,
                          ylim = c(0, max(production, imputed)),
                          col = is.na(productionSim) + 1))
with(romaniaMaize.df, lines(year, imputed, col = "steelblue", lwd = 2))
legend("bottomleft", legend = c("observed", "simulated missing",
                         "single imputation", "multiple imputation"),
       col = c("black", "red", "steelblue", "grey"), pch = 19, bty = "n")

## Fit the regression on single/multiple imputed data set
single.fit = lm(log(romaniaMaize.df$imputed) ~ year,
    data = romaniaMaize.df)
multiple.fit = lm.mi(log(productionSim) ~ year, romaniaMaize.mi)

## Compare the growth rate
(exp(coef(single.fit)[2]) - 1) * 100
((exp(coef(multiple.fit)[2]) - 1) * 100)

## Compare the P-value
2 * pt(single.fit$coef[2]/sqrt(diag(vcov(single.fit)))[2],
       df = single.fit$df.residual, lower.tail = FALSE)
2 * pt(multiple.fit@mi.pooled$coef[2]/multiple.fit@mi.pooled$se[2],
       df = single.fit$df.residual, lower.tail = FALSE)




########################################################################
## Title: Loess test
## Date: 2014-03-27
########################################################################

library(zoo)
x = c(arima.sim(100, model = list(ar = 0.8))) + 100
y = abs(rnorm(100, sd = 0.1)) * 1:100 + x
plot(x, ylim = c(range(c(x, y))), type = "l")
lines(y)

y.sim = y
y.sim[25:50] = NA
y.imp = na.approx(y.sim)
plot(x, ylim = c(range(c(x, y))), type = "l")
lines(y.imp)
plot(x, y)

plot(x, ylim = c(range(c(x, y))), type = "l")
lines(y.imp)
lines(y, col = "steelblue")
z = 1:100
lines(fitted(loess(y ~ x + z)), col = "red")
lines(fitted(loess(y ~ x + z, amily = "symmetric")), col = "green")





########################################################################
## Title: Wavelet test
## Date: 2013-04-06
########################################################################

## obtain the two series listed in Percival and Walden (2000), page 42
X1 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,.7,.9,0,.3)
X2 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,-.7,.9,0,.3)


## combine them and compute DWT
newX <- cbind(X1,X2)
wt <- dwt(newX, n.levels=3, boundary="reflection", fast=FALSE)

library(wmtsa)
test = c(arima.sim(model = list(ar = 0.8), n = 512))
result = wavDWT(test, wavelet="s8", n.levels=5, keep.series=TRUE)

wmtsa:::plot.wavTransform
myWav(result)

## plot the transform shifted for approximate zero
## phase alignment
plot(wavShift(result))

## plot summary
eda.plot(result)

## summarize the transform
summary(result)


testReconstruct = reconstruct(result)
plot(test, type = "l")
lines(testReconstruct, col = "steelblue")





########################################################################
## Title: data.table tricks
## Date: 2014-04-08
########################################################################
library(data.table)
n = 10000
DT = data.table(grp = 1:n,name = as.character(as.hexmode(1:n)),
    x = rnorm(10*n),y = rnorm(10*n))
setkey(DT,grp)
system.time(ans1<-DT[,lapply(.SD[,list(x,y)],sum),by=grp])
#   user  system elapsed
# 31.130   0.088  31.288    # bad
system.time(ans2<-DT[,lapply(list(x,y),sum),by=grp])
#   user  system elapsed
#  0.284   0.004   0.291    # good
setnames(ans2,names(ans1))
identical(ans1,ans2)
# [1] TRUE
system.time(ans3<-DT[,lapply(.SD,sum),by=grp,.SDcols=c("x","y")])
#   user  system elapsed
#  0.080   0.004   0.085    # even better (prior to v1.8.2 was slower and not recommended, no longer)
identical(ans1,ans3)
# [1] TRUE
tables()

library(data.table)
n=10000
DT = data.table(grp=1:n,x=rnorm(10*n),y=rnorm(10*n))
setkey(DT,grp)
system.time(ans1 <- DT[,transform(.SD,x2=x/sum(x),y2=y/sum(y)),by=grp])
# user  system elapsed
# 5.46    0.00    5.48  # slow
head(ans1,3)
#    grp          x           y        x2         y2
# 1:   1 -0.5848814 -0.41560829 0.6241268  0.5695575
# 2:   1 -0.6314059 -0.49076645 0.6737731  0.6725557
# 3:   1 -1.7694071  0.08860505 1.8881340 -0.1214260
system.time(tt <- DT[,list(x2=x/sum(x),y2=y/sum(y)),by=grp])
# user  system elapsed
# 0.02    0.00    0.02 (274 times faster!!!)

head(tt,3)
#    grp        x2         y2
# 1:   1 0.6241268  0.5695575
# 2:   1 0.6737731  0.6725557
# 3:   1 1.8881340 -0.1214260
system.time(ans2 <- cbind(DT,tt[,list(x2,y2)]))
# user  system elapsed
# 0.05    0.00    0.05    # very fast to add afterwards in bulk
head(ans2,3)
#     grp          x           y        x2         y2
# 1:   1 -0.5848814 -0.41560829 0.6241268  0.5695575
# 2:   1 -0.6314059 -0.49076645 0.6737731  0.6725557
# 3:   1 -1.7694071  0.08860505 1.8881340 -0.1214260
setkey(ans2,grp)
identical(ans1,ans2)
[1] TRUE
system.time(DT[, c('x2', 'y2') := list(x / sum(x), y / sum(y)), by = grp])
# user  system elapsed
# 0.07    0.00    0.07 # equivalent to cbind afterwards approach, but more memory efficient

# now DT has been updated
identical(ans1, DT)
# [1] TRUE
# remove new columns to show different approach
DT[, c('x2', 'y2') := NULL]

system.time(DT[, `:=`(x2=x / sum(x),y2= y / sum(y)), by = grp])
# user  system elapsed
#0.04    0.00    0.05 # this is slightly faster
identical(ans1, DT)

n=10000
DT = data.table(grp=1:n,x=rnorm(10*n),y=rnorm(10*n))
setkey(DT,grp)
system.time(ans1<-DT[,as.list(colSums(.SD)),by=grp])
system.time(ans2<-DT[,lapply(.SD,sum),by=grp])
identical(ans1,ans2)
system.time(ans3<-DT[,list(x=sum(x),y=sum(y)),by=grp])
identical(ans1,ans3)





########################################################################
## Title: examine the degree of freedom of splines
## Date: 2014-04-08
########################################################################

library(splines)
library(MASS)
plot(mcycle)
lines(predict(lm(accel ~ bs(times, df = 5), data = mcycle),
              newdata = mcycle), col = "red")
lines(predict(lm(accel ~ bs(times, df = 15), data = mcycle),
              newdata = mcycle), col = "green")
lines(predict(lm(accel ~ bs(times), data = mcycle),
              newdata = mcycle), col = "steelblue")




summary(lm(accel ~ bs(accel), data = mcycle))

tmp = bs(mcycle$accel)
plot(data.frame(mcycle$accel, tmp))



########################################################################
## Title: The Collatz Fractal
## Date: 2014-04-08
## Source: http://aschinchon.wordpress.com/2014/04/04/the-collatz-fractal/
########################################################################

library(ggplot2)
xrange <- seq(-8, 8, by = 0.01)
yrange <- seq(-3, 3, by = 0.01)
f  <- function (z) {1/4*(2+7*z-(2+5*z)*cos(pi*z))}
z <- outer(xrange, 1i*yrange,'+')
t <- mat.or.vec(nrow(z), ncol(z))
for (k in 1:10)
{
  z <- f(z)
  t <- t + (is.finite(z)+0)
}

## Supressing texts, titles, ticks, background and legend.
opt <- theme(legend.position="none",
             panel.background = element_blank(),
             axis.ticks=element_blank(),
             axis.title=element_blank(),
             axis.text =element_blank())
z = data.frame(expand.grid(x=xrange, y=yrange), z=as.vector(t))
ggplot(z, aes(x=x, y=y, color=z)) +
    geom_tile() + scale_colour_gradient(low="red", high="yellow") + opt


########################################################################
## Title: Extrapolation based on loess
## Date: 2014-04-09
########################################################################

library(MASS)
data(UKgas)
tmp = c(log(UKgas))
plot(tmp, type = "l")
tmp[1:30] = NA
T = 1:length(tmp)
tmp.fit = loess(tmp ~ T, control = loess.control(surface = "direct"))
tmp.predict = predict(tmp.fit, newdata = data.frame(T = T))
plot(c(log(UKgas)), type = "l", ylim = c(0, 8))
points(1:108, c(log(UKgas)),col = c(rep("red", 30), rep("black", 78)))
lines(tmp.predict)


########################################################################
## Title: Piecewise regression Date: 2014-04-09 Source:
## http://climateecology.wordpress.com/2012/08/19/r-for-ecologists-putting-together-a-piecewise-regression/
########################################################################


x <- c(1:10, 13:22)
y <- numeric(20)
## Create first segment
y[1:10] <- 20:11 + rnorm(10, 0, 1.5)
## Create second segment
y[11:20] <- seq(11, 15, len=10) + rnorm(10, 0, 1.5)
## Plot it
par(mar=c(4,4,1,1)+0.2)
plot(x,y, ylim=c(5, 20), pch=16)

breaks <- x[which(x >= 9 & x <= 17)]

mse <- numeric(length(breaks))
for(i in 1:length(breaks)){
 piecewise1 <- lm(y ~ x*(x < breaks[i]) + x*(x>=breaks[i]))
 mse[i] <- summary(piecewise1)[6]
}
mse <- as.numeric(mse)

breaks[which(mse==min(mse))]

piecewise2 <- lm(y ~ x*(x < 15) + x*(x > 15))
summary(piecewise2)

plot(x,y, ylim=c(5, 20), pch=16)
curve((3.3133 + 16.6352) + (0.5843-1.3025)*x, add=T, from=1, to=15)
curve((3.3133 - 0.9116) + 0.5843*x, add=T, from=15, to=max(x))
abline(v=15, lty=3)

library(segmented)
lin.mod <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=14)
plot(x,y, pch=16, ylim=c(5,20))
plot(segmented.mod, add = TRUE)



########################################################################
## Title: Examination of the FSI data
## Date: 2014-04-30
########################################################################

library(data.table)
library(FAOSTAT)
library(lattice)
load("C:/Users/kao/Dropbox/FSIproject/FSI2013/Inputs/final.Rdata")
pov.df = getWDI("SI.POV.DDAY", startDate = min(fsi.dt$Year),
    endDate = max(fsi.dt$Year))
final.df = mergeSYB(fsi.df, pov.df, all.x = TRUE)
final.df$Country = NULL
final.df = merge(final.df,
    FAOcountryProfile[!is.na(FAOcountryProfile$FAOST_CODE),
                      c("FAOST_CODE", "ABBR_FAO_NAME")],
    all.x = TRUE)
final.dt = data.table(final.df[final.df$FAOST_CODE %in%
    na.omit(FAOcountryProfile$FAOST_CODE), ])
setkeyv(final.dt, c("FAOST_CODE", "Year"))
final.dt[FAOST_CODE == 357, ABBR_FAO_NAME := "Taiwan and China"]
final.dt[FAOST_CODE == 107, ABBR_FAO_NAME := "Cote d'Ivoire"]
final.dt[FAOST_CODE == 284, ABBR_FAO_NAME := "Aland Islands"]
final.dt[FAOST_CODE == 279, ABBR_FAO_NAME := "Curacao"]
final.dt[FAOST_CODE == 182, ABBR_FAO_NAME := "Reunion"]
final.dt[FAOST_CODE == 282, ABBR_FAO_NAME := "Saint Barthelemy"]

xyplot(POU * 100 + SI.POV.DDAY ~ Year|ABBR_FAO_NAME,
       data = final.dt[Year %in% c(1985:2013), ], type = c("g", "l"))



xyplot(POU * 100 ~ SI.POV.DDAY|ABBR_FAO_NAME,
       data = final.dt[Year %in% c(1985:2013), ], type = c("g", "l"))


########################################################################
## Title: PISA visualization
## Date: 2014-05-29
########################################################################




library(ggplot2)
# read students data from PISA 2012
# directly from URL. This data is 210 mb
load("student2012.rda")


# calculate weighted mean from math / reading scores
# W_FSTUWT stands for final weights
mathScores <- unclass(by(student2012[,c("PV1MATH", "W_FSTUWT")],
                 student2012$CNT,
                 function(x) weighted.mean(x[,1], x[,2])) )
readScores <- unclass(by(student2012[,c("PV1READ", "W_FSTUWT")],
                 student2012$CNT,
                 function(x) weighted.mean(x[,1], x[,2])) )
# create a data.frame with scores and country names
# remove names with ( in name)
readMathScores <- data.frame(Country=names(readScores), readScores, mathScores)
readMathScores <- readMathScores[-grep(readMathScores$Country, pattern="(", fixed=TRUE),]

ggplot(readMathScores, aes(x=mathScores + readScores, y = mathScores - readScores, label = Country)) +
  geom_text() +
  theme_bw()

########################################################################
## Title: Test on the Numbeo api
## Date: 2014-06-05
########################################################################

link = "http://www.numbeo.com/api/"
key = "fao_1949"


Sys.setlocale("LC_ALL","C")
test = fromJSON("http://www.numbeo.com/api/country_prices?api_key=fao_1949&country=Italy")

cities = fromJSON("http://www.numbeo.com/api/cities?api_key=fao_1949")

sapply(cities$cities,
       function(x){
           tmp = try(x$coutntry)
           if(!inherits(tmp, "try-error")){
               if(tmp == "Italy")
                   try(print(x$city))
           }
       }
       )

########################################################################
## Title: box and beard plot
## Date: 2014-06-06
## Source: http://freakonometrics.hypotheses.org/14682
########################################################################

## Normal boxplot
set.seed(2)
x = rnorm(500)
boxplot(x,horizontal = TRUE,axes = FALSE)
axis(1)

## Add in beard
boxplot(x,horizontal = TRUE,xlim = c(-1,1.3),axes = FALSE)
axis(1)
Q = quantile(x,c(.25,.75))
y = cut(x[(x >= Q[1])&(x <= Q[2])],seq(Q[1],Q[2],length = 11))
tb = table(y)
u = seq(Q[1],Q[2],length = 11)
umid = (u[1:10]+u[2:11])/2
for(i in 1:10) segments(umid[i],1-.2,umid[i],1-.2-tb[i]/20,lwd = 3)

## Increase the band of the beard
rect(Q[1],1-.2,Q[2],1+.2)
segments(median(x),1-.2,median(x),1+.2,lwd = 2)
du = diff(umid)[1]
y = cut(x,seq(Q[1]-du*8,Q[2]+du*8,length = 11+16))
tb = table(y)
u = seq(Q[1]-du*8,Q[2]+du*8,length = 11+16)
umid = (u[1:26]+u[2:27])/2
for(i in 1:8) segments(umid[i],1,umid[i],1-.2-tb[i]/20,lwd = 3)
for(i in 19:26) segments(umid[i],1,umid[i],1-.2-tb[i]/20,lwd = 3)

## Add in smoother
vu = seq(-2.5,2.5,by = .02)
vv = dnorm(vu)
lines(vu,1-vv*4,col = "red",lty = 2)
d = density(x,bw = .1)
lines(d$x,1-d$y*4,col = "red")



########################################################################
## Title: Testing library ggivs
## Date: 2014-06-26
## Source: http://blog.revolutionanalytics.com/2014/06/interactive-web-ready-ggplot2-style-graphics-with-ggvis.html
########################################################################

library(ggvis)

mtcars %>%
  ggvis(~wt, ~mpg) %>%
  layer_smooths(span = input_slider(0.5, 1, value = 1)) %>%
  layer_points(size := input_slider(100, 1000, value = 100))



########################################################################
## Title: Checking the changes in arable land
## Date: 2014-06-03
########################################################################

library(FAOSTAT)
library(data.table)
library(lattice)
arableLand.df =
    getFAOtoSYB(name = "arableLand",
                domainCode = "RL",
                elementCode = 5110,
                itemCode = 6621)$entity

arableLand.dt = data.table(
    merge(arableLand.df,
          FAOcountryProfile[, c("FAOST_CODE", "FAO_TABLE_NAME")],
          by = "FAOST_CODE"
          )
    )
arableLand.dt = arableLand.dt[FAOST_CODE != 357, ]
setkeyv(arableLand.dt, c("FAOST_CODE", "Year"))

deltaArableLand.dt =
    arableLand.dt[, (exp(coef(lm(log(arableLand) ~ Year))[2]) - 1) * 100,
                  by = "FAO_TABLE_NAME"]
## Maximum 7% increase over 50 years, so very stable.

xyplot(sqrt(arableLand) ~ Year|FAO_TABLE_NAME,
       data = arableLand.dt[!FAOST_CODE %in% c(228, 185), ],
       type = c("g", "l"), ylim = c(0, 150))

xyplot(arableLand ~ Year|FAO_TABLE_NAME, data = arableLand.dt,
       scales = list(relation = "sliced", draw = FALSE),
       type = c("g", "l"))




########################################################################
## Title: simple bayesian regression
## Date: 2014-07-10
## source: http://stats.stackexchange.com/questions/47008/how-would-you-do-bayesian-anova-and-regression-in-r
########################################################################

library(bayesm)

podwt = structure(list(wt = c(1.76, 1.45, 1.03, 1.53, 2.34, 1.96,
1.79, 1.21, 0.49, 0.85, 1, 1.54, 1.01, 0.75, 2.11, 0.92), treat =
structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L,
2L, 2L), .Label = c("I", "U"), class = "factor"), mus = c(4.15, 2.76,
1.77, 3.11, 4.65, 3.46, 3.75, 2.04, 1.25, 2.39, 2.54, 3.41, 1.27,
1.26, 3.87, 1.01)), .Names = c("wt", "treat", "mus"), row.names =
c(NA, -16L), class = "data.frame")

# response
y1 <- podwt$wt

# First run a one-way anova

# Create the design matrix - need to insert a column of 1s
x1 <- cbind(matrix(1,nrow(podwt),1),podwt$treat)

# data for the Bayesian analysis
dt1 <- list(y=y1,X=x1)

# runiregGibbs uses a normal prior for the regression coefficients and
# an inverse chi-squared prior for va

# mean of the normal prior. We have 2 estimates - 1 intercept
# and 1 regression coefficient
betabar1 <- c(0,0)

# Pecision matrix for the normal prior. Again we have 2
A1 <- 0.01 * diag(2)
# note this is a very diffuse prior

# degrees of freedom for the inverse chi-square prior
n1 <- 3

# scale parameter for the inverse chi-square prior
ssq1 <- var(y1)

Prior1 <- list(betabar=betabar1, A=A1, nu=n1, ssq=ssq1)

# number of iterations of the Gibbs sampler
iter <- 10000

# thinning/slicing parameter. 1 means we keep all all values
slice <- 1

MCMC <- list(R=iter, keep=slice)

sim1 <- runiregGibbs(dt1, Prior1, MCMC)

plot(sim1$betadraw)
plot(sim1$sigmasqdraw)

summary(sim1$betadraw)
summary(sim1$sigmasqdraw)

# compare with maximum likelihood estimates:
fitpodwt <- lm(wt~treat, data=podwt)
summary(fitpodwt)
anova(fitpodwt)


# now for ordinary linear regression

x2 <- cbind(matrix(1,nrow(podwt),1),podwt$mus)

dt2 <- list(y=y1,X=x2)

sim2 <- runiregGibbs(dt1, Prior1, MCMC)

summary(sim1$betadraw)
    summary(sim1$sigmasqdraw)
plot(sim$betadraw)
    plot(sim$sigmasqdraw)

# compare with maximum likelihood estimates:
summary(lm(podwt$wt~mus,data=podwt))


# now with both variables

x3 <- cbind(matrix(1,nrow(podwt),1),podwt$treat,podwt$mus)

dt3 <- list(y=y1,X=x3)

# now we have an additional estimate so modify the prior accordingly

betabar1 <- c(0,0,0)
A1 <- 0.01 * diag(3)
Prior1 <- list(betabar=betabar1, A=A1, nu=n1, ssq=ssq1)

sim3 <- runiregGibbs(dt3, Prior1, MCMC)

plot(sim3$betadraw)
    plot(sim3$sigmasqdraw)
summary(sim3$betadraw)
    summary(sim3$sigmasqdraw)

# compare with maximum likelihood estimates:
summary(lm(podwt$wt~treat+mus,data=podwt))



########################################################################
## Title: Stan demo
## Date: 2014-07-15
## Source: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#prerequisites
########################################################################

library(rstan)
schools_code <- "
  data {
    int<lower=0> J; // number of schools
    real y[J]; // estimated treatment effects
    real<lower=0> sigma[J]; // s.e. of effect estimates
  }
  parameters {
    real mu;
    real<lower=0> tau;
    real eta[J];
  }
  transformed parameters {
    real theta[J];
    for (j in 1:J)
      theta[j] <- mu + tau * eta[j];
  }
  model {
    eta ~ normal(0, 1);
    y ~ normal(theta, sigma);
  }
"

schools_dat <- list(J = 8,
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit1 <- stan(model_code = schools_code, data = schools_dat,
            iter = 1000, chains = 4)

fit2 <- stan(fit = fit1, data = schools_dat, iter = 10000, chains = 4)
print(fit2)
plot(fit2)


########################################################################
## Title: Rice forecast
## Date: 2014-07-22
########################################################################

library(FAOSTAT)
library(data.table)
library(lattice)
library(lme4)

riceProduction.dt =
    data.table(
        getFAOtoSYB(domainCode = "QC",
                    itemCode = 27,
                    elementCode = 5510,
                    name = "riceProduction"
                    )$entity
        )

riceProduction.dt =
    merge(riceProduction.dt,
          FAOcountryProfile[, c("FAOST_CODE", "FAO_TABLE_NAME")],
          by = "FAOST_CODE"
          )

xyplot(log(riceProduction) ~ Year|FAO_TABLE_NAME,
       data = riceProduction.dt, type = "l")

xyplot(riceProduction ~ Year|FAO_TABLE_NAME,
       data = riceProduction.dt, type = "l")

## Remove Hong kong for production, since production does either not
## exist anymore or included in the statistic of China.
riceProduction.dt = riceProduction.dt[FAOST_CODE != 96, ]

aggregatedRiceProduction.dt =
    riceProduction.dt[, sum(riceProduction, na.rm = TRUE),
                      by = "Year"]

globalLinearPredict = predict(lm(V1 ~ Year,
    data = aggregatedRiceProduction.dt),
    newdata = data.frame(Year = 2013))

countryLinearPredict =
    sum(predict(lmer(riceProduction ~ Year|FAO_TABLE_NAME,
                 data = riceProduction.dt),
            newdata = data.frame(FAO_TABLE_NAME =
                unique(riceProduction.dt$FAO_TABLE_NAME),
                Year = 2013)))

finalLinearPredict = (globalLinearPredict + countryLinearPredict)/2

########################################################################
## Title: D3 network
## Date: 2014-08-07
## Source: http://christophergandrud.github.io/d3Network/
########################################################################

library(d3Network)

Source = c("A", "A", "A", "A", "B", "B", "C", "C", "D")
Target = c("B", "C", "D", "J", "E", "F", "G", "H", "I")
NetworkData <- data.frame(Source, Target)

d3SimpleNetwork(NetworkData, width = 400, height = 250,
                file = "simple.html")

# Load RCurl package for downloading the data
library(RCurl)

# Gather raw JSON formatted data
URL <- "https://raw.githubusercontent.com/christophergandrud/d3Network/master/JSONdata/miserables.json"
MisJson <- getURL(URL, ssl.verifypeer = FALSE)

# Convert JSON arrays into data frames
MisLinks <- JSONtoDF(jsonStr = MisJson, array = "links")
MisNodes <- JSONtoDF(jsonStr = MisJson, array = "nodes")

d3ForceNetwork(Links = MisLinks, Nodes = MisNodes,
               Source = "source", Target = "target",
               Value = "value", NodeID = "name",
               Group = "group", width = 800, height = 800,
               opacity = 0.9, file = "test.html")

URL <- "https://raw.githubusercontent.com/christophergandrud/d3Network/sankey/JSONdata/energy.json"
Energy <- getURL(URL, ssl.verifypeer = FALSE)
# Convert to data frame
EngLinks <- JSONtoDF(jsonStr = Energy, array = "links")
EngNodes <- JSONtoDF(jsonStr = Energy, array = "nodes")

# Plot
d3Sankey(Links = EngLinks, Nodes = EngNodes, Source = "source",
         Target = "target", Value = "value", NodeID = "name",
         fontsize = 12, nodeWidth = 30, width = 700,
         file = "sankeyTest.html")


URL <- "https://raw.githubusercontent.com/christophergandrud/d3Network/master/JSONdata/flare.json"
Flare <- getURL(URL)

# Convert to list format
Flare <- rjson::fromJSON(Flare)

# Create Graph
d3Tree(List = Flare, fontsize = 8, diameter = 800)


########################################################################
## Title: River plot
## Date: 2014-08-19
## Source: http://www.exegetic.biz/blog/2014/08/plotting-flows-with-riverplot/?utm_source=rss&utm_medium=rss&utm_campaign=plotting-flows-with-riverplot
########################################################################

## Set the edges
edges =
    data.frame(N1 = paste0(rep(LETTERS[1:4], each = 4),
                   rep(1:5, each = 16)),
               N2 = paste0(rep(LETTERS[1:4], 4), rep(2:6, each = 16)),
               Value = runif(80, min = 2, max = 5) *
               rep(c(1, 0.8, 0.6, 0.4, 0.3), each = 16),
               stringsAsFactors = FALSE
               )
edges = edges[sample(c(TRUE, FALSE), nrow(edges),
    replace = TRUE, prob = c(0.8, 0.2)),]
head(edges)

## Set the nodes
nodes = data.frame(ID = unique(c(edges$N1, edges$N2)),
    stringsAsFactors = FALSE)
nodes$x = as.integer(substr(nodes$ID, 2, 2))
nodes$y = as.integer(sapply(substr(nodes$ID, 1, 1), charToRaw)) - 65
rownames(nodes) = nodes$ID

## Set the colors and create the object style to associate with nodes
library(RColorBrewer)
palette = paste0(brewer.pal(4, "Set1"), "60")
styles = lapply(nodes$y, function(n) {
    list(col = palette[n+1], lty = 0, textcol = "black")
})
names(styles) = nodes$ID

## Plot the river plot
library(riverplot)
rp <- list(nodes = nodes, edges = edges, styles = styles)
class(rp) <- c(class(rp), "riverplot")
plot(rp, plot_area = 0.95, yscale=0.06)


riverExample = riverplot.example()
plot(riverExample)
plot(riverExample, srt=90, lty=1)


########################################################################
## Title: Introduction to reshape2 for Josef
## Date: 2014-08-20
########################################################################

## Install and Load library
install.packages("FAOSTAT")
library(FAOSTAT)
library(reshape2)

## Create query
query =
    data.frame(name = c("wheatProd", "riceProd", "maizeProd"),
               domainCode = rep("QC", 3),
               itemCode = c("15", "27", "56"),
               elementCode = rep("5510", 3))

## Raw data with codes
rawData.df = getFAOtoSYB(query = query, outputFormat = "long")$entity

## Convert 'n.a.' to empty string and makve value numeric. (I just
## found out this today, something changed in the API).
rawData.df$Value = as.numeric(gsub("n\\.a\\.", "", rawData.df$Value))


## Denormalize data
castedData.df = dcast(data = rawData.df,
    formula = FAOST_CODE + domainCode + itemCode +
    elementCode + name ~ Year,
    value.var = "Value")

## Normalize data, this becomes the same as the raw data.
meltedData.df =
    melt(castedData.df, variable.name = "Year",
         id.vars = c("FAOST_CODE", "domainCode", "itemCode",
             "elementCode", "name"))




########################################################################
## Title: Spatial analysis
## Date: 2014-08-25
########################################################################

x = sample(seq(0,1,0.001), replace=F)
y = sample(seq(0,1,0.001), replace=F)
z = runif(1001,min=0,max=1)
z[100] = 8
z[400] = 16
z[800] = 4

library(akima)
a = interp(x,y,z)
filled.contour(a$x,a$y,a$z)

library(fields)
test.spline = Tps(data.frame(x,y), z)
new.grid = predictSurface(test.spline, nx = 200, ny = 200)
image(new.grid)



########################################################################
## Title: Using bootMer to do model comparison in R
## Date: 2014-08-29
## Source: http://biologyforfun.wordpress.com/2014/07/13/using-bootmer-to-do-model-comparison-in-r/
########################################################################


library(lme4)
library(arm)
library(RColorBrewer)

##### work on model comparison using bootMer ##### simulate some data and fit a
##### random intercept model to them
x <- runif(100, 0, 10)
# the grouping variable
site <- gl(n = 10, k = 10)
# the random intercept effect, the simulated standard deviation around the
# intercept is 1
rnd <- rnorm(10, 0, 1)
# the simulated resposne variable, note that the fixed effect coefficient
# are 1 for the intercept and 3 for the slope. Also the simulated residuals
# will have a standard deviation of one
y <- rep(1 + rnd, each = 10) + 3 * x + rnorm(100, 0, 1)
# fit the model using Maximum Likelihood to be able to use the LRT
m1 <- lmer(y ~ x + (1 | site), REML = FALSE)

# simulate to generate credible intervals
simu <- sim(m1, n.sims = 1000)
# a new model matrix with ordered and equally spaced predictor values
new.x <- model.matrix(~x, data = data.frame(x = seq(0, 10, length.out = 100)))
new.y <- matrix(ncol = 1000, nrow = 100)
# get the predicted response values for each 1000 simulations of the fixed
# effect model parameters
new.y <- apply(simu@fixef, 1, function(x) new.x %*% x)
# compute the lower/upper quantile
lower <- apply(new.y, 1, function(x) quantile(x, prob = 0.025))
upper <- apply(new.y, 1, function(x) quantile(x, prob = 0.975))
median <- apply(new.y, 1, function(x) quantile(x, prob = 0.5))

# nice plot
pal <- brewer.pal(10, "RdYlBu")
plot(y ~ x, col = rep(pal, each = 10), pch = 16)
lines(new.x[, 2], median, col = "blue", lwd = 2)
lines(new.x[, 2], lower, col = "red", lwd = 2, lty = 2)
lines(new.x[, 2], upper, col = "red", lwd = 2, lty = 2)


# fit a second model with a random slope effect
m2 <- lmer(y ~ x + (x | site), REML = FALSE)

# using bootMer to compute 100 bootstrapped log-likelihood
b1 <- bootMer(m1, FUN = function(x) as.numeric(logLik(x)), nsim = 100)
b2 <- bootMer(m2, FUN = function(x) as.numeric(logLik(x)), nsim = 100)

# the observed LRT value
lrt <- as.numeric(-2 * logLik(m1) + 2 * logLik(m2))
# the 100 bootstrapped LRT
lrt.b <- -2 * b1$t + 2 * b2$t
# plot
quant <- quantile(lrt.b, probs = c(0.025, 0.975))
plot(1, lrt, xlab = "", ylab = "Likelihood ratio test", xaxt = "n",
     ylim = c(quant[1] +  1, quant[2] + 1))
abline(h = 0, lty = 2, lwd = 2, col = "red")
segments(1, quant[1], 1, quant[2], lend = 1)

# now simulate data from random intercept/ slope
rnd.slope <- rnorm(10, 0, 0.5)
y <- rep(1 + rnd, each = 10) + rep(3 + rnd.slope, each = 10) * x + rnorm(100,
    0, 1)

# the new models
m3 <- lmer(y ~ x + (x | site), REML = FALSE)
m4 <- lmer(y ~ x + (1 | site), REML = FALSE)

# LRT the observed values
lrt <- -2 * logLik(m4) + 2 * logLik(m3)
# the bootstrap
b3 <- bootMer(m3, FUN = function(x) as.numeric(logLik(x)), nsim = 100)
b4 <- bootMer(m4, FUN = function(x) as.numeric(logLik(x)), nsim = 100)

# the 100 bootstrapped LRT
lrt.b <- -2 * b4$t + 2 * b3$t

# the nice plot
quant <- quantile(lrt.b, probs = c(0.025, 0.975))
plot(1, lrt, xlab = "", ylab = "Likelihood ratio test", xaxt = "n",
     ylim = c(0,  quant[2] + 1))
abline(h = 0, lty = 2, lwd = 2, col = "red")
segments(1, quant[1], 1, quant[2], lend = 1)


########################################################################
## Title: Predictive Model Markup Language (PMML)
## Date:2014-08-31
########################################################################

fit = lm(Sepal.Length ~ ., data=iris)
pmml(fit)



########################################################################
## Title: Googlevis with interaction
## Date: 2014-09-02
## Source: http://lamages.blogspot.it/2014/09/zoom-zoom-googlevis.html
########################################################################

library(googleVis)
plot(
  gvisScatterChart(cars,
                   options=list(
                     explorer="{actions: ['dragToZoom',
                                          'rightClickToReset'],
                                maxZoomIn:0.05}",
                     chartArea="{width:'85%',height:'80%'}",
                     hAxis="{title: 'Speed (mph)',
                               titleTextStyle: {color: '#000000'}}",
                     vAxis="{title: 'Stopping distance (ft)',
                               titleTextStyle: {color: '#000000'}}",
                     title="Speed and stopping distances of cars in the 1920s",
                     width=550, height=500,
                     legend="none"),
                   chartid="ZoomZoom")
)


########################################################################
## Title: 2d histogram
## Date: 2014-09-02
########################################################################

## Color housekeeping
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

## Create normally distributed data for plotting
x <- rnorm(mean=1.5, 5000)
y <- rnorm(mean=1.6, 5000)
df <- data.frame(x,y)

## Plot
plot(df, pch=16, col='black', cex=0.5)

##### OPTION 1: hexbin from package 'hexbin' #######
library(hexbin)
## Create hexbin object and plot
h <- hexbin(df)
plot(h)
plot(h, colramp=rf)

## hexbinplot function allows greater flexibility
hexbinplot(y~x, data=df, colramp=rf)

## Setting max and mins
hexbinplot(y~x, data=df, colramp=rf, mincnt=2, maxcnt=60)

## Scaling of legend - must provide both trans and inv functions
hexbinplot(y~x, data=df, colramp=rf, trans=log, inv=exp)



##### OPTION 2: hist2d from package 'gplots' #######
library(gplots)

## Default call
h2 <- hist2d(df)

## Coarser binsizing and add colouring
h2 <- hist2d(df, nbins=25, col=r)

## Scaling with log as before
h2 <- hist2d(df, nbins=25, col=r, FUN=function(x) log(length(x)))


##### OPTION 4: kde2d from package 'MASS' #######
# Not a true heatmap as interpolated (kernel density estimation)
library(MASS)

## Default call
k <- kde2d(df$x, df$y)
image(k, col=r)

## Adjust binning (interpolate - can be computationally intensive for
## large datasets)
k <- kde2d(df$x, df$y, n=200)
image(k, col=r)



########################################################################
## Title: Testing similarity to calculate weights
## Date: 2014-09-04
########################################################################

a = rnorm(100, mean = 5)
b = rnorm(100, mean = 2)
c = rnorm(100, mean = 1.5)
d = rnorm(100, mean = 1.2)

test = matrix(c(a, b, c, d),nc = 4)
test2 = as.matrix(dist(t(test), diag = TRUE, upper = TRUE))
rowMeans(test2)




########################################################################
## Title: Testing prediction error for bootMer
## Date: 2014-09-09
########################################################################

library(lme4)

(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
summary(fm1)

(fm2 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy))
summary(fm2)


pe = function(x, y, newdata){
    sum((predict(x, newdata = newdata) - y)^2)
}

pe1 = bootMer(fm1,
        FUN = function(x){
            pe(x = x, y = sleepstudy$Reaction, newdata = sleepstudy)
        }, nsim = 10)


pe2 = bootMer(fm2,
        FUN = function(x){
            pe(x = x, y = sleepstudy$Reaction, newdata = sleepstudy)
        }, nsim = 10)




########################################################################
## Title: Testsing Benfords distribution
## Date: 2014-09-11
########################################################################


library(FAOSTAT)
library(data.table)
library(ggplot2)
library(BenfordTests)

dbenford = function(x){
    log(1 + 1/x, 10)
}

wheatProd.dt =
    data.table(getFAOtoSYB(name = "wheatProduction",
                           domainCode = "QC",
                           itemCode = 15,
                           elementCode = 5510)$entity
               )


num2FirstDigit = function(x){
    as.numeric(substr(as.character(x),1, 1))
}
wheatProd.dt[, wheatProduction :=
             as.numeric(gsub("n\\.a\\.", NA, wheatProduction))]
wheatProd.dt[, firstDigit := num2FirstDigit(as.numeric(wheatProduction))]

test.dt = wheatProd.dt[!is.na(wheatProduction) &
    wheatProduction != 0,
    meandigit.benftest(wheatProduction)$p.value,
    by = "Year"]
setkeyv(test.dt, "Year")

ggplot(data = test.dt, aes(x = Year, y = V1)) +
    geom_line() +
    geom_point() +
    ylab("P-value") +
    geom_abline(intercept = 0.05, slope = 0, col = "red")

leadingDigit = wheatProd.dt[!is.na(firstDigit) &
    firstDigit != 0, NROW(.SD), by = c("firstDigit")]
leadingDigit[, density := V1/sum(V1)]
leadingDigit[, benford := dbenford(firstDigit)]
with(leadingDigit, ks.test(density, benford))

ggplot(data = leadingDigit,
       aes(x = firstDigit, y = density)) +
    geom_histogram(binwidth = 1, stat = "identity") +
    geom_point(aes(x = firstDigit, y = dbenford(firstDigit)),
               col = "red") +
    geom_line(aes(x = firstDigit, y = dbenford(firstDigit)),
               col = "red") +
    scale_x_continuous(breaks = 1:9) +
    xlab("First Digit of Wheat Production") +
    ylab("Density")

leadingDigitByYear = wheatProd.dt[!is.na(firstDigit) &
    firstDigit != 0, NROW(.SD), by = c("Year", "firstDigit")]
leadingDigitByYear[, density := V1/sum(V1), by = "Year"]
leadingDigitByYear[, benford := dbenford(firstDigit)]

ggplot(data = leadingDigitByYear,
       aes(x = firstDigit, y = density)) +
    geom_histogram(binwidth = 1, stat = "identity") +
    facet_wrap(~Year) +
    geom_point(aes(x = firstDigit, y = dbenford(firstDigit)),
               col = "red") +
    geom_line(aes(x = firstDigit, y = dbenford(firstDigit)),
               col = "red") +
    scale_x_continuous(breaks = 1:9) +
    xlab("Leading Digit of Wheat Production") +
    ylab("Density")


########################################################################
## Title: Frequency table for burnt area for Simone
## Dates: 2014-09-14
########################################################################

library(raster)
test = raster("/home/mk/Downloads//burned_landtype.geotiff")
table(test[])


########################################################################
## Title: Just playing around with mathematics
## Date: 2014-09-16
########################################################################

library(animation)

## Create the square to start with
x = seq(-5, 5, length = 50)
y = seq(-5, 5, length = 50)
square = as.matrix(expand.grid(x, y))
plot(square)

## Create the rotation matrix
angle = pi/180
rotation =
    matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), ncol = 2)

## Plot
saveGIF(
    {
        init = square
        for(i in seq(0, 2 * pi, length = 360)){
            tmp = init
            distFromCenter = sqrt(tmp[, 1]^2 + tmp[, 2]^2)
            tmp[, 2] = tmp[, 2] + 10 * sin(i - distFromCenter)
            colIntensity = (tmp[, 2] + abs(min(tmp[, 2])))/
                max((tmp[, 2] + abs(min(tmp[, 2]))))
            plot(tmp[, c(1, 2)], xlim = c(-10, 10), ylim = c(-20, 20),
                 pch = ".", cex = 3, axes = FALSE, xlab = "", ylab = "",
                 col = rgb(colIntensity, 0, 0))
            init = init %*% rotation
        }
    },
    movie.name = "./wave.gif", interval = 0.005,
    nmax = 30, ani.width = 800,  ani.height = 800
    )




########################################################################
## Title: Test of using list to iterate the AUPUS
## Date: 2014-09-29
########################################################################

library(igraph)
wheat.df =
    data.frame(parent = c("wheat", "wheat", "wheat_flour"),
               item = c("wheat_flour", "wheat_bran", "bread"),
               init_existence = c(10000, 0, 0))

toTreeList = function(df){
    root = df[1, "parent"]
    baseTree = list(root = root)
    baseTree$children = as.list(df[df$parent == root, "item"])
    baseTree
}

toTreeList(wheat.df)


appendList <- function (x, val)
{
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    for (v in names(val)) {
        x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]]))
            appendList(x[[v]], val[[v]])
        else c(x[[v]], val[[v]])
    }
    x
}

wheat = list(name = "wheat")

wheat_flour = list(name = "wheat_flour")

appendList(wheat, wheat_flour)

mapply(function(x, y) x + y,
       c(a =  1, b = 2, c = 3),
       c(A = 10, B = 0, C = -10))


mapply(function(x, y) x$child = y,
       wheat.df[, 1],
       wheat.df[, 2],
       SIMPLIFY = FALSE)



base = list(wheat = list(init = 10000, production = 10000,
                share = 1, extractionRate = 1))
base$wheat$children =
    list(wheat_flour = list(init = 0, production = 0,
             share = 1, extractionRate = 0.8),
         wheat_bran = list(init = 0, production = 0,
             share = 1, extractionRate = 0.2))
base$wheat$children$wheat_flour$children =
    list(bread = list(init = 0, production = 0,
             share = 0.7, extractionRate = 1.2),
         bread = list(init = 0, production = 0,
             share = 0.3, extractionRate = 1.5))

iter = function(rootList){
    if(!is.null(rootList$children)){
        rootList =
            lapply(rootList,
                   FUN = function(x){
                   x$init = rootList$init * x$share
                   x$production = x$init * x$extractionRate
                   x
               }
               )
        if(length(rootList$children) > 1)
            rootList$children = lapply(rootList$children, iter)
    }
    rootList
}

########################################################################
## Title: binary search vs vector search
## Date: 2014-10-02
## Source: user2014.stat.ucla.edu/files/tutorial_Matt.pdf
########################################################################

pkgs <- c('data.table', 'rbenchmark')
lapply(pkgs, require, character.only = T)

load('2008.Rdata')
dt <- data.table(data)

benchmark(replications = 10, order = "elapsed",
  vector_search = {
    test1 <- dt[ArrTime == 1500 & Origin == 'ABE', ]
  },
  binary_search = {
    setkey(dt, ArrTime, Origin)
    test2 <- dt[.(1500, 'ABE'), ]
  }
)


########################################################################
## Title: Checking the entropy and KL function
## Date: 2014-10-16
########################################################################

library(entropy)
n = 1000
p = rnorm(n, sd = 5)
d = rnorm(n, sd = 10)
alphaRange = range(c(p, d), na.rm = TRUE)
p.density = hist(p, plot = FALSE,
    breaks = seq(alphaRange[1], alphaRange[2], length = n/10))$density
d.density = hist(d, plot = FALSE,
    breaks = seq(alphaRange[1], alphaRange[2], length = n/10))$density
entropy(p.density)
entropy(d.density)
KL.empirical(p.density, d.density)




########################################################################
## Title: Testing type in R
## Date: 2014-11-28
########################################################################

library(data.table)

DT = data.table(a = as.character(1:5), b = 1:5)
DT[3, a := as.numeric(5)]
DT[3, b := as.character(5)]


########################################################################
## Title: dplyr and pipe
## Date: 2014-12-02
########################################################################

pantheria <-
  "http://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR05_Aug2008.txt"
download.file(pantheria, destfile = "mammals.txt")

mammals <- read.table("mammals.txt", sep = "\t", header = TRUE,
  stringsAsFactors = FALSE)
names(mammals) <- sub("X[0-9._]+", "", names(mammals))
names(mammals) <- sub("MSW05_", "", names(mammals))
mammals <- dplyr::select(mammals, Order, Binomial, AdultBodyMass_g,
  AdultHeadBodyLen_mm, HomeRange_km2, LitterSize)
names(mammals) <- gsub("([A-Z])", "_\\L\\1", names(mammals), perl = TRUE)
names(mammals) <- gsub("^_", "", names(mammals), perl = TRUE)
mammals[mammals == -999] <- NA
names(mammals)[names(mammals) == "binomial"] <- "species"
mammals <- dplyr::tbl_df(mammals) # for prettier printing

library(dplyr)
glimpse(mammals)
select(mammals, adult_head_body_len_mm)
select(mammals, adult_head_body_len_mm, litter_size)
select(mammals, adult_head_body_len_mm:litter_size)
test = select(mammals, -adult_head_body_len_mm)
filter(mammals, species == "Balaena mysticetus")
glimpse(mutate(mammals, adult_body_mass_kg = adult_body_mass_g / 1000))

x <- rnorm(10)
x %>% max
mammals %>% arrange(adult_body_mass_g)

mammals %>%
  mutate(mass_to_length = adult_body_mass_g / adult_head_body_len_mm) %>%
  arrange(desc(mass_to_length)) %>%
  select(species, mass_to_length)

mammals %>% group_by(order) %>%
  summarise(median_litter = median(litter_size, na.rm = TRUE)) %>%
  filter(median_litter > 3) %>%
  arrange(desc(median_litter)) %>%
  select(order, median_litter)


library(data.table)
library(dplyr)
test.dt = data.table(a = 1:5, b = 5:1)
test <- test.dt %>% filter(a == 3) %>% select(a)

library(magrittr)
rnorm(200) %>%
    matrix(ncol = 2) %T>%
    plot %>% # plot usually does not return anything.
    colSums

iris$Sepal.Length %<>% sqrt



########################################################################
## Title: Checking the table in the total trade data set
## Date: 2014-12-08
########################################################################



library(faosws)
## Set up testing environments
if(Sys.getenv("USER") == "mk"){
    GetTestEnvironment(
        baseUrl = "https://hqlqasws1.hq.un.fao.org:8181/sws",
        token = "3f113726-f40e-44b3-b2af-d5f0de77c386"
        ## token = "7fe7cbec-2346-46de-9a3a-8437eca18e2a"
        )
}


allCountry = GetCodeList("trade", "total_trade", "geographicAreaM49")[, code]
allItem = GetCodeList("trade", "total_trade", "measuredItemHS")[, code]
allElement = GetCodeList("trade", "total_trade", "measuredElementTrade")[, code]
allYears = GetCodeList("trade", "total_trade", "timePointYears")[, code]

dimensions =
    list(Dimension(name = "geographicAreaM49", keys = allCountry),
         Dimension(name = "measuredItemHS", keys = allItem),
         Dimension(name = "measuredElementTrade", keys = allElement),
         Dimension(name = "timePointYears", keys = allYears))

newDataKey =
    DatasetKey(domain = "trade",
               dataset = "total_trade",
               dimensions = dimensions)
test = GetData(key = newDataKey)

########################################################################
## Title: Play around with caret
## Date: 2014-12-10
## Source: http://rforwork.info/2014/12/09/predictive-modelling-fun-with-the-caret-package/
########################################################################

library(caret)
leaf = read.csv("leaf.csv", header = FALSE)
colnames(leaf) = c("Class", "Specimen_Number", "Eccentricity", "Aspect_Ratio",
            "Elongation", "Solidity", "Stoch_Convexity",
            "Isoperimetric", "Max_Ind_Depth", "Lobedness", "Avg_Intensity",
            "Avg_Contrast", "Smoothness", "Third_Moment", "Uniformity",
            "Entropy")
leaf$Class = factor(leaf$Class)
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 5,
    selectionFunction = "oneSE")
in_train = createDataPartition(leaf$Class, p=.75, list=FALSE)

trf = train(Class ~ Eccentricity + Aspect_Ratio + Elongation +
            Solidity + Stoch_Convexity + Isoperimetric +
            Max_Ind_Depth + Lobedness + Avg_Intensity +
            Avg_Contrast + Smoothness + Third_Moment +
            Uniformity + Entropy, data=leaf, method="rf", metric="Kappa",
            trControl=ctrl, subset = in_train)

tgbm = train(Class ~ Eccentricity + Aspect_Ratio + Elongation +
             Solidity + Stoch_Convexity + Isoperimetric +
             Max_Ind_Depth + Lobedness + Avg_Intensity +
             Avg_Contrast + Smoothness + Third_Moment +
             Uniformity + Entropy, data=leaf, method="gbm", metric="Kappa",
             trControl=ctrl, subset = in_train, verbose=FALSE)



resampls = resamples(list(RF = trf,
                          GBM = tgbm))
difValues = diff(resampls)
summary(difValues)

test = leaf[-in_train,]
test$pred.leaf.rf = predict(trf, test, "raw")
confusionMatrix(test$pred.leaf.rf, test$Class)

varImp(trf, scale=FALSE)

## Lets see if the results change if we changed the sampling method
## from repeatedcv to boot.

ctrl = trainControl(method = "boot", number = 10, repeats = 5,
    selectionFunction = "oneSE")
in_train = createDataPartition(leaf$Class, p=.75, list=FALSE)

trf = train(Class ~ Eccentricity + Aspect_Ratio + Elongation +
            Solidity + Stoch_Convexity + Isoperimetric +
            Max_Ind_Depth + Lobedness + Avg_Intensity +
            Avg_Contrast + Smoothness + Third_Moment +
            Uniformity + Entropy, data=leaf, method="rf", metric="Kappa",
            trControl=ctrl, subset = in_train)

tgbm = train(Class ~ Eccentricity + Aspect_Ratio + Elongation +
             Solidity + Stoch_Convexity + Isoperimetric +
             Max_Ind_Depth + Lobedness + Avg_Intensity +
             Avg_Contrast + Smoothness + Third_Moment +
             Uniformity + Entropy, data=leaf, method="gbm", metric="Kappa",
             trControl=ctrl, subset = in_train, verbose=FALSE)



resampls = resamples(list(RF = trf,
                          GBM = tgbm))
difValues = diff(resampls)
summary(difValues)


########################################################################
## Title: mboost
## Data: 2015-01-23
########################################################################





########################################################################
## Title: Hack for Adam's trade mapping
## Date: 2015-02-12
########################################################################

tmp = read.csv("~/Github/fbs_compiler/fbs-compiler/Params/RawD2D3codes.csv")
tmp2 = melt(tmp, id.var = "varx")
tmp2$variable = NULL
tmp3 = tmp2[tmp2$value!= "", ]
tmp3$varx = gsub("[a-z|\\.]", "", tmp3$varx, ignore.case = TRUE)
tmp3$value = gsub("[a-z|\\.]", "", tmp3$value, ignore.case = TRUE)
tmp3 = tmp3[order(as.numeric(tmp3$varx)), ]
colnames(tmp3) = c("hs_parent", "hs_child")


write.csv(tmp3, file = "~/Github/sws_trade/tests/extension_tree_adam.csv", row.names = FALSE, na = "")


test = data.table(read.csv(file = "~/Downloads/csv_mapping_table.csv"))
test = test[, list(Source.Item, Target.Item)]
test[, Source.Item := gsub("[^0-9]", "", Source.Item)]
test[, Target.Item := gsub("[^0-9]", "", Target.Item)]
setnames(test, c("Source.Item", "Target.Item"), c("hs_child", "cpc_child"))

tmp4 = merge(data.table(tmp3), test, by = "hs_child", allow.cartesian = TRUE, all.x = TRUE)
setnames(test, c("hs_child", "cpc_child"), c("hs_parent", "cpc_parent"))
tmp5 = merge(tmp4, test, by = "hs_parent", allow.cartesian = TRUE, all.x = TRUE)



########################################################################
## Title: Example of ggvis scatterplot on over - display MGI gene
## symbol on click - open browser with Ensembl page
## Date :2015-02-16
########################################################################

library(ggvis)

# load dataset
dataset <- read.csv("http://dl.dropboxusercontent.com/u/232839/DO_liver_variability_sex.csv", as.is=TRUE)
head(dataset)

# what to do on hover
on_hover <- function(x) {
  if(is.null(x)) return(NULL)
  mgi_symbol <- dataset$Associated.Gene.Name[x$id]
  mgi_symbol
}

# what to do on click
on_click <- function(x) {
  if(is.null(x)) return(NULL)
  ensid <- dataset$Ensembl.Gene.ID[x$id]
  ensembl_url <- paste0("http://useast.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=", ensid)
  browseURL(ensembl_url)
  NULL
}

# start ggvis
point_size = 100 # if dots are too big/small, adjust this parameter
dataset %>%
  ggvis(~protein.Sex, ~mrna.Sex, key := ~id) %>%
  layer_points(size := point_size) %>%
  add_tooltip(on_hover, "hover") %>%
  add_tooltip(on_click, "click") %>% set_options(width=600, height=600)



library(FAOSTAT)
allFBS = getFAOtoSYB(query = .LastSearch, outputFormat = "long")


########################################################################
## Title: Entropic force
## Date: 2015-02-25
########################################################################



########################################################################
## Title: Sampling of Vatican area
## Date: 2015-03-02
########################################################################



require(geosphere)
require(ggmap)
require(plotGoogleMaps)
require(grDevices)

## Coordinates of Capella Sistina
capella = geocode("capella sistina, Vatican City, Roma")
## 20x20 grid of coordinates around the Capella
g = expand.grid(lon = seq(capella$lon-0.010, capella$lon+0.010, length.out = 20),
                lat = seq(capella$lat-0.005, capella$lat+0.005, length.out = 20))

## Hull Polygon containing coordinates
p = g[c(chull(g),chull(g)[1]),]

## Address of each coordinate of grid
a = apply(g, 1, revgeocode)

## Estimated area of the vatican city
length(grep("Vatican City", a))/length(a)*areaPolygon(p)/1000/1000
s = cbind(g, a)
s$InOut = apply(s, 1, function(x) grepl('Vatican City', x[3]))+0
coordinates(s) = ~lon+lat
proj4string(s) = CRS('+proj = longlat +datum = WGS84')
ic = iconlabels(s$InOut, height = 12)
plotGoogleMaps(s, iconMarker = ic, mapTypeId = "ROADMAP", legend = FALSE)


########################################################################
## Title: John snow and open street map
## Date: 2015-03-02
########################################################################

library(HistData)
data(Snow.deaths)
data(Snow.streets)

plot(Snow.deaths[,c("x","y")], col="red", pch=19, cex=.7,xlab="", ylab="",
     xlim=c(3,20), ylim=c(3,20))
slist <- split(Snow.streets[,c("x","y")],as.factor(Snow.streets[,"street"]))
invisible(lapply(slist, lines, col="grey"))

require(KernSmooth)
kde2d <- bkde2D(Snow.deaths[,2:3], bandwidth=c(0.5,0.5))
contour(x=kde2d$x1, y=kde2d$x2,z=kde2d$fhat, add=TRUE)

library(OpenStreetMap)
map = openmap(c(lat= 51.516,   lon= -.141),
              c(lat= 51.511,   lon= -.133))
map=openproj(map, projection = "+init=epsg:27700")
plot(map)


########################################################################
## Title: caret tutorial
## Date: 2015-03-09
########################################################################


library(caret)
library(caret)
library(mlbench)
data(Sonar)
set.seed(107)
inTrain <- createDataPartition(y = Sonar$Class,
                               ## the outcome data are needed
                               p = .75,
                               ## The percentage of data in the
                               ## training set
                               list = FALSE)
## The format of the results
## The output is a set of integers for the rows of Sonar
## that belong in the training set.
str(inTrain)

training <- Sonar[ inTrain,]
testing <- Sonar[-inTrain,]


plsFit <- train(Class ~ .,
                data = training,
                method = "pls",
                ## Center and scale the predictors for the training
                ## set and all future samples.
                preProc = c("center", "scale"))


## Change resapmling method from the default bootstrap to repeated cv
## and from 10 samples to 3. and also use summary function specific
## for binary classification.

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 3,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

plsFit <- train(Class ~ .,
                data = training,
                method = "pls",
                tuneLength = 15,
                trControl = ctrl,
                metric = "ROC",
                preProc = c("center", "scale"))

plot(plsFit)

plsClasses <- predict(plsFit, newdata = testing)

confusionMatrix(data = plsClasses, testing$Class)

rdaGrid = data.frame(gamma = (0:4)/4, lambda = 3/4)
set.seed(123)
rdaFit <- train(Class ~ .,
                data = training,
                method = "rda",
                tuneGrid = rdaGrid,
                trControl = ctrl,
                metric = "ROC")
rdaFit

rdaClasses <- predict(rdaFit, newdata = testing)
confusionMatrix(rdaClasses, testing$Class)


resamps <- resamples(list(pls = plsFit, rda = rdaFit))
summary(resamps)

xyplot(resamps, what = "BlandAltman")


########################################################################
## Title: Testing of trade reliability with eigenvalue and eigenvector
## Date: 2015-04-05
########################################################################


test.map =
    data.frame(children = c(1:5),
               parent = c(2, 3, 5, 5, 1),
               weights = runif(5, 0, 1))
test.graph = graph.data.frame(test.map, directed = FALSE)
plot(test.graph)
mymat = get.adjacency(test.graph, sparse = FALSE, attr = "weights")
f = function(x, extra = NULL){
    as.vector(mymat %*% x)
}

vals =
    arpack(f,
           options = list(n = vcount(test.graph), nev = 1, ncv = 4,
               which = "LA", maxiter = 200), sym = TRUE)


test.map =
    data.frame(children = c(1:5),
               parent = c(2, 3, 1, 5, 1),
               weights = runif(5, 0, 1))
test.graph = graph.data.frame(test.map, directed = FALSE)
page.rank(test.graph, weights = test.map$weights)
(reliability = evcent(test.graph, weights = test.map$weights)$vector)
plot(test.graph, edge.label=round(E(test.graph)$weights, 3),
     vertex.size = reliability * 20)


(test.mat = matrix(round(runif(4, 0, 10)), nc = 2))
eig = eigen(test.mat, symmetric = FALSE)
eig


test.mat = matrix(c(0.2, 0.3, 0.5, 0.7, 0.1, 0.2, 0.2, 0.4, 0.4), nc = 3,
    byrow = TRUE)
tmp = eigen(test.mat)

(test.mat - diag(rep(tmp$values[1], 3))) %*% tmp$vectors[, 1]
(test.mat - diag(rep(tmp$values[2], 3))) %*% tmp$vectors[, 2]
(test.mat - diag(rep(tmp$values[3], 3))) %*% tmp$vectors[, 3]

tmp$vectors


tmp$vectors * tmp$values[1]

tmp$values * test.mat


test.mat = matrix(c(0, 1, 0, 0, 0, 1, 4, -17, 8), nc = 3, byrow = TRUE)
tmp = eigen(test.mat)
(test.mat - diag(tmp$values)) %*% tmp$vectors

########################################################################
## Title: Time series validation based on state space model
## Date: 2015-04-16
########################################################################

library(dlm)

## Generate test
n = 50
set.seed(168)
test = arima.sim(model = list(ar = 0.8), n = n, rand.gen = rcauchy)
test[sample(n, n * 0.3)] = NA

## First order linear state-space model
testBuild = function(par) {
    dlmModPoly(1, dV = par[1], dW = par[2])
}
## Estimate the paramater (Variance)
testMLE = dlmMLE(test, rep(1, 2), testBuild)

## Build the model
testMod = testBuild(testMLE$par)

## Filter the time series
testFilt = dlmFilter(test, testMod)
## The variance of the time series, the mean is the filtered value. (testFilt$m)
testVar = unlist(dlmSvd2var(testFilt$U.C, testFilt$D.C))

## Plot the result
plot(cbind(test, testFilt$m[-1]), plot.type = "s",
     col = c("black","red"), ylab = "Level", main = "test", lwd = c(2, 1),
     lty = c(1, 2))

## Plot the bounds
pl = dropFirst(testFilt$m) + qnorm(0.05, sd = sqrt(testVar[-1]))
pu = dropFirst(testFilt$m) + qnorm(0.95, sd = sqrt(testVar[-1]))
lines(pl, lty = 2, col = "brown")
lines(pu, lty = 2, col = "brown")

## Potential invalid points
points(test, col = as.numeric(test < pl))
points(test, col = as.numeric(test > pu))





library("lme4")
library("ggplot2") # Plotting
data("Orthodont",package="MEMSS")
fm1 <- lmer(
    formula = distance ~ age*Sex + (age|Subject)
    , data = Orthodont
)

newdat <- expand.grid(
    age=c(8,10,12,14)
    , Sex=c("Female","Male")
    , distance = 0
)

mm <- model.matrix(terms(fm1),newdat)
newdat$distance <- predict(fm1,newdat,re.form=NA)
## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(fm1),mm))
tvar1 <- pvar1+VarCorr(fm1)$Subject[1]  ## must be adapted for more complex models
cmult <- 2 ## could use 1.96
newdat <- data.frame(
    newdat,
    plo = newdat$distance-cmult*sqrt(pvar1),
    phi = newdat$distance+cmult*sqrt(pvar1),
    tlo = newdat$distance-cmult*sqrt(tvar1),
    thi = newdat$distance+cmult*sqrt(tvar1)
    )

#plot confidence
g0 <- ggplot(newdat, aes(x=age, y=distance, colour=Sex))+geom_point()
g0 + geom_errorbar(aes(ymin = plo, ymax = phi))+
    labs(title="CI based on fixed-effects uncertainty ONLY")
#plot prediction
g0 + geom_errorbar(aes(ymin = tlo, ymax = thi))+
    labs(title="CI based on FE uncertainty + RE variance")



########################################################################
## Title: Predictin new levels in lmer
## Date: 2015-04-19
########################################################################

library(lme4)
time = 2010:2012
estimate.df =
    data.frame(country = rep(c("USA", "Italy"),
                   each = length(time)),
               year = rep(time, 2),
               region = rep(c("America", "Europe"),
                   each = length(time)),
               value = c(1, 1.5, 2, 2, 3, 4))

test = lmer(value ~ (-1 + year|region/country), data = estimate.df)

fixedPred = drop(model.matrix(test) %*% fixef(test))
tmp = lme4:::mkNewReTrms(object = test, newdata = estimate.df,
    re.form = lme4:::reOnly(formula(test)), allow.new.levels = TRUE)
randomPred = drop(as.matrix(tmp$b %*% tmp$Zt))
fixedPred + randomPred

predict.df =
    data.frame(country = rep(c("USA", "Italy", "France"),
                   each = length(time)),
               year = rep(time, 3),
               region = rep(c("America", "Europe", "Europe"),
                   each = length(time)),
               value = rep(NA, length(time) * 3))

predict.df$predicted = predict(test, predict.df, allow.new.level = TRUE)


RHS = formula(substitute(~R,
    list(R = lme4:::RHSForm(formula(test, fixed.only = TRUE)))))

newDesignMatrix = model.matrix(RHS, predict.df)
fixedPred = drop(newDesignMatrix %*% fixef(test))

tmp = lme4:::mkNewReTrms(object = test, newdata = predict.df,
    re.form = lme4:::reOnly(formula(test)), allow.new.levels = TRUE)
randomPred = drop(as.matrix(tmp$b %*% tmp$Zt))
fixedPred + randomPred

##This is for fixed effect
pvar1 = diag(newDesignMatrix %*% tcrossprod(vcov(test), newDesignMatrix))
## Add the standard deviation of the random effect and residual
tvar1 = pvar1 + sum(sapply(VarCorr(test), FUN = function(x) x[1]))



sleepstudy2 = sleepstudy[sleepstudy$Subject != "308", ]
sleepstudy2[sample(NROW(sleepstudy2), NROW(sleepstudy2) * 0.2), "Reaction"] = NA
fm1 = lmer(Reaction ~ Days + (Days | Subject), sleepstudy2)
sleepstudy2$predicted = predict(fm1, sleepstudy2, allow.new.levels = TRUE)
xyplot(Reaction + predicted ~ Days|Subject, sleepstudy2, auto.key = TRUE)

calculateSe = function(fit, newdata){
    RHS = formula(substitute(~R,
        list(R = lme4:::RHSForm(formula(fit, fixed.only = TRUE)))))
    ## Need to have the model.matrix for new levels, thats is, we need
    ## to take the design matrix for the upper level in the hierachy.
    ##
    ## Check how the function lme4 builds the new level.
    newmm = model.matrix(RHS, newdata)
    diag(newmm %*% tcrossprod(vcov(fit), newmm))
}

sleepstudy2$se = calculateSe(fm1, sleepstudy2)
sleepstudy2$ub = with(sleepstudy2, predicted + 1 * se)
sleepstudy2$lb = with(sleepstudy2, predicted - 1 * se)
xyplot(Reaction + predicted + ub + lb ~ Days|Subject, sleepstudy2, auto.key = TRUE)

sleepstudy2$pred = predict(fm1, sleepstudy2, re.form = NA)

xyplot(Reaction + predicted + pred~ Days|Subject, sleepstudy2, auto.key = TRUE)

with(sleepstudy2, sum((Reaction - pred)^2))

sleepstudy3 = sleepstudy
fm3 = lmer(Reaction ~ Days + (1| Subject), sleepstudy3, REML = TRUE)
sleepstudy3$pred = predict(fm3, sleepstudy3, re.form = NA)

with(sleepstudy3, sum((Reaction - pred)^2)/NROW(sleepstudy3))

z = model.matrix(~Days, sleepstudy3)

## This is the variance of the fixed effect
## VarF = diag(z %*% tcrossprod(vcov(fm3), z))
VarF <- var(as.vector(lme4::fixef(fm3) %*% t(fm3@pp$X)))

var(predict(fm3, newdata = sleepstudy3, re.form = NA))

VarRand <- sum(
    sapply(
        VarCorr(fm3)[!sapply(unique(unlist(strsplit(names(ranef(fm3)),":|/"))),
                             function(l) length(unique(fm3@frame[,l])) == nrow(fm3@frame))],
        function(Sigma) {
            X <- model.matrix(fm3)
            Z <- X[,rownames(Sigma)]
            sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X) } ) )

VarDisp <- unlist(VarCorr(fm3)[sapply(unique(unlist(strsplit(names(ranef(fm3)),":|/"))), function(l) length(unique(fm3@frame[,l])) == nrow(fm3@frame))])

VarResid <- attr(lme4::VarCorr(fm3), "sc")^2

Rm <- VarF/(VarF+VarRand+VarResid)
Rc <- (VarF+VarRand)/(VarF+VarRand+VarResid)

var(predict(fm3, newdata = sleepstudy3)) + sigma(fm3)^2


test = lm(Reaction ~ Days, data = sleepstudy)
mm = model.matrix(test)
var(fitted(test)) + var(test$residuals)

var(fitted(test))/var(sleepstudy$Reaction)
var(sleepstudy$Reaction)

var(sleepstudy$Reaction) - VarF - VarResid






fm4 = lmer(Reaction ~ Days + (1 + Days| Subject), sleepstudy3, REML = FALSE)

f = function(x, newdata = sleepstudy3){
    predict(x, newdata, allow.new.levels = TRUE)
}

bootstrapPrediction = bootMer(fm4, f, nsim = 300, use.u = TRUE,
    type = "parametric")

fixedVar = var(as.vector(lme4::fixef(fm4) %*% t(fm4@pp$X)))
residVar = attr(lme4::VarCorr(fm4), "sc")^2
totalVar = fixedVar + residVar


fixedPred = predict(fm4, sleepstudy3, allow.new.levels = TRUE)
predicted.df =
    data.frame(Days = sleepstudy3$Days,
               Subject = sleepstudy3$Subject,
               Reaction = sleepstudy3$Reaction,
               Prediction = c(t(bootstrapPrediction$t)),
               ub = fixedPred + 2 * sqrt(totalVar),
               lb = fixedPred - 2 * sqrt(totalVar),
               stringsAsFactors = FALSE)

library(ggplot2)
ggplot(predicted.df, aes(x = Days, y = Reaction)) +
    geom_point(aes(x = Days, y = Prediction), col = "red", size = 1) +
    geom_point() +
    geom_line(aes(x = Days, y = ub), col = "blue") +
    geom_line(aes(x = Days, y = lb), col = "blue") +
    facet_wrap(~Subject)






check = sleepstudy3
check$predicted = predict(fm4, check, allow.new.levels = TRUE)


ggplot(check, aes(x = Days, y = Reaction)) +
    geom_point() +
    geom_point(aes(x = Days, y = predicted), col = "red") +
    facet_wrap(~Subject)




ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
test.df = data.frame(weight = weight, group = group)
test.lm <- lm(weight ~ group, data = test.df)
predict(test.lm, test.df, se.fit = TRUE)

test.mm = model.matrix(test.lm)

sd(residuals(test.lm)) * solve(t(test.mm) %*% test.mm)


########################################################################
## Title: data.tree package test
## Date: 2015-04-21
########################################################################

library(data.tree)
wheatTree = Node$new("Wheat Tree")
wheatTree$q = 0
wheatTree$er = 1
wheatFlour = wheatTree$AddChild("Wheat Flour")
wheatFlour$er = 0.8
wheatFlour$q = 1000
wheatGerm = wheatTree$AddChild("Wheat Germ")
wheatGerm$er = 0.05
wheatGerm$q = 300
wheatGermOil = wheatGerm$AddChild("Wheat Germ Oil")
wheatGermOil$er = 0.9
wheatGermOil$q = 300
wheatBran = wheatTree$AddChild("Wheat Bran")
wheatBran$er = 0.15
wheatBran$q = 500
wheatBread = wheatFlour$AddChild("Wheat Bread")
wheatBread$er = 1
wheatBread$q = 10000
print(wheatTree, "er", "q")


calculateEffectiveRate = function(node){
    if(!node$isRoot){
        result = node$er * node$parent$er
    } else {
          result = node$er
      }
    if(length(result) == 0) {
        if (node$isLeaf){
            stop(paste("Cost for ", node$name, " not available!"))
        } else {
              result =
                  node$children %>%
                  sapply(standardize) %>%
                  sum
          }
    }
    return (result)
}

## Calculate effective extraction rate
wheatTree$Get(calculateEffectiveRate, assign = "eer")
print(wheatTree, "eer", "er", "q")

## Compute standardization
wheatTree$Get(function(x) x$eer * x$q, assign = "stdq")
print(wheatTree, "eer", "er", "q", "stdq")

## Somehow this isn't working, probably missing something
wheatTree$Get("Aggregate", "stdq", sum)

## Convert to data.frame and perform standardization.
wheatTree.df = wheatTree$ToDataFrame("name", "q", "eer", "stdq")
sum(wheatTree.df$stdq)





wheatTree$Get(standardize, assign = "check")

print(wheatTree, "check", "er", "q")

## Calculate standardized quantity
wheatTree$Get(function(x) sum(x$er * x$q), assign = "test")

## Standardization
sum(wheatTree$Get("stdq"))

## Sum up standardize quantity somehow it doesn't work
wheatTree$Get("Aggregate", "stdq", sum)
wheatTree$ToDataFrame()



library(data.tree)
acme = Node$new("Acme Inc.")
accounting = acme$AddChild("Accounting")
software = accounting$AddChild("New Software")
standards = accounting$AddChild("New Accounting Standards")
research = acme$AddChild("Research")
newProductLine = research$AddChild("New Product Line")
newLabs = research$AddChild("New Labs")
it = acme$AddChild("IT")
outsource = it$AddChild("Outsource")
agile = it$AddChild("Go agile")
goToR = it$AddChild("Switch to R")

print(acme)
acme$isRoot

software$cost = 1000000
standards$cost = 500000
newProductLine$cost = 2000000
newLabs$cost = 750000
outsource$cost = 400000
agile$cost = 250000
goToR$cost = 50000

software$p = 0.5
standards$p = 0.75
newProductLine$p = 0.25
newLabs$p = 0.9
outsource$p = 0.2
agile$p = 0.05
goToR$p = 1


acmedf = as.data.frame(acme)
acme$ToDataFrame()

acmedf$level <- acme$Get("level")
acmedf$cost <- acme$Get("cost")
acme$ToDataFrame("level", "cost")
print(acme, "level", "cost")

check = acme$ToDataFrame("level",
                 probability = acme$Get("p"),
                 cost = acme$Get("cost")
                )

f = function(x) sum(x^2)
acme$Get("Aggregate", "cost", sum)
acme$Sort("expectedCost", decreasing = TRUE)
print(acme, "expectedCost")



########################################################################
## Title: Compute density from a sample
## Date: 2015-04-22
########################################################################


sampleDensity = function(sample){
    sample.dens = density(sample, n = length(sample))
    ## predict density based on splines
    with(sample.dens, spline(x = x, y = y, xout = sample))$y
}

test = rnorm(10000)
test.den = sampleDensity(test)
plot(test.den, dnorm(test))
abline(0, 1, col = "red")

library(truncnorm)

test2 = rtruncnorm(100000, 0, 1000, 20, 30)
test2.den = sampleDensity(test2)
## plot(test2.den, dnorm(test2))
## abline(0, 1, col = "red")
hist(test2, freq = FALSE, breaks = length(test2)/100)
lines(density(test2, n = length(test2)), col = "red")
curve(dtruncnorm(x, 0, 1000, 20, 30), col = "green", add = TRUE, n = 1001)
lines(density(test2, n = length(test2), from = 0, kernel = "epanechnikov"),
      col = "blue")


test2 = rtruncnorm(10000, 0, 1000, 20, 30)
## check = density(test2, n = length(test2), from = 0, kernel = "epanechnikoav")
## check = density(test2, n = length(test2), from = 0, kernel = "rectangular")
hist(test2, freq = FALSE, breaks = log2(length(test2) * 100) + 1)
curve(dtruncnorm(x, 0, 1000, 20, 30), -10, 100, n = 1000, add = TRUE)
## lines(density(test2, n = length(test2), from = 0, kernel = "gaussian"),
##       col = "red")
lines(density(test2, n = length(test2), from = 0, kernel = "epanechnikov"),
      col = "orange")
## lines(density(test2, n = length(test2), from = 0, kernel = "rectangular"),
##       col = "yellow")
lines(density(test2, n = length(test2), from = 0, kernel = "biweight"),
      col = "green")
lines(density(test2, n = length(test2), from = 0, kernel = "cosine"),
      col = "blue")
lines(density(test2, n = length(test2), from = 0, kernel = "optcosine"),
      col = "violet")

########################################################################
## Title: Check on eigenvector centrality
## Date: 2015-04-27
########################################################################


library(igraph)
tf1 =
    data.frame(r = c("a", "a"),
               p = c("b", "c"),
               e = c(0.5, 0.5))
tf1.graph = graph.data.frame(tf1, directed = FALSE)
V(tf1.graph)$ind =
    evcent(tf1.graph, scale = FALSE,
           weights = E(tf1.graph)$e)$vector
plot(tf1.graph, edge.label = E(tf1.graph)$e,
     vertex.label = round(V(tf1.graph)$ind, 2))

tf2 =
    data.frame(r = c("a", "a", "a"),
               p = c("b", "c", "d"),
               e = c(0.5, 0.5, 0.5))
tf2.graph = graph.data.frame(tf2, directed = FALSE)
V(tf2.graph)$ind =
    evcent(tf2.graph, scale = FALSE,
           weights = E(tf2.graph)$e)$vector
plot(tf2.graph, edge.label = E(tf2.graph)$e,
     vertex.label = round(V(tf2.graph)$ind, 2))

## No gurantee that the index of bilateral partners is unique
tf3 =
    data.frame(r = c("a", "a", "b"),
               p = c("c", "b", "d"),
               e = c(0.5, 0.6, 0.5))
tf3.graph = graph.data.frame(tf3, directed = FALSE)
V(tf3.graph)$ind =
    evcent(tf3.graph, scale = FALSE,
           weights = E(tf3.graph)$e)$vector
plot(tf3.graph, edge.label = E(tf3.graph)$e,
     vertex.label = round(V(tf3.graph)$ind, 2))


tf4 =
    data.frame(r = c("a", "a", "b"),
               p = c("c", "b", "d"),
               e = c(0.2, 0.6, 0.5))

tf4.graph = graph.data.frame(tf4, directed = FALSE)
V(tf4.graph)$ind =
    evcent(tf4.graph, scale = FALSE,
           weights = E(tf4.graph)$e)$vector
plot(tf4.graph, edge.label = E(tf4.graph)$e,
     vertex.label = round(V(tf4.graph)$ind, 2))



########################################################################
## Title: Density estimatino for truncated distribution
## Date: 2015-05-01
########################################################################

library(truncnorm)
truncpoint = 15


loessDen = function(x){
    bin = hist(x, plot = FALSE, breaks = nclass.FD)
    den = with(bin,
        loess(density ~ mids, control = loess.control(surface = "direct"),
          enp.target = log(length(test)) - 3))
    ## newx = seq(truncpoint, max(test), length = 1000)
    ## predDen = predict(den, newdata = data.frame(mids = newx))
    x = seq(truncpoint, max(test), length = 1000)
    predDen = predict(den, newdata = data.frame(mids = x))
    predDen[predDen < 0] = 0
    list(x = x, y = predDen)
}

test = rtruncnorm(10000, truncpoint, 1000, 20, 10)
hist(test, breaks = length(test)/100, freq = FALSE)
with(loessDen(test), lines(x, y, col = "red", pch = "."))
lines(density(test, from = truncpoint), col = "green")
curve(dtruncnorm(x, truncpoint, 1000, 20, 10), add = TRUE, col = "blue", n = 1001)


## In order to implement the symmetric density estimation, we need to
## estimate the mode first.

mode = function(x){
    ## den = density(x, from = min(x), to = max(x), n = 512 * 2)
    ## ind = with(den, which.max(y))
    ## den$x[ind]
    den = hist(x, plot = FALSE, breaks = nclass.FD)
    ind = with(den, which.max(density))
    den$mids[ind]
}


check = replicate(10, sample(test, size = length(test) * 1000, replace = TRUE))
mean(apply(check, 2, mode))

mean(replicate(100, mode(rtruncnorm(10000, truncpoint, 1000, 20, 10))))




########################################################################
## Title: Minimise information discriminant balancing
## Date: 2015-05-09
########################################################################

## Need the lagrange multiplier

f = function(par){
    abs(1 * (log(par[1]) - log(10000)) + 0.6 * (log(par[2]) - log(2000)) +
                                                0.85 * (log(par[3]) - log(7000)))
}

f = function(par){
    1 * (log(10000/par[1])) + 0.6 * (log(2000/par[2])) +
                                                0.85 * (log(7000/par[3]))
}

optim(par = c(1, 1, 1), fn = f, lower = rep(1e-20, 3), method = "L-BFGS-B")
optim(par = c(1, 1, 1), fn = f)


f(c(1, 1, 1))
f(c(100000, 100000, 100000))

curve(-0.5 * log(0.5/x), 0, 1, n = 1001)
abline(v = 0.5)
abline(h = 0)


curve(1 * log(1/x), 0, 1, n = 1001)
curve(0.9 * log(0.9/x), add = TRUE, col = "red", n = 1001)
curve(0.8 * log(0.8/x), add = TRUE, col = "blue", n = 1001)


## minimize the function of x + y subject to x^2 + y^2 = 1
f = function(par){
    par[1]+ par[2] + 1/sqrt(2) * (par[1]^2 + par[2]^2 - 1)
}

## The minimum is -sqrt(2)
optim(par = c(0, 0), fn = f)


f = function(par){
    dnorm(par[1], mu = 0, sd = 1) + dnorm(par[2], mu = 0, sd = 2)
}

## Define the default distribution, then calculate the lambda. Then
## define the function and then pass into the optim for solving.

## Immplement the optimization then determine whether the standard
## deviation actually matters.


## This is the entropy balance!
cx = 0.9
cy = 0.6
mx = 1000
my = 100
lambda = (mx - my)/(1/(2 * cx) + 1/(2 * cy))

f = function(par){
    cx * (par[1] - mx)^2 + cy * (par[2] - my)^2 + lambda * (par[1] - par[2])
}
optim(par = c(mx, my), fn = f)

## Improvements:
## (1) change the cx and cy to the log
##
## (2) take absolute value rather than quadratic, so that we have an
##     interpretation of for every quantity change corresponds to
##     information gain in bit.



## Need to re-parametize the C in this case
cx = 1 - 0.9
cy = 1 - 0.6
mx = 1000
my = 900
lambda = (my - mx)/(1/(2 * log(cx)) + 1/(2 * log(cy)))

f = function(par){
    -log(cx) * (par[1] - mx)^2 - log(cy) * (par[2] - my)^2 +
        lambda * (par[1] - par[2])
}
optim(par = c(mx, my), fn = f)

## Can change the quadratic to the normal distribution which is the
## maximum entropy distribution of a continuous variable. Maybe use
## the standard normal distribution?

## Can use this to solve lambdas
test = D(expression(log(cx) * log(x)), "x")
f = function(x, cx) eval(test)

## Optimal implementation!
## This shows that the size of sigma does not matter, the solution is the same
cx = 1
cy = 0.1
mx = 1000
my = 1
sigma = 100
lambda = (mx - my)/(sigma^2/(cx) + sigma^2/(cy))

f = function(par){
    -cx * dnorm(par[1], mean = mx, sd = sigma, log = TRUE) -
        cy * dnorm(par[2], mean = my, sd = sigma, log = TRUE) +
            lambda * (par[1] - par[2])
}
optim(par = c(mx, my), fn = f)


test = rexp(100000, 5)
mean(test)
D(expression(log(lambda * exp(-lambda * (x - mu)))), "x")


curve(1/10 * exp(-1/10 * x), 0, 10)



## Optimal implementation!
## This shows that the size of sigma does not matter, the solution is the same
##
## This assumes the density is normal and f(x, y) is of the form cx *
## log(d(xbar)/d(x)).
cx = 1
cy = 0.1
mx = 100
my = 100
sigma = 1000
lambda = (my - mx)/(-sigma^2/(2 * cx) + -sigma^2/(2 * cy))

f = function(par){
    cx * ((par[1] - mx)/sigma)^2 + cy * ((par[2] - my)/sigma)^2 +
        lambda * (par[1] - par[2])
}
optim(par = c(mx, my), fn = f)


## Using the gamma distribution with the alpha, beta
## parametization. We set the beta as 1 for simplicity while
## calculating alpha to ensure the mode is equal to the provided
## value. This also makes the variance the same as the mean, which
## makes the balancing algorithm in favor of adjusting larger values.
##
## Maybe make the variance the same?
cx = 0.85
cy = 0.6
mx = 10000
ax = (mx + 1)
bx = 1
my = 7000
ay = my + 1
by = 1

## The formula for the lambda is really complicated. See if this can
## be simplified.
lambda = -(((ay - 1) * cy * bx * cx - (ax - 1) * cx * by * cy)/
               ((ax - 1) * cx + (ay -1 ) * cy))


f = function(par){
    cx * (dgamma(mx, shape = ax, rate = bx, log = TRUE) -
          dgamma(par[1], shape = ax, rate = bx, log = TRUE)) +
    cy * (dgamma(my, shape = ay, rate = by, log = TRUE) -
          dgamma(par[2], shape = ay, rate = by, log = TRUE)) +
    lambda * (par[1] - par[2])
}
optim(par = c(mx, my), fn = f)

curve(dgamma(x, shape = ay, scale = by), 0, 300, n = 1001)

## The variance is the same as the mode
test = rgamma(100000, shape = ay, scale = by)
mean(test)
var(test)

## Need to reparametize the gamma distribution according to the log of mx and my.

## Further extend the p(x) to x * p(x) and treat it as the posterior distribution




## Using the gamma distribution with the alpha, beta
## parametization. We set the beta as 1 for simplicity while
## calculating alpha to ensure the mode is equal to the provided
## value. This also makes the variance the same as the mean, which
## makes the balancing algorithm in favor of adjusting larger values.
##
## Maybe make the variance the same?
cx = 0.85
cy = 0.6
mx = 10000
ax = cx * mx^2
bx = cx * mx
my = 7000
ay = cy * my^2
by = cy * my


lambda = ((ay - 1) * bx - (ax - 1) * by)/((- ax - ay + 2))

f = function(par){
    -dgamma(par[1], shape = ax, rate = bx, log = TRUE) +
    -dgamma(par[2], shape = ay, rate = by, log = TRUE) +
    lambda * (par[1] - par[2])
}
optim(par = c(mx, my), fn = f)

## This means that the flag serve as a discrete prior distribution for
## the variability.


## Using the principle of minimising discriminant information for balancing.
cx = 0.85
cy = 0.6
mx = 10000
ax = cx * mx^2
bx = cx * mx
my = 7000
ay = cy * my^2
by = cy * my
kx = dgamma(mx, shape = ax, rate = bx)
ky = dgamma(my, shape = ay, rate = by)

lambda = ((ay - 1) * bx - (ax - 1) * by)/(-(ax - 1)/ky - (ay - 1)/kx)

f = function(par){
    kx * (dgamma(mx, shape = ax, rate = bx, log = TRUE) -
              dgamma(par[1], shape = ax, rate = bx, log = TRUE)) +
    ky * (dgamma(my, shape = ay, rate = by, log = TRUE) -
              dgamma(par[2], shape = ay, rate = by, log = TRUE)) +
    lambda * (par[1] - par[2])
}
optim(par = c(mx, my), fn = f)

curve(dgamma(x, shape = ax, rate = bx), n = 1001, 0, 10000)


curve(dgamma(x, shape = ay, rate = by), n = 1001, add = TRUE, col = "red")



## It seems the larger value are more in favor due to the fact that
## the distribution is non-symmetric.
cx = 0.85^2
cy = 0.8^2
mx = 1000
ax = cx * mx^2
bx = cx * mx
my = 7000
ay = cy * my^2
by = cy * my
kx = dgamma(mx, shape = ax, rate = bx)
ky = dgamma(my, shape = ay, rate = by)

lambda = ((ay - 1) * bx - (ax - 1) * by)/(-(ax - 1)/ky - (ay - 1)/kx)

f = function(par){
    kx * (dgamma(mx, shape = ax, rate = bx, log = TRUE) -
              dgamma(par[1], shape = ax, rate = bx, log = TRUE)) +
    ky * (dgamma(my, shape = ay, rate = by, log = TRUE) -
              dgamma(par[2], shape = ay, rate = by, log = TRUE)) +
    lambda * (par[1] - par[2])
}
optim(par = c(100, 200), fn = f)


## A different formulation would be to explicitly calculate a residual
## balance which represents the items unaccounted for. This residual
## balance is them distributed to other elements which maximises the
## probability. The distribution of the residual is a degenerate
## distribution at point zero. This would make the optimization an
## unconstrained one. But then this would not necessary balance.


test = matrix(rexp(1000000), nc = 1000)
hist(colMeans(test), breaks = 100, freq = FALSE)
curve(dnorm(x), add = TRUE)


########################################################################
## Title: Simulation of time series
## Date: 2015-06-04
########################################################################

library(forecast)
library(MASS)
test = arima.sim(model = list(ar = 0.8), n = 10000,
    rand.gen = function(n) rexp(n, rate = 0.2))


test = arima.sim(model = list(ar = 0.8), n = 10000,
    rand.gen = function(n) runif(n, min = -50, max = 50))

test = arima.sim(model = list(ar = 0.8), n = 10000,
    rand.gen = function(n) rlnorm(n, 0, 0.5))

hist(test, freq = FALSE, breaks = 100)
param = fitdistr(test, densfun = "normal")
curve(dnorm(x, mean = param$estimate[1], sd = param$estimate[2]),
      add = TRUE, col = "red", n = 1001)
plot(test)


########################################################################
## Title: Bertrand's paradox
## Date: 2015-06-14
########################################################################

## Initialize
n = 1e5
radius = 0.5

## length of the size of the equilateral triangle based on the radius
len = 2 * sqrt(radius^2 - (radius/2)^2)

jpeg(filename = "bertrand_illustration.jpeg", width = 960, height = 720,
     quality = 100)
par(mfrow = c(3, 4))

## Generate the circle
t = seq(0, 2 * pi, length=100)
coords = t(rbind(sin(t) * radius, cos(t) * radius))



## Random endpoints, approximately 33%
randomPointsLength = 2 * radius * sin(runif(n, 0, pi)/2)
sum(randomPointsLength > len)/n

hist(randomPointsLength, breaks = 1000,
     main = "Distribution of Chords Length", xlab = "", ylab = "Random End Points")
abline(v = len, col = "red")

theta = runif(2 * n, 0, 2 * pi)
x = radius * cos(theta)
y = radius * sin(theta)
x1 = x[1:n]
y1 = y[1:n]
x2 = x[(n + 1):(2 * n)]
y2 = y[(n + 1):(2 * n)]
plot.new()
plot.window(xlim = c(-radius, radius), ylim = c(-radius, radius))
lines(coords, col = "steelblue")
title(main = "Chord Distribution")
for(i in 1:800){
    lines(c(x1[i], x2[i]), c(y1[i], y2[i]), col = "#0000ff10")
}

plot.new()
plot.window(xlim = c(-radius, radius), ylim = c(-radius, radius))
title(main = "Mid Point Distribution")
t = seq(0, 2 * pi, length=100)
coords = t(rbind(sin(t) * radius, cos(t) * radius))
lines(coords, col = "steelblue")
points(c(x2[1:800] - x1[1:800])/2, c(y2[1:800] - y1[1:800])/2,
       col = "#0000ff22", pch = 19)

radiusLength = sqrt(((x2 - x1)/2)^2 + ((y2 - y1)/2)^2)
hist(radiusLength, breaks = 100, xlab = "", ylab = "",
     main = "Distribution of Radius")

## Random radius
radiusLength = runif(n, 0, radius)
randomRadiusLength = 2 * sqrt(radius^2 - radiusLength^2)
sum(randomRadiusLength > len)/n

hist(randomRadiusLength, breaks = 1000, xlab = "",
     ylab = "Random Radius", main = "")
abline(v = len, col = "red")

theta = runif(n, 0, 2 * pi)
r = runif(n, 0, radius)
thetachange = acos(r/radius)
theta1 = theta + thetachange
theta2 = theta - thetachange
x1 = radius * cos(theta1)
y1 = radius * sin(theta1)
x2 = radius * cos(theta2)
y2 = radius * sin(theta2)
plot.new()
plot.window(xlim = c(-radius, radius), ylim = c(-radius, radius))
t = seq(0, 2 * pi, length=100)
coords = t(rbind(sin(t) * radius, cos(t) * radius))
lines(coords, col = "steelblue")
## midpoints = c(cos(theta[1]), sin(theta[1])) * r[1]
## points(midpoints[1], midpoints[2], col = "red")
for(i in 1:800){
    lines(c(x1[i], x2[i]), c(y1[i], y2[i]),col = "#0000ff10")
}
## lines(radius * c(0, cos(theta)), radius * c(0, sin(theta)))

plot.new()
plot.window(xlim = c(-radius, radius), ylim = c(-radius, radius))
lines(coords, col = "steelblue")
points(r[1:800] * cos(theta[1:800]), r[1:800] * sin(theta[1:800]),
       col = "#0000ff22", pch = 19)

hist(radiusLength, breaks = 100, main = "", xlab = "", ylab = "")



## Random midpoints
##

t = 2 * pi * runif(n)
u = runif(n) + runif(n)
r = ifelse(u > 1, 2 - u, u)/2
x = r * cos(t)
y = r * sin(t)

midpointLength = sqrt(x^2 + y^2)
randomMidpointLength = 2 * sqrt(radius^2 - midpointLength^2)
sum(randomMidpointLength > len)/n

hist(randomMidpointLength, breaks = 100, main= "", ylab = "Random Mid Points")
abline(v = len, col = "red")

theta = runif(n, 0, 2 * pi)
r = sqrt(x^2 + y^2)
thetachange = acos(r/radius)
theta1 = theta + thetachange
theta2 = theta - thetachange
x1 = radius * cos(theta1)
y1 = radius * sin(theta1)
x2 = radius * cos(theta2)
y2 = radius * sin(theta2)
plot.new()
plot.window(xlim = c(-radius, radius), ylim = c(-radius, radius))
lines(coords, col = "steelblue")
for(i in 1:800){
    lines(c(x1[i], x2[i]), c(y1[i], y2[i]), col = "#0000ff10")
}


plot.new()
plot.window(xlim = c(-radius, radius), ylim = c(-radius, radius))
lines(coords, col = "steelblue")
points(c(x2[1:800] - x1[1:800])/2, c(y2[1:800] - y1[1:800])/2,
       col = "#0000ff22", pch = 19)

hist(midpointLength, breaks = 100, ylab = "", xlab = "", main = "")
dev.off()

########################################################################
## Title: Bertrand's Paradox
## Date: 2015-06-18
## Source: https://aschinchon.wordpress.com/2015/05/13/bertrand-or-the-importance-of-defining-problems-properly/
########################################################################

library(ggplot2)
n=1000
opt=theme(legend.position="none",
          panel.background = element_rect(fill="white"),
          panel.grid = element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          axis.text =element_blank())
#First approach
angle=runif(2*n, min = 0, max = 2*pi)
pt1=data.frame(x=cos(angle), y=sin(angle))
df1=cbind(pt1[1:n,], pt1[((n+1):(2*n)),])
colnames(df1)=c("x1", "y1", "x2", "y2")
df1$length=sqrt((df1$x1-df1$x2)^2+(df1$y1-df1$y2)^2)
p1=ggplot(df1) + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour=length>sqrt(3)), alpha=.4, lwd=.6)+
  scale_colour_manual(values = c("gray75", "red"))+opt
#Second approach
angle=2*pi*runif(n)
pt2=data.frame(aa=cos(angle), bb=sin(angle))
pt2$x0=pt2$aa*runif(n)
pt2$y0=pt2$x0*(pt2$bb/pt2$aa)
pt2$a=1+(pt2$x0^2/pt2$y0^2)
pt2$b=-2*(pt2$x0/pt2$y0)*(pt2$y0+(pt2$x0^2/pt2$y0))
pt2$c=(pt2$y0+(pt2$x0^2/pt2$y0))^2-1
pt2$x1=(-pt2$b+sqrt(pt2$b^2-4*pt2$a*pt2$c))/(2*pt2$a)
pt2$y1=-pt2$x0/pt2$y0*pt2$x1+(pt2$y0+(pt2$x0^2/pt2$y0))
pt2$x2=(-pt2$b-sqrt(pt2$b^2-4*pt2$a*pt2$c))/(2*pt2$a)
pt2$y2=-pt2$x0/pt2$y0*pt2$x2+(pt2$y0+(pt2$x0^2/pt2$y0))
df2=pt2[,c(8:11)]
df2$length=sqrt((df2$x1-df2$x2)^2+(df2$y1-df2$y2)^2)
p2=ggplot(df2) + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour=length>sqrt(3)), alpha=.4, lwd=.6)+
scale_colour_manual(values = c("gray75", "red"))+opt
#Third approach
angle=2*pi*runif(n)
radius=runif(n)
pt3=data.frame(x0=sqrt(radius)*cos(angle), y0=sqrt(radius)*sin(angle))
pt3$a=1+(pt3$x0^2/pt3$y0^2)
pt3$b=-2*(pt3$x0/pt3$y0)*(pt3$y0+(pt3$x0^2/pt3$y0))
pt3$c=(pt3$y0+(pt3$x0^2/pt3$y0))^2-1
pt3$x1=(-pt3$b+sqrt(pt3$b^2-4*pt3$a*pt3$c))/(2*pt3$a)
pt3$y1=-pt3$x0/pt3$y0*pt3$x1+(pt3$y0+(pt3$x0^2/pt3$y0))
pt3$x2=(-pt3$b-sqrt(pt3$b^2-4*pt3$a*pt3$c))/(2*pt3$a)
pt3$y2=-pt3$x0/pt3$y0*pt3$x2+(pt3$y0+(pt3$x0^2/pt3$y0))
df3=pt3[,c(6:9)]
df3$length=sqrt((df3$x1-df3$x2)^2+(df3$y1-df3$y2)^2)
p3=ggplot(df3) + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour=length>sqrt(3)), alpha=.4, lwd=.6)+scale_colour_manual(values = c("gray75", "red"))+opt


########################################################################
## Title: Distriutions
## Date:2015-06-29
########################################################################



obs = 22000
library(truncnorm)

curve(dnorm(x, 22000, 0.0915), 21999, 22001, n = 1001, ylim = c(0, 5))
curve(dcauchy(x, 22000, 0.09947), add = TRUE, col = "red", n = 1001)
curve(dlnorm(x, log(22000), 8.798935e-06), add = TRUE, col = "blue", n = 1001)
curve(dtruncnorm(x, a = 0, b = Inf, mean = 22000,
                 sd = exp(log(0.8))/(sqrt(2 * pi * exp(1)) * (1 - dnorm(0)))),
      add = TRUE, col = "green", n = 100)




curve(dnorm(x) - dnorm(x, mean = 3), -5, 5)
abline(h = 0, col = "red")

foo = function(param){
    dnorm(param) - dnorm(param, mean = -1)
}


## Monte carlo optimisation
n = 1e7
p = rnorm(n)
i = rnorm(n, mean = 2)
f = rnorm(n, mean = 3)
s = p + i - f
ll = -dnorm(p, log = TRUE) -
    dnorm(i, mean = 2, log = TRUE) -
    dnorm(f, mean = 3, log = TRUE) -
    dnorm(s, mean = 0, log = TRUE)
index = which.min(ll)
round(p[index], 1)
round(i[index], 1)
round(f[index], 1)
round(s[index], 1)
ll[index]



library(truncnorm)
library(faoswsFlag)
flagProduction = c("E", "E", "E")
valueProduction = 1000
flagImport = c("", "", "")
valueImport = 500
flagFood = c("E", "I", "I")
valueFood = 1000
flagSeed = c("E", "E", "I")
valueSeed = 600


faoswsFlagTable =
    data.frame(flagObservationStatus = c("", "T", "E", "I", "M"),
               flagObservationWeights = c(1, 0.8, 0.7, 0.75, 0.1))


selfInformation = function(prob){
    -log(prob)
}

siProduction = sum(selfInformation(flag2weight(flagProduction)))
siImport = sum(selfInformation(flag2weight(flagImport)))
siFood = sum(selfInformation(flag2weight(flagFood)))
siSeed = sum(selfInformation(flag2weight(flagSeed)))

ddegenerate = function(x, obsValue) ifelse(x == obsValue, 0, -1e5)

entropyTruncNormal = function(a = -Inf, b = Inf, mean = 0, sd = 1){
    alpha = (a - mean)/sd
    beta = (b - mean)/sd
    z = pnorm(beta) - pnorm(alpha)
    alphaTrans = ifelse(a == -Inf, 0, alpha * dnorm(alpha))
    betaTrans = ifelse(b == Inf, 0, beta * dnorm(beta))
    log(sqrt(2 * pi * exp(1)) * sd * z) + (alphaTrans - betaTrans)/(2 * z)
}

## Double check the entropy of the gamma distribution
entropyGamma = function(beta, mu){
    alpha = beta * mu
    alpha - log(beta) + lgamma(alpha) + (1 - alpha) * trigamma(alpha)
}


entropyLogNorm = function(obsValue, sd){
    log(obsValue) + sd^2 + log(sqrt(2 * pi * exp(1)) * sd)
}

test = uniroot(function(x) entropyLogNorm(obsValue = 22, sd = x) - log(0.8),
        interval = c(1e-50, 1e20))

log(22) + test$root^2



distributionise = function(obsValue, selfInformation, distribution){
    parameters = parameterise(obsValue = obsValue,
        selfInformation = selfInformation, distribution = distribution)
    switch(distribution,
           `normal` = {
               list(pdf = with(parameters,
                        function(x) dnorm(x, mean = mean, sd = sd)),
                    parameters = parameters)
           },
           `cauchy` = {
               list(pdf = with(parameters,
                        function(x) dcauchy(x, location = location,
                                            scale = scale)),
                    parameters = parameters)
           },
           `truncNorm` = {
               list(pdf = with(parameters,
                        function(x) dtruncnorm(x, a = 0, b = Inf, mean = mean,
                                               sd = sd)),
                    parameters = parameters)
           },
           `logNorm` = {
               list(pdf = with(parameters,
                        function(x) dlnorm(x, meanlog = meanlog, sdlog = sdlog)),
                    parameters = parameters)
           }
           )
}

par(mfrow = c(2, 1))
obsValue = 1
selfInformation = -log(0.1)
curve(distributionise(obsValue, selfInformation, "normal")$pdf(x),
      0, 10, n = 1001, ylim = c(0, 0.5))
curve(distributionise(obsValue, selfInformation, "cauchy")$pdf(x),
      0, 10, n = 1001, add = TRUE, col = "green")
curve(distributionise(obsValue, selfInformation, "truncNorm")$pdf(x),
      0, 10, n = 1001, add = TRUE, col = "red")
curve(distributionise(obsValue, selfInformation, "logNorm")$pdf(x),
      0, 10, n = 1001, add = TRUE, col = "orange")
selfInformation = -log(0.5)
curve(distributionise(obsValue, selfInformation, "normal")$pdf(x),
      0, 10, n = 1001, ylim = c(0, 2))
curve(distributionise(obsValue, selfInformation, "cauchy")$pdf(x),
      0, 10, n = 1001, add = TRUE, col = "green")
curve(distributionise(obsValue, selfInformation, "truncNorm")$pdf(x),
      0, 10, n = 1001, add = TRUE, col = "red")
curve(distributionise(obsValue, selfInformation, "logNorm")$pdf(x),
      0, 10, n = 1001, add = TRUE, col = "orange")
## selfInformation = -log(0.01)
## curve(distributionise(obsValue, selfInformation, "normal")$pdf(x),
##       0, 10, n = 1001, ylim = c(0, 0.05))
## curve(distributionise(obsValue, selfInformation, "cauchy")$pdf(x),
##       0, 10, n = 1001, add = TRUE, col = "green")
## curve(distributionise(obsValue, selfInformation, "truncNorm")$pdf(x),
##       0, 10, n = 1001, add = TRUE, col = "red")
## curve(distributionise(obsValue, selfInformation, "logNorm")$pdf(x),
##       0, 10, n = 1001, add = TRUE, col = "orange")




distributionise = function(obsValue, selfInformation, distribution){
    switch(distribution,
           `normal` = {
               mean = obsValue
               sd = sqrt(exp(2 * selfInformation)/(2 * pi * exp(1)))
               list(mean = mean, sd = sd)
           },
           `cauchy` = {
               location = obsValue
               scale = exp(selfInformation - log(4 * pi))
               list(location = location, scale = scale)
           },
           `truncNorm` = {
               ## The truncated normal here is defined as having lower
               ## bound of zero and upper bound of infinity. That is,
               ## on the positive real line.
               mean = obsValue
               sd = uniroot(
                   function(x){
                       entropyTruncNormal(a = 0, sd = x, mean = obsValue) -
                           selfInformation
                   },
                   interval = c(1e-20, 1e20))$root
               list(mean = mean, sd = sd)
           },
           `logNorm` = {
               sdlog = uniroot(
                   function(x){
                       entropyLogNorm(obsValue = obsValue, sd = x) -
                           selfInformation
                   },
                   interval = c(1e-50, 1e20))$root
               meanlog = log(obsValue) + sdlog^2
               list(meanlog = meanlog, sdlog = sdlog)
           },
           `exponential` = {
               rate = exp(-(selfInformation - 1))
               warning("Exponential distribution has mode at zero")
               list(rate = rate)
           }
           )
}




parameterise = function(obsValue, selfInformation, distribution){
    switch(distribution,
           `normal` = {
               sd = sqrt(exp(2 * selfInformation)/(2 * pi * exp(1)))
               list(rand = function(n)
                   rnorm(n, mean = obsValue, sd = sd),
                    pdf = function(x)
                        dnorm(x, mean = obsValue, sd = sd, log = TRUE))
           },
           `cauchy` = {
               gamma = exp(selfInformation - log(4 * pi))
               list(rand =
                        function(n) rcauchy(n, location = obsValue, scale = gamma),
                    pdf =
                        function(x) dcauchy(x, location = obsValue, scale = gamma,
                                            log = TRUE))
           },
           `truncNormal` = {
               f = function(x){
                   entropyTruncNormal(a = 0, sd = x, mean = obsValue) -
                       selfInformation
               }
               sd = uniroot(f, interval = c(1e-20, 1e20))$root
               list(rand =
                        function(n) rtruncnorm(n, a = 0, mean = obsValue, sd = sd),
                    pdf = {
                        if(selfInformation == 0){
                            function(x) dtruncnorm(x, a = 0, mean = obsValue,
                                                   sd = 1e-10)
                        } else {
                            function(x) dtruncnorm(x, a = 0, mean = obsValue,
                                                   sd = sd)
                        }}
                    )
           ## },
           ## `gamma` = {
           ##     f = function(x){
           ##         entropyGamma(beta = x, mu = obsValue) - selfInformation
           ##     }
           ##     beta = uniroot(f, interval = c(1e-20, 1e20))$root
           ##     alpha = beta * obsValue
           ##     list(rand = function(n) rgamma(n, shape = alpha, rate = beta),
           ##          pdf = function(x) dgamma(x, shape = alpha, rate = beta))
           }
           )
}

distProduction =
    parameterise(valueProduction, siProduction, "truncNormal")

distImport =
    parameterise(valueImport, siImport, "truncNormal")

distFood =
    parameterise(valueFood, siFood, "truncNormal")

distSeed =
    parameterise(valueSeed, siSeed, "truncNormal")

## Generate the variables
n = 1e6
randProduction = distProduction$rand(n)
randImport = distImport$rand(n)
randFood = distFood$rand(n)
randSeed = distSeed$rand(n)
randSeed = randProduction + randImport - randFood

ll = -distProduction$pdf(randProduction) - distImport$pdf(randImport) -
    distFood$pdf(randFood) - distSeed$pdf(randSeed)
index = which.min(ll)

## We can see that since seed is all imputed and also imputed value
## has the lowest level of uncertainty, so it is adjusted the
## most. Also since import were recorded without uncertainty, the
## value remained constant.
sapply(list(randProduction, randImport, randFood, randSeed),
       function(x) round(x[index]))

## library(GenSA)
## ll = function(x){
##     -distProduction$pdf(x[1]) - distImport$pdf(x[2]) -
##     distFood$pdf(x[3]) - distSeed$pdf(x[4])
## }
## test =
##     GenSA(par = c(valueProduction, valueImport, valueFood),
##           fn = ll, lower = rep(0, 3), upper = rep(100000, 3))

## ll(c(valueProduction, valueImport, valueFood, valueSeed))
## ll(c(valueProduction - 10, valueImport, valueFood, valueSeed))
## ll(c(test$par, test$par[1] + test$par[2] - test$par[3]))

## c(valueProduction, valueImport, valueFood, valueSeed)
## c(test$par, test$par[1] + test$par[2] - test$par[3])
## c(siProduction, siImport, siFood, siSeed)

ll = function(x){
    -log(distProduction$pdf(x[1])) - log(distImport$pdf(x[2])) -
    log(distFood$pdf(x[3])) - log(distSeed$pdf(x[4]))
}

## Need to check the specification of the distribution!!
ll = function(x){
    -log(dtruncnorm(x[1], a = 0, mean = valueProduction, sd = siProduction)) -
        log(dtruncnorm(x[2], a = 0, mean = valueImport, sd = siProduction)) -
            log(dtruncnorm(x[3], a = 0, mean = valueFood, sd = siFood)) -
                log(dtruncnorm(x[4], a = 0, mean = valueSeed, sd = siSeed))
}
equalCheck = function(x) x[1] + x[2] - x[3] - x[4]

## This seems promising, but you hae to play around with the parameters
test3 =
    solnp(pars = c(valueProduction, valueImport, valueFood, valueSeed),
          fun = ll,
          eqfun = equalCheck,
          eqB = 0,
          control = list(tol =1e-10))

test3$pars
c(valueProduction, valueImport, valueFood, valueSeed)
c(siProduction, siImport, siFood, siSeed)

ll(test3$pars)
ll(c(valueProduction, valueImport, valueFood, valueSeed))


curve(distProduction$pdf(x), 21800, 22200)

## One of the main problem is in the optimisation of the degenerate
## distribution when the uncertainty is zero, this also applies to
## structural zeroes.


entropyTruncNormal(sd = 1e-20, a = 0)

curve(entropyTruncNormal(sd = x, a = 0), 1e-20, 100, n = 10001)

f = function(x){
    ## entropyTruncNormal(a = 0, mean = 2200, sd = x) - log(0.8)
    entropyTruncNormal(a = 0, mean = 2200, sd = x)
}
uniroot(f, c(1e-50, 100))
curve(f(x), 1e-50, 10, n = 10001)

curve(dgamma(x, shape = 5, scale = 2), 0, 20)


entropyNormal = function(sd = 1){
    0.5 * log(2 * pi * exp(1) * sd^2)
}

entropyTruncNormal(a = 0, sd = 2)
entropyNormal(sd = 1)

uniroot(f = function(x) entropyNormal(sd = x) - log(1), interval = c(0, 10))
uniroot(f = function(x) entropyTruncNormal(a = 0, sd = x, mean = 22000) - log(0.8),
        interval = c(1e-10, 1e7))


sqrt(exp(2 * log(1))/(2 * pi * exp(1)))

entropyGamma(2, 1)


mymu = 0.2
mysi = log(0.8)
normal = parameterise(mymu, mysi, "normal")
## curve(exp(normal$pdf(x)), 21999, 22001, n = 1001, ylim = c(0, 5))
curve(exp(normal$pdf(x)), n = 1001)
cauchy = parameterise(mymu, mysi, "cauchy")
curve(exp(cauchy$pdf(x)), add = TRUE, col = "blue", n = 1001)
truncNorm = parameterise(mymu, mysi, "truncNormal")
curve(truncNorm$pdf(x), add = TRUE, col = "red", n = 1001)
mygamma = parameterise(mymu, mysi, "gamma")
curve(gamma$pdf(x), add = TRUE, col = "red", n = 1001)


########################################################################
## Title: A Simple Intro to Bayesian Change Point Analysis Date:
## 2015-07-19 Source:
## http://qualityandinnovation.com/2015/07/14/a-simple-intro-to-bayesian-change-point-analysis/
########################################################################


library(boot)
data(coal)
y <- tabulate(floor(coal[[1]]))
y <- y[1851:length(y)]
barplot(y,xlab="years", ylab="frequency of disasters")

# initialization
n <- length(y) # number of data elements to process
m <- 1000 # target length of the chain
L <- numeric(n) # likelihood fxn has one slot per year
k[1] <- sample(1:n,1) # pick 1 random year to start at
mu[1] <- 1
lambda[1] <- 1
b1 <- 1
b2 <- 1
# now set up blank 1000 element arrays for mu, lambda, and k
mu <- lambda <- k <- numeric(m)

# start at 2, so you can use initialization values as seeds
# and go through this process once for each of your m iterations
for (i in 2:m) {
 kt <- k[i-1] # start w/random year from initialization
 # set your shape parameter to pick mu from, based on the characteristics
 # of the early ("before") chunk of your data
 r <- .5 + sum(y[1:kt])
 # now use it to pick mu
 mu[i] <- rgamma(1,shape=r,rate=kt+b1)
 # if you're at the end of the time periods, set your shape parameter
 # to 0.5 + the sum of all the frequencies, otherwise, just set the shape
 # parameter that you will use to pick lambda based on the later ("after")
 # chunk of your data
 if (kt+1 > n) r <- 0.5 + sum(y) else r <- 0.5 + sum(y[(kt+1):n])
 lambda[i] <- rgamma(1,shape=r,rate=n-kt+b2)
 # now use the mu and lambda values that you got to set b1 and b2 for next iteration
 b1 <- rgamma(1,shape=.5,rate=mu[i]+1)
 b2 <- rgamma(1,shape=.5,rate=lambda[i]+1)
 # for each year, find value of LIKELIHOOD function which you will
 # then use to determine what year to hop to next
 for (j in 1:n) {
 L[j] <- exp((lambda[i]-mu[i])*j) * (mu[i]/lambda[i])^sum(y[1:j])
 }
 L <- L/sum(L)
 # determine which year to hop to next
 k[i] <- sample(1:n,prob=L,size=1)
}


b <- 201 # treats time until the 200th iteration as "burn-in"
mean(k[b:m])
mean(lambda[b:m])
mean(mu[b:m])

library(changepoint)
results <- cpt.mean(y,method="AMOC")
cpts(results)
param.est(results)
plot(results,cpt.col="blue",xlab="Index",cpt.width=4)


########################################################################
## TItle: Use Bayesian Search Theory to find best destination
## Date: 2015-08-05
########################################################################

library(ggplot2)


## Initialisation
gridSize = 30
td = 300
ed = 5

## Set up original traveller distance distribution
d = data.frame(x = rep(seq(-gridSize, gridSize), each = 61),
    y = rep(seq(-gridSize, gridSize), times = 61))
d$PrP = dnorm(d$x, 0, sqrt(td)) * dnorm(d$y, 0, sqrt(td))

ggplot(d, aes(x = x, y = y, z = PrP)) +
    geom_point(aes(alpha = PrP)) +
        stat_contour()

## Create the distance from search point
euclideanDistance = function(x, y, x0, y0){
    dx = (x - x0)^2
    dy = (y - y0)^2
    sqrt(dx + dy)
}

## Create the search distribution
searchDistribution = function(distance){
    1 - dnorm(distance/3, sd = 1)
}

## Update distribution assuming first destination is (20, 5)
d$searchedDist = with(d, searchDistribution(euclideanDistance(x, y, 20, 5)))

ggplot(d, aes(x = x, y = y, z = searchedDist)) +
    geom_point(aes(alpha = searchedDist)) +
        stat_contour() +
            geom_point(aes(x = 20, y = 5, color = "red"))

d$uPrP = d$PrP * d$searchedDist

ggplot(d, aes(x = x, y = y, z = uPrP)) +
    geom_point(aes(alpha = uPrP)) +
        stat_contour() +
            geom_point(aes(x = 20, y = 5, color = "red"))

## Update distribution assuming second destination is (20, 0)
d$searchedDist2 = with(d, searchDistribution(euclideanDistance(x, y, 20, 0)))

ggplot(d, aes(x = x, y = y, z = searchedDist2)) +
    geom_point(aes(alpha = searchedDist2)) +
        stat_contour() +
            geom_point(aes(x = 20, y = 0, color = "red"))

d$uPrP2 = d$uPrP * d$searchedDist2

ggplot(d, aes(x = x, y = y, z = uPrP2)) +
    geom_point(aes(alpha = uPrP2)) +
        stat_contour() +
            geom_point(aes(x = 20, y = 0, color = "red"))

## Function to compute multiple search and update probability
compoundSearchDist = function(data, iter){
    searched = dnorm(data$x, 0, sqrt(td)) * dnorm(data$y, 0, sqrt(td))
    dataCopy = data
    xvec = c()
    yvec = c()
    for(i in 1:iter){
        ## Search with update probability
        x0 = sample(data$x, size = 1, prob = searched)
        xvec = c(xvec, x0)
        y0 = sample(data$y, size = 1, prob = searched)
        yvec = c(yvec, y0)
        print(paste0("searching location (", x0, ", ", y0, ")"))
        searched =
            searched *
                searchDistribution(with(dataCopy, euclideanDistance(x, y, x0, y0)))
    }
    xvec <<- xvec
    yvec <<- yvec
    newDist = dataCopy$PrP * searched
    newDist/sum(newDist)
}

d$updatedSearch = compoundSearchDist(data = d, iter = 100)
d2 = data.frame(x = xvec, y = yvec)

ggplot(d, aes(x = x, y = y, z = updatedSearch)) +
    geom_point(aes(alpha = updatedSearch)) +
    stat_contour() +
    geom_point(data = d2, aes(x = x, y = y, z = 1, col = "red"))

## NOTE (Michael): Need to check the evolution of the density to
##                 determine the right parameters

########################################################################
## Title: Leaflet example
## 2015-09-17
########################################################################

library(leaflet)

m <- leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(lng=174.768, lat=-36.852, popup="The birthplace of R")
m


########################################################################
## Title: Class for model version control
## Date: 2015-09-24
########################################################################

setClass(Class = "sws_model",
         representation(domain = "character",
                        creator = "character",
                        creation_time = "POSIXt"),
         prototype(domain = "production",
                   creator = "mk",
                   creation_time = Sys.time()))
test = new("sws_model")

########################################################################
## Title: Play with altitude data
## Date: 2015-09-29
########################################################################
library(stringr)

## Download data
fileUrl = paste0("http://www.viewfinderpanoramas.org/DEM/TIF15/15-", LETTERS[1:24], ".zip")
fileDest = sapply(fileUrl, FUN = function(x) str_extract(x, "15-[A-X].zip"))

for(file in fileUrl){
    download.file(url = fileUrl, destfile = fileDest)
}


library(raster)
library(rasterVis)
str_name = "15-B.tif"
str_name = "15-H.tif"

alt.rast = raster(str_name)
plot3D(alt.rast)



########################################################################
## Title: Play with closure
## Date: 2015-09-29
########################################################################

test = function(transport){
    function(name){
        print(paste0("You are ", name, ", taking ", transport))
    }
}
test2 = test("submarine")
test2("Michael")

test3 = test("aeroplane")
test3("Tom")
## The same time of closure for javascript works in R


########################################################################
## Title: Fizz-Buzz test
## Date: 2015-09-30
########################################################################

fizzBuzz = function(n){
    for(i in 1:n){
        div3 = i%%3 == 0
        div5 = i%%5 == 0
        div3char = "Fizz"
        div5char = "Buzz"
        if(div3 & div5){
            print(paste0(div3char, "-", div5char))
        } else if(div3){
            print(div3char)
        } else if(div5){
            print(div5char)
        } else {
            print(i)
        }
    }
}
fizzBuzz(100)

computeFactorial = function(n){
    ## tail(cumprod(1:n), 1)
    k = 1
    for(i in 1:n){
        k = k * i
    }
    k
}


########################################################################
## Title: Spatial point analysis
## Date: 2015-10-28
########################################################################

## Set up the grid

createDensityMatrix = function(index, nrow = 100, ncol = 100){
    x = rep(1:nrow, ncol)
    y = rep(1:ncol, each = nrow)
    d = sqrt((x - index[1])^2 + (y - index[2])^2)
    dmat = matrix(d, nrow = nrow, ncol = ncol, byrow = FALSE)
    ds = sqrt((nrow/2)^2 + (ncol/2)^2)/10
    pmat = dnorm(dmat, mean = 0, sd = ds)
    pmat/sum(pmat)
}
test = createDensityMatrix(index = c(20, 50), 100, 200)

test1 = function(n, ncol = 100, nrow = 100){
    base = matrix(10000, nrow = nrow, ncol = ncol)
    for(i in 1:n){
        learning = createDensityMatrix(c(sample(nrow, 1), sample(ncol, 1)),
                                       nrow = nrow, ncol = ncol)
        base = base * (1 - learning)
    }
    base/sum(base)
}

system.time({r1 = test1(100, 100, 100)})
image(r1)


dimx = 100
dimy = 100
posx = runif(100, 0, 100)
posy = runif(100, 0, 100)
## posx = rnorm(100, 50, 10)
## posy = rnorm(100, 50, 10)
gridx = rep(1:dimx, dimy)
gridy = rep(1:dimy, each = dimx)

kernel2d = function(center, posx, posy, bandx = 1, bandy = 1){
    1/(2 * pi * bandx * bandy) *
        exp(-0.5 * ((posx - center[1])/bandx)^2 -0.5 * ((posy - center[2])/bandy)^2)
}



computeKernel = function(gridx, gridy, posx, posy){
    FUN = function(gx, gy){
        ## print(c(gx, gy))
        kernel2d(c(gx, gy), posx, posy)
    }
    matrix(mapply(FUN, gx = gridx, gy = gridy), nrow = max(gridx), ncol = max(gridy))
}

test = computeKernel(gridx, gridy, posx, posy)
image(test)
points(posx/100, posy/100, pch = 19)

plot(posx, posy, pch = 19)
image(test, add = TRUE)

library(MASS)

with(geyser,
     {
         plot(duration, waiting, xlim = c(0.5,6), ylim = c(40,100))
         f1 = kde2d(duration, waiting, n = 50, lims = c(0.5, 6, 40, 100))
         image(f1, zlim = c(0, 0.05))
         f2 = kde2d(duration, waiting, n = 50, lims = c(0.5, 6, 40, 100),
             h = c(width.SJ(duration), width.SJ(waiting)) )
         image(f2, zlim = c(0, 0.05))
     })

f3 = kde2d(posx, posy, n = 100, lims = c(0, 100, 0, 100))
image(f3)


final = r1 * f3$z/sum(f3$z)
## final = final/sum(final)
par(mfrow = c(1, 3))
image(r1)
image(f3)
image(final)

zr = range(c(r1, f3$z))

persp(r1, zlim = zr)
persp(f3, zlim = zr)
persp(final, zlim = zr)



########################################################################
## Title: Function to calculate bandwidth based on lat long for fix
##        distance
## Date: 2015-10-29
########################################################################


er = 3960
d2r = pi/180
r2d = 180/pi

calcBandLng = function(lat, d){
    r = er * cos(lat * d2r)
    d/r * r2d
}

calcBandLat = function(d){
    d/er * r2d
}

calcBandLng(0, 0.2)
calcBandLat(0.2)

calcBandLng(80, 0.2)
calcBandLat(0.2)

library(XML)
map = xmlParse("map.osm")
map.lst = xmlToList(map)

## Select only nodes
mapNodes.lst = map.lst[which(names(map.lst) == "node")]

##
hasList = sapply(mapNodes.lst, FUN = function(x) is.list(x))
mapTag.lst = mapNodes.lst[hasList]



## This is the successful implementation of the cafe density. Need to
## check how to download all the nodes
map = xmlParse("city.osm")
map.lst = xmlToList(map)
mapNodes.lst = map.lst[which(names(map.lst) == "node")]

hasList = sapply(mapNodes.lst, FUN = function(x) is.list(x))
mapTag.lst = mapNodes.lst[which(hasList)]
isCafe = sapply(mapTag.lst, FUN = function(x) x$tag["v"] == "cafe")
mapCafe.lst = mapTag.lst[which(isCafe)]

cafeLoc = sapply(mapCafe.lst, FUN = function(x) as.numeric(x$`.attrs`[c("lat", "lon")]))
cafeName = sapply(mapCafe.lst, FUN = function(x) x[[2]]["v"])
plot(cafeLoc[2, ], cafeLoc[1, ])
text(cafeLoc[2, ], cafeLoc[1, ], labels = cafeName)

latLonDim = as.numeric(map.lst[[1]])

## Get map center
mapCenter = c(latLonDim[2] + diff(latLonDim[c(2, 4)])/2,
    latLonDim[1] + diff(latLonDim[c(1, 3)])/2)

## set bandwidth as 100 meters, which is an assumed walking search distance
bandLng = calcBandLng(mapCenter[2], 100/1000)
bandLat = calcBandLat(100/1000)

cafeDen = kde2d(cafeLoc[2, ], cafeLoc[1, ], n = 300,
    lims = c(latLonDim[c(2, 4)], latLonDim[c(1, 3)]),
    h = c(bandLng, bandLat))
image(cafeDen)
points(cafeLoc[2, ], cafeLoc[1, ], pch = 19)


########################################################################
## Title: Open streetmap
## Date: 2015-10-30
########################################################################

library(osmar)
src = osmsource_api()
get_osm(node(18961430), source = src)
get_osm(way(3810479), source = src)
bb = center_bbox(174.76778, -36.85056, 700, 700)
ua = get_osm(bb, source = src)
summary(ua$nodes)
ts_ids = find(ua, node(tags(v == "traffic_signals")))
bs_ids = find(ua, node(tags(v %agrep% "busstop")))



########################################################################
## Title: Commodity price and foreign exchange rate
## Date: 2016-01-04
########################################################################

## Prepare the data for the FFPI
ffpiUrl = "http://www.fao.org/fileadmin/templates/worldfood/Reports_and_docs/Food_price_indices_data.csv"
ffpi = read.csv(file = ffpiUrl, skip = 2)
ffpi$Date = as.character(ffpi$Date)
ffpi$Month = as.numeric(sapply(strsplit(ffpi$Date, "/"), FUN = function(x) x[1]))
ffpi$Year = as.numeric(sapply(strsplit(ffpi$Date, "/"), FUN = function(x) x[2]))
ffpi$Date = NULL
redudantColumns = grep("^X", colnames(ffpi))
finalffpi = ffpi[-redudantColumns]

## Scrap usd exchange rate index
library(rvest)
usdFxUrl = "http://www.federalreserve.gov/releases/h10/summary/indexn_m.htm"

usdFx =
    usdFxUrl %>%
    read_html() %>%
    html_nodes(xpath = '//*[@id="printThis"]/table') %>%
    html_table() %>%
    .[[1]]



usdFx$MonthName = sapply(strsplit(usdFx$X1, " "), FUN = function(x) x[1])
usdFx$Year = as.numeric(sapply(strsplit(usdFx$X1, " "), FUN = function(x) x[2]))
usdFx$Month = match(usdFx$MonthName, toupper(month.abb))
usdFx$MonthName = NULL
usdFx$X1 = NULL
colnames(usdFx)[which(colnames(usdFx) == "X2")] = "usdFx"

## Merge the two datasets
all.df = merge(finalffpi, usdFx, by = c("Month", "Year"), all = FALSE)

## Correlation
##
## NOTE (Michael): It might be better off to predict each of the
##                 commodity index, rather than the aggregated
##                 fpi. This way, the effect can be isolated.
##
## NOTE (Michael): There seems to be very strong correlation between
##                 cereals and oil. Further, their correlation with
##                 the US exchange rate is very similar. This raise
##                 the question whether they are competing assets
##                 driven by the same factor.
cor(all.df)
pairs(~., data = all.df)


## Exploratory plots
par(mfrow = c(3, 1))
with(all.df,
{
    plot(Food.Price.Index, type = "l", ylim = c(0, 300))
    lines(usdFx, col = "red")
    plot(Food.Price.Index, usdFx)
    abline(lm(usdFx ~ Food.Price.Index), col = "red", lty = 2)
    ccf(Food.Price.Index, usdFx)
}
)

with(all.df,
{
    n = NROW(all.df)
    plot(Food.Price.Index[-1], type = "l", ylim = c(0, 300))
    lines(usdFx[-n], col = "red")
    plot(Food.Price.Index[-1], usdFx[-n])
    abline(lm(usdFx[-n] ~ Food.Price.Index[-1]), col = "red", lty = 2)
    ccf(Food.Price.Index[-1], usdFx[-n])
})

## It would appear that the us exchange rate has a 1 month lead ahead
## of the food price index, the correlation is negative.


########################################################################
## Title: Javascript V8
## Date: 2015-02-04
## Source: https://www.opencpu.org/posts/v8-release-0-10/
########################################################################

library(V8)

## This creates a javascript console
ctx <- V8::v8()
ctx$console()

## load the iris data into javascript
var iris = console.r.get("iris")

########################################################################
## Title: Hilbert's curve
## Date:2015-02-04
## Source: https://aschinchon.wordpress.com/2016/02/01/going-bananas-with-hilbert/
########################################################################

library(reshape2)
library(dplyr)
library(ggplot2)
opt=theme(legend.position="none",
          panel.background = element_rect(fill="white"),
          panel.grid=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank())

hilbert = function(m,n,r) {
  for (i in 1:n)
  {
    tmp=cbind(t(m), m+nrow(m)^2)
    m=rbind(tmp, (2*nrow(m))^r-tmp[nrow(m):1,]+1)
  }
  check <<- m
  melt(m) %>%
      plyr::rename(c("Var1" = "x", "Var2" = "y", "value"="order")) %>%
      arrange(order)
}

## Original
ggplot(hilbert(m=matrix(1), n=1, r=2), aes(x, y)) + geom_path()+ opt
ggplot(hilbert(m=matrix(1), n=2, r=2), aes(x, y)) + geom_path()+ opt
ggplot(hilbert(m=matrix(1), n=3, r=2), aes(x, y)) + geom_path()+ opt
ggplot(hilbert(m=matrix(1), n=4, r=2), aes(x, y)) + geom_path()+ opt
ggplot(hilbert(m=matrix(1), n=5, r=2), aes(x, y)) + geom_path()+ opt
ggplot(hilbert(m=matrix(1), n=7, r=2), aes(x, y)) + geom_path()+ opt

## Changing order
ggplot(hilbert(m=matrix(.5), n=5, r=2), aes(x, y)) + geom_path()+ opt
ggplot(hilbert(m=matrix(0), n=5, r=2), aes(x, y)) + geom_path()+ opt
ggplot(hilbert(m=matrix(tan(1)), n=5, r=2), aes(x, y)) + geom_path()+ opt
ggplot(hilbert(m=matrix(3), n=5, r=2), aes(x, y)) + geom_path()+ opt
ggplot(hilbert(m=matrix(-1), n=5, r=2), aes(x, y)) + geom_path()+ opt
ggplot(hilbert(m=matrix(log(.1)), n=5, r=2), aes(x, y)) + geom_path()+ opt
ggplot(hilbert(m=matrix(-15), n=5, r=2), aes(x, y)) + geom_path()+ opt
ggplot(hilbert(m=matrix(-0.001), n=5, r=2), aes(x, y)) + geom_path()+ opt

## Polygons
ggplot(hilbert(m=matrix(log(1)), n=4, r=2), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(.5), n=4, r=2), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(tan(1)), n=5, r=2), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(-15), n=4, r=2), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(-25), n=4, r=2), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(0), n=4, r=2), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(1000000), n=4, r=2), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(-1), n=4, r=2), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(-.00001), n=4, r=2), aes(x, y)) + geom_polygon()+ opt

## Changing exponent
gplot(hilbert(m=matrix(log(1)), n=4, r=-1), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(.5), n=4, r=-2), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(tan(1)), n=4, r=6), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(-15), n=3, r=sin(2)), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(-25), n=4, r=-.0001), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(0), n=4, r=200), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(1000000), n=3, r=.5), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(-1), n=4, r=sqrt(2)), aes(x, y)) + geom_polygon()+ opt
ggplot(hilbert(m=matrix(-.00001), n=4, r=52), aes(x, y)) + geom_polygon()+ opt

## Polar coordinates
ggplot(hilbert(m=matrix(1), n=4, r=2), aes(x, y)) + geom_polygon()+
    coord_polar()+opt
ggplot(hilbert(m=matrix(-1), n=5, r=2), aes(x, y)) + geom_polygon()+
    coord_polar()+opt
ggplot(hilbert(m=matrix(.1), n=2, r=.5), aes(x, y)) + geom_polygon()+
    coord_polar()+opt
ggplot(hilbert(m=matrix(1000000), n=2, r=.1), aes(x, y)) +
    geom_polygon()+ coord_polar()+opt
ggplot(hilbert(m=matrix(.25), n=3, r=3), aes(x, y)) +
    geom_polygon()+ coord_polar()+opt
ggplot(hilbert(m=matrix(tan(1)), n=5, r=1), aes(x, y)) +
    geom_polygon()+ coord_polar()+opt
ggplot(hilbert(m=matrix(1), n=4, r=1), aes(x, y)) + geom_polygon()+
    coord_polar()+opt
ggplot(hilbert(m=matrix(log(1)), n=3, r=sin(2)), aes(x, y)) + geom_polygon()+
    coord_polar()+opt
ggplot(hilbert(m=matrix(-.0001), n=4, r=25), aes(x, y)) +
    geom_polygon()+ coord_polar()+opt



########################################################################
## Title: bayesboot: An R package for doing the Bayesian bootstrap
## Date: 2016-02-21
## Source: http://www.sumsar.net/blog/2016/02/bayesboot-an-r-package/
########################################################################

library(bayesboot)

## Heights of the last ten American presidents in cm (Kennedy to Obama).
heights = c(183, 192, 182, 183, 177, 185, 188, 188, 182, 185)


b1 = bayesboot(heights, mean)
plot(b1)

b2 = bayesboot(heights, weighted.mean, use.weights = TRUE)
plot(b2)

## Given the model and the data, this is the probability that the mean
## heights of American presidents is above the mean heights of
## American males as given by
## www.cdc.gov/nchs/data/series/sr_11/sr11_252.pdf
mean( c(b2$V1 > 175.9, TRUE, FALSE) )


## The heights of opponents of American presidents (first time they were elected).
## From Richard Nixon to John McCain
heights_opponents = c(182, 180, 180, 183, 177, 173, 188, 185, 175)

## Running the Bayesian bootstrap for both datasets
b_presidents = bayesboot(heights, weighted.mean, use.weights = TRUE)
b_opponents  = bayesboot(heights_opponents, weighted.mean, use.weights = TRUE)

## Calculating the posterior difference and converting back to a
## bayesboot object for pretty plotting.
b_diff = as.bayesboot(b_presidents - b_opponents)
plot(b_diff)


## Bayesian bootstrap with LOESS

boot_fn <- function(cars, weights) {
  loess(dist ~ speed, cars, weights = weights)$fitted
}

bb_loess = bayesboot(cars, boot_fn, use.weights = TRUE)

## Plotting the data
plot(cars$speed, cars$dist, pch = 20, col = "tomato4", xlab = "Car speed in mph",
     ylab = "Stopping distance in ft",
     main = "Speed and Stopping distances of Cars")

## Plotting a scatter of Bootstrapped LOESS lines to represent the uncertainty.
for(i in sample(nrow(bb_loess), 20)) {
  lines(cars$speed, bb_loess[i,], col = "gray")
}
## Finally plotting the posterior mean LOESS line
lines(cars$speed, colMeans(bb_loess, na.rm = TRUE), type ="l",
      col = "tomato", lwd = 4)

########################################################################
## Title: Plotly interactive library
## Date: 2016-03-05
########################################################################

library(plotly)
p <- plot_ly(midwest, x = percollege, color = state, type = "box")
p

source("https://plot.ly/~Cat_Phish/8/halo-distribution.r")

########################################################################
## Title: Large scale eigenvalue decomposition and SVD with rARPACK
## Date: 2016-03-05
## Source: http://statr.me/2016/02/large-scale-eigen-and-svd-with-rarpack/
########################################################################

library(rARPACK)
set.seed(123)
## Some random data
x = matrix(rnorm(1000 * 100), 1000)
## If retvec == FALSE, we don't calculate eigenvectors
eigs_sym(cov(x), k = 5, which = "LM", opts = list(retvec = FALSE))

## eigs_sym is for symmetric matrices, otherwise use eigs


library(Matrix)
spmat = as(cov(x), "dgCMatrix")
eigs_sym(spmat, 2)

## Implicitly define the matrix by a function that calculates A %*% x
## Below represents a diagonal matrix diag(c(1:10))
fmat = function(x, args){
    return(x * (1:10))
}
eigs_sym(fmat, 3, n = 10, args = NULL)


## Benchmark between the three different PCA methods.
library(microbenchmark)
set.seed(123)
## Some random data
x = matrix(rnorm(2000 * 500), 2000)
pc = function(x, k){
    ## First center data
    xc = scale(x, center = TRUE, scale = FALSE)
    ## Partial SVD
    decomp = svds(xc, k, nu = 0, nv = k)
    return(list(loadings = decomp$v, scores = xc %*% decomp$v))
}
microbenchmark(princomp(x), prcomp(x), pc(x, 3), times = 5)


########################################################################
## Title Quick intro to Nonnegative Matrix Factorisation
## Date: 2015-03-08
## source: https://matloff.wordpress.com/2016/03/05/quick-intro-to-nmf-the-method-and-the-r-package/
########################################################################

library(NMF)
library(pixmap)
mtr = read.pnm('http://heather.cs.ucdavis.edu/MtRush.pgm')
a = mtr@grey
aout = nmf(a,50)
w = aout@fit@W
h = aout@fit@H
approxa = w %*% h
# brightness values must be in [0,1]
approxa = pmin(approxa,1)
mtrnew = mtr
mtrnew@grey = approxa
plot(mtrnew)


########################################################################
## Title: Plotting sunflower
## Date: 2016-03-14
## Source: https://aschinchon.wordpress.com/2016/03/14/sunflowers/
########################################################################

library(deldir)
library(ggplot2)
library(dplyr)
opt = theme(legend.position  = "none",
            panel.background = element_rect(fill="red4"),
            axis.ticks       = element_blank(),
            panel.grid       = element_blank(),
            axis.title       = element_blank(),
            axis.text        = element_blank())
CreateSunFlower <- function(nob=500, dx=0, dy=0) {
    data.frame(r = sqrt(1:nob),
               t = (1:nob)*(3-sqrt(5))*pi) %>%
        mutate(x = r * cos(t) + dx,
               y = r * sin(t) + dy)
}

g = seq(from=0, by = 45, length.out = 4)
jitter(g, amount=2) %>%
  expand.grid(jitter(g, amount=2)) %>%
  apply(1, function(x){
      CreateSunFlower(nob=round(jitter(220, factor=15)),
                      dx=x[1], dy=x[2])
      }) %>%
    do.call("rbind", .) %>%
    deldir() %>%
    .$dirsgs -> sunflowers
ggplot(sunflowers) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color="greenyellow") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  opt


########################################################################
## Title: rOpenSci geospatial libraries
## Date: 2016-03-21
## Source: http://ropensci.org/blog/2016/03/17/ropensci-geospatial-stack
########################################################################

## Create data from/to geojson
library("geojsonio")
geojson_json(c(-99.74, 32.45), pretty = TRUE)

library("wellknown")
point(data.frame(lon = -116.4, lat = 45.2))

## Need authentiation to create the gist
library("gistr")
cat(geojson_json(us_cities[1:100,], lat = 'lat', lon = 'long'),
    file = "map.geojson")
gist_create("map.geojson")

library("lawn")
lawn_hex_grid(c(-96,31,-84,40), 50, 'miles') %>% view

library("geoaxe")
library("rgeos")
wkt <- "POLYGON((-180 -20, -140 55, 10 0, -140 -60, -180 -20))"
poly <- rgeos::readWKT(wkt)
polys <- chop(x = poly)
plot(poly, lwd = 6, mar = c(0, 0, 0, 0))
plot(polys, lwd = 6, mar = c(0, 0, 0, 0))


########################################################################
## Title: testing invisible on lapply
## Date: 2016-04-15
########################################################################

f1 = function(x){
    x
}

f2 = function(n){
    invisible(lapply(1:n, f1))
}
f2(10)



########################################################################
## Title: Debugging test
## Date: 2016-50-05
########################################################################

SS <- function(mu, x) {
    d <- x - mu
    d2 <- d^2
    ss <- sum(d2)
    ss
}

## set the RNG seed so that the results are reproducible
set.seed(100)
x <- rnorm(100)



debug(SS)
SS(1, x)

nLL <- function(mu, x) {
    z <- mu * x
    lz <- log(z)
    L1 <- sum(lz)
    L2 <- mu/2
    LL <- -(L1 - L2)
    LL
}

f <- function(x) {
    r <- x - g(x)
    r
}

g <- function(y) {
    r <- y * h(y)
    r
}

h <- function(z) {
    r <- log(z)
    if (r < 10)
        r^2 else r^3
}

trace("h", quote( if(is.nan(r)) { recover() } ), at = 3, print = F)

flagObservationStatus = c("", "T", "E", "I", "M")
flagObservationWeights = c(1, 0.8, 0.75, 0.5, 0)
flagMethod = c("-", "q", "p", "h", "c", "b", "i", "s", "t", "e", "f", "n", "u")

completeFlag =
    expand.grid(flagObservationStatus = flagObservationStatus,
                flagMethod = flagMethod)

write.csv(completeFlag, row.names = FALSE, file = "~/Desktop/completeFlag.csv")


########################################################################
## Title: Ensurer package
## 2016-05-16
########################################################################


library(magrittr)
library(ensurer)

letters[1:5] %>%
    ## ensure_that(is.numeric(.)) %>%
    ensure_that(## checkMethodFlag(flag = .) ~ "Flag method incorrectly specified",
        is.character(.) ~ "values are not characters",
        is.numeric(.) ~ "values are not numeric")



test = list(a = "a", b = "b", c = "c")

with(test,{
    print(a)
    print(b)
    print(c)
})



foo = function(){
    list(.Internal(stop(TRUE, .makeMessage("test", domain = NULL))))
}

test = foo()

########################################################################
## Title: Plain vanilla recurrent neural networks in R: waves prediction
## Date: 2016-08-07
## Source: http://firsttimeprogrammer.blogspot.com/2016/08/plain-vanilla-recurrent-neural-networks.html
########################################################################

# Clear workspace
rm(list=ls())

# Load libraries
require(rnn)

# Set seed for reproducibility purposes
set.seed(10)

# Set frequency
f <- 5
w <- 2 * pi * f

# Create sequences
t <- seq(0.005,2,by = 0.005)
x <- sin(t * w) + rnorm(200, 0, 0.25)
y <- cos(t * w)

# Samples of 20 time series
X <- matrix(x, nrow = 40)
Y <- matrix(y, nrow = 40)

# Plot noisy waves
plot(as.vector(X), col = 'blue', type = 'l', ylab = "X,Y", main = "Noisy waves")
lines(as.vector(Y), col = "red")
legend("topright", c("X", "Y"), col = c("blue", "red"),
       lty = c(1, 1), lwd = c(1, 1))

# Standardize in the interval 0 - 1
X <- (X - min(X))/(max(X) - min(X))
Y <- (Y - min(Y))/(max(Y) - min(Y))

# Transpose
X <- t(X)
Y <- t(Y)

# Training-testing sets
train <- 1:8
test <- 9:10

# Train model. Keep out the last two sequences.
model <- trainr(Y = Y[train,],
                X = X[train,],
                learningrate = 0.05,
                hidden_dim = 16,
                numepochs = 1500)

# Predicted values
Yp <- predictr(model, X)

# Plot predicted vs actual. Training set + testing set
plot(as.vector(t(Y)), col = 'red', type = 'l',
     main = "Actual vs predicted", ylab = "Y,Yp")
lines(as.vector(t(Yp)), type = 'l', col = 'blue')
legend("topright", c("Predicted", "Real"),
       col = c("blue", "red"), lty = c(1, 1), lwd = c(1, 1))

# Plot predicted vs actual. Testing set only.
plot(as.vector(t(Y[test,])), col = 'red', type='l',
     main = "Actual vs predicted: testing set", ylab = "Y,Yp")
lines(as.vector(t(Yp[test,])), type = 'l', col = 'blue')
legend("topright", c("Predicted", "Real"), col = c("blue", "red"),
       lty = c(1, 1), lwd = c(1, 1))




########################################################################
## Title: Asynchronous and Distributed Programming in R with the
##        Future Package.
## Date: 2016-11-05 Source:
## https://alexioannides.com/2016/11/02/asynchronous-and-distributed-programming-in-r-with-the-future-package/
########################################################################

library(future)
 
f <- future({
  samples <- rnorm(10000000)
  mean(samples)
}) %plan% multiprocess
w <- value(f)
w

f_dots <- function() {
  f <- future({
    s <- rnorm(10000000)
    mean(s)
  }) %plan% multiprocess
 
  while (!resolved(f)) {
    cat("...")
  }
  cat("\n")
 
  value(f)
}
f_dots()

library(parallel)
 
sub_means =
    mclapply(
        X = 1:4,
        FUN = function(x) {
            samples <- rnorm(25000000);
            mean(samples)
        },
        mc.cores = 4)
 
final_mean = mean(unlist(sub_means))
final_mean

single_thread_mean <- function(n) {
  samples <- rnorm(n)
  mean(samples)
}
 
multi_thread_mean <- function(n) {
  f1 <- future({ single_thread_mean(n) }) %plan% multiprocess
  f2 <- future({ single_thread_mean(n) }) %plan% multiprocess
  f3 <- future({ single_thread_mean(n) }) %plan% multiprocess
  f4 <- future({ single_thread_mean(n) }) %plan% multiprocess
 
  mean(value(f1), value(f2), value(f3), value(f4))
}
 
multi_thread_mean()

library(microbenchmark)
n_sample = 1e8
microbenchmark(single_thread_mean(n_sample),
               multi_thread_mean(n_sample), times = 10)


########################################################################
## Title: Timeseries forecasting using extreme gradient boosting
## Date: 2016-11-07
## Source: http://ellisp.github.io/blog/2016/11/06/forecastxgb
########################################################################

## Univariate example
library(forecastxgb)
model = xgbts(gas)
gasForecast = forecast(model, h = 12)
plot(gasForecast)

## Here is a summary of the model. The feature are sorted in order of
## the second column which is gain. It is unclear how this gain is
## calculated.
summary(model)

## Multivariate example
library(fpp)
consumption <- usconsumption[ ,1]
income <- matrix(usconsumption[ ,2], dimnames = list(NULL, "Income"))
consumption_model <- xgbts(y = consumption, xreg = income)
summary(consumption_model)
plot(consumption_model)

## In order to forecast the response, we will need to forecast the
## regressors.
income_future <- matrix(forecast(xgbts(usconsumption[,2]), h = 10)$mean, 
                        dimnames = list(NULL, "Income"))
plot(forecast(consumption_model, xreg = income_future))

########################################################################
## Title: Naive Bayes: A Generative Model and Big Data Classifier
## Date: 2016-11-07
## Source: https://www.rstudio.com/rviews/2016/11/02/naive-bayes-a-generative-model-and-big-data-classifier/
########################################################################

library(mlbench)
library(e1071)
data("HouseVotes84")
model <- naiveBayes(Class ~ ., data = HouseVotes84)
print(model)


########################################################################
## Title: Introduction to map mate
## Date: 2016-11-07
## Source: https://leonawicz.github.io/mapmate/mapmate.html
##         https://leonawicz.github.io/mapmate/usage_and_limitations.html
########################################################################

## Load libraries and data
library(mapmate)
library(dplyr)
library(purrr)
data(monthlytemps)
monthlytemps

## Calculate moving average, but I am curious of how the season is
## determined.
get_ma(monthlytemps, type = "seasonal", season = "winter")
get_ma(monthlytemps, type = "annual", size = 20)

## Load annual data
data(annualtemps)

library(RColorBrewer)
pal <- rev(brewer.pal(11, "RdYlBu"))
temps <- mutate(annualtemps, frameID = Year - min(Year) + 1)
frame1 <- filter(temps, frameID == 1)  # subset to first frame
id <- "frameID"

## Create static maps
save_map(frame1, z.name = "z", id = id, ortho = FALSE, col = pal,
         type = "maptiles", save.plot = FALSE, return.plot = TRUE)
save_map(frame1, z.name = "z", id = id, col = pal, type = "maptiles",
         save.plot = FALSE, return.plot = TRUE)

## Create frames of the plot.
rng <- range(annualtemps$z, na.rm = TRUE)
n <- length(unique(annualtemps$Year))
suffix <- "annual_2D"
temps <- split(temps, temps$frameID)
walk(temps, ~save_map(.x, z.name = "z", id = id, ortho = FALSE,
                      col = pal, type = "maptiles", suffix = suffix,
                      z.range = rng))

########################################################################
## Title: Predictability in Network Models
## Date: 2016-11-07
## Source: http://jmbh.github.io//Predictability-in-network-models/
########################################################################


filePath ='http://psychosystems.org/wp-content/uploads/2014/10/Wenchuan.csv'
wenchuan.df <- read.csv(filePath)
wenchuan.df <- na.omit(wenchuan.df)
p <- ncol(wenchuan.df)

# Fit a mixed graphical model and calculate the predictabiility. The
# predictability here is basically the percentage of variance of a
# single node explained by all other nodes.
library(mgm)
fit_obj <- mgmfit(data = data, 
                  type = rep('g', p),
                  lev = rep(1, p),
                  rule.reg = 'OR')
pred_obj <- predict(fit_obj, data, 
                    error.continuous = 'VarExpl')


library(qgraph)
qgraph(fit_obj$wadj, # weighted adjacency matrix as input
       layout = 'spring', 
       pie = pred_obj$error$Error, # provide errors as input
       pieColor = rep('#377EB8',p),
       node.color = fit_obj$edgecolor,
       labels = colnames(data))

########################################################################
## Title: Easy Cross Validation in R with `modelr`
## Date: 2016-11-15
## Source: http://jacobsimmering.com/2016/11/11/CrossValidationInR/
########################################################################

library(dplyr)
library(modelr)
library(purrr)
library(tibble)
glimpse(cars)

## This createas a list of datasets for cross validated. Moreover,
## they are pointers to the data and lazy loaded.
cars_cv <- cars %>%
  select(-Doors, -Chevy, -sedan) %>%
  crossv_kfold(10)
cars_cv

## We don't see the data replicated 10 times
pryr::object_size(cars)
pryr::object_size(cars_cv)


## Fit the model to all the cv datasets then extract the RMSE.
cars_cv %>% 
    mutate(train = map(train, as_tibble)) %>%
    mutate(model = map(train, ~ lm(Price ~ Mileage + Cylinder +
                                       Cruise + Sound + Leather +
                                       Buick + Cadillac + Pontiac +
                                       Saab + Saturn + convertible +
                                       coupe +  hatchback + wagon, 
                                   data = .))) %>%
    mutate(rmse = map2_dbl(model, test, rmse)) %>%
    select(.id, rmse)

## Speed up the fitting
fitter <- function(training_data, testing_data){
    lm(Price ~ Mileage + Cylinder + Cruise + Sound + Leather +
           Buick + Cadillac + Pontiac + Saab + Saturn +
           convertible + coupe + hatchback + wagon, 
       data = training_data) %>% 
        rmse(testing_data)
}

ten_fold_2 <- function(x){
    x %>% 
        mutate(rmse = map2_dbl(train, test, fitter)) %>%
        select(.id, rmse)
}

microbenchmark::microbenchmark(ten_fold_2(cars_cv))

########################################################################
## Title: Sorry ARIMA, but I’m Going Bayesian
## Date: 2017-01-11
## Source: http://multithreaded.stitchfix.com/blog/2016/04/21/forget-arima/
########################################################################

## NOTE (Michael): This actually come from the official slide of the
##                 author of the bsts package.


library(lubridate)
library(bsts)
library(dplyr)
library(ggplot2)
library(forecast)

data("AirPassengers")

## Subset the data to create the training set
airPassengersTrain = 
    AirPassengers %>%
    window(., start = c(1949, 1), end = c(1959,12))

## Fitting ARIMA model
airPassengerArimaFit =
    arima(log10(airPassengersTrain), 
          order=c(0, 1, 1), 
          seasonal=list(order=c(0,1,1), period=12))

## Create the fitted data frame
airPassengerFit.df = 
    airPassengerArimaFit %>%
    {
        fitted = 10^as.numeric(fitted(.))
        predicted = 10^as.numeric(predict(., n.ahead = 12)$pred)
        data.frame(predicted = c(fitted, predicted),
                   actual = as.numeric(AirPassengers),
                   Date = as.Date(time(AirPassengers))
                   )
    }

## Calculate the MAPE
airPassengerMAPE =
    airPassengerFit.df %>%
    filter(., year(Date) > 1959) %>%
    ## summarise(check = mean(actual - predicted))
    summarise(MAPE = mean(abs(actual - predicted)/actual))

## Plot the result
ggplot(data = airPassengerFit.df, aes(x = Date)) +
    geom_line(aes(y = actual, colour = "actual"), size = 1.2) +
    geom_line(aes(y = predicted, colour = "predicted"),
              size = 1.2, linetype = 2) +
    theme_bw() +
    theme(legend.title = element_blank()) + 
    ylab("") +
    xlab("") +
    geom_vline(xintercept = as.numeric(as.Date("1959-12-01")),
               linetype = 2) +
    ggtitle(paste0("ARIMA -- Holdout MAPE = ",
                   round(100 * airPassengerMAPE, 2), "%")) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0))

## Fit the BSTS
logAirPassengers = log10(airPassengersTrain)
ss = AddLocalLinearTrend(list(), logAirPassengers)
ss = AddSeasonal(ss, logAirPassengers, nseasons = 12)
airPassengersBSTS =
    bsts(logAirPassengers, state.specification = ss,
         niter = 500, ping = 0, seed = 2016)

## Get a suggested number of burn-ins
burn = SuggestBurn(0.1, airPassengersBSTS)

## Create the plot data frame
airPassengerBSTS.df = 
    airPassengersBSTS %>%
    {
        onestopPE = airPassengersBSTS$one.step.prediction.errors
        fitted = 10^(-colMeans(onestopPE[-(1:burn), ]) +
                     logAirPassengers)
        predicted =
            predict.bsts(airPassengersBSTS, horizon = 12,
                             burn = burn, quantiles = c(0.025, 0.975))
        data.frame(predicted = c(fitted, 10^predicted$mean),
                   lowerBound =
                       c(rep(NA, length(fitted)),
                         10^as.numeric(predicted$interval[1, ])),
                   upperBound =
                       c(rep(NA, length(fitted)),
                         10^as.numeric(predicted$interval[2, ])), 
                   actual = as.numeric(AirPassengers),
                   Date = as.Date(time(AirPassengers)))
    }

## Calculate the MAPE, the MAPE is actually higher than that of ARIMA
airPassengerBSTSMAPE =
    filter(airPassengerBSTS.df, year(Date) > 1959) %>%
    summarise(MAPE = mean(abs(actual - predicted)/actual))

## Plot the result of BSTS.
ggplot(data = airPassengerBSTS.df, aes(x = Date)) +
    geom_line(aes(y = actual, colour = "actual"), size = 1.2) +
    geom_line(aes(y = predicted, colour = "predited"), size = 1.2,
              linetype = 2) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    ylab("") +
    xlab("") +
    geom_vline(xintercept = as.numeric(as.Date("1959-12-01")),
               linetype = 2) + 
    geom_ribbon(aes(ymin = lowerBound, ymax = upperBound),
                fill = "grey", alpha = 0.5) +
    ggtitle(paste0("BSTS -- Holdout MAPE = ",
                   round(100 * airPassengerBSTSMAPE,2), "%")) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0))


########################################################################
## Title: Entropy Based Image Binarization with imager and FSelectorRcpp
## Date: 2017-01-12
## Source: http://r-addict.com/2017/01/08/Entropy-Based-Image-Binarization.html
########################################################################

library(imager)
my_photo = load.image("~/Desktop/deepblu_profile.jpg")

## Create greyscale of photo
layout(t(c(1, 2)))
plot(my_photo, main = "My photo")
plot(grayscale(my_photo), main = "My photo in a grayscale")

# Scaling value.
layout(t(c(1,2,3)))
plot(my_photo)
plot(my_photo/2) # nothing happens on a plot 
plot(my_photo/2, rescale = FALSE)

layout(t(1))
# %>% operator means for
# f(a) %>% g(par1, par2) == g(f(a), par1, par2)
# and simplifies pipelines
grayscale(my_photo) %>% hist(main="Luminance values \nin my photo")

library(ggplot2)
mdf <- as.data.frame(my_photo)
head(mdf, 3)
mdf <- plyr::mutate(mdf, channel = factor(cc, labels = c('R', 'G', 'B')))
ggplot(mdf, aes(value, col = channel)) +
   geom_histogram(bins = 30) + 
   facet_wrap(~ channel)

f <- ecdf(grayscale(my_photo))
f(grayscale(my_photo)) %>% hist(main="Transformed luminance values")

layout(t(1:3))
plot(my_photo, main = "Original photo")
plot(grayscale(my_photo), main = "In a grayscale")
f(grayscale(my_photo)) %>% 
   as.cimg(dim=dim(grayscale(my_photo))) %>% 
    plot(main="With histogram equalisation")

f <- ecdf(grayscale(my_photo))
mdf <- f(grayscale(my_photo)) %>% 
         as.cimg(dim=dim(grayscale(my_photo))) %>%
         as.data.frame()

ggplot(mdf, aes(x, y)) + # create ggplot object
   geom_raster(aes(fill = value)) + # use raster geometry 
   scale_x_continuous(expand = c(0,0)) + # remove margins on X axis
   scale_y_continuous(expand = c(0,0), # remove margins on Y axis
      trans = scales::reverse_trans()) + # reverse Y axis
   scale_fill_gradient(low = "black", high = "white") # scale in grayscale - the default scale for ggplot is blue


my_photo_he <- f(grayscale(my_photo)) %>% 
         as.cimg(dim=dim(grayscale(my_photo)))
layout(t(matrix(1:6, ncol = 2, nrow = 3)))
plot(my_photo)
threshold(my_photo_he, thr = 0.2) %>%
    plot(main = "Determinant: \nover 0.2 value")
threshold(my_photo_he, thr = 0.4) %>%
    plot(main = "Determinant: \nover 0.4 value")
threshold(my_photo_he, thr = "25%") %>%
    plot(main = "Determinant: \n75% highest values")
threshold(my_photo_he, thr = "50%") %>%
    plot(main = "Determinant: \n50% highest values")
threshold(my_photo_he, thr = "auto") %>% plot(main = "auto")


########################################################################
## Title: When You Went too Far with Survival Plots During the
## survminer 1st Anniversary
## Date: 2017-01-20
## Source: http://r-addict.com/2017/01/15/Too-Far-With-Survival-Plots.html
########################################################################


devtools::install_github("kassambara/survminer", build_vignettes = TRUE)
library("survminer")
library("survival")
fit<- survfit(Surv(time, status) ~ sex, data = lung)
ggsurvplot(
   fit,                     # survfit object with calculated statistics.
   risk.table = TRUE,       # show risk table.
   pval = TRUE,             # show p-value of log-rank test.
   conf.int = TRUE,         # show confidence intervals for 
                            # point estimaes of survival curves.
   xlim = c(0,500),         # present narrower X axis, but not affect
                            # survival estimates.
   xlab = "Time in days",   # customize X axis label.
   break.time.by = 100,     # break X axis in time intervals by 500.
   ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
                            # in legend of risk table.
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Male", "Female"),    # change legend labels.
  palette = 
    c("#E7B800", "#2E9FDF"),# custom color palettes.
  main = "Survival curves",                       # specify the title of the plot
  submain = "Based on Kaplan-Meier estimates",    # the subtitle of the plot 
  caption = "created with survminer",             # the caption of the plot
  font.main = c(16, "bold", "darkblue"),          # font for titles of the plot, the table and censor part
  font.submain = c(15, "bold.italic", "purple"),  # font for subtitles in the plot, the table and censor part
  font.caption = c(14, "plain", "orange"),        # font for captions in the plot, the table and censor part
  font.x = c(14, "bold.italic", "red"),           # font for x axises in the plot, the table and censor part
  font.y = c(14, "bold.italic", "darkred"),       # font for y axises in the plot, the table and censor part
  font.tickslab = c(12, "plain", "darkgreen"),    # font for ticklabs in the plot, the table and censor part
  ########## risk table #########,
  risk.table.title = "Note the risk set sizes",          # the title of the risk table
  risk.table.subtitle = "and remember about censoring.", # the subtitle of the risk table
  risk.table.caption = "source code: website.com",       # the caption of the risk table
  risk.table.height = 0.35,                              # the height of the risk table
  ########## ncensor plot ######
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.title = "Number of censorings",           # as above but for the censoring plot
  ncensor.plot.subtitle = "over the time.",
  ncensor.plot.caption = "data available at data.com",
  ncensor.plot.height = 0.35
)


########################################################################
## Title: Library mlr
## Date: 2017-03-10
## Source: http://mlr-org.github.io/mlr-tutorial/release/html/
########################################################################

library(mlr)
data(iris)

## Define the task
task = makeClassifTask(id = "tutorial", data = iris, target = "Species")

## Define the learner
lrn = makeLearner("classif.lda")

## Define the resampling strategy
rdesc = makeResampleDesc(method = "CV", stratify = TRUE)

## Do the resampling
r = resample(learner = lrn, task = task, resampling = rdesc,
             show.info = FALSE)

## Get the mean misclassification error
r$aggr


########################################################################
## Title: M-value visualisation and analysis
## Date: 2017-09-14
########################################################################

library(dplyr)
compartmentCount = 16

## Workman formulation
halfTime = c(2.65, 7.94, 12.2, 18.5, 26.5, 37, 54, 79, 114, 146, 185,
             238, 304, 397, 503, 635)

## m0 = c(111.9, 89.1, 75.2, 68.8, 63.5, 57.3, 53.2, 51.9, 51.9, 50.2,
##        50.2, 47.3, 42.6, 42.6, 42.6, 42.6)

m0 = c(34.2, 27.2, 22.9, 21.0, 19.3, 17.4, 16.2, 15.8, 15.8, 15.3,
       15.3, 14.4, 12.9, 12.9, 12.9, 12.9)

m1 = c(1.2195, 1.2195, 1.2121, 1.1976, 1.1834, 1.1628, 1.1494, 1.1236,
       1.1236, 1.0707, 1.0707, 1.0593, 1.0395, 1.0395, 1.0395, 1.0395)

## Buhlman formulation
a = m0 - m1
b = 1/m1

pressure = seq(0, 10, length = 100)
depth = (pressure - 1) * 10

## Pressure plot
plot(pressure, pressure, type = "l")
for(i in 1:compartmentCount){
    mValue = (a[i] + (b[i] * pressure * 1000))/1000
    ## lines(pressure, mValue, col = rgb(i/16, 0, 0, ))
    lines(pressure, mValue, col = rgb(i/16, 0, 0, ))
}


## Depth plot
plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 10))
abline(0, 1, lwd = 5)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
box()
axis(1)
axis(2)

for(i in 1:compartmentCount){
    mValue = (a[i] + (b[i] * pressure * 1000))/1000
    ceiling = (mValue - 1) * 10
    lines(ceiling, depth, col = rgb(i/16, 0, 0, ))
}

bottomTime = 25
depth = 55

saturationLevel = function(depth, bottomTime){
    halfLifeCount = bottomTime/halfTime
    pressure = (depth/10) + 1
    pressureTissue = pressure * (1 - (0.5 ^ halfLifeCount))
    mValue = (a + b * pressure * 1000)/1000
    print(pressureTissue)
    print(mValue)
}
    
saturationLevel(depth = depth, bottomTime = bottomTime) %>% hist(., breaks = 16)


########################################################################
## Title: Model-based Boosting in R
## Date: 2018-03-03
## Source: https://cran.r-project.org/web/packages/mboost/vignettes/mboost_tutorial.pdf
########################################################################

library(mboost)
data("bodyfat", package = "TH.data")

## Reproduce formula of Garcia et al., 2005
lm1 <- lm(DEXfat ~ hipcirc + kneebreadth + anthro3a, data = bodyfat)
coef(lm1)

## Estimate same model by glmboost
glm1 <- glmboost(DEXfat ~ hipcirc + kneebreadth + anthro3a, data = bodyfat)
coef(glm1)

## Estimate with all variables
##
glm2 <- glmboost(DEXfat ~ ., data = bodyfat)
coef(glm2)
## Specifying the 'which' arguement in coef will list all the
## coefficient even if they are zero.
coef(glm2, which = "")

## Plot the coefficient path, this is similar to LARS
plot(glm2, off2int = TRUE)


## now change ylim to the range of the coefficients without intercept (zoom-in)
preds <- names(bodyfat[, names(bodyfat) != "DEXfat"])
plot(glm2, ylim = range(coef(glm2, which = preds)))


## Use a GAM model without assuming for linearity
gam1 <- gamboost(DEXfat ~ bbs(hipcirc) + bbs(kneebreadth) + bbs(anthro3a),
                 data = bodyfat)
par(mfrow = c(1,3)) ## 3 plots in one device
plot(gam1)
coef(gam1)


########################################################################
## Title: Get Started with XGBoost
## Date: 2018-03-06
## Source: https://xgboost.readthedocs.io/en/latest/get_started/
########################################################################

library(xgboost)
## load data
data(agaricus.train, package='xgboost')
data(agaricus.test, package='xgboost')
train <- agaricus.train
test <- agaricus.test
## fit model
bst <- xgboost(data = train$data, label = train$label, max.depth = 2, eta = 1, nround = 2,
               nthread = 2, objective = "binary:logistic")
## predict
pred <- predict(bst, test$data)

## Use the custom data type which supports initial prediction value,
## weight training and speed up the training.
dtrain <- xgb.DMatrix(train$data, label = train$label)
class(dtrain)


## Define custom objective and evaluation function
logregobj <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    preds <- 1/(1 + exp(-preds))
    grad <- preds - labels
    hess <- preds * (1 - preds)
    return(list(grad = grad, hess = hess))
}

evalerror <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    err <- sqrt(mean((preds-labels)^2))
    return(list(metric = "MSE", value = err))
}

dtest <- xgb.DMatrix(test$data, label = test$label)
watchlist <- list(eval = dtest, train = dtrain)
param <- list(max_depth = 2, eta = 1, silent = 1)
bst <- xgb.train(param, dtrain, nrounds = 2, watchlist, logregobj, evalerror, maximize = FALSE)


########################################################################
## Simulate
########################################################################

x1 = rnorm(100)
x2 = 2 + 1 * x1 + rnorm(100)

mod = lm(x2 ~ x1)
simulate(mod)


## The function simulate is really helpful, especially if you want to
## bootstrap from the distribution.
