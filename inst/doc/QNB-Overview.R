### R code from vignette source 'QNB-Overview.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=60)


###################################################
### code chunk number 2: input files
###################################################
library(QNB)
f1 <- system.file("extdata", "meth1.txt", package="QNB")
f2 <- system.file("extdata", "meth2.txt", package="QNB")
f3 <- system.file("extdata", "unmeth1.txt", package="QNB")
f4 <- system.file("extdata", "unmeth2.txt", package="QNB")

meth1 <- read.table(f1,header=TRUE)
meth2 <- read.table(f2,header=TRUE)
unmeth1 <- read.table(f3,header=TRUE)
unmeth2 <- read.table(f4,header=TRUE)
head(meth1)
head(unmeth1)



###################################################
### code chunk number 3: Differential methylation analysis
###################################################
result = qnbtest(meth1, meth2,unmeth1,unmeth2,mode="per-condition")
head(result)


###################################################
### code chunk number 4: <Differetial methylation analysis without replicates
###################################################
f1 <- system.file("extdata", "no_rep_meth1.txt", package="QNB")
f2 <- system.file("extdata", "no_rep_meth2.txt", package="QNB")
f3 <- system.file("extdata", "no_rep_unmeth1.txt", package="QNB")
f4 <- system.file("extdata", "no_rep_unmeth2.txt", package="QNB")

no_rep_meth1 <- read.table(f1,header=TRUE)
no_rep_meth2 <- read.table(f2,header=TRUE)
no_rep_unmeth1 <- read.table(f3,header=TRUE)
no_rep_unmeth2 <- read.table(f4,header=TRUE)
head(no_rep_meth1)
head(no_rep_unmeth1)
result = qnbtest(no_rep_meth1, 
                 no_rep_meth2,
                 no_rep_unmeth1,
                 no_rep_unmeth2,
                 mode="blind")


###################################################
### code chunk number 5: Differetial methylation analysis automatically
###################################################
result = qnbtest(meth1, meth2,unmeth1,unmeth2)