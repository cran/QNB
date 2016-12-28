### R code from vignette source 'QNB-Overview.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=60)


###################################################
### code chunk number 2: input files
###################################################
library(QNB)
f1 = system.file("extdata", "control_ip.txt", package="QNB")
f2 = system.file("extdata", "treated_ip.txt", package="QNB")
f3 = system.file("extdata", "control_input.txt", package="QNB")
f4 = system.file("extdata", "treated_input.txt", package="QNB")

meth1 = read.table(f1, header=TRUE)
meth2 = read.table(f2, header=TRUE)
unmeth1 = read.table(f3, header=TRUE)
unmeth2 = read.table(f4, header=TRUE)
head(meth1)
head(unmeth1)



###################################################
### code chunk number 3: Differential methylation analysis
###################################################
result = qnbtest(meth1, meth2,unmeth1,unmeth2,mode="per-condition")
head(result)


###################################################
### code chunk number 4: provide size factor
###################################################
size.factor = list(control_ip=c(0.9,0.9,0.9),
                   treated_ip=c(0.9,0.9,0.9), 
                   control_input=c(0.9,0.9,0.9),
                   treated_input=c(0.9,0.9,0.9))

result = qnbtest(meth1,meth2,unmeth1,unmeth2,
                 size.factor=size.factor,
                 mode="per-condition")

###################################################
### code chunk number 5: <Differetial methylation analysis without replicates
###################################################
f1 = system.file("extdata", "no_rep_controlip.txt", package="QNB")
f2 = system.file("extdata", "no_rep_treatedip.txt", package="QNB")
f3 = system.file("extdata", "no_rep_controlinput.txt", package="QNB")
f4 = system.file("extdata", "no_rep_treatedinput.txt", package="QNB")

no_rep_meth1 = read.table(f1, header=TRUE)
no_rep_meth2 = read.table(f2, header=TRUE)
no_rep_unmeth1 = read.table(f3, header=TRUE)
no_rep_unmeth2 = read.table(f4, header=TRUE)
head(no_rep_meth1)
head(no_rep_unmeth1)
result = qnbtest(no_rep_meth1, 
                 no_rep_meth2,
                 no_rep_unmeth1,
                 no_rep_unmeth2,
                 mode="blind")


###################################################
### code chunk number 6: Differetial methylation analysis automatically
###################################################
result = qnbtest(meth1, meth2, unmeth1, unmeth2)