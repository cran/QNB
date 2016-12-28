\name{qnbtest}
\alias{qnbtest}
\title{
This is the main function of differential RNA methylation analysis for MeRIP-seq.
}
\description{
This function tests for differential RNA methylation between two conditions.
}
\usage{
qnbtest(control_ip,treated_ip,control_input,treated_input,
size.factor=NA,mode="auto",plot.dispersion=TRUE,output.dir = NA)
}

\arguments{
  \item{control_ip}{
a  reads count data frame for IP samples of control
}
  \item{treated_ip}{
a  reads count data frame for IP samples of treated 
}
  \item{control_input}{
a  reads count data frame for input samples of control
}
  \item{treated_input}{
a  reads count data frame for input samples of treated
}
  \item{size.factor}{
A list of size factor. The size factor of the IP and input sample of the biological replicate and directly reflect their sequencing depth. If size.factor=NA, QNB will compute the size factor of IP and input sample of each replicate. If user could provide the size.factor, the name of each term must be control_ip, treated_ip, control_input, treated_input in list. 
}
  \item{mode}{
There are four ways how the empirical dispersion can be computed:

\itemize{
  \item{pooled -}{
Use the samples from all conditions with replicates to estimate a single pooled empirical dispersion value, called "pooled", and assign it to all samples. In this mode,
the number of replicates must be the same.
}
  \item{per-condition -}{
For each condition with replicates, compute an empirical dispersion value by considering the data from samples for this condition. 
}
  \item{blind -}{
Ignore the sample labels and compute an empirical dispersion value as if all samples were replicates of a single condition. This can be done even if there are no biological 
replicates.}
  \item{auto -}{
Select mode according to the size of samples automaticly.The default is auto. By default, QNB package implements the "per-conditon" mode for a more sensitive estimation of the raw variance parameter when biological replicates are provided; while the "blind" mode is implemented when biological replicates are not available.
}
}

}
  \item{plot.dispersion}{
The default is TRUE. If plot.dispersion=FALSE, it will not save the dispersion figure. 
}
  \item{output.dir}{
The saved file path. The default is NA. If output.dir=NA, the path is the current path.
}
}
\value{
Results will get a matrix including 7 columns (p.treated, p.control, log2.RR, log2.OR, pvalue, q, padj).
  \item{p.treated}{
  The percentage of methylation under treated condition.
  }
  \item{p.control}{
  The percentage of methylation under control condition.
  }
  \item{log2.RR}{
  The normalized risk ratio.
  }
   \item{log2.OR}{
  The normalized odds ratio.
  }
  \item{pvalue}{
Indicate the significance of the methylation site as an RNA differential methylation site
}
  \item{q}{
  The standardized feature abundance, which is proportional to the expression level of the RNA transcript.
  }
  \item{padj}{
  The FDR of the methylation site, indicating the significance of the peak as an RNA differential methylation site after multiple hypothesis correction usting BH method.  
  }
}
\author{
Lian Liu <liulian19860905@163.com>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## The function is currently defined as load library and specify the parameters

library("QNB")

f1 = system.file("extdata", "control_ip.txt", package="QNB")
f2 = system.file("extdata", "treated_ip.txt", package="QNB")
f3 = system.file("extdata", "control_input.txt", package="QNB")
f4 = system.file("extdata", "treated_input.txt", package="QNB")

meth1 = read.table(f1,header=TRUE)
meth2 = read.table(f2,header=TRUE)
unmeth1 = read.table(f3,header=TRUE)
unmeth2 = read.table(f4,header=TRUE)

# When there are replicates under two conditions, we could select 
# mode="per-condition" or mode="pooled" to estimate the dispersion. 
# The default is mode="auto".

result = qnbtest(meth1,meth2,unmeth1,unmeth2,mode="per-condition")

# When size.factor is not NA

## Not run:
size.factor = list(control_ip=c(0.9,0.9,0.9),
                   treated_ip=c(0.9,0.9,0.9), 
                   control_input=c(0.9,0.9,0.9),
                   treated_input=c(0.9,0.9,0.9))
                   
result = qnbtest(meth1,meth2,unmeth1,unmeth2,
                 size.factor=size.factor,
                 mode="per-condition")
## End(Not run)


# If you have replicates for one condition but not for the other, or without any 
# replicates for # two conditions, you can select mode="blind" to estimate 
# the dispersion. 

f1 = system.file("extdata", "no_rep_controlip.txt", package="QNB")
f2 = system.file("extdata", "no_rep_treatedip.txt", package="QNB")
f3 = system.file("extdata", "no_rep_controlinput.txt", package="QNB")
f4 = system.file("extdata", "no_rep_treatedinput.txt", package="QNB")

no_rep_meth1 = read.table(f1,header=TRUE)
no_rep_meth2 = read.table(f2,header=TRUE)
no_rep_unmeth1 = read.table(f3,header=TRUE)
no_rep_unmeth2 = read.table(f4,header=TRUE)

## Not run: 
result = qnbtest(no_rep_meth1, 
                 no_rep_meth2,
                 no_rep_unmeth1,
                 no_rep_unmeth2,
                 mode="blind")
## End(Not run)

# If you could not decide which mode to estimate dispersion, mode="auto" 
# will select suitable way to estimate dispersion according to the replicates.

## Not run: 
result = qnbtest(meth1, meth2,unmeth1,unmeth2)
## End(Not run)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ differential RNA methylation }
\keyword{ quad-negative-binomial model }% __ONLY ONE__ keyword per line