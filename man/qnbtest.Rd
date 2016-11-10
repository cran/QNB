\name{qnbtest}
\alias{qnbtest}
\title{
This is the main function of differential RNA methylation analysis for MeRIP-seq.
}
\description{
This function tests for differential RNA methylation between two conditions.
}
\usage{
qnbtest(meth1, meth2, unmeth1, unmeth2, 
        mode="auto", 
        plot.dispersion=TRUE,
        output.dir=NA)
}

\arguments{
  \item{meth1}{
        a  reads count matrix for IP samples of one condition
}
  \item{meth2}{
        a  reads count matrix for IP samples of two condition
}
  \item{unmeth1}{
        a  reads count matrix for input samples of one condition
}
  \item{unmeth2}{
        a  reads count matrix for input samples of two condition
}
  \item{mode}{
There are four ways how the empirical dispersion can be computed:

\itemize{
  \item{pooled -}{
Use the samples from all conditions with replicates to estimate a single pooled empirical dispersion value, called "pooled", and assign it to all samples.
}
  \item{per-condition -}{
For each condition with replicates, compute an empirical dispersion value by considering the data from samples for this condition. 
}
  \item{blind -}{
Ignore the sample labels and compute an empirical dispersion value as if all samples were replicates of a single condition. This can be done even if there are no biological 
replicates.}
  \item{auto -}{
Select mode according to the size of samples automaticly.The default is auto.
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
Results will get a matrix including 4 columns (pvalue, log2.fc, q, FDR).

  \item{pvalue}{
Indicate the significance of the methylation site as an RNA differential methylation site
}
  \item{log2.fc}{
  The log2 odds ratio between IP sample and Input sample.
  }
  \item{q}{
  The standardized feature abundance, which is proportional to the expression level of the RNA transcript.
  }
  \item{FDR}{
  FDR of the methylation site, indicating the significance of the peak as an RNA differential methylation site after multiple hypothesis correction usting BH method.  
  }
%%  ~Describe the value returned
}

\author{
Lian Liu <liulian19860905@163.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## The function is currently defined as load library and specify the parameters

library("QNB")

f1 = system.file("extdata", "meth1.txt", package="QNB")
f2 = system.file("extdata", "meth2.txt", package="QNB")
f3 = system.file("extdata", "unmeth1.txt", package="QNB")
f4 = system.file("extdata", "unmeth2.txt", package="QNB")

meth1 = read.table(f1,header=TRUE)
meth2 = read.table(f2,header=TRUE)
unmeth1 = read.table(f3,header=TRUE)
unmeth2 = read.table(f4,header=TRUE)

# When there are replicates under two conditions, we could select 
# \usage{mode="per-condition"} or \usage{mode="pooled"} to estimate the dispersion. 
# The default is \usage{mode="auto"}.

result = qnbtest(meth1,meth2,unmeth1,unmeth2,mode="per-condition")

# If you have replicates for one condition but not for the other, or without any 
# replicates for # two conditions, you can select \usage{mode="blind"} to estimate 
# the dispersion. 

f1 = system.file("extdata", "no_rep_meth1.txt", package="QNB")
f2 = system.file("extdata", "no_rep_meth2.txt", package="QNB")
f3 = system.file("extdata", "no_rep_unmeth1.txt", package="QNB")
f4 = system.file("extdata", "no_rep_unmeth2.txt", package="QNB")

no_rep_meth1 = read.table(f1,header=TRUE)
no_rep_meth2 = read.table(f2,header=TRUE)
no_rep_unmeth1 = read.table(f3,header=TRUE)
no_rep_unmeth2 = read.table(f4,header=TRUE)
head(no_rep_meth1)
head(no_rep_unmeth1)

## Not run: 
result = qnbtest(no_rep_meth1, 
                 no_rep_meth2,
                 no_rep_unmeth1,
                 no_rep_unmeth2,
                 mode="blind")
## End(Not run)

# If you could not decide which mode to estimate dispersion, \usage{mode="auto"} 
# will select suitable way to estimate dispersion according to the replicates.

## Not run: 
result = qnbtest(meth1, meth2,unmeth1,unmeth2)
## End(Not run)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ differential RNA methylation }
\keyword{ quad-negative-binomial model }% __ONLY ONE__ keyword per line