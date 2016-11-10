\name{QNB-package}
\alias{QNB-package}
\alias{QNB}
\docType{package}
\title{
Differential RNA Methylation Analysis with Quad-Negative-Binomial Model
}
\description{
With high throughput sequencing techniques, such as MeRIP-Seq, transcriptome-wide RNA methylation profile is now available as count-based measurements, with which it is often of special interests to study the dynamics in post-transcriptional regulation via differential RNA methylation analysis. In differential RNA methylation analysis, the sample size of RNA methylation experiment is usually very small due to its cost, making the inference very difficult. 
Meanwhile, there exist a large number of very lowly expressed genes whose methylation level 
cannot be accurately estimated due to limited sequencing depth.We present QNB, a statistical 
approach for differential RNA methylation analysis with count-based small-sample sequencing 
data. The method is based on 4 independent negative binomial dis-tributions, and with their 
variances and means linked by local regressions, respectively. And the algorithm can be in 
theory applied to data generated from other sequencing techniques related RNA modifications, 
including but not limited to RNA BS-Seq, Par-CLIP, RIP-Seq, etc.Please don't hesitate to contact <liulian19860905@163.com> if you have any questions.
}
\details{
\tabular{ll}{
Package: \tab QNB\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2016-11-5\cr
License: \tab GPL-2\cr
}
}
\author{
Lian Liu
Maintainer: <liulian19860905@163.com>
}
\keyword{ differential RNA methylation, quad-negative-binomial model }

\examples{
library(QNB)
?QNB
}