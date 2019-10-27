# Conformal Inference R Project
## Maintained by Ryan Tibshirani
### Based on work by Rina Barber, Emmanuel Candes, Max G'Sell, Jing Lei, Aaditya Ramdas, Alessandro Rinaldo, Ryan Tibshirani, Larry Wasserman

This repository contains R software tools for conformal inference. The current
emphasis is on conformal prediction in regression. We may eventually add tools
for density estimation and classification.

The folder "conformalInference" can be installed as an R package, providing
access to the software tools, and the file "conformalInference.pdf" contains 
documentation. 

The folder "lei2018" contains R code to reproduce all examples in the paper
[Distribution-Free Predictive Inference for Regression](http://www.stat.cmu.edu/~ryantibs/papers/conformal.pdf)
by Lei, G'Sell, Rinaldo, Tibshirani, Wasserman (2018).  The folder
"tibshirani2019" contains R code to reproduce all examples in the paper
[Conformal Prediction Under Covariate Shift](http://www.stat.cmu.edu/~ryantibs/papers/weightedcp.pdf)
by Tibshirani, Barber, Candes, Ramdas (2019).  This code all relies on the
"conformalInference" R package.

Relevant work (in reverse chronological order):

- [Conformal Prediction Under Covariate Shift](http://www.stat.cmu.edu/~ryantibs/papers/weightedcp.pdf)
  by Ryan Tibshirani, Rina Barber, Emmanuel Candes, Aaditya Ramdas (2019). 
- [Distribution-Free Predictive Inference for Regression](http://www.stat.cmu.edu/~ryantibs/papers/conformal.pdf) 
  by Jing Lei, Max G'Sell, Alessandro Rinaldo, Ryan Tibshirani, and Larry
  Wasserman, Journal of the American Statistical Association, 113(523),
  1094-1111, 2018.
- [Distribution-Free Prediction Bands for Non-parametric Regression](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12021)
  by Jing Lei and Larry Wasserman, Journal of the Royal Statistical Society:
  Series B, 76(1), 71-96, 2014. 
- [A Conformal Prediction Approach to Explore Functional Data](https://link.springer.com/article/10.1007%2Fs10472-013-9366-6)
  by Jing Lei, Alessandro Rinaldo, and Larry Wasserman, Annals of Mathematics
  and Artificial Intelligence, 74(4), 29-43, 2013.
- [Distribution Free Prediction Sets](https://amstat.tandfonline.com/doi/abs/10.1080/01621459.2012.751873)
  by Jing Lei, James Robins, and Larry Wasserman, Journal of the American
  Statistical Association, 108(501), 278-287, 2013.
- [On-line Predictive Linear Regression](http://www.alrw.net/articles/01.pdf) by
  Vladimir Vovk, Ilia Nouretdinov, and Alex Gammerman, Annals of Statistics,
  37(3), 1566-1590, 2009. 
- [Algorithmic Learning in a Random World](http://www.alrw.net) by Vladimir
  Vovk, Alex Gammerman, and Glenn Shafer, Springer, 2005.
  
### Install the R package

To install the conformalInference R package directly from github, run the
following in R:

```{r}
library(devtools)
install_github(repo="ryantibs/conformal", subdir="conformalInference")
```
