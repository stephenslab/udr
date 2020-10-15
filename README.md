# udr: Ultimate Deconvolution in R

R package implementing empirical Bayes methods for multivariate normal
means, building on the [extreme deconvolution][ed] method.

## Quick Start

Install the package and build the vignette:

```R
devtools::install_github("stephenslab/udr",build_vignettes = TRUE)
```

Note that installing the package will require a C++ compiler setup
that is appropriate for the version of R installed on your
computer. For details, refer to the [CRAN documentation][cran].

Once you have installed the package, load it into your environment:

```R
library(udr)
```

Work through the [introductory vignette][intro-vignette] illustrating
the methods on a small simulated data set:

```R
vignette("udr_intro")
```

## Credits

The udr R package was developed by [Yunqi Yang][yunqi] and
[Peter Carbonetto][peter] at the [University of Chicago][uchicago],
with contributions from [Matthew Stephens][matthew].

[ed]: https://github.com/jobovy/extreme-deconvolution
[cran]: https://cran.r-project.org
[uchicago]: https://www.uchicago.edu
[yunqi]: https://github.com/Nicholeyang0215
[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
[intro-vignette]: https://stephenslab.github.io/udr/articles/udr_intro.html
