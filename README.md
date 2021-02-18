
# NeuroNorm <img src="neuro_sticker.png" align="right" width="160" />

NeuroNorm is an R package that to preprocess structural magnetic resonance imaging (MRI) from multiple patients, diseases, scanners and sites. NeuroNorm transformed multiple raw T1-w images in the NIfTI format into preprocessed images comparable across patients, sites and diseases. Neuronorm performs inhomogeneity correction, spatial registration to a template, skull stripping, spatially informed MRI scan (brain segmentation) generation , intensity normalization and intensity adjustment. NeuroNorm comes up as a standard procedure to compare and analyze multiple T1-w scans of different neurodegenerative diseases. 

This package is an extension of the master thesis **Detection and Classification of Neurodegenerative Diseases: A Spatially Informed Bayesian neural Network** which conducts a population-level analysis of neurodegenerative patients.

## Installation

You can install NeuroNorm from github using `devtools`.

``` r
# install.packages("devtools")
devtools::install_github("DavidPayares/neuronorm")
```

`Neuronorm` relies on many neuroimaging packages: `fslr`, `ANTsr`,  `extrantsr`, `MNITemplate` and `RAVEL`.
The package `fslr` is available on CRAN, and requires FSL to be installed on
your machine; see the [FSL website](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) for installation. 
For `ANTsR`,`extrantsr` and `RAVEL`, it is recommended to install the latest stable version available at the [ANTsR GitHub page](https://github.com/stnava/ANTsR/releases/), 
[extrantsr GitHub page](https://github.com/muschellij2/extrantsr/releases/) and [RAVEL GitHub page](https://github.com/Jfortin1/RAVEL), respectively. 
For the template space, we use the MNI152 atlas with a isomorfic voxel size of 1mm included in the `MNITemplate` package, available on GitHub at <https://github.com/Jfortin1/MNITemplate>. 