
# NeuroNorm <img src="neuro_sticker.png" align="right" width="160" />

NeuroNorm is an R package that to preprocess structural magnetic resonance imaging (MRI) from multiple patients, diseases, scanners and sites. NeuroNorm transformed multiple raw T1-w images in the NIfTI format into preprocessed images comparable across patients, sites and diseases. Neuronorm performs inhomogeneity correction, spatial registration to a template, skull stripping, spatially informed MRI scan (brain segmentation) generation , intensity normalization and intensity adjustment. NeuroNorm comes up as a standard procedure to compare and analyze multiple T1-w scans of different neurodegenerative diseases. 

This package is an extension of the master thesis **Detection and Classification of Neurodegenerative Diseases: A Spatially Informed Bayesian neural Network** which conducts a population-level analysis of neurodegenerative patients.

## Installation

You can install NeuroNorm from github using `devtools`.

``` r
# install.packages("devtools")
devtools::install_github("DavidPayares/neuronorm")
```

`NeuroNorm` relies on many neuroimaging packages: `fslr`, `ANTsr`,  `extrantsr`, `MNITemplate` and `RAVEL`.
The package `fslr` is available on CRAN, and requires FSL to be installed on
your machine; see the [FSL website](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) for installation. 
For `ANTsR`,`extrantsr` and `RAVEL`, it is recommended to install the latest stable version available at the [ANTsR](https://github.com/stnava/ANTsR/releases/), 
[extrantsr](https://github.com/muschellij2/extrantsr/releases/) and [RAVEL](https://github.com/Jfortin1/RAVEL) GitHub pages, respectively. 
For the template space, the MNI152 atlas with a isomorfic voxel size of 1mm is used. It is included in the `MNITemplate` package, available on GitHub at <https://github.com/Jfortin1/MNITemplate>. 

## Usage

### Data extructure

For using `NeuroNorm`, data must follow a specific structure. This makes easier and more intuitive the loading of input MRI scans and organization of output MRI files. MRI images must be in `NiFTI` format. Currently, `NeuroNorm` only supports T1-w sequence scans. However, other modalities will be implemented in future versions. It is recommended to store your data in the following structure:

```r
├── General_folder              # main folder
│   ├── disease01_patient01     # patient-level folder
│   │   ├── T1-w                # image in NiFTI format
│   │   ├── T2-w
│   ├── disease01_patient02
│   │   ├── T1-w
│   │   ├── T2-w
│   ├── disease02_patient01
│   │   ├── T1-w
│   │   ├── T2-w
│   │   ├── FLAIR
│   ├── disease02_patient02
└── │   ├── T1-w
```

### NeuroNorm preprocessing

### Preprocessed images

## References