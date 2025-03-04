# r2tfdiffusion

**r2tfdiffusion** is an R package for modeling **network diffusiovity between receptors and transcription factors** using ex vivo data. to train the networks.


## 📥 Installation

You can install the package directly from GitHub using:
```r
install.packages("devtools")  # Install devtools if you haven't already
devtools::install_github("endreanderssen/r2tfdiffusion")

The package relies on paralell processing using the doParallel framework. Set up workers up to one per sample. If you are unfamilliar with doParallel setup is exlpained in the vignette. 

After setting up
```r
library(r2tfdiffusion)
data(r2tfData)
diff_result <- diffusionMap2(
    receptors = r2tfData$receptors,
    TFs = c('JUN','RELA'),
    feedbackNodes = r2tfData$feedbackNodes,
    M = r2tfDATA$Mexvivo,
    network = r2tfData$network,
    inputSignal = 0.99,
    nTicks = 2000
)
```

for a complete tutorial, check the vignette:
```r
browseVignettes("r2tfdiffusion")
```

This package is licensed under the MIT License.
