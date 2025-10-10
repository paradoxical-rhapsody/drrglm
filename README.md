
# drrglm

The `drrglm` package implements doubly regularized matrix-variate regression for generalized linear models. 
Additionally, it can be used for robust matrix factorization (such as robust PCA or factor models).



## Installation

The `drrglm` is currently on its way to CRAN.
For now, it can be installed from Github repository:
```R
library(remotes)
install_github("paradoxical-rhapsody/drrglm")
```

**Note**: On Windows platform, it requires [RTools](https://cran.r-project.org/bin/windows/Rtools/) to compile the source code.


## Usage

For detailed information on the usage and examples in the package, please refer to the help documentation:
```R
help(, drrglm)
```



## About `EEG` data

The `drrglm` provides preprocessed EEG data, which can be loaded via
```R
data(EEG, package="drrglm")
```


### Data Source
The original data are fully open-access and available at [UCI Machine Learning Repository](https://kdd.ics.uci.edu/databases/eeg/eeg.data.html).
The dataset includes two groups of subjects: 77 alcoholic and 45 control. 
In the original study, each subject was exposed to either a single stimulus (S1) or two stimuli (S1 and S2), which were pictures chosen from the 1980 Snodgrass and Vanderwart picture set. 
Each subject underwent 120 trials under each condition.


Here in this processed dataset, we provide the averages of 120 trials under S1 condition. 
The dataset is structured as a list containg two arrays,
- `EEG$alcholic`: Array of dimensions `256 x 64 x 77`.
- `EEG$control`:  Array of dimensions `256 x 64 x 45`.


### Data Citation

If you use this dataset, we would be very grateful if you could cite both the original data source and our work:

> Xu Z., Luo S., and Jiang B. [Doubly Regularized Matrix-Variate Regression](), Submitted.
