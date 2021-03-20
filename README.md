# pybragg
An implementation of the [Bortfeld bragg curve fit](https://pubmed.ncbi.nlm.nih.gov/9434986/) in python with scipy and numpy.

## Installation

You can install `pybragg` using `pip install <optional: '--user'> git+https://github.com/flome/pybragg.git`.

## Usage

If you have a data set `(z, D)` with two arrays, `z` being the depth in mm in a phantom and `D`being the deposited dose at that depth, you can fit the curve using the bortfeld fit with
```
from pybragg import fitBP
bortfeld_fit_params = fitBP(z, D)
```

The fit returns characteristic values like `D100` (max. Dose) or the proton beam range `R80D` (distal side of the beam, 80\% maximum dose).
