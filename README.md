# pybragg
An implementation of the [Bortfeld bragg curve fit](https://pubmed.ncbi.nlm.nih.gov/9434986/) in python with scipy and numpy.

## Installation

You can install `pybragg` using `pip install <optional: '--user'> git+https://github.com/flome/pybragg.git`.

## Usage

### Bragg peak fit
If you have a data set `(z, D)` with two arrays, `z` being the depth in mm in a phantom and `D`being the deposited dose at that depth, you can fit the curve using the bortfeld fit with
```
from pybragg import fitBP
bortfeld_fit_params = fitBP(z, D)
```

The fit returns characteristic values like `D100` (max. Dose) or the proton beam range `R80D` (distal side of the beam, 80\% maximum dose).

### The bortfeld equation
You can use the bortfeld equation to generate bragg peaks:
```
import numpy
from pybragg import bortfeld

z = np.linspace(0, 200, 200) # depth in mm
D = bortfeld(z, R0=130, sigma=1.3, Phi=600, epsilon=.3)

import matplotlib.pyplot as plt

plt.plot(z, D, label='bortfeld function')
plt.xlabel('$z$ in mm')
plt.ylabel('$D$ in a.u.')
plt.legend(loc='best')
plt.tight_layout()

plt.show()
```
