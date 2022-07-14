TranScan software presented in the paper entitled "Translocation Detection from Hi-C data via Scan Statistics" by Anthony Cheng, Disheng Mao, Yuping Zhang, Joseph Glaz, and Zhengqing Ouyang.

The TranScan R package detects chromosome translocation events in genome-wide proximity ligation data such as Hi-C. The accompanying translocation breakpoint finder is post-processing scripts that are written in Python.

# Software Requirements

## OS Requirements
This package has been tested on the following systems:
+ macOS Sierra (v 10.12.6)
+ Linux: CentOS 7

## Dependencies and Prerequisties
The following prerequisites are required to run TranScan
+ R (3.5.1)
+ ggplot2 (3.1.0)
+ dplyr (0.7.8)
+ ggpubr (0.2.5)
+ MASS (7.3-51.1)

The following prerequisites are required to run the breakpoint finder
+ Python (2.7.14)
+ matplotlib (2.1.0)
+ numpy (1.15.2)
+ seaborn (0.9.0)

### Optimized Performance (Optional)
For an optimized performance, consider installing OpenBLAS or ATLAS and export
the following environmental variable that points to the appropriate shared object.
There is an approximately 5 fold improvement in the runtime performance as compared
to the default R blas implementation: libRblas.

### OpenBLAS
```
export LD_PRELOAD=</path/to/libopenblas.so>
```

### ATLAS
```
export LD_PRELOAD=</path/to/libatlas.so>
```

# Installation

```
install.packages("devtools")
devtools::install("TranScan")
library(TranScan)
```

# Examples/Demos
Reproducible code examples can be found in the following HTML renders of the jupyter notebooks
+ Notebook01_Run_TranScan.html

# License
The implementations written for this project is covered by the GNU General Public License, version 3.0 (GPL-3.0).
