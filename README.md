## Granger Components Analysis
### Source code

This repository contains the source code for the paper "Granger Components Analysis: Unsupervised learning of latent temporal dependences". The paper is available on [openreview](https://openreview.net/forum?id=wqIm0Qsgy0).

### MATLAB

A demonstration of how to apply GCA to simulated VAR(3) data is provided:
- as a live script [here](code/matlab/demo_gca.mlx)
  - the output of the live script is [here](code/matlab/demo_gca.pdf)
- as a standard script [here](code/matlab/demo_gca.m)

The core function is [runGcaTrAlt.m](code/matlab/runGcaTrAlt.m).


### Python

_Disclaimer: this is a work-in-progress._ The MATLAB code has been extensively tested and is the recommended implementation.

A demonstration of how to apply GCA to simulated VAR(3) data is provided as a Jupyter notebook [here](code/python/notebooks/demo_gca.ipynb)

The core functions required to implement GCA in Python are [gca.py](code/python/gca.py).
