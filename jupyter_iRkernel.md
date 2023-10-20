## troubleshooting r package loading issue in jupyter-lab
### make sure both jupyterlab and jupyter are installed
conda install -c anaconda jupyter
conda install -c conda-forge jupyterlab

### within R, run:
install.packages("IRkernel")
IRkernel::installspec()

### check which R, and conda env
Sys.which("R")
Sys.getenv("CONDA_DEFAULT_ENV")
