## troubleshooting r package loading issue in jupyter-lab
### make sure both jupyterlab and jupyter are installed
```bash
conda create --name myjupyter python=3.9
conda install -c anaconda jupyter
conda install -c conda-forge jupyterlab
```
### within R, run:
```R
install.packages("IRkernel")
IRkernel::installspec()
```
### check which R, and conda env
```python
Sys.which("R")
Sys.getenv("CONDA_DEFAULT_ENV")
```
