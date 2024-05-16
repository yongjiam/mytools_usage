## conda init and conda.sh
```
conda init
. /scratch/pawsey0399/yjia/tools/miniconda3/etc/profile.d/conda.sh
```
## fixes mamba installation for conda
after installing mumba, the default conda envs_dirs changed from /data/tools/miniconda3/envs/ to /data/tools/miniforge3/envs/
Thus, environments located in /data/tools/miniconda3/envs/ could not activated by conda activate name, the fix is:
```bash
conda config --append envs_dirs /data/tools/miniconda3/envs
```
## install mamba to conda base env
```
conda install -n base mamba
```
## conda module not found error
https://github.com/conda/conda/issues/9672
```bash
# execute as root in /opt/conda/bin/
unlink python
ln -s python3.5 python
# you should now have python -> python3.5 
```
<img src="setonix_conda_bin.png" alt="miniconda3 python" width="600">
