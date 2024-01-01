## install alphafold_pulldown
https://github.com/KosinskiLab/AlphaPulldown
```bash
conda create -n AlphaPulldown -c omnia -c bioconda -c conda-forge python==3.10 openmm==8.0 pdbfixer==1.9 kalign2 cctbx-base pytest importlib_metadata hhsuite
source activate AlphaPulldown
conda install -c bioconda hmmer

source activate AlphaPulldown

python3 -m pip install alphapulldown==1.0.3
pip install jax==0.4.23 jaxlib==0.4.23+cuda11.cudnn86 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

## download alphafold database files

## trouble shooting

### 
pip install tensorflow-gpu==2.8.0
