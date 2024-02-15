## computer-based protein protein interaction simulation
https://github.com/patrickbryant1/SpeedPPI

## installation and database setup
```bash
## Python packages
git clone https://github.com/patrickbryant1/SpeedPPI
cd SpeedPPI
conda env create -f speed_ppi.yml
wait
conda activate speed_ppi
pip install --upgrade "jax[cuda12_local]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

conda activate speed_ppi
python3 ./src/test_gpu_avail.py

## HHblits
mkdir hh-suite
cd hh-suite
wget https://github.com/soedinglab/hh-suite/releases/download/v3.3.0/hhsuite-3.3.0-SSE2-Linux.tar.gz
tar xvfz hhsuite-3.3.0-SSE2-Linux.tar.gz
cd ..

## Uniclust30
wget http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz --no-check-certificate
tar -zxvf uniclust30_2018_08_hhsuite.tar.gz -C data/

## AlphaFold2 parameters
mkdir data/params
wget https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar
mv alphafold_params_2021-07-14.tar data/params
tar -xf data/params/alphafold_params_2021-07-14.tar
mv params_model_1.npz data/params

## Cleanup - remove unnecessary files
rm uniclust30_2018_08_hhsuite.tar.gz
rm data/params/alphafold_params_2021-07-14.tar
rm params_*.npz
rm hh-suite/hhsuite*.tar.gz
```
### Run the pipeline
```bash
## tests
conda activate speed_ppi
bash create_ppi_all_vs_all.sh ./data/test/test.fasta ./hh-suite/bin/hhblits 0.5 ./data/test/all_vs_all/

## some vs some
conda activate speed_ppi
bash create_ppi_some_vs_some.sh ./data/test/test1.fasta ./data/test/test2.fasta ./hh-suite/bin/hhblits 0.5 ./data/test/some_vs_some/
```
