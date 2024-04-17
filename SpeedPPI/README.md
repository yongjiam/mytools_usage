## SpeedPPI computer-based yeast-2-hybridization
https://github.com/patrickbryant1/SpeedPPI

## setup
```bash
#Get requirements
#conda install -c conda-forge -c bioconda hhsuite
#git clone https://github.com/soedinglab/hh-suite.git
#mkdir -p hh-suite/build && cd hh-suite/build
#cmake -DCMAKE_INSTALL_PREFIX=. ..
#make -j 4 && make install

#Download  Uniclust30
echo "Getting unlclust30..."
wget http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz --no-check-certificate
mkdir data/uniclust30
mv uniclust30_2018_08_hhsuite.tar.gz data/uniclust30
tar -zxvf data/uniclust30/uniclust30_2018_08_hhsuite.tar.gz

#Get the parameters for the Pfam annotation
# Get a TensorFlow SavedModel
wget -qN https://storage.googleapis.com/brain-genomics-public/research/proteins/pfam/models/single_domain_per_sequence_zipped_models/seed_random_32.0/5356760.tar.gz
# unzip
mv 5356760.tar.gz src/domain_mapping/
tar -xzf src/domain_mapping/5356760.tar.gz -C src/domain_mapping/

#Download AF2 parameters
echo "Getting AlphaFold parameters (v. 2021.07.14)..."
mkdir data/params
wget https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar
mv alphafold_params_2021-07-14.tar data/params
tar -xf data/params/alphafold_params_2021-07-14.tar
mv params_model_1.npz data/params

#Cleanup
echo "Cleaning up unnecessary files..."
rm data/uniclust30/uniclust30_2018_08_hhsuite.tar.gz
rm src/domain_mapping/5356760.tar.gz
rm data/params/alphafold_params_2021-07-14.tar
rm params_*.npz
```
## run
```
bash create_ppi_some_vs_some.sh  chunshen_query1.fasta chunshen_query2.fasta  ./hh-suite/bin/hhblits  0.5 ./chunshen_output2/
```
### Notes:
> nvidia-smi  ## check GPU memory \
> free -h ## check CPU memory \
> 2vs3 proteins use more memory than 1vs3
