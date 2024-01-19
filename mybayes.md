## Bayesian phylogeny using Mrbayes
#### resources and tutorials
http://treethinkers.org/tutorials/mrbayes/
https://github.com/fhcrc/seqmagick ## convert sequence file format
https://biopython.org/wiki/Concatenate_nexus ## concatenate multiple genes alignment

## installation
```bash
## build mrbayes
git clone https://github.com/NBISweden/MrBayes
cd MrBayes
./configure
make && sudo make install

## install mpi tools
sudo apt-get install -y mpich

```
## best substitution model selection using jmodeltest
```bash
wget https://github.com/ddarriba/jmodeltest2/files/157117/jmodeltest-2.1.10.tar.gz
tar -xzvf jmodeltest-2.1.10.tar.gz
cd jmodeltest-2.1.10/
java -jar jModelTest.jar -d mydata/combined50.nex -g 4 -i -f -AIC -BIC -a > log.txt
```

## run
## Step 1. prepare sequence alignment file in nexus format
```python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import nt_search
from Bio.Align.Applications import MafftCommandline
import os
import glob
import random

## loop through multiple fasta and perform codon-based alignment

def translate_cds_to_aa(cds_sequence):
    return cds_sequence.translate()

def map_gaps_to_cds(amino_acid_alignment, original_cds_sequences):
    aligned_cds_sequences = {}

    for aa_record in amino_acid_alignment:
        cds_record = next(cds for cds in original_cds_sequences if cds.id == aa_record.id)

        aligned_cds_seq = ""
        aa_pos = 0

        for aa in str(aa_record.seq):
            if aa == "-":
                aligned_cds_seq += "---"
            else:
                aligned_cds_seq += str(cds_record.seq[aa_pos * 3 : (aa_pos + 1) * 3])
                aa_pos += 1

        aligned_cds_record = SeqRecord(Seq(aligned_cds_seq), id=cds_record.id, description="")
        aligned_cds_sequences[aligned_cds_record.id] = aligned_cds_record

    return list(aligned_cds_sequences.values())

def perform_cds_alignment(input_cds_file, output_cds_alignment_file):
    # Step 1: Translate CDS to Amino Acid Sequences
    amino_acid_sequences = []

    for record in SeqIO.parse(input_cds_file, "fasta"):
        amino_acid_sequence = translate_cds_to_aa(record.seq)
        amino_acid_record = SeqRecord(amino_acid_sequence, id=record.id, description="")
        amino_acid_sequences.append(amino_acid_record)

    SeqIO.write(amino_acid_sequences, "amino_acid_sequences.fasta", "fasta")

    # Step 2: Amino Acid Sequence Alignment
    mafft_cline = MafftCommandline(input="amino_acid_sequences.fasta", auto=True)
    stdout, stderr = mafft_cline()

    with open("amino_acid_alignment.fasta", "w") as f:
        f.write(stdout)

    # Step 3: Convert Amino Acid Alignment to CDS Alignment
    amino_acid_alignment = list(SeqIO.parse("amino_acid_alignment.fasta", "fasta"))
    original_cds_sequences = list(SeqIO.parse(input_cds_file, "fasta"))

    aligned_cds_sequences = map_gaps_to_cds(amino_acid_alignment, original_cds_sequences)

    # Step 4: Write the final CDS alignment to a file
    SeqIO.write(aligned_cds_sequences, output_cds_alignment_file, "fasta")

if __name__ == "__main__":

    # List of FASTA files
    fasta_files = glob.glob("renamed_*.fasta")

    for input_cds_file in fasta_files:
        output_cds_alignment_file = input_cds_file[:-2]
        perform_cds_alignment(input_cds_file, output_cds_alignment_file)

## convert fasta format to nexus format
## install seqmagick
!pip install seqmagick
!ls *.fas |while read R;do seqmagick convert --output-format nexus --alphabet dna $R $R".nex";done
```
## Step 2. concatenate multiple gene alignment into a single nexus and partition on genes
```python
# Specify the desired directory path
new_directory = '/data/igenome/single-copy-OG/renamed_aligned_sequences50'

# Change the working directory
os.chdir(new_directory)

# the combine function takes a list of tuples [(name, nexus instance)...],
# if we provide the file names in a list we can use a list comprehension to
# create these tuples
from Bio.Nexus import Nexus
!rm combined.nex
nexi = []
file_list = glob.glob("*.nex")
file_list2 = file_list[0:50] ## random50 = random.sample(file_list, 50)
nexi = [(fname, Nexus.Nexus(fname)) for fname in file_list]

combined = Nexus.combine(nexi)
with open("combined.nex", "w") as f:
    combined.write_nexus_data(filename=f)
```
## Step 3. create mrbayes configuration file and run mrbayes
https://wiki.rc.usf.edu/index.php/MrBayes
#### sample.mb (create a mb file containing all commands)
```text
# My MrBayes input file

# Specify the data file
execute my_alignment.nex ## loads your Nexus format alignment file

# Set the outgroup (if needed)
outgroup taxon_name

# Set the number of generations and sample frequency
mcmc ngen=1000000 samplefreq=1000  

# Specify the substitution model
lset nst=6 rates=invgamma

# Set the number of chains and heating parameters
mcmc nchains=4 temp=0.2 0.2 0.2 0.2

# Run the analysis
mcmc
```

#### test.mb
```text
begin mrbayes;
  set autoclose=yes nowarn=yes;
  execute snp_annoted2.nex;
  lset nst=6 rates=invgamma;
  mcmc ngen=10000 samplefreq=100 printfreq=100 nchains=4 nruns=2 temp=0.2;
  mcmc;
  sump nruns=2;
  sumt nruns=2;
end;
```

## Step 4. batch mode: run mrbayes using mpirun
```bash
## non-mpi mode:
./mb your_file.mb

## batch mode
mpirun -np 30 mb test.mb > log.txt
```
