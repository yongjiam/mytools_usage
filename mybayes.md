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
## Step 1. prepare codon-based sequence alignment file in nexus format
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
#### single_copy_gene.mb (create a mb file containing all commands)
```text
begin mrbayes;
  set autoclose=yes nowarn=yes;
  execute Group2.nex;
  lset applyto=(all) nst=6 rates=invgamma;
  charset OG0015498 = 1-1044;
  charset OG0015994 = 1045-2910;
  charset OG0015012 = 2911-3921;
  charset OG0015559 = 3922-6024;
  charset OG0015340 = 6025-7671;
  charset OG0015641 = 7672-9297;
  charset OG0015336 = 9298-10617;
  charset OG0014911 = 10618-13833;
  charset OG0014703 = 13834-15648;
  charset OG0014957 = 15649-16839;
  charset OG0015189 = 16840-17994;
  charset OG0015919 = 17995-21063;
  charset OG0016312 = 21064-22914;
  charset OG0016172 = 22915-24195;
  charset OG0015468 = 24196-25893;
  charset OG0015178 = 25894-29079;
  charset OG0015300 = 29080-31719;
  charset OG0015822 = 31720-32937;
  charset OG0015838 = 32938-34311;
  charset OG0015617 = 34312-35247;
  charset OG0016226 = 35248-37494;
  charset OG0015146 = 37495-38850;
  charset OG0015277 = 38851-40617;
  charset OG0014827 = 40618-42759;
  charset OG0015981 = 42760-43569;
  charset OG0016104 = 43570-45231;
  charset OG0015954 = 45232-47193;
  charset OG0015053 = 47194-49944;
  charset OG0016024 = 49945-51771;
  charset OG0016350 = 51772-52824;
  charset OG0015711 = 52825-54078;
  charset OG0016222 = 54079-57021;
  charset OG0015633 = 57022-59232;
  charset OG0015083 = 59233-60537;
  charset OG0014713 = 60538-62349;
  charset OG0016180 = 62350-64179;
  charset OG0015792 = 64180-66696;
  charset OG0015853 = 66697-67521;
  charset OG0016255 = 67522-68463;
  charset OG0016050 = 68464-70629;
  charset OG0015934 = 70630-73248;
  charset OG0016174 = 73249-74607;
  charset OG0016184 = 74608-75333;
  charset OG0015903 = 75334-77241;
  charset OG0015602 = 77242-78255;
  charset OG0015255 = 78256-79632;
  charset OG0014690 = 79633-80415;
  charset OG0016234 = 80416-81711;
  charset OG0015852 = 81712-84153;
  charset OG0015298 = 84154-86643;
  charset OG0015486 = 86644-88230;
  charset OG0015282 = 88231-91380;
  charset OG0014815 = 91381-91941;
  charset OG0014809 = 91942-93453;
  partition Genes = 54:OG0015498,OG0015994,OG0015012,OG0015559,OG0015340,OG0015641,OG0015336,OG0014911,OG0014703,OG0014957,OG0015189,OG0015919,OG0016312,OG0016172,OG0015468,OG0015178,OG0015300,OG0015822,OG0015838,OG0015617,OG0016226,OG0015146,OG0015277,OG0014827,OG0015981,OG0016104,OG0015954,OG0015053,OG0016024,OG0016350,OG0015711,OG0016222,OG0015633,OG0015083,OG0014713,OG0016180,OG0015792,OG0015853,OG0016255,OG0016050,OG0015934,OG0016174,OG0016184,OG0015903,OG0015602,OG0015255,OG0014690,OG0016234,OG0015852,OG0015298,OG0015486,OG0015282,OG0014815,OG0014809;
  set partition=Genes;
  mcmc ngen=50000 samplefreq=100 printfreq=100 nchains=2 nruns=2 temp=0.2;
  sump nruns=2;
  sumt nruns=2;
end;
```

#### snp_annotated2.mb
```text
begin mrbayes;
  set autoclose=yes nowarn=yes;
  execute snp_annoted2.nex;
  lset nst=6 rates=invgamma;
  mcmc ngen=150000 samplefreq=100 printfreq=100 nchains=4 nruns=2 temp=0.2;
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

# Build docker image for mrbayes and run on sever
#### Dockerfile
```text
# Use the official Ubuntu 20.04 LTS image as the base image
FROM ubuntu:20.04

# Set the working directory to /app
WORKDIR /app

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Australia/Perth

# Update the package list and install any needed dependencies
RUN apt-get update && apt-get install -y \
    # Add your dependencies here, for example:
    --no-install-recommends tzdata \
    build-essential \
    libreadline8 \
    curl \
    mpich

    # other-dependency \

# Copy the current directory contents into the container at /app
COPY . /app

# Set environment variables if needed
ENV PATH="/app/src:${PATH}"

# Specify the default command to run when the container starts
CMD ["/bin/bash"]
```
### notes:
```bash
lsb_release -a ## get linux system version
ldd /data/tools/MrBayes/src/mb ## get requried packages for the mb command
```
### build, test, run docker image
```bash
docker build -t yongmrbayes .  ## build image
docker images  ## list images
docker run -it --rm yongmrbayes which mb ## check mb or mpirun available
docker run -it --rm -v /data/igenome/single-copy-OG/mrbayes:/data yongmrbayes /bin/bash -c "cd /data && mb test2.mb" ## can't run mb if binding data to $pwd
docker run -it --rm -v /data/igenome/single-copy-OG/mrbayes:/data yongmrbayes /bin/bash -c "cd /data && mpirun -np 3 mb test2.mb"
```
### build singularity container from docker image
```bash
docker images ## list images
singularity build yongmrbayes.sif docker-daemon://yongmrbayes:latest ## need to specify "docker-daemon" and the tag "lastest"
singularity run -B /data/igenome/single-copy-OG/mrbayes:/data yongmrbayes.sif /bin/bash -c "cd /data && mb test2.mb" ## test run using singularity
```
