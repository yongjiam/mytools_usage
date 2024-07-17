## change docker image storage location in ubuntu
https://www.guguweb.com/2019/02/07/how-to-move-docker-data-directory-to-another-location-on-ubuntu/
```
sudo service docker stop

mkdir /data/tools/docker
sudo nano /etc/docker/daemon.json
## /etc/docker/daemon.json
{
  "data-root": "/data/tools/docker"
}
##
sudo rsync -aP /var/lib/docker/ /path/to/your/docker
sudo mv /var/lib/docker /var/lib/docker.old
sudo service docker start
sudo rm -rf /var/lib/docker.old
```
## notes on docker usage
```bash
docker images\
docker ps\
docker ps -l\
docker rmi -f image_id\
docker build -t IMAGE_NAME .\
docker run -it --rm -v /data/igenome/single-copy-OG/mrbayes:/data yongmrbayes /bin/bash -c "cd /data && mpirun -np 3 mb test2.mb"\
docker stop container_id

## other notes
lsb_release -a ## get linux system version\
ldd /data/tools/MrBayes/src/mb ## get requried packages for the mb command
```
## modify a docker image
```
## pull juicer docker image
docker run aidenlab/juicer:latest
docker images

## run docker image interactively
docker run -it -v ${PWD}:/data --entrypoint=/bin/bash aidenlab/juicer:latest ## use entrypoint to stop automatic run

## make your changes
cd /data ## ln softlink not recoganized in docker mount, have to copy files
mkdir /aidenlab && cd /aidenlab
ln -s /data/tools/juicer/CPU scripts ## the program looks for /aidenlab/scripts/common/countligations.sh: No such file or directory

# Exit the container when you're done
exit

# Get the ID of the container you ran
docker ps -l ## get ID

# Commit the changes to a new Docker image
docker commit b16d1da08f93 aidenlab/juicer:yongjia

## build singularity image from docker
singularity build juicer.sif docker-daemon://aidenlab/juicer:yongjia
```
## push a local docker image to dockerhub
```
## push local docker image to docker hub
docker login
docker info|grep Username
docker tag yongmaker yongjia111/yongmaker:latest ### use docker hub account name yongjia111
docker push yongjia111/yongmaker:latest
singularity build yongmaker.sif docker://yongjia111/yongmaker:latest
```
### Dockerfile, mrbayes
```bash
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
### Dockerfile hic-pipeline
```
# Use a base image with Conda installed
FROM continuumio/miniconda3:latest

# Set the working directory inside the container
WORKDIR /app

# Copy the Conda environment YAML file into the container
COPY hic_pipeline.yml .

# Create the Conda environment based on the YAML file
RUN conda env create -f hic_pipeline.yml

# Activate the Conda environment
SHELL ["conda", "run", "-n", "hic-scaffolding-nf", "/bin/bash", "-c"]

# Any additional commands you need can go here
# For example, you might want to set up some environment variables or copy your application code
# Install additional packages into the Conda environment
RUN conda install -n hic-scaffolding-nf -c conda-forge nextflow

# Set the command to run your application
CMD ["/bin/bash"]
```
### Dockerfile old, miniconda
```bash
##Dockerfile
# Use a base image
FROM ubuntu:20.04

# Install dependencies
RUN apt-get update && \
    apt-get install -y wget && \
    rm -rf /var/lib/apt/lists/*

# Download and install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Set PATH environment variable to include Miniconda
ENV PATH="/opt/miniconda3/bin:${PATH}"
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Australia/Perth

# Create a Conda environment
COPY environment.yml .
#RUN conda env create -f environment.yml
#RUN conda env update --name base --file environment.yml

# Activate the Conda environment
#RUN echo "conda activate $(head -1 environment.yml | cut -d' ' -f2)" > ~/.bashrc

# Set the working directory
WORKDIR /app

# Copy the rest of your application code
COPY . .

# Set the default command to execute when the container starts
CMD ["/bin/bash"]
```
### Dockerfile, rMVP
```
# Use an official R base image from Docker Hub
FROM r-base:latest

# Install dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install rMVP package
RUN R -e "install.packages('rMVP', repos='http://cran.rstudio.com/')"

# Set the working directory
WORKDIR /app

# Start R
CMD ["R"]
```
### Dockerfile, ragtag
conda activate ragtag
conda env export > ragtag.yml
```
# Use the official Miniconda3 image as a base
FROM continuumio/miniconda3:latest

# Set the working directory inside the container
WORKDIR /app

# Copy the environment file to the container
COPY ragtag.yml .

# Create a conda environment with the specified dependencies
RUN conda env create -f ragtag.yml

# Activate the conda environment
SHELL ["conda", "run", "-n", "ragtag", "/bin/bash", "-c"]

# Install Ragtag package in the base environment
#RUN conda install -c conda-forge -c bioconda ragtag

# Set the default command to run when the container starts
CMD ["/bin/bash"]
```
### Dockerfile, earlGrey
```
# Use the Miniconda3 base image
FROM continuumio/miniconda3:latest

# Set the working directory inside the container
WORKDIR /data

# Update conda and install mamba in the base environment
RUN conda install -n base -c conda-forge mamba

# Install the earlgrey package using mamba
RUN mamba create -n earlgrey -c conda-forge -c bioconda earlgrey=4.2.4

# Activate the base environment by default
SHELL ["conda", "run", "-n", "earlgrey", "/bin/bash", "-c"]

# Optional: Set the entrypoint to bash
ENTRYPOINT ["/bin/bash"]
## Note: the above ENTRYPOINT is different from CMD, with ENTRYPOINT, you will run the job like below:
## srun --export=all -n 1 -c 64   singularity run -B ${PWD}:${PWD} docker://yongjia111/myearlgrey:latest earlgrey.sh
## but not: bash earlgrey.sh
```
### singularity, earlGrey
```
mkdir /data/singularity_tmp
export SINGULARITY_TMPDIR=/data/singularity_tmp

singularity build earlgreydfam38.sif docker://tobybaril/earlgrey_dfam3.8:latest

##setonix
module load singularity/3.11.4-slurm
export SINGULARITY_CACHEDIR=/scratch/pawsey0399/yjia/WBT/hifionly_run2/earlgrey/tmp
export IAMGE=/scratch/pawsey0399/yjia/huyou/containers/earlgreydfam38.sif
srun --export=all -n 1 -c 40   singularity exec -B ${PWD}:/data $IMAGE earlgrey.sh

## earlgrey.sh
#!/bin/bash
earlGrey -g /data/hap1.fasta \
        -s hap1 \
        -o /data/ \
        -r eukarya \
	      -d yes \
	      -m yes \
        -t 40
```
