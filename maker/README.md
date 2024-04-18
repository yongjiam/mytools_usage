## Build docker image for maker

### install maker and repeatmodeler in local conda
```
conda create --name maker bioconda::maker bioconda::repeatmodler
conda activate maker
conda env export > environment.yml
```
### create Dockerfile
```
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

# install maker
COPY environment.yml .
RUN conda env update --file environment.yml
#    conda install -y -c bioconda maker repeatmodeler

# Set the working directory
WORKDIR /data

# Set the default command to execute when the container starts
CMD ["/bin/bash"]
```
### build docker image
```
docker build -t yongmaker .
```
### build singularity image
```
docker images ## list images
export SINGULARITY_CACHEDIR=/path/to/tmp/ ##default to ~/.singularity/
singularity build yongmaker.sif docker-daemon://yongmaker:latest ## need to specify "docker-daemon" and the tag "lastest"
```
#### notes:
> docker images in linux are located at /var/lib/docker
> sudo apt clean ## remove some old packages
> 
