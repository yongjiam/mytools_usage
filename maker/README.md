## change docker image storage location in ubuntu
https://www.guguweb.com/2019/02/07/how-to-move-docker-data-directory-to-another-location-on-ubuntu/
```bash
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

## if /home drive out of space, build in another computer
docker save -o yongmaker_docker.tar yongmaker:latest
singularity build yongmaker.sif docker-archive://yongmaker_docker.tar

## push local docker image to docker hub
docker login
docker info|grep Username
docker tag yongmaker yongjia111/yongmaker:latest ### use docker hub account name yongjia111
docker push yongjia111/yongmaker:latest

SINGULARITY_CACHEDIR
singularity build yongmaker.sif docker://yongjia111/yongmaker:latest

## run docker image interactively
docker run -it -v ${PWD}:/data yongmaker:latest bash
source activate /path/to/conda/maker/env
exit
docker ps -l
docker commit b16d1da08f93 yongmaker:latest
docker commit b16d1da08f93 yongjia111/yongmaker:latest

```
#### notes:
> docker images in linux are located at /var/lib/docker
> sudo apt clean ## remove some old packages
> https://hub.docker.com/repositories/yongjia111
