## notes on docker usage
docker images\
docker ps\
docker ps -l\
docker rmi -f image_id\

### build docker image, mrbayes
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
### build docker image, miniconda
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