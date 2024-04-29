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
