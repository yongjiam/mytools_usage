## github website
https://github.com/tanghaibao/jcvi

#### install jcvi
# docker
docker pull tanghaibao/jcvi ## pull docker images
docker images ## list images

docker run your-image-name commands ## run a basic container
docker run --name my-container-name your-image-name ## assign a name to the container for easier management
docker run -it your-image-name commands ## run the container in interactive mode
docker run -d your-image-name ## run the container in detached model (in the background)
docker run -it -v host_dir:container_dir your-image-name ## bind a host directory to a container
docker run --rm your-image-name ## removing container on exit

docker ps ## once your container is running, list running containers
docker stop container_id ## stop a container

docker run -it your-image-name bash ## enter a docker image and start bash and run commands interactively
docker run -it -v ${PWD}:/mydata tanghaibao/jcvi bash ##bind current directory to /mydata new directory
