## docker
docker pull csicunam/get_homologues ##pull image form docker hub

docker run -v /tmp:/container/directory imagename ##add/mount shared folder to image

docker run -v /Volumes/Elements5T/Programs/docker_shared/you:/home/you csicunam/get_homologues
docker run --rm -it csicunam/get_homologues:latest /bin/bash ##run docker image

docker run -it -P --name GET_HOMO -v ~/container-data:/data csicunam/get_homologues
  
#list all images created
docker ps -a
docker attach ID(first four digit of container ID)

## quit docker
Exit

#### build docker image from docker file, then create singularity image sif from docker files
>download docker file Dockerfiler, then run:
docker build -t local/get_homologues .
>this will create docker image in docker /local
>you can also pull docker file form docker hub, but the image in docker hub does not contain the hmm and swissprot database, which need to downloaded during installation

> after build docker image at /local, then pull docker image to create singularity sif
singularity pull get_homologues.sif docker-daemon://local/get_homologues:latest
> copy sif file to anywhere to use

##
## samplot
https://github.com/ryanlayer/samplot
### pull samplot or build singularity image
```bash
docker pull dceoy/samplot ## docker
singularity build samplot.sif docker://dceoy/samplot ## singularity
```
### rum samplot 
```bash
##samplot.sh
samplot plot \
    -n HvBlp_mapping \
    -b W1.bam 1-2.sort.bam \
    -o updated3.png \
    -c chr1H \
    -s 518059036 \
    -e 518560254 \
    --coverage_only \
    --hide_annotation_labels \
    -t DEL \
    --dpi 600 \
    --long_read 1000 \
    -A GFF.bed.gz \
    --transcript_filename sorted_GFF.gff.gz

singularity exec --bind "$(pwd)":"$(pwd)" samplot.sif bash samplot.sh
```

    
