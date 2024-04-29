## ask ChatGPT
https://chat.openai.com/share/44e4cb88-2763-46e8-9a11-739b2fce9455

```bash
docker build -t myrmvp .

docker images

docker tag myrmvp yongjia111/myrmvp:latest

docker push yongjia111/myrmvp:latest

docker run -it myrmvp

singularity build myrmvp.sif docker://yongjia111/myrmvp:latest

singularity run docker://yongjia111/myrmvp:latest Rscript your.R
singularity run myrmvp.sif Rscript your.R
```
