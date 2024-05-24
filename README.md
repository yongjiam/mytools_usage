## mytools_usage
docker
singularity
github commands

## increase swap space in linux
https://askubuntu.com/questions/178712/how-to-increase-swap-space
```
dd if=/dev/zero of=/data/swapfile.img bs=1024 count=10M
mkswap /data/swapfile.img
sudo swapon /data/swapfile.img
```
