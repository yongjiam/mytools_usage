## aliyunpan cli use intructions
## installation
https://github.com/aliyun/aliyun-cli
https://github.com/tickstep/aliyunpan/blob/main/docs/manual.md

## configration
```bash
export ALIYUNPAN_CONFIG_DIR=/home/tickstep/tools/aliyunpan/config
aliyunpan config set -savedir <savedir>
```
## get login token
![get refreshtoken](get_token.png)

## login
```bash
aliyunpan login -RefreshToken=626a27b6193f4c5ca6ef0.......
#or
aliyunpan login
请输入RefreshToken, 回车键提交 > 626a27b6193f4c5ca6ef0.......
```
## navigation operations
```bash
aliyunpan update
aliyunpan loglist
aliyunpan who
aliyunpan help
```
