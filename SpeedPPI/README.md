## SpeedPPI computer-based yeast-2-hybridization
https://github.com/patrickbryant1/SpeedPPI

## setup

## run
```
bash create_ppi_some_vs_some.sh  chunshen_query1.fasta chunshen_query2.fasta  ./hh-suite/bin/hhblits  0.5 ./chunshen_output2/
```
### Notes:
> nvidia-smi  ## check GPU memory
> free -h ## check CPU memory
> 2vs3 proteins use more memory than 1vs3
