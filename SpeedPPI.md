## https://github.com/patrickbryant1/SpeedPPI
```bash
## input fasta file much have sequence in one line
bash create_ppi_all_vs_all.sh ./data/dev/test.fasta hh-suite/build/bin/hhblits 0.5 ./data/dev/all_vs_all/
bash create_ppi_some_vs_some.sh WT.fasta interacting_protein.fasta hh-suite/build/bin/hhblits 0.5 ./data/dev/some_vs_some/
```
