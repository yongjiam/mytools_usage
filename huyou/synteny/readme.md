## correct the chromosome name for hap1 and hap2
```bash
## sort by chromosome lengths
seqkit sort -lr hap1.fasta > tmp && mv tmp hap1.fasta
seqkit sort -lr hap2.fasta > tmp && mv tmp hap2.fasta

## long to short, named as chr1-9
### sed_id_hap1
s/h1tg000012l/H1_ch1/
s/h1tg000002l/H1_ch2/
s/h1tg000001l/H1_ch3/
s/h1tg000004l/H1_ch4/
s/h1tg000017l/H1_ch5/
s/h1tg000010l/H1_ch6/
s/h1tg000032l/H1_ch7/
s/h1tg000027l/H1_ch8/
s/h1tg000029l/H1_ch9/
### sed_id_hap2
s/h2tg000011l/H2_ch1/
s/h2tg000005l/H2_ch2/
s/h2tg000001l/H2_ch3/
s/h2tg000017l/H2_ch5/
s/h2tg000003l/H2_ch4/
s/h2tg000009l/H2_ch6/
s/h2tg000031l/H2_ch7/
s/h2tg000026l/H2_ch8/
s/h2tg000027l/H2_ch9/
```
