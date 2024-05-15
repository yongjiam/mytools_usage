## fast genome assembly alignment
### run
```
## approximate align -m
srun --export=all -n 1 -c 128 singularity exec $IMAGE /usr/software/wfmash/build/bin/wfmash \
       --threads 128 -n 2 -m /scratch/pawsey0399/yjia/shunlin/morexV3/genome.fasta \
       /scratch/pawsey0399/yjia/WBT/hifionly_run2/out_JBAT2_seded_sorted.FINAL.fa  > morexV3_align_out_JBAT2_seded_sorted.paf

## base pair align
srun --export=all -n 1 -c 128 singularity exec $IMAGE /usr/software/wfmash/build/bin/wfmash \
       --threads 128 -n 2 /scratch/pawsey0399/yjia/shunlin/morexV3/genome.fasta \
       /scratch/pawsey0399/yjia/WBT/hifionly_run2/out_JBAT2_seded_sorted.FINAL.fa  > morexV3_align_out_JBAT2_seded_sorted.paf

## output sam -a
srun --export=all -n 1 -c 128 singularity exec $IMAGE /usr/software/wfmash/build/bin/wfmash \
       --threads 128 -n 2 -a /scratch/pawsey0399/yjia/shunlin/morexV3/genome.fasta \
       /scratch/pawsey0399/yjia/WBT/hifionly_run2/out_JBAT2_seded_sorted.FINAL.fa  > morexV3_align_out_JBAT2_seded_sorted.sam
```
### usage
Singularity> /usr/software/wfmash/build/bin/wfmash
  /usr/software/wfmash/build/bin/wfmash [target] [queries...] {OPTIONS}

    wfmash: a pangenome-scale aligner, v0.11.0-13-g67ab187

  OPTIONS

      [ MANDATORY OPTIONS ]
        target                            alignment target/reference sequence file
      [ Files IO Options ]
        queries...                        query sequences file
        -Q[queries], --query-file-list=[queries]
                                          alignment queries files list
      [ Mapping Options ]
        **-p[%], --map-pct-id=[%]           percent identity in the mashmap step [default: 90]**
        **-s[N], --segment-length=[N]       segment seed length for mapping [default: 5k]**
        -l[N], --block-length=[N]         keep merged mappings supported by homologies of this total
                                          length [default: 5*segment-length]
        **-n[N], --num-mappings-for-segment=[N]**
                                          number of mappings to retain for each segment [default: 1]
        -S[N], --num-mappings-for-short-seq=[N]
                                          number of mappings to retain for each sequence shorter
                                          than segment length [default: 1]
        **-k[N], --kmer=[N]                 kmer size [default: 19]**
        -H[%], --kmer-threshold=[%]       ignore the top % most-frequent kmers [default: 0.001]
        -L, --lower-triangular            only map shorter sequences against longer
        -X, --skip-self                   skip self mappings when the query and target name is the
                                          same (for all-vs-all mode)
        -Y[C], --skip-prefix=[C]          skip mappings when the query and target have the same
                                          prefix before the last occurrence of the given character C
        -P[pfx], --target-prefix=[pfx]    use only targets whose name starts with this prefix
        -A[FILE], --target-list=[FILE]    file containing list of target sequence names to use
        -m, --approx-map                  skip base-level alignment, producing an approximate
                                          mapping in PAF
        -N, --no-split                    disable splitting of input sequences during mapping
                                          [default: enabled]
        -c[N], --chain-gap=[N]            chain mappings closer than this distance in query and
                                          target, sets approximate maximum variant length detectable
                                          in alignment [default: 20k]
        -K, --drop-low-map-id             drop mappings with estimated identity below --map-pct-id=%
        -f, --no-filter                   disable mapping filtering
        -x[FACTOR], --sparsify-mappings=[FACTOR]
                                          keep this fraction of mappings
        -w[N], --sketch-size=[N]          sketch size for sketching.
        -J[F], --kmer-complexity=[F]      Drop segments w/ predicted kmer complexity below this
                                          cutoff. Kmer complexity defined as #kmers / (s - k + 1)
        -1, --no-hg-filter                Don't use the hypergeometric filtering and instead use the
                                          MashMap2 first pass filtering.
        -2[%], --hg-filter-ani-diff=[%]   Filter out mappings unlikely to be this ANI less than the
                                          best mapping [default: 0.0]
        -3[%], --hg-filter-conf=[%]       Confidence value for the hypergeometric filtering
                                          [default: 99.9%]
        -M, --no-merge                    don't merge consecutive segment-level mappings
      [ Alignment Options ]
        -i[FILE], --input-paf=[FILE]      derive precise alignments for this input PAF
        -O, --invert-filtering            if an input PAF is specified, remove alignments with
                                          gap-compressed identity below --map-pct-id x 0.8, else
                                          keep all alignments [default: if an input PAF is
                                          specified, keep all alignments, else remove alignments
                                          with gap-compressed identity below --map-pct-id x 0.8]
        -W[N], --wflamda-segment=[N]      wflambda segment length: size (in bp) of segment mapped in
                                          hierarchical WFA problem [default: 256]
        -g[mismatch,gap1,ext1],
        --wfa-params=[mismatch,gap1,ext1] score parameters for the wfa alignment (affine); match
                                          score is fixed at 0 [default: 4,6,1]
        -G[mismatch,gap1,ext1],
        --wflign-params=[mismatch,gap1,ext1]
                                          score parameters for the wflign alignment (affine); match
                                          score is fixed at 0 [default: 4,6,1]
        -b[N], --max-mash-dist=[N]        maximum mash distance to perform the alignment in a
                                          wflambda segment [default: adaptive with respect to the
                                          estimated identity]
        -j[N], --wflign-min-wf-len=[N]    min wavefront length for heuristic WFlign [default: 1024]
        -q[N], --wflign-max-distance=[N]  max distance threshold for heuristic WFlign [default:
                                          2048/(estimated_identity^2)]
        -C[N], --max-patch-major=[N]      maximum length to patch in the major axis [default:
                                          512*segment-length]
        -F[N], --max-patch-minor=[N]      maximum length to patch in the minor axis [default:
                                          128*segment-length]
        -E[N], --erode-match-mismatch=[N] maximum length of match/mismatch islands to erode before
                                          patching [default: adaptive]
        --max-patching-score=[N]          maximum score allowed when patching [default: adaptive
                                          with respect to gap penalties and sequence length]
      [ Output Format Options ]
        -d, --md-tag                      output the MD tag
        **-a, --sam-format                  output in the SAM format (PAF by default)**
        -q, --no-seq-in-sam               do not fill the sequence field in the SAM format
      [ General Options ]
        -B[PATH], --tmp-base=[PATH]       base name for temporary files [default: `pwd`]
        -Z, --keep-temp                   keep intermediate files
      [ Threading ]
        **-t[N], --threads=[N]              use this many threads during parallel steps**
      [ Program Information ]
        -v, --version                     show version number and github commit hash
        -h, --help                        display this help menu
