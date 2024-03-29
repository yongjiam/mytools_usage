# guide on running paml natural selection

## install
```bash
conda install -c bioconda paml
```
## prepare input files
### 1. codon-based cds alignment file (I use MEGA, there may be other tools)\
cds.fas
### 2. phylogenetic tree, topology only, without branch length
tree.phy

(((((((((((TraesCS2A02G338600.1,TuG1812G0200003696.01.T01),(TRIDC2AG047520.1,TRITD2Av1G199330.1)),SECCE2Rv1G0105930.1),(TraesCS2B02G343300.1,TRITD2Bv1G162510.1)),(AET2Gv20736300.1,TraesCS2D02G324400.1)),Horvu_13821_2H01G457500.1),((TraesCS1B02G058900.1,TRIDC1BG007330.1),(Horvu_13821_1H01G536500.1,(TuG1812G0100000012.01.T02,(TRITD1Av1G000090.1,TuG1812G0200000917.01.T02))))),(Brame.09UG119000.1,(Brame.09PG123000.1,(Brasy9G166700.1,((Bhyb26.D05G134200.1,Bradi5g12390.1),(Bhyb26.S09G092600.1,Brast09G106000.1)))))),(LOC_Os04g38680.1,((TVU14179.1,TVU14241.1),((Pahal.7G177800.1,Seita.7G124400.1),(Sobic.006G106700.1,Zm00001d003422_T001))))),Ola019279.1),Joasc.15G002900.1);

### 3. specify foreground lineage using treeview software or manually in text editor
tree.phy

(((((((((((TraesCS2A02G338600.1,TuG1812G0200003696.01.T01),(TRIDC2AG047520.1,TRITD2Av1G199330.1)),SECCE2Rv1G0105930.1),(TraesCS2B02G343300.1,TRITD2Bv1G162510.1)),(AET2Gv20736300.1,TraesCS2D02G324400.1)),Horvu_13821_2H01G457500.1) #2,((TraesCS1B02G058900.1,TRIDC1BG007330.1),(Horvu_13821_1H01G536500.1,(TuG1812G0100000012.01.T02,(TRITD1Av1G000090.1,TuG1812G0200000917.01.T02)))) #1),(Brame.09UG119000.1,(Brame.09PG123000.1,(Brasy9G166700.1,((Bhyb26.D05G134200.1,Bradi5g12390.1),(Bhyb26.S09G092600.1,Brast09G106000.1))))) #3),(LOC_Os04g38680.1,((TVU14179.1,TVU14241.1),((Pahal.7G177800.1,Seita.7G124400.1),(Sobic.006G106700.1,Zm00001d003422_T001))))),Ola019279.1),Joasc.15G002900.1);

### 4. create codeml config file, specifying the input files and models
codeml.ctl

### branch specific model:
      seqfile = cds.fas
     treefile = tree.phy
      outfile = BLP_1branch

        noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        model = 0
                    * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = 0   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 2   * initial or fixed omega, for codons or codon-transltd AAs

    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = .0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 4   * # of categories in the dG or AdG models of rates

        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
       RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
       method = 0   * 0: simultaneous; 1: one branch at a time


      * Specifications for duplicating results for the small data set in table 1
      * of Yang (1998 MBE 15:568-573).
      * see the tree file lysozyme.trees for specification of node (branch) labels
### branch-site model
      seqfile = cds.fas
     treefile = tree.phy

      outfile = BLP_ModelA_BLP 
        noisy = 3
      verbose = 0
      runmode = 0

      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0:no clock, 1:global clock; 2:local clock; 3:TipDate
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
         aaRatefile = /../wag.dat
        model = 2  * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = 2  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
        kappa = 3  * initial or fixed kappa
    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1.5  * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 10  * # of categories in dG of NSsites models

        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
       RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

         Small_Diff = .5e-6
      *    cleandata = 1
      *       method = 1   * 0: simultaneous; 1: one branch at a time

### branch-site model null
      seqfile = cds.fas
     treefile = tree.phy

      outfile = BLP_modelAnull_BLP 
        noisy = 3
      verbose = 0
      runmode = 0

      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0:no clock, 1:global clock; 2:local clock; 3:TipDate
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
         aaRatefile = /../wag.dat
        model = 2  * models for codons: 
                    * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = 2  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
        kappa = 3  * initial or fixed kappa
    fix_omega = 1  * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1  * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 10  * # of categories in dG of NSsites models

        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
       RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

         Small_Diff = .5e-6
      *    cleandata = 1
      *       method = 1   * 0: simultaneous; 1: one branch at a time

### ancestral sequence construction
      seqfile = cds.fas * sequence data filename
     treefile = tree.phy      * tree structure file name
      outfile = BLP_ancestral           * main result file name

        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 3  * 0: concise; 1: detailed, 2: too much
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 2  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

            *        ndata = 10
        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
         aaRatefile = c:/users/a1611149/paml4.7/dat/wag.dat  * only used for aa seqs with model=empirical(_F)
                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

        model = 3  * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                         * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
        Mgene = 0
                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
                   * AA: 0:rates, 1:separate

    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2  * initial or fixed kappa
    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 
        omega = .4 * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 8  * # of categories in dG of NSsites models

        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
       RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

         Small_Diff = .5e-6
          cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
            *  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed
             method = 0  * Optimization method 0: simultaneous; 1: one branch a time
      
      * Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
      * 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
      * 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
      * 10: blepharisma nu.
      * These codes correspond to transl_table 1 to 11 of GENEBANK.

(((((((((((TraesCS2A02G338600.1,TuG1812G0200003696.01.T01),(TRIDC2AG047520.1,TRITD2Av1G199330.1)),SECCE2Rv1G0105930.1),(TraesCS2B02G343300.1,TRITD2Bv1G162510.1)),(AET2Gv20736300.1,TraesCS2D02G324400.1)),Horvu_13821_2H01G457500.1),((TraesCS1B02G058900.1,TRIDC1BG007330.1),(Horvu_13821_1H01G536500.1,(TuG1812G0100000012.01.T02,(TRITD1Av1G000090.1,TuG1812G0200000917.01.T02)))) #1),(Brame.09UG119000.1,(Brame.09PG123000.1,(Brasy9G166700.1,((Bhyb26.D05G134200.1,Bradi5g12390.1),(Bhyb26.S09G092600.1,Brast09G106000.1)))))),(LOC_Os04g38680.1,((TVU14179.1,TVU14241.1),((Pahal.7G177800.1,Seita.7G124400.1),(Sobic.006G106700.1,Zm00001d003422_T001))))),Ola019279.1),Joasc.15G002900.1);

### site specific model

## run the program:
### 1. cd to the directory containing the input cds.fas, tree.phy, and codeml.ctl
```bash
cd /target/directory/
```
### 2. run:
```bash
codeml
```
