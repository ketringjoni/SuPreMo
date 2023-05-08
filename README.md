Pipeline for scoring variants for disruption to genome folding using Akita (Fudenberg et. al. 2020).
  
Install Akita and its dependencies:  
https://github.com/calico/basenji/tree/master/manuscripts/akita  
  
Packages needed:  
	- Numpy    
	- Pandas  
	- Scipy   
	- Tensorflow  
	- Basenji and its dependencies  
	- Cooltools  
	- Pysam   
	- Math  
	- Collection  
	- Bioseq  
  
usage: Akita_variant_scoring [-h] --in IN_FILE [--fa FASTA]
[--chrom CHROM_LENGTHS]
                             [--centro CENTROMERE_COORDS]
                             [--score {mse,corr} [{mse,corr} ...]]
                             [--shift_by SHIFT_WINDOW [SHIFT_WINDOW ...]]
                             [--file OUT_FILE] [--dir OUT_DIR]
                             [--limit SVLEN_LIMIT]

Pipeline for scoring variants for disruption to genome folding using Akita
(Fudenberg et. al. 2020).

optional arguments:
  -h, --help            show this help message and exit
  --in IN_FILE          Input file with variants.
  --fa FASTA            hg38 reference genome fasta file.
  --chrom CHROM_LENGTHS
                        File with lengths of chromosomes in hg38. Columns:
                        chromosome (ex: 1), length; no header.
  --centro CENTROMERE_COORDS
                        Centromere coordinates for hg38. Columns: chromosome
                        (ex: chr1), start, end; no header.
  --score {mse,corr} [{mse,corr} ...]
                        Method(s) used to calculate disruption scores.
  --shift_by SHIFT_WINDOW [SHIFT_WINDOW ...]
                        Values for shifting prediciton windows (e.g. to
                        predict with no shift (variant centered) and shift by
                        1 bp to either side, use: -1 0 1). Values outside of
                        range -450000 ≤ x ≤ 450000 will be ignored. Prediction
                        windows at the edge of chromosome arms will only be
                        shifted in the direction that is possible (ex. for
                        window at chrom start, a -1 shift will be treated as a
                        1 shift since it is not possible to shift left.)
  --file OUT_FILE       Prefix for output files.
  --dir OUT_DIR         Output directory.
  --limit SVLEN_LIMIT   Maximum structural variant length to be scored (<=
                        700000).

  
Test run on a vcf of structural variants:  
> python score_var.py --in test/input/subset.somaticSV.vcf --shift_by -1 0 1 --file subset_SV --dir test/output
