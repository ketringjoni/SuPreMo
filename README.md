Pipeline for scoring variants for disruption to genome folding using Akita (Fudenberg et. al. 2020)  
  
Download and install Akita and its dependencies. Instructions here:  
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
  
usage: Akita_variant_scoring [-h] --in IN_FILE [--format {vcf,df}]  
                             [--type {simple,SV}] [--fa FASTA]  
                             [--chrom CHROM_LENGTHS]  
                             [--centro CENTROMERE_COORDS]  
                             [--score {mse,corr,both}]  
                             [--shift_by SHIFT_WINDOW [SHIFT_WINDOW ...]]  
                             [--out OUT_FILE] [--limit SVLEN_LIMIT]  
  
Score variants based on disruption to genome folding in the surround 1 Mb region using Akita (Fudenber et. al. 2020).

optional arguments:  
  -h, --help            show this help message and exit  
  --in IN_FILE          Input file with variants.  
  --format {vcf,df}     Format for input file.  
  --type {simple,SV}    Variant type: simple or SV.  
  --fa FASTA            hg38 reference genome fasta file.  
  --chrom CHROM_LENGTHS  
                        File with lengths of chromosomes in hg38. Columns:  
                        chromosome (ex: 1), length; no header.  
  --centro CENTROMERE_COORDS  
                        Centromere coordinates for hg38. Columns: chromosome  
                        (ex: chr1), start, end; no header.  
  --score {mse,corr,both}  
                        Method(s) used to calculate disruption scores.  
  --shift_by SHIFT_WINDOW [SHIFT_WINDOW ...]  
                        Values for shifting prediciton windows (e.g. to  
                        predict with no shift (variant centered) and shift by  
                        1 bp to either side, use: -1 0 1). Values outside of  
                        range -450000 ≤ x ≤ 450000 will be ignored.  
  --out OUT_FILE        Directory in which to save the output file(s).  
  --limit SVLEN_LIMIT   Maximum structural variant length to be scored.  
  
  
Test run on SVs and simple variants:  
python score_var.py --in test/input_file_types/subset.somaticSV.vcf --shift_by -1 0 1 --out subset_SV_vcf  
python score_var.py --in test/input/subset.consensus_somatic.public.vcf --type simple --shift_by -1 0 1 --out subset_simple_vcf  
  
The resulting outputs (saved as 'score_var_output' in cwd) are in /test/score_var_output.  

