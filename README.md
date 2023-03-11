# Akita_variant_scoring
Pipeline for scoring simple and structural variants using Akita.

Download and install Akita and its dependencies. Instructions here:
https://github.com/calico/basenji/tree/master/manuscripts/akita

You will need the following in a directory within cwd called “Akita_model”:
- model_best.h5  
- params.json

Input: Data frame with variant information. 
Required columns (data retrieved from vcf file):
- CHROM (string): Chromosome in form chr1.
- POS (integer): Variant position.
- REF (string): Reference allele, not used for SVs.
- ALT (string): Alternate allele for simple variants or coordinates for chromosomal rearrangements (SVTYPE = ‘BND’).
- END (string): SV only. End coordinate for SV (not used for BND)
- SVTYPE (string): SV only. Options: BND, DEL, DUP, INV, INS. Will not score INS without inserted sequence. Note: assumes duplications are tandem.

To run test:
python score_var.py test/simple_var_test /pollard/data/vertebrate_genomes/human/hg38/hg38/hg38.fa data/chrom_lengths data/centromere_coords both test/simple_var_out
python score_var.py test/SV_test /pollard/data/vertebrate_genomes/human/hg38/hg38/hg38.fa data/chrom_lengths data/centromere_coords both test/SV_out
