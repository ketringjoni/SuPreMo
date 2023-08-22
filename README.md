# SuPreMo: Sequence Mutator for Predictive Models

<ins>S</ins>equence M<ins>u</ins>tator for <ins>Pre</ins>dictive <ins>Mo</ins>dels (SuPreMo) is a pipeline for generating reference and mutated sequences from perturbations such as structural variants for input into predictive models and includes an application using the [Akita](https://www.nature.com/articles/s41592-020-0958-x) model to measure disruption to genome folding.

SuPreMo includes 2 modules. For both modules, the input is a set of perturbations. The possible parameters and outputs differ.

    1. **get_seq: Get sequences for predictive models**
    * This module generates reference and alternate sequences for each perturbation under each provided augmentation parameter. 
    * The sequences are accompanied by the relative position of the perturbation for each sequence. 
    * The perturbation is by default centered in the sequence if possible unless the user provides a shifting parameter. 
    * These sequences can be inputed into predictive models and the results can be compared between the reference and alternate predictions to evaluate the effect of the perturbation on the predicted outcome.

    2. **get_scores: Get disruption scores using the Akita model**
    * This module generated disruption scores by comparing predicted contact frequency maps from the reference and alternate sequences. 
    * Optional outputs are predicted maps or tracks of disruption scores along the map. 
    * The maps are accompanied by the relative position of the perturbation and the start coordinate that the predicted map corresponds to. 
    * The perturbation is by default centered in the map if possible unless the user provides a shifting parameter.
  
  
  
## Installation

To only use get_seq, create a conda environment with the packages outlined below under Overall or using [get_seq_env.yml](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/get_seq_env.yml).

To use get_scores, you'll need to install Akita, Basenji and their dependencies. Create a conda environment following the recommandation [here](https://github.com/calico/basenji/tree/master/manuscripts/akita). Compatible package versions shown in [get_scores_env.yml](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/get_scores_env.yml).
 
**Requirements:**
- Overall
    * Python
    * Pysam
    * Pandas
    * Numpy
    * Pathlib
    * Biopython
- Only for get_scores
    * Akita, Basenji, and their dependencies
- Only for walkthroughs
    * Jupyter
    * Matplotlib
    
    
Download SuPreMo by: 
```shell
git clone https://github.com/ketringjoni/Akita_variant_scoring.git
```



## Use

Both get_seq and get_scores are implemented under the same python script SuPreMo.py, used as follows:

```shell
python scripts/SuPreMo.py INPUT <args>
```

**Required arguments:**
- `INPUT`: Input file with perturbations. Accepted formats are: VCF, TSV (SVannot output format), BED, TXT. Can be gzipped. Coordinates should be 1-based and left open, except for SNPs (following vcf [4.1/4.2 specifications](https://samtools.github.io/hts-specs/VCFv4.1.pdf)).
- At least one of the following:
    * `--get_seq`: Get reference and mutated sequences. 
    * `--get_scores`: Get disruption scores using Akita.
    
**Optional arguments:**
- Overall arguments
    * `--help`: Print help page and exit.
    * `--fa FASTA_FILE`: path to hg38.fa file.
    * `--file FILE_PREFIX`: Prefix for output files.
    * `--dir DIRECTORY`: Path for output files.
    * `--limit LENGTH`: Limit for SV length.
    * `--nrows N_ROWS`: Number of rows to read at a time from input to save memory.
- get_seq arguments
    * `--shift_by SHIFT`: Values to shift sequence window by.
    * `--seq_len LENGTH`: Sequence length.
    * `--revcomp`: Take the reverse complement of the sequence. Options: no_revcomp,add_revcomp,only_revcomp
- get_scores arguments
    * `--scores SCORES`: Scores to be used to calculate 3D genome folding disruption. Options: mse, corr, ssi, scc, ins, di, dec, tri, and pca, from [Gunsalus and McArthur et al](https://www.biorxiv.org/content/10.1101/2023.04.04.535480v1.full.pdf).
    * `--augment`: Get scores for augmented predictions (Average with predictions from 1bp shifts and reverse complement). 
    * `--get_tracks`: Get disruption score tracks.
    * `--get_maps`: Get predicted contact frequency maps.
    
If multiple inputs are given for an argument, they should be space-separated.

For more details on how to use arguments, refer to help page printed at the top of [SuPreMo.py](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/scripts/SuPreMo.py).



## Tutorials

- Example application on structural variants from WGS in HCC1395 tumor cells.
    * [example_application.ipynb](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/walkthroughs/example_application.ipynb)
- Creating input files with custom perturbations
    * [custom_perturbations.ipynb](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/walkthroughs/custom_perturbations.ipynb)
- Using SuPreMo and interpreting the results
    * [get_seq_walkthrough.ipynb](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/walkthroughs/get_seq_walkthrough.ipynb)
    * [get_scores_walkthrough.ipynb](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/walkthroughs/get_scores_walkthrough.ipynb)



**Note:**
- SuPreMo currently only works with hg38
- Coordinates can be 0- or 1- based:
    * Input files: 0-based, left inclusive
    * Fasta file: 1-based (0-based once read into python)
    * [Chromosome lengths](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/data/chrom_lengths_hg38): 1-based, inclusive
    * [Centromere coordinates](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/data/centromere_coords_hg38): 0-based, left inclusive
    


## Test sets

We have generated two categories of test sets:

1. [Test sets for edge cases](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/test/test_set_edge_cases/)
These are sets of variants that are meant to include all edge cases that SuPreMo should handle and to ensure that SuPreMo is working appropriately and that only expected errors and warnings appear. These were run using the following commands and the outputs are saved in [test/test_set_edge_cases/](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/test/test_set_edge_cases/)

```shell
python scripts/SuPreMo.py test_data/test_set_edge_cases/test_set_edge_SV.bed \
                            --dir test_data/test_set_edge_cases \
                            --file test_set_edge_SV \
                            --shift_by -10000 0 10000 \
                            --get_scores
                            
python scripts/SuPreMo.py test_data/test_set_edge_cases/test_set_edge_simple.bed \
                            --dir test_data/test_set_edge_cases \
                            --file test_set_edge_simple \
                            --shift_by -10000 0 10000 \
                            --revcomp add_revcomp \
                            --get_scores
```

2. [Test set for sequences](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/test/test_set_sequences/)
This is a set of 10 variants for which sequences were generated manually to ensure that SuPreMo generates reference and mutated sequences accurately. The output sequences for this seq from the variant and the manual curation are not included in this repo due to size but they were tested to be exactly the same.

  

***
For any questions and/or feedback, please reach out to katie.gjoni at ucsf dot edu



