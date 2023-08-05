# SuPreMo: Sequence Mutator for Predictive Models

**S**equence M**u**tator for **Pre**dictive **Mo**dels (SuPreMo) is a pipeline for generating reference and mutated sequences from variants for input into predictive models and includes an application using the [Akita](https://www.nature.com/articles/s41592-020-0958-x) model to measure disruption to genome folding.
  
  
  
## Installation

To only use get_seq, create a conda environment with the packages outlined below under Overall. Compatible package versions shown in [get_seq_env.yml](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/get_seq_env.yml).

To use get_scores, you'll need to install Akita, Basenji and their dependencies. Create a conda environment following instructions [here](https://github.com/calico/basenji/tree/master/manuscripts/akita). Compatible package versions shown in [get_scores_env.yml](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/get_scores_env.yml).
 
**Requirements:**
- Overall
    * Python
    * Pysam
    * Pandas
    * Numpy
    * Pathlib
    * Biopython
- Only for get_scores (Akita, Basenji and their dependencies)
    * Pybedtools
    * PyBigWig
    * Cython
    * Astropy
    * Intervaltree
    * Natsort
    * Scipy
    * Cooltools
    * Basenji
    * Tensorflow
    * Jupyter
    * Protobuf
    * Patsy
    * Libclang
- Only for plotting maps
    * Matplotlib
    
    
Download tool by: 
```shell
git clone https://github.com/ketringjoni/Akita_variant_scoring.git
```



## Use

Both get_seq and get_scores are implemented under the same python script score_var.py, used as follows:

```shell
python scripts/score_var.py INPUT <args>
```
**Required argument:**
    * `INPUT`: Input file with variants. Accepted formats are: VCF, TSV (SVannot output format), BED, TXT. Can be gzipped. Coordinates should be 1-based and left open, except for SNPs (following vcf [4.1/4.2 specifications](https://samtools.github.io/hts-specs/VCFv4.1.pdf).
    
**Optional arguments:**
- Overall arguments
    * `--help`: Print help page and exit.
    * `--fa FASTA_FILE`: path to hg38.fa file.
    * `--file FILE_PREFIX`: Prefix for output files.
    * `--dir DIRECTORY`: Path for output files.
    * `--limit LENGTH`: Limit for SV length.
    * `--nrows N_ROWS`: Number of rows to read at a time from input to save memory.
- get_seq arguments
    * `--get_seq`: Get reference and mutated sequences.
    * `--shift_by SHIFT`: Values to shift sequence window by.
    * `--seq_len LENGTH`: Sequence length.
    * `--revcomp`: Take the reverse complement of the sequence. Options: no_revcomp,add_revcomp,only_revcomp
- get_scores arguments
    * `--get_scores`: Get disruption scores using Akita.
    * `--scores SCORES`: Scores to be used from [Gunsalus and McArthur et al](https://www.biorxiv.org/content/10.1101/2023.04.04.535480v1.full.pdf). Options: mse,corr,ssi,scc,ins,di,dec,tri,pca.
    * `--augment`: Get scores for augmented predictions (Average with predictions from 1bp shifts and reverse complement). 
    * `--get_tracks`: Get disruption score tracks.
    * `--get_maps`: Get predicted contact frequency maps.
    
For more details on how to use arguments, refer to help page printed at the top of [score_var.py](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/scripts/score_var.py).


### Test

Run the following command to test the tool in the tool directory. Expected output is in [test/output/](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/test/output/)

```shell
python scripts/score_var.py test/input/run_test_SV.bed --shift_by -10000 0 10000 --dir test/ouput --file run_test_SV --revcomp add_revcomp --get_scores
```



## Tutorials

- Creating input files with custom perturbations
    * [custom_perturbations.ipynb](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/walkthroughs/custom_perturbations.ipynb)
- Using the tool and interpreting the results
    * [get_seq_walkthrough.ipynb](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/walkthroughs/get_seq_walkthrough.ipynb)
    * [get_scores_walkthrough.ipynb](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/walkthroughs/get_scores_walkthrough.ipynb)



**Note:**
- This tool currently only works with hg38
- Coordinates can be 0- or 1- based:
    * Input files: 0-based, left inclusive
    * Fasta file: 1-based (0-based once read into python)
    * [Chromosome lengths](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/data/chrom_lengths_hg38): 1-based, inclusive
    * [Centromere coordinates](https://github.com/ketringjoni/Akita_variant_scoring/blob/main/data/centromere_coords_hg38): 0-based, left inclusive
    
    

For any questions and/or feedback, please reach out to katie.gjoni at ucsf dot edu



