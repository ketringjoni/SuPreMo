<img src="https://github.com/ketringjoni/SuPreMo/blob/main/SuPreMo_logo.png" alt="image" width="30%" height="auto">

<ins>S</ins>equence M<ins>u</ins>tator for <ins>Pre</ins>dictive <ins>Mo</ins>dels (SuPreMo) is a pipeline for generating reference and perturbed sequences for input into predictive models that is scalable, flexible, and comprehensive. SuPreMo-Akita applies the tool to an existing sequence-to-profile model, [Akita](https://www.nature.com/articles/s41592-020-0958-x), and generates scores that measure disruption to genome folding. 

**SuPreMo: Get sequences for predictive models**
- SuPreMo incorporates variants one at a time into the human reference genome and generates reference and alternate sequences for each perturbation under each provided augmentation parameter. 
- The sequences are accompanied by the relative position of the perturbation for each sequence. 
- The perturbation is by default centered in the sequence, if possible, unless the user provides a shifting parameter. 

**SuPreMo-Akita: Get 3D genome folding disruption scores using the Akita model**
- SuPreMo-Akita generates disruption scores by comparing contact frequency maps predicted from the reference and alternate sequences.
- The maps are accompanied by the relative position of the perturbation and the start coordinate that the predicted maps correspond to.
- The perturbation is by default centered in the map, if possible, unless the user provides a shifting parameter.
- Regions of interest can be inputted, such as accessible regions in the applicable cell type, to weight the disruption scores by.
- Optionally, the predicted maps and/or disruption score tracks along those maps can also be outputted. 
  
  
  
## Installation
  
### **Download repo**
```shell
git clone https://github.com/ketringjoni/SuPreMo.git
```


### **Install SuPreMo or SuPreMo-Akita**

**For SuPreMo:**

1. Create and activate a conda environment using the given yml file (this should typically be done in the conda "base" environment):
```shell
cd SuPreMo/
conda env create -f supremo_env.yml
```
*Make sure this runs with no errors.*

```shell
conda activate supremo_env
```

2. Test that all SuPreMo requirements are properly installed:
```shell
python scripts/test_install_SuPreMo.py
```
*This should print: SuPreMo packges successfully imported.*

You may also need to install the following: gcc, g++ and libz-dev.


**For SuPreMo-Akita:**

1. Create conda environment with tensorflow and check for its proper installation (this should typically be done in the conda "base" environment):
```shell
conda create -n supremo_akita_env python=3.10 numpy scipy pandas jupyter tensorflow
```
*You might be asked to confirm installation with a 'y'.*

```shell
conda activate supremo_akita_env
python -c "import tensorflow"
```
*This should not result in any errors (warnings are ok).*


2. Install [basenji](https://github.com/calico/basenji) with no dependencies and set environmental variables:
```shell
cd ~/
git clone https://github.com/calico/basenji.git
```
```shell
cd basenji/
python setup.py develop --no-deps
```
```shell
export BASENJIDIR=~/basenji
export PATH=$BASENJIDIR/bin:$PATH
export PYTHONPATH=$BASENJIDIR/bin:$PYTHONPATH
```

If downloading the Basenji repo in a directory other than home (~), adjust the path in the export command accordingly.


3. Install dependencies:
```shell
pip install astropy tensorflow-io-gcs-filesystem patsy libclang cooltools Cython biopython pathlib pysam natsort intervaltree pybedtools pybigwig qnorm seaborn statsmodels tabulate jax wrapt==1.14
```
*If pip installation does not work in your system, please troubleshoot the installation of all of the above pacakges.*

4. Test that all SuPreMo-Akita requirements are properly installed:
```shell
cd path_to_SuPreMo/
python scripts/test_install_SuPreMo-Akita.py
```
*This should print: SuPreMo packges successfully imported. SuPreMo-Akita packges successfully imported.*


For running walkthroughs/tutorials, Jupyter and Matplotlib are also required.


### **Download genome fasta file (optional\*)**
```shell
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/<genome>.fa.gz -P path_to_SuPreMo/data/
gunzip path_to_SuPreMo/data/<genome>.fa.gz 
```
Currently accepted genomes are hg19 and hg38. 

\*If you choose to opt out of this, you will need to specify the path to an existing fasta file using the --fa parameter.



## Use

### Input format

**Suggested input types**
- VCF file
    * following [vcf 4.1/4.2 specifications](https://samtools.github.io/hts-specs/VCFv4.1.pdf)
- TXT file
    * Columns required for simple variants: CHROM, POS, REF, ALT
    * Columns required for structural variants: CHROM, POS, REF, ALT, END, SVTYPE (SV type), SVLEN (SV length)
    
**Additional input types**
- BED-like file
    * Columns required for simple variants (exclude column names): chrom, pos, end, ref, alt
    * Columns required for structural variants(exclude column names): chrom, pos, end, ref, alt, SV type, SV length
- TSV file
    * following [AnnotSV output format](https://lbgi.fr/AnnotSV/Documentation/README.AnnotSV_latest.pdf)
    
**Optional additional inputs only for SuPreMo-Akita**
- Sequences outputted from SuPreMo in fasta format.
- Regions of interest to use for weighing. Text file with columns: *chrom, start, end* or *chrom, start1, end1, start2, end2* for paired regions.
- Amounts to shift prediction window by. Text file with column *shift_by* that contains the same number of rows as variants in the input.



### Running SuPreMo

Both SuPreMo and SuPreMo-Akita are implemented under the same python script SuPreMo.py, used as follows:

```shell
python scripts/SuPreMo.py INPUT <args>
```

**Required arguments**
- `INPUT`: Input file with perturbations. File name must end with one of the following suffixes: vcf, txt, bed, tsv. Can be gzipped. Coordinates should be ([1-based and left-open](https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/)), except for SNPs which are fully closed.
- At least one of the following:
    * `--get_seq`: Get reference and mutated sequences. 
    * `--get_Akita_scores`: Get disruption scores using Akita.
    
**Optional arguments**
- Overall arguments
    * `--help`: Print help page and exit.
    * `--fa FASTA_FILE`: Path to human fasta file. Default is data/{genome}.fa, where genome is hg38 or hg19.
    * `--file FILE_PREFIX`: Prefix for output files. Default: SuPreMo.
    * `--dir OUTPUT_DIRECTORY`: Path to directory for output files (created if not present). Default: SuPreMo_output.
    * `--nrows N_ROWS`: Number of rows to read at a time from input to save memory. Default: 1000.
    * `--genome GENOME`: Genome to be used: hg19 or hg38 (default).
- get_seq arguments
    * `--seq_len SEQ_LENGTH`: Sequence length. Default: 1048576.
    * `--limit VAR_LENGTH`: Limit for variant length. Default: 2/3 of seq_len.
    * `--shift_by SHIFT`: Values to shift sequence window by. Default: 0.
    * `--shifts_file FILE_PATH`: Path to file with column shift_by containing values to shift sequences by for each input variant. Default: None.
    * `--revcomp OPTION`: Take the reverse complement of the sequence. Options: no_revcomp (default),add_revcomp,only_revcomp.
- get_Akita_scores arguments
    * `--scores METHOD`: Method(s) to be used to calculate 3D genome folding disruption scores. Options: mse, corr, ssi, scc, ins, di, dec, tri, and pca, from [Gunsalus and McArthur et al](https://www.biorxiv.org/content/10.1101/2023.04.04.535480v1.full.pdf). Default: mse corr.
    * `--roi REGIONS`: Region of interest (roi) to weight. Provide path to BED-like file or 'genes' to weight gene transcription start sites. Default: None.
    * `--roi_scales SCALE`: Scale by which to weight roi. Default: None.
    * `--Akita_cell_types CELL_TYPE`: Cell type to make predictions for, from ones Akita was trained on: HFF (default), H1hESC, GM12878, IMR90, HCT116.
    * `--augment`: Get scores for augmented predictions (mean and median scores from predictions with shifts and reverse complement).
    * `--get_maps`: Get predicted contact frequency maps.
    * `--get_tracks`: Get disruption score tracks.

If multiple inputs are given for an argument, they should be space-separated.


For more details on how to use arguments, refer to help page printed at the top of [SuPreMo.py](https://github.com/ketringjoni/SuPreMo/blob/main/scripts/SuPreMo.py).


### Output format

**SuPreMo output**
- Sequences
    * Fasta file 
    * Sequence name format: {variant index<sup>1</sup>}_{shift<sup>2</sup>}_{reverse complement<sup>3</sup>}_{sequence index<sup>4</sup>}_{variant relative position<sup>5</sup>}
    * There are 2-3 entries per prediction (2 for non-BND variants and 3 for BND variants).

**SuPreMo-Akita output**
- 3D genome disruption scores
    * Tab-delimited text file
    * The first column name: var_index<sup>1</sup>
    * The rest of the column names: {method<sup>6</sup>}_{scale}-weighted<sup>7</sup>_{cell_type<sup>8</sup>}_{shift<sup>2</sup>}_{reverse complement<sup>3</sup>}
    * For variants with multiple alternate alleles, the score columns will contain a score for each allele, separated by a comma.    
- Contact frequency maps
    * Numpy file
    * Dictionary item name format: {variant index<sup>1</sup>}_maps_{cell_type<sup>8</sup>}_{shift<sup>2</sup>}_{reverse complement<sup>3</sup>}
    * There is 1 entry per prediction. Each entry contains the following: 2 (3 for chromosomal rearrangements) arrays that correspond to the upper right triangle of the predicted contact frequency maps, the relative variant position in the map, and the first coordinate of the sequence that the map corresponds to. 
- Disruption score tracks
    * Numpy file
    * Dictionary item name format: {variant index<sup>1</sup>}_{method<sup>6</sup>}_track_{cell_type<sup>8</sup>}_{shift<sup>2</sup>}_{reverse complement<sup>3</sup>}
    * There is 1 entry per prediction: a 448x1 array.
- IDs for weighted regions of interest (roi) in prediction window
    * Text file with columns: var_index	and roi_ids_{shift<sup>2</sup>}. It provides 0-based row indexes from roi input that overlap each prediction window.
    
    
**Superscript descriptions:**
1. Input row number. (For sequences, this is followed by _N for each allele of variants with multiple alternate alleles).
2. Integer that window is shifted by. If shifts_file is provided, it will just say 'shifted' instead;.
3. 'revcomp' present only if reverse complement of sequence was taken.
4. Index for sequences generated for that variant (0-1 for non-BND reference and alternate sequences and 0-2 for BND left and right reference sequence and alternate sequence). 
5. Relative position of variant in sequence (list of two for non-BND variant positions in reference and alternate sequence and an integer for BND breakend position in reference and alternate sequences).
6. Method used to calculate disruption score.
7. Scale by which regions of interests were weighted. Only present if --roi was specified.
8. Cell type used for Akita prediction.



## Tutorials

- Example application on structural variants from WGS in HCC1395 tumor cells.
    * [example_application.ipynb](https://github.com/ketringjoni/SuPreMo/blob/main/walkthroughs/example_application.ipynb)
- Creating input files with custom perturbations.
    * [custom_perturbations.ipynb](https://github.com/ketringjoni/SuPreMo/blob/main/walkthroughs/custom_perturbations.ipynb)
- Walkthrough for running SuPreMo with DeepSEA.
    * [SuPreMo_walkthrough.ipynb](https://github.com/ketringjoni/SuPreMo/blob/main/walkthroughs/SuPreMo_walkthrough.ipynb)
- Walkthrough for running SuPreMo-Akita, reading outputs, and interpreting results.
    * [SuPreMo-Akita_walkthrough.ipynb](https://github.com/ketringjoni/SuPreMo/blob/main/walkthroughs/SuPreMo-Akita_walkthrough.ipynb)
- Walkthrough for running SuPreMo-Akita with regions of interest (roi) to get weighted scores.
    * [SuPreMo-Akita_weighted_scores.ipynb](https://github.com/ketringjoni/SuPreMo/blob/main/walkthroughs/SuPreMo-Akita_weighted_scores.ipynb)



**To note**
- SuPreMo currently only works with hg38 and hg19
- Coordinates can be 0- or 1- based:
    * Input files: 1-based, left-open (same as 0-based, right-open)
    * Fasta file: 1-based (0-based once read into python)
    * [Chromosome lengths](https://github.com/ketringjoni/SuPreMo/blob/main/data/chrom_lengths_hg38): 1-based, fully closed
    * [Centromere coordinates](https://github.com/ketringjoni/SuPreMo/blob/main/data/centromere_coords_hg38): 0-based, right-open
    


## Test sets

We have generated two categories of test sets:

1. [Test sets for edge cases](https://github.com/ketringjoni/SuPreMo/blob/main/test_data/test_set_edge_cases/)

These are sets of 347 variants that are meant to include all edge cases that SuPreMo should handle and to ensure that SuPreMo is working appropriately with only expected errors and warnings appearing. These were run using the following commands and the outputs are saved in [test_data/test_set_edge_cases/output](https://github.com/ketringjoni/SuPreMo/blob/main/test_data/test_set_edge_cases/output/). If using regions of interest (roi) to weight, the simple variant test set was also run up-weighting genes and using random shifts for each variant, with outputs saved in [test_data/test_set_edge_cases/weighted_output](https://github.com/ketringjoni/SuPreMo/blob/main/test_data/test_set_edge_cases/weighted_output/).

```shell
python scripts/SuPreMo.py test_data/test_set_edge_cases/input/test_set_edge_SV.bed \
                            --dir test_data/test_set_edge_cases/output \
                            --file test_set_edge_SV \
                            --shift_by -10000 0 10000 \
                            --get_Akita_scores
```

```shell                            
python scripts/SuPreMo.py test_data/test_set_edge_cases/input/test_set_edge_simple.bed \
                            --dir test_data/test_set_edge_cases/output \
                            --file test_set_edge_simple \
                            --shift_by -10000 0 10000 \
                            --revcomp add_revcomp \
                            --get_Akita_scores
```

```shell
python scripts/SuPreMo.py test_data/test_set_edge_cases/input/test_set_edge_simple.bed \
                            --dir test_data/test_set_edge_cases/weighted_output \
                            --file test_set_edge_simple \
                            --shifts_file test_data/test_set_edge_cases/input/test_set_edge_simple_shift_df \
                            --get_Akita_scores \
                            --Akita_cell_types HFF H1hESC \
                            --roi genes \
                            --roi_scales 5 10
```


2. [Test set for sequences](https://github.com/ketringjoni/SuPreMo/blob/main/test_data/test_set_sequences/)

This is a set of 10 variants for which sequences were generated manually and used as ground truth to ensure that SuPreMo generates reference and mutated sequences accurately. The output sequences for this set from SuPreMo and from the manual curation are not included in this repo due to size, but were tested to be exactly the same.

  

***
*For any questions and/or feedback, please reach out to katie.gjoni at ucsf dot edu*



