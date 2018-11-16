# K562 allele-specific CRISPR target detection

Identification of allele-specific CRISPR targets in the K562 cell line, as described in Zhou et al., "Comprehensive, integrated, and phased whole-genome analysis of the primary ENCODE cell line K562".

The code is provided as-is. A more user-friendly version may be released in the future.



## Usage

1. Prepare list of off-target sites:

```bash
python prepareListOfftargetSites.py <dir> <genome> <variants>
```
where "dir" is the folder where data is saved, "genome" is file containing the chromosome sequences for the reference genome (one whole chromosome per line), and "variants" is the vcf file for the phased variants.


2. Extract the allelle-specific sequences:

```bash
python extractSequences.py <dir> <variants> <genome> <output>
```
where "dir", "genome" and "variants" are as above, and "output" is the file where these sequences will be saved.


3. Identify and filter CRISPR sites for these sequences:

```bash
python targetExtraction.py <dir> <sequences> <BowtieIndex>
```
where "dir" is as above, "sequences" is the file generated as output from step 2, and "BowtieIndex" is the index file for Bowtie2. In this analysis, it should therefore be `hg19`.


## Citing

If you use this software, please cite: 

```
B. Zhou et al. (2018). Comprehensive, integrated, and phased whole-genome analysis of the primary ENCODE cell line K562. BioRxiv 192344.
```
(paper under review; final citation details to be updated later)

## License

See `LICENSE`

