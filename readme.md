# PMC Detection in 3D confocal images

This is a repository to perform automatic segmentation and detection of primary mesenchyme cells (PMCs) in 3D confocal images. The workflow uses [snakemake](https://snakemake.readthedocs.io/en/stable/) and assumes the existence of a trained [Ilastik](https://www.ilastik.org/) model for PMC pixel classification.

The workflow begins by preprocessing immunolabeled PMC stains, performs semantic segmentation to identify individual cells, and then quantifies gene expression in each cell for any sm-FISH channel in the original confocal image.

FISH signal is quantified using two methods: first, for absolute quantification, single-molecules are counted using [Big-FISH](https://github.com/fish-quant/big-fish); second, relative quantification is performed by calculating the average intensity of the FISH signal within each PMC as Z-scores. While absolute quantification via spot counting is ideal, in practice we found this approach difficult if the given gene is expressed to a high enough degree that individual puncta are not easily discernable. The final output is a table of both PMC centroid coordinates, as well as gene expression measures for each cell and gene in the provided images. This will be written to a `.csv` file located at "<outdir>/final/counts.csv", where "<outdir>" is specified by the user.

## Installation

To install, clone this repository and ensure Snakemake, Ilastik, and [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) are installed. View the linked resources for how to best install each tool.

## Running the Pipeline

To run the pipeline, simply execute the command:

```
snakemake --use-conda --conda-frontend <conda/mamba> -j <desire number of jobs>
```

What files are processed and how they're processed is dictated by the configuration file. An example configuration file is found [here](files/config.yaml). Once cloned, you can either overwrite the values in this configuration file for your needs, or write your own configuration file. To make your configuration file accessible to snakemake, you will then either need to edit the `configfile` file name in the `Snakefile`, or supply your new configuration file at run time using the `--configfile` [parameter](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

## Configuration Parameters

The expected keys and values for workflow configuration are explained below:

```{yaml}
ilastik:
  model: # file path to ilastik model
  loc: # file path to run_ilastik.sh script 
input:
  logfile: # file path to log file containing meta data for each image file
  datadir: # directory containing '.nd2' 3D confocal images
output:
  dir: # desired output directed
quant:
  genes: # key-value pairs for each gene in the dataset that needs quantification 
    sm50: # gene name
      radius: [400, 350, 350] # expected radius of point spread function for spot counting.
    pks2: # gene name
      radius: [400, 350, 350] # expected radius of point spread function for spot counting.
```

## Image log file 

As noted above, the workflow expects a log file containing meta data fore ach image file in the dataset. This should be a tabular `.csv` file, and while the log file can contain any number of columns, the workflow **requries** the followings:

1. An index in the first column mapping file names for files in `datadir` to their meta data in the log file
2. `channel_order`: a column containing ";" delimited names for each channel in the image (e.g. "pmc;sm50;pks2")
3. `z_start`: z-slice where PMC signal starts, 0 and closed indexed.
4. `z_end`: z-slice where PMC signal ends, 0 and opened indexed (e.g. if slice 15 has the last PMC signal, `z_end` should be set to `16`)

## Implementation

For complete implementation details, see our paper on ICAT [paper](https://doi.org/10.1093/bioinformatics/btad278)

## Citation

If you used this pipeline in your work, please cite the original manuscript the workflow was developed for:

Dakota Y Hawkins, Daniel T Zuch, James Huth, Nahomie Rodriguez-Sastre, Kelley R McCutcheon, Abigail Glick, Alexandra T Lion, Christopher F Thomas, Abigail E Descoteaux, W Evan Johnson, Cynthia A Bradham, ICAT: A Novel Algorithm to Robustly Identify Cell States Following Perturbations in Single Cell Transcriptomes, Bioinformatics, 2023;, btad278, https://doi.org/10.1093/bioinformatics/btad278

### bibtex

```
@article {Hawkins2022.05.26.493603,
	author = {Dakota Y. Hawkins and Daniel T. Zuch and James Huth and Nahomie Rodriguez-Sastre and Kelley R. McCutcheon and Abigail Glick and Alexandra T. Lion and Christopher F. Thomas and Abigail E. Descoteaux and W. Evan Johnson and Cynthia A. Bradham},
	title = {ICAT: A Novel Algorithm to Robustly Identify Cell States Following Perturbations in Single Cell Transcriptomes},
	elocation-id = {2022.05.26.493603},
	year = {2023},
	doi = {10.1101/2022.05.26.493603},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Motivation The detection of distinct cellular identities is central to the analysis of single-cell RNA sequencing experiments. However, in perturbation experiments, current methods typically fail to correctly match cell states between conditions or erroneously remove population substructure. Here we present the novel, unsupervised algorithm ICAT that employs self-supervised feature weighting and control-guided clustering to accurately resolve cell states across heterogeneous conditions.Results Using simulated and real datasets, we show ICAT is superior in identifying and resolving cell states compared to current integration workflows. While requiring no a priori knowledge of extant cell states or discriminatory marker genes, ICAT is robust to low signal strength, high perturbation severity, and disparate cell type proportions. We empirically validate ICAT in a developmental model and find that only ICAT identifies a perturbation-unique cellular response. Taken together, our results demonstrate that ICAT offers a significant improvement in defining cellular responses to perturbation in single-cell RNA sequencing data.},
	URL = {https://www.biorxiv.org/content/early/2023/03/04/2022.05.26.493603},
	eprint = {https://www.biorxiv.org/content/early/2023/03/04/2022.05.26.493603.full.pdf},
	journal = {bioRxiv}
}
```
