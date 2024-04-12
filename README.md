# InfUSER
## About
`InfUSER` proposes tools for the analysis of Hi-C data. It currently only contains one fucnction, `singletree`, that aims to infer the state of progenitor cells, given a differentition tree and data of the corresponding data type at the leaves. It can take multiple data types as input, but will need all 
inputs of a single run to be of the same data type.

## Quick Start
[TODO]
### CLI
### API

### Cloning the repository
Use the following command to clone the repository: 
   
    git clone https://github.com/AudreyBaguette/InfUSER.git

### Installing `InfUSER` and its dependencies
[TODO]


## Usage
### info
Prints the version of the package and its source (this Git page).

### singletree
Run InfUSER with a single data type.
#### Required inputs
- Tree topology file (API: tree_path, CLI: TREEPATH)
    The tree topology file is a file that records the topology of the differentiation tree to use. The first row contains only one field, the name of the root. The following rows contain two fields, separated by tabs. The first one is the name of a new node, the second is the name of its parent. The tree is contrsucted from the root down, so the specification of parents must be written before its children. The tree is not binary, a parent can have more than two children. Each node name must be unique.
    Example:  
    root  
    n1 &nbsp; &nbsp;root  
    n2 &nbsp; &nbsp;root  
    l1 &nbsp; &nbsp;n1  
    l2 &nbsp; &nbsp;n1  
    l3 &nbsp; &nbsp;n1  
    l4 &nbsp; &nbsp;n2  
    l5 &nbsp; &nbsp;n2  
    Gives the following tree:  
    root  
    ├── n1  
    │ &nbsp; &nbsp; ├── l1  
    │ &nbsp; &nbsp; ├── l2  
    │ &nbsp; &nbsp; └── l3  
    └── n2  
    &nbsp; &nbsp; &nbsp; &nbsp;├── l4  
    &nbsp; &nbsp; &nbsp; &nbsp;└── l5
    The file corresponding to the example file above is provided at `Examples/tree_file.tsv`.
- Samples description file (API: sample_file, CLI: SAMPLEFILE)
    The samples files contains two fileds, separated by tabs. The first field is the name of the sample. Samples names must be the same as the one in the topology file. The second field is the path to the corresponding data file (.mcool or .bed). All leaves (terminal nodes) present in the topology file must be present in the samples file with a valid path. If extra samples are present in the samples files and not in the topology, they will be ignored.
    - For Hi-C:
    The samples file corresponding to tree topology file above is provided at `Examples/HiC_samples_file.tsv`.
    - For other data types:
    The data files need to have the same number of rows, in the same order. The samples file corresponding to tree topology file above is provided at `Examples/ChIP_samples_file.tsv`.
- Output directory (API: output_dir, CLI: OUTDIR)
    The path to the output directory (see Outputs section for a description of the created files and the file structure)
- File of chromosome sizes (API: chrom_sizes, CLI: CHROMSIZES)
    Optional, the path to the file containing the size (in bp) of each chromosome (default "data/hg38.chrom.sizes")
- Chromosomes to process (API: chromlist, CLI: CHROMLIST)
    Optional, the names of the chromosomes to consider. This list is ignored if subset is not null. 

#### Optional inputs
- res/r/resolution : int
    The resolution to consider. The input files need to have been generated with that resolution. (Hi-C only, default 10000)
- Subset file (specified by subset)
    The subset file is only used when the input data is Hi-C. It must contain three columns, separated by tabs. The first column is the chromosome name. The second column is the start of the region to consider. The third column is the end of the region to consider. The file must contain a header, with the following names:
	seqnames	start	end
	Notes: the start and end coordinates are rounded down to the nearest bin, relative to the resolution. The end bin is excluded. If the start and end region fall within the exact same bin, the region is considered too small and is ignored. An example is provided at `Examples/subset_file.tsv`.
- dist : int
    Optional, The distance to consider. All interactions beyond that distance will be ignored. If set to 0, all interactions are kept (Hi-C only, default 0)
- n_values : int
    The number of values that need to be stored in the nodes. In other words, how many values should be considered to discretize the continuous data. (default 9)
- min : float
    The minimal Z-score value to consider (default -4)
- max : float
    The maximal Z-score value to consider (default 4)
- column : int
    Optional, the column conting the score to consider. The first column is column 1 (1D data only, default 4)
- transform
- balance
- n_jobs

#### Outputs
The output folder will contain one file and three sub-folders:

- tree_structure.txt
This files records the tree topology constructed from the input file.

- internal_nodes
This folder contains the infered data for the internal nodes of the differentiation tree. It contains one sub-folder per node. The names will match the names in the input files.
For Hi-C:
Each sample sub-folder will contain one file per (computed) chomosome and one .cool file merging them all.
For other data types:
Each sample sub-folder will ontain one tsv file. This only contains one unnamed column, containing the inferred values in the same order as the input file.

- leaves
This folder contains the transformed data for the leaves (terminal nodes) of the differentiation tree. It contains one sub-folder per node. The names will match the names in the input files.
For Hi-C:
Each sample sub-folder will contain one file per (computed) chomosome and one .cool file merging them all. The .cool files are different from the input files in that they are limited to one resolution, may not contain data genomewide (depending on the chromlist and subset_file), and the contact values went through the same Z-transformed, discretization and de-transformation as the internal nodes.
For other data types:
Each sample sub-folder will ontain one tsv file. This only contains one unnamed column, containing the inferred values in the same order as the input file. The .tsv files are different from the input files in that their contact values went through the same Z-transformed, discretization and de-transformation as the internal nodes.

- edges
This folder contains the changes between each parent-child pair. Those changes are computed as differences between Z-scores, before de-transformation. The folder contains one sub-folder per edge in the topology. The names will match the names in the input files and record the direction that was considered (from parent to child).
For Hi-C:
Each sub-folder will contain one file per chomosome and one .cool file merging them all. The files do not contain a contact frequency, but differences in contact frequencies.
For other data types:
Each sub-folder will contain one tsv file. The files do not contain a signal value, but differences in signal values.


## Contributing
### Contributors
- Audrey Baguette
- Tunde Lapohos

## References

## Citing `InfUSER_single_tree`
[TODO]
