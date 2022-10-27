# mitowrap
```
----------------------------------------------
           _ _                                
 _ __ ___ (_) |_ _____      ___ __ __ _ _ __  
| '_ ` _ \| | __/ _ \ \ /\ / / '__/ _` | '_ \ 
| | | | | | | || (_) \ V  V /| | | (_| | |_) |
|_| |_| |_|_|\__\___/ \_/\_/ |_|  \__,_| .__/ 
                                       |_|        
----------------------------------------------
```

A snakemake pipeline wrapping MitoZ and getOrganelle for de novo mitogenome assembly using short reads and subsequent QC. 

## Installation
Firstly, clone this repository:

```
git clone https://github.com/charlesfoster/mitowrap.git

cd mitowrap
```

Create a new conda environment:

```
# install mamba if not already installed
conda install -c conda-forge mamba
mamba env create -f environment.yml
```

If you are using Linux and want to run programs using Singularity containers (see sections below), you will also need to take the following step:

```
conda activate mitowrap
mamba install -c conda-forge singularity==3.7.1
```

## Usage
Each time you wish to run the program, make sure to activate the conda environment first:

```
conda activate mitowrap
```

Then, update the paths and parameters in `config.yaml` as you see fit. At minimum, you will likely want to update:

* reference: the path to a file with reference mitochondrial genome(s) in fasta format. Used to extract mitochondrial reads from your input data.
* suffix: the suffix in your reads filenames used to determine sample names properly. For example, if your forward reads are named "Sample1_R1_001.fastq.gz" and "Sample_R1_001.fastq.gz", respectively, then a suffix of "_R1_001.fastq.gz" will allow `snakemake` to recognise the sample name as "Sample1".
* reads_dir: the path to the directory containing all sequencing reads to be used in the analysis
* outdir: the path to the directory where you want your results to be saved

Minimal instructions to run:

```
snakemake -j 16
```

Note 1: Update the value after `-j` to reflect how many threads/cores you would like to allow `snakemake` to use.

Note 2: This method assumes all dependencies are installed correctly in your path. An easier option is to let the `snakemake` pipeline take care of all software dependencies for you by creating internal conda environments or using containers with Singularity.

Note 3: This method assumes all paths etc. for your run are defined in `config.yaml`. You can choose to have several config files instead, each with different names. If you take this path, you will always need to append the config file's name (e.g., 'other_config.yaml') to your `snakemake` command, e.g. `--configfile other_config.yaml`.

### Using conda
Running the pipeline with conda is currently more stable than with Singularity, and is supported on both Mac OS and Linux. Usage:

```
snakemake -j 16 --configfile config.yaml --use-conda 
```

### Using singularity
*Note* currently experiencing some issues with the `getOrganelle` container. Try using `conda` instead.

Running the pipeline with Singularity is only supported on Linux. Usage:

```
snakemake -j 16 --use-singularity --singularity-args
```

Sometimes if you have necessary data stored on mounted drives, you need to tell Singularity to bind those drives. Example:

```
snakemake -j 16 ---use-singularity --singularity-args '--bind /home:/home,/data/data'
```

### Other snakemake options
There are plenty of other options available within `snakemake` that I haven't touched on here. For example, you can save an image of the directed acyclical graph describing the workflow:

```
snakemake -c 1 --dag | dot  -Tpdf > dag.pdf
```

Note: requires external `dot` program

You can also choose to print all commands being run by `snakemake` during the workflow by appending one simple flag to your command:

```
snakemake -j 16 --use-conda --printshellcmds
```

For all other options, look at `snakemake -h` or the `snakemake` website.

## What does the pipeline do?
For each sample:

1. Raw reads are quality/adapter trimmed with `fastp`
2. Clean reads are mapped to the specified reference mitochondrial genome, and mitochondrial reads are extracted
3. Reads are de novo assembled into mitogenomes using (a) `MitoZ`, and (b) `getOrganelle` (then annotated using `MitoZ`)
4. QC metrics are calculated for each sample, and saved in an overall file

## What do you need to run the program?
Check the config file. At minimum:
* Paired-end short sequencing reads for one or more samples
* File with mitochondrial genomes in fasta format
* An up to date `config.yaml` file

## Citations
If you use this program and find it useful, I'd appreciate some kind of attribution, such as a link to this GitHub repo. Please also cite the programs used within this pipeline.
