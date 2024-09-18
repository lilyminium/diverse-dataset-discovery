# diverse-dataset-discovery

This repo contains tools for evaluating datasets and selecting molecules that contain chemistries
where we are seeking to improve our data coverage.

:rocket: We would love to get sent any selected molecules. These molecules represent areas of chemistry for which we currently have low data coverage in both training and benchmarking our force fields. Any contributions would directly help us improve our force fields! :rocket:

**Contents**

<!-- vscode-markdown-toc -->
* 1. [How to use](#Howtouse)
	* 1.1. [Installation (optional)](#Installationoptional)
	* 1.2. [Download the tool](#Downloadthetool)
	* 1.3. [Run the tool](#Runthetool)
		* 1.3.1. [ More options](#Moreoptions)
	* 1.4. [Contributing data back](#Contributingdataback)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

##  1. <a name='Howtouse'></a>How to use

The tools in [bin](bin/) will work in any Python 3.9+ environment with the [OpenFF Toolkit](https://docs.openforcefield.org/projects/toolkit/en/stable/) v0.10+ and *either* a working installation of [RDKit](https://www.rdkit.org/) or [OpenEye](https://www.eyesopen.com/). (However, we generally
recommend using up-to-date software as much as possible).

###  1.1. <a name='Installationoptional'></a>Installation (optional)

If you do not have an OpenFF environment already set up, create one with the below command. (If you prefer the conda or micromamba package managers, replace `mamba` with the correct executable.)

```
mamba create --name openff-toolkit -c conda-forge "openff-toolkit>0.10.0"
mamba activate openff-toolkit
```

###  1.2. <a name='Downloadthetool'></a>Download the tool

Download the single-file script either by cloning the GitHub repo:

```
git clone https://github.com/openforcefield/diverse-dataset-discovery.git
cd diverse-dataset-discovery/bin
```

or just download the file:

```
curl https://raw.githubusercontent.com/openforcefield/diverse-dataset-discovery/main/bin/sort-by-rare-groups_2.2.0_v1.py > sort-by-rare-groups_2.2.0_v1.py
```

###  1.3. <a name='Runthetool'></a>Run the tool

The tool expects an input list of SMILES and outputs a CSV. More options are documented in the help string (below).


For example:
```
python select-interesting-molecules_2.2.0_v1 -i my-input-smiles.smi -o selected.smi -np 16
```

An example run, example inputs, and outputs are shown in [example/](example/).

The molecules in `selected.smi` match at least one environment for which we are
seeking to increase coverage in our datasets. We search for these matches by either functional group, or force field parameters.


####  1.3.1. <a name='Moreoptions'></a> More options

The tool comes with a number of options to control the output of the script and multiprocessing. For example, `-n 100` returns only the top 100 molecules.

These are listed in the help string as below.

```
$ python bin/sort-by-rare-groups_2.2.0_v1.py --help
usage: sort-by-rare-groups_2.2.0_v1.py [-h] [-i INPUT] [-o OUTPUT] [-np NPROC] [--no-write-all] [-c COUNT_THRESHOLD] [-n ONLY_TOP_N]

Sort a list of SMILES strings by the number of matches to rare environments. This takes in a list of SMILES strings and outputs a CSV file with the SMILES strings and the number of matches to rare environments.
The rare environments include a list of functional groups and openff-2.2.0.offxml parameters for which there is low coverage in our already available data.
The output CSV file will contain columns for SMILES strings (column: 'SMILES', type: str) and the number of matches to rare environments (column: 'Count', type: int). If the flag --no-write-all is set, only these two columns will be written. If the flag is not set, all the rare environments that were searched for will be written as well, with boolean values to indicate if they are present or not.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to a file containing SMILES strings, with one on each line.
  -o OUTPUT, --output OUTPUT
                        Path to the output CSV file.
  -np NPROC, --nproc NPROC
                        Number of processes to use.
  --no-write-all        Write all columns to the CSV file. If False, only columns for SMILES and the number matches to rare groups will be
                        written. If True, all the rare groups that were searched for will be written.
  -c COUNT_THRESHOLD, --count-threshold COUNT_THRESHOLD
                        Threshold for considering a parameter 'rare'. Only groups with a count greater than or equal to this threshold will be
                        written.
  -n ONLY_TOP_N, --only-top-n ONLY_TOP_N
                        Only write the top N entries to the output file.
```

###  1.4. <a name='Contributingdataback'></a>Contributing data back

OpenFF would love to receive molecules that are matches, as they represent areas of chemistry where we currently have limited data coverage for training and benchmarking our force fields. Your contributions would directly help us improve our force fields by filling these gaps and improving their accuracy. 

Feel free to reach out to us on the Open Force Field Slack or by email.

Please note that all molecules that are sent to us will get released openly on the [MolSSI QCArchive](https://qcarchive.molssi.org/),
with a permissive license (typically [CC-BY](https://creativecommons.org/licenses/by/4.0/deed.en)).