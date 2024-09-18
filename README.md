# diverse-dataset-discovery

This repo contains tools for evaluating datasets and selecting molecules that contain chemistries
where we are seeking to improve our data coverage.

:rocket: We would love to get sent any selected molecules. These molecules represent areas of chemistry for which we currently have low data coverage in both training and benchmarking our force fields. Any contributions would directly help us improve our force fields! :rocket:

Otherwise, if you let us know of any matched groups and the number of matches
(which can be saved with the `-oc` flag), we would appreciate that in guiding our
curation of public databases.

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
curl https://raw.githubusercontent.com/openforcefield/diverse-dataset-discovery/main/bin/select-interesting-molecules_2.2.0_v1.py > select-interesting-molecules_2.2.0_v1.py
```

###  1.3. <a name='Runthetool'></a>Run the tool

The tool expects an input list of SMILES and outputs a list of SMILES. More options are documented in the help string (below).


For example:
```
python select-interesting-molecules_2.2.0_v1.py -i my-input-smiles.smi -o selected.smi -np 16 -oc counts.csv
```

An example run, example inputs, and outputs are shown in [example/](example/).

The molecules in `selected.smi` match at least one environment for which we are
seeking to increase coverage in our datasets. We search for these matches by either functional group, or force field parameters.

The `-oc` flag writes an output CSV of how many matches to each group were found.
This is also very useful information for us as it allows us to prioritize groups
for improving data coverage.

The `-np` flag controls how many processes to use.


####  1.3.1. <a name='Moreoptions'></a> More options

The tool comes with a number of options to control the output of the script and multiprocessing. For example, `-n 100` returns only the top 100 molecules.

These are listed in the help string as below.

```
$ python bin/select-interesting-molecules_2.2.0_v1.py --help
usage: select-interesting-molecules_2.2.0_v1.py [-h] -i INPUT [-o OUTPUT]
                                                [-n ONLY_TOP_N] [-np NPROC]
                                                [-oc OUTPUT_COUNT]
                                                [-of OUTPUT_FULL]
                                                [-c COUNT_THRESHOLD]

Select molecules with chemistries where we are seeking to improve data coverage.
This takes in a multi-molecule SMILES file and outputs a multi-molecule SMILES
file, where each molecule is on a separate line.
The chemistries we are selecting for include a list of functional groups and
openff-2.2.0.offxml parameters for which there is low coverage in our already
available data.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to a file containing SMILES strings, with one on
                        each line.
  -o OUTPUT, --output OUTPUT
                        Path to the output SMILES file. (Default: output.smi)
  -n ONLY_TOP_N, --only-top-n ONLY_TOP_N
                        Only write the top N molecules to the output file. If
                        not specified, write all molecules.
  -np NPROC, --nproc NPROC
                        Number of processes to use. (Default: 1)
  -oc OUTPUT_COUNT, --output-count OUTPUT_COUNT
                        If specified, write the counts of each group as a CSV
                        to the given path.
  -of OUTPUT_FULL, --output-full OUTPUT_FULL
                        If specified, write matches to low coverage groups as a
                        CSV to the given path. Each group will be a column,
                        with boolean values to indicate if this group is
                        present in the molecule. Each row will correspond to a
                        molecule in the input file. A column 'Count' will be
                        included to indicate the total number of matches. If
                        not specified, this file will not be written.
  -c COUNT_THRESHOLD, --count-threshold COUNT_THRESHOLD
                        Number of matches to groups with low existing data
                        coverage. Only molecules with a count greater than or
                        equal to this threshold will be written as output.
                        (Default: 1)
```

###  1.4. <a name='Contributingdataback'></a>Contributing data back

OpenFF would love to receive molecules that are matches, as they represent areas of chemistry where we currently have limited data coverage for training and benchmarking our force fields. Your contributions would directly help us improve our force fields by filling these gaps and improving their accuracy. 

Feel free to reach out to us on the Open Force Field Slack or by email.

Please note that all molecules that are sent to us will get released openly on the [MolSSI QCArchive](https://qcarchive.molssi.org/),
with a permissive license (typically [CC-BY](https://creativecommons.org/licenses/by/4.0/deed.en)).