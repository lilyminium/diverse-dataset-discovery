# Regenerating the scripts in this repo

This repo contains tools to re-generate the script.

## Installation

First create a development environment:

```
mamba env create -f devtools/conda-envs/dev.yaml
conda activate discoset-dev
pip install -e .
```

## Running

Then use the cli to generate a new file:

```
discoset generate -ff "openff-2.2.0.offxml" -o diverse-dataset-discovery/bin -v -np 16 -fft 10 -fgt 10 -vn 1
```

This will:

* download all available datasets from the QCArchive
* locate all unique molecules within the datasets
* label each molecule with functional groups (checkmol) or force field parameters
* find parameters or functional groups with low coverage
* save those to a new script file


## Testing

GitHub CI includes a broad matrix of tests. These currently range from Python 3.8 to 3.11,
cover OpenFF toolkits 0.10 onwards, and include either RDKit or OpenEye.
These have a `workflow_dispatch` trigger, i.e. can be run manually.
