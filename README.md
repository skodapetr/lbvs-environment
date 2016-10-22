# lbvs-environment
Ligand-based virtual screening (LBVS) environment (LBVS-environment) is designed to support researchers in task of development and benchmarking of chemical representations. The process of ligand-based virtual screening is divided into several phases: data preparation, screening and evaluation. This environment consists of a set of scripts enabling easy implementation of the above mentioned steps.

### Datasets preparation
Benchmarking datasets and selections can be downloaded using the `preparation/download.py` script.

### Screening
Given the method and specification of datasets performs a simulated LBVS.

### Evaluation
Evaluates quality of the individual methods used in the screening phase.

## Usage example
```
# Download datasets - Press '3' to download all, exit the script by pressing '4' when the download is finished.
python preparation/download.py
# Perform screening.
cd screening
python python_local.py -i ../data/datasets/10.1016%2Fj.ymeth.2014.11.015/selections/random_00_05_100_20_4900 -m ap/*
cd ..
# Evaluate data.
python evaluation/evaluation_local.py.py
```

