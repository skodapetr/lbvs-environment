# lbvs-environment
Ligand-based virtual screening (LBVS) environment (LBVS-environment) is 
designed to support researchers in task of development and benchmarking of 
chemical representations. The process of LBVS is divided into several phases: 
* data preparation
* screening 
* evaluation. 
This environment consists of a set of scripts enabling easy implementation 
of the above mentioned steps.

## Usage examples 

### Ligand-based virtual screening
Benchmark performance of `tt_tanimoto` on the `01` collection.  

```
cd scripts
# Download datasets - Choose download all, when done exit the script.
python download_data.py
# Perform screening.
python screen_localhost.py -c 01 -m tt/tt_tanimoto
# Evaluate data.
python evaluate_screening.py
python export_results.py
```

### Molecule export
Export molecules (test, train, validation set) used in collection `01`.

```
cd scripts
# Download datasets - Choose download all, when done exit the script.
python download_data.py
python export_molecules.py -c 01 -o ../molecule-export --flat --info
```
