# lbvs-environment
Ligand-based virtual screening (LBVS) environment (LBVS-environment) is designed to support researchers in task of development and benchmark of chemical representations. The process of ligand based virtual screening is divided into sevaral phases: data preparation, screening, evaluation and analysis. 

### Datasets preparation
Datasets and selections can be downloaded using the `preparation/download.py` script. 

### Screeening
Given the method and specification of data perform a simulated LBVS.

### Evalution
Compute evaluation methods of the results of the screening pahse.

### Analysis
Provide aggregation of evaluated results into a single file.

## Usage example
```
# Download datasets - Press '3' to download all, exit the script by pressing '4' when the download is finished.
python preparation/download.py
# Perform screening.
python screening/screening.py -f ./collections/auc_80-85.dat -M tt_tanimoto
# Evaluate data.
python evaluation/evaluation.py
# Write averega AUC perforamnce to the CSV table.
python analysis/saveAsCsv.py -f ./collections/single.dat -o performance.csv -n AUC -t avg
```


