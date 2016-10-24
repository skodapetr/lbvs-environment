## Datasets
Datasets are not stored in the GitHub repository due to their size. 
They can be downloaded using ```preparation/download.py``` Python 
script.

### Structure
Once downloaded each dataset is store in a directory in this folder. 
All datasets follow the same structure.

A dataset consists of molecules and selections. The molecules are 
stored in the SDF format in ```molecules/sdf directory```.

Selection are stored in the ```selection``` directory sub-tree. The 
first level corresponds to the method used to create the selection. 
For example the ```random_00_05_100_20_4900``` state for selection
created by random selection with 5 actives and 100 inactives in 
the train set and 20 actives and 4900 inactives in the test set.
The directory of selection methods contains ```targets``` 
(```groups```). In common case the dataset contains collections of 
actives and inactives against different targets. There collections 
are represented by ```groups```. In each ```groups``` are two types
of files group definition file and selection files. The group 
definition fil (```groups.json```) holds definition of activity for 
molecules used in the test set of given groups. The selection
files (test instances, ```s_*.json```) contains definition of train 
and test set.

