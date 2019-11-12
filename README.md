# OLDRIM: A Method for Overlapping Driver Module Identification in Cancer

### This is the original repository for the OLDRIM paper


**Installing the dependencies**

```
pip install -r requirements.txt
```

## **Input**

1. The PPI network file:

The file is located at data/intAct_PPI.txt

```
gene1 gene2 confidence value
MDM2  TP53  0.99
APP APP 0.99
MYC MAX 0.98
...
```
Note: this network contains the edges with confidence value > 0.35

2. Mutation data
The file contains the list of mutated genes. Each line consists of a gene and the list of patients where that gene is mutated
The file is located at data/gene_patients.txt

```
RNF14 TCGA-BH-A0BR TCGA-A5-A0GJ ...
UBE2Q1 TCGA-B0-5709 TCGA-AA-3549 ...
RNF17 TCGA-CN-4731 TCGA-D1-A0ZQ ...
....
```
This file is needed for calculating the mutex and coverage scores

3. The Ranked seed list

The file is located at data/intAct_seeds.txt
```
TP53
CUL9
NMT1
...
```
Note: this list is generated based on the network provided above and the pancancer mutation data


## **Seed list generation**
To generate the ranked seed list for a given network and patients mutation data,
the following script should be run first
```
cd src/seed_generation/
python generate_seed_file.py
```
Modify the appropriate parameters in the script to change the network,
the mutation data, and the output files path

## **Run**

To run OLDRIM on the given input files:

```
cd src
python oldrim.py
```
`src/config.py` contains different parameters for the inout files path and the functions used for scoring
