# RoDTIS: Robust DNA Targeting *in silico*
#### A script that returns the gene sequences of interest in an electronic PCR-like manner. 
This script was used in Nielsen *et al.* (2018). Working title: Association of Outer Membrane Protein A.

<br>

## Assumptions

 * Python 3 is installed on machine. [Download python](https://www.anaconda.com/download/#macos).
 * Biopython(https://biopython.org/wiki/Download) is installed on machine. 
 * User has bash terminal. 

<br>

## Input
##### Standard: 
1. Create an Excel sheet with the gene name, primers, and size range of product as follows:
E.g.: ![Input](panel1.png)

2. ```RoDTIS.py Panel1.xlsx fasta_files/* > Output.fasta```


##### Advanced: 
RoDTIS allows for variable matching of the last half of the primer used, which sets it apart for simply matching the primer in the sequence. 

 * To use the advanced option, use the ```-s``` flag. 
 	* The ```-s``` flag must be used with an interval between 0.0-1.0.
 		* The higher the number:
 		* The lower the number: 
 	* e. g. ```RoDTIS.py -s 0.2 Panel1.xlsx fasta_files/* > Output.fasta```

<br>

## Output
Open FASTA file in the program of your choice, and view the genes within each fasta file.


<br>
