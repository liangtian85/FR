# Functional Redundancy in the Human Microbiome

A MATLAB package to quantify the functional redundancy of human microbial community

Liang Tian1,2, Xu-Wen Wang1, Ang-Kun Wu1,3, Yuhang Fan1,4, Jonathan Friedman5, Amber Dahlin1, Matthew K. Waldor6,7, George M. Weinstock8, Scott T. Weiss1 & Yang-Yu Liu1
1Channing Division of Network Medicine, Brigham and Women’s Hospital and Harvard Medical School, Boston, Massachusetts 02115, USA
2Department of Physics, Hong Kong Baptist University, Hong Kong, China
3Department of Physics and Astronomy, Rutgers University, Piscataway, New Jersey 08854, USA
4Department of Bioengineering, Stanford University, Stanford, California 94305, USA
5Department of Plant Pathology and Microbiology, Faculty of Agriculture, Food and Environment, The Hebrew University of Jerusalem, Jerusalem, Israel
6Division of Infectious Diseases, Brigham and Women’s Hospital and Harvard Medical School, Boston, Massachusetts 02115, USA
7Howard Hughes Medical Institute, Boston, Massachusetts 02115, USA
8The Jackson Laboratory for Genomic Medicine, Farmington, Connecticut 06117, USA 

FR_tutorial_1.0.zip  (1.6 MB, version 1.0)
Last update: Sep 06, 2020

The file FR_tutorial_1.0.zip contains: 
1) Abundance table and genomic content network for Stool sample of HMP (HMP_Stool.mat): Table of abundance at stain level of 553 stool samples from the HMP (http://www.hmpdacc.org/ ).
2) Matlab scripts (script_FR_GCN_randomization.m & script_FR_otu_randomization): The matlab code imports the abundance table and genomic content network, calculates the functional redundancy for the real abundance and genomic content network and also for the randomized counterpart, and then plot the figures.
3) Matlab functions used by the scripts.
Running the tutorial:
1) Extract the content of the enclosed FR_tutorial_1.0.zip file to a local directory.
2) Run the Matlab file. Running times for both of the scripts are less than 1 min. 
(The code was written on MATLAB R2016b)

Output:
1) Matlab figure showing normalized functional redundancy of real abundance and genomic content network compared with that for randomized genomic content network (left) and for randomized abundance table (right).
2) Mann-Whitney U test p values and adjusted q values (Benjamini-Hochberg).
 
![FR](https://github.com/liangtian85/FR/blob/master/FR_figure.png)
