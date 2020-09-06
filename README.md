# Functional Redundancy in the Human Microbiome

A MATLAB package to quantify the functional redundancy of human microbial community

Liang Tian<sup>1,2</sup>, Xu-Wen Wang<sup>1</sup>, Ang-Kun Wu<sup>1,3</sup>, Yuhang Fan<sup>1,4</sup>, Jonathan Friedman<sup>5</sup>, Amber Dahlin<sup>1</sup>, Matthew K. Waldor<sup>6,7</sup>,  George M. Weinstock<sup>8</sup>, Scott T. Weiss<sup>1</sup> & Yang-Yu Liu<sup>1</sup>

<sup>1</sup>Channing Division of Network Medicine, Brigham and Women’s Hospital and Harvard Medical School, Boston, Massachusetts 02115, USA  
<sup>2</sup>Department of Physics, Hong Kong Baptist University, Hong Kong, China  
<sup>3</sup>Department of Physics and Astronomy, Rutgers University, Piscataway, New Jersey 08854, USA  
<sup>4</sup>Department of Bioengineering, Stanford University, Stanford, California 94305, USA  
<sup>5</sup>Department of Plant Pathology and Microbiology, Faculty of Agriculture, Food and Environment, The Hebrew University of Jerusalem, Jerusalem, Israel  
<sup>6</sup>Division of Infectious Diseases, Brigham and Women’s Hospital and Harvard Medical School, Boston, Massachusetts 02115, USA  
<sup>7</sup>Howard Hughes Medical Institute, Boston, Massachusetts 02115, USA  
<sup>8</sup>The Jackson Laboratory for Genomic Medicine, Farmington, Connecticut 06117, USA 

Last update: Sep 06, 2020

The files contain: 
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
