blot
====
### Branch Lengths for Open Tree ###

Tool to semi-automate incorporating raw branch lengths into subtrees from the opentree of life project. 

---
#### Requirements ####
The following are the package requirements to run the scripts, and how to install them. (1. is for Debian / Ubuntu, 2. is for Fedora / Redhat)
```
1. sudo apt-get install
2. sudo yum -y install
```

#### Example Walkthrough ####
This assumes some type of Linux environment with access to a package manager like apt-get or yum.

The first step is to download the scripts and the example data. You can acheive this with the following command:
``` 
git clone https://github.com/lcoghill/blot
```
In this case the example data is a tree of Centarchidae downloaded from the [tree.opentreeoflife.org](https://tree.opentreeoflife.org) web application. You can replicate this by navigating to the page, typing Centarchidae in the search box on the upper right side of the page, and then clicking on the link that says "[Download subtree as Newick string](https://tree.opentreeoflife.org/opentree/argus/ottol@437610/Centrarchidae)".

You should end up with a newick tree that has a format like:
```
((Acantharchus_pomotis_ott476356)Acantharchus_ott476355,...
```
In this case, the format is Genus_species_opentreeIDnumber.

Next, we need to set up a corrections file. This is just a simple file in csv format that makes the script aware of any changes that need to be applied to taxon names. The OTT (open tree taxonomy) is a amalgamation of many other taxonomies into a single source. This can create a problem with fetching data from Genbank as it will only recognize scientific names that are present in the NCBI taxonomy. There is a script to help with this, called validate-data.py. You can run this script with the tree file and target gene, and it will let you know if there are taxa that aren't found, or any taxa for which your target gene has no data present in Genbank. 

```
python validate-data.py
```
If there is data present for all your taxa on the target gene, you are set to move on. If not, you can run this script multiple times and use it to help guide building your corrections.csv file. 



Now that we have the data, we need to set some basic parameters in the script: 
```
gene = ""
email = ""
seqs_per_taxon = X
tree_file_name = ""
corrections_file_name = ""
raxml_dir = ""
raxml_constraint = ""
max_muscle_iters = ""
max_muscle_trees = ""
```
