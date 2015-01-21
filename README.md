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

##### Downloading and setting up the pipeline #####
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

##### Building a RAxML constraint tree #####

##### Validating data and building a corrections file #####

Next, we need to set up a corrections file. This is just a simple file in csv format that makes the script aware of any changes that need to be applied to taxon names. The OTT (open tree taxonomy) is a amalgamation of many other taxonomies into a single source. This can create a problem with fetching data from Genbank as it will only recognize scientific names that are present in the NCBI taxonomy. There is a script to help with this, called validate-data.py. You can run this script with the tree file and target gene, and it will let you know if there are taxa that aren't found, or any taxa for which your target gene has no data present in Genbank. 

```
python validate-data.py
```
If there is data present for all your taxa on the target gene, you are set to move on. If not, you can run this script multiple times and use it to help guide building your corrections.csv file. If you run this on the raw tree we downloaded without any corrections, your output should resemble something like:
```
Looking for newick tree file... Found. 
Parsing newick file centarchidae.tre... Done.
Checking for a corrections file... 
No corrections file found. Proceeding without it.
Checking that NCBI has records for all target taxa...
Some terminal taxa don't have data available in Genbank for the gene ND2.
Please check the gene name for errors or double check that data is available in Genbank.
----------------------------------------------------------------------------------------------------
Taxon                              	 Taxon Present 	 Data Present 	 # of Seqs
-------------                       	 ------------- 	 ------------- 	 -------------
Ambloplites constellatus            	 True         	 True              	 1
Lepomis miniatus                    	 True         	 True              	 2
Chaenobryttus gulosus               	 True         	 True              	 2
Micropterus punctulatus punctulatus 	 False         	 False               0
Micropterus cf coosae MRB-2014      	 True         	 True              	 16
Lepomis cyanellus                   	 True         	 True              	 2
Micropterus dolomieu                	 True         	 True              	 18
Acantharchus pomotis                	 True         	 True              	 1
Pomoxis nigromaculatus              	 True         	 True              	 3
Ambloplites cavifrons               	 True         	 True              	 1
Pomoxis annularis                   	 True         	 True              	 2
Micropterus cf coosae REB-2012a     	 True         	 True              	 5
Micropterus cf coosae REB-2012b     	 True         	 True              	 2
Micropterus warriorensis            	 True         	 True              	 4
Micropterus cahabae                 	 True         	 True              	 4
Micropterus cataractae              	 True         	 True              	 5
Micropterus coosae                  	 True         	 True              	 13
Micropterus floridanus              	 True         	 True              	 20
Micropterus pallidus                	 False         	 False               0
Poxomis nigromaculatus              	 False         	 False               0
Micropterus henshalli               	 True         	 True              	 5
Lepomis pallidus                    	 False         	 False               0
Lepomis punctatus                   	 True         	 True              	 1
Centrarchus macropterus             	 True         	 True              	 1
Lepomis auritus                     	 True         	 True              	 1
Xenotis megalotis                   	 False         	 False               0
Enneacanthus chaetodon              	 True         	 True              	 1
Micropterus salmoides salmoides     	 True         	 True              	 2
Micropterus chattahoochae           	 True         	 True              	 3
Lepomis marginatus                  	 True         	 True              	 1
Micropterus tallapoosae             	 True         	 True              	 4
Micropterus treculii                	 True         	 True              	 2
Micropterus notius                  	 True         	 True              	 4
Lepomis macrochirus                 	 True         	 True              	 5
Ambloplites ariommus                	 True         	 True              	 1
Ambloplites rupestris               	 True         	 True              	 1
Enneacanthus obesus                 	 True         	 True              	 1
Lepomis microlophus                 	 True         	 True              	 2
Lepomis megalotis                   	 True         	 True              	 20
Archoplites interruptus             	 True         	 True              	 1
Lepomis humilis                     	 True         	 True              	 2
Enneacanthus gloriosus              	 True         	 True              	 1
Lepomis peltastes                   	 False         	 False               0
Lepomis gibbosus                    	 True         	 True              	 1
Lepomis symmetricus                 	 True         	 True              	 1
----------------------------------------------------------------------------------------------------

```
The output can be helpful for us to design our corrections file. For example, we can see that several of the taxa are not in the NCBI taxonomy, and subsequently don't have sequences available in Genbank. We can develop a corrections file that will tell the script how to proceed in these cases. For example, "Lepomis pallidus" is also and more commonly known as Lepomis macrochirus (bluegill). Chances are there are sequences in Genbank for this taxon, but under L. macrochirus instead of L. pallidus. In order to correct this, we can create a csv file and add a line with the following: 
```
Lepomis pallidus, Lepomis machrochirus, organism
```
what this is telling the script is that, where it would normally search for L. pallidus, replace that search looking for data for L. machrochirus. The third parameter, "organism" is telling the script to search genbank with the organism keyword. This will help it retreive any hits that exists for this species. We could also pass it an accession number if we want to use a specific sequence for a given taxa:
```
Lepomis pallidus, KF571702.1, accession
```
this will force the script to only validate that the given accession number is present in genbank. If it is, that taxon will pass validation. 

This will likely be the most time consuming part of the pipeline. You'll likely have to do some literature research, or have an extensive knowledge of your group, in order to build a helpful and accurate corrections file. In the end, some calls may be subjective in nature as well. If you are concerned about that, you can always tweak the corrections file and run it multiple times to inspect your results.

The final item to check in the validation output is the estimated number of sequences in genbank for a given taxon and gene. We can see in our example output above that many of a taxa only have 1 or 2 sequences availabe for the ND2 gene. Some have 0, but it is likely once we build the corrections file, some data will be present for all taxa. However, keep this in mind when suggestion the number of sequences per taxa to use when running blot. The script will run, even if your desired number of sequences per taxa aren't present, however, it will throw some errors up along the way.

##### Building the final tree #####

Once your data (raw tree, corrections file and constraint tree) are ready, we can actually proceed to building the tree with branch lengths.

First, we need to set some basic parameters in the script as before.
Below are the suggestions for the example data.


This is the target gene to fetch data for:

```gene = "ND2"```

Your email, used for the NCBI Entrez calls in case there are any problems.

```email = "YOU@someplace.edu"```

Desired number of sequences to use per taxon. Ultimately these will be combined into a single consensus sequence. If you have a desired sequence for a given taxa, you could specific the accession in the corrections file.

```seqs_per_taxon = 3```

Name / location of the raw newick tree file.

```tree_file_name = "examples/centarchidae.tre"```

Name / location of the corrections file.

```corrections_file_name = "examples/corrections.csv"```

Name / location of a directory to save the RAxML output trees.

```raxml_dir = "raxml/"```

Name / location of the RAxML constraint tree we generated earlier.

```raxml_constraint = "examples/constraint_tree.tre"```

Number of iterations for the Muscle alignment. 

```max_muscle_iters = "3"```

Number of trees to build per Muscle iteration.

```max_muscle_trees = "2"```
