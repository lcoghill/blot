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
