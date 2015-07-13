import re
import os
import sys
import getopt
from ete2 import Tree
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO as SIO
from termcolor import colored



def check_data(gene, taxa_names, email, seqs_per_taxon, corrections, max_length) :
## Check that NCBI has data for the target gene for all leafs in tree
    print "Checking that NCBI has records for all target taxa..."
    ncbi_data = {}
    taxon_data = {}
    no_data = False
    count = 1
    est_seq_count = {}
    for t in taxa_names :
        seqs_present = True
        result, est_seqs = check_seqs(t, gene, email, seqs_per_taxon, ott_ids[t], corrections, max_length)
        taxon_result = check_taxon(t, email, corrections)
        est_seq_count[t] = est_seqs
        if result == 0 :
            seqs_present = False
            no_data = True
        ncbi_data[t] = seqs_present
        if taxon_result == 0 :
            taxon_data[t] = False

    if no_data == False :
        print "All terminal taxa appear to have data available in Genbank for gene %s." %gene
    elif no_data == True :  
        print "Some terminal taxa don't have data available in Genbank for the gene %s." %gene
        print "Please check the gene name for errors or double check that data is available in Genbank."
        print "-"*100
        table_fmt = [5, 13, 12, 9]
        long_key = 0
        for key in ncbi_data.keys() :
            if len(key) > long_key :
                long_key = len(key)
            
        
        print "Taxon%s\t Taxon Present \t Data Present \t # of Seqs" %(" "*(long_key-len("Taxon")))
        print "%s \t %s \t %s \t %s" %("-"*13+" "*(long_key-13),"-"*13,"-"*13,"-"*13)
        for key, val in ncbi_data.items() :
            chars = len(key)
            fmt_spaces = " " * (long_key - chars)

            if val == False and ncbi_data[key] == False :
               print "%s%s \t %s%s \t %s%s \t %s" %(key, fmt_spaces, colored(ncbi_data[key], 'red'), " "*8, colored(val, 'red'), " "*13,  est_seq_count[key])    
            elif ncbi_data[key] == False :
               print "%s%s \t %s%s \t %s%s \t %s" %(key, fmt_spaces, colored(ncbi_data[key], 'red'), " "*8, val, " "*13,  est_seq_count[key])
            elif val == False :
                print "%s%s \t %s%s \t %s%s \t %s" %(key, fmt_spaces, ncbi_data[key], " "*8, colored(val, 'red'), " "*13,  est_seq_count[key])
            else :
               print "%s%s \t %s%s \t %s%s \t %s" %(key, fmt_spaces, colored(ncbi_data[key], 'green'), " "*8, colored(val, 'green'), " "*13, est_seq_count[key])
        
        print "-"*100
        
def check_seqs(t, gene, email, seqs_per_taxon, ott_id, corrections, max_length) :
    seq_present = 0 
    Entrez.email = email
    ## get list of gi values that match organism and gene
    if t in corrections.keys() :
        if corrections[t].split(",")[1] == "organism" :
             term = corrections[t].split(",")[0] + "[Organism] AND " + gene + "[Gene]"
        elif corrections[t].split(",")[1] == "accession" :
             term = corrections[t].split(",")[0] + "[Accession]" 
    else :
        term = t + "[Organism] AND " + gene + "[Gene]"

    handle = Entrez.esearch(db="nucleotide", term=term)
    record = Entrez.read(handle)
    
    if len(record['IdList']) > 0 :
        ## get seq length here
        #for each id, get seq length here if there is at least one present < 2000, set value to 1. 
        for gi in record['IdList'] :
            if seq_present == 0 :
                handle = Entrez.efetch(db="nucleotide", id=gi ,rettype="gb", retmode="text")
                gbrecord = SeqIO.read(handle, "genbank")
                if len(gbrecord.seq) <= max_length :
                    seq_present = 1   
        
    else :
        seq_present = 0

    est_seqs = len(record['IdList'])
    return seq_present, est_seqs

def check_taxon(t, email, corrections) :
    seq_present = 0 
    Entrez.email = email
    ## get list of gi values that match organism and gene
    if t in corrections.keys() :
        if corrections[t].split(",")[1] == "organism" :
             term = corrections[t].split(",")[0] + "[Organism]"
        elif corrections[t].split(",")[1] == "accession" :
             term = corrections[t].split(",")[0] + "[Accession]" 
    else :
        term = t + "[Organism]"

    handle = Entrez.esearch(db="taxonomy", term=term)
    record = Entrez.read(handle)
    
    if len(record['IdList']) > 0 :
        seq_present = 1   
        
    else :
        seq_present = 0

    return seq_present



## parameters
gene = "ND2"
email = "lcoghill@fieldmuseum.org"
seqs_per_taxon = 5
tree_file_name = 'examples/centarchidae.tre'
corrections_file_name = 'examples/corrections.csv'
max_length = 2000




## get leaf names from newick
print "Looking for newick tree file...",
if os.path.isfile(tree_file_name) :
    tree_handle = open(tree_file_name, 'r')
    print "Found. \nParsing newick file %s..." %tree_file_name,
    for line in tree_handle :
        newick = line.strip('\n')+";"
    print "Done."
else :
    print "\nNo tree file found. Cannot continue."
    sys.exit()

ott_ids = {}
taxa_names = []
tree = Tree(newick, format=9)
for leaf in tree :
    name = re.sub(r'_ott.*?$', '', leaf.name)
    name = name.replace("_", " ")
    ott_id = leaf.name.split("_ott")[1]
    ott_ids[name] = ott_id
    taxa_names.append(name.replace("'", ""))


## check for, and process a corrections .csv file here if present
corrections = {}
print "Checking for a corrections file...",
if os.path.isfile(corrections_file_name) :
    corrections_file = open(corrections_file_name, 'r')
    print "Found. \nLoading corrections corrections file...",
    for line in corrections_file :
        rec = line.strip("\n").split(",")
        corrections[rec[0]] = rec[1] + "," + rec[2].replace(" ", "").strip()
    print "Done."
else :
    print "\nNo corrections file found. Proceeding without it."


check_data(gene, taxa_names, email, seqs_per_taxon, corrections, max_length)