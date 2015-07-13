import re
import os
import sys
import getopt
import subprocess
from ete2 import Tree
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo
from Bio.Phylo.Consensus import *
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO as SIO



def get_seqs(t, gene, email, seqs_per_taxon, corrections, max_length, names_and_numbers) :

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
    
    if len(record['IdList']) < seqs_per_taxon :
       print "Warning: there are fewer than %s %s sequences available for %s. Using all %s available." %(seqs_per_taxon, gene, t, len(record['IdList']))

    ## get sequence data
    seq_recs = []
    for gi in record['IdList'] :
        handle = Entrez.efetch(db="nucleotide", id=gi ,rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        if len(seq_recs) <= seqs_per_taxon and len(record.seq) <= max_length :
            seq_recs.append(record)
        # TODO: filter sequences for quality in some way

    # build file object to pass to muscle for alignment
    output = open("fasta/"+t.replace(" ", "_")+"_"+gene+"_seqs.fas", 'w')
    SeqIO.write(seq_recs, output, "fasta")

    return len(seq_recs)



## parameters
gene = "ND2"
email = "test@example.com"
max_length = 2000
seqs_per_taxon = 3
tree_file_name = 'examples/centarchidae.tre'
corrections_file_name = 'examples/corrections.csv'
max_muscle_iters = "3"
max_muscle_trees = "2"




## get leaf names from newick
os.mkdir("fasta")
tree_handle = open(tree_file_name, 'r')
print "Parsing newick file %s..." %tree_file_name
for line in tree_handle :
    newick = line.strip('\n')+";"

ott_ids = {}
names_and_numbers = {}
taxa_names = []
tree = Tree(newick, format=8)
for leaf in tree :
    name = re.sub(r'_ott.*?$', '', leaf.name)
    name = name.replace("_", " ")
    ott_id = leaf.name.split("_ott")[1]
    ott_ids[name] = ott_id
    taxa_names.append(name.replace("'", ""))
    
tree.write(format=9, outfile="examples/constraint_tree.tre")
## check for, and process a corrections .csv file here if present
corrections_file = open(corrections_file_name, 'r')
corrections = {}
print "Loading corrections corrections file..."
for line in corrections_file :
    rec = line.strip("\n").split(",")
    corrections[rec[0]] = rec[1] + "," + rec[2].replace(" ", "").strip()


## get sequence data for X number of records for gene X
## and build consensus sequence of those sequences for each taxa
leaf_con_seqs = {}
con_seqs = []
for t in taxa_names :
    # get sequences from genbank for all taxa, build fasta style file object and return
    num_res = get_seqs(t, gene, email, seqs_per_taxon, corrections, max_length, names_and_numbers)
    ## build alignment for those sequences
    if num_res > 1 :
        seqs_f = "fasta/"+t.replace(" ", "_") + "_" + gene + "_seqs.fas"
        align_f = "fasta/"+t.replace(" ", "_") + "_" + gene + "_alignment.fas"
        subprocess.call(["muscle", "-in", seqs_f, "-out", align_f, "-maxiters", max_muscle_iters, "-quiet"])
        msa = AlignIO.read(align_f, "fasta")
        msa_summary = AlignInfo.SummaryInfo(msa)
        consensus = msa_summary.dumb_consensus()
        # build seqrec object here and add to con_seqs
        record = SeqRecord(Seq(str(consensus), IUPAC.ambiguous_dna), id=t.replace(" ", "_")+"_ott"+ott_ids[t], description="")
        con_seqs.append(record)

    ## get all the singles
    else :
        fasta_f = "fasta/"+t.replace(" ", "_") + "_" + gene + "_seqs.fas"
        handle = open(fasta_f, "rU")
        for record in SeqIO.parse(handle, "fasta") :
            record.id = t.replace(" ", "_")+"_ott"+ott_ids[t]
            record.description = ""
            con_seqs.append(record)

## write all species seqs to fasta
out_fasta_file = tree_file_name.split(".")[0] + "_" + gene + "_sequences.fas"
out_fasta = open(out_fasta_file, 'w')
print "Writing %s raw sequences fasta file to %s..." %(len(con_seqs), out_fasta_file)
SeqIO.write(con_seqs, out_fasta,  "fasta")
out_fasta.close()
 
## build alignment of new fasta
print "Building alignment of %s..." %out_fasta_file
final_msa_file = tree_file_name.split("/")[1].replace(".tre", "") + "_" + gene + "_alignment.fas"
subprocess.call(["muscle", "-in", out_fasta_file, "-diags", "-maxiters", max_muscle_iters,
                "-maxtrees", max_muscle_trees, "-out", final_msa_file]) 

## convert to relaxed phylip here
temp_alignment = AlignIO.read(final_msa_file, "fasta")
AlignIO.write(temp_alignment, final_msa_file.replace(".fas", ".phy"), "phylip-relaxed")

## build a ML tree with new alignment
print "Building final tree with constrained topology to estimate branch lengths..."
subprocess.call(["phyml", "-i", final_msa_file.replace(".fas", ".phy"), "-u", "examples/constraint_tree.tre", "-d" "nt", "-m", "GTR", "-o", "l"])


## write code to replace short names with original names in tree and alignment files
print "Complete."
