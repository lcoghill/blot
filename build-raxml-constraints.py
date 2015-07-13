from ete2 import Tree



## params
tree_file = 'examples/centarchidae.tre'

## read the original tree into a ete2 tree object
handle = open(tree_file, 'r')
for line in handle :
    tree_string = line.strip("\n")
tree_string = tree_string.replace("''", "") + ";" # used to parse out escaped quotes present in some ot trees
tree = Tree(tree_string, format=1)


## collect nodes 1 step above leaves
leaf_parents = []
for leaf in tree :
    node = leaf.up
    if node not in leaf_parents and node.up not in leaf_parents :
        leaf_parents.append(node)

## for each parent node, get all the leaves under it, and build a clade list
clades = []
single_taxa = []
for node in leaf_parents :
    clade = []
    clade = node.get_leaf_names()
    if len(clade) > 1 :
        clades.append(clade)
    else :
        clades.append(clade[0])

## remove any single taxon clades if that taxa is present in another clade
for clade in clades :
    clades = [x for x in single_taxa if x not in clade]

## build the newick string
newick = "("
for s in single_taxa :
    newick = newick + s + ","
newick = newick[:-1]
for group in groups :
    constraint_clade = ",".join(group)
    newick = newick + ",(" + constraint_clade + ")"
newick  = newick + ");"

## dump newick to a file
handle = open('constraint_tree.tre', 'w')
handle.write(newick)
handle.close()