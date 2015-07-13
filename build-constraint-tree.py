from ete2 import Tree



## params
tree_file = 'examples/centarchidae.tre'

## read the original tree into a ete2 tree object
handle = open(tree_file, 'r')
for line in handle :
    tree_string = line.strip("\n")
tree_string = tree_string.replace("''", "") + ";" # used to parse out escaped quotes present in some ot trees
tree = Tree(tree_string, format=9)

for leaf in tree :
    print leaf.get_sisters()