Created by: 
Paul Lanctot
lanct020

Matlab files in this directory:
rRnaAlign.m
rRnaDist.m
Sankoff.m
treeParse.m
toys.m

Extra files in this directory
multialign_results.txt
pairwise_distances.txt
parsimony_results.txt
tree.jpg
tree_nodes.jpg
Toy1.pdf
Toy2.pdf
Toy3.pdf
README.txt

Matlab files
-------------
Note: for Problems 1 & 2, the fasta sequences are assumed to be in a directory named 'sequences'.
-----------

-------------------------
rRnaAlign.m (Problem 1)
-----------------------
Function rRnaAlign aligns multiple rRNA Sequences
     Input: File Directory containing fasta files
    Output: Multiple Sequence Alignment

You can call the function in Matlab by:
rRnaAlign('path_to_directory')

-------------------------
rRnaDist.m (Problem 2)
-----------------------
Function rRnaDist computes the pairwise distances of sequences 
 with function seqpdist and applies seqneighjoin to infer a phenological tree

     Input: File Directory containing fasta files
    Output: Pairwise distances of sequences, phylogenetic tree

You can call the function in Matlab by:
rRnaDist('path_to_directory')

-----------------------
Sankoff.m (Problem 3)
----------------------
Function Sankoff that performs Sankoff's algorithm on two phylogenetic nodes
         Input: Left and right node of a branch
        Output: The parsimony values for each nucleotide at the node 

You can call the function in Matlab by:
Sankoff(left_node_of_tree_branch, right_node_of_tree_branch)


--------------------------
treeParse.m (Problem 4)
-------------------------
Function treeParse that performs Sankoff's algorithm on a phylogenetic tree.
     Input: phylogenetic tree, multiple sequence alignment
    Output: the parsimony score of the input tree, the sequence alignment of each internal node/branch

You can call the function in Matlab by:
treeParse(phylogenetic_tree)

--------------------
toys.m (Problem 3)
-------------------
Initialization and testing of three toy examples (taken from the slides, pgs 82, 90, 83, respectively), for use with Sankoff.m.
 Generates figure with parsimony score at each node.

Not meant to be called in Matlab. Instead highlight each individual toy example, right-click and select 'evalutate selection' (or press f9).


---------
Results
--------
multialign_results.txt: The output result of rRnaAlign - a multiple sequence alignment.
pairwise_distances.txt: An output result of rRnaDist - the pairwise distances of all the input sequences.
parsimony_results.txt:  An output result of treeParse - the parsimony score of the input tree and the sequence alignment of each internal node/branch.
tree.jpg:        An output result of rRnaDist - the phylogenetic tree of all the input sequences (with branch labels).
tree_nodes.jpg:  An output result of rRnaDist - the phylogenetic tree of all the input sequences (with branch and leaf node labels).
Toy1.pdf: pdf file of the figure generated in Matlab by running Toy example 1, from pg 82, with Sankoff.m
Toy2.pdf: pdf file of the figure generated in Matlab by running Toy example 2, from pg 90, with Sankoff.m
Toy3.pdf: pdf file of the figure generated in Matlab by running Toy example 3, from pg 83, with Sankoff.m
