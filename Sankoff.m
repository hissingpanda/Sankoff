function out = Sankoff(node_l, node_r)

%{
    Function Sankoff that performs Sankoff's algorithm on two phylogenetic
     nodes
         Input: Left and right node of a branch
        Output: The parsimony values for each nucleotide at the node 
%}

% Default scoring matrix
    %	A	U	G	C	-
scor = [0,	3,	4,	9,	8 ;  %A
        3,	0,	2,	4,	8 ;  %U
        4,	2,	0,	4,	8 ;  %G
        9,	4,	4,	0,	8 ;  %C
        8,	8,	8,	8,	8 ;]; %- (gap)

% Initializes output array
out = [inf, inf, inf, inf, inf];

%{
    Loops through 5 possible nucleotides at each node (includes dash character)
     keeping track of the minimum parsimony at each position by using the 
     scoring matrix and current position.

     Tries all 25 possibilities (for the 5 possible nucleotides).
%}

len = length(node_l);
for i=1:len
    
    left = inf; 
    for j=1:len      
        cur = scor(i,j) + node_l(j);
        left = min(left, cur);
    end
    
    right = inf;
    for j=1:len
        cur = scor(i,j) + node_r(j);
        right = min(right, cur);
    end
    out(i) = left + right; 
    
end
%display(out);
