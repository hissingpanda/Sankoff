function [parsimony, seqs] = treeParse(tree, alignment)

%{
    Function treeParse that performs Sankoff's algorithm on a phylogenetic
    tree.
         Input: phylogenetic tree, multiple sequence alignment
        Output: parsimony score for the input tree
%}

% Initializes individual nucleotides
a = [0, inf, inf, inf, inf];
u = [inf, 0, inf, inf, inf];
g = [inf, inf, 0, inf, inf];
c = [inf, inf, inf, 0, inf];
dash = [inf, inf, inf, inf, 0];
parsimony = 0;

names    = get(tree,'LeafNames');
pointers = get(tree,'Pointers');
numLeaves   = get(tree,'NumLeaves');
branches   = get(tree,'NumBranches');
branchNames = get (tree, 'BranchNames');

len = branches;

% Gets the sequences at each leaf
seqlen = length(alignment);
for i=1:seqlen
    aligned_seqs{i} = alignment(i).Sequence;
end
% display(aligned_seqs{1}(1));

%{
    Loops through entire sequence length, translating the lowest parsimony
    value to a nucleotide to be used as each internal nodes sequence
    alignment
%}
seqlen = length(aligned_seqs{1});
for z=1:seqlen

    for i=1:len
        % Pointers are the branch/leaf connectivity array
        left = pointers(i,1);
        right = pointers(i,2);
        
        % At a leaf node
        if (left <= numLeaves)
            if (right <= numLeaves)
                leafSeq(1) = aligned_seqs{left}(z);
                leafSeq(2) = aligned_seqs{right}(z);
                
                %bool, true if at leaf, false if at internal node
                atLeaf(1) = 1;
                atLeaf(2) = 1;
        % Second node is at ancestor node
            else
                leafSeq(1) = aligned_seqs{left}(z);
                internalSeq{2} = scores((right - numLeaves),:);
                atLeaf(1) = 1;
                atLeaf(2) = 0;
            end
        % First node is ancestor node
        elseif (left > numLeaves)
            if (right <= numLeaves)
                internalSeq{1} = scores((left - numLeaves),:);
                leafSeq(2) = aligned_seqs{right}(z);
                atLeaf(1) = 0;
                atLeaf(2) = 1;
            % Both nodes are ancestors
            else
                internalSeq{1} = scores((left - numLeaves),:);
                internalSeq{2} = scores((right - numLeaves),:);
                atLeaf(1) = 0;
                atLeaf(2) = 0;
            end
        end
        
        for x=1:2
            %converts to ascii value for easy switch case 65 = 'A' for leaf
            %sequences, this is to give the default parsimony score for
            %each nucleoted A = [0, inf, inf, inf, inf]
            if (atLeaf(x))
                switch uint8(leafSeq(x))
                    % A, U, G, C, - in order
                    case 65
                        node{x} = a;
                    case 85
                        node{x} = u;
                    case 71
                        node{x} = g;
                    case 67
                        node{x} = c;
                    otherwise
                        node{x} = dash;
                end            
            end
        end
        
        % Performs Sankoff's algorithm on each node by determining if it is
        % an internal node or leaf node
        if (atLeaf(1) && atLeaf(2))
            scores(i,:) = Sankoff(node{1}, node{2});
        elseif (atLeaf(1) && atLeaf(2) == 0)
            scores(i,:) = Sankoff(node{1}, internalSeq{2});
        elseif (atLeaf(1) == 0 && atLeaf(2))
            scores(i,:) = Sankoff(internalSeq{1}, node{2});
        else
            scores(i,:) = Sankoff(internalSeq{1}, internalSeq{2});
        end
    end
    
    % Gets the minimum parsimony for each branch at the sequence position
    minimum_val = min(scores,[],2);
    % 37 = the root node, what we will use for the parsimony score
    tmp_parsimony = minimum_val(37);
    parsimony = parsimony + tmp_parsimony;
    %display(scores);
    
    % Switch case used to generate a nucleotide sequence at each position
    % for the internal nodes
    for i=1:len
        for k=1:5
            if (scores(i,k) == minimum_val(i))
                switch k
                    % A, U, G, C, -
                    case 1
                        seqs(i,z) = 'A';
                    case 2
                        seqs(i,z) = 'U';
                    case 3
                        seqs(i,z) = 'G';
                    case 4
                        seqs(i,z) = 'C';
                    otherwise
                        seqs(i,z) = '-';
                end            
            end
        end
    end
    %display(seqs);
    %display(branchNames);
end
display(seqs);
%display(minimum_val);
