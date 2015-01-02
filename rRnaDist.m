function [alignment_distances, tree] = rRnaDist(directory)

%{
   Function rRnaDist computes the pairwise distances of sequences 
    with function seqpdist and applies seqneighjoin to infer a phenological tree

         Input: File Directory containing fasta files
        Output: Pairwise distances of sequences, phylogenetic tree
%}

% Traverses input directory and stores all the filenames (should be fasta files)
Files=dir(directory);
for k=1:length(Files)
   FileNames=Files(k).name;
end

% For loop that reads in the fasta files and stores them in seqs array.
%  Skips first two because they are for current and previous directory
for i=3:length(Files)
   temp = 'sequences/';
   temp = strcat(temp, Files(i).name);
   [head{i-2}, seqs{i-2}] = fastaread(temp); 
end

%computes pairwise distances
alignment_distances = seqpdist(seqs, 'Alphabet', 'NT');

%{
% Alternative
alignment_distances = seqpdist(seqs,'Method','jukes-cantor',...
                'Alphabet','NT',...
                'Indels','score',...
                'Scoringmatrix','pam250',...
                'PairwiseAlignment',true);
%}

% Generates and view phylogenetic tree
tree = seqneighjoin(alignment_distances, 'equivar', head);
view(tree);
