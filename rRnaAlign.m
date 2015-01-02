function alignment = rRnaAlign(directory)

%{
  Function rRnaAlign aligns multiple rRNA Sequences
         Input: File Directory containing fasta files
        Output: Multiple Sequence Alignment
%}

% Traverses input directory and stores all the filenames (should be fasta files)
Files=dir(directory);
for k=1:length(Files)
   FileNames=Files(k).name;
end

% Hardcoded first sequence for the purpose of concatenating the Sequence matrix
temp = 'sequences/acinetobacter_sp.fasta';
seqs = fastaread(temp);

%{
  For loop that reads in the fasta files and stores them in seqs array.
   Skips first two integers because they are for current and previous
   directory, 3rd is for 1st file that was previously stored as temp
%}
for i=4:length(Files)
   % The name of the directory inside the current folder
   temp = 'sequences/';
   temp = strcat(temp, Files(i).name);
   seqs = [seqs, fastaread(temp)]; 
end

% Aligns sequences
alignment = multialign(seqs);
%alignment = multialign(seqs, 'terminalGapAdjust', true);
