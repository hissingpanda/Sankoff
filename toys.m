%{
    Initialization and testing of three toy examples (taken from the slides,
     pgs 82, 90, 83, respectively).

    Generates figure with parsimony score at each node.
%}

%%%%%%%%%%%%%%%%%%%%%%
%Toy1 example pg 82
%%%%%%%%%%%%%%%

% Initializes individual nucleotides
a = [0, inf, inf, inf, inf];
u = [inf, 0, inf, inf, inf];
g = [inf, inf, 0, inf, inf];
c = [inf, inf, inf, 0, inf];
%dash = [inf, inf, inf, inf, 0];

% Stores nodes in leafs and branches in proper order
%  for ease of use with final output figure
Leaf{1} = a;
Leaf{3} = u;
Leaf{4} = g;
Leaf{2} = c;

% Computes Sankoff's algorithm nodes up to the root (Branch1)
Branch{2} = Sankoff(a, c);
Branch{3} = Sankoff(u, g);
Branch{1} = Sankoff(Branch{2}, Branch{3});

display(Branch{2});
display(Branch{3});
display(Branch{1});

% Creates treeplot with parsimony values at each node
%  for better visualization.
nodes = [0 1 1 2 2 3 3];
treeplot(nodes, '');
count = size(nodes,2);
[x,y] = treelayout(nodes);
x = x';
y = y';

% Loop that generates text boxes to be displayed at each node of the
% treeplot with the appropriate parsimony values
for i=1:count
    if (i <= 3)
        values = cellstr(mat2str(Branch{i}));
        text(x(i,1), y(i,1), values, 'VerticalAlignment','bottom','HorizontalAlignment','right')
    else
        values = cellstr(mat2str(Leaf{i-3}));
        text(x(i,1), y(i,1), values, 'VerticalAlignment','bottom','HorizontalAlignment','right')
    end
end

title({'Toy1 Phylogenic Tree'},'FontSize',12,'FontName','Times New Roman');



%%%%%%%%%%%%%%%%%%%%%
%Toy2 example pg 90
%%%%%%%%%%%%%%

% Initializes individual nucleotides
a = [0, inf, inf, inf, inf];
u = [inf, 0, inf, inf, inf];
g = [inf, inf, 0, inf, inf];
c = [inf, inf, inf, 0, inf];
%dash = [inf, inf, inf, inf, 0];

% Stores nodes in leafs and branches in proper order
%  for ease of use with final output figure
Leaf{1} = a;
Leaf{2} = c;
Leaf{3} = g;
Leaf{4} = g;

% Computes Sankoff's algorithm nodes up to the root (Branch1)
Branch{2} = Sankoff(a, c);
Branch{3} = Sankoff(g, g);
Branch{1} = Sankoff(Branch{2}, Branch{3});

display(Branch{2});
display(Branch{3});
display(Branch{1});

% Creates treeplot with parsimony values at each node
%  for better visualization.
nodes = [0 1 1 2 2 3 3];
treeplot(nodes, '');
count = size(nodes,2);
[x,y] = treelayout(nodes);
x = x';
y = y';

% Loop that generates text boxes to be displayed at each node of the
% treeplot with the appropriate parsimony values
for i=1:count
    if (i <= 3)
        values = cellstr(mat2str(Branch{i}));
        text(x(i,1), y(i,1), values, 'VerticalAlignment','bottom','HorizontalAlignment','right')
    else
        values = cellstr(mat2str(Leaf{i-3}));
        text(x(i,1), y(i,1), values, 'VerticalAlignment','bottom','HorizontalAlignment','right')
    end
end

title({'Toy2 Phylogenic Tree'},'FontSize',12,'FontName','Times New Roman');




%%%%%%%%%%%%%%%%%%%%%
%Toy3 example pg 83
%%%%%%%%%%%%%%


a = [0, inf, inf, inf, inf];
u = [inf, 0, inf, inf, inf];
g = [inf, inf, 0, inf, inf];
c = [inf, inf, inf, 0, inf];
%dash = [inf, inf, inf, inf, 0];


% Stores nodes in leafs and branches in proper order
%  for ease of use with final output figure
Leaf{1} = u;
Leaf{2} = c;
Leaf{3} = g;

% Computes Sankoff's algorithm nodes up to the root (Branch1)
Branch{2} = Sankoff(c, g);
Branch{1} = Sankoff(Branch{2}, u);
display(Branch{2});
display(Branch{1});

% Creates treeplot with parsimony values at each node
%  for better visualization.
nodes = [0 1 1 2 2];
treeplot(nodes, '');
count = size(nodes,2);
[x,y] = treelayout(nodes);
x = x';
y = y';

% Loop that generates text boxes to be displayed at each node of the
% treeplot with the appropriate parsimony values
for i=1:count
    if (i <= 2)
        values = cellstr(mat2str(Branch{i}));
        text(x(i,1), y(i,1), values, 'VerticalAlignment','bottom','HorizontalAlignment','right')
    else
        values = cellstr(mat2str(Leaf{i-2}));
        text(x(i,1), y(i,1), values, 'VerticalAlignment','bottom','HorizontalAlignment','right')
    end
end

title({'Toy3 Phylogenic Tree'},'FontSize',12,'FontName','Times New Roman');






%{
data = {'u' 'u';
        'a' 'a';
        'g' 'g';
        'c' 'c';};

for ind = 1:4
    seqs(ind).Header   = data{ind,1};
   % seqs(ind).Sequence = s( round(rand(1,sLength)*numRands) );
end
seqs(1).Sequence = 'AG'
seqs(2).Sequence = 'CG'

seqs(3).Sequence = 'AU'
seqs(4).Sequence = 'CC'


%alignment = multialign(seqs);
distances = seqpdist(seqs,'Method','Jukes-Cantor','Alpha','DNA');
tree = seqlinkage(distances,'UPGMA',seqs);
names = get(tree,'LeafNames')
pointers=get(tree,'Pointers')


h = plot(tree,'orient','top');
ylabel('Evolutionary distance');
text( h(1,1), h(1,2), 'Label first node');
%}
%set(h.terminalNodeLabels,'Rotation',65)
%annotation('textbox', [0.2,0.4,0.1,0.1],...
%           'String', out);
%result = treeParse(tree);
