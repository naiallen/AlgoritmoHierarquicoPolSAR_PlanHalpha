function plotTree(tree)
% Get the leafs IDs
[level, col] = size(tree);
leafs = zeros(col,1);
for ii = 1:col
    leafs(ii) = tree{level, ii}.ID;
end
maxLeaf = max(leafs);
nodes = zeros(1,maxLeaf)*nan;
ent = zeros(1,maxLeaf)*nan;
for lev = 1:level;
    for c =1:lev
        nodes(tree{lev, c}.ID) = tree{lev, c}.parentID;          
        ent(tree{lev, c}.ID) = round(tree{lev, c}.EntropyGain,2);  
    end
end
treeplot(nodes,'ob', '-b')
[x,y] = treelayout(nodes);
for i=1:length(x)
    text( x(i),y(i),...
        char({ strcat('ID: ', num2str(i)), strcat( 'H:', num2str(ent(i)) )}) );
  % {strcat(strcat(strcat('ID: ', num2str(i)), ' H Gain', num2str(ent(i)) )) );
end
grid on
ylabel('level')
end







