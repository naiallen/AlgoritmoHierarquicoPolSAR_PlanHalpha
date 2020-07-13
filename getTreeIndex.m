function index = getTreeIndex( tree,level, ID)
%% look for ID
for ii = 1:(size(tree, 2))
    if (~isempty(tree{level, ii}))
        if (tree{level, ii}.ID == ID)
            index = ii;
            break;
        end
    end
end
end

