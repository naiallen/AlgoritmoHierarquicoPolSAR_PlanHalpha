function tree = setChildID(tree, parent_id, level)
 %look for max ID
    ID = 0;
    for index =1:(size(tree, 2))
        if (~isempty(tree{level, index}))
            maxID = tree{level, index}.ID;
            if (maxID > ID)
                ID = maxID;
            end
        end
    end
    %% look for ID
    parent_index = getTreeIndex( tree,level, parent_id);
    
    %% Update Tree
    tree{level, parent_index}.ChildsID(1) = ID+1;    
    tree{level, parent_index}.ChildsID(2) = ID+2;
end