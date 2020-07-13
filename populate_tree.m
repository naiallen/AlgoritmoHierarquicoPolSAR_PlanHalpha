%%=========================================================================
% Brief: Populate hirarchical tree
%   Input: 
%       tree: hierarchical clustering tree
%       parent_id:  Parebt cluster ID
%       level: tree level
%   Output
%       tree: hierarchical clustering tree
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================
function tree = populate_tree(tree, parent_id, level)

if isempty(tree)
    tree{level} = initTreeParameter();
    %Update actual level---------------------------------------------------
    tree{level}.ID = parent_id;
    tree{level}.parentID = 0;
    tree{level}.Color = getcolors();
    tree{level}.Parameter = nan; 
    tree{level}.ChildsID = [parent_id+1 parent_id+2];
    tree{level}.Seeds = nan; 
    tree{level}.Entropy = nan; 
    tree{level}.EntropyGain = nan; 
    %----------------------------------------------------------------------
    
    next_level = level+1;
    
    %Create the childs----------------------------------------------------- 
    tree{next_level, 1} = initTreeParameter();
    tree{next_level, 1}.ID = tree{level}.ChildsID(1);
    tree{next_level, 1}.parentID = tree{level}.ID ;
    tree{next_level, 1}.Color = getcolors();
    tree{next_level, 1}.Parameter = nan;
    tree{next_level, 1}.ChildsID = [nan nan];
    tree{next_level, 1}.Seeds = nan; 
    tree{next_level, 1}.Entropy = nan; 
    tree{next_level, 1}.EntropyGain = nan; 
    
    tree{next_level, 2} = initTreeParameter();
    tree{next_level, 2}.ID = tree{level}.ChildsID(2);
    tree{next_level, 2}.parentID = tree{level}.ID ;
    tree{next_level, 2}.Color = getcolors();
    tree{next_level, 2}.Parameter = nan;
    tree{next_level, 2}.ChildsID = [nan nan];    
    tree{next_level, 2}.Seeds = nan; 
    tree{next_level, 2}.Entropy = nan; 
    tree{next_level, 2}.EntropyGain = nan; 
    %----------------------------------------------------------------------
else
    %% look for ID
    parent_index = getTreeIndex( tree,level, parent_id);
   
    %% look for max ID
    ID = 0;
    for index =1:(size(tree, 2))
        if (~isempty(tree{level, index}))
            maxID = tree{level, index}.ID;
            if (maxID > ID)
                ID = maxID;
            end
        end
    end
 
    %Create the childs
    next_level = level+1;
    tree{next_level, 1} = initTreeParameter();
    tree{next_level, 1}.ID = tree{level, parent_index}.ChildsID(1);
    tree{next_level, 1}.parentID = tree{level, parent_index}.ID;
    tree{next_level, 1}.Color = getcolors();
    tree{next_level, 1}.Parameter = nan;
    tree{next_level, 1}.ChildsID = [nan nan];
    tree{next_level, 1}.Seeds = nan;
    tree{next_level, 1}.Entropy = nan;
    tree{next_level, 1}.EntropyGain = nan;
    
    tree{next_level, 2} = initTreeParameter();
    tree{next_level, 2}.ID = tree{level, parent_index}.ChildsID(2);
    tree{next_level, 2}.parentID =  tree{level, parent_index}.ID ;
    tree{next_level, 2}.Color = getcolors();
    tree{next_level, 2}.Parameter = nan;
    tree{next_level, 2}.ChildsID = [nan nan];
    tree{next_level, 2}.Seeds = nan;
    tree{next_level, 2}.Entropy = nan;
    tree{next_level, 2}.EntropyGain = nan;
    
    
    %Popula os nos subsequentes
    index = 3;
    for ii =1:(size(tree, 2))
        if (~isempty(tree{level, ii}))
            if (ii ~= parent_index)
                tree{next_level, index} = tree{level, ii};
                index=index+1;
            end
        end
    end
end
end