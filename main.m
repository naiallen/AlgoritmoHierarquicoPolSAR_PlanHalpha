
clear
close all
clc

folder = 'C:\Users\ncarval1\Documents\PersonalCodes\Matlab\Doutorado\';
L = 3;
q = 3;

%% Load covariance matrix image
im_filename = 'polSAR2.dat_cov.dat';
hdr_filename = 'polSAR2.dat_cov.hdr';
imcov = openPolSARimage(folder, im_filename, hdr_filename);
figure;
 plotPolSARimage( imcov )
[n_row, n_col, n_bands] = size(imcov);


%% Algoritimo hierárquico--------------------------------------------------
tree = [];
parent_id = 1;
im_dummy = ones(n_row, n_col);
numinte = zeros(20,1);
for level = 1:30
    
    next_level = level+1;
    %% Split levels
    if level ==  1
        %% Populate tree
        tree = populate_tree(tree, parent_id, level);
        parent_index = 1;
                
        %% Vectorize image
        data_cov = im2vec( imcov, n_bands);
%         dataPolSAR = im2vec( polSAR_im, 3);
        %% Get parent parameter
        covariance = average_covariance(data_cov);
        par.L = L;
        par.q =  q;
        par.cov = covariance;
        tree{parent_index}.Parameter = par;
        
        %% Get seeds
%         [seeds, data1, data2] = get_seeds_pca(data);
        [seeds, data1, data2] = cov_pca_intrisic_mean(data_cov);
%         [seeds, data1, data2] =  get_seeds_pca_complex(dataPolSAR, data);        
%         [seeds, data1, data2] =  get_seeds_pca_real(dataPolSAR, data);
%         [indexEM, seeds, data1, data2] = get_seed_EM(data_cov,  L, q);
%         [seeds, data1, data2] =  get_seeds_random(data);

        tree{parent_index}.Seeds = seeds;
        tree{parent_index}.SeedsData = {data1 data2};
        
        %% Get entropy parent
        tree{parent_index}.Entropy = whisart_entropy( par, n_bands);
                  
        %% Populate dummy image
        im_dummy = ones(n_row, n_col)*parent_id;
        
        %% Get child ID
        child1_id = tree{parent_index}.ChildsID(1);
        child2_id = tree{parent_index}.ChildsID(2);
        
        child1_index = getTreeIndex( tree, next_level, child1_id);
        child2_index = getTreeIndex( tree, next_level, child2_id);
    else
        
        %% Choose child to open
        n_leaf = size(tree, 2);
        HG = zeros(n_leaf, 1);

        for ii = 1:n_leaf        
            HG(ii) = tree{level, ii}.EntropyGain;
        end
        [~, pos] = max(HG);
        parent_index = pos;
        
        %% Populate tree
        parent_id = tree{level, parent_index}.ID;
        tree = setChildID(tree, parent_id, level);
        tree = populate_tree(tree, parent_id, level);  
        
        child1_id = tree{level, parent_index}.ChildsID(1);
        child2_id = tree{level, parent_index}.ChildsID(2);
        
        child1_index = getTreeIndex( tree, next_level, child1_id);
        child2_index = getTreeIndex( tree, next_level, child2_id);
        %% Get seeds
        seeds = tree{level, parent_index}.Seeds;        
    end
    
    im_aux = im_dummy;
    im_aux(im_aux ~= parent_id) = NaN;
    im_test = ones(n_row, n_col, n_bands);
    for nb = 1:n_bands
        im_test(:,:,nb) = im_aux.*(imcov(:,:,nb));
    end

    %% Run Stochastic Clustering-------------------------------------------
    tic
    [numinte(level), output_im] = stochasticClustering(im_test, 'Bhattacharyya', 2, seeds, 20, L);
    toc
    %----------------------------------------------------------------------   
    
    %% Get Image ----------------------------------------------------------
    im1_cov = zeros(n_row, n_col, n_bands);
    d = output_im(:,:,1);
    d1 = d(:);
    d = output_im(:,:,2);
    d2 = d(:);
    for nb = 1:n_bands
        im1_cov(:,:,nb) = output_im(:,:,1).*(imcov(:,:,nb));
    end
    
    im2_cov = zeros(n_row, n_col, n_bands);
    for nb = 1:n_bands
        im2_cov(:,:,nb) = output_im(:,:,2).*(imcov(:,:,nb));
    end
    %----------------------------------------------------------------------
    
    %% Get vector data ----------------------------------------------------
    data1_cov = im2vec( im1_cov, n_bands);
    data2_cov = im2vec( im2_cov, n_bands);
    %----------------------------------------------------------------------
    
    %% update image dummy
    im_dummy(output_im(:,:,1)== 1) = child1_id;
    im_dummy(output_im(:,:,2)== 1) = child2_id;
    
    %% Get childs indexes
    child_id = tree{level, parent_index}.ChildsID(1);
    child1_index = getTreeIndex( tree, level+1, child_id);
    child_id = tree{level, parent_index}.ChildsID(2);
    child2_index = getTreeIndex( tree, level+1, child_id);
    
    %% Compute the child 1 seeds and parameters----------------------------
    w1_d1 = size(data1_cov,1)/(n_row*n_col);
    w1_d2 = size(data2_cov,1)/(n_row*n_col);
    w2_d1 = size(data1_cov,1)/(size(data1_cov,1)+size(data2_cov,1));
    w2_d2 = size(data2_cov,1)/(size(data1_cov,1)+size(data2_cov,1));
    
    next_level = level+1;
    
    
    if (isempty(data1_cov) || isempty(data2_cov) ||...
       (w2_d1 < 0.01) || (w2_d2 < 0.01))
        %% If there is no division, the branch stops
        tree{next_level, child1_index}.EntropyGain = 0;
        tree{next_level, child2_index}.EntropyGain = 0;
    else
        %% Update child node
        if (w1_d1 < 0.005)
            %% If there is no division, the branch stops
            tree{next_level, child1_index}.EntropyGain = 0;
        else
            tree = updateChild(tree, next_level, child1_index, parent_index, data1_cov, L, q, n_row, n_col);
        end
        if (w1_d2 < 0.005)
            tree{next_level, child2_index}.EntropyGain = 0;
        else
            tree = updateChild(tree, next_level, child2_index, parent_index, data2_cov, L, q,  n_row, n_col);
        end
    end
    
    HG = 0;
    for index =1:(size(tree, 2))
        if (~isempty(tree{next_level, index}))
            HG = HG + tree{next_level, index}.EntropyGain;
        end
    end
    if (HG <= 0)
        break;
    end   
end
imwrite(im2uint8(mat2gray(im_dummy)), strcat(folder, 'im_classifgray.png'))
