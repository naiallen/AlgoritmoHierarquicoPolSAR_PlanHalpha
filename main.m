
clear
close all
clc

folder = 'C:\Users\naial\Documents\Naiallen\Doutorado\Data\Oberpffafenhofe_Gilching\';
L = 3;
q = 3;

%% Load coherence matrix image
im_filename = 'ESAR97VV_small_coh.dat';
hdr_filename = 'ESAR97VV_small_coh.hdr';
im_coh = openPolSARimage(folder, im_filename, hdr_filename);
imcoh = im_coh;
figure;
plotPolSARimage( im_coh );

% imcoh2 = fmedia(im_coh,L);
% imcoh = filter_imean (im_coh, L);
% figure;
% subplot(121)
% plotPolSARimage( imcoh2 );
% title('Box Filter')
% subplot(122)
% plotPolSARimage( imcoh );
% title('Riemann Box Filter')

%% Load covariance matrix image
im_filename = 'ESAR97VV_small_cov.dat';
hdr_filename = 'ESAR97VV_small_cov.hdr';
im_cov = openPolSARimage(folder, im_filename, hdr_filename);
imcov = im_cov;
figure;
plotPolSARimage( imcov );

% imcov2 = fmedia(im_cov,L);
% imcov = filter_imean (im_cov, L);
% figure;
% subplot(121)
% plotPolSARimage( imcov2 );
% title('Box Filter')
% subplot(122)
% plotPolSARimage( imcov );
% title('Riemann Box Filter')

%% Compute spacial number of looks-----------------------------------------
% tic
% 
% % imcov = fmedia(polSAR_im, L);
% toc

% % Compute Covariance------------------------------------------------------
% tic
% imcov = getPolSARcovariance(nlook_polSAR_im);
[n_row, n_col, n_bands] = size(imcov);
% figure;
% plotPolSARimage( imcov )
% toc

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
        data_coh = im2vec( imcoh, n_bands);
        [entropia, alpha] = planHalpha( data_coh);
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
      
    %% Plot parent image---------------------------------------------------
    h1 = figure('Position', [10 10 1400 700]);
    subplot(242)
    plotPolSARimage(im_test)
    t1 = '';
    title({strcat(strcat(strcat('Parent - Level: ', num2str(level)), ' - ID:'), num2str(tree{level, parent_index}.ID)), t1}, 'color', tree{level, parent_index}.Color);
%     data1 = [entropia(~isnan(tree{level, parent_index}.SeedsData{1})) alpha(~isnan(tree{level, parent_index}.SeedsData{1}))];
%     data2 = [entropia(~isnan(tree{level, parent_index}.SeedsData{2})) alpha(~isnan(tree{level, parent_index}.SeedsData{2}))];
%     subplot(2, 3, 2)
%     plot( data1(:,1), data1(:,2), '.', 'color', tree{next_level, child1_index}.Color );
%     hold on
%     plot( data2(:,1), data2(:,2), '.', 'color', tree{next_level, child2_index}.Color );
%     hold on
%     plotPlanHalpha
%     
%     plot3(sqrt(data1(:, 1)), sqrt(data1(:, 5)), sqrt(data1(:, 9)), '.', 'color', tree{next_level, child1_index}.Color);
%     hold on
%     plot3(sqrt(data2(:, 1)), sqrt(data2(:, 5)), sqrt(data2(:, 9)), '.', 'color', tree{next_level, child2_index}.Color);
%     grid on;
    %----------------------------------------------------------------------

    
    %% Run Stochastic Clustering-------------------------------------------
    tic
    [numinte(level), output_im] = stochasticClustering(im_test, 'Bhattacharyya', 2, seeds, 20, L);
%     numinte(level) = numinte(level)+indexEM;
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
    
    subplot(246)
    plotPolSARimage(im1_cov);
    t1 = '';%sprintf('Parameter %5.2f:' , im_tree{next_level, 1}.Parameter);
    title({strcat(strcat(strcat('Child 1 - Level: ', num2str(next_level)), ' - ID:'), num2str(tree{next_level, child1_index}.ID)), t1}, 'color', tree{next_level, child1_index}.Color);

    subplot(247)
    plotPolSARimage( im2_cov )
    t1 = '';%sprintf('Parameter %5.2f:' , im_tree{next_level, 2}.Parameter);
    title({strcat(strcat(strcat('Child 2 - Level: ', num2str(next_level)), ' - ID:'), num2str(tree{next_level, child2_index}.ID)), t1}, 'color', tree{next_level, child2_index}.Color);
    
    
    data1 = [entropia(~isnan(d1)) alpha(~isnan(d1))];
    data2 = [entropia(~isnan(d2)) alpha(~isnan(d2))];
    data = [data1; data2];
    
    subplot(243)
%     plot( data(:,1), data(:,2), '.', 'color', tree{level, parent_index}.Color );
%     hold on
    plotPlanHalpha( data(:,1), data(:,2));
    
    ent1 = 0;
    ent2 = 0;
    alp1 = 0;
    alp2 = 0;
    if (isempty(data1) || isempty(data2))   
        error_entropy = 0;
        error_alpha = 0;
    else
        ent1 = mean(data1(:,1));
        alp1 = mean(data1(:,2));
        
        ent2 = mean(data2(:,1));
        alp2 = mean(data2(:,2));
        
        error_entropy = abs(ent1 - ent2);
        error_alpha = abs(alp1 - alp2);
    end
    
    subplot(244)
    plotHalphaOnly;
    hold on
    plot(  ent1,  alp1, 'o', 'LineWidth',2, 'color', tree{next_level, child1_index}.Color );    
    hold on
    plot(  ent2,  alp2, 'o', 'LineWidth',2, 'color', tree{next_level, child2_index}.Color );    
    plot(  mean(data(:,1)),  mean(data(:,2)), 'o', 'LineWidth',2, 'color', tree{level, parent_index}.Color );    
    title({strcat('Entropy Error: ', num2str(error_entropy)),strcat('Alpha Angle Error: ', num2str(error_alpha))})
    
    subplot(245)
%     plot( data1(:,1), data1(:,2), '.', 'color', tree{next_level, child1_index}.Color );
%     hold on
    plotPlanHalpha( data1(:,1), data1(:,2) )
%     temp = [mean(data1); mean(data2)];
%     temp(:,2) = temp(:,2)/90;
%     error = sqrt(mean(abs(temp(1,:).^2 - temp(2,:).^2)));    
%     title(strcat('Error: ', num2str(error)));
    
    subplot(248)
%     plot( data2(:,1), data2(:,2), '.', 'color', tree{next_level, child2_index}.Color );
%     hold on
    plotPlanHalpha(data2(:,1), data2(:,2))
%     temp = [mean(data1); mean(data2)];
%     temp(:,2) = temp(:,2)/90;
%     error = sqrt(mean(abs(temp(1,:).^2 - temp(2,:).^2)));    
%     title(strcat('Error: ', num2str(error)));
    
    
%     %% PolSAR image
%     imPolSAR1 = zeros(n_row, n_col, n_bands);
%     for nb = 1:3
%         imPolSAR1(:,:,nb) = output_im(:,:,1).*(polSAR_im(:,:,nb));
%     end
%     
%     imPolSAR2 = zeros(n_row, n_col, n_bands);
%     for nb = 1:3
%         imPolSAR2(:,:,nb) = output_im(:,:,2).*(polSAR_im(:,:,nb));
%     end
%     
%     dataPolSAR1 = im2vec( imPolSAR1, 3);
%     dataPolSAR2 = im2vec( imPolSAR2, 3);
    
    
%     
%     plot3(sqrt(data1(:,1)), sqrt(data1(:, 5)), sqrt(data1(:, 9)), '.', 'color', tree{next_level, child1_index}.Color);
%     hold on
%     plot3(sqrt(data2(:, 1)), sqrt(data2(:,5)), sqrt(data2(:, 9)), '.', 'color', tree{next_level, child2_index}.Color);
%     grid on;
    
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
            tree = updateChild(tree, next_level, child1_index, data1_cov, L, q, n_row, n_col);
        end
        if (w1_d2 < 0.005)
            tree{next_level, child2_index}.EntropyGain = 0;
        else
            tree = updateChild(tree, next_level, child2_index, data2_cov, L, q,  n_row, n_col);
        end
    end
    
    saveas(h1, strcat(folder, strcat(strcat('level_', num2str(level)), '.png')));
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
h4  = figure('Position', [10 10 2000 700]);
plotTree(tree);
saveas(h4, strcat(folder, 'tree.png'));
% level = level+1;
% n = size( tree, 2);
% im_classif = zeros(n_row, n_col, 3);
% for ii = 1:n_row
%     for jj=1:n_col    
%         ID = im_dummy(ii, jj);
%         index = getTreeIndex( tree,level, ID);
%         im_classif(ii, jj, :) = tree{level, index}.Color;      
%     end
% end

% h2  = figure();
% imagesc (im_classif);
% camroll(-90);
% set(gca, 'xdir', 'reverse');
% axis equal
% saveas(h2,strcat(folder, 'im_classif.png'));
% % imwrite(im_classif, strcat(folder, 'im_classif.png'))


h2  = figure();
imagesc (im_dummy);
colormap gray
camroll(-90);
set(gca, 'xdir', 'reverse');
axis equal
imwrite(im2uint8(mat2gray(im_dummy)), strcat(folder, 'im_classifgray.png'))


h3 = figure;
plot(numinte)
xlabel('Level')
ylabel('Number of iteration');
grid on
saveas(h3, strcat(folder, 'numberofiteration.png'));