function tree = updateChild(tree, next_level, child_index, data_cov, L, q,  n_row, n_col)

% [ch_seeds, data1, data2] = get_seeds_pca_complex(dataPolSAR, data);
% [ch_seeds, data1, data2] = get_seeds_pca_real(dataPolSAR, data_cov);
% [~, ch_seeds, data1, data2] = get_seed_EM(data_cov, L, q);
% [ch_seeds, data1, data2] = get_seeds_pca(data);
[ch_seeds, data1, data2] = cov_pca_intrisic_mean(data_cov);
% [ch_seeds, data1, data2] = get_seeds_random(data);

tree{next_level, child_index}.Seeds = ch_seeds;
tree{next_level, child_index}.SeedsData = {data1 data2};

%% Get parent parameter
m_cov = intrisic_mean( data_cov );
covariance = reshape(m_cov,3,3)';% average_covariance(data);
par.L = L;
par.q =  q;
par.cov = covariance;
tree{next_level, child_index}.Parameter = par;

%% Get entropy
n_bands = size(data_cov, 2);
tree{next_level, child_index}.Entropy = whisart_entropy(par, n_bands);

%% Compute entropy gain - child 1
par.L = L;
par.q =  q;
par.cov = ch_seeds{1};
entropy1 = whisart_entropy(par, n_bands);

%% Get size weigth - child 1
sw1 = 1;%size(data1, 1)/(n_row*n_col);

%% Get entropy gain - child 1
HG1 = abs((tree{next_level, child_index}.Entropy-entropy1)*sw1);

%% Compute entropy gain - child 2
par.L = L;
par.q =  q;
par.cov = ch_seeds{2};
entropy2 = whisart_entropy(par, n_bands);

%% Get size weigth - child 2
sw2 = 1;%size(data2, 1)/(n_row*n_col);

%% Get entropy gain - child 2
HG2 = abs((tree{next_level, child_index}.Entropy-entropy2)*sw2);

HG = mean([HG1 HG2]);
if (HG < 1e-1)
    HG = 0;
end
tree{next_level, child_index}.EntropyGain = HG;
end

