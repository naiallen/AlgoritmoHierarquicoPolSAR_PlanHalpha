function [index, seed, data1, data2] = get_seed_EM(imacov_vec,  L, q)

[index, data1, data2] = EM( imacov_vec, L, q );

if (isempty(data1))  
    %% Get Seed 2----------------------------------------------------------
    seed{2}= average_covariance(data2);
    
    %% Get Seed 1----------------------------------------------------------
    seed{1} = seed{2};
    data1 = data2(1:3,:);
elseif (isempty(data2))
    %% Get Seed 1----------------------------------------------------------
    seed{1} = average_covariance(data1);
    
    %% Get Seed 2----------------------------------------------------------
    seed{2}= seed{1};
    data2 = data1(1:3,:);
else
    %% Get Seed 1----------------------------------------------------------
    seed{1} = average_covariance(data1);
    
    %% Get Seed 2----------------------------------------------------------
    seed{2}= average_covariance(data2);
end

end

