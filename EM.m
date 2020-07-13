function [index, data1, data2] = EM( imacov_vec, L, q )
s = size(imacov_vec, 1);

[ seed, ~, ~ ] = get_seeds_pca( imacov_vec );
 
%% Get parameter 1---------------------------------------------------------
p{1}.weight = 0.5;
p{1}.L = L;
p{1}.q = q;
p{1}.covariance = seed{1};

%% Get parameter 2---------------------------------------------------------
p{2}.weight = 0.5;
p{2}.L = L;
p{2}.q = q;
p{2}.covariance = seed{2};
   

epsilon = 1e-2; % precision

error = 2;
index = 0;
step = ceil(size(imacov_vec,1)*0.0005);
while (error > epsilon)
    %% Expectation---------------------------------------------------------
    data1 = [];
    data2 = [];
    for ii = 1:step:size(imacov_vec,1)
        
        cov = reshape(imacov_vec(ii, :),3,3)';
        p_cluster1 = log(p{1}.weight) - GetWhishart( cov, p{1} );
        p_cluster2 = log(p{2}.weight) - GetWhishart( cov, p{2} );
                
        if abs(p_cluster1) > abs(p_cluster2)
            data1 = [data1; imacov_vec(ii, :)];
        else
            data2 = [data2; imacov_vec(ii, :)];
        end
        
    end
    
    if ( isempty(data1) || isempty(data2) )
        break
    end
    %----------------------------------------------------------------------
    
    %% Maximization--------------------------------------------------------
    %weigths
    p{1}.weight = 0.5;%size(data1,1) / (size(imacov_vec,1)/10);
    p{2}.weight = 1 - p{1}.weight;
    
    prev_temp = [diag(p{1}.covariance) diag(p{2}.covariance)]';
    if((~isempty(data1)) && (size(data1,1)>1))
        m_cov = intrisic_mean( data1 );
        p{1}.covariance = reshape(m_cov,3,3)';%average_covariance(data1);
    end
    
    if((~isempty(data2)) && (size(data2,1)>1))
        m_cov = intrisic_mean( data2 );
        p{2}.covariance = reshape(m_cov,3,3)';%average_covariance(data2);
    end
    
    actual_temp = [diag(p{1}.covariance) diag(p{2}.covariance)]';
    temp_num = sum( ((prev_temp - actual_temp) .^ 2)' );
    tem_den = sum( (prev_temp.^2)' );
    error_v = temp_num./tem_den;
    error = max(error_v);
    index = index+1;
    if(index > 10)
        break;
    end
end

end

