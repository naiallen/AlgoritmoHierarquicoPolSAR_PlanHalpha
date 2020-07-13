function m_cov = intrisic_mean( imacov_vec )
%% Intrisic mean
n_bands = size(imacov_vec, 2);
N = size(imacov_vec, 1);
m_cov = sum(imacov_vec)/N;
for jj =1:10
    
    delta_m = zeros(1,n_bands);
    for kk =1:n_bands
        delta_m(kk) = sum( log(imacov_vec(:,kk))/ log(m_cov(:,kk)) )/N;
    end
    m_cov = m_cov.^(delta_m);   
end

end

