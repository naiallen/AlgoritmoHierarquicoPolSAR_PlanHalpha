function m_cov = intrisic_mean( imacov_vec )
%% Intrisic mean
n_bands = size(imacov_vec, 2);
N = size(imacov_vec, 1);
m_cov_temp = sum(imacov_vec)/N;
m_cov = reshape(m_cov_temp,3,3);
error = 100;
index = 0;
while(error > 0.01)
    delta_m = zeros(sqrt(n_bands), sqrt(n_bands));
    for ii =1:N
        A = reshape(imacov_vec(ii,:),3,3);
        delta_m = delta_m + (logm(A));
    end
    delta_m = delta_m/N;
    m_cov_prev = m_cov;
    m_cov = expm(delta_m);    
    error = abs( trace(m_cov_prev) - trace(m_cov) );
    index = index+1;
    if (index > 5)
        break;
    end
%     delta_m = zeros(1,n_bands);
%     for kk =1:n_bands
%         delta_m(kk) = sum( log(imacov_vec(:,kk))/ log(m_cov(:,kk)) )/N;
%     end
%     m_cov = m_cov.^(delta_m);   
end
m_cov = reshape(m_cov, 1, n_bands);

end

