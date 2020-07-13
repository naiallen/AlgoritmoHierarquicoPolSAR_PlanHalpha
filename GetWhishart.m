function output = GetWhishart( cov1, p )

output = p.L*log(det(p.covariance)) + (p.L-p.q)*log(det(cov1)) - p.L*trace((p.covariance^(-1))*cov1);

end

