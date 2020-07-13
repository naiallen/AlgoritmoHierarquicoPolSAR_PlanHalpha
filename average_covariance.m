function covariance = average_covariance(imacov_vec)

%--------------------------------------------------------------------------
c11 = mean(imacov_vec(:, 1));
c21 = mean(imacov_vec(:, 2));
c31 = mean(imacov_vec(:, 3));

c12 = mean(imacov_vec(:, 4));
c22 = mean(imacov_vec(:, 5));
c32 = mean(imacov_vec(:, 6));

c13 = mean(imacov_vec(:, 7));
c23 = mean(imacov_vec(:, 8));
c33 = mean(imacov_vec(:, 9));

covariance = [c11 c12 c13; c21 c22 c23; c31 c32 c33];
%================================================

end

