function value = stochastic_distance(distance_name, cov1, cov2, L)

cov1 = cov1+eye(3)*0.001;
cov2 = cov2+eye(3)*0.001;

q = size(cov1,1);
det1 = (det( cov1));
det2 = (det( cov2));

%% Check if determinant is zero -> matrix has no inverse
if ((abs(det1) == 0) || (abs(det2) == 0) || isnan(det1) || isnan(det2) )
    value = 0;
    return
end

inv1 = (cov1)^(-1);
inv2 = (cov2)^(-1);

%% Bhattacharyya
if ( strcmp(distance_name, 'Bhattacharyya') )
    value = L*( ( ( log(det1) + log(det2) )/2 ) - log(det( ((inv1+inv2)/2)^(-1) ) ) );
    
%% Kullback_Leibler
elseif(strcmp(distance_name, 'Kullback-Leibler') )
    value = L*( trace( (inv1*cov2)+(inv2*cov1) )/2 - q);

%% Hellinger
elseif(strcmp(distance_name, 'Hellinger') )
    value = 1 - ( det( 2.0^(-1)*((inv1+inv2))^(-1))/(sqrt(det1*det2)) )^L;
    
    %% Renyi
elseif(strcmp(distance_name, 'Renyi') )
    bt = 0.9;
    fir = ( ( det1^(-bt) )*( det2^(bt-1) )*det( inv( bt*inv1 + (1-bt)*inv2 ) ) )^L;
    sec = ( ( det1^(bt-1) )*( det2^(-bt) )*det( inv( bt*inv2 + (1-bt)*inv1 ) ) )^L;
    value = ( log10(2)/(1-bt) )+ ( (1/(bt-1))*(fir+sec) );
    
    %% Chi-quadrado
elseif(strcmp(distance_name, 'Qui-quadrado') )
    fir = ( (det1/(det2^2))*abs(det(inv(2*inv2 - inv1))) )^L;
    sec = ( (det2/(det1^2))*abs(det(inv(2*inv1 - inv2))) )^L;
    value = 0.25*(fir + sec - 2);
    
    %% Euclidean Distance
elseif strcmp(distance_name, 'Euclidean')
    mu_1 = diag(cov1);
    mu_2 = diag(cov2);
    value = sqrt(sum((mu_1 - mu_2 )));
end

end