function H = whisart_entropy( parameters, n_bands)
g = 0;
for ii=1:parameters.q
    g = g+(gamma((parameters.L-(ii-1))/2));
end
H = n_bands/2 + (n_bands/2)*log(pi) + 0.5*log(det(parameters.cov)) - 0.5*g;
end

