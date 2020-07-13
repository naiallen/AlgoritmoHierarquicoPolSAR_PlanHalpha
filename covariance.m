function sigma = covariance(data, n_bands)
    if (n_bands > 1 && size(data, 1) > 1)
        m =  mean( data );
        for ii = 1:n_bands
            data(:, ii) = data(:, ii) - m(ii);
        end
        sigma = ( 1/(size(data, 1)-1) )*(data'*data);
    else
        sigma = nan;
    end
end 