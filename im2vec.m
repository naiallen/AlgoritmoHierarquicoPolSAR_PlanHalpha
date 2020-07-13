function [ vector ] = im2vec( image, n_bands)

vector = [];
for ii =1:n_bands
    clear aux
    aux =  image(:, :, ii); aux=(aux(~isnan(aux)));
    vector = [vector aux(:)];
end

end

