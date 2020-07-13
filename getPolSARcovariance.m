%%=========================================================================
% Brief: Compute the PolSAR covariance
%   Input: 
%       PolSAR image: 3 bands image
%   Output
%       Covariance Image: 9 bands image
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================

%%
function [ cov ] = getPolSARcovariance( polSAR_im )
Shh = polSAR_im(:,:,1);
Shv = polSAR_im(:,:,2);
Svv = polSAR_im(:,:,3);

cov = zeros(size(Shh,1), size(Shh,2), 9);
raiz = sqrt(2);
cov(:, :, 1) = abs(Shh.*Shh);
cov(:, :, 4) = Shh.*conj(Shv)*raiz;
cov(:, :, 7) = Shh.*conj(Svv);

cov(:, :, 2) = Shv.*conj(Shh)*raiz;
cov(:, :, 5) = 2*abs(Shv.*Shv);
cov(:, :, 8) = Shv.*conj(Svv)*raiz;

cov(:, :, 3) = Svv.*conj(Shh);
cov(:, :, 6) = Svv.*conj(Shv)*raiz;
cov(:, :, 9) = abs(Svv.*Svv);
end

