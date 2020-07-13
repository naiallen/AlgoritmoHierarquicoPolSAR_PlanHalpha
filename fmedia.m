%%=========================================================================
% Brief: Filter PolSAR images - Get a N-look image
%   Input: 
%       input_img: 3 bands image
%       L: Number of Looks
%   Output
%       output_image: 3 bands image
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================

%%
function output_image = fmedia (input_img, L)

% Ceate output image
output_image = input_img;

[n_row, n_col, n_bands] = size(input_img);
% Contex
step = floor(L/2);

for row = 1+step:n_row-step
    for col = 1+step:n_col-step
        for bands = 1:n_bands
            output_image(row, col, bands) = mean(mean(input_img(row-step:row+step, col-step:col+step, bands)));
        end
    end
end
end