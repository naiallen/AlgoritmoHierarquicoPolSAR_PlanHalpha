%%=========================================================================
% Brief: Read PolSAR covariance image data
%   Input: 
%       folder: Image SAR data folder address
%       im_filename:  Image SAR data file name
%       hdr_filename: Image SAR data header file name
%   Output
%       PolSAR Image: n bands image
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================
function [ output ] = openCOVimage(folder, im_filename, hdr_filename)
%% Read header
address = fullfile(folder, hdr_filename);
[n_row, n_col, n_bands, format] = read_hdr( address );

%% Read image
address = fullfile(folder, im_filename);
fileID = fopen(address);%, 'r', 'ieee-le');
if fileID == -1, error('Cannot open file: %s', filename); end
Data = fread(fileID, Inf, format);
fclose(fileID);

offset = size(Data,1) - (n_col*n_row*n_bands*2);
Data(1:offset) = [];


%% Mouting Image
ImTemp = zeros(n_row, n_col, n_bands);
for ii = 1:n_bands    
    clear temp_real temp_img real imag
    temp_real = Data( ((ii-1)*2)+1 : ( n_bands*2) : end ,1);
    temp_img =  Data( ((ii-1)*2)+2 : ( n_bands*2) : end, 1);

    real = reshape(temp_real, n_row, n_col);
    imag = reshape(temp_img, n_row, n_col);
    ImTemp(:,:,ii) = real + 1i*imag; 
end
output = ImTemp(600:end, 1:400, :);
end

