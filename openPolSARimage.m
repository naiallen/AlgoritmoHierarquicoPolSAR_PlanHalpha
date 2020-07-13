%%=========================================================================
% Brief: Read PolSAR image data
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
function [ output ] = openPolSARimage(folder, im_filename, hdr_filename)
%% Read header
address = fullfile(folder, hdr_filename);
[ n_row, n_col, n_bands, format, interleave, byte_order] = read_hdr( address );

%% Read image
address = fullfile(folder, im_filename);
fileID = fopen(address);%, 'r', 'ieee-le');
if fileID == -1, error('Cannot open file: %s', filename); end
if (byte_order>0)
    Data = fread(fileID, Inf, format, 'ieee-be');
else
    Data = fread(fileID, Inf, format);
end
fclose(fileID);

offset = size(Data,1) - (n_col*n_row*n_bands*2);
Data(1:offset) = [];


%% Mouting Image
ImTemp = zeros(n_row, n_col, n_bands);
for ii = 1:n_bands
    if(strcmp(interleave, ' bsq'))
        ImSize = n_col*n_row*2; %Real and Imag
        clear temp real imag
        temp = Data(1:ImSize,1);
        Data(1:ImSize) = [];
        
        temp_real = reshape(temp(1:2:end-1), n_row, n_col);
        temp_imag = reshape(temp(2:2:end), n_row, n_col);
        ImTemp(:,:,ii) = temp_real + 1i*temp_imag; 
    elseif(strcmp(interleave,' bip'))
        temp_real = reshape(Data( ((ii-1)*2)+1 : ( n_bands*2) : end ,1), n_row, n_col);
        temp_imag = reshape(Data( ((ii-1)*2)+2 : ( n_bands*2) : end, 1), n_row, n_col);
        ImTemp(:,:,ii) = temp_real + 1i*temp_imag; 
    end
end
output = ImTemp(950:1200, 500:900, :);%(950:960, 500:510, :);%
end

