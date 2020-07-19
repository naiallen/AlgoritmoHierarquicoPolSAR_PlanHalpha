function [dataz] = zscore(data)

[~, ncol] = size(data);
dataz = zeros(size(data));
for ii = 1:ncol
    dataz(:, ii) = (data(:,ii) - mean(data(:,ii)))./std(data(:,ii));
end
end

