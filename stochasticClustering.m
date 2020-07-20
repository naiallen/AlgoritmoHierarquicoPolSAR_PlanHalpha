function [index, output_im] = stochasticClustering(input_im, distance, n_classes, seeds, iteration, L)

ampl_max = [max(input_im(:,1)) max(input_im(:,5)) max(input_im(:,9))];
ampl_min = [min(input_im(:,1)) min(input_im(:,5)) min(input_im(:,9))];
error_th = mean(abs(ampl_max - ampl_min))*0.05;
[n_row, n_col, n_bands] = size(input_im);
index =1;
counter = 1;
data_total = im2vec( input_im, n_bands);
error_prev = 0;
while (index <= iteration)
    output_im = ones(n_row, n_col, n_classes)*nan;
    for ii = 1:n_row          %linha
        for jj = 1:n_col      %coluna
            if isnan(input_im(ii, jj,:))
                continue;
            end

            %% Get covariance
            cov1 = reshape(input_im(ii, jj,:), 3,3)';
            dist = zeros(2,1);
            for kk = 1:n_classes
                cov2 = seeds{kk};
                dist(kk, 1) = abs( stochastic_distance(distance, cov1, cov2, L) );
            end
            [~, pos] = min(dist);
            output_im(ii, jj, pos) = 1;
            
        end
    end
    
    %% Geet new seeds
    ampl_prev = []; 
    ampl_actual = []; 
    
    for kk=1:n_classes
        ampl_prev = [ampl_prev seeds{kk}(1,1) seeds{kk}(2,2) seeds{kk}(3,3)];
        im = zeros(n_row, n_col, n_bands);
        for nb = 1:n_bands
            im(:, :, nb) = output_im(:,:, kk).*(input_im(:,:,nb));
        end
        clear data
        data = im2vec( im, n_bands);
        if((~isempty(data)) && (size(data,1)>1))
            m_cov = intrisic_mean( data );
            seeds{kk} = reshape(m_cov,3,3)';
%             seeds{kk} = average_covariance(data);
        else
            s = size(data_total, 1);
            index1 = randi(s);
%             data = [ reshape(seeds{kk}, 1, 9); data_total(index1, :)];
            m_cov = data_total(index1, :);%intrisic_mean( data ); %
            seeds{kk} = reshape(m_cov,3,3)';
%             seeds{kk} = average_covariance(data);
            counter = counter+1;
        end
        ampl_actual = [ampl_actual seeds{kk}(1,1) seeds{kk}(2,2) seeds{kk}(3,3)];
%         seeds{kk} = average_covariance(data);
    end 
    %% Compute RSME
    error = sqrt( mean((abs(ampl_prev.^2 - ampl_actual.^2))));
    if((error < error_th)||...
        (index >= 20)||...
        (counter >= 3) )
        break;
    end
%     error_prev = error;
    index=index+1;
end
index=index-1;
end

