function [entropia, alpha] = planHalpha( data_coh)
%% Calcula a entropia e angulo alpha
entropia     = zeros(size(data_coh,1), 1);
alpha        = zeros(size(data_coh,1), 1);
% anisotropia  = zeros(size(data_coh,1), 1);

for ii = 1:size(data_coh,1)
    m  = conj(reshape(data_coh(ii,:),3,3));
    [V, D] = eig(m);
       
    %Ordena os autovalores do maior pro menor
    [autovalor, b] = sort(diag(D), 'descend');
    autovetor = zeros(size(V));
    for ev = 1:size(autovalor, 1)
        autovetor(:, ev) = V(:, b(ev));
    end
    
    %Calcula a entropia
    soma_lambda = sum(autovalor);
    PP = zeros(3,1);
    PP(1) = single( autovalor(1)/soma_lambda);
    PP(2) = single( autovalor(2)/soma_lambda);
    PP(3) = single( autovalor(3)/soma_lambda);
    
    for ev = 1:size(autovalor, 1)
        if(PP(ev)>0)
            entropia(ii) = entropia(ii) - ( PP(ev)*(log10(PP(ev))/log10(3)) );
        end
    end
    entropia(ii) = abs(entropia(ii));
    
    %Angulo alpha
    a = zeros(3,1);
    a(1) = acos(abs(autovetor(1, 1)));
    a(2) = acos(abs(autovetor(1, 2)));
    a(3) = acos(abs(autovetor(1, 3)));
    alpha(ii) = abs(a(1)*PP(1) + a(2)*PP(2) + a(3)*PP(3));
    
%     %Anisotropia
%     anisotropia(ii) = abs((autovalor(2) - autovalor(3)) / (autovalor(2) + autovalor(3)));
    
end
alpha = abs(rad2deg(alpha));

% %Plota dados
% hold on
% plot(entropia(:),  alpha(:), '.', 'color', cl)
end

