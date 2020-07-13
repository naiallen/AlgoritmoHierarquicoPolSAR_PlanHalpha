%%=========================================================================
% Brief: Get the initial seeds using PCA
%   Input:
%       Covariance Image: 9 bands image
%   Output
%       seed, data1, data2
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================

%%
function [ seed, data1, data2 ] = get_seeds_random( imacov_vec )

s = size(imacov_vec, 1);

index1 = randi(s);
index2 = randi(s);

data1 = imacov_vec(index1,:);
data2 = imacov_vec(index2,:);


%% Get Seed 1--------------------------------------------------------------
seed{1} = reshape(imacov_vec(index1,:), 3, 3)';

%% Get Seed 2--------------------------------------------------------------
seed{2} =  reshape(imacov_vec(index2,:), 3, 3)';

% %Plot data---------------------------------------------------------------
% h2 = figure('Position', [10 10 900 600]);
%     
% chisquare_val = 2.7955; %sqrt(7.815)
% 
% a = chisquare_val*sqrt(eigenval(1));
% b = chisquare_val*sqrt(eigenval(2));
% c = chisquare_val*sqrt(eigenval(3));
% 
% %Transform points to x,y,z
% x = [m(1) m(1)+b*eigenvec(1,2) m(1)+c*eigenvec(1,3)];
% y = [m(2) m(2)+b*eigenvec(2,2) m(2)+c*eigenvec(2,3)];
% z = [m(3) m(3)+b*eigenvec(3,2) m(3)+c*eigenvec(3,3)];
% 
% %Initial Points
% X0 = m(1);
% Y0 = m(2);
% Z0 = m(3);
% 
% %End Points
% xLim = [min(x) max(x)];
% zLim = [min(z) max(z)];
% [X,Z] = meshgrid(xLim,zLim);
% Y = (A * X + C * Z + D)/ (-B);
% reOrder = [ 4 3 1 2];
% X = X(reOrder);
% Y = Y(reOrder);
% Z = Z(reOrder);
% 
% plot3( real(data1(1:1:end, 1)),  real(data1(1:1:end, 5)), real(data1(1:1:end, 9)), '.m')
% hold on
% plot3( real(data2(1:1:end, 1)),  real(data2(1:1:end, 5)), real(data2(1:1:end, 9)), '.c')
% 
% quiver3(X0, Y0, Z0, real(eigenvec(1,1))*sqrt(eigenval(1)),  real(eigenvec(2,1))*sqrt(eigenval(1)), real(eigenvec(3,1))*sqrt(eigenval(1)), '-k', 'LineWidth',2);
% hold on;
% quiver3(X0, Y0, Z0, real(eigenvec(1,2))*sqrt(eigenval(2)),  real(eigenvec(2,2))*sqrt(eigenval(2)), real(eigenvec(3,2))*sqrt(eigenval(2)), '-k', 'LineWidth',2);
% hold on;
% quiver3(X0, Y0, Z0, real(eigenvec(1,3))*sqrt(eigenval(3)),  real(eigenvec(2,3))*sqrt(eigenval(3)), real(eigenvec(3,3))*sqrt(eigenval(3)), '-k', 'LineWidth',2);
% hold on
% fill3(X, Y, Z,'b');
% grid on;
% alpha(0.3);

end

