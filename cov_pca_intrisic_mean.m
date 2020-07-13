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
function [ seed, data1, data2] = cov_pca_intrisic_mean( imacov_vec)

%% Compute parameters -------------------------------------------------
data = [sqrt(imacov_vec(:,1)) sqrt(imacov_vec(:,5)) sqrt(imacov_vec(:,9))];
dataz = abs(zscore(data));
index = zeros(size(data,1),1);
index( (dataz(:,1)>2.8) | (dataz(:,2)>2.8) | (dataz(:,3)>2.8),:) = 1;
imacov_vec(index==1,:) = [];
m_cov = intrisic_mean( imacov_vec );
m = [sqrt(m_cov(1)) sqrt(m_cov(5)) sqrt(m_cov(9))];
Sigma = reshape(m_cov,3,3)';%average_covariance(imacov_vec);

%% Compue eigenvalue and eigenvector---------------------------------
[V, D] = eig(Sigma);

%% Sort eigenvalue---------------------------------------------------
[eigenval, b] = sort(diag(D), 'descend');

%% Sort eigenvector--------------------------------------------------
eigenvec = zeros(size(V));
for ev = 1:size(eigenval,1)
    eigenvec(:, ev) = V(:,b(ev));
end

%% Compute the normal plane to the biggest eigenvalue----------------
plan = cross(eigenvec(:, 2), eigenvec(:, 3));

%% Find plan equation coefficients-----------------------------------
A = plan(1);
B = plan(2);
C = plan(3);
D = -dot(plan,m);

%% Divide data given the plan----------------------------------------------
d = A*sqrt(imacov_vec(:, 1)) + B*sqrt(imacov_vec(:, 5)) + C*sqrt(imacov_vec(:, 9)) + D;
data1 = zeros(size(imacov_vec,1), 1)*NaN;
data2 = zeros(size(imacov_vec,1), 1)*NaN;
% data1(d>0,:) = 1;
% data2(d<0,:) = 1;

data1 = imacov_vec(d>0,:);
data2 = imacov_vec(d<0,:);

% data1 = imacoh_vec(d>0,:);
% data2 = imacoh_vec(d<0,:);

%% Get Seed 1--------------------------------------------------------------
m_cov = intrisic_mean(  imacov_vec(d>0,:) );
seed{1} = reshape(m_cov,3,3)';

%% Get Seed 2--------------------------------------------------------------
m_cov = intrisic_mean( imacov_vec(d<0,:) );
seed{2} =  reshape(m_cov,3,3)';

% %% Plot data---------------------------------------------------------------
% h2 = figure('Position', [10 10 900 600]);
% 
% subplot(224)
% chisquare_val = 2.7955; %sqrt(7.815)
% a = chisquare_val*sqrt(eigenval(1));
% b = chisquare_val*sqrt(eigenval(2));
% c = chisquare_val*sqrt(eigenval(3));
% 
% %Transform points to x,y,z
% x = [m_cov(1) m_cov(1)+b*eigenvec(1,2) m_cov(1)+c*eigenvec(1,3)];
% y = [m_cov(2) m_cov(2)+b*eigenvec(2,2) m_cov(2)+c*eigenvec(2,3)];
% z = [m_cov(3) m_cov(3)+b*eigenvec(3,2) m_cov(3)+c*eigenvec(3,3)];
% 
% %Initial Points
% X0 = m_cov(1);
% Y0 = m_cov(2);
% Z0 = m_cov(3);
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
% plot3( (data1(1:1:end, 1)),  (data1(1:1:end, 5)), (data1(1:1:end, 9)), '.m')
% hold on
% plot3( (data2(1:1:end, 1)),  (data2(1:1:end, 5)), (data2(1:1:end, 9)), '.c')
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
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% 
% subplot(221)
% plot3( (data1(1:1:end, 1)),  (data1(1:1:end, 5)), (data1(1:1:end, 9)), '.m')
% hold on
% plot3( (data2(1:1:end, 1)),  (data2(1:1:end, 5)), (data2(1:1:end, 9)), '.c')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% view(2)
% grid on;
% 
% subplot(222)
% plot3( (data1(1:1:end, 1)),  (data1(1:1:end, 5)), (data1(1:1:end, 9)), '.m')
% hold on
% plot3( (data2(1:1:end, 1)),  (data2(1:1:end, 5)), (data2(1:1:end, 9)), '.c')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% view(90,0)
% grid on;
% 
% %Plot ZY plan
% subplot(223)
% plot3( (data1(1:1:end, 1)),  (data1(1:1:end, 5)), (data1(1:1:end, 9)), '.m')
% hold on
% plot3( (data2(1:1:end, 1)),  (data2(1:1:end, 5)), (data2(1:1:end, 9)), '.c')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% view(0,90)
% grid on;

end

