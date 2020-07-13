function plotPlanHalpha(entropy, alpha)

ent_alpha = [ entropy, alpha ];
ent = linspace(0, 1, 100);
alp = linspace(0, 90, 100);
[x,y]=meshgrid(ent,alp);
z = zeros(size(x));
for ii = 2:size(ent,2)
    temp = ent_alpha( (ent_alpha(:,1) > ent(ii-1)) &...
                                    (ent_alpha(:,1) <= ent(ii) ), :);
   
    if(~isempty(temp))
        for jj = 2:size(alp,2)
            num = numel( temp( (temp(:,2) > alp(jj-1)) & (temp(:,2) <= alp(jj)),1));        
            if(num>0)
                z(jj, ii) = num;
            end
        end
    end
end
% th = 10*mean(z(:));
temp = z(:);
% temp(temp>th) = temp(temp>th)*0.1;
th =  max(z(:))*0.005;
temp(temp<th & temp ~=0) = max(z(:))*0.1;
z = reshape(temp, size(z));

s = surf(x,y,z);
s.EdgeColor = 'none';
% colormap(flipud(hot))
colormap(mycolormap());
view(2)
c = colorbar;
step = max(z(:))/10;
num_range = 0:step:max(z(:));

perc = zeros(size(num_range));
label = {'0%', '0%', '0%', '0%', '0%', '0%', '0%', '0%', '0%', '0%', '0%'};;
temp1 = z(:);
temp1(temp1 == 0) = [];
div = size(temp1,1);
for ii = 2:size(num_range,2)    
    temp = ( (temp1>num_range(ii-1))&(temp1<=num_range(ii)) );
    perc(ii) = round( (sum(temp)/div)*100,2);
    label{ii} = strcat(num2str( perc(ii) ), '%');
end

c.Ticks = num_range;
c.TickLabels = label;
hold on
plot(ent_alpha(:,1),ent_alpha(:,2), '.', 'color', [0 0 0.5313])

clear x y
hold on
m = 0:1/512:1;
log_3 =  log10(3.0);

max_z = max(z(:));
clear z
%curva I
a1 = m .* 180 ./ (1 + 2 .* m);
h1 = -( 1 ./(1+2 .* m) ) .* log10( (m .^ (2 .* m)) ./ ( (1 + 2 .* m).^(2 .* m+1) ) ) ./ log_3;
z = ones(size(a1))*max_z;
plot3 (h1, a1, z, 'k', 'LineWidth',2);

%curva II
a2 = zeros(1,513);
h2 = zeros(1,513);
ss = 256;
ss2 = ss+1;
a2(1:ss)  = 90;
a2(ss2:513) = 180 ./ (1 + 2 .* m(ss2:513));
h2(1:ss)    = - ( 1 ./ (1 + 2 .* m(1:ss))) .* log10( (2 .*m(1:ss).^(2 .*m(1:ss))) ./ ( (1 + 2 .*m(1:ss)) .^ (2 .* m(1:ss)+1) ) )./log_3;
h2(ss2:513) = - ( 1 ./ (1 + 2 .* m(ss2:513)) ).*log10( (2.*m(ss2:513)-1).^(2.*m(ss2:513)-1) ./ ((1+2 .* m(ss2:513)) .^(2 .* m(ss2:513)+1)) )./log_3;
z = ones(size(a2))*max_z;
hold on;
plot3 (h2,a2, z, 'k', 'LineWidth',2);

%linha 1
x = m .* 0.5;
y(1:513) = 42.5;
z = ones(size(y))*max_z;
hold on;
plot3 (x, y, z, ':k', 'LineWidth',1);

%linha 2
y(1:513) = 47.5;
hold on;
plot3 (x, y, z, ':k', 'LineWidth',1);

%linha 3
x = 0.5+m.*(0.4);
y(1:513)= 50;
hold on;
plot3 (x, y, z, ':k', 'LineWidth',1);

%linha 4
y(1:513)= 40;
hold on;
plot3 (x, y, z, ':k', 'LineWidth',1);

%linha 5
% x = 0.9 + m.*(0.1);
% hold on;
% plot (x, y, 'k', 'LineWidth',2);

%linha 6
x = 0.9 + m.*(0.1);
y(1:513) = 55;
hold on;
plot3 (x, y,z, ':k', 'LineWidth',1);

%coluna 1
x(1:513) = 0.5;
y = m.*90;
y(y<14) = 14;
hold on;
plot3 (x, y,z, ':k', 'LineWidth',1);

%coluna 2
x(1:513) = 0.9;
y = m.*90;
y(y<38) = 38;
y(y>78) = 78;
hold on;
plot3 (x, y,z, ':k', 'LineWidth',1);

axis([0 1 0 90])

end

