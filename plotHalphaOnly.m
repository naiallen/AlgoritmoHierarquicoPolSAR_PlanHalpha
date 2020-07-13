function plotHalphaOnly
clear x y
hold on
m = 0:1/512:1;
log_3 =  log10(3.0);

%curva I
a1 = m .* 180 ./ (1 + 2 .* m);
h1 = -( 1 ./(1+2 .* m) ) .* log10( (m .^ (2 .* m)) ./ ( (1 + 2 .* m).^(2 .* m+1) ) ) ./ log_3;
plot(h1, a1, 'k', 'LineWidth',2);

%curva II
a2 = zeros(1,513);
h2 = zeros(1,513);
ss = 256;
ss2 = ss+1;
a2(1:ss)  = 90;
a2(ss2:513) = 180 ./ (1 + 2 .* m(ss2:513));
h2(1:ss)    = - ( 1 ./ (1 + 2 .* m(1:ss))) .* log10( (2 .*m(1:ss).^(2 .*m(1:ss))) ./ ( (1 + 2 .*m(1:ss)) .^ (2 .* m(1:ss)+1) ) )./log_3;
h2(ss2:513) = - ( 1 ./ (1 + 2 .* m(ss2:513)) ).*log10( (2.*m(ss2:513)-1).^(2.*m(ss2:513)-1) ./ ((1+2 .* m(ss2:513)) .^(2 .* m(ss2:513)+1)) )./log_3;
hold on;
plot (h2,a2, 'k', 'LineWidth',2);

%linha 1
x = m .* 0.5;
y(1:513) = 42.5;
hold on;
plot (x, y, ':k', 'LineWidth',1);

%linha 2
y(1:513) = 47.5;
hold on;
plot (x, y, ':k', 'LineWidth',1);

%linha 3
x = 0.5+m.*(0.4);
y(1:513)= 50;
hold on;
plot(x, y, ':k', 'LineWidth',1);

%linha 4
y(1:513)= 40;
hold on;
plot (x, y, ':k', 'LineWidth',1);

%linha 5
% x = 0.9 + m.*(0.1);
% hold on;
% plot (x, y, 'k', 'LineWidth',2);

%linha 6
x = 0.9 + m.*(0.1);
y(1:513) = 55;
hold on;
plot (x, y, ':k', 'LineWidth',1);

%coluna 1
x(1:513) = 0.5;
y = m.*90;
y(y<14) = 14;
hold on;
plot (x, y, ':k', 'LineWidth',1);

%coluna 2
x(1:513) = 0.9;
y = m.*90;
y(y<38) = 38;
y(y>78) = 78;
hold on;
plot (x, y, ':k', 'LineWidth',1);

axis([0 1 0 90])
end