angles = [0:0.2:2*pi]';
%angles = angles(randn(numel(angles),1) > 0);
N = numel(angles);

r = ones(N,1);
%r = r + randn(N,1)*0.03;
rin = 0.1;
rout = 0.1;

Xon = [r.*cos(angles),r.*sin(angles)];
Von = zeros(N,1);

Xin = [(r+rin).*cos(angles),(r+rin).*sin(angles)];
Vin = ones(N,1)*100;

Xout = [(r-rout).*cos(angles),(r-rout).*sin(angles)];
Vout = -ones(N,1)*10;

X = [Xon;Xin;Xout];
V = [Von;Vin;Vout];


[xi,yi,F] = thinPlateSpline(X,V, 1);
%surf(xi,yi,abs(F).^(1/4)); colormap gray;

surf(xi,yi,abs(F).^(1/2)); colormap hot;


newplot;
hold on;
contour(xi,yi,F,[0 0],'r-');
plot(X(:,1),X(:,2), 'go'); 
axis equal;
hold off;


Xout2 = [1.5*cos(angles),1.5*sin(angles)];
Vout2 = ones(N,1)*40;
X = [Xon;Xin;Xout;Xout2];
V = [Von;Vin;Vout;Vout2];




hold on
plot3(X(:,1),X(:,2),zeros(3*N,1), '*r')
grid
stem3(X(:,1),X(:,2),V,'^','fill')
hold off