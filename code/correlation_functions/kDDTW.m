function v = kDDTW(x, y)
xp1 = x(3:end);
xm1 = x(1:end-2);
x = x(2:end-1);
yp1 = y(3:end);
ym1 = y(1:end-2);
y = y(2:end-1);
v = mean(((xp1+2*x-3*xm1)/4-(yp1+2*y-3*ym1)/4).^2);
end

