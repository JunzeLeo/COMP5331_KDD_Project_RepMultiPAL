function v = IPD(x, y)
    p = ceil((length(x)-1)/2);
    v = 0;
    for i=1:p
        lami = 2*pi*i/length(x);
        sx = 0;
        sy = 0;
        for t=1:length(x)
            sx = sx + x(t)*exp(-1j*lami*t);
            sy = sy + y(t)*exp(-1j*lami*t);
        end
        sx = abs(sx);
        sy = abs(sy);
        v = v + log(sx)/sx - log(sy)/sy;
    end
    v = sqrt(abs(v))/p;
end

