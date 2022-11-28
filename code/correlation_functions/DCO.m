function v = DCO(x, y)
    CORR = sum(diff(x).*diff(y));
    CORR = CORR/(sqrt(sum(diff(x).^2))*sqrt(sum(diff(y).^2)));
    CORT = sum((x-mean(x)).*(y-mean(y)));
    CORT = CORT/(sqrt(sum((x-mean(x)).^2))*sqrt(sum((y-mean(y)).^2)));
    v = CORR*CORT;
end

