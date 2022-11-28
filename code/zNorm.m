function y = zNorm(x)
    y = (x-mean(x))/std(x, 1);
end