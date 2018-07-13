function [thepsi,thegrad, thehess] = mypsi(b, beta, sigma, x, y)
ld = -0.5 * real(logdet(sigma));
bd = b - beta;
dp = bd * (sigma\bd);
sz = size(x);
nr = sz(1);
auxlogexp = zeros(1, nr);
auxlog = zeros(1, nr)
for s=1:nr
    tmp = 0;
    for j = 1:nr
        xx = x(s, j, :);
        tmp = tmp + exp(xx * b);
    end
    auxlog(s) = tmp;
    auxlogexp(s) = -log(tmp);
end
tmp = zeros(1, nr);
for j=1:nr
    tmpj =0;
    aux2 = zeros(1, nr);
    for s = 1:nr
        xx = x(s, j, :);
        yy = y(s, j);
        aux2(s) = yy*(xx * b + auxlogexp(s));
        tmpj = tmpj + aux2(s);
    end
    tmp(j) = tmpj;
end
thesum = sum(tmp);
thepsi = ld + dp + thesum;
if nargout == 1
    return
end

end
    
function themat = tenscomm(x, i, j, k)
x1 = x(i, j, :);
x2 = x(i, k, :);
themat = x1 ( x1' - x1 * x2';
end

function em = expmat(x, b)
em = cellfun(@(xx) exp(xx * b), x);
end