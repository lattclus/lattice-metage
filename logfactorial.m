function logmat = logfactorial(mat)
% mat is a 1/2/3 dimensional array

nd = ndims(mat);
mat(find(mat==0))=1; % since factorial of zero = 1

if(nd==1)
    a = length(mat);
    n = max(mat);
    logmat = zeros(a);
else if(nd ==2)
    [a,b]=size(mat);
    n = max(max(mat));
    logmat = zeros(a,b);
    else
        [a,b,c] = size(mat);
        n = max(max(max(mat)));
        logmat = zeros(a,b,c);
    end
end


logfac = zeros(1,n);
i=2;

while i<=n
    logfac(i) = logfac(i-1)+log(i);
    i=i+1;
end

logmat = logfac(mat);
