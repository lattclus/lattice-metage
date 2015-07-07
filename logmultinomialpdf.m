function y = logmultinomialpdf(X, theta)

[~,k] = size(theta);

minlen = min(sum(X));

if k==1
    c = bsxfun(@times,X,log(theta));
    c = sum(c - log(factorial(X)));
    y = bsxfun(@plus,c,logfacbyfac(sum(X),minlen)); 
else
    error('Parameters mismatched.');
end

