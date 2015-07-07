function [y] = logcondmultinomialpdf(X,theta_mod)
% output: sum over c_k and w_j (n_i(w_jc_k)*log(theta))
% changes: input X = count of n_i(c_k |l) 

if size(X,2)==size(theta_mod,2)
    theta_mod(theta_mod==0) = 1;
    c = bsxfun(@times,X,log(theta_mod)); % len(c_k) 4^k read_i. %% instead of 4^k it should be #groups
    c = sum(c,1); % sum over c_k so 1 4^k read_i
    y = sum(squeeze(c)); % sum over w_j so a #reads vector %% summed over groups
else
    error('Parameters mismatched.');
end



       