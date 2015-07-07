function [z] = logcondmultinomialpdf(X,theta_mod,Im, theta_b)
% output: sum over c_k and w_j (n_i(w_jc_k)*log(theta))
% changes: input X = count of n_i(c_k |l) 
% theta_mod is [c_k,l], X is [c_k,l,i], I is [l,c_k]

if size(X,2)==size(theta_mod,2)
    %size(theta_b)
    %size(Im')
    %size(theta_mod)
    %size(theta_b)
    %size((1-Im)')
    %size(X)
    %size(log(theta_mod))
    c = bsxfun(@times,X,(1-Im)'.*log(theta_mod));
    theta_b(theta_b==0)=1;
        c = c + bsxfun(@times,X,Im'.*log(theta_b) ); % len(c_k) 4^k read_i. %% instead of 4^k it should be #groups
    
    c = sum(c,1); % sum over c_k so 1 4^k read_i
    y = sum(squeeze(c)); % sum over w_j so a #reads vector %% summed over groups
    z=y';
else
    error('Parameters mismatched.');
end



       