function theta = multidimprod(X,Q)

[b,d,~] = size(X);
[~,k]= size(Q);
theta = zeros(k,d,b);

for m=1:k,
    for pa_j=1:d,
        for c_k=1:b,
            theta(m,pa_j,c_k) = Q(:,m)'*squeeze(X(c_k,pa_j,:));
        end
        theta(m,pa_j,:) = theta(m,pa_j,:)/sum(squeeze(theta(m,pa_j,:)));
    end
end