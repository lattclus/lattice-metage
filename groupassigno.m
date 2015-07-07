function  [C,groupcount] = groupassign(Q,theta,M,L,wordcountfile)

load(wordcountfile);
[b,d,n] = size(X); % b=#c_k, d=#words, n=#reads

c = zeros(1,L);
for m=1:M,
    for j=1:d,
        for l=1:L
            y = squeeze(X(:,j,:));
            y = bsxfun(@times,y,log(squeeze(theta(m,l,:))));
            c(l) = sum(y)*Q(:,m); % sum over all c_k: n_i(c_k|w_j)*log(theta)
        end
        [xx,C(m,j)] = max(c);
    end
end
clear y X;
groupcount = getgroupcount(C,L,wordcountfile);