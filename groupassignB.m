function  [C,groupcount] = groupassignB(Q,theta,L,wordcount2)

%load(wordcountfile);
[b,d,n,M] = size(wordcount2); % b=#c_k, d=#words, n=#reads, M=#species

c = zeros(1,L);

for m=1:M
    for j=1:d
        for l=1:L
            y = squeeze(wordcount2(:,j,:,m));
            y = bsxfun(@times,y,log(squeeze(theta(m,l,:))));
            c(l) = sum(y)*Q(:,m); % sum over all c_k: n_i(c_k|w_j)*log(theta)
        end
        [xx,C(m,j)] = max(c);
    end
end
clear y;
groupcount = getgroupcountB(C,L,wordcount2);