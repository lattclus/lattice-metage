function  [C,groupcount] = groupassign(Q,theta,I,theta_back,M,L,wordcountfile)

load(wordcountfile);
[b,d,n] = size(X); % b=#c_k, d=#words, n=#reads

c = zeros(1,L);

%theta_final = theta;

for b1=1:b
    [t1,t2] = find(squeeze(I(:,:,b))==1);
    theta(t1,t2,b1) = theta_back(t1,t2,b1);
end
% for m=1:M
%     for l=1:L
%         for b1=1:b
%             if I(m,l,b1) == 1
%                 theta_final = theta_back;
%             else 
%                 theta_final = theta;
%         end
%     end
% end
for m=1:M
    for j=1:d
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