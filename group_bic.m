[b,d,n] = size(X);
X1 = squeeze(sum(X,1));% word counts for each read
L = 0;
for i=1:10
L=L+20;
C = zeros(M,d);%repmat(randsample(L,d,'true')',M,1); %words assigned to a group
% For a species m
clear ll v f;
for m=1:M;
    X2 = X1(:,find(Q(:,m)==1));%all reads in species
    f = [squeeze(sum(X2,2)) squeeze(var(X2'))'];% some read assignment done.
    C(m,:) = kmeans(f,L,'EmptyAction','singleton');
    vall = [var(f(:,1),1); var(f(:,2),1)];
    clear D W Nc X3 int
    for l=1:L;
        X3 = X2(find(C(m,:)==l),:);
        Nc(m,l) = length(find(C(m,:)==l));%#elements in cluster l
        if isempty(find(C(m,:)==l));
            v(:,l)= [0;0];
            D(l) = 0;
            int(l) = 0;
        else
            v(:,l) = [var(f(find(C(m,:)==l),1),1) ; var(f(find(C(m,:)==l),2),1)];
            D(l) = norm(X3-ones(size(X3,1),1)*mean(X3));
            int(l) = D(l)/(2*Nc(m,l));
        end
        ll(l) = -Nc(m,l)*sum(log(v(:,l)+vall)/2);
    end
    BIC(m) = -2 * sum(ll) + 2*L * log(n);
    %int = D./(2*squeeze(Nc(m,:)));
    W(m) = sum( int);
    %AIC(m) = -2 * sum(ll) + 2*L;
end
BICa(i) = mean(BIC)
Wmean(i) = mean(W)
end