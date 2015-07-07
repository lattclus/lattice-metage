% Script to generate data to check the model.
M=5; L = 50; b=4;
n = 10000;
d=4^3; % LEt number of words be d.theta is word freq.
  
alpha=ones(M,1)/M; % ALPHA
r = rand(n,1);
for i=1:n;
    if r(i) < 0.2
        Q(i,:) = [ 1 0 0 0 0];
        s(i) = 1;
    elseif r(i) < 0.4
        Q(i,:) = [ 0 1 0 0 0];
        s(i) = 2;
    elseif r(i) < 0.6
        Q(i,:) = [ 0 0 1 0 0];
        s(i) = 3;
    elseif r(i) < 0.8
        Q(i,:) = [ 0 0 0 1 0];
        s(i) = 4;
    else
        Q(i,:) = [ 0 0 0 0 1];
        s(i) = 5;
    end
end
C=randi(L,[M,d]); % C
theta=rand([M,L,b]); % theta
gc=zeros(b,L,n,M); %gc
wc=zeros(b,d,n); %wc
for i=1:n;
    for m=1:M;
        for b1=1:b;
            for l=1:L;
                if s(i) == m;
                    gc(b1,l,i,m) = theta(m,l,b1);
                    wc(b1,:,i) = gc(b1,l,i,m) /16*6*ones(1,d);
                end
            end
        end
    end
end
