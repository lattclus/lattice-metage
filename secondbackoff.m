% Naive Bayes: Mixture of Multinomials - K species
function [theta, C] = mixturecondMultinomials2(K, d, n, L, wordcount)
% fnames contains the cell array of all files to be read. All files contain
% a number of reads.
% species is a vector containing each species. Only length of species used
% in code.
% nmer is the length of the k-mers used.
% L is the number of word groups used.
% theta_b1 and theta_b2 are the two coarser/ smaller word distributions.

% label is the final predicted species to which each read belongs. 
% model contains alpha, theta(#species *#groups *4)  and c (the word group to which a word wj in species m belongs to).
% llh is the array of likelihood values.
% CM is the final confusion matrix.

%K = length(species);

%% EM algorithm:
[Q, C,groupcount] = initializationB(d,n, K, L,wordcount);
%fprintf('EM for Multinomial Conditional mixture: running ... \n');
tol = 1e-15;
maxiter = 30;
llh = -inf(1,maxiter);
converged = false;
t = 1;

while ~converged && t < maxiter
    t = t+1;
    model = maximizationB(groupcount,Q,C,wordcount);
    [Q,llh(t)] = expectationB(groupcount, model,Q);
    [xx,label] = max(Q,[],2);
	%temp = confusionmat(class,label);
	%CM(:,:,t) = temp;
    converged = llh(t)-llh(t-1) < tol*abs(llh(t));
	%acc(t)=sum(max(temp))/sum(sum(temp));
end
theta = model.theta;
C = model.C;
%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q, C,groupcount] = initializationB(d,n, init, L,wordcount)
% Initialize the mixing coefficients alpha's and species parameters lambda
    %[nb,d,n] = size(X);
if(length(init)==1) % Initialize randomly
    k = init;   % number of classes
%     idx = randsample(n,k);
%     A = reshape(X,nb*d,n);
%     m = A(:,idx); % take #species instances of reads
%     [xx,label] = max(bsxfun(@minus,m'*A,sum(m.^2,1)'/2)); % each read is assigned a species
%     while k ~= unique(label) % all species are represented.
%         idx = randsample(n,k);
%         m = A(:,idx);
%         [xx,label] = max(bsxfun(@minus,m'*A,sum(m.^2,1)'/2));
%     end
%     Q = full(sparse(1:n,label,1,n,k,n));
    tmp = randi(k,[n,1]);
    Q = zeros(n,k);
    for i=1:k
         ind = find(tmp==i);
         Q(ind,i) = 1;
    end
    C = repmat(randsample(L,d,'true')',k,1);% randomly assign groups
    while L~=length(unique(C(1,:)))
    C = repmat(randsample(L,d,'true')',k,1);% randomly assign groups
    end
    groupcount = getgroupcountB(C,L,wordcount);
elseif size(init,1) == 1 && size(init,2) == n  % initialize with labels
    label = init;
    k = max(label);
    Q = full(sparse(1:n,label,1,n,k,n));
%elseif size(init,1) == d && size(init,2) > 1  %initialize with only centers
  %  k = size(init,2);
   % m = init;
    %[xx,label] = max(bsxfun(@minus,m'*X,sum(m.^2,1)'/2));
%    Q = full(sparse(1:n,label,1,n,k,n));
else
    error('ERROR: init is not valid.');
end

%%%%%%%%%%%%% EXPECTATION STEP %%%%%%%%%%%%%%%%%%%%%%%%
function [Q,llh] = expectationB(groupcount, model,Q)
theta = model.theta;
alpha = model.alpha;
C = model.C;

n = size(groupcount,3);
k = size(theta,1); % number of species
R = zeros(n,k);

for i = 1:k
    [R(:,i)] = logcondmultinomialpdfB(squeeze(groupcount(:,:,:,i)),squeeze(theta(i,:,:))'); % returns sum over c_k AND l n_i(c_k|l)*log(theta) % input X=n_i(c_k|l)
end
R = bsxfun(@plus,R,log(alpha)); % returns sum( log(alpha*theta^n_i(w_jc_k)) )
llh = sum(sum(bsxfun(@times,Q,R))); % log likelihood using the previous Q
T = logsumexp(R,2);
Q = bsxfun(@minus,R,T);
Q = exp(Q);

%%%%%%%%%%%%%% MAXIMIZATION STEP %%%%%%%%%%%%%%%%%%%%%%%%
function  model = maximizationB(groupcount,Q,C,wordcount) % get X as an output

% function  model = maximization(X,Q,C,L)
%[b,d,n] = size(X); % b=#c_k, d=#words, n=#reads
[b,L,n,M] = size(groupcount); % b=#c_k, L=#groups, n=#reads
M = size(Q,2); % M=#species % could be left out
s = sum(Q,1); % sum(Q(:,m))
alpha = s/n;

theta = zeros(M,L,b);
for m=1:M,
    for l=1:L,
        ind = (C(m,:)==l); % in species m, ind=[ 1 0 0 1] implies words 1 and 4 belong to group l
        for k=1:b,
            tmp = squeeze((groupcount(k,:,:,m))); % fixing c_k, tmp(w_j,r_i) becomes tmp(l,r_i)
            % theta(m,l,k) = sum(tmp(ind,:),1)*Q(:,m)/(s(m)*sum(ind)); 
            theta(m,l,k) = tmp(l,:)*Q(:,m)/(s(m)*sum(ind));
		end
        theta(m,l,:) = theta(m,l,:)/sum(squeeze(theta(m,l,:)));% m,l,c_k
    end
end

[C,groupcount] = groupassignB(Q,theta,L,wordcount);
model.alpha = alpha;
model.theta = theta;
model.C = C;

