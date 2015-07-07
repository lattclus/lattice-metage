% Naive Bayes: Mixture of Multinomials - K species
function [label,model,llh] = mixturecondMultinomials(species, d, n, L, wordcountfile, theta_back, C_back, theta2, L2)
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

K = length(species);

%% EM algorithm:
[Q, C,groupcount] = initialization(d,n, K, L,wordcountfile);
fprintf('EM for Multinomial Conditional mixture: running ... \n');
tol = 1e-30;
maxiter = 50;
llh = -inf(1,maxiter);
converged = false;
t = 1;

while ~converged && t < maxiter
    t = t+1;
    model = maximization(groupcount,Q,C,wordcountfile);
    [Q,llh(t)] = expectation(groupcount, model,Q);
    [xx,label] = max(Q,[],2);
	converged = llh(t)-llh(t-1) < tol*abs(llh(t));
end

%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q, C,groupcount] = initialization(d,n, init, L,wordcountfile)
% Initialize the mixing coefficients alpha's and species parameters lambda
if(length(init)==1) % Initialize randomly
    k = init;   % number of classes
    tmp = randi(k,[n,1]);
    Q = zeros(n,k);
    for i=1:k
         ind = find(tmp==i);
         Q(ind,i) = 1;
    end
    C = repmat(randsample(L,d,'true')',k,1);% randomly assign groups
    groupcount = getgroupcount(C,L,wordcountfile);
elseif size(init,1) == 1 && size(init,2) == n  % initialize with labels
    label = init;
    k = max(label);
    Q = full(sparse(1:n,label,1,n,k,n));
else
    error('ERROR: init is not valid.');
end

%%%%%%%%%%%%% EXPECTATION STEP %%%%%%%%%%%%%%%%%%%%%%%%
function [Q,llh] = expectation(groupcount, model,Q)
% using groupcount updated after getting C
theta = model.theta;
alpha = model.alpha;
C = model.C;

n = size(groupcount,3);
k = size(theta,1); % number of species
R = zeros(n,k);

for i = 1:k
    [R(:,i)] = logcondmultinomialpdf(squeeze(groupcount(:,:,:,i)),squeeze(theta(i,:,:))'); % returns sum over c_k AND l n_i(c_k|l)*log(theta) % input X=n_i(c_k|l)
end
R = bsxfun(@plus,R,log(alpha)); % returns sum( log(alpha*theta^n_i(w_jc_k)) )
llh = sum(sum(bsxfun(@times,Q,R))); % log likelihood using the previous Q
T = logsumexp(R,2);
Q = bsxfun(@minus,R,T);
Q = exp(Q);

%%%%%%%%%%%%%% MAXIMIZATION STEP %%%%%%%%%%%%%%%%%%%%%%%%
function  model = maximization(groupcount,Q,C,wordcountfile) % get X as an output
% function  model = maximization(X,Q,C,L)
[b,L,n,M] = size(groupcount); % b=#c_k, L=#groups, n=#reads
M = size(Q,2); % M=#species % could be left out
s = sum(Q,1); % sum(Q(:,m))
alpha = s/n;
tau = 10; 

theta = zeros(M,L,b);
I=zeros(m,l,k); 
for m=1:M,
    for l=1:L,
        ind = (C(m,:)==l); % in species m, ind=[ 1 0 0 1] implies words 1 and 4 belong to group l
        for k=1:b,
            tmp = squeeze((groupcount(k,:,:,m))); % fixing c_k, tmp(w_j,r_i) becomes tmp(l,r_i)
			count(m,l,k) = tmp(l,:)*Q(:,m);
            theta(m,l,k) = count(m,l,k)/(s(m)*sum(ind));
			if count(m,l,k) >= tau
				I(m,l,k) = 1; % This one will be used for the computation of lambda (l1 and l2).
			end
		end
        theta(m,l,:) = theta(m,l,:)/sum(squeeze(theta(m,l,:)));% m,l,c_k
    end
end
% Equation for lambda:
for m=1:M
	for l=1:L
		% theta1(m,l) is incorrect
		temp = squeeze(theta1(m,l,:)) + squeeze(theta2(m,l,:));
		temp1 = ( squeeze(theta1(m,l,:))*.squeeze(I(m,l,:)) )./temp;
		temp2 = ( squeeze(theta2(m,l,:))*.squeeze(I(m,l,:)) )./temp;
		a1 = ( temp1' * squeeze(groupcount(:,l,:,m)) ) * Q(:,m);
		a2 = squeeze(theta1(m,l,:))'*squeeze(I(m,l,:));
		a3 = ( squeeze(I(m,l,:))'*squeeze(groupcount(:,l,:,m)) ) * Q(:,m);
		a4 = 1 - squeeze(theta(m,l,:))'*(1-squeeze(I(m,l,:)));
		a5 = ( temp2' * squeeze(groupcount(:,l,:,m)) ) * Q(:,m);
		a6 = squeeze(theta2(m,l,:))'*squeeze(I(m,l,:));
		l1(m,l) = a1/(a2*a3/a4);
		l2(m,l) = a5/(a6*a3/a4);
		lam1(m,l) = l1(m,l) / ( l1(m,l) + l2(m,l) );
		lam2(m,l) = l1(m,l) / ( l1(m,l) + l2(m,l) );
		%let's clear stuff here
	end
end
% Equation for delta:
for m=1:M
	for l=1:L
		s1 = squeeze(theta(m,l,:))'*(1-squeeze(I(m,l,:)));
		s2 = lam1*(squeeze(theta1(m,l,:))'*squeeze(I(m,l,:))) + lam2*(squeeze(theta2(m,l,:))'*squeeze(I(m,l,:)));
		delta(m,l) = (1-s1)/s2;
	end
end
[C,groupcount] = groupassign(Q,theta,M,L,wordcountfile);  % update groupcount after getting C
% Backing off to a smaller word model: finding the parameter to be used.
% NEW I: using new C, to be used for Qs next.

for m=1:M,
    for l=1:L,
        ind = (C(m,:)==l); % in species m, ind=[ 1 0 0 1] implies words 1 and 4 belong to group l
        for k=1:b,
			for j=1:length(ind)
				if ( ind(j) >= 1+(k-1)*b^(d-1) && ind(j) <= b^(d-1) ) 
					ind2(j) = ind(j) - i*b^(d-1); % word index of smaller size words
				end
			end
            tmp = squeeze((groupcount(k,:,:,m))); % fixing c_k, tmp(w_j,r_i) becomes tmp(l,r_i)
            count(m,l,k) = tmp(l,:)*Q(:,m);
			if count(m,l,k) >= tau
				I(m,l,k) = 1;
			end
		end
		theta1(m,l,:) = theta_back( m, mode(C_back(m,ind2)), : ); % what group do these smaller words belong to?
    end
end

theta2=secondbackoff( species, L, n, L2, groupcount,wordcountfile); % total number of words is L
model.alpha = alpha;
model.theta = theta;
model.C = C;

% Ask for new EM results for the backoff models. Give inputs n_{im} {c_k|w_j}
% smaller words model will be found outside, and the members of word group will be polled to see which group in the smaller word model lies closest.
% for word groups, n_im(c_k|l) will be used.lc_k is the new word.
