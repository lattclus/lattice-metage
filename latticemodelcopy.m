% For back off, combine results from all words together. lambda(m,l,j)
% Naive Bayes: Mixture of Multinomials - K species
function [Q,label,model,llh] = latticemodelcopy(M, nmer, nmer2, n, L, wordcountfile,L2,theta_back1, C_back1,tau,fill3)
% THRESHOLD KEPT AT TAU PERCENT
% fnames contains the cell array of all files to be read. All files contain
% a number of reads.
% species is a vector containing each species. Only length of species used
% in code.
% nmer is the length of the k-mers used.
% L is the number of word groups used.
% theta_b1 and theta_b2 are the two coarser/ smaller word distributions.

% label is the final predicted species to which each read belongs. 
% model contains alpha, theta(#species *#groups *4)  and c (the word group 
% to which a word wj in species m belongs to).
% llh is the array of likelihood values.
% CM is the final confusion matrix.

%% EM algorithm:
d = 4^nmer;
flag_bo = 1; thresh = 0;
[Q, C,groupcount] = initialization(d,n, M, L,wordcountfile,tau,fill3);
fprintf('EM for Multinomial Conditional mixture: running ... \n');
tol = 1e-15;
maxiter = 10;
llh = -inf(1,maxiter);
converged = false;
t = 1;

while ~converged && t < maxiter
    t = t+1;    
    [model,thresh] = maximization(groupcount,Q,C,wordcountfile,nmer, nmer2, L2,theta_back1, C_back1,tau,thresh,flag_bo);
    fprintf('E for Multinomial Conditional mixture: running ... \n');
    flag_bo = 1;
    groupcount = model.gc;
    [Q,llh(t)] = expectation(model,Q);
	if (llh(t)-llh(t-1) > 0 && t~=2)
		[xx,label] = max(Q,[],2);
		clear xx;
		
	else
		fprintf('convergence issue ... \n');		
	end
	converged = llh(t)-llh(t-1) < tol*abs(llh(t));
end

%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q, C,groupcount] = initialization(d,n, M, L,wordcountfile, tau,fill3)
% Initialize the mixing coefficients alpha's and species parameters lambda
if(length(M)==1) % Initialize randomly
    k = M;   % number of classes
    tmp2 = fill3; %randi(k,[n,1]);
    Q = zeros(n,k);
    for i=1:k
         ind = find(tmp2==i);
         Q(ind,i) = 1; % for m=1:M; Q(10000*(m-1)+1:10000*m,m) = 1;end
    end
    C = repmat(randsample(L,d,'true')',k,1);% randomly assign groups
    groupcount = getgroupcount(C,L,wordcountfile);
%    [b,L,n,M] = size(groupcount); % b=#c_k, L=#groups, n=#reads, M=#species
% 
%     for m=1:M,
%         for l=1:L
%             for k=1:b
%                 tmp = squeeze((groupcount(k,:,:,m))); % fixing c_k, tmp(w_j,r_i) becomes tmp(l,r_i)
%                 count(m,l,k) = tmp(l,:)*Q(:,m);
%             end
%         end
% 	end
%             
%     thresh = min(min(min(count)));
% 	thresh = thresh + tau*( max(max(max(count))) - thresh )
elseif size(M,1) == 1 && size(M,2) == n  % initialize with labels
    label = M;
    k = max(label);
    Q = full(sparse(1:n,label,1,n,k,n));
else
    error('ERROR: init is not valid.');
end

%%%%%%%%%%%%% EXPECTATION STEP %%%%%%%%%%%%%%%%%%%%%%%%
function [Q,llh] = expectation(model,Q)
% using groupcount updated after getting C
%theta = model.theta;
alpha = model.alpha;
%theta_back = model.theta_back;
%C = model.C;
%I = model.I;
groupcount = model.gc;
thetaf= model.thetaf;
n = size(groupcount,3);
M = size(thetaf,1); % number of species
R = zeros(n,M);

for m = 1:M
    %[R(:,m)] = logcondmultinomialpdf(squeeze(groupcount(:,:,:,m)),squeeze(theta(m,:,:))',squeeze(I(m,:,:)), squeeze(theta_back(m,:,:))' ); % returns sum over c_k AND l n_i(c_k|l)*log(theta) % input X=n_i(c_k|l)
    [R(:,m)] = logcondmultinomialpdfo(squeeze(groupcount(:,:,:,m)),squeeze(thetaf(m,:,:))' ); % returns sum over c_k AND l n_i(c_k|l)*log(theta) % input X=n_i(c_k|l)

end
R = bsxfun(@plus,R,log(alpha)); % returns sum( log(alpha*theta^n_i(w_jc_k)) )
llh = sum(sum(bsxfun(@times,Q,R))); % log likelihood using the previous Q
T = logsumexp(R,2);
Q = bsxfun(@minus,R,T); % normalization
Q = exp(Q);


%%%%%%%%%%%%%% MAXIMIZATION STEP %%%%%%%%%%%%%%%%%%%%%%%%
function  [model,thresh] = maximization(groupcount,Q,C,wordcountfile,nmer,nmer2,L2,theta_back1, C_back1,tau,thresh,flag_bo) % get X as an output
[b,L,n,M] = size(groupcount); % b=#c_k, L=#groups, n=#reads, M=#species
if flag_bo==1;
    for m=1:M,
        for l=1:L
            for k=1:b
                tmp = squeeze((groupcount(k,:,:,m))); % fixing c_k, tmp(w_j,r_i) becomes tmp(l,r_i)
                count(m,l,k) = tmp(l,:)*Q(:,m);
            end
        end
    end
    clear tmp           
    thresh = min(min(min(count)));
    thresh = thresh + tau*( max(max(max(count))) - thresh )

end
s = sum(Q,1); % sum(Q(:,m))
alpha = s/n;
[theta_back2,C_back2] = secondbackoff( M, L, n, L2, groupcount); % total number of words is L
theta_back2(isnan(theta_back2)) = 0; % ISNAN TO 0
theta = zeros(M,L,b);
I = zeros(M,L,b); 
%load(wordcountfile);
for m=1:M,
    for l=1:L
        ind = (C(m,:)==l); % in species m, ind=[ 1 0 0 1] implies words 1 and 4 belong to group l
        ind1 = find(C(m,:)==l); % collect all word index in the current group
        ind2 =zeros(1,length(ind1));
            for k=1:b^(nmer-nmer2),
			    for j=1:length(ind1)
				    if ( ind1(j) >= 1+(k-1)*4^(nmer2) && ind1(j) <= k*4^(nmer2) ) 
					    ind2(j) = ind1(j) - (k-1)*4^(nmer2); % word index of smaller size words
                        passonind{m,l}(1,j) = ind2(j);
                    end
                end
            end
		    for k=1:b,
                tmp = squeeze((groupcount(k,:,:,m))); % fixing c_k, tmp(w_j,r_i) becomes tmp(l,r_i)
                count(m,l,k) = tmp(l,:)*Q(:,m);
				x = count(m,l,k)/(s(m)*sum(ind));
				x(isnan(x)) = 0; % ISNAN TO 0
                xx(k) = x;
				if count(m,l,k) <= thresh
					I(m,l,k) = 1;
				end
            end
        xxx = theta_back2(m, C_back2(m,l), :);
        xx = xx/sum(xx);% m,l,c_k
        xx(isnan(xx)) = 0; % ISNAN TO 0
        theta(m,l,:) = xx;
        xxx(isnan(xxx)) = 0; % ISNAN TO 0
        theta2(m,l,:) = xxx;
        clear xx xxx;
% b=#c_k, d=#words, n=#reads
% Equation for lambda:
%lam1 = ones(M,L); lam2 = ones(M,L);
% for m=1:M
% 	for l=1:L
       
        for j=1:length(ind2)
            tb1 = squeeze(theta_back1(m, C_back1(m,ind2(j)),:)); % theta1
            temp = tb1 + squeeze(theta2(m,l,:)); % theta1+theta2
            temp1 = ( tb1.*squeeze(I(m,l,:)) )./temp; % theta1*I/(theta1+theta2)
            temp2 = ( squeeze(theta2(m,l,:)).*squeeze(I(m,l,:)) )./temp; % theta2*I/(theta1+theta2)
            %a1 = ( temp1' * squeeze(X(:,j,:)) ) * Q(:,m);%c_k,L,i,M % numerator of lam1
			a1 = ( temp1' * squeeze(groupcount(:,l,:,m)) ) * Q(:,m);%c_k,L,i,M % numerator of lam1
            a2 = tb1'*squeeze(I(m,l,:)); % sum_c I*theta1
            %a3 = ( squeeze(I(m,l,:))'*squeeze(X(:,j,:)) ) * Q(:,m);
			a3 = ( squeeze(I(m,l,:))'*squeeze(groupcount(:,l,:,m)) ) * Q(:,m);
            a4 = 1 - squeeze(theta(m,l,:))'*(1-squeeze(I(m,l,:))); % 1-\sum_c(1-I)*theta
            %a5 = ( temp2' * squeeze(X(:,j,:)) ) * Q(:,m);
			a5 = ( temp2' * squeeze(groupcount(:,l,:,m)) ) * Q(:,m); % numerator of lam2
            a6 = squeeze(theta2(m,l,:))'*squeeze(I(m,l,:)); % sum_c I*theta2
            %l1(m,l) = a1/(a2*a3/a4);
            %l2(m,l) = a5/(a6*a3/a4);
            if a1*a6 + a2*a5 ~= 0
                lam1{m,l,j} = (a1*a6)/(a1*a6+a2*a5);%l1(m,l) / ( l1(m,l) + l2(m,l) );
                lam2{m,l,j} = (a5*a2)/(a5*a2+a6*a1);%l1(m,l) / ( l1(m,l) + l2(m,l) );
            else
                lam1{m,l,j} = 0;
                lam2{m,l,j} = 0;
            end
            %lam1{m,l,j} = 0;
            %lam2{m,l,j} = 1;
        end
		clear a1 a2 a3 a4 a5 a6 temp temp1 temp2 %let's clear stuff here
	end
end
%clear X;
% Equation for delta:
% lam1(cellfun(@isempty,lam1)) = {0}; lam2(cellfun(@isempty,lam2)) = {0};
delta = zeros(M,L);
s3 = zeros(M,L,b);
for m=1:M
	for l=1:L
		s1 = squeeze(theta(m,l,:))'*(1-squeeze(I(m,l,:)));
        s2 = 0;
        
        for j=1:length(find(C(m,:)==l))
            tb1 = squeeze( theta_back1(m, C_back1(m,passonind{m,l}(1,j)),:) );
            s2 = s2 + lam1{m,l,j}*(tb1'*squeeze(I(m,l,:))) + lam2{m,l,j}*(squeeze(theta2(m,l,:))'*squeeze(I(m,l,:)));
            % theta1 is {m,l,j}(c_k)
            s3(m,l,:) = squeeze(s3(m,l,:)) + lam1{m,l,j}*tb1 + lam2{m,l,j}*squeeze(theta2(m,l,:));
        end
        if s2~=0
            delta(m,l) = (1-s1)/s2;
        end
	end
end
delta(isnan(delta))=0;
for m=1:M
    for l=1:L
        theta_back(m,l,:) = delta(m,l)*squeeze(s3(m,l,:));
        
    end
end
thetaf = theta;
for b1= 1:b;
    [i1, i2] = find(I(:,:,b1)==1);
    for ic=1:length(i1);
        thetaf(i1(ic),i2(ic),b1) = theta_back(i1(ic),i2(ic),b1);
    end
end

 % replacing zeros with ones.
%[C,groupcount] = groupassign(Q,theta,I,theta_back,M,L,wordcountfile);  % update groupcount after getting C
[C,groupcount] = groupassigno(Q,thetaf,M,L,wordcountfile);  % update groupcount after getting C

   % model.thresh=thresh;
model.alpha = alpha;
model.theta = theta;
model.thetaf = thetaf;
model.C = C;
model.I = I;
model.theta_back = theta_back;
model.gc = groupcount;
model.lam = lam1;
%% Combining the results for theta


% Ask for new EM results for the backoff models. Give inputs n_{im} {c_k|w_j}
% smaller words model will be found outside, and the members of word group will be polled to see which group in the smaller word model lies closest.
% for word groups, n_im(c_k|l) will be used.lc_k is the new word.
