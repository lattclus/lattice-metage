function gc = getgroupcountB(C,L,wordcount)
% groupcount gc is (c_k,l,i,m)
%load(wordcountfile);
[M,wj] = size(C);
[b,wj,n,M]=size(wordcount);
gc=zeros(b,L,n,M);
for m=1:M
	%for k=1:b
		%for i=1:n
			for l=1:L
				ind=find(C(m,:)==l);
				gc(:,l,:,m)=sum(wordcount(:,ind,:,m),2); % ??
			end
		%end
	%end
end