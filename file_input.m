fnames = {'Agrobacterium', 'Bordetella', 'Campylobacter', 'Escherichia', 'Haemophilus'}
nmer = 5
nmer2 = 3
L = 400
L2 = 100
M = 5
%%
if(length(fnames) ~= M)
   display('Error');
end
[header,sequence]= fastaread(fnames{1});
for i=2:length(fnames)
    [tmpheader,tmpsequence]= fastaread(fnames{i});
    sequence = [sequence(:); tmpsequence(:)]; % all reads collected
end
clear tmpheader tmpsequence header;
X = condfreqComp(sequence,nmer); % 4 X 4^nmer_len X #reads sized matrix (number of words w_j^c_k in read i )
wordcountfile=sprintf('wordcount%d.mat',nmer);
save(wordcountfile,'X');
clear sequence X;

%%
mfile=sprintf('model%d_L%d.mat',nmer2,L);
load(mfile);
n= M*10000; % #reads
tau = .01
%%
[Q,label,model,llh] = latticemodelcopy( M, nmer, nmer2, n, L, wordcountfile,L2,model.theta, model.C,tau );
	