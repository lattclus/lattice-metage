function condfreq = condfreqComp(reads, nmer)

%==== Determine the frequency of various nmers in reads{i}
% 'reads': sequences of reads (nreads of cell type)
% 'nmer': Size of nmers(trimer, tetramer)
% 'freq': Returns a 4 X 4^nmer_len sized matrix * nreads 

%--------------------------------------------------------------------------

nreads = max(size(reads));
% nmers
nmers = pick('ACTG',nmer,'or');
base = ['A','C','T','G'];

ncol = 4^nmer;
condfreq = zeros(4,ncol,nreads);
for i=1:nreads, % index for reads
    tmpseq = reads{i};
    len = length(tmpseq);
    for j=1:ncol,   % index for parents pa_j
        ind = findstr(tmpseq,nmers(j,:));
        for l=1:length(ind),
            for k=1:4,   % index for c_k
                if(ind(l)+nmer<= len && (tmpseq(ind(l)+nmer) == base(k)))
                    condfreq(k,j,i) = condfreq(k,j,i)+1;
                end
            end
        end
    end
end
   