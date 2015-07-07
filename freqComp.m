function freq = freqComp(reads, nmer_len)

%==== Determine the frequency of various nmers in reads{i}
% 'reads': sequences of reads (nreads of cell type)
% 'nmer_len': Size of nmers(trimer, tetramer)
% 'freq': Returns a nreads X 4^nmer_len sized matrix 

%--------------------------------------------------------------------------

nreads = max(size(reads));
% nmers
nmers = pick('ACTG',nmer_len,'or');

ncol = 4^nmer_len;
freq = zeros(nreads,ncol);
for i=1:nreads,
    for j=1:ncol,
        freq(i,j) = size(findstr(reads{i},nmers(j,1:nmer_len)),2);
    end
end


    