clear; clc;

dirname=('../species');
nmer = 4;
L=10;

fnames = dir(dirname);
k=1;
for i=1:length(fnames),
    if(~strcmp(fnames(i).name,'.') && ~strcmp(fnames(i).name,'..'))
        species{k} = fnames(i).name;
        k = k+1;
    end
end

K = length(species);
k=1;
for i=6:10,
    for j=1:K,
        if(i~=j)
            sprintf('%s vs %s .......\n',species{i},species{j})
            fi = strcat('../species/',species{i});  
            fj = strcat('../species/',species{j});
            [label,model,llh,CM] = mixturecondMultinomials({fi fj},{species{i} species{j}}, nmer, L);
            sp = {species{i} species{j}};
            save(strcat([species{i} species{j}],'.mat'),'label','model','CM');
            k = k+1;
        end
    end
end


