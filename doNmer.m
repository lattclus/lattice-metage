clear; clc;

dirname=('../species');
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
species = {'Escherichia','Haemophilus','Pseudomonas','Aeruginosa','Azotobacter','Bordetella','Neisseria','Meningitidis','Agrobacterium','Rhodobacter','Sphaeroides','Rickettsia','Helicobacter','Campylobacter','Bacillus','Thuringiensis','Lactococcus','Streptococcus','Saccharopolyspora','Mycobacterium','Mycoplasma','Synechococcus','Synechocystis','Borrelia','Treponema','Thermus','Halobacterium','Methanocaldococcus','Pyrococcus','Sulfolobus'};
i=1;nmer = 5;
%for i=1:4,
    for j=2:4,
  %      for nmer=2:6,
            if(i~=j)
                sprintf('%s vs %s for nmer %d.......\n',species{i},species{j},nmer)
                fi = strcat('../species/',species{i});  
                fj = strcat('../species/',species{j});
                [label,model,llh,CM] = mixturecondMultinomials({fi fj},{species{i} species{j}}, nmer,L);
                sp = {species{i} species{j}};
                save(strcat('results/',[species{i} species{j}],num2str(nmer),'.mat'),'label','model','CM');
            end
   %     end
    end
%end



