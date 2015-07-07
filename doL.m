clear; clc;

species = {'Escherichia','Haemophilus','Pseudomonas','Aeruginosa','Azotobacter','Bordetella','Neisseria','Meningitidis','Agrobacterium','Rhodobacter','Sphaeroides','Rickettsia','Helicobacter','Campylobacter','Bacillus','Thuringiensis','Lactococcus','Streptococcus','Saccharopolyspora','Mycobacterium','Mycoplasma','Synechococcus','Synechocystis','Borrelia','Treponema','Thermus','Halobacterium','Methanocaldococcus','Pyrococcus','Sulfolobus'};

nmer = 5;
K = length(species);
L = [10, 20, 30, 50, 70, 100];
for i=3:4,
    for j=1:4,
        for l=1:length(L)
            if(i~=j)
                sprintf('%s vs %s .......\n',species{i},species{j})
                fi = strcat('../species/',species{i});  
                fj = strcat('../species/',species{j});
                [label,model,llh,CM] = mixturecondMultinomials({fi fj},{species{i} species{j}}, nmer, L(l));
                sp = {species{i} species{j}};
                save(strcat('results/',[species{i} species{j}],num2str(L(l)),'.mat'),'label','model','CM','llh');
            end
        end
    end
end


