opts = detectImportOptions('bipartite_trait_gene.csv');
T = readtable('bipartite_trait_gene.csv',opts);

uTraits = unique(T(:,2)); 
uTraits=uTraits.Trait;
ugenes = unique(T(:,3)); 
ugenes=ugenes.Significant_Genes;

%read in data and put in a bipartite Trait x Gene matrix
TraitGeneMatrix=zeros(size(uTraits,1), size(ugenes,1));
for ii= 1: size(uTraits,1)
    t = uTraits{ii};
    g = T.Significant_Genes(find(strcmp(t, T.Trait(:,1))));
    g2=zeros(size(g,1),1);
    for jj = 1:size(g,1)
        g2(jj) = find(strcmp(g{jj}, ugenes));
    end
    TraitGeneMatrix(ii,g2)=1;
end

GGmatrix=TraitGeneMatrix' * TraitGeneMatrix;
GGmatrix=GGmatrix-diag(diag(GGmatrix));


%lets do the shuffle 100 times(music time)
for i = 1:100
    [R, eff1] = randmio_und_bipartite(TraitGeneMatrix, 10);
    [R, eff2] = randmio_und_bipartite(R', 10);
    Random_GGmatrix=R';
    Random_GGmatrix_table_{i} = array2table(Random_GGmatrix,'RowNames',uTraits,'VariableNames',ugenes);
    filename = ['random_GTM_', num2str(i), '.csv'];
    writetable(Random_GGmatrix_table_{i}, filename, 'WriteRowNames', true);
end




