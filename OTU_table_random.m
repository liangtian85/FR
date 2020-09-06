function [ otu_table_NULL ] = OTU_table_random( otu_table, r_scheme )
%randomize otu abundance profile
%r_scheme=0 no randomization
%r_scheme=1 randomization along column (permutation)
%r_scheme=2 randomization along row (permutation)
%r_shceme=3 randomization (permutation of whole otu table)
%r_scheme=4 randomization along column (only for non-zero elements)
%r_scheme=5 randomization along row (only for non-zero elements)
%r_scheme=6 randomization (permutation of whole otu table for non-zero elements)
%r_scheme=7 uniform distribution
%r_scheme=8 uniform distribution (only for non-zero elements)
%r_scheme=9 gaussian distribution
%r_scheme=10 gaussian distribution (only for non-zero elements)

[Num_spe, Num_samp]=size(otu_table);
otu_table_NULL=zeros(size(otu_table));

switch r_scheme
    case 0
        otu_table_NULL=otu_table;
        
    case 1
        for i=1:Num_samp
            otu_table_NULL(:,i)=otu_table(randperm(Num_spe),i);
        end
        
    case 2
        for i=1:Num_spe
            otu_table_NULL(i,:)=otu_table(i,randperm(Num_samp));
        end
        
    case 3
        otu_table_vector=otu_table(:);
        otu_table_vector=otu_table_vector(randperm(length(otu_table_vector)));
        otu_table_NULL=reshape(otu_table_vector,[Num_spe,Num_samp]);
    
    case 4
        for i=1:Num_samp
            Inonzeros=find(otu_table(:,i)>0);
            otu_table_NULL(Inonzeros,i)=otu_table(Inonzeros(randperm(length(Inonzeros))),i);
        end
    
    case 5
        for i=1:Num_spe
            Inonzeros=find(otu_table(i,:)>0);
            otu_table_NULL(i,Inonzeros)=otu_table(i,Inonzeros(randperm(length(Inonzeros))));
        end
        
    case 6
        otu_table_vector=otu_table(:);
        Inonzeros=find(otu_table_vector>0);
        otu_table_vector(Inonzeros)=otu_table_vector(Inonzeros(randperm(length(Inonzeros))));
        otu_table_NULL=reshape(otu_table_vector,[Num_spe,Num_samp]);
        
    case 7
        otu_table_vector=otu_table(:);
        otu_table_vector=rand(length(otu_table_vector),1);
        otu_table_NULL=reshape(otu_table_vector,[Num_spe,Num_samp]);
        
    case 8
        otu_table_vector=otu_table(:);
        Inonzeros=find(otu_table_vector>0);
        otu_table_vector(Inonzeros)=rand(length(Inonzeros),1);
        otu_table_NULL=reshape(otu_table_vector,[Num_spe,Num_samp]);
        
    case 9
        otu_table_vector=otu_table(:);
        otu_table_vector=normrnd(0.5,0.5/3,length(otu_table_vector),1);
        otu_table_vector(otu_table_vector<0)=0;
        otu_table_vector(otu_table_vector>1)=1;
        otu_table_NULL=reshape(otu_table_vector,[Num_spe,Num_samp]);
        
    case 10
        otu_table_vector=otu_table(:);
        Inonzeros=find(otu_table_vector>0);
        otu_table_vector(Inonzeros)=normrnd(0.5,0.5/3,length(Inonzeros),1);
        otu_table_vector(otu_table_vector<0)=0;
        otu_table_vector(otu_table_vector>1)=1;
        otu_table_NULL=reshape(otu_table_vector,[Num_spe,Num_samp]);
               
end

end

