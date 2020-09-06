function [ SD, FD, FR ] = FDFR_Rao_q ( otu_table, dij, q )
% SD = \sum_i \sum_j (p_i)^q (p_j)^q 
% FD = \sum_i \sum_j d_ij (p_i)^q (p_j)^q 
% FR = SD(q) - FD(q)

[Num_spe, Num_samp] = size (otu_table);

otu_table=otu_table./repmat(sum(otu_table),Num_spe,1);

dij=squareform(dij);

for i=1:Num_samp
    otu_vector=otu_table(:,i);
    otu_matrix=otu_vector*otu_vector';
    otu_matrix(1:Num_spe+1:Num_spe*Num_spe)=0;
    if q==0
        SD(i)=sum(sum(otu_matrix>0));
        FD(i)=sum(sum((otu_matrix>0).*dij));
        FR(i)=SD(i)-FD(i);
        continue;
    end
    SD(i)=sum(sum(otu_matrix.^q));
    FD(i)=sum(sum(otu_matrix.^q.*dij));
    FR(i)=SD(i)-FD(i);
end

end

