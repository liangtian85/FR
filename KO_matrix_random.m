function [ ko_matrix_null ] = KO_matrix_random ( ko_matrix, r_scheme )
%randomly rewiring the Gene Content Bipartite network
%r_scheme=0 no randomization
%r_scheme=1 complete randomization
%r_scheme=2 species-degree preserved randomization
%r_scheme=3 function-degree preserved randomization
%r_scheme=4 both species- and function- degrees preserved randomization 

[Num_spe,Num_ko]=size(ko_matrix);
ko_matrix_null=zeros(Num_spe,Num_ko);

switch r_scheme
    case 0
        ko_matrix_null=ko_matrix;
    
    case 1
        Num_link=sum(sum(ko_matrix));
        rew_spe=randi(Num_spe,Num_link,1);
        rew_ko=randi(Num_ko,Num_link,1);
        ko_matrix_null=histcounts2(rew_spe,rew_ko,1:(Num_spe+1),1:(Num_ko+1));
        
    case 2
        for i=1:Num_spe
            Num_link=sum(ko_matrix(i,:));
            rew_link=randi(Num_ko,Num_link,1);
            ko_matrix_null(i,:)=histcounts(rew_link,1:(Num_ko+1));
        end
        
    case 3
        for i=1:Num_ko
            Num_link=sum(ko_matrix(:,i));
            rew_link=randi(Num_spe,Num_link,1);
            ko_matrix_null(:,i)=histcounts(rew_link,1:(Num_spe+1));
        end
        
        
    case 4
        Num_link=sum(sum(ko_matrix));
        spe_pre_deg=sum(ko_matrix,2);
        ko_pre_deg=sum(ko_matrix,1);
        spe_deg=zeros(size(spe_pre_deg));
        ko_deg=zeros(size(ko_pre_deg));
        spe_pool=1:Num_spe;
        ko_pool=1:Num_ko;
        for i=1:Num_link
            spe_index=randi(length(spe_pool));
            spe_node=spe_pool(spe_index);
            ko_index=randi(length(ko_pool));
            ko_node=ko_pool(ko_index);
            ko_matrix_null(spe_node,ko_node)=ko_matrix_null(spe_node,ko_node)+1;
            spe_deg(spe_node)=spe_deg(spe_node)+1;
            ko_deg(ko_node)=ko_deg(ko_node)+1;
            if spe_deg(spe_node)==spe_pre_deg(spe_node)
                spe_pool(spe_index)=[];
            end
            if ko_deg(ko_node)==ko_pre_deg(ko_node)
                ko_pool(ko_index)=[];
            end
        end
          
end%r_scheme

end

