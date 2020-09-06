clear
clc
% probability per step/HGT between species
% Event genome modication occurs 
%      subevent I (HGT): with probability q_hgt
%      subevent II (gene loss): with probability q_gl
%      subevent III (gene gain): with probablity q_gg
g=genpath('./bimat');
addpath(g);
p_gg=0.005;
n_step=500000;
p_hgt=0.795;

%initial network
rng(0)
ns0=500;
ng0=200;
GCN=zeros(ns0,ng0);
GCN(1:ns0,1:ng0)=rand(ns0,ng0)<0.8;

for k=0:n_step
    rng(k)
    probability=rand();
    %% randomly select one species
    [ns,ng]=size(GCN);
    species_degree=sum(GCN,2);
    if sum(species_degree)>0
        pro_spec=species_degree.^2/sum(species_degree.^2);
        pro_spec=[0;cumsum(pro_spec)];
        n1=find(rand<pro_spec);
        index_spec=n1(1)-1;
    else
        index_spec=datasample(1:ns,1);
    end
    %% Event I: gene loss
    p_gl=1-p_hgt-p_gg;
    if probability<p_gl
        m1=find(GCN(index_spec,:)==1);
        if ~isempty(m1)
            gene_index=datasample(m1,1);
            GCN(index_spec,gene_index)=0;
        end
        %% event II gene gain
    elseif probability>=p_gl&&probability<(p_gl+p_gg)
        new_gene=zeros(size(GCN,1),1);
        new_gene(index_spec)=1;
        GCN=[GCN,new_gene];
    else
        %% event III: HGT
        donor_spec=randi(ns-1);
        if donor_spec==index_spec
            donor_spec=donor_spec+1;
        end
        HGT_gene_pool=GCN(donor_spec,:)-GCN(index_spec,:);
        HGT_gene_slot=find(HGT_gene_pool==1);% gene pool that can be uptaken by this species
        if ~isempty(HGT_gene_slot)
            HGT_gene=datasample(HGT_gene_slot,1);
            GCN(index_spec,HGT_gene)=1;
        end
    end
    %% remove null gene
    m1=find(sum(GCN,1)==0);
    GCN(:,m1(1:end))=[];
    %% plot the results at different steps
    if mod(k,100000)==0
        dis=pdist(GCN,'jaccard');
        %plot GCN
        figure
        subplot(2,2,1)
        bp=Bipartite(logical(GCN));
        bp.ntc.CalculateNestedness();
        NTC=bp.ntc.N;
        NODF=bp.nodf.nodf;
        
        imagesc(1:ng,1:ns,bp.ntc.MatrixMinimal);
        set(gca,'ylim',[0,size(GCN,1)]);
        set(gca,'YDir','reverse');
        set(gca,'xlim',[0,size(GCN,2)]);
        set(gca,'fontsize',8);
        box on;
        axis square;
        
        %plot distance
        d1=sum(GCN,1);
        d2=sum(GCN,2);
        
        subplot(2,2,2)
        edge=0:0.02:1;
        x=edge(1:end-1)+0.01;
        prob = histcounts(dis,edge, 'Normalization', 'probability');
        plot(x,prob,'color',[1.0,0.0,0.0],'LineWidth',2);
        set(gca,'fontsize',8);
        set(gca,'xtick',[0,0.2,0.4,0.6,0.8,1]);
        set(gca,'TickLength',[0.03,0.03]);
        box on;
        axis square;
        
        %plot degree of species
        subplot(2,2,3)
        h=histogram(d2,25,'Normalization', 'probability');
        h.EdgeColor = [1,1,1];
        box on;
        axis square;
        set(gca,'fontsize',8);
        set(gca,'TickLength',[0.03,0.03]);
        
        %plot degree of genes
        subplot(2,2,4)
        edge=sort(unique(d1));
        edge=[edge,edge(end)+1];
        Fun_degdis=histcounts(d1,edge);
        Fun_degdis=Fun_degdis/sum(Fun_degdis);
        plot(edge(1:end-1),Fun_degdis,'o','color',[0,0,1],'MarkerSize',6);
        
        box on;
        axis square;
        set(gca,'XScale','log');
        set(gca,'YScale','log');
        set(gca,'fontsize',8);
        set(gca,'TickLength',[0.03,0.03]);
    end
    GCN=sparse(GCN);
end
