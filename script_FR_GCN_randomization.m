clear;

%% load data
load HMP_Stool.mat;% taxonomic abundance and gcn

ab_table_real=ab_table_Stool_real;
ko_table_real=ko_table_Stool_real;

%% Taxonomic diversity, Functional diversity, and Functional redundancy calculation
% filtering
ab_table_real=ab_table_real(:,sum(ab_table_real>0,1)>5);

% Real taxonomic profile
dij_real=pdist(ko_table_real,@distfun_WeightedJaccard);
[SD_real, FD_real, FR_real]=FDFR_Rao_q (ab_table_real, dij_real);

% Null-GCN-1
ko_table_NULL_1=KO_matrix_random(ko_table_real,1);
dij_NULL_1=pdist(ko_table_NULL_1,@distfun_WeightedJaccard);
[SD_NULL_1, FD_NULL_1, FR_NULL_1]=FDFR_Rao_q (ab_table_real, dij_NULL_1);

% Null-GCN-2
ko_table_NULL_2=KO_matrix_random(ko_table_real,2);
dij_NULL_2=pdist(ko_table_NULL_2,@distfun_WeightedJaccard);
[SD_NULL_2, FD_NULL_2, FR_NULL_2]=FDFR_Rao_q (ab_table_real, dij_NULL_2);

% Null-GCN-3
ko_table_NULL_3=KO_matrix_random(ko_table_real,3);
dij_NULL_3=pdist(ko_table_NULL_3,@distfun_WeightedJaccard);
[SD_NULL_3, FD_NULL_3, FR_NULL_3]=FDFR_Rao_q (ab_table_real, dij_NULL_3);

% Null-GCN-4
ko_table_NULL_4=KO_matrix_random(ko_table_real,4);
dij_NULL_4=pdist(ko_table_NULL_4,@distfun_WeightedJaccard);
[SD_NULL_4, FD_NULL_4, FR_NULL_4]=FDFR_Rao_q (ab_table_real, dij_NULL_4);


%% Mann-Whitney U test
p_values_1=signrank(FR_real./SD_real,FR_NULL_1./SD_NULL_1);
p_values_2=signrank(FR_real./SD_real,FR_NULL_2./SD_NULL_2);
p_values_3=signrank(FR_real./SD_real,FR_NULL_3./SD_NULL_3);
p_values_4=signrank(FR_real./SD_real,FR_NULL_4./SD_NULL_4);

q_values=mafdr([p_values_1,p_values_2,p_values_3,p_values_4],'BHFDR', true);

%% Figure
real_posi=1;
NULL_1_posi=2;
NULL_2_posi=3;
NULL_3_posi=4;
NULL_4_posi=5;

FR_real=[FR_real./SD_real];
g_FR_real=[ones(1,length(FR_real))*real_posi(1)];

FR_NULL_1=[FR_NULL_1./SD_NULL_1];
g_FR_NULL_1=[ones(1,length(FR_NULL_1))*NULL_1_posi(1)];

FR_NULL_2=[FR_NULL_2./SD_NULL_2];
g_FR_NULL_2=[ones(1,length(FR_NULL_2))*NULL_2_posi(1)];

FR_NULL_3=[FR_NULL_3./SD_NULL_3];
g_FR_NULL_3=[ones(1,length(FR_NULL_3))*NULL_3_posi(1)];

FR_NULL_4=[FR_NULL_4./SD_NULL_4];
g_FR_NULL_4=[ones(1,length(FR_NULL_4))*NULL_4_posi(1)];

real_color=[59/255,59/255,59/255];
NULL_1_color=[0,0.45,0.74];
NULL_2_color=[0.85,0.33,0.1];
NULL_3_color=[0.47,0.67,0.19];
NULL_4_color=[0.93,0.69,0.13];

figure('position',[537 713 977/5*3/4 420*2/3]);
hold on;
boxplot(FR_real,g_FR_real,'color',real_color,'positions',real_posi,'width',0.5,'Symbol','.','OutlierSize',0.5);
boxplot(FR_NULL_1,g_FR_NULL_1,'color',NULL_1_color,'positions',NULL_1_posi,'width',0.5,'Symbol','.','OutlierSize',0.5);
boxplot(FR_NULL_2,g_FR_NULL_2,'color',NULL_2_color,'positions',NULL_2_posi,'width',0.5,'Symbol','.','OutlierSize',0.5);
boxplot(FR_NULL_3,g_FR_NULL_3,'color',NULL_3_color,'positions',NULL_3_posi,'width',0.5,'Symbol','.','OutlierSize',0.5);
boxplot(FR_NULL_4,g_FR_NULL_4,'color',NULL_4_color,'positions',NULL_4_posi,'width',0.5,'Symbol','.','OutlierSize',0.5);

set(gca,'fontsize',10)
set(gca,'XTickLabel',{' '});
set(gca,'ylim',[0,0.8]);
set(gca,'xlim',[-2,8]);
set(gca,'xtick',[]);









