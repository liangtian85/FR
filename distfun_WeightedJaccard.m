function [ d2 ] = distfun_WeightedJaccard ( XI, XJ )
% Distance function for weighted Jaccard distance
%
%XI=repmat(XI,size(XJ,1),1);
%

d2=sum(min(XI,XJ),2)./sum(max(XI,XJ),2);
d2(isnan(d2))=0;
d2=1-d2;

end

