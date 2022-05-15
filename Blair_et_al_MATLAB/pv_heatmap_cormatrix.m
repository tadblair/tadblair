function [mapR, mapP] = pv_heatmap_cormatrix(traincurves,testcurves)

if size(traincurves,1)>1
    
        for i=1:size(traincurves,2)
            for j=1:size(testcurves,2)
                temp1=squeeze(traincurves(:,i));
                temp2=squeeze(testcurves(:,j));
                good=find(~isnan(temp1) & ~isnan(temp2));
                [R,P] = corrcoef(temp1(good),temp2(good));
                middleR(i,j)=R(1,2); 
                middleP(i,j)=P(1,2);
                mapR(i,j)=R(1,2); 
                mapP(i,j)=-log10(P(1,2))*sign(R(1,2));
                %mapR(i,j)=1-pdist2(temp1(good)',temp2(good)','cosine');
            end
        end      
    
else
    
    mapR=NaN(23);
    mapP=NaN(23);
        
end