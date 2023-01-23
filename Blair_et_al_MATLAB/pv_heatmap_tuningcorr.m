function [resultR, resultP] = pv_heatmap_tuningcorr(traincurves,testcurves)

        for i=1:size(traincurves,1)
                temp1=squeeze(traincurves(i,:));
                temp2=squeeze(testcurves(i,:));
                good=find(~isnan(temp1) & ~isnan(temp2));
                [R,P] = corrcoef(temp1(good),temp2(good));
                resultR(i)=R(1,2); 
                %resultP(i)=P(1,2);
               resultP(i)=-log10(P(1,2))*sign(R(1,2));
               if isinf(resultP(i))
                   resultP(i)=NaN;
               end
                %mapR(i,j)=1-pdist2(temp1(good)',temp2(good)','cosine');
        end      
        