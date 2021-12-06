%% STRIATUM CELL TYPES

%load rat data
%load 'H1'; load 'H5'; load 'H7';

%load 'gaussmixturemodel' allobj; %gaussian mixtre model for cell type classification
adbv=0.000000004577636718750000; %voltage conversion factor

figure(1); clf;

%% NON RESPONSIVE TO SWRs 
subplot(2,6,1); hold off;

SW_NR = [H7.hippo.spkwid(H7.hippo.NR_index)'; H5.hippo.spkwid(H5.hippo.NR_index)'; H1.hippo.spkwid(H1.hippo.NR_index)'];
SW_exc = [H7.hippo.spkwid(H7.hippo.exc_index)'; H5.hippo.spkwid(H5.hippo.exc_index)'; H1.hippo.spkwid(H1.hippo.exc_index)'];
SW_inh = [H7.hippo.spkwid(H7.hippo.inh_index)'; H5.hippo.spkwid(H5.hippo.inh_index)'; H1.hippo.spkwid(H1.hippo.inh_index)'];
SW_all = [SW_NR; SW_exc; SW_inh];
SW_swr = [SW_exc; SW_inh];

FR_NR = [H7.hippo.FR(H7.hippo.NR_index)'; H5.hippo.FR(H5.hippo.NR_index)'; H1.hippo.FR(H1.hippo.NR_index)'];
FR_exc = [H7.hippo.FR(H7.hippo.exc_index)'; H5.hippo.FR(H5.hippo.exc_index)'; H1.hippo.FR(H1.hippo.exc_index)'];
FR_inh = [H7.hippo.FR(H7.hippo.inh_index)'; H5.hippo.FR(H5.hippo.inh_index)'; H1.hippo.FR(H1.hippo.inh_index)'];
FR_all = [FR_NR; FR_exc; FR_inh];
FR_swr = [FR_exc; FR_inh];

POP2_NR = [H7.hippo.bucketPOP2(H7.hippo.NR_index)'; H5.hippo.bucketPOP2(H5.hippo.NR_index)'; H1.hippo.bucketPOP2(H1.hippo.NR_index)'];
POP2_exc = [H7.hippo.bucketPOP2(H7.hippo.exc_index)'; H5.hippo.bucketPOP2(H5.hippo.exc_index)'; H1.hippo.bucketPOP2(H1.hippo.exc_index)'];
POP2_inh = [H7.hippo.bucketPOP2(H7.hippo.inh_index)'; H5.hippo.bucketPOP2(H5.hippo.inh_index)'; H1.hippo.bucketPOP2(H1.hippo.inh_index)'];
POP2_all = [POP2_NR; POP2_exc; POP2_inh];
POP2_swr = [POP2_exc; POP2_inh];

p_NR = [H7.hippo.SWR_pvalue(H7.hippo.NR_index); H5.hippo.SWR_pvalue(H5.hippo.NR_index); H1.hippo.SWR_pvalue(H1.hippo.NR_index)];
p_exc = [H7.hippo.SWR_pvalue(H7.hippo.exc_index); H5.hippo.SWR_pvalue(H5.hippo.exc_index); H1.hippo.SWR_pvalue(H1.hippo.exc_index)];
p_inh = [H7.hippo.SWR_pvalue(H7.hippo.inh_index); H5.hippo.SWR_pvalue(H5.hippo.inh_index); H1.hippo.SWR_pvalue(H1.hippo.inh_index)];
p_all = [p_NR; p_exc; p_inh];
p_swr = [p_exc; p_inh];
%SWR_pvalue_hipNR = [H7.hippo.SWR_pvalue(H7.hippo.NR_index); H5.hippo.SWR_pvalue(H5.hippo.NR_index); H1.hippo.SWR_pvalue(H1.hippo.NR_index)];
hipNRwave = [H7.hippo.spkwave(H7.hippo.NR_index,:); H5.hippo.spkwave(H5.hippo.NR_index,:); H1.hippo.spkwave(H1.hippo.NR_index,:)];
hipexcwave = [H7.hippo.spkwave(H7.hippo.exc_index,:); H5.hippo.spkwave(H5.hippo.exc_index,:); H1.hippo.spkwave(H1.hippo.exc_index,:)];
hipinhwave = [H7.hippo.spkwave(H7.hippo.inh_index,:); H5.hippo.spkwave(H5.hippo.inh_index,:); H1.hippo.spkwave(H1.hippo.inh_index,:)];
hipallwave = [hipNRwave; hipexcwave; hipinhwave];

 %ezcontour(@(y,x)pdf(allobj,[x y]),[0 1],[200 800]); title(''); xlabel(''); ylabel('');
subplot(3,2,3); hold off;
scatter(log10(FR_NR),SW_NR,'k.'); hold on;
scatter(log10(FR_exc),SW_exc,'b.');
scatter(log10(FR_inh),SW_inh,'r.');
set(gca,'YLim',[100 900]);

subplot(3,2,1); hold off;
temp=histogram(log10(FR_inh(find(SW_inh>470))),-2:.25:2); hold on;
Y=temp.Values;
temp=histogram(log10(FR_exc(find(SW_exc>470))),-2:.25:2); hold on;
Y=[Y; temp.Values];
temp=histogram(log10(FR_NR(find(SW_NR>470))),-2:.25:2); hold on;
Y=[Y; temp.Values];
bar(temp.BinEdges(1:end-1)+.25/2,Y',1,'stacked');
title([num2str(median(FR_all(find(SW_all>500))))]);
axis tight; set(gca,'XLim',[-2 2]);

for i=1:8
    YY(:,i)=Y(:,(i-1)*2+1)+Y(:,(i-1)*2+2);
end
Y=YY;
clear YY;
for i=1:4
    YY(:,i)=Y(:,(i-1)*2+1)+Y(:,(i-1)*2+2);
end
Y=YY;
percSWR_PYR=[sum(Y(1:2,:))./sum(Y); Y(3,:)./sum(Y)];

subplot(3,2,5); hold off;
[length(find(SW_exc<=470)) length(find(SW_inh<=470)) length(find(SW_NR<=470))]
[length(find(SW_exc>470)) length(find(SW_inh>470)) length(find(SW_NR>470))]
p=ranksum(log10(FR_all(find(SW_all<=470))),log10(FR_all(find(SW_all>470)))); %pyr vs int
p=ranksum(log10(FR_all(find(SW_swr<=470))),log10(FR_all(find(SW_NR<=470))))
p=ranksum(log10(FR_all(find(SW_swr>470))),log10(FR_all(find(SW_NR>470))))
p=ranksum(log10(FR_all(find(SW_exc<=470))),log10(FR_all(find(SW_inh<=470))))
p=ranksum(log10(FR_all(find(SW_exc>470))),log10(FR_all(find(SW_inh>470))))
temp=histogram(log10(FR_inh(find(SW_inh<=470))),-2:.25:2); 
Y=temp.Values;
temp=histogram(log10(FR_exc(find(SW_exc<=470))),-2:.25:2); 
Y=[Y; temp.Values];
temp=histogram(log10(FR_NR(find(SW_NR<=470))),-2:.25:2); 
Y=[Y; temp.Values];
bar(temp.BinEdges(1:end-1)+.25/2,Y',1,'stacked');
axis tight; set(gca,'XLim',[-2 2]);

for i=1:8
    YY(:,i)=Y(:,(i-1)*2+1)+Y(:,(i-1)*2+2);
end
Y=YY;
clear YY;
for i=1:4
    YY(:,i)=Y(:,(i-1)*2+1)+Y(:,(i-1)*2+2);
end
Y=YY;
percSWR_INT=[sum(Y(1:2,:))./sum(Y); Y(3,:)./sum(Y)];

subplot(3,2,4); hold off;
temp=histogram(SW_inh,31.25*2:2*31.25:31*31.25);
Y=temp.Values;
temp=histogram(SW_exc,31.25*2:2*31.25:31*31.25);
Y=[Y; temp.Values];
temp=histogram(SW_NR,31.25*2:2*31.25:31*31.25);
Y=[Y; temp.Values];
barh(temp.BinEdges(1:end-1)+31.25,Y',1,'stacked');
%barh(temp.BinEdges(1:end-1)+31.25,temp.Values,1);
axis tight; set(gca,'YLim',[100 900]);

subplot(3,2,2); hold off;
 m=nanmean(hipallwave(find(SW_all>530),:)); s=nanstd(hipallwave(find(SW_all>530),:))/sqrt(length(find(SW_all>530)));
 m=m*adbv*2000000;  s=s*adbv*2000000; N1=length(find(SW_all>530));
 shadedErrorBar(-(1000/32)+[1:32]*(1000/32),m,s,'lineProps',{'k'}); hold on;
 set(gca,'YLim',[-100 150]);

 subplot(3,2,6); hold off;
  m=nanmean(hipallwave(find(SW_all<=340),:)); s=nanstd(hipallwave(find(SW_all<=340),:))/sqrt(length(find(SW_all<=340)));
 m=m*adbv*2000000;  s=s*adbv*2000000; N1=length(find(SW_all<=340));
 shadedErrorBar(-(1000/32)+[1:32]*(1000/32),m,s,'lineProps',{'b'}); hold on;
 set(gca,'YLim',[-100 150]);

 figure(200); subplot(1,3,1); hold off; bar([percSWR_PYR(1,1) percSWR_INT(1,1); percSWR_PYR(1,2) percSWR_INT(1,2); percSWR_PYR(1,3) percSWR_INT(1,3); percSWR_PYR(1,4) percSWR_INT(1,4)],1);% hold on; plot(percSWR_INT(1,:)); set(gca,'YLim',[0 1]);
set(gca,'YLim',[0 1]);

