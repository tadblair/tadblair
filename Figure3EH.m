%% STRIATUM CELL TYPES

%load rat data
%load 'H1'; load 'H5'; load 'H7';

%load 'gaussmixturemodel' allobj; %gaussian mixtre model for cell type classification
adbv=0.000000004577636718750000; %voltage conversion factor

%% NON RESPONSIVE TO SWRs 

SW_NR = [H7.septum.spkwid(H7.septum.NR_index)'; H5.septum.spkwid(H5.septum.NR_index)'];
SW_exc = [H7.septum.spkwid(H7.septum.exc_index)'; H5.septum.spkwid(H5.septum.exc_index)'];
SW_inh = [H7.septum.spkwid(H7.septum.inh_index)'; H5.septum.spkwid(H5.septum.inh_index)'];
SW_all = [SW_NR; SW_exc; SW_inh];
SW_swr = [SW_exc; SW_inh];

FR_NR = [H7.septum.FR(H7.septum.NR_index)'; H5.septum.FR(H5.septum.NR_index)'];
FR_exc = [H7.septum.FR(H7.septum.exc_index)'; H5.septum.FR(H5.septum.exc_index)'];
FR_inh = [H7.septum.FR(H7.septum.inh_index)'; H5.septum.FR(H5.septum.inh_index)'];
FR_all = [FR_NR; FR_exc; FR_inh];
FR_swr = [FR_exc; FR_inh];

POP2_NR = [H7.septum.bucketPOP2(H7.septum.NR_index)'; H5.septum.bucketPOP2(H5.septum.NR_index)'];
POP2_exc = [H7.septum.bucketPOP2(H7.septum.exc_index)'; H5.septum.bucketPOP2(H5.septum.exc_index)'];
POP2_inh = [H7.septum.bucketPOP2(H7.septum.inh_index)'; H5.septum.bucketPOP2(H5.septum.inh_index)'];
POP2_all = [POP2_NR; POP2_exc; POP2_inh];
POP2_swr = [POP2_exc; POP2_inh];

p_NR = [H7.septum.SWR_pvalue(H7.septum.NR_index)'; H5.septum.SWR_pvalue(H5.septum.NR_index)'];
p_exc = [H7.septum.SWR_pvalue(H7.septum.exc_index)'; H5.septum.SWR_pvalue(H5.septum.exc_index)'];
p_inh = [H7.septum.SWR_pvalue(H7.septum.inh_index)'; H5.septum.SWR_pvalue(H5.septum.inh_index)'];
p_all = [p_NR; p_exc; p_inh];
p_swr = [p_exc; p_inh];

excpeth = [H7.septum.bucketpeth(:,H7.septum.exc_index)'; H5.septum.bucketpeth(:,H5.septum.exc_index)'];
inhpeth = [H7.septum.bucketpeth(:,H7.septum.inh_index)'; H5.septum.bucketpeth(:,H5.septum.inh_index)'];

figure(1); clf;
for i=1:length(find(SW_exc>350))
  subplot(7,7,i);
  bar(-495:10:495,excpeth(i,:),1); axis tight;
end
figure(3); clf;
for i=1:length(find(SW_exc<=350))
  subplot(4,4,i);
  bar(-495:10:495,excpeth(i,:),1); axis tight;
end
%SWR_pvalue_sepNR = [H7.septum.SWR_pvalue(H7.septum.NR_index); H5.septum.SWR_pvalue(H5.septum.NR_index); H1.septum.SWR_pvalue(H1.septum.NR_index)];
sepNRwave = [H7.septum.spkwave(H7.septum.NR_index,:); H5.septum.spkwave(H5.septum.NR_index,:)];
sepexcwave = [H7.septum.spkwave(H7.septum.exc_index,:); H5.septum.spkwave(H5.septum.exc_index,:)];
sepinhwave = [H7.septum.spkwave(H7.septum.inh_index,:); H5.septum.spkwave(H5.septum.inh_index,:)];
sepallwave = [sepNRwave; sepexcwave; sepinhwave];


figure(2); clf;
subplot(2,6,1); hold off;

subplot(3,2,3); hold off;
scatter(log10(FR_NR),SW_NR,'k.'); hold on;
scatter(log10(FR_exc),SW_exc,'b.');
scatter(log10(FR_inh),SW_inh,'r.');
set(gca,'YLim',[100 900],'XLim',[-2 2]);

subplot(3,2,1); hold off;
temp=histogram(log10(FR_inh(find(SW_inh>350))),-2:.25:2); hold on;
Y=temp.Values;
temp=histogram(log10(FR_exc(find(SW_exc>350))),-2:.25:2); hold on;
Y=[Y; temp.Values];
temp=histogram(log10(FR_NR(find(SW_NR>350))),-2:.25:2); hold on;
Y=[Y; temp.Values];
bar(temp.BinEdges(1:end-1)+.25/2,Y',1,'stacked');
title([num2str(median(FR_all(find(SW_all>350))))]);
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
percSWR_MSN=[sum(Y(1:2,:))./sum(Y); Y(3,:)./sum(Y)];

subplot(3,2,5); hold off;
temp=histogram(log10(FR_inh(find(SW_inh<=350))),-2:.25:2); 
Y=temp.Values;
temp=histogram(log10(FR_exc(find(SW_exc<=350))),-2:.25:2); 
Y=[Y; temp.Values];
temp=histogram(log10(FR_NR(find(SW_NR<=350))),-2:.25:2); 
Y=[Y; temp.Values];
bar(temp.BinEdges(1:end-1)+.25/2,Y',1,'stacked');
title([num2str(median(FR_all(find(SW_all<350))))]);
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
 m=nanmean(sepallwave(find(SW_all>530),:)); s=nanstd(sepallwave(find(SW_all>530),:))/sqrt(length(find(SW_all>530)));
 m=m*adbv*2000000;  s=s*adbv*2000000; N1=length(find(SW_all>530));
 shadedErrorBar(-(1000/32)+[1:32]*(1000/32),m,s,'lineProps',{'k'}); hold on;

 subplot(3,2,6); hold off;

  m=nanmean(sepallwave(find(SW_all<=340 & log10(FR_all)<Inf),:)); s=nanstd(sepallwave(find(SW_all<350 & log10(FR_all)<Inf),:))/sqrt(length(find(SW_all<350 & log10(FR_all)<Inf)));
 m=m*adbv*2000000;  s=s*adbv*2000000; N1=length(find(SW_all<350 & log10(FR_all)<Inf));
 shadedErrorBar(-(1000/32)+[1:32]*(1000/32),m,s,'lineProps',{'b'}); hold on;
 set(gca,'YLim',[-100 150]);
 
[length(find(SW_all>350)) length(find(SW_all<=350))]
[length(find(SW_all>350)) length(find(SW_all<=350))]/226
 [length(find(SW_exc<=350)) length(find(SW_inh<=350)) length(find(SW_NR<=350))]
[length(find(SW_exc>350)) length(find(SW_inh>350)) length(find(SW_NR>350))]

p=ranksum(log10(FR_all(find(SW_all<=350))),log10(FR_all(find(SW_all>350)))); %pyr vs int
p=ranksum(log10(FR_all(find(SW_swr>350))),log10(FR_all(find(SW_NR>350))))
p=ranksum(log10(FR_all(find(SW_swr<=350))),log10(FR_all(find(SW_NR<=350))))
p=ranksum(log10(FR_all(find(SW_exc>350))),log10(FR_all(find(SW_inh>350))))
p=ranksum(log10(FR_all(find(SW_exc<=350))),log10(FR_all(find(SW_inh<=350))))


 figure(200); subplot(1,3,2); hold off; bar([percSWR_MSN(1,1) percSWR_INT(1,1); percSWR_MSN(1,2) percSWR_INT(1,2); percSWR_MSN(1,3) percSWR_INT(1,3); percSWR_MSN(1,4) percSWR_INT(1,4)],1);% hold on; plot(percSWR_INT(1,:)); set(gca,'YLim',[0 1]);
set(gca,'YLim',[0 1]);
