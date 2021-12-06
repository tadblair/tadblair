%% STRIATUM CELL TYPES

H1.hippo.spkwid(find(H1.hippo.spkwid==500 & log10(H1.hippo.spkwid)>.9))=H1.hippo.spkwid(find(H1.hippo.spkwid==500 & log10(H1.hippo.spkwid)>.9))-31.25*5;
H5.hippo.spkwid(find(H5.hippo.spkwid==500 & log10(H5.hippo.spkwid)>.9))=H5.hippo.spkwid(find(H5.hippo.spkwid==500 & log10(H5.hippo.spkwid)>.9))-31.25*5;
H7.hippo.spkwid(find(H7.hippo.spkwid==500 & log10(H7.hippo.spkwid)>.9))=H7.hippo.spkwid(find(H7.hippo.spkwid==500 & log10(H7.hippo.spkwid)>.9))-31.25*5;
%load rat data
%load 'H1'; load 'H5'; load 'H7';

%load 'gaussmixturemodel' allobj; %gaussian mixtre model for cell type classification
adbv=0.000000004577636718750000; %voltage conversion factor

figure(3); clf;
%histogram([SW_strNR; SW_strexc; SW_strinh],31.25:2*31.25:31*31.25);

%allobj = fitgmdist([[SW_strNR; SW_strexc; SW_strinh] [bucketPOP2_strNR; bucketPOP2_strexc; bucketPOP2_strinh]],4)

%% NON RESPONSIVE TO SWRs 
subplot(2,6,1); hold off;

SW_NR = [H7.striatum.spkwid(H7.striatum.NR_index)'; H5.striatum.spkwid(H5.striatum.NR_index)'; H1.striatum.spkwid(H1.striatum.NR_index)'];
SW_exc = [H7.striatum.spkwid(H7.striatum.exc_index)'; H5.striatum.spkwid(H5.striatum.exc_index)'; H1.striatum.spkwid(H1.striatum.exc_index)'];
SW_inh = [H7.striatum.spkwid(H7.striatum.inh_index)'; H5.striatum.spkwid(H5.striatum.inh_index)'; H1.striatum.spkwid(H1.striatum.inh_index)'];
SW_all = [SW_NR; SW_exc; SW_inh];
SW_swr = [SW_exc; SW_inh];


FR_NR = [H7.striatum.FR(H7.striatum.NR_index)'; H5.striatum.FR(H5.striatum.NR_index)'; H1.striatum.FR(H1.striatum.NR_index)'];
FR_exc = [H7.striatum.FR(H7.striatum.exc_index)'; H5.striatum.FR(H5.striatum.exc_index)'; H1.striatum.FR(H1.striatum.exc_index)'];
FR_inh = [H7.striatum.FR(H7.striatum.inh_index)'; H5.striatum.FR(H5.striatum.inh_index)'; H1.striatum.FR(H1.striatum.inh_index)'];
FR_all = [FR_NR; FR_exc; FR_inh];
FR_swr = [FR_exc; FR_inh];

idstr_NR = [H7.striatum.idstring(H7.striatum.NR_index); H5.striatum.idstring(H5.striatum.NR_index); H1.striatum.idstring(H1.striatum.NR_index)];
lowrate=find(FR_NR<1);
lowid=idstr_NR(lowrate,:);

POP2_NR = [H7.striatum.bucketPOP2(H7.striatum.NR_index)'; H5.striatum.bucketPOP2(H5.striatum.NR_index)'; H1.striatum.bucketPOP2(H1.striatum.NR_index)'];
POP2_exc = [H7.striatum.bucketPOP2(H7.striatum.exc_index)'; H5.striatum.bucketPOP2(H5.striatum.exc_index)'; H1.striatum.bucketPOP2(H1.striatum.exc_index)'];
POP2_inh = [H7.striatum.bucketPOP2(H7.striatum.inh_index)'; H5.striatum.bucketPOP2(H5.striatum.inh_index)'; H1.striatum.bucketPOP2(H1.striatum.inh_index)'];
POP2_all = [POP2_NR; POP2_exc; POP2_inh];

p_NR = [H7.striatum.SWR_pvalue(H7.striatum.NR_index)'; H5.striatum.SWR_pvalue(H5.striatum.NR_index)'; H1.striatum.SWR_pvalue(H1.striatum.NR_index)'];
p_exc = [H7.striatum.SWR_pvalue(H7.striatum.exc_index)'; H5.striatum.SWR_pvalue(H5.striatum.exc_index)'; H1.striatum.SWR_pvalue(H1.striatum.exc_index)'];
p_inh = [H7.striatum.SWR_pvalue(H7.striatum.inh_index)'; H5.striatum.SWR_pvalue(H5.striatum.inh_index)'; H1.striatum.SWR_pvalue(H1.striatum.inh_index)'];
p_all = [p_NR; p_exc; p_inh];

%SWR_pvalue_strNR = [H7.striatum.SWR_pvalue(H7.striatum.NR_index); H5.striatum.SWR_pvalue(H5.striatum.NR_index); H1.striatum.SWR_pvalue(H1.striatum.NR_index)];
strNRwave = [H7.striatum.spkwave(H7.striatum.NR_index,:); H5.striatum.spkwave(H5.striatum.NR_index,:); H1.striatum.spkwave(H1.striatum.NR_index,:)];
strexcwave = [H7.striatum.spkwave(H7.striatum.exc_index,:); H5.striatum.spkwave(H5.striatum.exc_index,:); H1.striatum.spkwave(H1.striatum.exc_index,:)];
strinhwave = [H7.striatum.spkwave(H7.striatum.inh_index,:); H5.striatum.spkwave(H5.striatum.inh_index,:); H1.striatum.spkwave(H1.striatum.inh_index,:)];
strallwave = [strNRwave; strexcwave; strinhwave];



% idx=cluster(allobj,[SW_strNR SWR_pvalue_strNR]);
%  d=1; scatter(SWR_pvalue_strNR(find(idx==d)), SW_strNR(find(idx==d)),'g.'); hold on;
%  d=2; scatter(SWR_pvalue_strNR(find(idx==d)), SW_strNR(find(idx==d)),'r.');
%  d=3; scatter(SWR_pvalue_strNR(find(idx==d)), SW_strNR(find(idx==d)),'k.');
%  d=4; scatter(SWR_pvalue_strNR(find(idx==d)), SW_strNR(find(idx==d)),'b.');

 %ezcontour(@(y,x)pdf(allobj,[x y]),[0 1],[200 800]); title(''); xlabel(''); ylabel('');
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

clear percSWR*;

for i=1:8
    YY(:,i)=Y(:,(i-1)*2+1)+Y(:,(i-1)*2+2);
end
Y=YY; clear YY;
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
Y=YY; clear YY;
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

% subplot(3,2,2); hold off;
%  m=nanmean(strallwave(find(SW_all>530),:)); s=nanstd(strallwave(find(SW_all>530),:))/sqrt(length(find(SW_all>530)));
%  m=m*adbv*2000000;  s=s*adbv*2000000; N1=length(find(SW_all>530));
%  shadedErrorBar(-(1000/32)+[1:32]*(1000/32),m,s,'lineProps',{'k'}); hold on;
%  set(gca,'YLim',[-100 150]);

 subplot(3,2,6); hold off;
%   m=nanmean(strallwave(find(SW_all<=340 & log10(FR_all)<=0),:)); s=nanstd(strallwave(find(SW_all<350 & log10(FR_all)<=0),:))/sqrt(length(find(SW_all<350 & log10(FR_all)<=0)));
%  m=m*adbv*2000000;  s=s*adbv*2000000; N1=length(find(SW_all<350 & log10(FR_all)<=0));
%  shadedErrorBar(-(1000/32)+[1:32]*(1000/32),m,s,'lineProps',{'g'}); hold on;
% 
%    m=nanmean(strallwave(find(SW_all<=340 & log10(FR_all)>0),:)); s=nanstd(strallwave(find(SW_all<350 & log10(FR_all)>0),:))/sqrt(length(find(SW_all<350 & log10(FR_all)>0)));
%  m=m*adbv*2000000;  s=s*adbv*2000000; N1=length(find(SW_all<350 & log10(FR_all)>0));
%  shadedErrorBar(-(1000/32)+[1:32]*(1000/32),m,s,'lineProps',{'b'}); hold on;

   m=nanmean(strallwave(find(SW_all<=340 & log10(FR_all)<0),:)); s=nanstd(strallwave(find(SW_all<350 & log10(FR_all)<0),:))/sqrt(length(find(SW_all<350 & log10(FR_all)<0)));
 m=m*adbv*2000000;  s=s*adbv*2000000; N1=length(find(SW_all<350 & log10(FR_all)<0));
 shadedErrorBar(-(1000/32)+[1:32]*(1000/32),m,s,'lineProps',{'b'}); hold on;
 set(gca,'YLim',[-100 150]);
 subplot(3,2,2); hold off;
   m=nanmean(strallwave(find(SW_all<=340 & log10(FR_all)>0),:)); s=nanstd(strallwave(find(SW_all<350 & log10(FR_all)>0),:))/sqrt(length(find(SW_all<350 & log10(FR_all)>0)));
 m=m*adbv*2000000;  s=s*adbv*2000000; N1=length(find(SW_all<350 & log10(FR_all)>0));
 shadedErrorBar(-(1000/32)+[1:32]*(1000/32),m,s,'lineProps',{'b'}); hold on;
 set(gca,'YLim',[-100 150]);
 
 
 [length(find(SW_all>350)) length(find(SW_all<=350))]
[length(find(SW_all>350)) length(find(SW_all<=350))]/366
 [length(find(SW_exc<=350)) length(find(SW_inh<=350)) length(find(SW_NR<=350))]
[length(find(SW_exc>350)) length(find(SW_inh>350)) length(find(SW_NR>350))]


p=ranksum(log10(FR_all(find(SW_all<=350))),log10(FR_all(find(SW_all>350)))); %pyr vs int
p=ranksum(log10(FR_swr(find(SW_swr>350))),log10(FR_NR(find(SW_NR>350))))
p=ranksum(log10(FR_swr(find(SW_swr<=350))),log10(FR_NR(find(SW_NR<=350))))
p=ranksum(log10(FR_exc(find(SW_exc>350))),log10(FR_inh(find(SW_inh>350))))
p=ranksum(log10(FR_exc(find(SW_exc<=350))),log10(FR_inh(find(SW_inh<=350))))

p=ranksum(log10(FR_all(find(SW_swr>350))),log10(FR_all(find(SW_NR>350))))

[r, p]=corr(log10(p_all(find(SW_swr>350))),log10(FR_all(find(SW_swr>350))))
[r, p]=corr((p_all(find(SW_swr>350))),(FR_all(find(SW_swr>350))))

[r, p]=corr(log10(p_all(find(SW_swr<=350))),log10(FR_all(find(SW_swr<=350))))
[r, p]=corr((p_all(find(SW_swr<=350))),(FR_all(find(SW_swr<=350))))

 figure(200); subplot(1,3,3); hold off; bar([percSWR_MSN(1,1) percSWR_INT(1,1); percSWR_MSN(1,2) percSWR_INT(1,2); percSWR_MSN(1,3) percSWR_INT(1,3); percSWR_MSN(1,4) percSWR_INT(1,4)],1);% hold on; plot(percSWR_INT(1,:)); set(gca,'YLim',[0 1]);
set(gca,'YLim',[0 1]);

 %  d=2; m=mean(strNRwave(find(idx==d),:)); s=std(strNRwave(find(idx==d),:))/sqrt(length(find(idx==d)));
%  m=m*adbv*2000000-200; s=s*adbv*2000000; N2=length(find(idx==d));
%  shadedErrorBar(1:32,m,s,'lineProps',{'r'}); hold on;
%  d=3; m=mean(strNRwave(find(idx==d),:)); s=std(strNRwave(find(idx==d),:))/sqrt(length(find(idx==d)));
%  m=m*adbv*2000000; s=s*adbv*2000000; N3=length(find(idx==d));
%  shadedErrorBar(33:64,m,s,'lineProps',{'k'}); hold on;
%  d=4; m=mean(strNRwave(find(idx==d & SW_strNR<350 & SWR_pvalue_strNR<.6),:)); s=std(strNRwave(find(idx==d & SW_strNR<350 & SWR_pvalue_strNR<.6),:))/sqrt(length(find(idx==d)));
%  m=m*adbv*2000000-200; s=s*adbv*2000000; N4=length(find(idx==d));
%  shadedErrorBar(33:64,m,s,'lineProps',{'b'}); hold on;
%  axis tight;
% subplot(2,6,2+6); hold off;
% labels={'G','K','R','B'};
% pie([N1 N3 N2 N4],labels);
% [N1 N3 N2 N4]
% %% EXCITED BY SWRs 
% subplot(2,6,3); hold off;
% 
% SW_strexc = [H7.striatum.spkwid(H7.striatum.exc_index); H5.striatum.spkwid(H5.striatum.exc_index); H1.striatum.spkwid(H1.striatum.exc_index)];
% SWR_pvalue_strexc = [H7.striatum.SWR_pvalue(H7.striatum.exc_index); H5.striatum.SWR_pvalue(H5.striatum.exc_index); H1.striatum.SWR_pvalue(H1.striatum.exc_index)];
% strexcwave = [H7.striatum.spkwave(H7.striatum.exc_index,:); H5.striatum.spkwave(H5.striatum.exc_index,:); H1.striatum.spkwave(H1.striatum.exc_index,:)];
% 
% idx=cluster(allobj,[SW_strexc SWR_pvalue_strexc]);
%  d=1; scatter(SWR_pvalue_strexc(find(idx==d)), SW_strexc(find(idx==d)),'g.'); hold on;
%  d=2; scatter(SWR_pvalue_strexc(find(idx==d)), SW_strexc(find(idx==d)),'r.');
%  d=3; scatter(SWR_pvalue_strexc(find(idx==d)), SW_strexc(find(idx==d)),'k.');
%  d=4; scatter(SWR_pvalue_strexc(find(idx==d)), SW_strexc(find(idx==d)),'b.');
% %ezcontour(@(y,x)pdf(allobj,[x y]),[0 1],[200 800]); title(''); xlabel(''); ylabel('');
% subplot(2,6,4); hold off;
%  d=1; m=mean(strexcwave(find(idx==d),:)); s=std(strexcwave(find(idx==d),:))/sqrt(length(find(idx==d)));
%  m=m*adbv*2000000;  s=s*adbv*2000000; N1=length(find(idx==d));
%  shadedErrorBar(1:32,m,s,'lineProps',{'g'}); hold on;
%  d=2; m=mean(strexcwave(find(idx==d),:)); s=std(strexcwave(find(idx==d),:))/sqrt(length(find(idx==d)));
%  m=m*adbv*2000000-200; s=s*adbv*2000000; N2=length(find(idx==d));
%  shadedErrorBar(1:32,m,s,'lineProps',{'r'}); hold on;
%  d=3; m=mean(strexcwave(find(idx==d),:)); s=std(strexcwave(find(idx==d),:))/sqrt(length(find(idx==d)));
%  m=m*adbv*2000000; s=s*adbv*2000000; N3=length(find(idx==d));
%  shadedErrorBar(33:64,m,s,'lineProps',{'k'}); hold on;
%  d=4; m=mean(strexcwave(find(idx==d),:)); s=std(strexcwave(find(idx==d),:))/sqrt(length(find(idx==d)));
%  m=m*adbv*2000000-200; s=s*adbv*2000000; N4=length(find(idx==d));
%  shadedErrorBar(33:64,m,s,'lineProps',{'b'}); hold on;
%  axis tight;
%  subplot(2,6,4+6); hold off;
% labels={'G','K','R','B'};
% pie([N1 N3 N2 N4],labels);
% [N1 N3 N2 N4]
% %% INHIBITED BY SWRs 
% subplot(2,6,5); hold off;
% 
% SW_strinh = [H7.striatum.spkwid(H7.striatum.inh_index); H5.striatum.spkwid(H5.striatum.inh_index); H1.striatum.spkwid(H1.striatum.inh_index)];
% SWR_pvalue_strinh = [H7.striatum.SWR_pvalue(H7.striatum.inh_index); H5.striatum.SWR_pvalue(H5.striatum.inh_index); H1.striatum.SWR_pvalue(H1.striatum.inh_index)];
% strinhwave = [H7.striatum.spkwave(H7.striatum.inh_index,:); H5.striatum.spkwave(H5.striatum.inh_index,:); H1.striatum.spkwave(H1.striatum.inh_index,:)];
% 
% idx=cluster(allobj,[SW_strinh SWR_pvalue_strinh]);
%  d=1; scatter(SWR_pvalue_strinh(find(idx==d)), SW_strinh(find(idx==d)),'g.'); hold on;
%  d=2; scatter(SWR_pvalue_strinh(find(idx==d)), SW_strinh(find(idx==d)),'r.');
%  d=3; scatter(SWR_pvalue_strinh(find(idx==d)), SW_strinh(find(idx==d)),'k.');
%  d=4; scatter(SWR_pvalue_strinh(find(idx==d)), SW_strinh(find(idx==d)),'b.');
% %ezcontour(@(y,x)pdf(allobj,[x y]),[0 1],[200 800]); title(''); xlabel(''); ylabel('');
% subplot(2,6,6); hold off;
%  d=1; m=mean(strinhwave(find(idx==d),:)); s=std(strinhwave(find(idx==d),:))/sqrt(length(find(idx==d)));
%  m=m*adbv*2000000;  s=s*adbv*2000000; N1=length(find(idx==d));
%  shadedErrorBar(1:32,m,s,'lineProps',{'g'}); hold on;
%  d=2; m=mean(strinhwave(find(idx==d),:)); s=std(strinhwave(find(idx==d),:))/sqrt(length(find(idx==d)));
%  m=m*adbv*2000000-200; s=s*adbv*2000000; N2=length(find(idx==d));
%  shadedErrorBar(1:32,m,s,'lineProps',{'r'}); hold on;
%  d=3; m=mean(strinhwave(find(idx==d),:)); s=std(strinhwave(find(idx==d),:))/sqrt(length(find(idx==d)));
%  m=m*adbv*2000000; s=s*adbv*2000000; N3=length(find(idx==d));
%  shadedErrorBar(33:64,m,s,'lineProps',{'k'}); hold on;
%  d=4; m=strinhwave(find(idx==d),:);% s=std(strinhwave(find(idx==d),:))/sqrt(length(find(idx==d)));
%  m=m*adbv*2000000-200; N4=length(find(idx==d));
%  plot(33:64,m,'b'); hold on;
%  axis tight;
% subplot(2,6,6+6); hold off;
% labels={'G','K','R','B'};
% pie([N1 N3 N2 N4],labels);
% [N1 N3 N2 N4]
% % dexes(tt1)
% 
% figure(100); clf; 
% 
% strwids=[SW_strexc; SW_strinh; SW_strNR];
% strwids(find(strwids<200))=250;
% strwids(find(strwids>900))=750;
% subplot(3,2,3); histogram(strwids,31.25:2*31.25:31*31.25);
% set(gca,'XLim',[200 900]);
% 
% strwids=[SW_strexc; SW_strinh; SW_strNR];
% strwids(find(strwids<200))=250;
% strwids(find(strwids>900))=750;
% subplot(3,2,5); histogram(strwids,31.25:2*31.25:31*31.25);
% set(gca,'XLim',[200 900]);
% 
% strpops=[SWR_pvalue_strexc; SWR_pvalue_strinh; SWR_pvalue_strNR];
% subplot(3,2,4); histogram(strpops,0:.1:1);
% %set(gca,'XLim',[200 900]);
% 
% strpops=[SWR_pvalue_strexc; SWR_pvalue_strinh; SWR_pvalue_strNR];
% subplot(3,2,6); histogram(strpops,0:.1:1);
% %set(gca,'XLim',[200 900]);