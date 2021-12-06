clear;

load 'C:\Users\Blair Lab\Documents\MATLAB\AndrewData\revisions\H7r';
load 'C:\Users\Blair Lab\Documents\MATLAB\AndrewData\revisions\H5r';
load 'C:\Users\Blair Lab\Documents\MATLAB\AndrewData\revisions\H1r';



%hipcells;

hip_exc_r=[];
hip_exc_m=[];
hip_inh_r=[];
hip_inh_m=[];
hip_NR_r=[];
hip_NR_m=[];

clear  H1_hip_* H5_hip_* H7_hip_*; 
%%%% ----------------------------------------------------------------- H1

% load h1results;
% for i=1:length(results_idStr)
%     temp=results_idStr{i};
%     kk=strfind(temp,'tt');
%     ll=strfind(temp(kk:end),'_');
%     ll=ll(1)+kk-2;
%     tt=str2num(temp((kk+2):ll));
%     results(8,i)=tt;
% end
% r=1;
% hipdex=[];
% for i=1:length(rat(r).tt)
%     if ~isempty(rat(r).tt(i))
%         for j=1:length(rat(r).tt(i).cell)
%             temp=find(results(10,:)==rat(r).tt(i).cell(1,j).sess_clust(1) & results(9,:)==rat(r).tt(i).cell(1,j).sess_clust(2)  & results(8,:)==i);
%             hipdex=[hipdex; temp];
%         end
%     end
% end

H7_hip_exc = H7.hippo.exc_index;%[10 11 64:66 104:107 139 141 174 313 234 281 311 312 67 108 109 146 111 144 316 179 181 317];
H7_hip_inh = H7.hippo.inh_index;%[314];
H7_hip_NR = H7.hippo.NR_index;%setdiff(hipdex,[H1_hip_exc H1_hip_inh]);


phdists=H7.hippo.phasedist(:,H7_hip_exc);
for i=1:length(H7_hip_exc)
   H7_hip_exc_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H7_hip_exc_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H7_hip_exc_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
figure(100);
goodp=find(H7_hip_exc_p<.01);
subplot(3,4,2); %histogram(H7_hip_exc_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H7_hip_exc_m(goodp),[-pi:(2*pi/12):pi]+pi/12);
%subplot(3,2,2); histogram(H7_hip_exc_r,0:.1:1); title('H7 hip exc');
[H7_stats_hip_exc.Zp H7_stats_hip_exc.Z] = circ_rtest(H7_hip_exc_m(goodp)');
[H7_stats_hip_exc.Vp H7_stats_hip_exc.V] = circ_vtest(H7_hip_exc_m(goodp)',0);
H7_stats_hip_exc.df=length(goodp);
hip_exc_r=[hip_exc_r H7_hip_exc_r];
H7_hip_exc_mean = circ_mean(H7_hip_exc_m(goodp)');
hip_exc_m=[hip_exc_m circ_dist(H7_hip_exc_m(goodp),H7_hip_exc_mean)];

phdists=H7.hippo.phasedist(:,H7_hip_inh);
for i=1:length(H7_hip_inh)
   H7_hip_inh_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H7_hip_inh_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H7_hip_inh_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H7_hip_inh_p<.01);
subplot(3,4,6); %histogram(H7_hip_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H7_hip_inh_m(goodp),[-pi:(2*pi/12):pi]+pi/12);
H7_hip_inh_mean = circ_mean(H7_hip_inh_m(goodp)');
%subplot(3,2,4); histogram(H7_hip_inh_r,0:.1:1); title('H7 hip inh');
[H7_stats_hip_inh.Zp H7_stats_hip_inh.Z] = circ_rtest(H7_hip_inh_m(goodp)');
[H7_stats_hip_inh.Vp H7_stats_hip_inh.V] = circ_vtest(H7_hip_inh_m(goodp)',0);
H7_stats_hip_inh.df=length(goodp);
hip_inh_r=[hip_inh_r H7_hip_inh_r];
hip_inh_m=[hip_inh_m circ_dist(H7_hip_inh_m(goodp),H7_hip_exc_mean)];

phdists=H7.hippo.phasedist(:,H7_hip_NR);
for i=1:length(H7_hip_NR)
   H7_hip_NR_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H7_hip_NR_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H7_hip_NR_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H7_hip_NR_p<.01);
subplot(3,4,10); %histogram(H7_hip_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H7_hip_NR_m(goodp),[-pi:(2*pi/12):pi]+pi/12);
%subplot(3,2,6); histogram(H7_hip_NR_r,0:.1:1); title('H7 hip NR');
hip_NR_r=[hip_NR_r H7_hip_NR_r];
[H7_stats_hip_NR.Zp H7_stats_hip_NR.Z] = circ_rtest(H7_hip_NR_m(goodp)');
[H7_stats_hip_NR.Vp H7_stats_hip_NR.V] = circ_vtest(H7_hip_NR_m(goodp)',0);
H7_hip_NR_mean = circ_mean(H7_hip_NR_m(goodp)');
H7_stats_hip_NR.df=length(goodp);
hip_NR_m=[hip_NR_m circ_dist(H7_hip_NR_m(goodp),H7_hip_exc_mean)];


% phdists=H1.hippo.phasedist(:,H1_hip_NR);
% for i=1:length(H1_hip_NR)
%    H1_hip_NR_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
%    H1_hip_NR_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
%    H1_hip_NR_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
% end
% goodp=find(H1_hip_NR_p<.01);
% H1_hip_NR_mean = circ_mean(H1_hip_NR_m(goodp)');
% subplot(3,4,11); %histogram(H1_hip_inh_m(goodp),-pi:(2*pi/12):pi);
% polarhistogram(H1_hip_NR_m(goodp),[-pi:(2*pi/12):pi-(2*pi/12)]+(2*pi/12));
% %subplot(3,2,6); histogram(H1_hip_NR_r,0:.1:1); title('H1 hip NR');
% [H1_stats_hip_NR.Zp H1_stats_hip_NR.Z] = circ_rtest(H1_hip_NR_m(goodp)');
% [H1_stats_hip_NR.Vp H1_stats_hip_NR.V] = circ_vtest(H1_hip_NR_m(goodp)',pi);
% H1_stats_hip_NR.df=length(goodp);
% hip_NR_r=[hip_NR_r H1_hip_NR_r];
% hip_NR_m=[hip_NR_m circ_dist(H1_hip_NR_m(goodp),H1_hip_exc_mean)];

%% --------------------------------------------------------------------------------------------------------



% load h5results;
% for i=1:length(results_idStr)
%     temp=results_idStr{i};
%     kk=strfind(temp,'tt');
%     ll=strfind(temp(kk:end),'_');
%     ll=ll(1)+kk-2;
%     tt=str2num(temp((kk+2):ll));
%     results(8,i)=tt;
% end
% r=5;
% hipdex=[];
% for i=1:length(rat(r).tt)
%     if ~isempty(rat(r).tt(i))
%         for j=1:length(rat(r).tt(i).cell)
%             temp=find(results(10,:)==rat(r).tt(i).cell(1,j).sess_clust(1) & results(9,:)==rat(r).tt(i).cell(1,j).sess_clust(2)  & results(8,:)==i);
%             hipdex=[hipdex; temp];
%         end
%     end
% end

% H5_hip_exc = [439 168 446 447 448 540 750 301 376 174 176 303 379 304 547 549 2 3 4 459 761 460 461 462 463 464 465 764 5 6 177 766];
% H5_hip_inh = [71];
% H5_hip_NR = setdiff(hipdex,[H5_hip_exc H5_hip_inh]);

H5_hip_exc = H5.hippo.exc_index;%[10 11 64:66 104:107 139 141 174 313 234 281 311 312 67 108 109 146 111 144 316 179 181 317];
H5_hip_inh = H5.hippo.inh_index;%[314];
H5_hip_NR = H5.hippo.NR_index;%setdiff(hipdex,[H1_hip_exc H1_hip_inh]);

phdists=H5.hippo.phasedist(:,H5_hip_exc);
for i=1:length(H5_hip_exc)
   H5_hip_exc_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H5_hip_exc_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H5_hip_exc_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
figure(100);  
goodp=find(H5_hip_exc_p<.01);
subplot(3,4,3); %histogram(H5_hip_exc_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H5_hip_exc_m(goodp),[-pi:(2*pi/12):pi]+pi/12);
%subplot(3,2,2); histogram(H5_hip_exc_r,0:.1:1); title('H5 hip exc');
[H5_stats_hip_exc.Zp H5_stats_hip_exc.Z] = circ_rtest(H5_hip_exc_m(goodp)');
[H5_stats_hip_exc.Vp H5_stats_hip_exc.V] = circ_vtest(H5_hip_exc_m(goodp)',pi);
H5_stats_hip_exc.df=length(goodp);
hip_exc_r=[hip_exc_r H5_hip_exc_r];
H5_hip_exc_mean = circ_mean(H5_hip_exc_m(goodp)');
hip_exc_m=[hip_exc_m circ_dist(H5_hip_exc_m(goodp),H5_hip_exc_mean)];

phdists=H5.hippo.phasedist(:,H5_hip_inh);
for i=1:length(H5_hip_inh)
   H5_hip_inh_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H5_hip_inh_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H5_hip_inh_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H5_hip_inh_p<.01);
subplot(3,4,7); %histogram(H5_hip_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H5_hip_inh_m(goodp),[-pi:(2*pi/12):pi]+pi/12);
%subplot(3,2,4); histogram(H5_hip_inh_r,0:.1:1); title('H5 hip inh');
H5_hip_inh_mean = circ_mean(H5_hip_inh_m(goodp)');
[H5_stats_hip_inh.Zp H5_stats_hip_inh.Z] = circ_rtest(H5_hip_inh_m(goodp)');
[H5_stats_hip_inh.Vp H5_stats_hip_inh.V] = circ_vtest(H5_hip_inh_m(goodp)',pi);
H5_stats_hip_inh.df=length(goodp);
hip_inh_r=[hip_inh_r H5_hip_inh_r];
hip_inh_m=[hip_inh_m circ_dist(H5_hip_inh_m(goodp),H5_hip_exc_mean)];

phdists=H5.hippo.phasedist(:,H5_hip_NR);
for i=1:length(H5_hip_NR)
   H5_hip_NR_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H5_hip_NR_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H5_hip_NR_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H5_hip_NR_p<.01);
subplot(3,4,11); %histogram(H5_hip_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H5_hip_NR_m(goodp),[-pi:(2*pi/12):pi]+pi/12);
H5_hip_NR_mean = circ_mean(H5_hip_NR_m(goodp)');
%subplot(3,2,6); histogram(H5_hip_NR_r,0:.1:1); title('H5 hip NR');
[H5_stats_hip_NR.Zp H5_stats_hip_NR.Z] = circ_rtest(H5_hip_NR_m(goodp)');
[H5_stats_hip_NR.Vp H5_stats_hip_NR.V] = circ_vtest(H5_hip_NR_m(goodp)',pi);
H5_stats_hip_NR.df=length(goodp);
hip_NR_r=[hip_NR_r H5_hip_NR_r];
hip_NR_m=[hip_NR_m circ_dist(H5_hip_NR_m(goodp),H5_hip_exc_mean)];


H1_hip_exc = H1.hippo.exc_index;%[10 11 64:66 104:107 139 141 174 313 234 281 311 312 67 108 109 146 111 144 316 179 181 317];
H1_hip_inh = H1.hippo.inh_index;%[314];
H1_hip_NR = H1.hippo.NR_index;%setdiff(hipdex,[H1_hip_exc H1_hip_inh]);

phdists=H1.hippo.phasedist(:,H1_hip_exc);
for i=1:length(H1_hip_exc)
   H1_hip_exc_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H1_hip_exc_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H1_hip_exc_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end

figure(100); 
goodp=find(H1_hip_exc_p<.01);
subplot(3,4,4); %histogram(H1_hip_exc_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H1_hip_exc_m(goodp),[-pi:(2*pi/12):pi]+pi/12);
%subplot(3,4,7); histogram(H1_hip_exc_r,0:.1:1); title('H1 hip exc');
hip_exc_r=[hip_exc_r H1_hip_exc_r];
H1_hip_exc_mean = circ_mean(H1_hip_exc_m(goodp)');
[H1_stats_hip_exc.Zp H1_stats_hip_exc.Z] = circ_rtest(H1_hip_exc_m(goodp)');
[H1_stats_hip_exc.Vp H1_stats_hip_exc.V] = circ_vtest(H1_hip_exc_m(goodp)',pi);
H1_stats_hip_exc.df=length(goodp);
hip_exc_m=[hip_exc_m circ_dist(H1_hip_exc_m(goodp),H1_hip_exc_mean)];

phdists=H1.hippo.phasedist(:,H1_hip_inh);
for i=1:length(H1_hip_inh)
   H1_hip_inh_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H1_hip_inh_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H1_hip_inh_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H1_hip_inh_p<.01);
%subplot(3,4,7); %histogram(H1_hip_inh_m(goodp),-pi:(2*pi/12):pi);
%[stats_hip_inh.Zp stats_hip_inh.Z] = circ_rtest(H1_hip_inh_m(goodp)');
%[stats_hip_inh.Vp stats_hip_inh.V] = circ_vtest(H1_hip_inh_m(goodp)',pi);
%stats_hip_inh.df=length(goodp);
%subplot(3,2,4); histogram(H1_hip_inh_r,0:.1:1); title('H1 hip inh');
hip_inh_r=[hip_inh_r H1_hip_inh_r];
%hip_inh_m=[hip_inh_m circ_dist(H1_hip_inh_m(goodp),H1_hip_exc_mean)];

phdists=H1.hippo.phasedist(:,H1_hip_NR);
for i=1:length(H1_hip_NR)
   H1_hip_NR_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H1_hip_NR_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H1_hip_NR_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H1_hip_NR_p<.01);
%subplot(3,4,7); %histogram(H1_hip_NR_m(goodp),-pi:(2*pi/12):pi);
%[stats_hip_NR.Zp stats_hip_NR.Z] = circ_rtest(H1_hip_NR_m(goodp)');
%[stats_hip_NR.Vp stats_hip_NR.V] = circ_vtest(H1_hip_NR_m(goodp)',pi);
%stats_hip_NR.df=length(goodp);
%subplot(3,2,4); histogram(H1_hip_NR_r,0:.1:1); title('H1 hip NR');
hip_NR_r=[hip_NR_r H1_hip_NR_r];
%hip_NR_m=[hip_NR_m circ_dist(H1_hip_NR_m(goodp),H1_hip_exc_mean)];

%% --------------------------------------------------------------------------------------------------------
% load h7results;
% load newidstring7;
% for i=1:length(newidstring)
%     temp=newidstring{i};
%     if ~isempty(temp);
%     kk=strfind(temp,'tt');
%     ll=strfind(temp(kk:end),'_');
%     ll=ll(1)+kk-2;
%     tt=str2num(temp((kk+2):ll));
%     c=find(results(1,:)==i);
%     results(8,c)=tt;
%     end
% end
% r=7;
% hipdex=[];
% H7_idstr=[];
% for i=1:length(rat(r).tt)
%     if ~isempty(rat(r).tt(i))
%         for j=1:length(rat(r).tt(i).cell)
%             if ~isempty(rat(r).tt(i).cell(j).sess_clust)
%                 found=0; crow=1;
%                 while ~found & crow<=size(rat(r).tt(i).cell(j).sess_clust,1)
%                 temp=find(results(10,:)==rat(r).tt(i).cell(j).sess_clust(crow,1) & results(9,:)==rat(r).tt(i).cell(j).sess_clust(crow,2)  & results(2,:)==i);
%                 if ~isempty(temp)
%                     if ~isempty(newidstring(results(1,temp)))
%                         H7_idstr=[H7_idstr; newidstring(results(1,temp))];
%                     end
%                     hipdex=[hipdex; temp(1)];
%                     found=1;
%                 else
%                     crow=crow+1;
%                 end
%                 end
%             else
%                 [i j]
%             end
%         end
%     end
% end

% H7_hip_exc = [104 223 158 308 360 681 684 1345 1472 1473 1619 1618 234 594 925 1761:1763 1023 1027 1026 1029 1414 1478 1477 440 514 441 513 515 597 598 752 753 931 934 1092 1186 1416 1479 1546 1485 1420 1488 1480 1690 1481 1482 1483 1484 1490 1545 1624 1622 1623 1625 755 1037 1096 1035 1188 1189 1554 1627 1423 1555 843 939 1039 1038 1192 1558 1684 1560 1624 1683 1630 1684];
% H7_hip_inh = [156 748 233 1024 510 756 1557 1631];
% H7_hip_NR = setdiff(hipdex,[H7_hip_exc H7_hip_inh]);
% H7_hip_NR=H7_hip_NR(1:42);




%%%% ----------------------------------------------------------------- all rats


figure(11); clf;
subplot(3,1,3);
histogram(log10(hip_NR_r),-4:.25:1);%;,-3:.2:0); 
title('hip NR r');
set(gca,'XLim',[-4 1]); title('hip exc r');
subplot(3,1,1);
histogram(log10(hip_exc_r),[-4:.25:1]-100);%,-2:.5:2)%,-3:.2:0); 
hold on;
histogram(log10(hip_exc_r),-4:.25:1);%,-2:.5:2)%,-3:.2:0); 
my=get(gca,'YLim');
set(gca,'YLim',[0 my(2)+1],'XLim',[-4 1]); title('hip exc r');
subplot(3,1,2);
histogram(log10(hip_exc_r),[-4:.25:1]-100);%,-2:.5:2)%,-3:.2:0); 
hold on;
histogram(log10(hip_exc_r),[-4:.25:1]-100);%,-2:.5:2)%,-3:.2:0); 
histogram(log10(hip_inh_r),-4:.25:1);%,-2:.5:2)%,-3:.2:0);
axis tight;
my=get(gca,'YLim');
set(gca,'YLim',[0 my(2)+1],'XLim',[-4 1]); title('hip inh r');

figure(12); clf;
subplot(2,1,1);
errorbar([1 2 3], [mean(log10(hip_exc_r)) mean(log10(hip_inh_r)) mean(log10(hip_NR_r))],[std(log10(hip_exc_r))/sqrt(length(hip_exc_r)) std(log10(hip_inh_r))/sqrt(length(hip_inh_r)) std(log10(hip_NR_r))/sqrt(length(hip_NR_r))],'o-');
set(gca,'YLim',[-1.5 -.5]);
subplot(2,1,2);
errorbar([1 2 3], [mean(log10(H7_hip_exc_r)) mean(log10(H7_hip_inh_r)) mean(log10(H7_hip_NR_r))],[std(log10(H7_hip_exc_r))/sqrt(length(H7_hip_exc_r)) std(log10(H7_hip_inh_r))/sqrt(length(H7_hip_inh_r)) std(log10(H7_hip_NR_r))/sqrt(length(H7_hip_NR_r))],'o-');
hold on;
errorbar([1 2 3], [mean(log10(H5_hip_exc_r)) mean(log10(H5_hip_inh_r)) mean(log10(H5_hip_NR_r))],[std(log10(H5_hip_exc_r))/sqrt(length(H5_hip_exc_r)) std(log10(H5_hip_inh_r))/sqrt(length(H5_hip_inh_r)) std(log10(H5_hip_NR_r))/sqrt(length(H5_hip_NR_r))],'o-');
set(gca,'YLim',[-1.2 -.5]);


anovan([log10(hip_NR_r) log10(hip_exc_r) log10(hip_inh_r)]',[hip_NR_r*0+1 hip_exc_r*0+2 hip_inh_r*0+3]');
[h,p,ci,stats_hip_exc_v_inh]=ttest2(log10(hip_exc_r), log10(hip_inh_r)); stats_hip_exc_v_inh.p=p;
[h,p,ci,stats_hip_exc_v_NR]=ttest2(log10(hip_exc_r), log10(hip_NR_r)); stats_hip_exc_v_NR.p=p;
[h,p,ci,stats_hip_inh_v_NR]=ttest2(log10(hip_inh_r), log10(hip_NR_r)); stats_hip_inh_v_NR.p=p;

[length(hip_exc_m) length(hip_exc_r)]
length(hip_exc_m)/length(hip_exc_r)

[length(hip_inh_m) length(hip_inh_r)]
length(hip_inh_m)/length(hip_inh_r)

[length(hip_NR_m) length(hip_NR_r)]
length(hip_NR_m)/length(hip_NR_r)

[length(hip_exc_m) length(hip_exc_r)]
length(hip_exc_m)/length(hip_exc_r)

[length(hip_inh_m) length(hip_inh_r)]
length(hip_inh_m)/length(hip_inh_r)

[length(hip_NR_m) length(hip_NR_r)]
length(hip_NR_m)/length(hip_NR_r)


% figure(102); clf;
% 
% subplot(3,4,1); polarhistogram(hip_exc_m,[-pi:(2*pi/12):pi]+pi/12); 
% subplot(3,4,2); polarhistogram(hip_exc_m(find(SW_exctheta>350))+pi/12,[-pi:(2*pi/12):pi]); 
% subplot(3,4,3); polarhistogram(hip_exc_m(find(SW_exctheta<=350))+pi/12,[-pi:(2*pi/12):pi]);
[x_all y_all]=pol2cart(circ_mean(hip_exc_m'),circ_r(hip_exc_m'));
% [x_msn y_msn]=pol2cart(circ_mean(hip_exc_m(find(SW_exctheta>350))'),circ_r(hip_exc_m(find(SW_exctheta>350))'));
% [x_int y_int]=pol2cart(circ_mean(hip_exc_m(find(SW_exctheta<=350))'),circ_r(hip_exc_m(find(SW_exctheta<=350))'));

% subplot(3,4,4); compass([x_all x_msn x_int],[y_all y_msn y_int]);
% [circ_mean(hip_exc_m') circ_r(hip_exc_m'); circ_mean(hip_exc_m(find(SW_exctheta>350))') circ_r(hip_exc_m(find(SW_exctheta>350))'); circ_mean(hip_exc_m(find(SW_exctheta<=350))') circ_r(hip_exc_m(find(SW_exctheta<=350))')]

%subplot(3,4,4+1); polarhistogram(hip_inh_m,[-pi:(2*pi/12):pi]+pi/12); 
% subplot(3,4,5+1); polarhistogram(hip_inh_m(find(SW_inhtheta>350))+pi/12,[-pi:(2*pi/12):pi]); 
% subplot(3,4,6+1); polarhistogram(hip_inh_m(find(SW_inhtheta<=350))+pi/12,[-pi:(2*pi/12):pi]);
[x_all y_all]=pol2cart(circ_mean(hip_inh_m'),circ_r(hip_inh_m'));
% [x_msn y_msn]=pol2cart(circ_mean(hip_inh_m(find(SW_inhtheta>350))'),circ_r(hip_inh_m(find(SW_inhtheta>350))'));
% [x_int y_int]=pol2cart(circ_mean(hip_inh_m(find(SW_inhtheta<=350))'),circ_r(hip_inh_m(find(SW_inhtheta<=350))'));

% subplot(3,4,7+1); compass([x_all x_msn x_int],[y_all y_msn y_int]);
% [circ_mean(hip_inh_m') circ_r(hip_inh_m'); circ_mean(hip_inh_m(find(SW_inhtheta>350))') circ_r(hip_inh_m(find(SW_inhtheta>350))'); circ_mean(hip_inh_m(find(SW_inhtheta<=350))') circ_r(hip_inh_m(find(SW_inhtheta<=350))')]

%subplot(3,4,7+2); polarhistogram(hip_NR_m,[-pi:(2*pi/12):pi]+pi/12); 
% subplot(3,4,8+2); polarhistogram(hip_NR_m(find(SW_NRtheta>350))+pi/12,[-pi:(2*pi/12):pi]); 
% subplot(3,4,9+2); polarhistogram(hip_NR_m(find(SW_NRtheta<=350))+pi/12,[-pi:(2*pi/12):pi]); 
% [x_all y_all]=pol2cart(circ_mean(hip_NR_m'),circ_r(hip_NR_m'));
% [x_msn y_msn]=pol2cart(circ_mean(hip_NR_m(find(SW_NRtheta>350))'),circ_r(hip_NR_m(find(SW_NRtheta>350))'));
% [x_int y_int]=pol2cart(circ_mean(hip_NR_m(find(SW_NRtheta<=350))'),circ_r(hip_NR_m(find(SW_NRtheta<=350))'));
% subplot(3,4,10+2); compass([x_all x_msn x_int],[y_all y_msn y_int]);
% [circ_mean(hip_NR_m') circ_r(hip_NR_m'); circ_mean(hip_NR_m(find(SW_NRtheta>350))') circ_r(hip_NR_m(find(SW_NRtheta>350))'); circ_mean(hip_NR_m(find(SW_NRtheta<=350))') circ_r(hip_NR_m(find(SW_NRtheta<=350))')]
% 
% [length(hip_exc_m) length(hip_inh_m) length(hip_NR_m)]
% [length(hip_NR_m(find(SW_NRtheta>350))) length(hip_NR_m(find(SW_NRtheta<=350)))]
% [Vp V] = circ_rtest(hip_NR_m(find(SW_NRtheta>350))')
% 
% adat=hip_NR_m(find(SW_NRtheta<=350))';
% rad2deg(pi+pi/12+circ_mean(adat))
% [Vp V] = circ_vtest(pi+pi/12+adat,pi/2)
% 
% [Vp V] = circ_vtest(adat,circ_mean(adat))

% [pval table] = circ_wwtest(hip_NR_m(find(SW_NRtheta>350)),hip_NR_m(find(SW_NRtheta<=350)))
% [pval table] = circ_wwtest(hip_exc_m(find(SW_exctheta>350)),hip_exc_m(find(SW_exctheta<=350)))
% [pval table] = circ_wwtest(hip_NR_m,hip_exc_m)
% 
% rad2deg(pi+pi/12+circ_mean(hip_NR_m(find(SW_NRtheta>350))'))
% 
% [pval z] = circ_mtest(hip_inh_m(find(SW_inhtheta>-Inf)), 85.1)
%polarhistogram(hip_exc_m,[-pi:(2*pi/12):pi]+pi/12); 
%polarhistogram(H5_hip_NR_m(goodp),[-pi:(2*pi/12):pi]+pi/12);
figure(100); 
subplot(3,4,1); polarhistogram(circ_dist(hip_exc_m,pi),[-pi:(2*pi/12):pi]+pi/12); 
subplot(3,4,5); polarhistogram(circ_dist(hip_inh_m,pi),[-pi:(2*pi/12):pi]+pi/12);
subplot(3,4,9); polarhistogram(circ_dist(hip_NR_m,pi),[-pi:(2*pi/12):pi]+pi/12);
%subplot(3,4,7);polar([H7_hip_exc_mean 0 H5_hip_exc_mean 0 H1_hip_exc_mean],[1 0 1 0 1]);