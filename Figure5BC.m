% clear;
% 
%  load 'hip_means_revised';
% % load 'H1';
% % load 'H5';
% % load 'H7';
% % load 'H7bad';
% load 'C:\Users\Blair Lab\Documents\MATLAB\AndrewData\revisions\H7r';
% load 'C:\Users\Blair Lab\Documents\MATLAB\AndrewData\revisions\H5r';
% load 'C:\Users\Blair Lab\Documents\MATLAB\AndrewData\revisions\H1r';

septum_exc_r=[];
septum_exc_m=[];
striatum_exc_r=[];
striatum_exc_m=[];
septum_inh_r=[];
septum_inh_m=[];
striatum_inh_r=[];
striatum_inh_m=[];
septum_NR_r=[];
septum_NR_m=[];
striatum_NR_r=[];
striatum_NR_m=[];

% - INDICES
% badseptum=[];
% for i=1:length(H7bad.septum.idstring)
%     for j=1:length(H7.septum.idstring)
%         if strcmp(H7bad.septum.idstring{i},H7.septum.idstring{j})
%             badseptum=[badseptum j];
%         end
%     end
% end
% badstriatum=[];
% for i=1:length(H7bad.striatum.idstring)
%     for j=1:length(H7.striatum.idstring)
%         if strcmp(H7bad.striatum.idstring{i},H7.striatum.idstring{j})
%             badstriatum=[badstriatum j];
%         end
%     end
% end
H7_septum_exc = H7.septum.exc_index;
H7_septum_inh = H7.septum.inh_index;
H7_septum_NR = H7.septum.NR_index;

H5_septum_exc = H5.septum.exc_index;
H5_septum_inh = H5.septum.inh_index;
H5_septum_NR = H5.septum.NR_index;

H7_striatum_exc = H7.striatum.exc_index;
H7_striatum_inh = H7.striatum.inh_index;
H7_striatum_NR = H7.striatum.NR_index;

H5_striatum_exc = H5.striatum.exc_index;
H5_striatum_inh = H5.striatum.inh_index;
H5_striatum_NR = H5.striatum.NR_index;

H1_striatum_exc = H1.striatum.exc_index;
H1_striatum_inh = H1.striatum.inh_index;
H1_striatum_NR = H1.striatum.NR_index;
%%%% ----------------------------------------------------------------- striatum H7

phdists=H7.striatum.phasedist(:,H7_striatum_exc);
for i=1:length(H7_striatum_exc)
   H7_striatum_exc_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H7_striatum_exc_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H7_striatum_exc_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
figure(100); clf; 
goodp=find(H7_striatum_exc_p<.01);
subplot(3,4,2); %histogram(H7_striatum_exc_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H7_striatum_exc_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_striatum_exc.H7_Zp stats_striatum_exc.H7_Z] = circ_rtest(H7_striatum_exc_m(goodp)');
[stats_striatum_exc.H7_Vp stats_striatum_exc.H7_V] = circ_vtest(H7_striatum_exc_m(goodp)',0);
stats_striatum_exc.H7_df=length(goodp);
H7_striatum_exc_mean = circ_mean(H7_striatum_exc_m(goodp)');
%subplot(3,4,2); histogram(H7_striatum_exc_r,0:.1:1); title('H7 striatum exc');
striatum_exc_r=[striatum_exc_r H7_striatum_exc_r];
striatum_exc_m=[striatum_exc_m circ_dist(H7_striatum_exc_m(goodp),H7_hip_exc_mean)];

phdists=H7.striatum.phasedist(:,H7_striatum_inh);
for i=1:length(H7_striatum_inh)
   H7_striatum_inh_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H7_striatum_inh_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H7_striatum_inh_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H7_striatum_inh_p<.01);
subplot(3,4,6); %histogram(H7_striatum_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H7_striatum_inh_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_striatum_inh.H7_Zp stats_striatum_inh.H7_Z] = circ_rtest(H7_striatum_inh_m(goodp)');
[stats_striatum_inh.H7_Vp stats_striatum_inh.H7_V] = circ_vtest(H7_striatum_inh_m(goodp)',pi);
stats_striatum_inh.H7_df=length(goodp);
H7_striatum_inh_mean = circ_mean(H7_striatum_inh_m(goodp)');
%subplot(3,4,4); histogram(H7_striatum_inh_r,0:.1:1); title('H7 striatum inh');
striatum_inh_r=[striatum_inh_r H7_striatum_inh_r];
striatum_inh_m=[striatum_inh_m circ_dist(H7_striatum_inh_m(goodp),H7_hip_exc_mean)];

phdists=H7.striatum.phasedist(:,H7_striatum_NR);
for i=1:length(H7_striatum_NR)
   H7_striatum_NR_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H7_striatum_NR_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H7_striatum_NR_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H7_striatum_NR_p<.01);
subplot(3,4,10); %histogram(H7_striatum_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H7_striatum_NR_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_striatum_NR.H7_Zp stats_striatum_NR.H7_Z] = circ_rtest(H7_striatum_NR_m(goodp)');
[stats_striatum_NR.H7_Vp stats_striatum_NR.H7_V] = circ_vtest(H7_striatum_NR_m(goodp)',0);
stats_striatum_NR.H7_df=length(goodp);
H7_striatum_NR_mean = circ_mean(H7_striatum_NR_m(goodp)');
%subplot(3,4,6); histogram(H7_striatum_NR_r,0:.1:1); title('H7 striatum NR');
striatum_NR_r=[striatum_NR_r H7_striatum_NR_r];
striatum_NR_m=[striatum_NR_m circ_dist(H7_striatum_NR_m(goodp),H7_hip_exc_mean)];

%%%% ----------------------------------------------------------------- striatum H5

figure(100);

phdists=H5.striatum.phasedist(:,H5_striatum_exc);
for i=1:length(H5_striatum_exc)
   H5_striatum_exc_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H5_striatum_exc_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H5_striatum_exc_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H5_striatum_exc_p<.01);
subplot(3,4,3); %histogram(H5_striatum_exc_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H5_striatum_exc_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_striatum_exc.H5_Zp stats_striatum_exc.H5_Z] = circ_rtest(H5_striatum_exc_m(goodp)');
[stats_striatum_exc.H5_Vp stats_striatum_exc.H5_V] = circ_vtest(H5_striatum_exc_m(goodp)',pi);
stats_striatum_exc.H5_df=length(goodp);
H5_striatum_exc_mean = circ_mean(H5_striatum_exc_m(goodp)');
%subplot(3,4,2); histogram(H5_striatum_exc_r,0:.1:1); title('H5 striatum exc');
striatum_exc_r=[striatum_exc_r H5_striatum_exc_r];
striatum_exc_m=[striatum_exc_m circ_dist(H5_striatum_exc_m(goodp),H5_hip_exc_mean)];

phdists=H5.striatum.phasedist(:,H5_striatum_inh);
for i=1:length(H5_striatum_inh)
   H5_striatum_inh_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H5_striatum_inh_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H5_striatum_inh_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H5_striatum_inh_p<.01);
subplot(3,4,7); %histogram(H5_striatum_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H5_striatum_inh_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_striatum_inh.H5_Zp stats_striatum_inh.H5_Z] = circ_rtest(H5_striatum_inh_m(goodp)');
[stats_striatum_inh.H5_Vp stats_striatum_inh.H5_V] = circ_vtest(H5_striatum_inh_m(goodp)',0);
stats_striatum_inh.H5_df=length(goodp);
H5_striatum_inh_mean = circ_mean(H5_striatum_inh_m(goodp)');
%subplot(3,4,4); histogram(H5_striatum_inh_r,0:.1:1); title('H5 striatum inh');
striatum_inh_r=[striatum_inh_r H5_striatum_inh_r];
striatum_inh_m=[striatum_inh_m circ_dist(H5_striatum_inh_m(goodp),H5_hip_exc_mean)];

phdists=H5.striatum.phasedist(:,H5_striatum_NR);
for i=1:length(H5_striatum_NR)
   H5_striatum_NR_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H5_striatum_NR_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H5_striatum_NR_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H5_striatum_NR_p<.01);
subplot(3,4,11); %histogram(H5_striatum_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H5_striatum_NR_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_striatum_NR.H5_Zp stats_striatum_NR.H5_Z] = circ_rtest(H5_striatum_NR_m(goodp)');
[stats_striatum_NR.H5_Vp stats_striatum_NR.H5_V] = circ_vtest(H5_striatum_NR_m(goodp)',pi);
stats_striatum_NR.H5_df=length(goodp);
H5_striatum_NR_mean = circ_mean(H5_striatum_NR_m(goodp)');
%subplot(3,4,6); histogram(H5_striatum_NR_r,0:.1:1); title('H5 striatum NR');
striatum_NR_r=[striatum_NR_r H5_striatum_NR_r];
striatum_NR_m=[striatum_NR_m circ_dist(H5_striatum_NR_m(goodp),H5_hip_exc_mean)];

%%%% ----------------------------------------------------------------- striatum H1

figure(100); 

phdists=H1.striatum.phasedist(:,H1_striatum_exc);
for i=1:length(H1_striatum_exc)
   H1_striatum_exc_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H1_striatum_exc_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H1_striatum_exc_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H1_striatum_exc_p<.01);
subplot(3,4,4); %histogram(H1_striatum_exc_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H1_striatum_exc_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_striatum_exc.H1_Zp stats_striatum_exc.H1_Z] = circ_rtest(H1_striatum_exc_m(goodp)');
[stats_striatum_exc.H1_Vp stats_striatum_exc.H1_V] = circ_vtest(H1_striatum_exc_m(goodp)',pi);
stats_striatum_exc.H1_df=length(goodp);
H1_striatum_exc_mean = circ_mean(H1_striatum_exc_m(goodp)');
%subplot(3,4,2); histogram(H1_striatum_exc_r,0:.1:1); title('H1 striatum exc');
striatum_exc_r=[striatum_exc_r H1_striatum_exc_r];
striatum_exc_m=[striatum_exc_m circ_dist(H1_striatum_exc_m(goodp),H1_hip_exc_mean)];

phdists=H1.striatum.phasedist(:,H1_striatum_NR);
for i=1:length(H1_striatum_NR)
   H1_striatum_NR_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H1_striatum_NR_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H1_striatum_NR_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H1_striatum_NR_p<.01);
subplot(3,4,12); %histogram(H1_striatum_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H1_striatum_NR_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_striatum_NR.H1_Zp stats_striatum_NR.H1_Z] = circ_rtest(H1_striatum_NR_m(goodp)');
[stats_striatum_NR.H1_Vp stats_striatum_NR.H1_V] = circ_vtest(H1_striatum_NR_m(goodp)',pi);
stats_striatum_NR.H1_df=length(goodp);
H1_striatum_NR_mean = circ_mean(H1_striatum_NR_m(goodp)');
%subplot(3,4,6); histogram(H1_striatum_NR_r,0:.1:1); title('H1 striatum NR');
striatum_NR_r=[striatum_NR_r H1_striatum_NR_r];
striatum_NR_m=[striatum_NR_m circ_dist(H1_striatum_NR_m(goodp),H1_hip_exc_mean)];

figure(10); clf;
subplot(3,1,3);
histogram(log10(striatum_NR_r),-4:.25:1);%;,-3:.2:0); 
title('striatum NR r');
subplot(3,1,1);
histogram(log10(striatum_exc_r),[-4:.25:1]-100);%,-2:.5:2)%,-3:.2:0); 
hold on;
histogram(log10(striatum_exc_r),-4:.25:1);%,-2:.5:2)%,-3:.2:0); 
my=get(gca,'YLim');
set(gca,'YLim',[0 my(2)+1],'XLim',[-4 1]); title('striatum exc r');
subplot(3,1,2);
histogram(log10(striatum_exc_r),[-4:.25:1]-100);%,-2:.5:2)%,-3:.2:0); 
hold on;
histogram(log10(striatum_exc_r),[-4:.25:1]-100);%,-2:.5:2)%,-3:.2:0); 
histogram(log10(striatum_inh_r),-4:.25:1);%,-2:.5:2)%,-3:.2:0);
axis tight;
my=get(gca,'YLim');
set(gca,'YLim',[0 my(2)+1],'XLim',[-4 1]); title('striatum inh r');

anovan([log10(striatum_NR_r) log10(striatum_exc_r) log10(striatum_inh_r)]',[striatum_NR_r*0+1 striatum_exc_r*0+2 striatum_inh_r*0+3]');
[h,p,ci,stats_striatum_exc_v_inh]=ttest2(log10(striatum_exc_r), log10(striatum_inh_r)); stats_striatum_exc_v_inh.p=p;
[h,p,ci,stats_striatum_exc_v_NR]=ttest2(log10(striatum_exc_r), log10(striatum_NR_r)); stats_striatum_exc_v_NR.p=p;
[h,p,ci,stats_striatum_inh_v_NR]=ttest2(log10(striatum_inh_r), log10(striatum_NR_r)); stats_striatum_inh_v_NR.p=p;


%%%% ----------------------------------------------------------------- septum H7

figure(101); clf; 

phdists=H7.septum.phasedist(:,H7_septum_exc);
for i=1:length(H7_septum_exc)
   H7_septum_exc_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H7_septum_exc_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H7_septum_exc_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H7_septum_exc_p<.01);
subplot(3,4,2); %histogram(H7_septum_exc_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H7_septum_exc_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_septum_exc.H7_Zp stats_septum_exc.H7_Z] = circ_rtest(H7_septum_exc_m(goodp)');
[stats_septum_exc.H7_Vp stats_septum_exc.H7_V] = circ_vtest(H7_septum_exc_m(goodp)',pi/2);
stats_septum_exc.H7_df=length(goodp);
H7_septum_exc_mean = circ_mean(H7_septum_exc_m(goodp)');
%subplot(3,4,2); histogram(H7_septum_exc_r,0:.1:1); title('H7 septum exc');
septum_exc_r=[septum_exc_r H7_septum_exc_r];
septum_exc_m=[septum_exc_m circ_dist(H7_septum_exc_m(goodp),H7_hip_exc_mean)];

phdists=H7.septum.phasedist(:,H7_septum_inh);
for i=1:length(H7_septum_inh)
   H7_septum_inh_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H7_septum_inh_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H7_septum_inh_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H7_septum_inh_p<.01);
subplot(3,4,6); %histogram(H7_septum_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H7_septum_inh_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_septum_inh.H7_Zp stats_septum_inh.H7_Z] = circ_rtest(H7_septum_inh_m(goodp)');
[stats_septum_inh.H7_Vp stats_septum_inh.H7_V] = circ_vtest(H7_septum_inh_m(goodp)',-pi/2);
stats_septum_inh.H7_df=length(goodp);
H7_septum_inh_mean = circ_mean(H7_septum_inh_m(goodp)');
%subplot(3,4,4); histogram(H7_septum_inh_r,0:.1:1); title('H7 septum inh');
septum_inh_r=[septum_inh_r H7_septum_inh_r];
septum_inh_m=[septum_inh_m circ_dist(H7_septum_inh_m(goodp),H7_hip_exc_mean)];

phdists=H7.septum.phasedist(:,H7_septum_NR);
for i=1:length(H7_septum_NR)
   H7_septum_NR_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H7_septum_NR_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H7_septum_NR_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H7_septum_NR_p<.01);
subplot(3,4,10); %histogram(H7_septum_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H7_septum_NR_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_septum_NR.H7_Zp stats_septum_NR.H7_Z] = circ_rtest(H7_septum_NR_m(goodp)');
[stats_septum_NR.H7_Vp stats_septum_NR.H7_V] = circ_vtest(H7_septum_NR_m(goodp)',pi/2);
stats_septum_NR.H7_df=length(goodp);
H7_septum_NR_mean = circ_mean(H7_septum_NR_m(goodp)');
%subplot(3,4,6); histogram(H7_septum_NR_r,0:.1:1); title('H7 septum NR');
septum_NR_r=[septum_NR_r H7_septum_NR_r];
septum_NR_m=[septum_NR_m circ_dist(H7_septum_NR_m(goodp),H7_hip_exc_mean)];

%%%% ----------------------------------------------------------------- septum H5

%figure(55); clf; 

phdists=H5.septum.phasedist(:,H5_septum_exc);
for i=1:length(H5_septum_exc)
   H5_septum_exc_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H5_septum_exc_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H5_septum_exc_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H5_septum_exc_p<.01);
subplot(3,4,3); %histogram(H5_septum_exc_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H5_septum_exc_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_septum_exc.H5_Zp stats_septum_exc.H5_Z] = circ_rtest(H5_septum_exc_m(goodp)');
[stats_septum_exc.H5_Vp stats_septum_exc.H5_V] = circ_vtest(H5_septum_exc_m(goodp)',-pi/2+.6);
stats_septum_exc.H5_df=length(goodp);
H5_septum_exc_mean = circ_mean(H5_septum_exc_m(goodp)');
%subplot(3,4,2); histogram(H5_septum_exc_r,0:.1:1); title('H5 septum exc');
septum_exc_r=[septum_exc_r H5_septum_exc_r];
septum_exc_m=[septum_exc_m circ_dist(H5_septum_exc_m(goodp),H5_hip_exc_mean)];

phdists=H5.septum.phasedist(:,H5_septum_inh);
for i=1:length(H5_septum_inh)
   H5_septum_inh_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H5_septum_inh_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H5_septum_inh_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H5_septum_inh_p<.01);
subplot(3,4,7); %histogram(H5_septum_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H5_septum_inh_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_septum_inh.H5_Zp stats_septum_inh.H5_Z] = circ_rtest(H5_septum_inh_m(goodp)');
[stats_septum_inh.H5_Vp stats_septum_inh.H5_V] = circ_vtest(H5_septum_inh_m(goodp)',pi/2);
stats_septum_inh.H5_df=length(goodp);
H5_septum_inh_mean = circ_mean(H5_septum_inh_m(goodp)');
%subplot(3,4,4); histogram(H5_septum_inh_r,0:.1:1); title('H5 septum inh');
septum_inh_r=[septum_inh_r H5_septum_inh_r];
septum_inh_m=[septum_inh_m circ_dist(H5_septum_inh_m(goodp),H5_hip_exc_mean)];

phdists=H5.septum.phasedist(:,H5_septum_NR);
for i=1:length(H5_septum_NR)
   H5_septum_NR_r(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
   H5_septum_NR_p(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,phdists(:,i));
   H5_septum_NR_m(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',phdists(:,i));
end
goodp=find(H5_septum_NR_p<.01);
subplot(3,4,11); %histogram(H5_septum_inh_m(goodp),-pi:(2*pi/12):pi);
polarhistogram(H5_septum_NR_m(goodp),[-pi:(2*pi/12):pi]+pi/12); title([num2str(length(goodp)) '-' num2str(size(phdists,2))]);
[stats_septum_NR.H5_Zp stats_septum_NR.H5_Z] = circ_rtest(H5_septum_NR_m(goodp)');
[stats_septum_NR.H5_Vp stats_septum_NR.H5_V] = circ_vtest(H5_septum_NR_m(goodp)',-pi/2);
stats_septum_NR.H5_df=length(goodp);
H5_septum_NR_mean = circ_mean(H5_septum_NR_m(goodp)');
%subplot(3,4,6); histogram(H5_septum_NR_r,0:.1:1); title('H5 septum NR');
septum_NR_r=[septum_NR_r H5_septum_NR_r];
septum_NR_m=[septum_NR_m circ_dist(H5_septum_NR_m(goodp),H5_hip_exc_mean)];


figure(11); clf;
subplot(3,1,3);
histogram(log10(septum_NR_r),-4:.25:1);%;,-3:.2:0); 
title('septum NR r');
subplot(3,1,1);
histogram(log10(septum_exc_r),[-4:.25:1]-100);%,-2:.5:2)%,-3:.2:0); 
hold on;
histogram(log10(septum_exc_r),-4:.25:1);%,-2:.5:2)%,-3:.2:0); 
my=get(gca,'YLim');
set(gca,'YLim',[0 my(2)+1],'XLim',[-4 1]); title('septum exc r');
subplot(3,1,2);
histogram(log10(septum_exc_r),[-4:.25:1]-100);%,-2:.5:2)%,-3:.2:0); 
hold on;
histogram(log10(septum_exc_r),[-4:.25:1]-100);%,-2:.5:2)%,-3:.2:0); 
histogram(log10(septum_inh_r),-4:.25:1);%,-2:.5:2)%,-3:.2:0);
axis tight;
my=get(gca,'YLim');
set(gca,'YLim',[0 my(2)+1],'XLim',[-4 1]); title('septum inh r');

figure(12); clf;
subplot(2,1,1);
errorbar([1 2 3], [mean(log10(striatum_exc_r)) mean(log10(striatum_inh_r)) nanmean(log10(striatum_NR_r))],[std(log10(striatum_exc_r))/sqrt(length(striatum_exc_r)) std(log10(striatum_inh_r))/sqrt(length(striatum_inh_r)) nanstd(log10(striatum_NR_r))/sqrt(length(striatum_NR_r))],'o-');
set(gca,'YLim',[-1.5 -.5]);
subplot(2,1,2);
errorbar([1 2 3], [mean(log10(H7_striatum_exc_r)) mean(log10(H7_striatum_inh_r)) mean(log10(H7_striatum_NR_r))],[std(log10(H7_striatum_exc_r))/sqrt(length(H7_striatum_exc_r)) std(log10(H7_striatum_inh_r))/sqrt(length(H7_striatum_inh_r)) std(log10(H7_striatum_NR_r))/sqrt(length(H7_striatum_NR_r))],'o-');
hold on;
errorbar([1 2 3], [mean(log10(H5_striatum_exc_r)) mean(log10(H5_striatum_inh_r)) mean(log10(H5_striatum_NR_r))],[std(log10(H5_striatum_exc_r))/sqrt(length(H5_striatum_exc_r)) std(log10(H5_striatum_inh_r))/sqrt(length(H5_striatum_inh_r)) std(log10(H5_striatum_NR_r))/sqrt(length(H5_striatum_NR_r))],'o-');
errorbar([1 2 3], [mean(log10(H1_striatum_exc_r)) NaN mean(log10(H1_striatum_NR_r))],[std(log10(H1_striatum_exc_r))/sqrt(length(H1_striatum_exc_r)) NaN std(log10(H1_striatum_NR_r))/sqrt(length(H1_striatum_NR_r))],'o-');
set(gca,'YLim',[-1.5 -.5]);

figure(13); clf;
subplot(2,1,1);
errorbar([1 2 3], [mean(log10(septum_exc_r)) mean(log10(septum_inh_r)) mean(log10(septum_NR_r))],[std(log10(septum_exc_r))/sqrt(length(septum_exc_r)) std(log10(septum_inh_r))/sqrt(length(septum_inh_r)) std(log10(septum_NR_r))/sqrt(length(septum_NR_r))],'o-');
set(gca,'YLim',[-1.5 -.5]);
subplot(2,1,2);
errorbar([1 2 3], [mean(log10(H7_septum_exc_r)) mean(log10(H7_septum_inh_r)) mean(log10(H7_septum_NR_r))],[std(log10(H7_septum_exc_r))/sqrt(length(H7_septum_exc_r)) std(log10(H7_septum_inh_r))/sqrt(length(H7_septum_inh_r)) std(log10(H7_septum_NR_r))/sqrt(length(H7_septum_NR_r))],'o-');
hold on;
errorbar([1 2 3], [mean(log10(H5_septum_exc_r)) mean(log10(H5_septum_inh_r)) mean(log10(H5_septum_NR_r))],[std(log10(H5_septum_exc_r))/sqrt(length(H5_septum_exc_r)) std(log10(H5_septum_inh_r))/sqrt(length(H5_septum_inh_r)) std(log10(H5_septum_NR_r))/sqrt(length(H5_septum_NR_r))],'o-');
set(gca,'YLim',[-1.5 -.5]);

anovan([log10(septum_NR_r) log10(septum_exc_r) log10(septum_inh_r)]',[septum_NR_r*0+1 septum_exc_r*0+2 septum_inh_r*0+3]');
[h,p,ci,stats_septum_exc_v_inh]=ttest2(log10(septum_exc_r), log10(septum_inh_r)); stats_septum_exc_v_inh.p=p;
[h,p,ci,stats_septum_exc_v_NR]=ttest2(log10(septum_exc_r), log10(septum_NR_r)); stats_septum_exc_v_NR.p=p;
[h,p,ci,stats_septum_inh_v_NR]=ttest2(log10(septum_inh_r), log10(septum_NR_r)); stats_septum_inh_v_NR.p=p;

[length(striatum_exc_m) length(striatum_exc_r)]
length(striatum_exc_m)/length(striatum_exc_r)

[length(striatum_inh_m) length(striatum_inh_r)]
length(striatum_inh_m)/length(striatum_inh_r)

[length(striatum_NR_m) length(striatum_NR_r)]
length(striatum_NR_m)/length(striatum_NR_r)

[length(septum_exc_m) length(septum_exc_r)]
length(septum_exc_m)/length(septum_exc_r)

[length(septum_inh_m) length(septum_inh_r)]
length(septum_inh_m)/length(septum_inh_r)

[length(septum_NR_m) length(septum_NR_r)]
length(septum_NR_m)/length(septum_NR_r)

septum_all_r=[septum_exc_r septum_inh_r septum_NR_r];
striatum_all_r=[striatum_exc_r striatum_inh_r striatum_NR_r];
%load 'hip_r' hip_all_r hip_exc_r hip_inh_r hip_NR_r;
%anovan([log10(hip_all_r) log10(septum_all_r) log10(striatum_all_r)]',[hip_all_r*0+1 septum_all_r*0+2 striatum_all_r*0+3]');
%[h,p,ci,stats_hip_v_septum]=ttest2(log10(hip_all_r), log10(septum_all_r)); stats_hip_v_septum.p=p;
%[h,p,ci,stats_hip_v_striatum]=ttest2(log10(hip_all_r), log10(striatum_all_r)); stats_hip_v_striatum.p=p;
%[h,p,ci,stats_septum_v_striatum]=ttest2(log10(septum_all_r), log10(striatum_all_r)); stats_septum_v_striatum.p=p;

figure(100); 
subplot(3,4,1); polarhistogram(circ_dist(striatum_exc_m+pi/12,pi),[-pi:(2*pi/12):pi]+pi/12); 
subplot(3,4,5); polarhistogram(circ_dist(striatum_inh_m+pi/12,pi),[-pi:(2*pi/12):pi]+pi/12); 
subplot(3,4,9); polarhistogram(circ_dist(striatum_NR_m+pi/12,pi),[-pi:(2*pi/12):pi]+pi/12);
figure(101); 
subplot(3,4,1); polarhistogram(circ_dist(septum_exc_m+pi/12,pi),[-pi:(2*pi/12):pi]+pi/12); 
subplot(3,4,5); polarhistogram(circ_dist(septum_inh_m+pi/12,pi),[-pi:(2*pi/12):pi]+pi/12); 
subplot(3,4,9); polarhistogram(circ_dist(septum_NR_m+pi/12,pi),[-pi:(2*pi/12):pi]+pi/12);

