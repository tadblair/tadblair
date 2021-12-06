%% STRIATUM CELL TYPES


adbv=0.000000004577636718750000; %voltage conversion factor

%% NON RESPONSIVE TO SWRs 

RN_NR = [1+0*H7.septum.spkwid(H7.septum.NR_index)'; 2+0*H5.septum.spkwid(H5.septum.NR_index)'];
RN_exc = [1+0*H7.septum.spkwid(H7.septum.exc_index)'; 2+0*H5.septum.spkwid(H5.septum.exc_index)'];
RN_inh = [1+0*H7.septum.spkwid(H7.septum.inh_index)'; 2+0*H5.septum.spkwid(H5.septum.inh_index)'];
RN_all = [RN_NR; RN_exc; RN_inh];
RN_swr = [RN_exc; RN_inh];

SW_NR = [H7.septum.spkwid(H7.septum.NR_index)'; H5.septum.spkwid(H5.septum.NR_index)'];
SW_exc = [H7.septum.spkwid(H7.septum.exc_index)'; H5.septum.spkwid(H5.septum.exc_index)'];
SW_inh = [H7.septum.spkwid(H7.septum.inh_index)'; H5.septum.spkwid(H5.septum.inh_index)'];
SW_all = [SW_NR; SW_exc; SW_inh];
SW_swr = [SW_exc; SW_inh];

FR_NR = [H7.septum.FR(H7.septum.NR_index)'; H5.septum.FR(H5.septum.NR_index)']; %H1.septum.FR(H1.septum.NR_index)'];
FR_exc = [H7.septum.FR(H7.septum.exc_index)'; H5.septum.FR(H5.septum.exc_index)']; %H1.septum.FR(H1.septum.exc_index)'];
FR_inh = [H7.septum.FR(H7.septum.inh_index)'; H5.septum.FR(H5.septum.inh_index)']; %H1.septum.FR(H1.septum.inh_index)'];
FR_all = [FR_NR; FR_exc; FR_inh];
FR_swr = [FR_exc; FR_inh];

POP2_NR = [H7.septum.bucketPOP2(H7.septum.NR_index)'; H5.septum.bucketPOP2(H5.septum.NR_index)']; %H1.septum.bucketPOP2(H1.septum.NR_index)'];
POP2_exc = [H7.septum.bucketPOP2(H7.septum.exc_index)'; H5.septum.bucketPOP2(H5.septum.exc_index)']; %H1.septum.bucketPOP2(H1.septum.exc_index)'];
POP2_inh = [H7.septum.bucketPOP2(H7.septum.inh_index)'; H5.septum.bucketPOP2(H5.septum.inh_index)']; %H1.septum.bucketPOP2(H1.septum.inh_index)'];
POP2_all = [POP2_NR; POP2_exc; POP2_inh];

p_NR = [H7.septum.SWR_pvalue(H7.septum.NR_index)'; H5.septum.SWR_pvalue(H5.septum.NR_index)']; %H1.septum.SWR_pvalue(H1.septum.NR_index)'];
p_exc = [H7.septum.SWR_pvalue(H7.septum.exc_index)'; H5.septum.SWR_pvalue(H5.septum.exc_index)']; %H1.septum.SWR_pvalue(H1.septum.exc_index)'];
p_inh = [H7.septum.SWR_pvalue(H7.septum.inh_index)'; H5.septum.SWR_pvalue(H5.septum.inh_index)']; %H1.septum.SWR_pvalue(H1.septum.inh_index)'];
p_all = [p_NR; p_exc; p_inh];

%SWR_pvalue_strNR = [H7.septum.SWR_pvalue(H7.septum.NR_index); H5.septum.SWR_pvalue(H5.septum.NR_index)]; %H1.septum.SWR_pvalue(H1.septum.NR_index)];
strNRwave = [H7.septum.spkwave(H7.septum.NR_index,:); H5.septum.spkwave(H5.septum.NR_index,:)]; %H1.septum.spkwave(H1.septum.NR_index,:)];
strexcwave = [H7.septum.spkwave(H7.septum.exc_index,:); H5.septum.spkwave(H5.septum.exc_index,:)]; %H1.septum.spkwave(H1.septum.exc_index,:)];
strinhwave = [H7.septum.spkwave(H7.septum.inh_index,:); H5.septum.spkwave(H5.septum.inh_index,:)]; %H1.septum.spkwave(H1.septum.inh_index,:)];
strallwave = [strNRwave; strexcwave; strinhwave];

NRphase = [H7.septum.phasedist(:,H7.septum.NR_index)'; H5.septum.phasedist(:,H5.septum.NR_index)']; %H1.septum.phasedist(:,H1.septum.NR_index)'];
excphase = [H7.septum.phasedist(:,H7.septum.exc_index)'; H5.septum.phasedist(:,H5.septum.exc_index)']; %H1.septum.phasedist(:,H1.septum.exc_index)'];
inhphase = [H7.septum.phasedist(:,H7.septum.inh_index)'; H5.septum.phasedist(:,H5.septum.inh_index)']; %H1.septum.phasedist(:,H1.septum.inh_index)'];
allphase = [H7.septum.phasedist(:,H7.septum.ls_index)'; H5.septum.phasedist(:,H5.septum.ls_index)']; %H1.septum.phasedist(:,:)'];

clear NR_theta* exc_theta* inh_theta*;
for i=1:size(NRphase,1)
    NR_thetar(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',NRphase(i,:)');
    NR_thetam(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',NRphase(i,:)');
    NR_thetap(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,NRphase(i,:));
end
for i=1:size(excphase,1)
    exc_thetar(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',excphase(i,:)');
    exc_thetam(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',excphase(i,:)');
    exc_thetap(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,excphase(i,:));
end
for i=1:size(inhphase,1)
    inh_thetar(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',inhphase(i,:)');
    inh_thetam(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',inhphase(i,:)');
    inh_thetap(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,inhphase(i,:));
end
for i=1:size(allphase,1)
    all_thetar(i)=circ_r([(2*pi/36):(2*pi/36):2*pi]',allphase(i,:)');
    all_thetam(i)=circ_mean([(2*pi/36):(2*pi/36):2*pi]',allphase(i,:)');
    all_thetap(i)=circ_rtest((2*pi/36):(2*pi/36):2*pi,allphase(i,:));
end

figure(4); clf;
subplot(2,3,2); hold off;
RCI_dmsint_inh=log10(-log10(inh_thetap(find(SW_inh<=350))));
RCI_dmsint_inh(find(isinf(RCI_dmsint_inh)))=2.99;
temp=histogram(RCI_dmsint_inh,-3:.3:3); 
Y=temp.Values;
RCI_dmsint_exc=log10(-log10(exc_thetap(find(SW_exc<=350))));
RCI_dmsint_exc(find(isinf(RCI_dmsint_exc)))=2.99;
temp=histogram(RCI_dmsint_exc,-3:.3:3); 
Y=[Y; temp.Values];
RCI_dmsint_NR=log10(-log10(NR_thetap(find(SW_NR<=350))));
RCI_dmsint_NR(find(isinf(RCI_dmsint_NR)))=2.99;
temp=histogram(RCI_dmsint_NR,-3:.3:3); 
Y=[Y; temp.Values];
bar(temp.BinEdges(1:end-1)+.25/2,Y',1,'stacked'); axis tight;
title('INT');

subplot(2,3,5); hold off;
RCI_dmsmsn_inh=log10(-log10(inh_thetap(find(SW_inh>350))));
RCI_dmsmsn_inh(find(isinf(RCI_dmsmsn_inh)))=2.99;
temp=histogram(RCI_dmsmsn_inh,-3:.3:3); 
Y=temp.Values;
RCI_dmsmsn_exc=log10(-log10(exc_thetap(find(SW_exc>350))));
RCI_dmsmsn_exc(find(isinf(RCI_dmsmsn_exc)))=2.99;
temp=histogram(RCI_dmsmsn_exc,-3:.3:3); 
Y=[Y; temp.Values];
RCI_dmsmsn_NR=log10(-log10(NR_thetap(find(SW_NR>350))));
RCI_dmsmsn_NR(find(isinf(RCI_dmsmsn_NR)))=2.99;
temp=histogram(RCI_dmsmsn_NR,-3:.3:3); 
Y=[Y; temp.Values];
bar(temp.BinEdges(1:end-1)+.25/2,Y',1,'stacked'); axis tight;
title('MSN');
%title([num2str(median(FR_all(find(SW_all<350))))]);
%axis tight; set(gca,'XLim',[-2 2]);

RN_lsint_all=RN_all(find(SW_all<=350));
RN_lsmsn_all=RN_all(find(SW_all>350));

[sum(RN_lsint_all==1) sum(RN_lsint_all==2) sum(RN_lsint_all==3); sum(RN_lsmsn_all==1) sum(RN_lsmsn_all==2) sum(RN_lsmsn_all==3)]

%ordered by rat
RCI_dmsexc_all=log10(-log10(exc_thetap));%[RCI_dmsint_NR RCI_dmsint_exc RCI_dmsint_inh];
RCI_dmsinh_all=log10(-log10(inh_thetap));%[RCI_dmsint_NR RCI_dmsint_exc RCI_dmsint_inh];
RCI_dmsNR_all=log10(-log10(NR_thetap));%[RCI_dmsint_NR RCI_dmsint_exc RCI_dmsint_inh];
SW_exctheta=SW_exc(find(RCI_dmsexc_all>.3));
SW_inhtheta=SW_inh(find(RCI_dmsinh_all>.3));
SW_NRtheta=SW_NR(find(RCI_dmsNR_all>.3));


RCI_dmsint_all=[RCI_dmsint_NR RCI_dmsint_exc RCI_dmsint_inh];
RCI_dmsmsn_all=[RCI_dmsmsn_NR RCI_dmsmsn_exc RCI_dmsmsn_inh];
RCI_dms_all=[RCI_dmsint_all RCI_dmsmsn_all];

RCI_dmsinh_all=[RCI_dmsmsn_inh RCI_dmsint_inh];
RCI_dmsexc_all=[RCI_dmsmsn_exc RCI_dmsint_exc];
RCI_dmsNR_all=[RCI_dmsmsn_NR RCI_dmsint_NR];

RCI_dmsint_all=log10(-log10(all_thetap(find(SW_all<=350))));%[RCI_dmsint_NR RCI_dmsint_exc RCI_dmsint_inh];
logFR_dmsint_all=(log10(FR_all(find(SW_all<=350))));
[h,p,ci,stats]=ttest2(logFR_dmsint_all(find(RCI_dmsint_all>.3)), logFR_dmsint_all(find(RCI_dmsint_all<=.3)))

RCI_dmsmsn_all=log10(-log10(all_thetap(find(SW_all>350))));%[RCI_dmsint_NR RCI_dmsint_exc RCI_dmsint_inh];
logFR_dmsmsn_all=(log10(FR_all(find(SW_all>350))));
[h,p,ci,stats]=ttest2(logFR_dmsmsn_all(find(RCI_dmsmsn_all>.3)), logFR_dmsmsn_all(find(RCI_dmsmsn_all<=.3)))

