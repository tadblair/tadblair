%% STRIATUM CELL TYPES

%load rat data
load 'H1'; load 'H5'; load 'H7';

adbv=0.000000004577636718750000; %voltage conversion factor


%% NON RESPONSIVE TO SWRs 

RN_NR = [1+0*H7.hippo.spkwid(H7.hippo.NR_index)'; 2+0*H5.hippo.spkwid(H5.hippo.NR_index)'; 3+0*H1.hippo.spkwid(H1.hippo.NR_index)'];
RN_exc = [1+0*H7.hippo.spkwid(H7.hippo.exc_index)'; 2+0*H5.hippo.spkwid(H5.hippo.exc_index)'; 3+0*H1.hippo.spkwid(H1.hippo.exc_index)'];
RN_inh = [1+0*H7.hippo.spkwid(H7.hippo.inh_index)'; 2+0*H5.hippo.spkwid(H5.hippo.inh_index)'; 3+0*H1.hippo.spkwid(H1.hippo.inh_index)'];
RN_all = [RN_NR; RN_exc; RN_inh];
RN_swr = [RN_exc; RN_inh];

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

NRphase = [H7.hippo.phasedist(:,H7.hippo.NR_index)'; H5.hippo.phasedist(:,H5.hippo.NR_index)'; H1.hippo.phasedist(:,H1.hippo.NR_index)'];
excphase = [H7.hippo.phasedist(:,H7.hippo.exc_index)'; H5.hippo.phasedist(:,H5.hippo.exc_index)'; H1.hippo.phasedist(:,H1.hippo.exc_index)'];
inhphase = [H7.hippo.phasedist(:,H7.hippo.inh_index)'; H5.hippo.phasedist(:,H5.hippo.inh_index)'; H1.hippo.phasedist(:,H1.hippo.inh_index)'];
allphase = [H7.hippo.phasedist(:,:)'; H5.hippo.phasedist(:,:)'; H1.hippo.phasedist(:,:)'];

clear NR_theta* exc_theta* inh_theta* all_theta*;
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

figure(4);

subplot(2,3,1); hold off;
RCI_hipint_inh=log10(-log10(inh_thetap(find(SW_inh<=350))));
RCI_hipint_inh(find(isinf(RCI_hipint_inh)))=2.99;
temp=histogram(RCI_hipint_inh,-3:.3:3); 
Y=temp.Values;
RCI_hipint_exc=log10(-log10(exc_thetap(find(SW_exc<=350))));
RCI_hipint_exc(find(isinf(RCI_hipint_exc)))=2.99;
temp=histogram(RCI_hipint_exc,-3:.3:3); 
Y=[Y; temp.Values];
RCI_hipint_NR=log10(-log10(NR_thetap(find(SW_NR<=350))));
RCI_hipint_NR(find(isinf(RCI_hipint_NR)))=2.99;
temp=histogram(RCI_hipint_NR,-3:.3:3); 
Y=[Y; temp.Values];
bar(temp.BinEdges(1:end-1)+.25/2,Y',1,'stacked'); axis tight;
title('INT');

subplot(2,3,4); hold off;
RCI_hippyr_inh=log10(-log10(inh_thetap(find(SW_inh>350))));
RCI_hippyr_inh(find(isinf(RCI_hippyr_inh)))=2.99;
temp=histogram(RCI_hippyr_inh,-3:.3:3); 
Y=temp.Values;
RCI_hippyr_exc=log10(-log10(exc_thetap(find(SW_exc>350))));
RCI_hippyr_exc(find(isinf(RCI_hippyr_exc)))=2.99;
temp=histogram(RCI_hippyr_exc,-3:.3:3); 
Y=[Y; temp.Values];
RCI_hippyr_NR=log10(-log10(NR_thetap(find(SW_NR>350))));
RCI_hippyr_NR(find(isinf(RCI_hippyr_NR)))=2.99;
temp=histogram(RCI_hippyr_NR,-3:.3:3); 
Y=[Y; temp.Values];
bar(temp.BinEdges(1:end-1)+.25/2,Y',1,'stacked'); axis tight;
title('PYR');

%ordered by rat
RCI_dmsexc_all=log10(-log10(exc_thetap));%[RCI_dmsint_NR RCI_dmsint_exc RCI_dmsint_inh];
RCI_dmsinh_all=log10(-log10(inh_thetap));%[RCI_dmsint_NR RCI_dmsint_exc RCI_dmsint_inh];
RCI_dmsNR_all=log10(-log10(NR_thetap));%[RCI_dmsint_NR RCI_dmsint_exc RCI_dmsint_inh];
SW_exctheta=SW_exc(find(RCI_dmsexc_all>.3));
SW_inhtheta=SW_inh(find(RCI_dmsinh_all>.3));
SW_NRtheta=SW_NR(find(RCI_dmsNR_all>.3));

RCI_hipint_all=[RCI_hipint_NR RCI_hipint_exc RCI_hipint_inh];
RCI_hippyr_all=[RCI_hippyr_NR RCI_hippyr_exc RCI_hippyr_inh];
RCI_hip_all=[RCI_hipint_all RCI_hippyr_all];

%%%% THETA RESPONSIVENESS OF INTERNEURONS VS PYRAMIDAL CELLS

RCI_hipall=[log10(-log10(NR_thetap)) log10(-log10(exc_thetap)) log10(-log10(inh_thetap))];
RCI_hipall(find(isinf(RCI_hipall)))=2.99;
[h,p,ci,stats]=ttest2(RCI_hipall(find(SW_all<=350)), RCI_hipall(find(SW_all>350)))
[h,p,ci,stats]=ttest2(RCI_hipall(find(SW_all<=350 & RN_all==1)), RCI_hipall(find(SW_all>350 & RN_all==1)))
[h,p,ci,stats]=ttest2(RCI_hipall(find(SW_all<=350 & RN_all==2)), RCI_hipall(find(SW_all>350 & RN_all==2)))
[h,p,ci,stats]=ttest2(RCI_hipall(find(SW_all<=350 & RN_all==3)), RCI_hipall(find(SW_all>350 & RN_all==3)))

[length(RCI_hipall(find(SW_all'<=350 & RCI_hipall>.3 & RN_all'==1))) length(RCI_hipall(find(SW_all'<=350 & RCI_hipall<=.3 & RN_all'==1)))]
[length(RCI_hipall(find(SW_all'>350 & RCI_hipall>.3 & RN_all'==1))) length(RCI_hipall(find(SW_all'>350 & RCI_hipall<=.3 & RN_all'==1)))]

[length(RCI_hipall(find(SW_all'<=350 & RCI_hipall>.3 & RN_all'==2))) length(RCI_hipall(find(SW_all'<=350 & RCI_hipall<=.3 & RN_all'==2)))]
[length(RCI_hipall(find(SW_all'>350 & RCI_hipall>.3 & RN_all'==2))) length(RCI_hipall(find(SW_all'>350 & RCI_hipall<=.3 & RN_all'==2)))]

[length(RCI_hipall(find(SW_all'<=350 & RCI_hipall>.3 & RN_all'==3))) length(RCI_hipall(find(SW_all'<=350 & RCI_hipall<=.3 & RN_all'==3)))]
[length(RCI_hipall(find(SW_all'>350 & RCI_hipall>.3 & RN_all'==3))) length(RCI_hipall(find(SW_all'>350 & RCI_hipall<=.3 & RN_all'==3)))]

RCI_hipinh_all=[RCI_hippyr_inh RCI_hipint_inh];
RCI_hipexc_all=[RCI_hippyr_exc RCI_hipint_exc];
RCI_hipNR_all=[RCI_hippyr_NR RCI_hipint_NR];

RCI_hipint_all=log10(-log10(all_thetap(find(SW_all<=350))));%[RCI_hipint_NR RCI_hipint_exc RCI_hipint_inh];
logFR_hipint_all=(log10(FR_all(find(SW_all<=350))));
[h,p,ci,stats]=ttest2(logFR_hipint_all(find(RCI_hipint_all>.3)), logFR_hipint_all(find(RCI_hipint_all<=.3)))

RCI_hippyr_all=log10(-log10(all_thetap(find(SW_all>350))));%[RCI_hipint_NR RCI_hipint_exc RCI_hipint_inh];
logFR_hippyr_all=(log10(FR_all(find(SW_all>350))));
[h,p,ci,stats]=ttest2(logFR_hippyr_all(find(RCI_hippyr_all>.3)), logFR_hippyr_all(find(RCI_hippyr_all<=.3)))

