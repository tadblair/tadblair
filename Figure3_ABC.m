figure(); clf;

%load rat data
load 'H1r'; load 'H5r'; load 'H7r';

labels={'exc','inh','NR'};

%pie chart for CA1
subplot(3,1,1);
hippo_exc = length(H7.hippo.exc_index) + length(H5.hippo.exc_index) + length(H1.hippo.exc_index);
hippo_inh = length(H7.hippo.inh_index) + length(H5.hippo.inh_index) + length(H1.hippo.inh_index);
hippo_NR = length(H7.hippo.NR_index) + length(H5.hippo.NR_index) + length(H1.hippo.NR_index);
pie([hippo_exc hippo_inh hippo_NR],labels);
N=hippo_exc+hippo_inh+hippo_NR;
[hippo_exc hippo_inh hippo_NR N]
[hippo_exc/N hippo_inh/N hippo_NR/N]

[length(H7.hippo.exc_index)  length(H5.hippo.exc_index)  length(H1.hippo.exc_index)]
[length(H7.hippo.inh_index)  length(H5.hippo.inh_index)  length(H1.hippo.inh_index)]
[length(H7.hippo.NR_index)  length(H5.hippo.NR_index)  length(H1.hippo.NR_index)]


%pie chart for LS
subplot(3,1,2);
septum_exc = length(H7.septum.exc_index) + length(H5.septum.exc_index);
septum_inh = length(H7.septum.inh_index) + length(H5.septum.inh_index);
septum_NR = length(H7.septum.NR_index) + length(H5.septum.NR_index);
pie([septum_exc septum_inh septum_NR],labels);
N=septum_exc+septum_inh+septum_NR;
[septum_exc septum_inh septum_NR N]
[septum_exc/N septum_inh/N septum_NR/N]

[length(H7.septum.exc_index)  length(H5.septum.exc_index)]
[length(H7.septum.inh_index)  length(H5.septum.inh_index)]
[length(H7.septum.NR_index)  length(H5.septum.NR_index)]

subplot(3,1,3);
striatum_exc = length(H7.striatum.exc_index) + length(H5.striatum.exc_index) + length(H1.striatum.exc_index);
striatum_inh = length(H7.striatum.inh_index) + length(H5.striatum.inh_index) + length(H1.striatum.inh_index);
striatum_NR = length(H7.striatum.NR_index) + length(H5.striatum.NR_index) + length(H1.striatum.NR_index);
pie([striatum_exc striatum_inh striatum_NR],labels);
N=striatum_exc+striatum_inh+striatum_NR;
[striatum_exc striatum_inh striatum_NR N]
[striatum_exc/N striatum_inh/N striatum_NR/N]

figure(100); clf; 

%H7.hippo.SWR_psign=sign(on_diff)';

subplot(2,2,1); hold off;
pvals_sep=[log10([H7.septum.SWR_pvalue(ismember([1:133],H7.septum.ls_index) & H7.septum.SWR_psign<0) H5.septum.SWR_pvalue(H5.septum.SWR_psign<0) ]) -log10([H7.septum.SWR_pvalue(ismember([1:133],H7.septum.ls_index) & H7.septum.SWR_psign>=0) H5.septum.SWR_pvalue(H5.septum.SWR_psign>=0) ])];
pd=fitdist(pvals_sep','Normal');
pvals_sep(find(pvals_sep<-5))=-5;
pvals_sep(find(pvals_sep>5))=5;
histogram(pvals_sep(find(pvals_sep>-2 & pvals_sep<2)),-10:.5:10); hold on;
histogram(pvals_sep(find(pvals_sep<=-2)),-10:.5:10); 
histogram(pvals_sep(find(pvals_sep>=2)),-10:.5:10); 
ci = paramci(pd,'Alpha',.01)
[h1_sep, p1_sep]= lillietest(pvals_sep(find(abs(pvals_sep)<5)));
[h2_sep, p2_sep]= jbtest(pvals_sep(find(abs(pvals_sep)<5)));
[h3_sep, p3_sep]= adtest(pvals_sep(find(abs(pvals_sep)<5)));

subplot(2,2,2);
pvals_str=[log10([H7.striatum.SWR_pvalue(ismember([1:236],H7.striatum.dms_index) & H7.striatum.SWR_psign<0) H5.striatum.SWR_pvalue(H5.striatum.SWR_psign<0) H1.striatum.SWR_pvalue(H1.striatum.SWR_psign<0)]) -log10([ H7.striatum.SWR_pvalue(ismember([1:236],H7.striatum.dms_index) & H7.striatum.SWR_psign>=0) H5.striatum.SWR_pvalue(H5.striatum.SWR_psign>=0) H1.striatum.SWR_pvalue(H1.striatum.SWR_psign>=0)])];
pd=fitdist(pvals_str','Normal')
pvals_str(find(pvals_str<-5))=-5;
pvals_str(find(pvals_str>5))=5;
histogram(pvals_str(find(pvals_str>-2 & pvals_str<2)),-10:.5:10); hold on;
histogram(pvals_str(find(pvals_str<=-2)),-10:.5:10); 
histogram(pvals_str(find(pvals_str>=2)),-10:.5:10); 
ci = paramci(pd,'Alpha',.01)
h1_str = lillietest(pvals_str);
h2_str = jbtest(pvals_str);
h3_str = adtest(pvals_str);

subplot(2,2,3);
pvals_hip=[log10([H7.hippo.SWR_pvalue(H7.hippo.SWR_psign<0)' H5.hippo.SWR_pvalue(H5.hippo.SWR_psign<0)' H1.hippo.SWR_pvalue(H1.hippo.SWR_psign<0)']) -log10([H7.hippo.SWR_pvalue(H7.hippo.SWR_psign>=0)' H5.hippo.SWR_pvalue(H5.hippo.SWR_psign>=0)' H1.hippo.SWR_pvalue(H1.hippo.SWR_psign>=0)'])];
pd=fitdist(pvals_hip','Normal')
pvals_hip(find(pvals_hip<-5))=-5;
pvals_hip(find(pvals_hip>5))=5;
histogram(pvals_hip(find(pvals_hip>-2 & pvals_hip<2)),-10:.5:10); hold on;
histogram(pvals_hip(find(pvals_hip<=-2)),-10:.5:10); 
histogram(pvals_hip(find(pvals_hip>=2)),-10:.5:10); 
ci = paramci(pd,'Alpha',.01)
h1_hip = lillietest(pvals_hip);
h2_hip = jbtest(pvals_hip);
h3_hip = adtest(pvals_hip);