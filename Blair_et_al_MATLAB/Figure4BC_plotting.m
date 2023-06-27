rfigoff=0;
distoff=2;

clim=[-1 1]*.75;

load 'Behav_effects';

behavshkgroup=[Effects.num_short.shock([1 2 3 4 10 5 6 8 13 ],6); Effects.num_short.shock2([10 13],6)]./[Effects.num_short.shock([1 2 3 4 10 5 6 8 13 ],3); Effects.num_short.shock([10 13],3)];
behavscopgroup=Effects.num_short.scoposhk([9 11 12 14 5:8 10],6)./Effects.num_short.scoposhk([9 11 12 14 5:8 10],3);
behavbargroup=Effects.num_short.barrier([5 6 8 13 9 12],6)./Effects.num_short.barrier([5 6 8 13 9 12 ],3);
%behavbargroup=behavbargroup(3:6);

grpshk1=[1 2 3 4 6];
grpbar=7:12;
grpshk2=13:16;
grpshk2_2=26:27;
grpscpshk=17:25;
grpshock=[grpshk1 grpshk2 26 27];

shkdex=[1:11];
shkdex2=[1:22];
%uncomment for to select only shk rats with blocked avoidance (n=9)
shkdex=[1:7 9 11];
shkdex2=[1:14 17 18 21:22];
grpshock=grpshock(shkdex);
behavshkgroup=behavshkgroup(shkdex);

scpdex=[1:9];
scpdex2=[1:18];
%uncomment for to select only scop rats with blocked avoidance (n=7)
scpdex=[1 2 3 4:9];
scpdex2=[1:18];
grpscpshk=grpscpshk(scpdex);
behavscopgroup=behavscopgroup(scpdex);

%grpscp=21:22;
[]
[1 2 3 4 5 6 7 8 9]
grpall=[grpshk1 grpbar grpshk2 grpscpshk ];

%analysis_field_Nrecur./(analysis_field_Nrecur+analysis_field_Nnonrecur)

grp.shockbeh=behavshkgroup;
grp.scpshkbeh=behavscopgroup;
grp.barbeh=behavbargroup;

grp.shock=grpshock;
grp.shock1=[1:4 6 13:16];
grp.shock2=[26:27];
grp.scpshk=grpscpshk;
grp.bar=grpbar;

figure(800);
clf;
subplot(1,3,1); histogram([3 4 4 4 4 4 4 4 5],0:6); set(gca,'YLim',[0 8]);
subplot(1,3,2); histogram([4 4 4 4 4 4 4],0:6); set(gca,'YLim',[0 8]);
subplot(1,3,3); histogram([4 4 4 4 4 3],0:6); set(gca,'YLim',[0 8]);


% figure(801); clf; 
% 
% subplot(2,2,1);
% [pre_v_cross_peakshift] = between_session_analysis_df(34, 35, 'pshkft', analysis_results, grp, [0 100]);
% 
% analysis_results(:,36:37) = analysis_field_Nrecur./(analysis_field_Nrecur+analysis_field_Nnonrecur);
% 
% subplot(2,2,2);
% [pre_v_cross_recurperc] = between_session_analysis_df(36, 37, '% recur', analysis_results, grp, [0 1]);
% 
% subplot(2,2,3);
% %[pre_v_cross_simR] = between_session_analysis_df(24, 25, 'simR', analysis_results, grp, [0 1]);
% [safe_v_unsafe_pvc] = between_session_analysis(19, 23, 'pvc', analysis_results, grp, [0 1]);
% 
% %subplot(2,2,4);
% %[pre_v_cross_simS] = between_session_analysis_df(38, 39, 'simS', analysis_results, grp, [0 5]);
% 
% figure(803); clf;
% subplot(2,3,1);
% [pre_v_cross_shortpath] = between_session_analysis_df(32, 33, 'all pvc', analysis_results, grp, [-.2 1]);
% subplot(2,3,3);
% [pre_v_cross_safe] = between_session_analysis_df(18, 19, 'safe pvc', analysis_results, grp, [-.2 1]);
% subplot(2,3,2);
% [pre_v_cross_unsafe] = between_session_analysis_df(20, 23, 'unsafe pvc', analysis_results, grp, [-.2 1]);
% 


% subplot(2,3,4);
% scatter(analysis_results(grpshock,23),behavshkgroup,'k');
% subplot(2,3,5);
% scatter(analysis_results(grpscpshk,23),behavscopgroup,'r');
% subplot(2,3,6);
% scatter(analysis_results(grpbar,23),behavbargroup,'b');

%------------------------ PLOT HEATMAPS-------------------------


figure(figoff+702); clf;

hm_shockpost_shk_place=[hm_shockpost_shk(find(hm_pflags_post_shk(:,1) & ismember(hm_post_shk_group,grpshock)),1:23); hm_shockpost_shk(find(hm_pflags_post_shk(:,2) & ismember(hm_post_shk_group,grpshock)),24:46)];
[temp peakpre]=max(hm_shockpost_shk_place');
clear approach;
for kk=1:length(peakpre)
    if peakpre(kk)<13
        approach(kk)=true;
    elseif peakpre(kk)>13
        approach(kk)=false;
    else
        if sum(hm_shockpost_shk_place(kk,1:12))<sum(hm_shockpost_shk_place(kk,14:end))
            approach(kk)=true;
        else
            approach(kk)=false;
        end
    end
end
[sum(approach) sum(~approach)]/length(approach)
[sum(approach) sum(~approach)]
hm_shockpost_shk_sort=sortrows([peakpre' hm_shockpost_shk_place],1);
hm_shockpost_shk_sort=hm_shockpost_shk_sort(:,2:end);
subplot(6,6,[1 7 13]+0); imagesc(1000*conv2(1,1,hm_shockpost_shk_sort,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_shk_sort,1));

hm_post_shk_place=[hm_post_shk(find(hm_pflags_post_shk(:,1) & ismember(hm_post_shk_group,grpshock)),1:23); hm_post_shk(find(hm_pflags_post_shk(:,2) & ismember(hm_post_shk_group,grpshock)),24:46)];
[temp peak]=max(hm_post_shk_place');
clear approach;
for kk=1:length(peak)
    if peak(kk)<13
        approach(kk)=true;
    elseif peak(kk)>13
        approach(kk)=false;
    else
        if sum(hm_shockpost_shk_place(kk,1:12))<sum(hm_shockpost_shk_place(kk,14:end))
            approach(kk)=true;
        else
            approach(kk)=false;
        end
    end
end
[sum(approach) sum(~approach)]/length(approach)
[sum(approach) sum(~approach)]
hm_post_shk_sort=sortrows([peakpre' hm_post_shk_place],1);
hm_post_shk_sort=hm_post_shk_sort(:,2:end);
subplot(6,6,[1 7 13]+1); imagesc(1000*conv2(1,1,hm_post_shk_sort,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_shk_sort,1));

% [temp, pkd]=max(hm_pre_shk_sort');
% subplot(6,6,19); histogram(pkd,'Normalization','cdf'); 


%binedges=[-.5 4.5:3:19.5 23.5];
clear postnorm prenorm;
binedges=[-.5 1.5:3:22.5 23.5];
binradius=1;
figure(1000+distoff); clf;
subplot(3,3,1); 
temp=histogram2(peak,peakpre,binedges,binedges,'DisplayStyle','tile');
peakdist=temp.Values;
clear postnonorm postnorm pstpre_*;
for i=1:(length(binedges)-1)
    postnonorm(i,:)=peakdist(i,:);
    postnorm(i,:)=peakdist(i,:)/sum(peakdist(i,:));
    pstpre_shk(i)=3*mean(postnorm(i,max([1 i-binradius]):min([(length(binedges)-1) i+binradius])));
    pstpre_shk_in(i)=sum(postnonorm(i,max([1 i-binradius]):min([(length(binedges)-1) i+binradius])));
    pstpre_shk_out(i)=sum(postnonorm(i,:))-pstpre_shk_in(i);
end
imagesc(postnorm(:,end:-1:1)); 
xlabel('post sess'); ylabel('pre sess'); 
colormap hot; colorbar;
subplot(3,3,4); 
%imagesc(peakdist((length(binedges)-1):-1:1,:)); colorbar;
temp=histogram(peakpre,binedges,'Normalization','probability');
prepeakdist=temp.Values;
temp=histogram(peak,binedges,'Normalization','probability');
postpeakdist=temp.Values;
hold off; plot(prepeakdist); hold on; plot(postpeakdist); 
axis square; set(gca,'YLim',[0 .25]);
subplot(3,3,7);
clear prenonorm prenorm prepst_*;
for i=1:(length(binedges)-1)
    prenonorm(i,:)=peakdist(:,i);
    prenorm(i,:)=peakdist(:,i)/sum(peakdist(:,i));
    prepst_shk(i)=3*mean(prenorm(i,max([1 i-binradius]):min([(length(binedges)-1) i+binradius])));
    prepst_shk_in(i)=sum(prenonorm(i,max([1 i-binradius]):min([(length(binedges)-1) i+binradius])));
    prepst_shk_out(i)=sum(prenonorm(i,:))-prepst_shk_in(i);
end
imagesc(prenorm(:,end:-1:1)'); 
xlabel('pre sess'); ylabel('post sess'); 
colormap hot; colorbar;
caxis([0 .5]);

figure(figoff+702);

[pvcR_df Rdf resultP] = pv_heatmap_cormatrix_pdR(hm_post_shk_sort,hm_shockpost_shk_sort);


[simRL result_ss] = pv_heatmap_tuningcorr(conv2(1,1,hm_post_shk_place,'same'),conv2(1,1,hm_shockpost_shk_place,'same'));
subplot(6,6,32-7); thist=histogram(result_ss,-5:.5:20,'Normalization','cdf');
subplot(6,6,32); plot(thist.BinEdges(1:end-1)+.25,thist.Values); hold on;  set(gca,'XLim',[-5 20],'YLim',[0 1]);
group_ss=[hm_post_shk_group(find(hm_pflags_post_shk(:,1))); hm_post_shk_group(find(hm_pflags_post_shk(:,2)))];
%result_ss=result_ss(ismember(group_ss,grpshock));
group_ss=group_ss(ismember(group_ss,grpshock));
subplot(6,6,[20 26]);
boxplot(result_ss,group_ss,'ori','horizontal','Symbol',''); 
analysis_results(grp.shock,39) = groupsummary(result_ss(:),group_ss(:),'median');
set(gca,'XLim',[-5 20],'YLim',[.5 11.5]); title(sum(result_ss>-log10(.05))/length(result_ss)); 
subplot(6,6,31);
scatter(analysis_results(grp.shock,39),behavshkgroup);
set(gca,'XLim',[-.5 4],'YLim',[0 2]);

% [temp, pkd]=max(hm_pre_scpshk_sort');
% subplot(6,6,19+2); histogram(pkd,'Normalization','cdf'); 

hm_shockpost_scpshk_place=[hm_shockpost_scpshk(find(hm_pflags_post_scpshk(:,1) & ismember(hm_post_scpshk_group,grpscpshk)),1:23); hm_shockpost_scpshk(find(hm_pflags_post_scpshk(:,2) & ismember(hm_post_scpshk_group,grpscpshk)),24:46)];
for i=1:length(grpscpshk)
%    [grpshock(i) sum(group_ss==grpshock(i))]
    [grpscpshk(i) size([hm_shockpost_scpshk(find(hm_pflags_post_scpshk(:,1) & ismember(hm_post_scpshk_group,grpscpshk(i))),1:23); hm_shockpost_scpshk(find(hm_pflags_post_scpshk(:,2) & ismember(hm_post_scpshk_group,grpscpshk(i))),24:46)],1)]
end
[temp peakpre]=max(hm_shockpost_scpshk_place');
clear approach;
for kk=1:length(peakpre)
    if peakpre(kk)<13
        approach(kk)=true;
    elseif peakpre(kk)>13
        approach(kk)=false;
    else
        if sum(hm_shockpost_shk_place(kk,1:12))<sum(hm_shockpost_shk_place(kk,14:end))
            approach(kk)=true;
        else
            approach(kk)=false;
        end
    end
end
[sum(approach) sum(~approach)]/length(approach)
[sum(approach) sum(~approach)]
hm_shockpost_scpshk_sort=sortrows([peakpre' hm_shockpost_scpshk_place],1);
hm_shockpost_scpshk_sort=hm_shockpost_scpshk_sort(:,2:end);
subplot(6,6,[1 7 13]+2); imagesc(1000*conv2(1,1,hm_shockpost_scpshk_sort,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_scpshk_sort,1));

hm_post_scpshk_place=[hm_post_scpshk(find(hm_pflags_post_scpshk(:,1) & ismember(hm_post_scpshk_group,grpscpshk)),1:23); hm_post_scpshk(find(hm_pflags_post_scpshk(:,2) & ismember(hm_post_scpshk_group,grpscpshk)),24:46)];
[temp peak]=max(hm_post_scpshk_place');
clear approach;
for kk=1:length(peak)
    if peakpre(kk)<13
        approach(kk)=true;
    elseif peak(kk)>13
        approach(kk)=false;
    else
        if sum(hm_shockpost_shk_place(kk,1:12))<sum(hm_shockpost_shk_place(kk,14:end))
            approach(kk)=true;
        else
            approach(kk)=false;
        end
    end
end
[sum(approach) sum(~approach)]/length(approach)
[sum(approach) sum(~approach)]
hm_post_scpshk_sort=sortrows([peakpre' hm_post_scpshk_place],1);
hm_post_scpshk_sort=hm_post_scpshk_sort(:,2:end);
subplot(6,6,[1 7 13]+3); imagesc(1000*conv2(1,1,hm_post_scpshk_sort,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_scpshk_sort,1));

figure(1000+distoff); 
subplot(3,3,2); 
temp=histogram2(peak,peakpre,binedges,binedges,'DisplayStyle','tile');
peakdist=temp.Values;
for i=1:(length(binedges)-1)
    postnonorm(i,:)=peakdist(i,:);
    postnorm(i,:)=peakdist(i,:)/sum(peakdist(i,:));
    pstpre_scp(i)=3*mean(postnorm(i,max([1 i-binradius]):min([(length(binedges)-1) i+binradius])));
    pstpre_scp_in(i)=sum(postnonorm(i,max([1 i-binradius]):min([(length(binedges)-1) i+binradius])));
    pstpre_scp_out(i)=sum(postnonorm(i,:))-pstpre_scp_in(i);
end
imagesc(postnorm(:,end:-1:1)'); 
xlabel('post sess'); ylabel('pre sess'); 
colormap hot; colorbar;
subplot(3,3,5); 
%imagesc(peakdist((length(binedges)-1):-1:1,:)); colorbar;
temp=histogram(peakpre,binedges,'Normalization','probability');
prepeakdist=temp.Values;
temp=histogram(peak,binedges,'Normalization','probability');
postpeakdist=temp.Values;
hold off; plot(prepeakdist); hold on; plot(postpeakdist); 
axis square; set(gca,'YLim',[0 .25]);
subplot(3,3,8);
for i=1:(length(binedges)-1)
    prenonorm(i,:)=peakdist(:,i);
    prenorm(i,:)=peakdist(:,i)/sum(peakdist(:,i));
    prepst_scp(i)=3*mean(prenorm(i,max([1 i-binradius]):min([(length(binedges)-1) i+binradius])));
    prepst_scp_in(i)=sum(prenonorm(i,max([1 i-binradius]):min([(length(binedges)-1) i+binradius])));
    prepst_scp_out(i)=sum(prenonorm(i,:))-prepst_scp_in(i);
end
imagesc(prenorm(:,end:-1:1)'); 
xlabel('pre sess'); ylabel('post sess'); 
colormap hot; colorbar;
caxis([0 .5]);

figure(2000); clf; temp=histogram2(peak,peakpre,[.5:23.5],[.5:23.5],'DisplayStyle','tile');
peakdist=temp.Values;
for i=1:23
    shknm(i,:)=peakdist(i,:)/sum(peakdist(i,:));
    postpre_scp(i)=3*mean(shknm(i,max([1 i-2]):min([23 i+2])));
    prenm(i,:)=peakdist(:,i)/sum(peakdist(:,i));
    prepost_scp(i)=3*mean(prenm(i,max([1 i-2]):min([23 i+2])));
end

figure(figoff+702);

[pvcR_sc Rsc resultP] = pv_heatmap_cormatrix_pdR(hm_post_scpshk_sort,hm_shockpost_scpshk_sort);

[simRL result_ss] = pv_heatmap_tuningcorr(conv2(1,1,hm_post_scpshk_place,'same'),conv2(1,1,hm_shockpost_scpshk_place,'same'));
subplot(6,6,32+2-7); thist=histogram(result_ss,-5:.5:20,'Normalization','cdf');
subplot(6,6,32); plot(thist.BinEdges(1:end-1)+.25,thist.Values);  set(gca,'XLim',[-5 20],'YLim',[0 1]);
group_ss=[hm_post_scpshk_group(find(hm_pflags_post_scpshk(:,1))); hm_post_scpshk_group(find(hm_pflags_post_scpshk(:,2)))];
subplot(6,6,[20 26]+2);
%result_ss=result_ss(ismember(group_ss,grpscpshk));
group_ss=group_ss(ismember(group_ss,grpscpshk));
boxplot(result_ss,group_ss,'ori','horizontal','Symbol',''); 
analysis_results(grp.scpshk,39) = groupsummary(result_ss(:),group_ss(:),'median');
set(gca,'XLim',[-5 20],'YLim',[.5 11.5]); title(sum(result_ss>-log10(.05))/length(result_ss)); 
subplot(6,6,33);
scatter(analysis_results(grp.scpshk,39),behavscopgroup);
temp=polyfit(analysis_results(grp.scpshk,39),behavscopgroup,1);
hold on; plot([-.5 4],temp(2)+[[-.5 4]]*temp(1),'r');
set(gca,'XLim',[-.5 4],'YLim',[0 2]);

hm_post_bar_place=[hm_post_bar(find(hm_pflags_post_bar(:,1)),1:23); hm_post_bar(find(hm_pflags_post_bar(:,2)),24:46)];
[temp peak]=max(hm_post_bar_place');
hm_post_bar_sort=sortrows([peak' hm_post_bar_place],1);
hm_post_bar_sort=hm_post_bar_sort(:,2:end);
subplot(6,6,[1 7 13]+4); imagesc(1000*conv2(1,1,hm_post_bar_sort,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_bar_sort,1));

hm_shockpost_bar_place=[hm_shockpost_bar(find(hm_pflags_post_bar(:,1)),1:23); hm_shockpost_bar(find(hm_pflags_post_bar(:,2)),24:46)];
[temp peakpre]=max(hm_shockpost_bar_place');
hm_shockpost_bar_sort=sortrows([peak' hm_shockpost_bar_place],1);
hm_shockpost_bar_sort=hm_shockpost_bar_sort(:,2:end);
subplot(6,6,[1 7 13]+5); imagesc(1000*conv2(1,1,hm_shockpost_bar_sort,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_bar_sort,1));

figure(1000+distoff); 
subplot(3,3,3); 
%temp=histogram2(peak,peakpre,[-.5 1.5 4.5 7.5 10.5 12.5 15.5 19.5 22.5 23.5],[-.5 1.5 4.5 7.5 10.5 12.5 15.5 19.5 22.5 23.5],'DisplayStyle','tile');
temp=histogram2(peak,peakpre,binedges,binedges,'DisplayStyle','tile');
peakdist=temp.Values;
for i=1:(length(binedges)-1)
    postnorm(i,:)=peakdist(i,:)/sum(peakdist(i,:));
    pstpre_bar(i)=3*mean(postnorm(i,max([1 i-binradius]):min([(length(binedges)-1) i+binradius])));
end
imagesc(postnorm((length(binedges)-1):-1:1,:)'); 
xlabel('post sess'); ylabel('pre sess'); 
colormap hot; colorbar;
subplot(3,3,6); 
%imagesc(peakdist((length(binedges)-1):-1:1,:)); colorbar;
temp=histogram(peakpre,binedges,'Normalization','probability');
prepeakdist=temp.Values;
temp=histogram(peak,binedges,'Normalization','probability');
postpeakdist=temp.Values;
hold off; plot(prepeakdist); hold on; plot(postpeakdist); 
axis square; set(gca,'YLim',[0 .25]);
subplot(3,3,9);
for i=1:(length(binedges)-1)
    prenorm(i,:)=peakdist(:,i)/sum(peakdist(:,i));
    prepst_bar(i)=3*mean(prenorm(i,max([1 i-binradius]):min([(length(binedges)-1) i+binradius])));
end
imagesc(prenorm((length(binedges)-1):-1:1,:)'); 
xlabel('pre sess'); ylabel('post sess'); 
colormap hot; colorbar;
caxis([0 .5]);

figure(2000); clf; temp=histogram2(peak,peakpre,[.5:23.5],[.5:23.5],'DisplayStyle','tile');
peakdist=temp.Values;
for i=1:23
    shknm(i,:)=peakdist(i,:)/sum(peakdist(i,:));
    postpre_bar(i)=3*mean(shknm(i,max([1 i-2]):min([23 i+2])));
    prenm(i,:)=peakdist(:,i)/sum(peakdist(:,i));
    prepost_bar(i)=3*mean(prenm(i,max([1 i-2]):min([23 i+2])));
end

figure(figoff+702);

[simRL result_ss] = pv_heatmap_tuningcorr(conv2(1,1,hm_post_bar_place,'same'),conv2(1,1,hm_shockpost_bar_place,'same'));
subplot(6,6,32+4-7); thist=histogram(result_ss,-5:.5:20,'Normalization','cdf');
subplot(6,6,32); plot(thist.BinEdges(1:end-1)+.25,thist.Values); set(gca,'XLim',[-5 20],'YLim',[0 1]);
group_ss=[hm_post_bar_group(find(hm_pflags_post_bar(:,1))); hm_post_bar_group(find(hm_pflags_post_bar(:,2)))];
subplot(6,6,[20 26]+4);
boxplot(result_ss,group_ss,'ori','horizontal','Symbol',''); 
analysis_results(grp.bar,39) = groupsummary(result_ss(:),group_ss(:),'median');
set(gca,'XLim',[-5 20],'YLim',[.5 11.5]); title(sum(result_ss>-log10(.05))/length(result_ss));
subplot(6,6,35);
scatter(analysis_results(grp.bar,39),behavbargroup);
set(gca,'XLim',[-.5 4],'YLim',[0 2]);

% figure(2000); clf;
% subplot(3,2,1); 
% plot(postpre_shk); hold on; 
% plot(postpre_scp);  
% plot(postpre_bar);  
% subplot(3,2,2); set(gca,'YLim',[0 1]);
% plot(prepost_shk); hold on;
% plot(prepost_scp);
% plot(prepost_bar); set(gca,'YLim',[0 1]);

figure(2001); %clf;
subplot(3,2,1+distoff); hold off;
plot(pstpre_shk); hold on;
if distoff==4
    title('shock'); 
else
    title('non shock');
end
plot(pstpre_scp);  
%plot(pstpre_bar);  
set(gca,'YLim',[0 1]);
subplot(3,2,2+distoff); hold off;
plot(prepst_shk); hold on;
plot(prepst_scp);
%plot(prepst_bar);  
set(gca,'YLim',[0 1]);

clim=[-.5 .5];

figure(3001); %clf;

if distoff==2
    distoff=3;
end
subplot(3,2,1+distoff-3);
imagesc(Rdf); colormap jet;
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;

subplot(3,2,3+distoff-3);
imagesc(Rsc); colormap jet;
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;

subplot(3,2,5+distoff-3);

for i=2:22
    dfval=[];
    scval=[];
    for j=(i-1):(i+1)
        for k=(i-1):(i+1)
            dfval=[dfval Rdf(j,k)];
            scval=[scval Rsc(j,k)];
        end
    end
    [h p]=ttest2(dfval,scval);
    pp(i)=-log10(p)*sign(mean(dfval)-mean(scval));
end
plot(pp); hold on;
plot([0 22],[-2 -2]);
set(gca,'YLim',[-15 2]);

%----------------------- DRUG FREE --------------------------

figure(651);
clf; 

rvalx=-.975:.05:.975;

cmin=-1; cmax=1;

clear p_*;
for r=1:length(grpshock)
    rr=shkdex(r);
    try
    [p_signrank(r), h]=signrank(hm_pre_shk_shortpath(rr,:),hm_post_shk_shortpath(rr,:));
    catch
        p_signrank(r)=NaN;
    end
%    [p_ranksum(r), h]=ranksum(hm_pre_shk_shortpath(r,:),hm_post_shk_shortpath(r,:));
end
shk_shortpath_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_shk_shortpath(shkdex,:)')-nanmedian(hm_post_shk_shortpath(shkdex,:)'));
shk_shortpath_dex=-log10(p_signrank) .* sign(nanmedian(hm_pre_shk_shortpath(shkdex,:)')-nanmedian(hm_post_shk_shortpath(shkdex,:)'));
subplot(9,6,32-6*4-1); 
temp1=histogram(hm_pre_shk_shortpath(shkdex,:),-1:.05:1,'Normalization','cdf'); title(['inc ' num2str(sum(shk_shortpath_sig>0))]);
subplot(9,6,32-6*2-1); 
temp2=histogram(hm_post_shk_shortpath(shkdex,:),-1:.05:1,'Normalization','cdf'); title(['dec ' num2str(sum(shk_shortpath_sig>0))]);
subplot(9,6,32-6*3); 
%plot(rvalx,temp1.Values);
hold on; plot(rvalx,temp2.Values); 
set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+1-6*3);
% boxplot(hm_pre_shk_shortpath','ori','horizontal','Symbol',''); title('shortpath');
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
subplot(9,6,[19 25]+1);
boxplot(hm_post_shk_shortpath(shkdex,:)','ori','horizontal','Symbol',''); 
set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% subplot(9,6,[31 32 37 38]);
% imagesc(squeeze(nanmean(Rmatrix_place_shock1_pre,1))); colormap jet;
% set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
% subplot(9,6,[31 32 37 38]+12);
% imagesc(squeeze(nanmean(Rmatrix_place_shock1_post(shkdex2,:,:),1))); colormap jet;
% set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;

subplot(9,6,[31 32 37 38]+12);
imagesc(Rdf); colormap jet;
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;

clear p_*;
for r=1:length(grpshock)
    rr=shkdex(r);
    try
    [p_signrank(r), h]=signrank(hm_pre_shk_safe(rr,:),hm_post_shk_safe(rr,:));
    catch
        p_signrank(r)=NaN;
    end
%    [p_ranksum(r), h]=ranksum(hm_pre_shk_safe(r,:),hm_post_shk_safe(r,:));
end
shk_safe_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_shk_safe(shkdex,:)')-nanmedian(hm_post_shk_safe(shkdex,:)'));
shk_safe_dex=-log10(p_signrank) .* sign(nanmedian(hm_pre_shk_safe(shkdex,:)')-nanmedian(hm_post_shk_safe(shkdex,:)'));
subplot(9,6,[32 38]+1);
scatter(shk_safe_dex,behavshkgroup);
[R_shk_safe,P]=corrcoef(shk_safe_dex,behavshkgroup)
subplot(9,6,32-6*4+1); 
temp1=histogram(hm_pre_shk_safe(shkdex,:),-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*2+1); 
temp2=histogram(hm_post_shk_safe(shkdex,:),-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*3+2); 
%plot(rvalx,temp1.Values);
hold on; plot(rvalx,temp2.Values); 
set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+3-6*3);
% boxplot(hm_pre_shk_safe','ori','horizontal','Symbol',''); title('safe');
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
subplot(9,6,[19 25]+3);
boxplot(hm_post_shk_safe(shkdex,:)','ori','horizontal','Symbol','');
set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);


clear p_*;
for r=1:length(grpshock)
    rr=shkdex(r);
    try
    [p_signrank(r), h]=signrank(hm_pre_shk_unsafe(rr,:),hm_post_shk_unsafe(rr,:));
    catch
        p_signrank(r)=NaN;
    end
%    [p_ranksum(r), h]=ranksum(hm_pre_shk_unsafe(r,:),hm_post_shk_unsafe(r,:));
end
shk_unsafe_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_shk_unsafe(shkdex,:)')-nanmedian(hm_post_shk_unsafe(shkdex,:)'));
shk_unsafe_dex=-log10(p_signrank) .* sign(nanmedian(hm_pre_shk_unsafe(shkdex,:)')-nanmedian(hm_post_shk_unsafe(shkdex,:)'));
subplot(9,6,[32 38]+2+1);
scatter(shk_unsafe_dex,behavshkgroup);
[R_shk_unsafe,P]=corrcoef(shk_unsafe_dex,behavshkgroup)
subplot(9,6,32+2-6*4+1); 
temp1=histogram(hm_pre_shk_unsafe(shkdex,:),-1:.05:1,'Normalization','cdf');
subplot(9,6,32+2-6*2+1); 
temp2=histogram(hm_post_shk_unsafe(shkdex,:),-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*3+4); 
%plot(rvalx,temp1.Values);
hold on; plot(rvalx,temp2.Values); 
set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+5-6*3);
% boxplot(hm_pre_shk_unsafe','ori','horizontal','Symbol',''); title('unsafe');
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
subplot(9,6,[19 25]+5);
boxplot(hm_post_shk_unsafe(shkdex,:)','ori','horizontal','Symbol',''); 
set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);

%----------------------- SCOPOLAMINE --------------------------

% figure(602);
% clf; 

rvalx=-.975:.05:.975;

cmin=-1; cmax=1;

clear p_*;
for r=1:length(grpscpshk)
    rr=scpdex(r);
    [p_signrank(r), h]=signrank(hm_pre_scp_shortpath(rr,:),hm_post_scp_shortpath(rr,:));
%    [p_ranksum(r), h]=ranksum(hm_pre_scp_shortpath(r,:),hm_post_scp_shortpath(r,:));
end
scp_shortpath_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_scp_shortpath(scpdex,:)')-nanmedian(hm_post_scp_shortpath(scpdex,:)'));
subplot(9,6,32-6*4-1); 
temp1=histogram(hm_pre_scp_shortpath(scpdex,:),-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*2-1); 
temp2=histogram(hm_post_scp_shortpath(scpdex,:),-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*3); 
%plot(rvalx,temp1.Values);
hold on; plot(rvalx,temp2.Values); 
set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+1-6*3);
% boxplot(hm_pre_scp_shortpath','ori','horizontal','Symbol',''); title('shortpath');
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
%subplot(9,6,[19 25]+1);
subplot(9,6,[19 25]+1-6*3);
boxplot(hm_post_scp_shortpath(scpdex,:)','ori','horizontal','Symbol',''); 
set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% subplot(9,6,[31 32 37 38]+2);
% imagesc(squeeze(nanmean(Rmatrix_place_scopshock_pre,1))); colormap jet;
% set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
% subplot(9,6,[31 32 37 38]+12+2);
% imagesc(squeeze(nanmean(Rmatrix_place_scopshock_post(scpdex2,:,:),1))); colormap jet;
% set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;

subplot(9,6,[31 32 37 38]+12);
imagesc(Rsc); colormap jet;
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;

clear p_*;
for r=1:length(grpscpshk)
    rr=scpdex(r);
    [p_signrank(r), h]=signrank(hm_pre_scp_safe(rr,:),hm_post_scp_safe(rr,:));
%    [p_ranksum(r), h]=ranksum(hm_pre_scp_safe(r,:),hm_post_scp_safe(r,:));
end
scp_safe_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_scp_safe(scpdex,:)')-nanmedian(hm_post_scp_safe(scpdex,:)'));
subplot(9,6,32-6*4+1); 
temp1=histogram(hm_pre_scp_safe(scpdex,:),-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*2+1); 
temp2=histogram(hm_post_scp_safe(scpdex,:),-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*3+2); 
%plot(rvalx,temp1.Values);
hold on; plot(rvalx,temp2.Values); 
set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+3-6*3);
% boxplot(hm_pre_scp_safe','ori','horizontal','Symbol',''); title('safe');
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
subplot(9,6,[19 25]+3-6*3);
boxplot(hm_post_scp_safe(scpdex,:)','ori','horizontal','Symbol','');
set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);

clear p_*;
for r=1:length(grpscpshk)
    rr=scpdex(r);
    [p_signrank(r), h]=signrank(hm_pre_scp_unsafe(rr,:),hm_post_scp_unsafe(rr,:));
%    [p_ranksum(r), h]=ranksum(hm_pre_scp_unsafe(r,:),hm_post_scp_unsafe(r,:));
end
scp_unsafe_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_scp_unsafe(scpdex,:)')-nanmedian(hm_post_scp_unsafe(scpdex,:)'));
subplot(9,6,32+2-6*4+1); 
temp1=histogram(hm_pre_scp_unsafe(scpdex,:),-1:.05:1,'Normalization','cdf');
subplot(9,6,32+2-6*2+1); 
temp2=histogram(hm_post_scp_unsafe(scpdex,:),-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*3+4); 
%plot(rvalx,temp1.Values);
hold on; plot(rvalx,temp2.Values); 
set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+5-6*3);
% boxplot(hm_pre_scp_unsafe','ori','horizontal','Symbol',''); title('unsafe');
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
subplot(9,6,[19 25]+5-6*3);
boxplot(hm_post_scp_unsafe(scpdex,:)','ori','horizontal','Symbol',''); 
set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);

%----------------------- BARRIER --------------------------

% figure(603);
% clf; 

rvalx=-.975:.05:.975;

cmin=-1; cmax=1;

clear p_*;
for r=1:6
    [p_signrank(r), h]=signrank(hm_pre_bar_shortpath(r,:),hm_post_bar_shortpath(r,:));
%    [p_ranksum(r), h]=ranksum(hm_pre_bar_shortpath(r,:),hm_post_bar_shortpath(r,:));
end
bar_shortpath_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_bar_shortpath')-nanmedian(hm_post_bar_shortpath'));
subplot(9,6,32-6*4-1); 
temp1=histogram(hm_pre_bar_shortpath,-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*2-1); 
temp2=histogram(hm_post_bar_shortpath,-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*3); 
%plot(rvalx,temp1.Values);
hold on; plot(rvalx,temp2.Values); 
set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+1-6*3);
% boxplot(hm_pre_bar_shortpath','ori','horizontal','Symbol',''); title('shortpath');
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
%subplot(9,6,[19 25]+1);
subplot(9,6,[32 38]+0);
boxplot(hm_post_bar_shortpath','ori','horizontal','Symbol',''); 
set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% subplot(9,6,[31 32 37 38]+4);
% imagesc(squeeze(nanmean(Rmatrix_place_barrier_pre,1))); colormap jet;
% set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
subplot(9,6,[31 32 37 38]+12+4);
imagesc(squeeze(nanmean(Rmatrix_place_barrier_post,1))); colormap jet;
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;


clear p_*;
for r=1:6
    [p_signrank(r), h]=signrank(hm_pre_bar_safe(r,:),hm_post_bar_safe(r,:));
%    [p_ranksum(r), h]=ranksum(hm_pre_bar_safe(r,:),hm_post_bar_safe(r,:));
end
bar_safe_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_bar_safe')-nanmedian(hm_post_bar_safe'));
subplot(9,6,32-6*4+1); 
temp1=histogram(hm_pre_bar_safe,-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*2+1); 
temp2=histogram(hm_post_bar_safe,-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*3+2); 
%plot(rvalx,temp1.Values);
hold on; plot(rvalx,temp2.Values); 
set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+3-6*3);
% boxplot(hm_pre_bar_safe','ori','horizontal','Symbol',''); title('safe');
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% subplot(9,6,[19 25]+3);
subplot(9,6,[32 38]+2);
boxplot(hm_post_bar_safe','ori','horizontal','Symbol','');
set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);

clear p_*;
for r=1:6
    [p_signrank(r), h]=signrank(hm_pre_bar_unsafe(r,:),hm_post_bar_unsafe(r,:));
%    [p_ranksum(r), h]=ranksum(hm_pre_bar_unsafe(r,:),hm_post_bar_unsafe(r,:));
end
bar_unsafe_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_bar_unsafe')-nanmedian(hm_post_bar_unsafe'));
subplot(9,6,32+2-6*4+1); 
temp1=histogram(hm_pre_bar_unsafe,-1:.05:1,'Normalization','cdf');
subplot(9,6,32+2-6*2+1); 
temp2=histogram(hm_post_bar_unsafe,-1:.05:1,'Normalization','cdf');
subplot(9,6,32-6*3+4); 
%plot(rvalx,temp1.Values);
hold on; plot(rvalx,temp2.Values); 
set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+5-6*3);
% boxplot(hm_pre_bar_unsafe','ori','horizontal','Symbol',''); title('unsafe');
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% subplot(9,6,[19 25]+5);
subplot(9,6,[32 38]+4);
boxplot(hm_post_bar_unsafe','ori','horizontal','Symbol',''); 
set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);

clim=[-.5 .5];

figure(1001); %clf;

subplot(2,2,1+distoff-3);
imagesc(Rdf); colormap jet;
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;

subplot(2,2,3+distoff-3);
imagesc(Rsc); colormap jet;
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;

% for i=1:12
%     subplot(3,4,i);
%     imagesc(squeeze(Rmatrix_place_barrier_post(i,:,:)));
%     caxis([-1 1]); colormap jet;
% end
% 
% temp1=hm_post_shk_safe(shkdex,:);
% temp2=hm_post_scp_safe(scpdex,:);
% temp3=hm_post_bar_safe;
% temp4=hm_post_shk_unsafe(shkdex,:);
% temp5=hm_post_scp_unsafe(scpdex,:);
% temp6=hm_post_bar_unsafe;

%     between_factors=[ones(1,length(temp1(:))) 2*ones(1,length(temp2(:))) 3*ones(1,length(temp3(:)))]';
%     
%     datamat=[[temp1(:); temp2(:); temp3(:)] ...
%              [temp4(:); temp5(:); temp6(:)]];
%              
%     [safe_vs_unsafe.tbl,rm] = simple_mixed_anova(datamat, between_factors, {'safeunsafe'}, {'traintype'});
