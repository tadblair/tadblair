%to be executed after running the Hipp_analysis_between_dfdf_shock_revRL

figure(3); clf; 

load 'Behav_effects';

behavshkgroup=[Effects.num_short.shock([1 2 3 4 10 5 6 8 13 ],6); Effects.num_short.shock2([10 13],6)]./[Effects.num_short.shock([1 2 3 4 10 5 6 8 13 ],3); Effects.num_short.shock([10 13],3)];
behavscopgroup=Effects.num_short.scoposhk([9 11 12 14 5:8 10],6)./Effects.num_short.scoposhk([9 11 12 14 5:8 10],3);
behavbargroup=Effects.num_short.barrier([5 6 8 13 9 12],6)./Effects.num_short.barrier([5 6 8 13 9 12 ],3);

grpshk1=[1 2 3 4 6];
grpbar=7:12;
grpshk2=13:16;
grpshk2_2=26:27;
grpscpshk=17:25;
grpshock=[grpshk1 grpshk2 26 27];

shkdex=[1:11];
shkdex2=[1:22];
%uncomment for to select only shk rats with blocked avoidance (n=9)
% shkdex=[1:7 9 11];
% shkdex2=[1:14 17 18 21 22];
% grpshock=grpshock(shkdex);
% behavshkgroup=behavshkgroup(shkdex);

scpdex=[1:9];
scpdex2=[1:18];
%uncomment for to select only scop rats with blocked avoidance (n=7)
% scpdex=[1 3 4:9];
% scpdex2=[1:2 7:18];
% grpscpshk=grpscpshk(scpdex);
% behavscopgroup=behavscopgroup(scpdex);

%grpscp=21:22;
[]
[1 2 3 4 5 6 7 8 9]
grpall=[grpshk1 grpbar grpshk2 grpscpshk ];

analysis_field_Nrecur./(analysis_field_Nrecur+analysis_field_Nnonrecur)

grp.shockbeh=behavshkgroup;
grp.scpshkbeh=behavscopgroup;
grp.barbeh=behavbargroup;

grp.shock=grpshock;
grp.shock1=[1:4 6 13:16];
grp.shock2=[26:27];
grp.scpshk=grpscpshk;
grp.bar=grpbar;

grp2=grp;
grp2.shock=grp.shock([1:7 9 10]);
grp2.scpshk=grp.scpshk([1 4:9]);

subplot(3,3,2);
mshk=11.4*(mean(allshockcells_df));
sdshk=11.4*(std(allshockcells_df))/sqrt(size(allshockcells_df,1));
% hold on; plot(mean(allnonshockcells_df));
% hold on; plot(mean(allshockinhcells_df));
hold on; 
mnon=11.4*(mean([allnonshockcells_df]));
sdnon=11.4*(std([allnonshockcells_df]))/sqrt(size([allnonshockcells_df],1));
shadedErrorBar((1/11.4)*[-22:22],mshk,sdshk); hold on;
shadedErrorBar((1/11.4)*[-22:22],mnon,sdnon); 
set(gca,'YLim',[0 4]);

prewin_df=mean(allshockcells_df(:,1:22)')';
postwin_df=mean(allshockcells_df(:,24:45)')';
prewin_sc=mean(allshockcells_scop(:,1:22)')';
postwin_sc=mean(allshockcells_scop(:,24:45)')';

[mean([prewin_df]) mean([postwin_df])]
[mean([prewin_sc]) mean([postwin_sc])]

    between_factors=[ones(1,length([prewin_df])) 2*ones(1,length([prewin_sc]))]';
    
    datamat=[[prewin_df(:); prewin_sc(:)] ...
             [postwin_df(:); postwin_sc(:)]];
             
    [shockcells.tbl,rm] = simple_mixed_anova(datamat, between_factors, {'prepost'}, {'condition'});

    
prewin_dfs=prewin_df;    
prewin_scs=prewin_sc;   
postwin_dfs=postwin_df;    
postwin_scs=postwin_sc; 
    
prewin_df=mean(allnonshockcells_df(:,1:22)')';
postwin_df=mean(allnonshockcells_df(:,24:45)')';
prewin_sc=mean(allnonshockcells_scop(:,1:22)')';
postwin_sc=mean(allnonshockcells_scop(:,24:45)')';

[mean([prewin_df]) mean([postwin_df])]
[mean([prewin_sc]) mean([postwin_sc])]

    between_factors=[ones(1,length([prewin_df])) 2*ones(1,length([prewin_sc]))]';
    
    datamat=[[prewin_df(:); prewin_sc(:)] ...
             [postwin_df(:); postwin_sc(:)]];
             
    [nonshockcells.tbl,rm] = simple_mixed_anova(datamat, between_factors, {'prepost'}, {'condition'});
        
    [h BLdf, ci, statsBLdf]=ttest2(prewin_dfs(:),prewin_df(:));
    [h BLsc, ci, statsBLsc]=ttest2(prewin_scs(:),prewin_sc(:));

subplot(3,3,3);
temp=histogram(allshockcells_peakdf(:),'Normalization','pdf');
shkdist=temp.Values;
% temp=histogram(allnonshockcells_peakdf,'Normalization','pdf');
% nondist=temp.Values;
% temp=histogram(allshockinhcells_peakdf,'Normalization','pdf');
temp=histogram([allnonshockcells_peakdf(:)],'Normalization','pdf');
inhdist=temp.Values;
plot(shkdist); hold on;
%plot(nondist);
plot(inhdist);
set(gca,'YLim',[0 .15]);
title([num2str(median(allshockcells_peakdf)) ' ' num2str(median(allnonshockcells_peakdf))]);

subplot(3,3,1);
histogram(analyze_percplace_shockresp(grp2.shock),0:.1:1);
title(num2str(mean(analyze_percplace_shockresp(grp2.shock))));

% subplot(2,4,4);
% pie([length(allshockcells_peakdf(:)) length([allnonshockcells_peakdf(:)])]);

subplot(3,3,5);
mshk=11.4*(mean(allshockcells_scop));
sdshk=11.4*(std(allshockcells_scop))/sqrt(size(allshockcells_scop,1));
% hold on; plot(mean(allnonshockcells_scop));
% hold on; plot(mean(allshockinhcells_scop));
hold on; 
mnon=11.4*(mean([allnonshockcells_scop]));
sdnon=11.4*(std([allnonshockcells_scop]))/sqrt(size([allnonshockcells_scop],1));
shadedErrorBar((1/11.4)*[-22:22],mshk,sdshk); hold on;
shadedErrorBar((1/11.4)*[-22:22],mnon,sdnon); 
set(gca,'YLim',[0 4]);

[h,p_shvnon_dfpre,ci,stats_shvnon_dfpre]=ttest2(prewin_dfs,prewin_df);
[h,p_shvnon_scpre,ci,stats_shvnon_scpre]=ttest2(prewin_scs,prewin_sc);
[h,p_dfvsc_pre,ci,stats_dfvsc_pre]=ttest2(prewin_dfs,prewin_scs);
[h,p_dfvsc_post,ci,stats_dfvsc_post]=ttest2(postwin_dfs,postwin_scs);

subplot(3,3,8); 
text(1,1,['dfpre p=' num2str(p_shvnon_dfpre)]);
text(1,-2,['scpre p=' num2str(p_shvnon_scpre)]);
text(1,-5,['dfvsc pre ' num2str(p_dfvsc_pre)]);
text(1,-8,['dfvsc post ' num2str(p_dfvsc_post)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);

subplot(3,3,6);
temp=histogram(allshockcells_peakscop(:),'Normalization','pdf');
shkdist=temp.Values;
% temp=histogram(allnonshockcells_peakscop,'Normalization','pdf');
% nondist=temp.Values;
% temp=histogram(allshockinhcells_peakscop,'Normalization','pdf');
temp=histogram([allnonshockcells_peakscop(:)],'Normalization','pdf');
inhdist=temp.Values;
plot(shkdist); hold on;
%plot(nondist);
plot(inhdist);
set(gca,'YLim',[0 .15]);
title([num2str(median(allshockcells_peakscop)) ' ' num2str(median(allnonshockcells_peakscop))]);


subplot(3,3,9); 
text(1,1,['df shk ' num2str(signtest(allshockcells_peakdf-12))]);
text(1,-2,['df non ' num2str(signtest(allnonshockcells_peakdf-12))]);
text(1,-5,['sc shk ' num2str(signtest(allshockcells_peakscop-12))]);
text(1,-8,['sc non ' num2str(signtest(allnonshockcells_peakscop-12))]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);
xlabel('sign tests');

subplot(3,3,4);
histogram(analyze_percplace_shockresp(grp2.scpshk),0:.1:1);
title(num2str(mean(analyze_percplace_shockresp(grp2.scpshk))));

[h,p_perc,ci,stats]=ttest2(analyze_percplace_shockresp(grp2.shock),analyze_percplace_shockresp(grp2.scpshk));

subplot(3,3,7); 
text(1,1,'df vs scop');
text(1,-2,['p=' num2str(p_perc)]);
text(1,-5,['t(' num2str(stats.df) ')=' num2str(stats.tstat)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);

% subplot(2,4,8);
% pie([length(allshockcells_peakscop(:)) length([allnonshockcells_peakscop(:)])]);

mean(analyze_percplace_shockresp(grp2.shock))
mean(analyze_percplace_shockresp(grp2.scpshk))
std(analyze_percplace_shockresp(grp2.shock))/sqrt(8)
std(analyze_percplace_shockresp(grp2.scpshk))/sqrt(6)

[[sum(allshockcells_approachdf) sum(~allshockcells_approachdf)]/(sum(allshockcells_approachdf)+sum(~allshockcells_approachdf)); [sum(allnonshockcells_approachdf) sum(~allnonshockcells_approachdf)]/(sum(allnonshockcells_approachdf)+sum(~allnonshockcells_approachdf))] 

[[sum(allshockcells_approachscop) sum(~allshockcells_approachscop)]/(sum(allshockcells_approachscop)+sum(~allshockcells_approachscop)); [sum(allnonshockcells_approachscop) sum(~allnonshockcells_approachscop)]/(sum(allnonshockcells_approachscop)+sum(~allnonshockcells_approachscop))] 
