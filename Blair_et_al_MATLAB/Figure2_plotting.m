% must run 'Figure2_analysis' to executing this script, so
% that relevant variables will exist in the workspace


figoff=0;

clim=[-1 1];

load 'Behav_effects';

behavshkgroup=[Effects.num_short.shock([1 2 3 4 10 5 6 8 13 ],6); Effects.num_short.shock2([9 11],6)]./[Effects.num_short.shock([1 2 3 4 10 5 6 8 13 ],3); Effects.num_short.shock2([9 11],3)];
behavscopgroup=Effects.num_short.scoposhk([9 11 12 14 5:8 10],6)./Effects.num_short.scoposhk([9 11 12 14 5:8 10],3);
behavbargroup=Effects.num_short.barrier([5 6 8 13 9 12],6)./Effects.num_short.barrier([5 6 8 13 9 12 ],3);

grp.shockbeh=behavshkgroup;
grp.scpshkbeh=behavscopgroup;
grp.barbeh=behavbargroup;

grpshk1=[1 2 3 4 6];
grpbar=7:12;
grpshk2=13:16;
grpshk2_2=26:27;
grpscpshk=17:25;
%grpscp=21:22;
grpshock=[grpshk1 grpshk2 26 27];
[]
[1 2 3 4 5 6 7 8 9]
grpall=[grpshk1 grpbar grpshk2 grpscpshk ];

recur=analysis_field_Nrecur./(analysis_field_Nrecur+analysis_field_Nnonrecur)

grp.shock=grpshock;%([1:7 9:10]);
grp.shock1=[1:4 6 13:16];
grp.shock2=[26:27];
grp.scpshk=grpscpshk;%([1 4:9]);
grp.bar=grpbar;

shkmat_slices=1:22;%[1:14 17:20 ];
scpmat_slices=1:18;%[1:2 7:18];
barmat_slices=1:12;

shkmat_slices=[1:14 17:20 ];
scpmat_slices=[1:2 7:18];
barmat_slices=1:12;
scopmat_slices=1:6;

df_retain_rats=grp.shock([1:7 9:10]); %
scop_forget_rats=grp.scpshk([1 4:end]);

%---------------------------------------------Fig 2C

figure(20); clf;

%------left panel: number of beelines

shk_pre_Nbee=analysis_Nbee(df_retain_rats,[1 ])' + analysis_Nbee(df_retain_rats,[4 ])';
shk_train_Nbee=sum(analysis_Nbee(df_retain_rats,[2 3])') + sum(analysis_Nbee(df_retain_rats,[5 6])');

scp_pre_Nbee=analysis_Nbee(scop_forget_rats,[1 ])' + analysis_Nbee(scop_forget_rats,[4 ])';
scp_train_Nbee=sum(analysis_Nbee(scop_forget_rats,[2 3])') + sum(analysis_Nbee(scop_forget_rats,[5 6])');

subplot(9,6,[31 32 37 38]-24); hold off;
%mean bars
bar(1:2,nanmean([shk_pre_Nbee; shk_train_Nbee]')); hold on;
bar(3:4,nanmean([scp_pre_Nbee; scp_train_Nbee]')); hold on;
%symbol+line
plot([shk_pre_Nbee*0+1; shk_train_Nbee*0+2],[shk_pre_Nbee; shk_train_Nbee],'o-c');
plot([scp_pre_Nbee*0+3; scp_train_Nbee*0+4],[scp_pre_Nbee; scp_train_Nbee],'o-m');
%recolor lines for females
plot([shk_pre_Nbee([1 4 8])*0+1; shk_train_Nbee([1 4 8])*0+2],[shk_pre_Nbee([1 4 8]); shk_train_Nbee([1 4 8])],'-b');
plot([scp_pre_Nbee([2 5 6])*0+3; scp_train_Nbee([2 5 6])*0+4],[scp_pre_Nbee([2 5 6]); scp_train_Nbee([2 5 6])],'-r');
set(gca,'XLim',[0 5],'YLim',[0 100],'XTick',[]);
ylabel('Nbee');    

    between_factors=[ones(1,length(shk_pre_Nbee)) 2*ones(1,length(scp_pre_Nbee)) ]';
    
    datamat=[[shk_pre_Nbee'; scp_pre_Nbee' ] ...
             [shk_train_Nbee'; scp_train_Nbee' ]];
             
    [Nbee_tbl,rm] = simple_mixed_anova(datamat, between_factors, {'session'}, {'dfscop'});
        
    [h,p_df,ci,stats]=ttest(shk_pre_Nbee,shk_train_Nbee)
    [h,p_sc,ci,stats]=ttest(scp_pre_Nbee,scp_train_Nbee)
    
subplot(9,6,[31 44]-24+12); hold off;
text(1,1,['2x2 (session): ' num2str(Nbee_tbl{4,5})]);
text(1,-2,['2x2 (drug): ' num2str(Nbee_tbl{2,5})]);
text(1,-5,['2x2 (sess x drug): ' num2str(Nbee_tbl{5,5})]);
text(1,-8,['DF pre v trn: ' num2str(p_df)]);
text(1,-11,['SC pre v trn: ' num2str(p_sc)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);

%---------right panel: mean running speed

shk_pre_spd=analysis_results(df_retain_rats,[29 ])';
shk_train_spd=mean(analysis_results(df_retain_rats,[30 31])');
scp_pre_spd=analysis_results(scop_forget_rats,[29 ])';
scp_train_spd=mean(analysis_results(scop_forget_rats,[30 31])');

subplot(9,6,[31 38]-20); hold off;
bar(1:2,nanmean([shk_pre_spd; shk_train_spd]')); hold on;
bar(3:4,nanmean([scp_pre_spd; scp_train_spd]')); hold on;
plot([shk_pre_spd*0+1; shk_train_spd*0+2],[shk_pre_spd; shk_train_spd],'o-c');
plot([scp_pre_spd*0+3; scp_train_spd*0+4],[scp_pre_spd; scp_train_spd],'o-m');
plot([shk_pre_spd([1 4 8])*0+1; shk_train_spd([1 4 8])*0+2],[shk_pre_spd([1 4 8]); shk_train_spd([1 4 8])],'-b');
plot([scp_pre_spd([2 5 6])*0+3; scp_train_spd([2 5 6])*0+4],[scp_pre_spd([2 5 6]); scp_train_spd([2 5 6])],'-r');
set(gca,'XLim',[0 5],'YLim',[0 80],'XTick',[]);
ylabel('spd');    

    datamat=[[shk_pre_spd'; scp_pre_spd' ] ...
             [shk_train_spd'; scp_train_spd' ]];
             
    [spd_tbl,rm] = simple_mixed_anova(datamat, between_factors, {'session'}, {'dfscop'});

    [h,p_df,ci,stats]=ttest(shk_pre_Nbee,shk_train_spd)
    [h,p_sc,ci,stats]=ttest(scp_pre_Nbee,scp_train_spd)

subplot(9,6,[31 44]-20+12); hold off;
text(1,1,['2x2 (session): ' num2str(spd_tbl{4,5})]);
text(1,-2,['2x2 (drug): ' num2str(spd_tbl{2,5})]);
text(1,-5,['2x2 (sess x drug): ' num2str(spd_tbl{5,5})]);
text(1,-8,['DF pre v trn: ' num2str(p_df)]);
text(1,-11,['SC pre v trn: ' num2str(p_sc)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[],'XTick',[]);

set(gcf,'Name','Figure 2C');

%---------------------------------------------Fig 2D

figure(21); clf;

subplot(9,6,[31 32 37 38]-24); hold off;

%mean bars
bar(1:2,nanmean([analysis_Nplacecells(df_retain_rats,1)'; analysis_Nplacecells(df_retain_rats,2)']')); hold on;
bar(3:4,nanmean([analysis_Nplacecells(scop_forget_rats,1)'; analysis_Nplacecells(scop_forget_rats,2)']')); hold on;
%line and symbol
plot([df_retain_rats*0+1; df_retain_rats*0+2],[analysis_Nplacecells(df_retain_rats,1)'; analysis_Nplacecells(df_retain_rats,2)'],'o-c');
plot([scop_forget_rats*0+3; scop_forget_rats*0+4],[analysis_Nplacecells(scop_forget_rats,1)'; analysis_Nplacecells(scop_forget_rats,2)'],'o-m');
%females
plot([grp.shock([1 4 9])*0+1; grp.shock([1 4 9])*0+2],[analysis_Nplacecells(grp.shock([1 4 9]),1)'; analysis_Nplacecells(grp.shock([1 4 9]),2)'],'-b');
plot([grp.scpshk([4 7 8])*0+3; grp.scpshk([4 7 8])*0+4],[analysis_Nplacecells(grp.scpshk([4 7 8]),1)'; analysis_Nplacecells(grp.scpshk([4 7 8]),2)'],'-r');
set(gca,'XLim',[0 5],'YLim',[0 400]);

    datamat=[[analysis_Nplacecells(df_retain_rats,1); analysis_Nplacecells(scop_forget_rats,1) ] ...
             [analysis_Nplacecells(df_retain_rats,2); analysis_Nplacecells(scop_forget_rats,2) ]];
             
    [nplace_tbl,rm] = simple_mixed_anova(datamat, between_factors, {'session'}, {'dfscop'}); 

    [h,p_df,ci,stats]=ttest(analysis_Nplacecells(df_retain_rats,1),analysis_Nplacecells(df_retain_rats,2))
    [h,p_sc,ci,stats]=ttest(analysis_Nplacecells(scop_forget_rats,1),analysis_Nplacecells(scop_forget_rats,2))
    
subplot(9,6,[31 44]-24+12); hold off;
text(1,1,['2x2 (session): ' num2str(nplace_tbl{4,5})]);
text(1,-2,['2x2 (drug): ' num2str(nplace_tbl{2,5})]);
text(1,-5,['2x2 (sess x drug): ' num2str(nplace_tbl{5,5})]);
text(1,-8,['DF pre v trn: ' num2str(p_df)]);
text(1,-11,['SC pre v trn: ' num2str(p_sc)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);

set(gcf,'Name','Figure 2D');

subplot(9,6,[31 32 37 38]-20); hold off;

%mean bars
bar(1:2,nanmean([analysis_meanplaceperc(df_retain_rats,1)'; analysis_meanplaceperc(df_retain_rats,2)']')); hold on;
bar(3:4,nanmean([analysis_meanplaceperc(scop_forget_rats,1)'; analysis_meanplaceperc(scop_forget_rats,2)']')); hold on;
%line and symbol
plot([df_retain_rats*0+1; df_retain_rats*0+2],[analysis_meanplaceperc(df_retain_rats,1)'; analysis_meanplaceperc(df_retain_rats,2)'],'o-c');
plot([scop_forget_rats*0+3; scop_forget_rats*0+4],[analysis_meanplaceperc(scop_forget_rats,1)'; analysis_meanplaceperc(scop_forget_rats,2)'],'o-m');
%females
plot([grp.shock([1 4 9])*0+1; grp.shock([1 4 9])*0+2],[analysis_meanplaceperc(grp.shock([1 4 9]),1)'; analysis_meanplaceperc(grp.shock([1 4 9]),2)'],'-b');
plot([grp.scpshk([4 7 8])*0+3; grp.scpshk([4 7 8])*0+4],[analysis_meanplaceperc(grp.scpshk([4 7 8]),1)'; analysis_meanplaceperc(grp.scpshk([4 7 8]),2)'],'-r');
set(gca,'XLim',[0 5],'YLim',[0 1]);

    datamat=[[analysis_meanplaceperc(df_retain_rats,1); analysis_meanplaceperc(scop_forget_rats,1) ] ...
             [analysis_meanplaceperc(df_retain_rats,2); analysis_meanplaceperc(scop_forget_rats,2) ]];
             
    [placeperc_tbl,rm] = simple_mixed_anova(datamat, between_factors, {'session'}, {'dfscop'}); 

    [h,p_df,ci,stats]=ttest(analysis_meanplaceperc(df_retain_rats,1),analysis_meanplaceperc(df_retain_rats,2))
    [h,p_sc,ci,stats]=ttest(analysis_meanplaceperc(scop_forget_rats,1),analysis_meanplaceperc(scop_forget_rats,2))
    
subplot(9,6,[31 44]-20+12); hold off;
text(1,1,['2x2 (session): ' num2str(placeperc_tbl{4,5})]);
text(1,-2,['2x2 (drug): ' num2str(placeperc_tbl{2,5})]);
text(1,-5,['2x2 (sess x drug): ' num2str(placeperc_tbl{5,5})]);
text(1,-8,['DF pre v trn: ' num2str(p_df)]);
text(1,-11,['SC pre v trn: ' num2str(p_sc)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);

%---------------------------------------------Fig 2E


figure(22); clf;

subplot(2,2,1); pie(sum(analysis_onlyboth(df_retain_rats,[2 1 3])));
title(sum(sum(analysis_onlyboth(df_retain_rats,[2 1 3]))));
subplot(2,2,3); pie(sum(analysis_onlyboth(df_retain_rats,[5 4 6])));
title(sum(sum(analysis_onlyboth(df_retain_rats,[5 4 6]))));
subplot(2,2,2); pie(sum(analysis_onlyboth(scop_forget_rats,[2 1 3])));
title(sum(sum(analysis_onlyboth(scop_forget_rats,[2 1 3]))));
subplot(2,2,4); pie(sum(analysis_onlyboth(scop_forget_rats,[5 4 6])));
title(sum(sum(analysis_onlyboth(scop_forget_rats,[5 4 6]))));


set(gcf,'Name','Figure 2E');

%---------------------------------------------Fig 2F

figure(23); clf; 

subplot(9,6,[31 32 37 38]-30); hold off;
%mean bars
bar(1:2,nanmean(1000*[analysis_results(df_retain_rats,15)'; analysis_results(df_retain_rats,16)']')); hold on;
bar(3:4,nanmean(1000*[analysis_results(scop_forget_rats,15)'; analysis_results(scop_forget_rats,16)']')); hold on;
%lines & symbols
plot([[1:7 9:10]*0+1; [1:7 9:10]*0+2],1000*[analysis_results(df_retain_rats,15)'; analysis_results(df_retain_rats,16)'],'o-c');
plot([[1 4:9]*0+3; [1 4:9]*0+4],1000*[analysis_results(scop_forget_rats,15)'; analysis_results(scop_forget_rats,16)'],'o-m');
%females
plot([[1 4 9]*0+1; [1 4 9]*0+2],1000*[analysis_results(grp.shock([1 4 9]),15)'; analysis_results(grp.shock([1 4 9]),16)'],'-b');
plot([[4 7 8]*0+3; [4 7 8]*0+4],1000*[analysis_results(grp.scpshk([4 7 8]),15)'; analysis_results(grp.scpshk([4 7 8]),16)'],'-r');
set(gca,'XLim',[0 5],'YLim',[0 6]);
ylabel('peak');
    
    datamat=[[analysis_results(df_retain_rats,15); analysis_results(scop_forget_rats,15) ] ...
             [analysis_results(df_retain_rats,16); analysis_results(scop_forget_rats,16) ]];
             
    [peakrate_tbl,rm] = simple_mixed_anova(datamat, between_factors, {'session'}, {'dfscop'});
    
    [h,p_df,ci,stats]=ttest(analysis_results(df_retain_rats,15),analysis_results(df_retain_rats,16))
    [h,p_sc,ci,stats]=ttest(analysis_results(scop_forget_rats,15),analysis_results(scop_forget_rats,16))

subplot(9,6,[31 38]-24+6); hold off;
text(1,1,['2x2 (session): ' num2str(peakrate_tbl{4,5})]);
text(1,-2,['2x2 (drug): ' num2str(peakrate_tbl{2,5})]);
text(1,-5,['2x2 (sess x drug): ' num2str(peakrate_tbl{5,5})]);
text(1,-8,['DF pre v trn: ' num2str(p_df)]);
text(1,-11,['SC pre v trn: ' num2str(p_sc)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);
    
subplot(9,6,[31 38]-26); hold off;
%mean bars
bar(1:2,nanmean([analysis_results(df_retain_rats,12)'; analysis_results(df_retain_rats,13)']')); hold on;
bar(3:4,nanmean([analysis_results(scop_forget_rats,12)'; analysis_results(scop_forget_rats,13)']')); hold on;
%lines & symbols
plot([[1:7 9:10]*0+1; [1:7 9:10]*0+2],[analysis_results(df_retain_rats,12)'; analysis_results(df_retain_rats,13)'],'o-c');
plot([[1 4:9]*0+3; [1 4:9]*0+4],[analysis_results(scop_forget_rats,12)'; analysis_results(scop_forget_rats,13)'],'o-m');
%females
plot([[1 4 9]*0+1; [1 4 9]*0+2],[analysis_results(grp.shock([1 4 9]),12)'; analysis_results(grp.shock([1 4 9]),13)'],'-b');
plot([[4 7 8]*0+3; [4 7 8]*0+4],[analysis_results(grp.scpshk([4 7 8]),12)'; analysis_results(grp.scpshk([4 7 8]),13)'],'-r');
set(gca,'XLim',[0 5],'YLim',[0 2]);
ylabel('blrate');

    datamat=[[analysis_results(df_retain_rats,12); analysis_results(scop_forget_rats,12) ] ...
             [analysis_results(df_retain_rats,13); analysis_results(scop_forget_rats,13) ]];
             
    [bl_tbl,rm] = simple_mixed_anova(datamat, between_factors, {'session'}, {'dfscop'});

    
    [h,p_df,ci,stats]=ttest(analysis_results(df_retain_rats,12),analysis_results(df_retain_rats,13))
    [h,p_sc,ci,stats]=ttest(analysis_results(scop_forget_rats,12),analysis_results(scop_forget_rats,13))

subplot(9,6,[31 38]-20+6); hold off;
text(1,1,['2x2 (session): ' num2str(bl_tbl{4,5})]);
text(1,-2,['2x2 (drug): ' num2str(bl_tbl{2,5})]);
text(1,-5,['2x2 (sess x drug): ' num2str(bl_tbl{5,5})]);
text(1,-8,['DF pre v trn: ' num2str(p_df)]);
text(1,-11,['SC pre v trn: ' num2str(p_sc)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[],'XTick',[]);    
    
subplot(9,6,[31 38]); hold off;
%mean bars
bar(1:2,nanmean([analysis_results(df_retain_rats,26)'; analysis_results(df_retain_rats,27)']')); hold on;
bar(3:4,nanmean([analysis_results(scop_forget_rats,26)'; analysis_results(scop_forget_rats,27)']')); hold on;
%lines & symbols
plot([[1:7 9:10]*0+1; [1:7 9:10]*0+2],[analysis_results(df_retain_rats,26)'; analysis_results(df_retain_rats,27)'],'o-c');
plot([[1 4:9]*0+3; [1 4:9]*0+4],[analysis_results(scop_forget_rats,26)'; analysis_results(scop_forget_rats,27)'],'o-m');
%females
plot([[1 4 9]*0+1; [1 4 9]*0+2],[analysis_results(grp.shock([1 4 9]),26)'; analysis_results(grp.shock([1 4 9]),27)'],'-b');
plot([[4 7 8]*0+3; [4 7 8]*0+4],[analysis_results(grp.scpshk([4 7 8]),26)'; analysis_results(grp.scpshk([4 7 8]),27)'],'-r');
%set(gca,'XLim',[0 5],'YLim',[0 40]);
ylabel('width');

    datamat=[[analysis_results(df_retain_rats,26); analysis_results(scop_forget_rats,26) ] ...
             [analysis_results(df_retain_rats,27); analysis_results(scop_forget_rats,27) ]];
             
    [fs_tbl,rm] = simple_mixed_anova(datamat, between_factors, {'session'}, {'dfscop'});

    [h,p_df,ci,stats]=ttest(analysis_results(df_retain_rats,26),analysis_results(df_retain_rats,27))
    [h,p_sc,ci,stats]=ttest(analysis_results(scop_forget_rats,26),analysis_results(scop_forget_rats,27))

subplot(9,6,[31 38]+12); hold off;
text(1,1,['2x2 (session): ' num2str(fs_tbl{4,5})]);
text(1,-2,['2x2 (drug): ' num2str(fs_tbl{2,5})]);
text(1,-5,['2x2 (sess x drug): ' num2str(fs_tbl{5,5})]);
text(1,-8,['DF pre v trn: ' num2str(p_df)]);
text(1,-11,['SC pre v trn: ' num2str(p_sc)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[],'XTick',[]);
    
    
subplot(9,6,[31 38]+4); hold off;
%mean bars
bar(1:2,nanmean([analysis_results(df_retain_rats,8)'; analysis_results(df_retain_rats,9)']')); hold on;
bar(3:4,nanmean([analysis_results(scop_forget_rats,8)'; analysis_results(scop_forget_rats,9)']')); hold on;
%lines & symbols
plot([[1:7 9:10]*0+1; [1:7 9:10]*0+2],[analysis_results(df_retain_rats,8)'; analysis_results(df_retain_rats,9)'],'o-c');
plot([[1 4:9]*0+3; [1 4:9]*0+4],[analysis_results(scop_forget_rats,8)'; analysis_results(scop_forget_rats,9)'],'o-m');
%females
plot([[1 4 9]*0+1; [1 4 9]*0+2],[analysis_results(grp.shock([1 4 9]),8)'; analysis_results(grp.shock([1 4 9]),9)'],'-b');
plot([[4 7 8]*0+3; [4 7 8]*0+4],[analysis_results(grp.scpshk([4 7 8]),8)'; analysis_results(grp.scpshk([4 7 8]),9)'],'-r');
set(gca,'XLim',[0 5],'YLim',[0 .6]);
ylabel('info');

    datamat=[[analysis_results(df_retain_rats,8); analysis_results(scop_forget_rats,8) ] ...
             [analysis_results(df_retain_rats,9); analysis_results(scop_forget_rats,9) ]];
             
    [bps_tbl,rm] = simple_mixed_anova(datamat, between_factors, {'session'}, {'dfscop'});

    [h,p_df,ci,stats]=ttest(analysis_results(df_retain_rats,8),analysis_results(df_retain_rats,9))
    [h,p_sc,ci,stats]=ttest(analysis_results(scop_forget_rats,8),analysis_results(scop_forget_rats,9))
    

subplot(9,6,[31 38]+16); hold off;
text(1,1,['2x2 (session): ' num2str(bps_tbl{4,5})]);
text(1,-2,['2x2 (drug): ' num2str(bps_tbl{2,5})]);
text(1,-5,['2x2 (sess x drug): ' num2str(bps_tbl{5,5})]);
text(1,-8,['DF pre v trn: ' num2str(p_df)]);
text(1,-11,['SC pre v trn: ' num2str(p_sc)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[],'XTick',[]);

set(gcf,'Name','Figure 2F');

%---------------------------------------------Fig 2G

figure(24); clf; 

subplot(9,6,[31 32 37 38]-30); hold off;

recur=analysis_field_Nrecur./(analysis_field_Nrecur+analysis_field_Nnonrecur);

bar(1:2,[mean(recur(df_retain_rats)) nanmean(recur(scop_forget_rats))]); hold on;
scatter(1+df_retain_rats*0,recur(df_retain_rats),'oc');
scatter(2+scop_forget_rats*0,recur(scop_forget_rats),'oc');
set(gca,'XLim',[0 3],'YLim',[0 1]);

[h,p,ci,stats]=ttest2(recur(df_retain_rats),recur(scop_forget_rats(~isnan(recur(scop_forget_rats)))))

subplot(9,6,[31 38]-24+6); hold off;
text(1,1,['DF v SCP: ' num2str(p)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);

set(gcf,'Name','Figure 2G');

%---------------------------------------------Fig 2H

figure(25); clf; 

hm_shkpost_shk_placeLR=[hm_shockpost_shk(find(hm_pflags_post_shk(:,1) & ismember(hm_post_shk_group,df_retain_rats)),1:23)];
[temp peakLR]=max(hm_shkpost_shk_placeLR');
hm_shkpost_shk_sortLR=sortrows([peakLR' hm_shkpost_shk_placeLR],1);
hm_shkpost_shk_sortLR=hm_shkpost_shk_sortLR(:,2:end);

hm_shkpost_shk_placeRL=[hm_shockpost_shk(find(hm_pflags_post_shk(:,2) & ismember(hm_post_shk_group,df_retain_rats)),24:46)];
[temp peakRL]=max(hm_shkpost_shk_placeRL');
hm_shkpost_shk_sortRL=sortrows([peakRL' hm_shkpost_shk_placeRL],1);
hm_shkpost_shk_sortRL=hm_shkpost_shk_sortRL(:,2:end);

subplot(6,6,[1 7 13]+2-1); imagesc(1000*conv2(1,1,hm_shkpost_shk_sortLR,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); %ylabel(size(hm_shkpost_shk_sortLR,1));
title('DF train1');
subplot(6,6,[1 7 13]+20-1); imagesc(1000*conv2(1,1,hm_shkpost_shk_sortRL,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); %ylabel(size(hm_shkpost_shk_sortRL,1));

hm_post_shk_placeLR=[hm_post_shk(find(hm_pflags_post_shk(:,1) & ismember(hm_post_shk_group,df_retain_rats)),1:23)];
%[temp peakLR]=max(hm_post_shk_placeLR');
hm_post_shk_sortLR=sortrows([peakLR' hm_post_shk_placeLR],1);
hm_post_shk_sortLR=hm_post_shk_sortLR(:,2:end);

hm_post_shk_placeRL=[hm_post_shk(find(hm_pflags_post_shk(:,2) & ismember(hm_post_shk_group,df_retain_rats)),24:46)];
%[temp peakRL]=max(hm_post_shk_placeRL');
hm_post_shk_sortRL=sortrows([peakRL' hm_post_shk_placeRL],1);
hm_post_shk_sortRL=hm_post_shk_sortRL(:,2:end);

subplot(6,6,[1 7 13]+3-1); imagesc(1000*conv2(1,1,hm_post_shk_sortLR,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); %ylabel(size(hm_post_shk_sortLR,1));
title('DF train2');
subplot(6,6,[1 7 13]+21-1); imagesc(1000*conv2(1,1,hm_post_shk_sortRL,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); %ylabel(size(hm_post_shk_sortRL,1));

hm_pre_shk_placeLR=[hm_pre_shk(find(hm_pflags_post_shk(:,1) & ismember(hm_pre_shk_group,df_retain_rats)),1:23)];
%[temp peakLR]=max(hm_pre_shk_placeLR');
hm_pre_shk_sortLR=sortrows([peakLR' hm_pre_shk_placeLR],1);
hm_pre_shk_sortLR=hm_pre_shk_sortLR(:,2:end);

hm_pre_shk_placeRL=[hm_pre_shk(find(hm_pflags_post_shk(:,2) & ismember(hm_pre_shk_group,df_retain_rats)),24:46)];
%[temp peakRL]=max(hm_pre_shk_placeRL');
hm_pre_shk_sortRL=sortrows([peakRL' hm_pre_shk_placeRL],1);
hm_pre_shk_sortRL=hm_pre_shk_sortRL(:,2:end);

subplot(6,6,[1 7 13]); imagesc(1000*conv2(1,1,hm_pre_shk_sortLR,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(['LR ' size(hm_pre_shk_sortLR,1)]);
title('DF pre48');
subplot(6,6,[1 7 13]+18); imagesc(1000*conv2(1,1,hm_pre_shk_sortRL,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(['RL ' size(hm_pre_shk_sortRL,1)]);

hm_shkpost_scpshk_placeLR=[hm_shockpost_scpshk(find(hm_pflags_post_scpshk(:,1) & ismember(hm_post_scpshk_group,scop_forget_rats)),1:23)];
[temp peakLR]=max(hm_shkpost_scpshk_placeLR');
hm_shkpost_scpshk_sortLR=sortrows([peakLR' hm_shkpost_scpshk_placeLR],1);
hm_shkpost_scpshk_sortLR=hm_shkpost_scpshk_sortLR(:,2:end);

hm_shkpost_scpshk_placeRL=[hm_shockpost_scpshk(find(hm_pflags_post_scpshk(:,2) & ismember(hm_post_scpshk_group,scop_forget_rats)),24:46)];
[temp peakRL]=max(hm_shkpost_scpshk_placeRL');
hm_shkpost_scpshk_sortRL=sortrows([peakRL' hm_shkpost_scpshk_placeRL],1);
hm_shkpost_scpshk_sortRL=hm_shkpost_scpshk_sortRL(:,2:end);

subplot(6,6,[1 7 13]+4); imagesc(1000*conv2(1,1,hm_shkpost_scpshk_sortLR,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); %ylabel(size(hm_shkpost_scpshk_sortLR,1));
title('SCOP train1');
subplot(6,6,[1 7 13]+22); imagesc(1000*conv2(1,1,hm_shkpost_scpshk_sortRL,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); %ylabel(size(hm_shkpost_scpshk_sortRL,1));

hm_post_scpshk_placeLR=[hm_post_scpshk(find(hm_pflags_post_scpshk(:,1) & ismember(hm_post_scpshk_group,scop_forget_rats)),1:23)];
%[temp peakLR]=max(hm_post_scpshk_placeLR');
hm_post_scpshk_sortLR=sortrows([peakLR' hm_post_scpshk_placeLR],1);
hm_post_scpshk_sortLR=hm_post_scpshk_sortLR(:,2:end);

hm_post_scpshk_placeRL=[hm_post_scpshk(find(hm_pflags_post_scpshk(:,2) & ismember(hm_post_scpshk_group,scop_forget_rats)),24:46)];
%[temp peakRL]=max(hm_post_scpshk_placeRL');
hm_post_scpshk_sortRL=sortrows([peakRL' hm_post_scpshk_placeRL],1);
hm_post_scpshk_sortRL=hm_post_scpshk_sortRL(:,2:end);

subplot(6,6,[1 7 13]+5); imagesc(1000*conv2(1,1,hm_post_scpshk_sortLR,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); %ylabel(size(hm_post_scpshk_sortLR,1));
title('SCOP train2');
subplot(6,6,[1 7 13]+23); imagesc(1000*conv2(1,1,hm_post_scpshk_sortRL,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); %ylabel(size(hm_post_scpshk_sortRL,1));

hm_pre_scpshk_placeLR=[hm_pre_scpshk(find(hm_pflags_post_scpshk(:,1) & ismember(hm_pre_scpshk_group,scop_forget_rats)),1:23)];
%[temp peakLR]=max(hm_pre_scpshk_placeLR');
hm_pre_scpshk_sortLR=sortrows([peakLR' hm_pre_scpshk_placeLR],1);
hm_pre_scpshk_sortLR=hm_pre_scpshk_sortLR(:,2:end);

hm_pre_scpshk_placeRL=[hm_pre_scpshk(find(hm_pflags_post_scpshk(:,2) & ismember(hm_pre_scpshk_group,scop_forget_rats)),24:46)];
%[temp peakRL]=max(hm_pre_scpshk_placeRL');
hm_pre_scpshk_sortRL=sortrows([peakRL' hm_pre_scpshk_placeRL],1);
hm_pre_scpshk_sortRL=hm_pre_scpshk_sortRL(:,2:end);

subplot(6,6,[1 7 13]+3); imagesc(1000*conv2(1,1,hm_pre_scpshk_sortLR,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(['LR ' size(hm_pre_scpshk_sortLR,1)]);
title('SCOP pre48');
subplot(6,6,[1 7 13]+21); imagesc(1000*conv2(1,1,hm_pre_scpshk_sortRL,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(['RL ' size(hm_pre_scpshk_sortRL,1)]);

set(gcf,'Name','Figure 2H');

figure(26); clf; 

subplot(2,2,1); pie(sum(analysis_onlyboth_recur(df_retain_rats,[2 1 3])));
subplot(2,2,2); pie(sum(analysis_onlyboth_recur(scop_forget_rats,[2 1 3])));

set(gcf,'Name','Figure 2H pies');

%---------------------------------------------Fig 2I

figure(27); clf;

btwn_df=nanmedian(hm_pre_shk_shortpath');
wthn_df=nanmedian(hm_post_shk_shortpath');
subplot(2,2,1);
imagesc(squeeze(nanmean(Rmatrix_place_shock1_pre(shkmat_slices,:,:),1))); colormap jet;
title('DF between');
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
subplot(2,2,2);
imagesc(squeeze(nanmean(Rmatrix_place_shock1_post(shkmat_slices,:,:),1))); 
title('DF within');
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;

btwn_sc=nanmedian(hm_pre_scp_shortpath');
wthn_sc=nanmedian(hm_post_scp_shortpath');
subplot(2,2,3);
imagesc(squeeze(nanmean(Rmatrix_place_scopshock_pre(scpmat_slices,:,:),1))); colormap jet;
title('SCOP between');
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
subplot(2,2,4);
imagesc(squeeze(nanmean(Rmatrix_place_scopshock_post(scpmat_slices,:,:),1))); 
title('SCOP within');
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;

set(gcf,'Name','Figure 2I matrices');

figure(28); clf;

subplot(9,6,[31 32 37 38]-30); hold off;

bar(1:2,nanmean([btwn_df([1:7 9:10]); wthn_df([1:7 9:10])]')); hold on;
bar(3:4,nanmean([btwn_sc([1 4:9]); wthn_sc([1 4:9])]')); hold on;
plot([btwn_df([1:7 9:10])*0+1; wthn_df([1:7 9:10])*0+2],[btwn_df([1:7 9:10]); wthn_df([1:7 9:10])],'o-c');
plot([btwn_sc([1 4:9])*0+3; wthn_sc([1 4:9])*0+4],[btwn_sc([1 4:9]); wthn_sc([1 4:9])],'o-m');
plot([btwn_df([1 4 9])*0+1; wthn_df([1 4 9])*0+2],[btwn_df([1 4 9]); wthn_df([1 4 9])],'o-g');
plot([btwn_sc([4 7 8])*0+3; wthn_sc([4 7 8])*0+4],[btwn_sc([4 7 8]); wthn_sc([4 7 8])],'o-g');
set(gca,'XLim',[0 5],'YLim',[-.5 1]);

    datamat=[[btwn_df([1:7 9:10])'; btwn_sc([1 4:9])' ] ...
             [wthn_df([1:7 9:10])'; wthn_sc([1 4:9])' ]];
             
    [popvec_tbl,rm] = simple_mixed_anova(datamat, between_factors, {'session'}, {'dfscop'}); 
    [h,p_btwn,ci,stats]=ttest2(btwn_df([1:7 9:10])', btwn_sc([1 4:9])')
    [h,p_wthn,ci,stats]=ttest2(wthn_df([1:7 9:10])', wthn_sc([1 4:9])')

subplot(9,6,[31 38]-24+6); hold off;
text(1,1,['2x2 (session): ' num2str(popvec_tbl{4,5})]);
text(1,-2,['2x2 (drug): ' num2str(popvec_tbl{2,5})]);
text(1,-5,['2x2 (sess x drug): ' num2str(popvec_tbl{5,5})]);
text(1,-8,['BTWN df v sc: ' num2str(p_btwn)]);
text(1,-11,['WTHN df v sc: ' num2str(p_wthn)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);

set(gcf,'Name','Figure 2I 2x2 bars');
