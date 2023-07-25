figoff=0;

shkdex=[1:7 9 10];
shkdex2=[1:14 17:20];
scpdex=[1 4:9];
scpdex2=[1:2 7:18];

clim=[-1 1]*.75;

load 'Behav_effects';

% F M M F M  M  M  F  F  M  M
% 1 2 3 4 6  13 14 15 16 26 27
% 6 7 8 9 31 12 13 18 35 30 32 df shock
%                   ^        ^

% M  M  F  F  M  M  F  F  M
% 30 32 34 36 12 13 15 18 31 scop shock
%     ^  ^

% M  M  F  F  M  F  
% 12 13 18 35 30 34 barrier
%  
behavshkgroup=[Effects.num_short.shock([1 2 3 4 10 5 6 8 13 ],6); Effects.num_short.shock2([9 11],6)]./[Effects.num_short.shock([1 2 3 4 10 5 6 8 13 ],3); Effects.num_short.shock([9 11],3)];
behavscopgroup=Effects.num_short.scoposhk([9 11 12 14 5:8 10],6)./Effects.num_short.scoposhk([9 11 12 14 5:8 10],3);
behavbargroup=Effects.num_short.barrier([5 6 8 13 9 12],6)./Effects.num_short.barrier([5 6 8 13 9 12 ],3);
%behavbargroup=behavbargroup(3:6);

grpshk1=[1 2 3 4 6];
grpbar=7:12;
grpshk2=13:16;
grpshk2_2=26:27;
grpscpshk=17:25;
grpshock=[grpshk1 grpshk2 26 27];

%analysis_field_Nrecur./(analysis_field_Nrecur+analysis_field_Nnonrecur)

grp.shockbeh=behavshkgroup(shkdex);
grp.scpshkbeh=behavscopgroup([1 4:end]);
grp.barbeh=behavbargroup;

grp.shock=grpshock(shkdex);
grp.shock1=grpshock([1:7 9]);
grp.shock2=grpshock(shkdex(end));
grp.shockF_retain=grpshock([1 4 9]);
grp.scpshk=grpscpshk([1 4:end]);
grp.scpshkF_forget=grpscpshk([4 7 8]);
grp.bar=grpbar;
grp.barF=grpbar([3 4 6]);

%--------------- Figure 3A

figure(30);
clf;
subplot(1,3,1); histogram(analysis_sessint(grp.shock),0:6); set(gca,'YLim',[0 8]);
subplot(1,3,2); histogram(analysis_sessint(grp.scpshk),0:6); set(gca,'YLim',[0 8]);
subplot(1,3,3); histogram(analysis_sessint(grp.bar),0:6); set(gca,'YLim',[0 8]);

set(gcf,'Name','Figure 3A');

%--------------- Figure 3B

figure(31); clf;

% LEFT PANEL

    subplot(2,2,1); hold off; a=3; b=4;
    bar(1:2,mean([sum(analysis_beelines([grp.shock1 grp.shock2],[a b])') sum(analysis_beelines([grp.shock1 grp.shock2],[a b]+2)')])./[1  1]); hold on;
    for i=grp.shock1
        plot(1:2,[sum(analysis_beelines(i,[a b])) sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-ok');     hold on;
[sum(analysis_beelines(i,[a b])) sum(analysis_beelines(i,[a b]+2))]
    end
    for i=grp.shock2
        plot(1:2,[sum(analysis_beelines(i,[a b])) sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-sk'); 
    end
    for i=grp.shockF_retain
        plot(1:2,[sum(analysis_beelines(i,[a b])) sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-g'); 
    end
    
    bar(5:6,mean([sum(analysis_beelines([7:12],[a b])') sum(analysis_beelines([7:12],[a b]+2)')])./[1  1]); 
    hold on;
    for i=7:12
        plot(5:6,[sum(analysis_beelines(i,[a b])') sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-db'); 
    end    
    for i=[9 10 12]
        plot(5:6,[sum(analysis_beelines(i,[a b])') sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-g'); 
    end    
    
    bar(3:4,mean([sum(analysis_beelines([grp.scpshk],[a b])')])./[1  1]); 
    hold on;
    for i=grp.scpshk(1:2)
        plot(3:4,[sum(analysis_beelines(i,[a b])) sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-or'); 
    end    
    for i=grp.scpshk(3:7)
        plot(3:4,[sum(analysis_beelines(i,[a b])) sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-sr'); 
    end    %     for i=17:25
    for i=grp.scpshkF_forget
        plot(3:4,[sum(analysis_beelines(i,[a b])) sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-g'); 
    end    %     for i=17:25
    set(gca,'XLim',[0 7],'YLim',[0 100]); 
    

      between_factors=[ones(1,length(grp.shock)) 2*ones(1,length(grp.scpshk)) 3*ones(1,length(grp.bar))]';
      datamat=[sum(analysis_beelines([grp.shock1 grp.shock2],[a b])')'; sum(analysis_beelines(grp.scpshk,[a b])')'; sum(analysis_beelines(grp.bar,[a b])')'];
      [p,tbl_bee,stats] = anova1(datamat,between_factors,'off')    

    [h,p_dfsc]=ttest2(sum(analysis_beelines([grp.shock1 grp.shock2],[a b])')',sum(analysis_beelines(grp.scpshk,[a b])')');
    [h,p_dfbar]=ttest2(sum(analysis_beelines([grp.shock1 grp.shock2],[a b])')',sum(analysis_beelines(grp.bar,[a b])')');
    [h,p_barsc]=ttest2(sum(analysis_beelines(grp.bar,[a b])')',sum(analysis_beelines(grp.scpshk,[a b])')');

    subplot(2,2,3);
text(1,1,['1x3 (# bee): ' num2str(tbl_bee{2,6})]);
text(1,-2,['DF v SC: ' num2str(p_dfsc)]);
text(1,-5,['DF v BAR: ' num2str(p_dfbar)]);
text(1,-8,['BAR v SC: ' num2str(p_barsc)]);
set(gca,'XLim',[.8 5],'YLim',[-9 3],'XTick',[]);    

% RIGHT PANEL

subplot(2,2,2);
hold off; a=2; b=3;
    bar(1:2,mean(analysis_meanspeed([grp.shock1 grp.shock2],[a b]))./[1  1]); hold on;
    for i=grp.shock1
        plot(1:2,analysis_meanspeed(i,[a b])./[1 1 ],'-ok');     hold on;

    end
    for i=grp.shock2
        plot(1:2,analysis_meanspeed(i,[a b])./[1 1 ],'-sk'); 
    end
    for i=grp.shockF_retain
        plot(1:2,analysis_meanspeed(i,[a b])./[1 1 ],'-g'); 
    end
    
    bar(5:6,mean(analysis_meanspeed([7:12],[a b]))./[1  1]); 
    hold on;
    for i=7:12
        plot(5:6,analysis_meanspeed(i,[a b])./[1 1 ],'-db'); 
    end    
    for i=[9 10 12]
        plot(5:6,analysis_meanspeed(i,[a b])./[1 1 ],'-g'); 
    end    
    
    bar(3:4,mean(analysis_meanspeed([grp.scpshk],[a b]))./[1  1]); 
    hold on;
    for i=grp.scpshk(1:2)
        plot(3:4,analysis_meanspeed(i,[a b])./[1 1 ],'-or'); 
    end    
    for i=grp.scpshk(3:7)
        plot(3:4,analysis_meanspeed(i,[a b])./[1 1 ],'-sr'); 
    end    %     for i=17:25
    for i=grp.scpshkF_forget
        plot(3:4,analysis_meanspeed(i,[a b])./[1 1 ],'-g'); 
    end    %     for i=17:25
    set(gca,'XLim',[0 7],'YLim',[0 100]); 
    
    between_factors=[ones(1,length(grp.shock)) 2*ones(1,length(grp.scpshk)) 3*ones(1,length(grp.bar))]';
    
    datamat=[[analysis_meanspeed(grp.shock,a); analysis_meanspeed(grp.scpshk,a); analysis_meanspeed(grp.bar,a)] ...
             [analysis_meanspeed(grp.shock,b); analysis_meanspeed(grp.scpshk,b); analysis_meanspeed(grp.bar,b)]];
             i
    [pre48_v_post_meanspd.tbl,rm] = simple_mixed_anova(datamat, between_factors, {'sesspair'}, {'traintype'});

    [h,p_df]=ttest(analysis_meanspeed(grp.shock,a),analysis_meanspeed(grp.shock,b));
    [h,p_sc]=ttest(analysis_meanspeed(grp.scpshk,a),analysis_meanspeed(grp.scpshk,b));
    [h,p_bar]=ttest(analysis_meanspeed(grp.bar,a),analysis_meanspeed(grp.bar,b));

    subplot(2,2,4);
text(1,1,['3x3 (train): ' num2str(pre48_v_post_meanspd.tbl{2,5})]);
text(1,-2,['3x3 (sess): ' num2str(pre48_v_post_meanspd.tbl{4,5})]);
text(1,-5,['train x sess: ' num2str(pre48_v_post_meanspd.tbl{5,5})]);
text(1,-8,['DF pre v post: ' num2str(p_df)]);
text(1,-11,['SC pre v post: ' num2str(p_sc)]);
text(1,-14,['BAR pre v post: ' num2str(p_bar)]);
set(gca,'XLim',[.8 5],'YLim',[-15 3],'XTick',[]);

        set(gcf,'Name','Figure 3B');  
    
%--------------- Figure 3C
figure(39); clf; 

% subplot(2,2,1);
% [pre_v_cross_peakshift] = between_session_analysis_df(34, 35, 'pshkft', analysis_results, grp, [0 100]);

% LEFT PANEL

subplot(2,2,1);
hold off; a=2; b=3;
    bar(1:2,mean(analysis_meanplaceperc([grp.shock1 grp.shock2],[a b]))./[1  1]); hold on;
    for i=grp.shock1
        plot(1:2,analysis_meanplaceperc(i,[a b])./[1 1 ],'-ok');     hold on;
    end
    for i=grp.shock2
        plot(1:2,analysis_meanplaceperc(i,[a b])./[1 1 ],'-sk'); 
    end
    for i=grp.shockF_retain
        plot(1:2,analysis_meanplaceperc(i,[a b])./[1 1 ],'-g'); 
    end
    
    bar(5:6,mean(analysis_meanplaceperc([7:12],[a b]))./[1  1]); 
    hold on;
    for i=7:12
        plot(5:6,analysis_meanplaceperc(i,[a b])./[1 1 ],'-db'); 
    end    
    for i=[9 10 12]
        plot(5:6,analysis_meanplaceperc(i,[a b])./[1 1 ],'-g'); 
    end    
    
    bar(3:4,mean(analysis_meanplaceperc([grp.scpshk],[a b]))./[1  1]); 
    hold on;
    for i=grp.scpshk(1:2)
        plot(3:4,analysis_meanplaceperc(i,[a b])./[1 1 ],'-or'); 
    end    
    for i=grp.scpshk(3:7)
        plot(3:4,analysis_meanplaceperc(i,[a b])./[1 1 ],'-sr'); 
    end    %     for i=17:25
    for i=grp.scpshkF_forget
        plot(3:4,analysis_meanplaceperc(i,[a b])./[1 1 ],'-g'); 
    end    %     for i=17:25
    set(gca,'XLim',[0 7],'YLim',[0 1]); 
    
    between_factors=[ones(1,length(grp.shock)) 2*ones(1,length(grp.scpshk)) 3*ones(1,length(grp.bar))]';
    
    datamat=[[analysis_meanplaceperc(grp.shock,a); analysis_meanplaceperc(grp.scpshk,a); analysis_meanplaceperc(grp.bar,a)] ...
             [analysis_meanplaceperc(grp.shock,b); analysis_meanplaceperc(grp.scpshk,b); analysis_meanplaceperc(grp.bar,b)]];
             
    [placeperc_tbl,rm] = simple_mixed_anova(datamat, between_factors, {'sesspair'}, {'traintype'});

    [h,p_df]=ttest(analysis_meanplaceperc(grp.shock,a),analysis_meanplaceperc(grp.shock,b));
    [h,p_sc]=ttest(analysis_meanplaceperc(grp.scpshk,a),analysis_meanplaceperc(grp.scpshk,b));
    [h,p_bar]=ttest(analysis_meanplaceperc(grp.bar,a),analysis_meanplaceperc(grp.bar,b));

    subplot(2,2,3);
text(1,1,['3x2 (train): ' num2str(placeperc_tbl{2,5})]);
text(1,-2,['3x2 (sess): ' num2str(placeperc_tbl{4,5})]);
text(1,-5,['train x sess: ' num2str(placeperc_tbl{5,5})]);
text(1,-8,['DF pre v post: ' num2str(p_df)]);
text(1,-11,['SC pre v post: ' num2str(p_sc)]);
text(1,-14,['BAR pre v post: ' num2str(p_bar)]);
set(gca,'XLim',[.8 5],'YLim',[-15 3],'XTick',[]);



figure(32); clf; 

% subplot(2,2,1);
% [pre_v_cross_peakshift] = between_session_analysis_df(34, 35, 'pshkft', analysis_results, grp, [0 100]);

% LEFT PANEL

subplot(2,2,1);
hold off; a=2; b=3;
    bar(1:2,mean(analysis_Nplacecells([grp.shock1 grp.shock2],[a b]))./[1  1]); hold on;
    for i=grp.shock1
        plot(1:2,analysis_Nplacecells(i,[a b])./[1 1 ],'-ok');     hold on;
    end
    for i=grp.shock2
        plot(1:2,analysis_Nplacecells(i,[a b])./[1 1 ],'-sk'); 
    end
    for i=grp.shockF_retain
        plot(1:2,analysis_Nplacecells(i,[a b])./[1 1 ],'-g'); 
    end
    
    bar(5:6,mean(analysis_Nplacecells([7:12],[a b]))./[1  1]); 
    hold on;
    for i=7:12
        plot(5:6,analysis_Nplacecells(i,[a b])./[1 1 ],'-db'); 
    end    
    for i=[9 10 12]
        plot(5:6,analysis_Nplacecells(i,[a b])./[1 1 ],'-g'); 
    end    
    
    bar(3:4,mean(analysis_Nplacecells([grp.scpshk],[a b]))./[1  1]); 
    hold on;
    for i=grp.scpshk(1:2)
        plot(3:4,analysis_Nplacecells(i,[a b])./[1 1 ],'-or'); 
    end    
    for i=grp.scpshk(3:7)
        plot(3:4,analysis_Nplacecells(i,[a b])./[1 1 ],'-sr'); 
    end    %     for i=17:25
    for i=grp.scpshkF_forget
        plot(3:4,analysis_Nplacecells(i,[a b])./[1 1 ],'-g'); 
    end    %     for i=17:25
    set(gca,'XLim',[0 7],'YLim',[0 420]); 
    
    between_factors=[ones(1,length(grp.shock)) 2*ones(1,length(grp.scpshk)) 3*ones(1,length(grp.bar))]';
    
    datamat=[[analysis_Nplacecells(grp.shock,a); analysis_Nplacecells(grp.scpshk,a); analysis_Nplacecells(grp.bar,a)] ...
             [analysis_Nplacecells(grp.shock,b); analysis_Nplacecells(grp.scpshk,b); analysis_Nplacecells(grp.bar,b)]];
             
    [nplace_tbl,rm] = simple_mixed_anova(datamat, between_factors, {'sesspair'}, {'traintype'});

    [h,p_df]=ttest(analysis_Nplacecells(grp.shock,a),analysis_Nplacecells(grp.shock,b));
    [h,p_sc]=ttest(analysis_Nplacecells(grp.scpshk,a),analysis_Nplacecells(grp.scpshk,b));
    [h,p_bar]=ttest(analysis_Nplacecells(grp.bar,a),analysis_Nplacecells(grp.bar,b));

    subplot(2,2,3);
text(1,1,['3x2 (train): ' num2str(nplace_tbl{2,5})]);
text(1,-2,['3x2 (sess): ' num2str(nplace_tbl{4,5})]);
text(1,-5,['train x sess: ' num2str(nplace_tbl{5,5})]);
text(1,-8,['DF pre v post: ' num2str(p_df)]);
text(1,-11,['SC pre v post: ' num2str(p_sc)]);
text(1,-14,['BAR pre v post: ' num2str(p_bar)]);
set(gca,'XLim',[.8 5],'YLim',[-15 3],'XTick',[]);

% RIGHT PANEL

analysis_results(:,36) = analysis_field_Nrecur./(analysis_field_Nrecur+analysis_field_Nnonrecur);
analysis_results(:,37) = analysis_field_Nrecur./(analysis_field_Nrecur+analysis_field_Nnonrecur);

subplot(2,2,2);
[pre_v_cross_recurperc] = between_session_analysis_df(36, 37, '% recur', analysis_results, grp, [0 1]);
hold off;
    bar(1,mean(analysis_results([grp.shock1 grp.shock2],37))); hold on;
    for i=grp.shock1
        scatter(1,analysis_results(i,37),'ok');    
    end
    for i=grp.shock2
        scatter(1,analysis_results(i,37),'sk');    
    end
    
    bar(3,mean(analysis_results([7:12],37))); 
    hold on;
    for i=7:12
        scatter(3,analysis_results(i,37),'db');    
    end    
    
    bar(2,mean(analysis_results([grp.scpshk],37))); 
    hold on;
    for i=grp.scpshk(1:2)
        scatter(2,analysis_results(i,37),'or');    
    end    
    for i=grp.scpshk(3:7)
        scatter(2,analysis_results(i,37),'sr');    
    end    %     for i=17:25
    set(gca,'XLim',[0 4],'YLim',[0 1]); 

    subplot(2,2,4);
text(1,1,['1x3 (traintype): ' num2str(pre_v_cross_recurperc.tbl{2,6})]);
text(1,-2,['DF vs SC: ' num2str(pre_v_cross_recurperc.dfsc_unpair_cross)]);
text(1,-5,['DF vs BAR: ' num2str(pre_v_cross_recurperc.dfbar_upair_cross)]);
text(1,-8,['SC vs BAR: ' num2str(pre_v_cross_recurperc.scbar_unpair_cross)]);
%text(1,-11,['WTHN df v sc: ' num2str(p_wthn)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);    
    
    set(gcf,'Name','Figure 3C');    
    
%--------------- Figure 3D
figure(33); clf;

% HEATMAPS

hm_post_shk_placeLR=[hm_post_shk(find(hm_pflags_post_shk(:,1) & ismember(hm_post_shk_group,grp.shock)),1:23)];
[temp peakLR]=max(hm_post_shk_placeLR');
hm_post_shk_sortLR=sortrows([peakLR' hm_post_shk_placeLR],1);
hm_post_shk_sortLR=hm_post_shk_sortLR(:,2:end);

hm_post_shk_placeRL=[hm_post_shk(find(hm_pflags_post_shk(:,2) & ismember(hm_post_shk_group,grp.shock)),24:46)];
[temp peakRL]=max(hm_post_shk_placeRL');
hm_post_shk_sortRL=sortrows([peakRL' hm_post_shk_placeRL],1);
hm_post_shk_sortRL=hm_post_shk_sortRL(:,2:end);

subplot(6,6,[1 7 13]+0); imagesc(1000*conv2(1,1,hm_post_shk_sortLR,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_shk_sortLR,1));
subplot(6,6,[1 7 13]+18); imagesc(1000*conv2(1,1,hm_post_shk_sortRL,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_shk_sortRL,1));

hm_shockpost_shk_placeLR=[hm_shockpost_shk(find(hm_pflags_post_shk(:,1) & ismember(hm_post_shk_group,grp.shock)),1:23)];
hm_shockpost_shk_sortLR=sortrows([peakLR' hm_shockpost_shk_placeLR],1);
hm_shockpost_shk_sortLR=hm_shockpost_shk_sortLR(:,2:end);

hm_shockpost_shk_placeRL=[hm_shockpost_shk(find(hm_pflags_post_shk(:,2) & ismember(hm_post_shk_group,grp.shock)),24:46)];
hm_shockpost_shk_sortRL=sortrows([peakRL' hm_shockpost_shk_placeRL],1);
hm_shockpost_shk_sortRL=hm_shockpost_shk_sortRL(:,2:end);

subplot(6,6,[1 7 13]+1); imagesc(1000*conv2(1,1,hm_shockpost_shk_sortLR,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_shk_sortLR,1));
subplot(6,6,[1 7 13]+19); imagesc(1000*conv2(1,1,hm_shockpost_shk_sortRL,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_shk_sortRL,1));

[temp postpeaks_df]=max([hm_post_shk_sortLR; hm_post_shk_sortRL]');
[temp prepeaks_df]=max([hm_shockpost_shk_sortLR; hm_shockpost_shk_sortRL]');

hm_post_scpshk_placeLR=[hm_post_scpshk(find(hm_pflags_post_scpshk(:,1) & ismember(hm_post_scpshk_group,grp.scpshk)),1:23)];
[temp peakLR]=max(hm_post_scpshk_placeLR');
hm_post_scpshk_sortLR=sortrows([peakLR' hm_post_scpshk_placeLR],1);
hm_post_scpshk_sortLR=hm_post_scpshk_sortLR(:,2:end);

hm_post_scpshk_placeRL=[hm_post_scpshk(find(hm_pflags_post_scpshk(:,2) & ismember(hm_post_scpshk_group,grp.scpshk)),24:46)];
[temp peakRL]=max(hm_post_scpshk_placeRL');
hm_post_scpshk_sortRL=sortrows([peakRL' hm_post_scpshk_placeRL],1);
hm_post_scpshk_sortRL=hm_post_scpshk_sortRL(:,2:end);

subplot(6,6,[1 7 13]+2); imagesc(1000*conv2(1,1,hm_post_scpshk_sortLR,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_scpshk_sortLR,1));
subplot(6,6,[1 7 13]+20); imagesc(1000*conv2(1,1,hm_post_scpshk_sortRL,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_scpshk_sortRL,1));

hm_shockpost_scpshk_placeLR=[hm_shockpost_scpshk(find(hm_pflags_post_scpshk(:,1) & ismember(hm_post_scpshk_group,grp.scpshk)),1:23)];
hm_shockpost_scpshk_sortLR=sortrows([peakLR' hm_shockpost_scpshk_placeLR],1);
hm_shockpost_scpshk_sortLR=hm_shockpost_scpshk_sortLR(:,2:end);

hm_shockpost_scpshk_placeRL=[hm_shockpost_scpshk(find(hm_pflags_post_scpshk(:,2) & ismember(hm_post_scpshk_group,grp.scpshk)),24:46)];
hm_shockpost_scpshk_sortRL=sortrows([peakRL' hm_shockpost_scpshk_placeRL],1);
hm_shockpost_scpshk_sortRL=hm_shockpost_scpshk_sortRL(:,2:end);

subplot(6,6,[1 7 13]+3); imagesc(1000*conv2(1,1,hm_shockpost_scpshk_sortLR,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_scpshk_sortLR,1));
subplot(6,6,[1 7 13]+21); imagesc(1000*conv2(1,1,hm_shockpost_scpshk_sortRL,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_scpshk_sortRL,1));

[temp postpeaks_sc]=max([hm_post_scpshk_sortLR; hm_post_scpshk_sortRL]');
[temp prepeaks_sc]=max([hm_shockpost_scpshk_sortLR; hm_shockpost_scpshk_sortRL]');

hm_post_bar_placeLR=[hm_post_bar(find(hm_pflags_post_bar(:,1) & ismember(hm_post_bar_group,grp.bar)),1:23)];
[temp peakLR]=max(hm_post_bar_placeLR');
hm_post_bar_sortLR=sortrows([peakLR' hm_post_bar_placeLR],1);
hm_post_bar_sortLR=hm_post_bar_sortLR(:,2:end);

hm_post_bar_placeRL=[hm_post_bar(find(hm_pflags_post_bar(:,2) & ismember(hm_post_bar_group,grp.bar)),24:46)];
[temp peakRL]=max(hm_post_bar_placeRL');
hm_post_bar_sortRL=sortrows([peakRL' hm_post_bar_placeRL],1);
hm_post_bar_sortRL=hm_post_bar_sortRL(:,2:end);

subplot(6,6,[1 7 13]+4); imagesc(1000*conv2(1,1,hm_post_bar_sortLR,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_bar_sortLR,1));
subplot(6,6,[1 7 13]+22); imagesc(1000*conv2(1,1,hm_post_bar_sortRL,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_bar_sortRL,1));

hm_shockpost_bar_placeLR=[hm_shockpost_bar(find(hm_pflags_post_bar(:,1) & ismember(hm_post_bar_group,grp.bar)),1:23)];
hm_shockpost_bar_sortLR=sortrows([peakLR' hm_shockpost_bar_placeLR],1);
hm_shockpost_bar_sortLR=hm_shockpost_bar_sortLR(:,2:end);

hm_shockpost_bar_placeRL=[hm_shockpost_bar(find(hm_pflags_post_bar(:,2) & ismember(hm_post_bar_group,grp.bar)),24:46)];
hm_shockpost_bar_sortRL=sortrows([peakRL' hm_shockpost_bar_placeRL],1);
hm_shockpost_bar_sortRL=hm_shockpost_bar_sortRL(:,2:end);

subplot(6,6,[1 7 13]+5); imagesc(1000*conv2(1,1,hm_shockpost_bar_sortLR,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_bar_sortLR,1));
subplot(6,6,[1 7 13]+23); imagesc(1000*conv2(1,1,hm_shockpost_bar_sortRL,'same')); colormap hot;
set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_bar_sortRL,1));

[temp postpeaks_bar]=max([hm_post_bar_sortLR; hm_post_bar_sortRL]');
[temp prepeaks_bar]=max([hm_shockpost_bar_sortLR; hm_shockpost_bar_sortRL]');

    set(gcf,'Name','Figure 3D - heatmaps');   
    
% PIES

figure(34); clf;
subplot(2,3,1); pie(sum(analysis_onlyboth_recur(grp.shock,[2 1 3])));
title(num2str(sum(sum(analysis_onlyboth_recur(grp.shock,[2 1 3])))));
subplot(2,3,2); pie(sum(analysis_onlyboth_recur(grp.scpshk,[2 1 3])));
title(num2str(sum(sum(analysis_onlyboth_recur(grp.scpshk,[2 1 3])))));
subplot(2,3,3); pie(sum(analysis_onlyboth_recur(grp.bar,[2 1 3])));
title(num2str(sum(sum(analysis_onlyboth_recur(grp.bar,[2 1 3])))));

    set(gcf,'Name','Figure 3D - pies');    

%--------------- Figure 3E

figure(35); clf;

subplot(2,3,1);
histogram(prepeaks_df,.5:1:23.5); hold on;
histogram(postpeaks_df,.5:1:23.5); 
[h,p_df,KS_df]=kstest2(prepeaks_df,postpeaks_df); title(p_df);
subplot(2,3,4);
text(1,1,['DF KS: ' num2str(KS_df)]);
text(1,-2,['DF pval: ' num2str(p_df)]);
%text(1,-11,['WTHN df v sc: ' num2str(p_wthn)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);    

subplot(2,3,2);
histogram(prepeaks_sc,.5:1:23.5); hold on;
histogram(postpeaks_sc,.5:1:23.5); 
[h,p_sc,KS_sc]=kstest2(prepeaks_sc,postpeaks_sc); title(p_sc);
subplot(2,3,5);
text(1,1,['SCOP KS: ' num2str(KS_sc)]);
text(1,-2,['SCOP pval: ' num2str(p_sc)]);
%text(1,-11,['WTHN df v sc: ' num2str(p_wthn)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);    

subplot(2,3,3);
histogram(prepeaks_bar,.5:1:23.5); hold on;
histogram(postpeaks_bar,.5:1:23.5); 
[h,p_bar,KS_bar]=kstest2(prepeaks_bar,postpeaks_bar); title(p_bar);
subplot(2,3,6);
text(1,1,['BAR KS: ' num2str(KS_bar)]);
text(1,-2,['BAR pval: ' num2str(p_bar)]);
%text(1,-11,['WTHN df v sc: ' num2str(p_wthn)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);    

    set(gcf,'Name','Figure 3E');
    
%--------------- Figure 3F

figure(36);
clf; 

rvalx=-.975:.05:.975;

cmin=-1; cmax=1;
    
subplot(2,3,1);
imagesc(squeeze(nanmean(Rmatrix_place_shock1_post(shkdex2,:,:),1))); colormap jet;
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
xlabel('drug free shock');

subplot(2,3,2);
imagesc(squeeze(nanmean(Rmatrix_place_scopshock_post(scpdex2,:,:),1))); colormap jet;
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
xlabel('scopolamine shock');

subplot(2,3,3);
imagesc(squeeze(nanmean(Rmatrix_place_barrier_post,1))); colormap jet;
set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
xlabel('barrier no shock');

    set(gcf,'Name','Figure 3F');
    

%--------------- Figure 3H

figure(37);
clf; 

subplot(3,9,2);
[pre_v_cross_safeL] = between_session_analysis_df(20, 33, 'safeL pvc', analysis_results, grp, [-1 1]);
set(gca,'YLim',[-.5 1],'XTick','');
subplot(3,9,1);
text(1,1,['1x3 (traintype): ' num2str(pre_v_cross_safeL.tbl{2,6})]);
text(1,-2,['DF vs SC: ' num2str(pre_v_cross_safeL.dfsc_unpair_cross)]);
text(1,-5,['DF vs BAR: ' num2str(pre_v_cross_safeL.dfbar_upair_cross)]);
text(1,-8,['SC vs BAR: ' num2str(pre_v_cross_safeL.scbar_unpair_cross)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);    

subplot(3,9,11)
[pre_v_cross_centerz] = between_session_analysis_df(20, 23, 'centerz pvc', analysis_results, grp, [-1 1]);
set(gca,'YLim',[-.5 1],'XTick','');
subplot(3,9,10);
text(1,1,['1x3 (traintype): ' num2str(pre_v_cross_centerz.tbl{2,6})]);
text(1,-2,['DF vs SC: ' num2str(pre_v_cross_centerz.dfsc_unpair_cross)]);
text(1,-5,['DF vs BAR: ' num2str(pre_v_cross_centerz.dfbar_upair_cross)]);
text(1,-8,['SC vs BAR: ' num2str(pre_v_cross_centerz.scbar_unpair_cross)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);  

subplot(3,9,20);
[pre_v_cross_safeR] = between_session_analysis_df(18, 19, 'safeR pvc', analysis_results, grp, [-1 1]);
set(gca,'YLim',[-.5 1],'XTick','');
subplot(3,9,19);
text(1,1,['1x3 (traintype): ' num2str(pre_v_cross_safeR.tbl{2,6})]);
text(1,-2,['DF vs SC: ' num2str(pre_v_cross_safeR.dfsc_unpair_cross)]);
text(1,-5,['DF vs BAR: ' num2str(pre_v_cross_safeR.dfbar_upair_cross)]);
text(1,-8,['SC vs BAR: ' num2str(pre_v_cross_safeR.scbar_unpair_cross)]);
set(gca,'XLim',[.8 5],'YLim',[-12 3],'XTick',[]);  

t=table([1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3]',[pre_v_cross_safeL.df_data; pre_v_cross_safeL.sc_data; pre_v_cross_safeL.bar_data],[pre_v_cross_centerz.df_data; pre_v_cross_centerz.sc_data; pre_v_cross_centerz.bar_data],[pre_v_cross_safeR.df_data; pre_v_cross_safeR.sc_data; pre_v_cross_safeR.bar_data],'VariableNames',{'condition','z1','z2','z3'});
Meas = table([1 2 3]','VariableNames',{'Zones'});
rm = fitrm(t,'z1-z3~condition','WithinDesign',Meas);
ranovatbl = ranova(rm)

sex_order=[1 4 8 2 5 6 7 3 9]; %-drug free shock
subplot(3,9,[2 3 4]+1);
try
boxplot(hm_post_shk_safeL(shkdex,:)','Symbol',''); 
set(gca,'XLim',[.5 10.5],'YLim',[-.5 1],'XTick','');
end
subplot(3,9,[10 11 12]+2);
boxplot(hm_post_shk_centerz(shkdex(sex_order),:)','Symbol',''); 
set(gca,'XLim',[.5 10.5],'YLim',[-.5 1],'XTick','');
subplot(3,9,[18 19 20]+3);
try
boxplot(hm_post_shk_safeR(shkdex(sex_order),:)','Symbol',''); 
set(gca,'XLim',[.5 10.5],'YLim',[-.5 1],'XTick','');
end

sex_order=[6 2 5 4 3 1 7]; %--scopoalmine shock

subplot(3,9,[5 6]+1);
try
boxplot(hm_post_scp_safeL(scpdex(sex_order),:)','Symbol',''); 
set(gca,'XLim',[.5 7.5],'YLim',[-.5 1],'XTick','');
end
subplot(3,9,[13 14]+2);
boxplot(hm_post_scp_centerz(scpdex(sex_order),:)','Symbol',''); 
set(gca,'XLim',[.5 7.5],'YLim',[-.5 1],'XTick','');
subplot(3,9,[21 22]+3);
try
boxplot(hm_post_scp_safeR(scpdex(sex_order),:)','Symbol',''); 
set(gca,'XLim',[.5 7.5],'YLim',[-.5 1],'XTick','');
end

sex_order=[5 6 3 4 1 2]; %--scopoalmine shock

subplot(3,9,[5 6]+3);
try
boxplot(hm_post_bar_safeL(sex_order,:)','Symbol',''); 
set(gca,'XLim',[-.5 6.5],'YLim',[-.5 1],'XTick','');
end
subplot(3,9,[13 14]+4);
boxplot(hm_post_bar_centerz(sex_order,:)','Symbol',''); 
set(gca,'XLim',[-.5 6.5],'YLim',[-.5 1],'XTick','');
subplot(3,9,[21 22]+5);
try
boxplot(hm_post_bar_safeR(sex_order,:)','Symbol',''); 
set(gca,'XLim',[-.5 6.5],'YLim',[-.5 1],'XTick','');
end

    set(gcf,'Name','Figure 3H');    

    
    
    
    
    
    
    
% clear p_*;
% for r=1:length(grp.shock)
%     rr=shkdex(r);
%     [p_signrank(r), h]=signrank(hm_pre_shk_centerz(rr,:),hm_post_shk_centerz(rr,:));
% %    [p_ranksum(r), h]=ranksum(hm_pre_shk_centerz(r,:),hm_post_shk_centerz(r,:));
% end
% shk_centerz_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_shk_centerz(shkdex,:)')-nanmedian(hm_post_shk_centerz(shkdex,:)'));
% shk_centerz_dex=-log10(p_signrank) .* sign(nanmedian(hm_pre_shk_centerz(shkdex,:)')-nanmedian(hm_post_shk_centerz(shkdex,:)'));
% subplot(9,6,[32 38]+2+1);
% scatter(shk_centerz_dex,grp.shockbeh);
% [R_shk_centerz,P]=corrcoef(shk_centerz_dex,grp.shockbeh)
% subplot(9,6,32+2-6*4+1); 
% temp1=histogram(hm_pre_shk_centerz(shkdex,:),-1:.05:1,'Normalization','cdf');
% subplot(9,6,32+2-6*2+1); 
% temp2=histogram(hm_post_shk_centerz(shkdex,:),-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*3+4); 
% %plot(rvalx,temp1.Values);
% hold on; plot(rvalx,temp2.Values); 
% set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+5-6*3);
% % boxplot(hm_pre_shk_centerz','ori','horizontal','Symbol',''); title('centerz');
% % set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% % subplot(9,6,[19 25]+5);
% boxplot(hm_post_shk_centerz(shkdex,:)','ori','horizontal','Symbol',''); 
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% violin_shk_centerz=hm_post_shk_centerz(shkdex,:);
% violin_shk_centerz=violin_shk_centerz(:);
% 
% 
% 
% 
% subplot(2,2,3);
% %[pre_v_cross_simR] = between_session_analysis_df(24, 25, 'simR', analysis_results, grp, [0 1]);
% [safeL_v_centerz_pvc] = between_session_analysis_df(19, 23, 'pvc', analysis_results, grp, [0 1]);
% 
% %subplot(2,2,4);
% 
% figure(803); clf;
% % subplot(2,3,1);
% % [pre_v_cross_safeR] = between_session_analysis_df(32, 33, 'all pvc', analysis_results, grp, [-.2 1]);
% %[pre_v_cross_simS] = between_session_analysis_df(24, 39, 'simS', analysis_results, grp, [-5 20]);
% subplot(2,3,1);
% [pre_v_cross_safeL] = between_session_analysis_df(18, 19, 'safeL pvc', analysis_results, grp, [-1 1]);
% subplot(2,3,2)
% [pre_v_cross_centerz] = between_session_analysis_df(20, 23, 'centerz pvc', analysis_results, grp, [-1 1]);
% subplot(2,3,3);
% [pre_v_cross_safeR] = between_session_analysis_df(20, 33, 'safeR pvc', analysis_results, grp, [-1 1]);
% 
% % subplot(2,3,5); 
% a=19; b=23;
% %     bar(1:2,mean(analysis_results([grp.shock1 grp.shock2],[a b]))./[1  1]); hold on;
% %     for i=grp.shock1
% %         plot(1:2,analysis_results(i,[a b])./[1 1 ],'-ok'); hold on;
% %     end
% %     for i=grp.shock2
% %         plot(1:2,analysis_results(i,[a b])./[1 1 ],'-sk'); hold on;
% %     end
% %     
% %     bar(5:6,mean(analysis_results([7:12],[a b]))./[1  1]); hold on;
% %     for i=7:12
% %         plot(5:6,analysis_results(i,[a b])./[1 1 ],'-db'); hold on;
% %     end    
% %     
% %     bar(3:4,mean(analysis_results([17:25],[a b]))./[1  1]); hold on;
% %     for i=grp.scpshk(1:2)
% %         plot(3:4,analysis_results(i,[a b])./[1 1 ],'-or'); hold on;
% %     end    
% %     for i=grp.scpshk(3:7)
% %         plot(3:4,analysis_results(i,[a b])./[1 1 ],'-sr'); hold on;
% %     end    %     for i=17:25
% %     set(gca,'XLim',[0 7],'YLim',[-.2 1]); 
% 
%     between_factors=[ones(1,length(grp.shock)) 2*ones(1,length(grp.scpshk)) 3*ones(1,length(grp.bar))]';
%     
%     datamat=[[analysis_results(grp.shock,a); analysis_results(grp.scpshk,a); analysis_results(grp.bar,a)] ...
%              [analysis_results(grp.shock,b); analysis_results(grp.scpshk,b); analysis_results(grp.bar,b)]];
%              i
%     [safeL_v_centerz.tbl,rm] = simple_mixed_anova(datamat, between_factors, {'sesspair'}, {'traintype'});
% 
%     [h,p]=ttest2(analysis_results(grp.shock,a),analysis_results(grp.shock,b))
%     [h,p]=ttest2(analysis_results(grp.scpshk,a),analysis_results(grp.scpshk,b))
%     [h,p]=ttest2(analysis_results(grp.bar,a),analysis_results(grp.bar,b))
% 
% 
%     subplot(2,3,6); hold off; a=3; b=4;
%     bar(1:2,mean([sum(analysis_beelines([grp.shock1 grp.shock2],[a b])') sum(analysis_beelines([grp.shock1 grp.shock2],[a b]+2)')])./[1  1]); hold on;
%     for i=grp.shock1
%         plot(1:2,[sum(analysis_beelines(i,[a b])) sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-ok');     hold on;
% [sum(analysis_beelines(i,[a b])) sum(analysis_beelines(i,[a b]+2))]
%     end
%     for i=grp.shock2
%         plot(1:2,[sum(analysis_beelines(i,[a b])) sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-sk'); 
%     end
%     
%     bar(5:6,mean([sum(analysis_beelines([7:12],[a b])') sum(analysis_beelines([7:12],[a b]+2)')])./[1  1]); 
%     hold on;
%     for i=7:12
%         plot(5:6,[sum(analysis_beelines(i,[a b])') sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-db'); 
%     end    
%     
%     bar(3:4,mean([sum(analysis_beelines([grp.scpshk],[a b])')])./[1  1]); 
%     hold on;
%     for i=grp.scpshk(1:2)
%         plot(3:4,[sum(analysis_beelines(i,[a b])) sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-or'); 
%     end    
%     for i=grp.scpshk(3:7)
%         plot(3:4,[sum(analysis_beelines(i,[a b])) sum(analysis_beelines(i,[a b]+2))]./[1 1 ],'-sr'); 
%     end    %     for i=17:25
%     set(gca,'XLim',[0 7],'YLim',[0 100]); 
%     
% 
%       between_factors=[ones(1,length(grp.shock)) 2*ones(1,length(grp.scpshk)) 3*ones(1,length(grp.bar))]';
%       datamat=[sum(analysis_beelines([grp.shock1 grp.shock2],[a b])')'; sum(analysis_beelines(grp.scpshk,[a b])')'; sum(analysis_beelines(grp.bar,[a b])')'];
%       [p,tbl,stats] = anova1(datamat,between_factors,'off')    
% 
%     [h,p]=ttest2(sum(analysis_beelines([grp.shock1 grp.shock2],[a b])')',sum(analysis_beelines(grp.scpshk,[a b])')')
%     [h,p]=ttest2(sum(analysis_beelines([grp.shock1 grp.shock2],[a b])')',sum(analysis_beelines(grp.bar,[a b])')')
%     [h,p]=ttest2(sum(analysis_beelines(grp.bar,[a b])')',sum(analysis_beelines(grp.scpshk,[a b])')')
% 
% subplot(2,3,5);
% hold off; a=2; b=3;
%     bar(1:2,mean(analysis_meanspeed([grp.shock1 grp.shock2],[a b]))./[1  1]); hold on;
%     for i=grp.shock1
%         plot(1:2,analysis_meanspeed(i,[a b])./[1 1 ],'-ok');     hold on;
% 
%     end
%     for i=grp.shock2
%         plot(1:2,analysis_meanspeed(i,[a b])./[1 1 ],'-sk'); 
%     end
%     for i=grp.shock2
%         plot(1:2,analysis_meanspeed(i,[a b])./[1 1 ],'-sk'); 
%     end
%     
%     bar(5:6,mean(analysis_meanspeed([7:12],[a b]))./[1  1]); 
%     hold on;
%     for i=7:12
%         plot(5:6,analysis_meanspeed(i,[a b])./[1 1 ],'-db'); 
%     end    
%     
%     bar(3:4,mean(analysis_meanspeed([grp.scpshk],[a b]))./[1  1]); 
%     hold on;
%     for i=grp.scpshk(1:2)
%         plot(3:4,analysis_meanspeed(i,[a b])./[1 1 ],'-or'); 
%     end    
%     for i=grp.scpshk(3:7)
%         plot(3:4,analysis_meanspeed(i,[a b])./[1 1 ],'-sr'); 
%     end    %     for i=17:25
%     set(gca,'XLim',[0 7],'YLim',[0 100]); 
%     
%     between_factors=[ones(1,length(grp.shock)) 2*ones(1,length(grp.scpshk)) 3*ones(1,length(grp.bar))]';
%     
%     datamat=[[analysis_meanspeed(grp.shock,a); analysis_meanspeed(grp.scpshk,a); analysis_meanspeed(grp.bar,a)] ...
%              [analysis_meanspeed(grp.shock,b); analysis_meanspeed(grp.scpshk,b); analysis_meanspeed(grp.bar,b)]];
%              i
%     [pre48_v_post_meanspd.tbl,rm] = simple_mixed_anova(datamat, between_factors, {'sesspair'}, {'traintype'});
% 
%     [h,p]=ttest(analysis_meanspeed(grp.shock,a),analysis_meanspeed(grp.shock,b))
%     [h,p]=ttest(analysis_meanspeed(grp.scpshk,a),analysis_meanspeed(grp.scpshk,b))
%     [h,p]=ttest(analysis_meanspeed(grp.bar,a),analysis_meanspeed(grp.bar,b))
% 
% 
%     
%     
%     
%     
% figure(804); clf;
% 
% subplot(4,3,1);
% scatter(analysis_results([15 27],25),behavshkgroup([8 11]),'*k'); hold on;
% scatter(analysis_results(grp.shock1,25),behavshkgroup([1:7 9]),'ok'); hold on;
% scatter(analysis_results(grp.shock2,25),behavshkgroup(10),'sk');
% set(gca,'XLim',[-1 3.5],'YLim',[-.25 1.75]); axis square;
% subplot(4,3,2);
% scatter(analysis_results(18:19,25),behavscopgroup([2 3]),'*r'); hold on;
% scatter(analysis_results(grp.scpshk(1:2),25),behavscopgroup([1 4]),'r'); hold on;
% scatter(analysis_results(grp.scpshk(3:7),25),behavscopgroup([5:end]),'sr');
% temp=polyfit(analysis_results(grp.scpshk,25),behavscopgroup([1 4 5:end]),1);
% plot([0 3],[temp(2) temp(2)+temp(1)*3]);
% set(gca,'XLim',[-1 3.5],'YLim',[-.25 1.75]); axis square;
% subplot(4,3,3);
% scatter(analysis_results(grpbar,25),behavbargroup,'db');
% set(gca,'XLim',[-1 3.5],'YLim',[-.25 1.75]); axis square;
% 
% 
% subplot(4,3,4);
% scatter(analysis_results([15 27],19),behavshkgroup([8 11]),'*k'); hold on;
% scatter(analysis_results(grp.shock1,19),behavshkgroup([1:7 9]),'ok'); hold on;
% scatter(analysis_results(grp.shock2,19),behavshkgroup(10),'sk');
% set(gca,'XLim',[-.2 1],'YLim',[-.25 1.75]); axis square;
% subplot(4,3,5);
% scatter(analysis_results(18:19,19),behavscopgroup([2 3]),'*r'); hold on;
% scatter(analysis_results(grp.scpshk(1:2),19),behavscopgroup([1 4]),'r'); hold on;
% scatter(analysis_results(grp.scpshk(3:7),19),behavscopgroup([5:end]),'sr');
% set(gca,'XLim',[-.2 1],'YLim',[-.25 1.75]); axis square;
% subplot(4,3,6);
% scatter(analysis_results(grpbar,19),behavbargroup,'db');
% set(gca,'XLim',[-.2 1],'YLim',[-.25 1.75]); axis square;
% 
% subplot(4,3,7);
% scatter(analysis_results([15 27],23),behavshkgroup([8 11]),'*k'); hold on;
% scatter(analysis_results(grp.shock1,23),behavshkgroup([1:7 9]),'ok'); hold on;
% scatter(analysis_results(grp.shock2,23),behavshkgroup(10),'sk');
% set(gca,'XLim',[-.2 1],'YLim',[-.25 1.75]); axis square;
% subplot(4,3,8);
% scatter(analysis_results(18:19,23),behavscopgroup([2 3]),'*r'); hold on;
% scatter(analysis_results(grp.scpshk(1:2),23),behavscopgroup([1 4]),'r'); hold on;
% scatter(analysis_results(grp.scpshk(3:7),23),behavscopgroup(5:end),'sr');
% set(gca,'XLim',[-.2 1],'YLim',[-.25 1.75]); axis square;
% subplot(4,3,9);
% scatter(analysis_results(grpbar,23),behavbargroup,'db');
% set(gca,'XLim',[-.2 1],'YLim',[-.25 1.75]); axis square;
% 
% subplot(4,3,10);
% scatter(analysis_results([15 27],33),behavshkgroup([8 11]),'*k'); hold on;
% scatter(analysis_results(grp.shock1,33),behavshkgroup([1:7 9]),'ok'); hold on;
% scatter(analysis_results(grp.shock2,33),behavshkgroup(10),'sk');
% set(gca,'XLim',[-.2 1],'YLim',[-.25 1.75]); axis square;
% subplot(4,3,11);
% scatter(analysis_results(18:19,33),behavscopgroup([2 3]),'*r'); hold on;
% scatter(analysis_results(grp.scpshk(1:2),33),behavscopgroup([1 4]),'r'); hold on;
% scatter(analysis_results(grp.scpshk(3:7),33),behavscopgroup([5:end]),'sr');
% temp=polyfit(analysis_results(grp.scpshk,33),behavscopgroup([1 4 5:end]),1);
% plot([0 .75],[temp(2) temp(2)+temp(1)*.75]);
% set(gca,'XLim',[-.2 1],'YLim',[-.25 1.75]); axis square;
% subplot(4,3,12);
% scatter(analysis_results(grpbar,33),behavbargroup,'db');
% set(gca,'XLim',[-.2 1],'YLim',[-.25 1.75]); axis square;
% 
% %------------------------ PLOT HEATMAPS-------------------------
% 
% 
% figure(figoff+702); clf;
% 
% hm_post_shk_place=[hm_post_shk(find(hm_pflags_post_shk(:,1) & ismember(hm_post_shk_group,grp.shock)),1:23); hm_post_shk(find(hm_pflags_post_shk(:,2) & ismember(hm_post_shk_group,grp.shock)),24:46)];
% [temp peak]=max(hm_post_shk_place');
% hm_post_shk_sort=sortrows([peak' hm_post_shk_place],1);
% hm_post_shk_sort=hm_post_shk_sort(:,2:end);
% subplot(6,6,[1 7 13]+0); imagesc(1000*conv2(1,[.5 1 .5]/2,hm_post_shk_sort,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_shk_sort,1));
% 
% hm_shockpost_shk_place=[hm_shockpost_shk(find(hm_pflags_post_shk(:,1) & ismember(hm_post_shk_group,grp.shock)),1:23); hm_shockpost_shk(find(hm_pflags_post_shk(:,2) & ismember(hm_post_shk_group,grp.shock)),24:46)];
% hm_shockpost_shk_sort=sortrows([peak' hm_shockpost_shk_place],1);
% hm_shockpost_shk_sort=hm_shockpost_shk_sort(:,2:end);
% subplot(6,6,[1 7 13]+1); imagesc(1000*conv2(1,[.5 1 .5]/2,hm_shockpost_shk_sort,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_shk_sort,1));
% 
% figure(563); clf;
% 
% hm_post_shk_placeLR=[hm_post_shk(find(hm_pflags_post_shk(:,1) & ismember(hm_post_shk_group,grp.shock)),1:23)];
% [temp peakLR]=max(hm_post_shk_placeLR');
% hm_post_shk_sortLR=sortrows([peakLR' hm_post_shk_placeLR],1);
% hm_post_shk_sortLR=hm_post_shk_sortLR(:,2:end);
% 
% hm_post_shk_placeRL=[hm_post_shk(find(hm_pflags_post_shk(:,2) & ismember(hm_post_shk_group,grp.shock)),24:46)];
% [temp peakRL]=max(hm_post_shk_placeRL');
% hm_post_shk_sortRL=sortrows([peakRL' hm_post_shk_placeRL],1);
% hm_post_shk_sortRL=hm_post_shk_sortRL(:,2:end);
% 
% 
% subplot(6,6,[1 7 13]+0); imagesc(1000*conv2(1,1,hm_post_shk_sortLR,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_shk_sortLR,1));
% subplot(6,6,[1 7 13]+18); imagesc(1000*conv2(1,1,hm_post_shk_sortRL,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_shk_sortRL,1));
% 
% hm_shockpost_shk_placeLR=[hm_shockpost_shk(find(hm_pflags_post_shk(:,1) & ismember(hm_post_shk_group,grp.shock)),1:23)];
% hm_shockpost_shk_sortLR=sortrows([peakLR' hm_shockpost_shk_placeLR],1);
% hm_shockpost_shk_sortLR=hm_shockpost_shk_sortLR(:,2:end);
% 
% hm_shockpost_shk_placeRL=[hm_shockpost_shk(find(hm_pflags_post_shk(:,2) & ismember(hm_post_shk_group,grp.shock)),24:46)];
% hm_shockpost_shk_sortRL=sortrows([peakRL' hm_shockpost_shk_placeRL],1);
% hm_shockpost_shk_sortRL=hm_shockpost_shk_sortRL(:,2:end);
% 
% subplot(6,6,[1 7 13]+1); imagesc(1000*conv2(1,1,hm_shockpost_shk_sortLR,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_shk_sortLR,1));
% subplot(6,6,[1 7 13]+19); imagesc(1000*conv2(1,1,hm_shockpost_shk_sortRL,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_shk_sortRL,1));
% 
% [temp postpeaks_df]=max([hm_post_shk_sortLR; hm_post_shk_sortRL]');
% [temp prepeaks_df]=max([hm_shockpost_shk_sortLR; hm_shockpost_shk_sortRL]');
% 
% figure(figoff+702);
% 
% [simRL result_ss] = pv_heatmap_tuningcorr(conv2(1,[.5 1 .5]/2,hm_post_shk_place,'same'),conv2(1,[.5 1 .5]/2,hm_shockpost_shk_place,'same'));
% subplot(6,6,32-7); thist=histogram(result_ss,-5:.5:20,'Normalization','cdf');
% subplot(6,6,32); plot(thist.BinEdges(1:end-1)+.25,thist.Values); hold on;  set(gca,'XLim',[-5 20],'YLim',[0 1]);
% group_ss=[hm_post_shk_group(find(hm_pflags_post_shk(:,1))); hm_post_shk_group(find(hm_pflags_post_shk(:,2)))];
% %result_ss=result_ss(ismember(group_ss,grp.shock));
% group_ss=group_ss(ismember(group_ss,grp.shock));
% subplot(6,6,[20 26]);
% boxplot(result_ss,group_ss,'ori','horizontal','Symbol',''); 
% analysis_results(grp.shock,39) = groupsummary(result_ss(:),group_ss(:),'median');
% set(gca,'XLim',[-5 20],'YLim',[.5 11.5]); title(sum(result_ss>-log10(.05))/length(result_ss)); 
% subplot(6,6,31);
% scatter(analysis_results(grp.shock,39),grp.shockbeh);
% set(gca,'XLim',[-.5 4],'YLim',[0 2]);
% violin_shk_simS=result_ss;
% 
% hm_post_scpshk_place=[hm_post_scpshk(find(hm_pflags_post_scpshk(:,1) & ismember(hm_post_scpshk_group,grp.scpshk)),1:23); hm_post_scpshk(find(hm_pflags_post_scpshk(:,2) & ismember(hm_post_scpshk_group,grp.scpshk)),24:46)];
% [temp peak]=max(hm_post_scpshk_place');
% hm_post_scpshk_sort=sortrows([peak' hm_post_scpshk_place],1);
% hm_post_scpshk_sort=hm_post_scpshk_sort(:,2:end);
% subplot(6,6,[1 7 13]+2); imagesc(1000*conv2(1,[.5 1 .5]/2,hm_post_scpshk_sort,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_scpshk_sort,1));
% 
% hm_shockpost_scpshk_place=[hm_shockpost_scpshk(find(hm_pflags_post_scpshk(:,1) & ismember(hm_post_scpshk_group,grp.scpshk)),1:23); hm_shockpost_scpshk(find(hm_pflags_post_scpshk(:,2) & ismember(hm_post_scpshk_group,grp.scpshk)),24:46)];
% %hm_shockpost_scpshk_place=hm_shockpost_scpshk(find(sum(hm_pflags_post_scpshk')),:);
% hm_shockpost_scpshk_sort=sortrows([peak' hm_shockpost_scpshk_place],1);
% hm_shockpost_scpshk_sort=hm_shockpost_scpshk_sort(:,2:end);
% subplot(6,6,[1 7 13]+3); imagesc(1000*conv2(1,[.5 1 .5]/2,hm_shockpost_scpshk_sort,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_scpshk_sort,1));
% 
% figure(563); 
% 
% hm_post_scpshk_placeLR=[hm_post_scpshk(find(hm_pflags_post_scpshk(:,1) & ismember(hm_post_scpshk_group,grp.scpshk)),1:23)];
% [temp peakLR]=max(hm_post_scpshk_placeLR');
% hm_post_scpshk_sortLR=sortrows([peakLR' hm_post_scpshk_placeLR],1);
% hm_post_scpshk_sortLR=hm_post_scpshk_sortLR(:,2:end);
% 
% hm_post_scpshk_placeRL=[hm_post_scpshk(find(hm_pflags_post_scpshk(:,2) & ismember(hm_post_scpshk_group,grp.scpshk)),24:46)];
% [temp peakRL]=max(hm_post_scpshk_placeRL');
% hm_post_scpshk_sortRL=sortrows([peakRL' hm_post_scpshk_placeRL],1);
% hm_post_scpshk_sortRL=hm_post_scpshk_sortRL(:,2:end);
% 
% subplot(6,6,[1 7 13]+2); imagesc(1000*conv2(1,1,hm_post_scpshk_sortLR,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_scpshk_sortLR,1));
% subplot(6,6,[1 7 13]+20); imagesc(1000*conv2(1,1,hm_post_scpshk_sortRL,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_scpshk_sortRL,1));
% 
% hm_shockpost_scpshk_placeLR=[hm_shockpost_scpshk(find(hm_pflags_post_scpshk(:,1) & ismember(hm_post_scpshk_group,grp.scpshk)),1:23)];
% hm_shockpost_scpshk_sortLR=sortrows([peakLR' hm_shockpost_scpshk_placeLR],1);
% hm_shockpost_scpshk_sortLR=hm_shockpost_scpshk_sortLR(:,2:end);
% 
% hm_shockpost_scpshk_placeRL=[hm_shockpost_scpshk(find(hm_pflags_post_scpshk(:,2) & ismember(hm_post_scpshk_group,grp.scpshk)),24:46)];
% hm_shockpost_scpshk_sortRL=sortrows([peakRL' hm_shockpost_scpshk_placeRL],1);
% hm_shockpost_scpshk_sortRL=hm_shockpost_scpshk_sortRL(:,2:end);
% 
% subplot(6,6,[1 7 13]+3); imagesc(1000*conv2(1,1,hm_shockpost_scpshk_sortLR,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_scpshk_sortLR,1));
% subplot(6,6,[1 7 13]+21); imagesc(1000*conv2(1,1,hm_shockpost_scpshk_sortRL,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_scpshk_sortRL,1));
% 
% [temp postpeaks_sc]=max([hm_post_scpshk_sortLR; hm_post_scpshk_sortRL]');
% [temp prepeaks_sc]=max([hm_shockpost_scpshk_sortLR; hm_shockpost_scpshk_sortRL]');
% 
% 
% figure(figoff+702);
% 
% 
% [simRL result_ss] = pv_heatmap_tuningcorr(conv2(1,[.5 1 .5]/2,hm_post_scpshk_place,'same'),conv2(1,[.5 1 .5]/2,hm_shockpost_scpshk_place,'same'));
% subplot(6,6,32+2-7); thist=histogram(result_ss,-5:.5:20,'Normalization','cdf');
% subplot(6,6,32); plot(thist.BinEdges(1:end-1)+.25,thist.Values);  set(gca,'XLim',[-5 20],'YLim',[0 1]);
% group_ss=[hm_post_scpshk_group(find(hm_pflags_post_scpshk(:,1))); hm_post_scpshk_group(find(hm_pflags_post_scpshk(:,2)))];
% subplot(6,6,[20 26]+2);
% %result_ss=result_ss(ismember(group_ss,grp.scpshk));
% group_ss=group_ss(ismember(group_ss,grp.scpshk));
% boxplot(result_ss,group_ss,'ori','horizontal','Symbol',''); 
% analysis_results(grp.scpshk,39) = groupsummary(result_ss(:),group_ss(:),'median');
% set(gca,'XLim',[-5 20],'YLim',[.5 11.5]); title(sum(result_ss>-log10(.05))/length(result_ss)); 
% subplot(6,6,33);
% scatter(analysis_results(grp.scpshk,39),grp.scpshkbeh);
% temp=polyfit(analysis_results(grp.scpshk,39),grp.scpshkbeh,1);
% hold on; plot([-.5 4],temp(2)+[[-.5 4]]*temp(1),'r');
% set(gca,'XLim',[-.5 4],'YLim',[0 2]);
% violin_scp_simS=result_ss;
% 
% 
% hm_post_bar_place=[hm_post_bar(find(hm_pflags_post_bar(:,1)),1:23); hm_post_bar(find(hm_pflags_post_bar(:,2)),24:46)];
% [temp peak]=max(hm_post_bar_place');
% hm_post_bar_sort=sortrows([peak' hm_post_bar_place],1);
% hm_post_bar_sort=hm_post_bar_sort(:,2:end);
% subplot(6,6,[1 7 13]+4); imagesc(1000*conv2(1,[.5 1 .5]/2,hm_post_bar_sort,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_bar_sort,1));
% 
% hm_shockpost_bar_place=[hm_shockpost_bar(find(hm_pflags_post_bar(:,1)),1:23); hm_shockpost_bar(find(hm_pflags_post_bar(:,2)),24:46)];
% hm_shockpost_bar_sort=sortrows([peak' hm_shockpost_bar_place],1);
% hm_shockpost_bar_sort=hm_shockpost_bar_sort(:,2:end);
% subplot(6,6,[1 7 13]+5); imagesc(1000*conv2(1,[.5 1 .5]/2,hm_shockpost_bar_sort,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_bar_sort,1));
% 
% figure(563); 
% 
% hm_post_bar_placeLR=[hm_post_bar(find(hm_pflags_post_bar(:,1) & ismember(hm_post_bar_group,grp.bar)),1:23)];
% [temp peakLR]=max(hm_post_bar_placeLR');
% hm_post_bar_sortLR=sortrows([peakLR' hm_post_bar_placeLR],1);
% hm_post_bar_sortLR=hm_post_bar_sortLR(:,2:end);
% 
% hm_post_bar_placeRL=[hm_post_bar(find(hm_pflags_post_bar(:,2) & ismember(hm_post_bar_group,grp.bar)),24:46)];
% [temp peakRL]=max(hm_post_bar_placeRL');
% hm_post_bar_sortRL=sortrows([peakRL' hm_post_bar_placeRL],1);
% hm_post_bar_sortRL=hm_post_bar_sortRL(:,2:end);
% 
% subplot(6,6,[1 7 13]+4); imagesc(1000*conv2(1,1,hm_post_bar_sortLR,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_bar_sortLR,1));
% subplot(6,6,[1 7 13]+22); imagesc(1000*conv2(1,1,hm_post_bar_sortRL,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_post_bar_sortRL,1));
% 
% hm_shockpost_bar_placeLR=[hm_shockpost_bar(find(hm_pflags_post_bar(:,1) & ismember(hm_post_bar_group,grp.bar)),1:23)];
% hm_shockpost_bar_sortLR=sortrows([peakLR' hm_shockpost_bar_placeLR],1);
% hm_shockpost_bar_sortLR=hm_shockpost_bar_sortLR(:,2:end);
% 
% hm_shockpost_bar_placeRL=[hm_shockpost_bar(find(hm_pflags_post_bar(:,2) & ismember(hm_post_bar_group,grp.bar)),24:46)];
% hm_shockpost_bar_sortRL=sortrows([peakRL' hm_shockpost_bar_placeRL],1);
% hm_shockpost_bar_sortRL=hm_shockpost_bar_sortRL(:,2:end);
% 
% subplot(6,6,[1 7 13]+5); imagesc(1000*conv2(1,1,hm_shockpost_bar_sortLR,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_bar_sortLR,1));
% subplot(6,6,[1 7 13]+23); imagesc(1000*conv2(1,1,hm_shockpost_bar_sortRL,'same')); colormap hot;
% set(gca,'XTick',[],'YTick',[]); caxis([0 11]); ylabel(size(hm_shockpost_bar_sortRL,1));
% 
% [temp postpeaks_bar]=max([hm_post_bar_sortLR; hm_post_bar_sortRL]');
% [temp prepeaks_bar]=max([hm_shockpost_bar_sortLR; hm_shockpost_bar_sortRL]');
% 
% figure(figoff+702);
% 
% [simRL result_ss] = pv_heatmap_tuningcorr(conv2(1,[.5 1 .5]/2,hm_post_bar_place,'same'),conv2(1,[.5 1 .5]/2,hm_shockpost_bar_place,'same'));
% subplot(6,6,32+4-7); thist=histogram(result_ss,-5:.5:20,'Normalization','cdf');
% subplot(6,6,32); plot(thist.BinEdges(1:end-1)+.25,thist.Values); set(gca,'XLim',[-5 20],'YLim',[0 1]);
% group_ss=[hm_post_bar_group(find(hm_pflags_post_bar(:,1))); hm_post_bar_group(find(hm_pflags_post_bar(:,2)))];
% subplot(6,6,[20 26]+4);
% boxplot(result_ss,group_ss,'ori','horizontal','Symbol',''); 
% analysis_results(grp.bar,39) = groupsummary(result_ss(:),group_ss(:),'median');
% set(gca,'XLim',[-5 20],'YLim',[.5 11.5]); title(sum(result_ss>-log10(.05))/length(result_ss));
% subplot(6,6,35);
% scatter(analysis_results(grp.bar,39),behavbargroup);
% set(gca,'XLim',[-.5 4],'YLim',[0 2]);
% violin_bar_simS=result_ss;
% 
% cd ..
% figure(100); clf; 
% subplot(1,4,1);
% violinplot([violin_shk_simS'; violin_scp_simS'; violin_bar_simS'],[1+0*violin_shk_simS'; 2+0*violin_scp_simS'; 3+0*violin_bar_simS']);
% set(gca,'YLim',[-10 25]);
% 
% [kwp, kwtbl, kwstats]=kruskalwallis([violin_shk_simS'; violin_scp_simS'; violin_bar_simS'],[1+0*violin_shk_simS'; 2+0*violin_scp_simS'; 3+0*violin_bar_simS'],'off');
% [mup,muh,mustats] = ranksum(violin_shk_simS', violin_scp_simS')
% [mup,muh,mustats] = ranksum(violin_shk_simS', violin_bar_simS')
% [mup,muh,mustats] = ranksum(violin_scp_simS', violin_bar_simS')
% 
% %----------------------- DRUG FREE --------------------------
% 
% figure(601);
% clf; 
% 
% rvalx=-.975:.05:.975;
% 
% cmin=-1; cmax=1;
%     
% clear p_*;
% for r=1:length(grp.shock)
%     rr=shkdex(r);
%     [p_signrank(r), h]=signrank(hm_pre_shk_safeR(rr,:),hm_post_shk_safeR(rr,:));
% %    [p_ranksum(r), h]=ranksum(hm_pre_shk_safeR(r,:),hm_post_shk_safeR(r,:));
% end
% shk_safeR_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_shk_safeR(shkdex,:)')-nanmedian(hm_post_shk_safeR(shkdex,:)'));
% shk_safeR_dex=-log10(p_signrank) .* sign(nanmedian(hm_pre_shk_safeR(shkdex,:)')-nanmedian(hm_post_shk_safeR(shkdex,:)'));
% subplot(9,6,32-6*4-1); 
% temp1=histogram(hm_pre_shk_safeR(shkdex,:),-1:.05:1,'Normalization','cdf'); title(['inc ' num2str(sum(shk_safeR_sig>0))]);
% subplot(9,6,32-6*2-1); 
% temp2=histogram(hm_post_shk_safeR(shkdex,:),-1:.05:1,'Normalization','cdf'); title(['dec ' num2str(sum(shk_safeR_sig>0))]);
% subplot(9,6,32-6*3); 
% %plot(rvalx,temp1.Values);
% hold on; plot(rvalx,temp2.Values); 
% set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+1-6*3);
% % boxplot(hm_pre_shk_safeR','ori','horizontal','Symbol',''); title('safeR');
% % set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% % subplot(9,6,[19 25]+1);
% boxplot(hm_post_shk_safeR(shkdex,:)','ori','horizontal','Symbol',''); 
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% subplot(9,6,[31 32 37 38]);
% imagesc(squeeze(nanmean(Rmatrix_place_shock1_pre,1))); colormap jet;
% set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
% subplot(9,6,[31 32 37 38]+12);
% imagesc(squeeze(nanmean(Rmatrix_place_shock1_post(shkdex2,:,:),1))); colormap jet;
% set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
% violin_shk_safeR=hm_post_shk_safeR(shkdex,:);
% violin_shk_safeR=violin_shk_safeR(:);
% 
% clear p_*;
% for r=1:length(grp.shock)
%     rr=shkdex(r);
%     [p_signrank(r), h]=signrank(hm_pre_shk_safeL(rr,:),hm_post_shk_safeL(rr,:));
% %    [p_ranksum(r), h]=ranksum(hm_pre_shk_safeL(r,:),hm_post_shk_safeL(r,:));
% end
% shk_safeL_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_shk_safeL(shkdex,:)')-nanmedian(hm_post_shk_safeL(shkdex,:)'));
% shk_safeL_dex=-log10(p_signrank) .* sign(nanmedian(hm_pre_shk_safeL(shkdex,:)')-nanmedian(hm_post_shk_safeL(shkdex,:)'));
% subplot(9,6,[32 38]+1);
% scatter(shk_safeL_dex,grp.shockbeh);
% [R_shk_safeL,P]=corrcoef(shk_safeL_dex,grp.shockbeh)
% subplot(9,6,32-6*4+1); 
% temp1=histogram(hm_pre_shk_safeL(shkdex,:),-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*2+1); 
% temp2=histogram(hm_post_shk_safeL(shkdex,:),-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*3+2); 
% %plot(rvalx,temp1.Values);
% hold on; plot(rvalx,temp2.Values); 
% set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+3-6*3);
% % boxplot(hm_pre_shk_safeL','ori','horizontal','Symbol',''); title('safeL');
% % set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% % subplot(9,6,[19 25]+3);
% boxplot(hm_post_shk_safeL(shkdex,:)','ori','horizontal','Symbol','');
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% violin_shk_safeL=hm_post_shk_safeL(shkdex,:);
% violin_shk_safeL=violin_shk_safeL(:);
% 
% 
% clear p_*;
% for r=1:length(grp.shock)
%     rr=shkdex(r);
%     [p_signrank(r), h]=signrank(hm_pre_shk_centerz(rr,:),hm_post_shk_centerz(rr,:));
% %    [p_ranksum(r), h]=ranksum(hm_pre_shk_centerz(r,:),hm_post_shk_centerz(r,:));
% end
% shk_centerz_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_shk_centerz(shkdex,:)')-nanmedian(hm_post_shk_centerz(shkdex,:)'));
% shk_centerz_dex=-log10(p_signrank) .* sign(nanmedian(hm_pre_shk_centerz(shkdex,:)')-nanmedian(hm_post_shk_centerz(shkdex,:)'));
% subplot(9,6,[32 38]+2+1);
% scatter(shk_centerz_dex,grp.shockbeh);
% [R_shk_centerz,P]=corrcoef(shk_centerz_dex,grp.shockbeh)
% subplot(9,6,32+2-6*4+1); 
% temp1=histogram(hm_pre_shk_centerz(shkdex,:),-1:.05:1,'Normalization','cdf');
% subplot(9,6,32+2-6*2+1); 
% temp2=histogram(hm_post_shk_centerz(shkdex,:),-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*3+4); 
% %plot(rvalx,temp1.Values);
% hold on; plot(rvalx,temp2.Values); 
% set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+5-6*3);
% % boxplot(hm_pre_shk_centerz','ori','horizontal','Symbol',''); title('centerz');
% % set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% % subplot(9,6,[19 25]+5);
% boxplot(hm_post_shk_centerz(shkdex,:)','ori','horizontal','Symbol',''); 
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% violin_shk_centerz=hm_post_shk_centerz(shkdex,:);
% violin_shk_centerz=violin_shk_centerz(:);
% 
% %----------------------- SCOPOLAMINE --------------------------
% 
% figure(602);
% clf; 
% 
% rvalx=-.975:.05:.975;
% 
% cmin=-1; cmax=1;
% 
% clear p_*;
% for r=1:length(grp.scpshk)
%     rr=scpdex(r);
%     [p_signrank(r), h]=signrank(hm_pre_scp_safeR(rr,:),hm_post_scp_safeR(rr,:));
% %    [p_ranksum(r), h]=ranksum(hm_pre_scp_safeR(r,:),hm_post_scp_safeR(r,:));
% end
% scp_safeR_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_scp_safeR(scpdex,:)')-nanmedian(hm_post_scp_safeR(scpdex,:)'));
% subplot(9,6,32-6*4-1); 
% temp1=histogram(hm_pre_scp_safeR(scpdex,:),-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*2-1); 
% temp2=histogram(hm_post_scp_safeR(scpdex,:),-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*3); 
% %plot(rvalx,temp1.Values);
% hold on; plot(rvalx,temp2.Values); 
% set(gca,'XLim',[-1 1]);
% % subplot(9,6,[19 25]+1-6*3);
% % boxplot(hm_pre_scp_safeR','ori','horizontal','Symbol',''); title('safeR');
% % set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% % %subplot(9,6,[19 25]+1);
% subplot(9,6,[19 25]+1-6*3);
% boxplot(hm_post_scp_safeR(scpdex,:)','ori','horizontal','Symbol',''); 
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% subplot(9,6,[31 32 37 38]+2);
% imagesc(squeeze(nanmean(Rmatrix_place_scopshock_pre,1))); colormap jet;
% set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
% subplot(9,6,[31 32 37 38]+12+2);
% imagesc(squeeze(nanmean(Rmatrix_place_scopshock_post(scpdex2,:,:),1))); colormap jet;
% set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
% violin_scp_safeR=hm_post_scp_safeR(scpdex,:);
% violin_scp_safeR=violin_scp_safeR(:);
% 
% 
% clear p_*;
% for r=1:length(grp.scpshk)
%     rr=scpdex(r);
%     [p_signrank(r), h]=signrank(hm_pre_scp_safeL(rr,:),hm_post_scp_safeL(rr,:));
% %    [p_ranksum(r), h]=ranksum(hm_pre_scp_safeL(r,:),hm_post_scp_safeL(r,:));
% end
% scp_safeL_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_scp_safeL(scpdex,:)')-nanmedian(hm_post_scp_safeL(scpdex,:)'));
% subplot(9,6,32-6*4+1); 
% temp1=histogram(hm_pre_scp_safeL(scpdex,:),-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*2+1); 
% temp2=histogram(hm_post_scp_safeL(scpdex,:),-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*3+2); 
% %plot(rvalx,temp1.Values);
% hold on; plot(rvalx,temp2.Values); 
% set(gca,'XLim',[-1 1]);
% % subplot(9,6,[19 25]+3-6*3);
% % boxplot(hm_pre_scp_safeL','ori','horizontal','Symbol',''); title('safeL');
% % set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% subplot(9,6,[19 25]+3-6*3);
% boxplot(hm_post_scp_safeL(scpdex,:)','ori','horizontal','Symbol','');
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% violin_scp_safeL=hm_post_scp_safeL(scpdex,:);
% violin_scp_safeL=violin_scp_safeL(:);
% 
% clear p_*;
% for r=1:length(grp.scpshk)
%     rr=scpdex(r);
%     [p_signrank(r), h]=signrank(hm_pre_scp_centerz(rr,:),hm_post_scp_centerz(rr,:));
% %    [p_ranksum(r), h]=ranksum(hm_pre_scp_centerz(r,:),hm_post_scp_centerz(r,:));
% end
% scp_centerz_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_scp_centerz(scpdex,:)')-nanmedian(hm_post_scp_centerz(scpdex,:)'));
% subplot(9,6,32+2-6*4+1); 
% temp1=histogram(hm_pre_scp_centerz(scpdex,:),-1:.05:1,'Normalization','cdf');
% subplot(9,6,32+2-6*2+1); 
% temp2=histogram(hm_post_scp_centerz(scpdex,:),-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*3+4); 
% %plot(rvalx,temp1.Values);
% hold on; plot(rvalx,temp2.Values); 
% set(gca,'XLim',[-1 1]);
% subplot(9,6,[19 25]+5-6*3);
% % boxplot(hm_pre_scp_centerz','ori','horizontal','Symbol',''); title('centerz');
% % set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% % subplot(9,6,[19 25]+5-6*3);
% boxplot(hm_post_scp_centerz(scpdex,:)','ori','horizontal','Symbol',''); 
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% violin_scp_centerz=hm_post_scp_centerz(scpdex,:);
% violin_scp_centerz=violin_scp_centerz(:);
% 
% %----------------------- BARRIER --------------------------
% 
% figure(603);
% clf; 
% 
% rvalx=-.975:.05:.975;
% 
% cmin=-1; cmax=1;
% 
% clear p_*;
% for r=1:6
%     [p_signrank(r), h]=signrank(hm_pre_bar_safeR(r,:),hm_post_bar_safeR(r,:));
% %    [p_ranksum(r), h]=ranksum(hm_pre_bar_safeR(r,:),hm_post_bar_safeR(r,:));
% end
% bar_safeR_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_bar_safeR')-nanmedian(hm_post_bar_safeR'));
% subplot(9,6,32-6*4-1); 
% temp1=histogram(hm_pre_bar_safeR,-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*2-1); 
% temp2=histogram(hm_post_bar_safeR,-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*3); 
% %plot(rvalx,temp1.Values);
% hold on; plot(rvalx,temp2.Values); 
% set(gca,'XLim',[-1 1]);
% % subplot(9,6,[19 25]+1-6*3);
% % boxplot(hm_pre_bar_safeR','ori','horizontal','Symbol',''); 
% % set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% subplot(9,6,[19 25]+1);
% subplot(9,6,[32 38]+0);
% boxplot(hm_post_bar_safeR','ori','horizontal','Symbol',''); 
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]); title('safeR');
% subplot(9,6,[31 32 37 38]+4);
% imagesc(squeeze(nanmean(Rmatrix_place_barrier_pre,1))); colormap jet;
% set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
% subplot(9,6,[31 32 37 38]+12+4);
% imagesc(squeeze(nanmean(Rmatrix_place_barrier_post,1))); colormap jet;
% set(gca,'XTick',[],'YTick',[]); caxis(clim); axis square;
% violin_bar_safeR=hm_post_bar_safeR(:);
% 
% 
% clear p_*;
% for r=1:6
%     [p_signrank(r), h]=signrank(hm_pre_bar_safeL(r,:),hm_post_bar_safeL(r,:));
% %    [p_ranksum(r), h]=ranksum(hm_pre_bar_safeL(r,:),hm_post_bar_safeL(r,:));
% end
% bar_safeL_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_bar_safeL')-nanmedian(hm_post_bar_safeL'));
% subplot(9,6,32-6*4+1); 
% temp1=histogram(hm_pre_bar_safeL,-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*2+1); 
% temp2=histogram(hm_post_bar_safeL,-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*3+2); 
% %plot(rvalx,temp1.Values);
% hold on; plot(rvalx,temp2.Values); 
% set(gca,'XLim',[-1 1]);
% % subplot(9,6,[19 25]+3-6*3);
% % boxplot(hm_pre_bar_safeL','ori','horizontal','Symbol',''); 
% % set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% subplot(9,6,[19 25]+3);
% subplot(9,6,[32 38]+2);
% boxplot(hm_post_bar_safeL','ori','horizontal','Symbol','');
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]); title('safeL');
% violin_bar_safeL=hm_post_bar_safeL(:);
% 
% clear p_*;
% for r=1:6
%     [p_signrank(r), h]=signrank(hm_pre_bar_centerz(r,:),hm_post_bar_centerz(r,:));
% %    [p_ranksum(r), h]=ranksum(hm_pre_bar_centerz(r,:),hm_post_bar_centerz(r,:));
% end
% bar_centerz_sig=-(p_signrank<.01) .* sign(nanmedian(hm_pre_bar_centerz')-nanmedian(hm_post_bar_centerz'));
% subplot(9,6,32+2-6*4+1); 
% temp1=histogram(hm_pre_bar_centerz,-1:.05:1,'Normalization','cdf');
% subplot(9,6,32+2-6*2+1); 
% temp2=histogram(hm_post_bar_centerz,-1:.05:1,'Normalization','cdf');
% subplot(9,6,32-6*3+4); 
% %plot(rvalx,temp1.Values);
% hold on; plot(rvalx,temp2.Values); 
% set(gca,'XLim',[-1 1]);
% % subplot(9,6,[19 25]+5-6*3);
% % boxplot(hm_pre_bar_centerz','ori','horizontal','Symbol',''); 
% % set(gca,'XLim',[-1 1],'YLim',[.5 11.5]);
% subplot(9,6,[19 25]+5);
% subplot(9,6,[32 38]+4);
% boxplot(hm_post_bar_centerz','ori','horizontal','Symbol',''); 
% set(gca,'XLim',[-1 1],'YLim',[.5 11.5]); title('centerz');
% violin_bar_centerz=hm_post_bar_centerz(:);
% 
% figure(1000); clf;
% for i=1:12
%     subplot(3,4,i);
%     imagesc(squeeze(Rmatrix_place_barrier_post(i,:,:)));
%     caxis([-1 1]); colormap jet;
% end
% 
% temp1=hm_post_shk_safeL(shkdex,:);
% temp2=hm_post_scp_safeL(scpdex,:);
% temp3=hm_post_bar_safeL;
% temp4=hm_post_shk_centerz(shkdex,:);
% temp5=hm_post_scp_centerz(scpdex,:);
% temp6=hm_post_bar_centerz;
% 
% figure(100); 
% subplot(1,4,2);
% violinplot([violin_shk_safeR; violin_scp_safeR; violin_bar_safeR],[1+0*violin_shk_safeR; 2+0*violin_scp_safeR; 3+0*violin_bar_safeR]);
% 
% figure(100); 
% subplot(1,4,3);
% violinplot([violin_shk_safeL; violin_scp_safeL; violin_bar_safeL],[1+0*violin_shk_safeL; 2+0*violin_scp_safeL; 3+0*violin_bar_safeL]);
% 
% figure(100); 
% subplot(1,4,4);
% violinplot([violin_shk_centerz; violin_scp_centerz; violin_bar_centerz],[1+0*violin_shk_centerz; 2+0*violin_scp_centerz; 3+0*violin_bar_centerz]);
% 
% % goodshk=find(~isnan(violin_shk_safeL));
% % goodscp=find(~isnan(violin_scp_safeL));
% % goodbar=find(~isnan(violin_bar_safeL));
% 
% 
% 
% %      between_factors=[ones(1,length(violin_shk_safeL)) 2*ones(1,length(violin_scp_safeL)) 3*ones(1,length(violin_bar_safeL))]';
% %      datamat=[violin_shk_safeL; violin_scp_safeL; violin_bar_safeL];
% %      [p,tbl,stats] = anova1(datamat,between_factors)
% % 
% %      [h,p,ci,stats] = ttest2(violin_shk_safeL,violin_scp_safeL)
% %      [h,p,ci,stats] = ttest2(violin_shk_safeL,violin_bar_safeL)
% %      [h,p,ci,stats] = ttest2(violin_scp_safeL,violin_bar_safeL)
%      
% %     
% %     datamat=[[violin_shk_safeL; violin_scp_safeL; violin_bar_safeL] ...
% %              [violin_shk_centerz; violin_scp_centerz; violin_bar_centerz]];
%              
%  %   [tbl,rm] = simple_mixed_anova(datamat, between_factors, {'sesspair'}, {'traintype'});
% 
% set(gca,'YLim',[-10 25]);
% 
% %     between_factors=[ones(1,length(temp1(:))) 2*ones(1,length(temp2(:))) 3*ones(1,length(temp3(:)))]';
% %     
% %     datamat=[[temp1(:); temp2(:); temp3(:)] ...
% %              [temp4(:); temp5(:); temp6(:)]];
% %              
% %     [safeL_vs_centerz.tbl,rm] = simple_mixed_anova(datamat, between_factors, {'safeLcenterz'}, {'traintype'});
