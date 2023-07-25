function [pre_vs_cross] = between_session_analysis_df(a,b,lbl,analysis_results, grp, ylim)

shklist=[];
    bar(1:2,mean(analysis_results([grp.shock1 grp.shock2],[a b]))./[1  1]); hold on;
    for i=grp.shock1
        scatter(2,analysis_results(i,[ b])/[1  ],'ok'); hold on;
        shklist=[shklist analysis_results(i,[ b])];
    end
    for i=grp.shock2
        scatter(2,analysis_results(i,[ b])/[1  ],'sk'); hold on;
        shklist=[shklist analysis_results(i,[ b])];
    end

    barlist=[];
    bar(5:6,mean(analysis_results([7:12],[a b]))./[1  1]); hold on;
    for i=7:12
        scatter(6,analysis_results(i,[ b])/[1  ],'db'); hold on;
        barlist=[barlist analysis_results(i,[ b])];
    end    
    
    scoplist=[];
%     for i=[17:20]
%         scatter(4,analysis_results(i,[ b])/[1  ],'or'); hold on;
%         scoplist=[scoplist analysis_results(i,[ b])];
%     end    
    for i=grp.scpshk
        scatter(4,analysis_results(i,[ b])/[1  ],'sr'); hold on;
        scoplist=[scoplist analysis_results(i,[ b])];
    end    %     for i=17:25
    bar([1 3 5],[mean(shklist) mean(scoplist) mean(barlist)]);
    set(gca,'XLim',[0 7],'YLim',ylim); title(lbl);
    
%     clf;
%     scatter(grp.shockbeh,analysis_results(grp.shock,a));
    
     pre_vs_cross.df_means=nanmean(analysis_results(grp.shock,[b])); 
     pre_vs_cross.sc_means=nanmean(analysis_results(grp.scpshk,[b]));
     pre_vs_cross.bar_means=nanmean(analysis_results(grp.bar,[b]));

     pre_vs_cross.df_data=(analysis_results(grp.shock,[b])); 
     pre_vs_cross.sc_data=(analysis_results(grp.scpshk,[b]));
     pre_vs_cross.bar_data=(analysis_results(grp.bar,[b]));

%     [p1, h]=signrank(analysis_results(grp.shock,a),analysis_results(grp.shock,b));
%     [p2, h]=signrank(analysis_results(grp.shock,c),analysis_results(grp.shock,b));
%     [p3, h]=signrank(analysis_results(grp.shock,c),analysis_results(grp.shock,a));
%     [p1 p2 p3]
%     [p1, h]=signrank(analysis_results(grp.scpshk,a),analysis_results(grp.scpshk,b));
%     [p2, h]=signrank(analysis_results(grp.scpshk,c),analysis_results(grp.scpshk,b));
%     [p3, h]=signrank(analysis_results(grp.scpshk,c),analysis_results(grp.scpshk,a));
%     [p1 p2 p3]
%     [p1, h]=signrank(analysis_results(grp.bar,a),analysis_results(grp.bar,b));
%     [p2, h]=signrank(analysis_results(grp.bar,c),analysis_results(grp.bar,b));
%     [p3, h]=signrank(analysis_results(grp.bar,c),analysis_results(grp.bar,a));
%     [p1 p2 p3]    
%     [p1, h]=signrank(mean(analysis_results(grp.shock,[a c])'),analysis_results(grp.shock,b))
%     [p1, h]=signrank(mean(analysis_results(grp.scpshk,[a c])'),analysis_results(grp.scpshk,b))
%     [p1, h]=signrank(mean(analysis_results(grp.bar,[a c])'),analysis_results(grp.bar,b))

%pre-training vs cross-training
%     [pre_vs_cross.df_p, h]=signrank(analysis_results(grp.shock,a),analysis_results(grp.shock,b));
%     [pre_vs_cross.sc_p, h]=signrank(analysis_results(grp.scpshk,a),analysis_results(grp.scpshk,b));
%     [pre_vs_cross.bar_p, h]=signrank(analysis_results(grp.bar,a),analysis_results(grp.bar,b));

% [R,P]= corrcoef(grp.shockbeh,analysis_results(grp.shock,a)); pre_vs_cross.dfpre_corrp=P(1,2); pre_vs_cross.dfpre_corrR=R(1,2);
% [R,P]= corrcoef(grp.scpshkbeh,analysis_results(grp.scpshk,a)); pre_vs_cross.scpre_corrp=P(1,2); pre_vs_cross.scpre_corrR=R(1,2);
% [R,P]= corrcoef(grp.barbeh,analysis_results(grp.bar,a)); pre_vs_cross.barpre_corrp=P(1,2); pre_vs_cross.barpre_corrR=R(1,2);

[R,P]= corrcoef(grp.shockbeh,analysis_results(grp.shock,b)); pre_vs_cross.dfcross_corrp=P(1,2); pre_vs_cross.dfcross_corrR=R(1,2);
[R,P]= corrcoef(grp.scpshkbeh,analysis_results(grp.scpshk,b)); pre_vs_cross.sccross_corrp=P(1,2); pre_vs_cross.sccross_corrR=R(1,2);
[R,P]= corrcoef(grp.barbeh,analysis_results(grp.bar,b)); pre_vs_cross.barcross_corrp=P(1,2); pre_vs_cross.barcross_corrR=R(1,2);

% [R,P]= corrcoef(grp.shockbeh,analysis_results(grp.shock,b)-analysis_results(grp.shock,a)); pre_vs_cross.dfdiff_corrp=P(1,2); pre_vs_cross.dfdiff_corrR=R(1,2);
% [R,P]= corrcoef(grp.scpshkbeh,analysis_results(grp.scpshk,b)-analysis_results(grp.scpshk,a)); pre_vs_cross.scdiff_corrp=P(1,2); pre_vs_cross.scdiff_corrR=R(1,2);
% [R,P]= corrcoef(grp.barbeh,analysis_results(grp.bar,b)-analysis_results(grp.bar,a)); pre_vs_cross.bardiff_corrp=P(1,2); pre_vs_cross.bardiff_corrR=R(1,2);


     between_factors=[ones(1,length(grp.shock)) 2*ones(1,length(grp.scpshk)) 3*ones(1,length(grp.bar))]';
%     
     datamat=[analysis_results(grp.shock,b); analysis_results(grp.scpshk,b); analysis_results(grp.bar,b)];% ...
%              [analysis_results(grp.shock,b); analysis_results(grp.scpshk,b); analysis_results(grp.bar,b)]];
             
    [pre_vs_cross.p, pre_vs_cross.tbl] = anova1(datamat, between_factors, 'off');

%     [h, pre_vs_cross.df_pairtp, ci, pre_vs_cross.df_pairt_stats]=ttest(analysis_results(grp.shock,a),analysis_results(grp.shock,b));
%     [h, pre_vs_cross.sc_pairtp, ci, pre_vs_cross.sc_pairt_stats]=ttest(analysis_results(grp.scpshk,a),analysis_results(grp.scpshk,b));
%     [h, pre_vs_cross.bar_pairtp, ci, pre_vs_cross.bar_pairt_stats]=ttest(analysis_results(grp.bar,a),analysis_results(grp.bar,b));
% 
     [h, pre_vs_cross.dfsc_unpair_cross, ci, pre_vs_cross.dfsc_unpair_crossstats]=ttest2(analysis_results(grp.shock,b),analysis_results(grp.scpshk,b));
     [h, pre_vs_cross.dfbar_upair_cross, ci, pre_vs_cross.dfbar_unpair_crossstats]=ttest2(analysis_results(grp.shock,b),analysis_results(grp.bar,b));
     [h, pre_vs_cross.scbar_unpair_cross, ci, pre_vs_cross.scbar_unpair_crossstats]=ttest2(analysis_results(grp.scpshk,b),analysis_results(grp.bar,b));
% 
%     [h, pre_vs_cross.df_prev0, ci, pre_vs_cross.df_prev0_stats]=ttest(analysis_results(grp.shock,a));
%     [h, pre_vs_cross.sc_prev0, ci, pre_vs_cross.sc_prev0_stats]=ttest(analysis_results(grp.scpshk,a));
%     [h, pre_vs_cross.bar_prev0, ci, pre_vs_cross.bar_prev0_stats]=ttest(analysis_results(grp.bar,a));
% 
%     [h, pre_vs_cross.df_crossv0, ci, pre_vs_cross.df_crossv0_stats]=ttest(analysis_results(grp.shock,b));
%     [h, pre_vs_cross.sc_crossv0, ci, pre_vs_cross.sc_crossv0_stats]=ttest(analysis_results(grp.scpshk,b));
%     [h, pre_vs_cross.bar_crossv0, ci, pre_vs_cross.bar_crossv0_stats]=ttest(analysis_results(grp.bar,b));
