clear;

rats_to_analyze=[1:4 6:27];
shockflag=0;

placemap_shock=[];
placemap_barrier=[];
placemap_scopshock=[];

load 'Behav_effects_old';

resubsample=false;
plotbeh=false;

maxhighbins=12;
rundir=0;
goodpeakbins = 1:23; %bins of allowed place peaks
minplacescore = 2; %must beat this to count as a place cell (place score is -log10 of median p value from 200 split correlations)
minspikes = 1; %must beat this to count as a place cell (mean spikes fired per beeline )
L_zone=2:8; M_zone=9:15; R_zone=16:23; M_wid=3;

% session 15 of Hipp35

smear=1;
%convkernel=[zeros(1,smear-1) ones(1,smear)]/smear;
convkernel=ones(1,smear)/smear;

analyze_shock = true;

analysis_shift_dfpre=[];
analysis_shift_scpre=[];
analysis_shift_barpre=[];
analysis_shift_dfpost=[];
analysis_shift_scpost=[];
analysis_shift_barpost=[];

dfpre_placecellmap=[];
dfpost_placecellmap=[];
scpre_placecellmap=[];
scpost_placecellmap=[];
barpre_placecellmap=[];
barpost_placecellmap=[];

finaldir='D:\MiniScopeData\Hipp_LFOV\'; %folder where the goodies are stored

rootdir='D:\MiniScopeData\pvAnalysis\'; %folder where the goodies are stored
lineardir='D:\MiniScopeData\Shock_LFOV\'; %folder where the goodies are stored

shockdir='D:\MiniScopeData\Blair_et_al\shocktimes\'; %folder where the goodies are stored
cellmapdir='D:\MiniScopeData\Blair_et_al\cellmaps\'; %folder where the goodies are stored
sessiondir='D:\MiniScopeData\Blair_et_al\sessiondata\';
datadir='D:\MiniScopeData\Blair_et_al\prepost\';

pethwid = 22; %total width of PETH window is 2*pethwid+1

spbinsize=13.33;                  % width of each spatial bin in cm
bincenters=-(spbinsize*11):spbinsize:(spbinsize*11);
binedges=[-500 [-(10*spbinsize):spbinsize:(spbinsize*11)]-spbinsize/2 500];

cellmat_shockcols = [...
    4 5 8;...  %r=1  Hipp6 ext ---- sep3
    4 6 11;... %r=2 Hipp7 ---- sep4
    1 3 7;...  %r=3 Hipp8 ---- sep4
    2 4 8;...  %r=4 Hipp9 ---- sep4
    6 8 12;... %r=5 Hipp15
    3 5 9;...  %r=6 Hipp31 ---- sep4
    3 5 9;... %r=7 Hipp12 (barrier)
    3 5 9;... %r=8 Hipp13 (barrier)
    6 6 9;... %r=9 Hipp18 (barrier) ext
    5 7 11;... %r=10 Hipp35 (barrier)
    5 5 9;... %r=11 Hipp30 (barrier) ext
    5 7 11;... %r=12 Hipp34 (barrier)
    
    7 9 13;... %r=13 Hipp12 (shock) ---- sep4
    7 9 13;... %r=14 Hipp13 (shock) ---- sep4
    12 13 16;... %r=15 Hipp18 (shock) ext ---- sep2
    9 11 15;...   %r=16 Hipp35 (shock) ---- sep4
    
    7 9 13;... %r=17 Hipp30 (shock+scop)
    10 10 14;... %r=18 Hipp32 (shock+scop)
    11 11 15;... %r=19 Hipp34 (shock+scop)
    9 11 15;...   %r=20 Hipp36 (shock+scop)
    
    20 22 26;... %r=21 Hipp12 (shock+scop)
    20 22 26;... %r=22 Hipp13 (shock+scop)
    20 22 26;... %r=23 Hipp15 (shock+scop)
    20 22 26;...   %r=24 Hipp18 (shock+scop)
    19 19 23;...   %r=25 Hipp31 (shock+scop)
    
    13 13 17;... %r=26 Hipp30 (drug free shock2) ---- sep2
    16 19 21;... %r=27 Hipp32 (drug free shock2) ---- sep5
    
    14 14 16;... %Hipp31 (scopolamine alone)
    5 5 7;... %Hipp32 (scopolamine alone)
    6 6 8]; %Hipp36 (scopolamine alone)

%shocked first
rat(1).name='Hipp6';
rat(2).name='Hipp7';
rat(3).name='Hipp8';
rat(4).name='Hipp9';
rat(5).name='Hipp15'; %F:\ScaryData2\Hipp12to18\small_contours\Hipp15\submatching\shock added 6 7 8 9 10 12 13
rat(6).name='Hipp31';
%barrier first
rat(7).name='Hipp12';
rat(8).name='Hipp13';
rat(9).name='Hipp18';
rat(10).name='Hipp35';
rat(11).name='Hipp30';
rat(12).name='Hipp34';
%shock after barrier (nor prior shock)
rat(13).name='Hipp12';
rat(14).name='Hipp13';
rat(15).name='Hipp18';
rat(16).name='Hipp35';
%shock on scopolamine (no prior shock, all use -6,0,+6 spacing)
rat(17).name='Hipp30';
rat(18).name='Hipp32';
rat(19).name='Hipp34';
rat(20).name='Hipp36';
%shock on scopolamine (prior shock, all use -6,0,+6 spacing)
rat(21).name='Hipp12';
rat(22).name='Hipp13';
rat(23).name='Hipp15';
rat(24).name='Hipp18';
rat(25).name='Hipp31';
%drug free shock 2
rat(26).name='Hipp30';
rat(27).name='Hipp32';
%scopolamine alone
rat(28).name='Hipp31';
rat(29).name='Hipp32';
rat(30).name='Hipp36';

rat(1).shockrow=1;
rat(2).shockrow=2;
rat(3).shockrow=3;
rat(4).shockrow=4;
rat(5).shockrow=7; %F:\ScaryData2\Hipp12to18\small_contours\Hipp15\submatching\shock added 6 7 8 9 10 12 13
rat(6).shockrow=10;

rat(7).shockrow=5;
rat(8).shockrow=6;
rat(9).shockrow=8;
rat(10).shockrow=13;
rat(11).shockrow=9;
rat(12).shockrow=12;

rat(13).shockrow=5;
rat(14).shockrow=6;
rat(15).shockrow=8;
rat(16).shockrow=13;

rat(17).shockrow=9;
rat(18).shockrow=11;
rat(19).shockrow=12;
rat(20).shockrow=14;

rat(21).shockrow=5;
rat(22).shockrow=6;
rat(23).shockrow=7;
rat(24).shockrow=8;
rat(25).shockrow=10;

rat(26).shockrow=1;
rat(27).shockrow=2;

clear heatmaps* all_chmat Loss* Decode* TCcorr*;

heatmaps_shocksess = []; heatmaps_shockflags=[]; heatmaps_shockrows=[];
heatmaps_presess = []; heatmaps_preflags=[]; heatmaps_prerows=[];
heatmaps_postsess = []; heatmaps_postflags=[]; heatmaps_postrows=[];

sessions_to_analyze = [2 1 3]; %dont use

hm_pre_shk=[];
hm_shockpre_shk=[];
hm_shockpost_shk=[];
hm_post_shk=[];
hm_pflags_pre_shk=[];
hm_pflags_post_shk=[];

hm_pflags_pre_shk=[];
hm_pflags_post_shk=[];

hm_post_shk_safeR=[];
hm_post_bar_safeR=[];
hm_post_scp_safeR=[];

hm_post_shk_group=[];

hm_post_shk_safeL=[];
hm_post_bar_safeL=[];
hm_post_scp_safeL=[];

hm_post_shk_centerz=[];
hm_post_bar_centerz=[];
hm_post_scp_centerz=[];

hm_post_bar=[];
hm_pflags_post_bar=[];

hm_shockpost_bar=[];
hm_shockpost_bar_group=[];

hm_post_bar_group=[];
hm_post_scpshk_group=[];

hm_shockpost_scpshk=[];
hm_post_scpshk=[];
hm_pflags_post_scpshk=[];

hm_post_scp=[];
hm_pflags_post_scp=[];

analysis_results=[];

for r=rats_to_analyze %rats_to_analyze %--- loop through rats
    
    r
    %need to keep some variables when new rat analysis starts
    clearvars -except hm_* plotbeh resubsample placemap_* *_placecellmap event_place_cells* shockflag global_p Effects rundir *_zone M_wid maxhighbins goodpeakbins minplacescore minspikes *_segment* Rmatrix* Pmatrix* barriertimes allshockcells* allnonshockcells* SHKexc recur_shock* predata postdata shockdata precol postcol analyze_* data otherdata smear convkernel Decode* TCcorr* Loss* chi* MdlSpec rat* session* analysis* r goodsess s *dir *_result *shocktimes* cellmat* heatmaps* pethwid all_chmat binedges;
    
    clear sessionNums;
    
    %read in matching matrix for this rat
    if r<=6 | (r>=13 & r<=16) %drug-free shocks
        
        load([cellmapdir rat(r).name '_shock_cmap']);
        
    elseif r>=17 & r<=25 %scopolamine shocks
        
        load([cellmapdir rat(r).name '_scopshock_cmap']);
        
    elseif (r>=7 & r<=12) %barrier
        
        load([cellmapdir rat(r).name '_barrier_cmap']);
        
    elseif (r<=27) %second drug free shock
        
        load([cellmapdir rat(r).name '_shock2_cmap']);
        
    else %second drug free shock
        
        load([cellmapdir rat(r).name '_scopo_cmap']);
        
    end
    
    
    %cmap column indices for the sessions of interest
    pre_column=find(sessionNums==cellmat_shockcols(r,2));
    post_column=find(sessionNums==cellmat_shockcols(r,3));

    
    analysis_sessint(r)=cellmat_shockcols(r,3)-cellmat_shockcols(r,2);
    
    precol=2; %pre session
    load([sessiondir rat(r).name '_linear' num2str(cellmat_shockcols(r,precol)) '_sess']);
    eval(['preframe = frame' num2str(cellmat_shockcols(r,precol)) ';']);
    load([datadir rat(r).name '_linear' num2str(cellmat_shockcols(r,precol)) '_predata']);
    
    %load shock session data
    postcol=3;
    load([sessiondir rat(r).name '_linear' num2str(cellmat_shockcols(r,postcol)) '_sess']);
    eval(['postframe = frame' num2str(cellmat_shockcols(r,postcol)) ';']);
    load([datadir rat(r).name '_linear' num2str(cellmat_shockcols(r,postcol)) '_postdata']);
    
    analysis_beelines(r,:)=[0 0 size(predata.LRint,1) size(predata.RLint,1) size(postdata.LRint,1) size(postdata.RLint,1)];
    
    analysis_beelines(:,1)=analysis_beelines(:,3)+analysis_beelines(:,4);
    analysis_beelines(:,2)=analysis_beelines(:,5)+analysis_beelines(:,6);
        
    [pre_peak.valLR, pre_peak.LRd]=max(squeeze(predata.dcurve_LR)'); [pre_peak.valRL, pre_peak.RLd]=max(squeeze(predata.dcurve_RL)');
    [post_peak.valLR, post_peak.LRd]=max(squeeze(postdata.dcurve_LR)'); [post_peak.valRL, post_peak.RLd]=max(squeeze(postdata.dcurve_RL)');
    
    clear premaxrateLR premaxrateRL prehighbinsLR prehighbinsRL;
    for i=1:size(preframe.S,1)
        premaxrateLR(i)=max(predata.dcurve_LR(i,:)*1000);
        premaxrateRL(i)=max(predata.dcurve_RL(i,:)*1000);
        prehighbinsLR(i)=sum(predata.dcurve_LR(i,:)>max(predata.dcurve_LR(i,:))*.7);
        prehighbinsRL(i)=sum(predata.dcurve_RL(i,:)>max(predata.dcurve_RL(i,:))*.7);
        predata.infoLR(i) = predata.dcurve_LR_bps(i);
        predata.infoRL(i) = predata.dcurve_RL_bps(i);
    end
    clear postmaxrateLR postmaxrateRL posthighbinsLR posthighbinsRL;
    for i=1:size(postframe.deconv,1)
        postmaxrateLR(i)=max(postdata.dcurve_LR(i,:)*1000);
        postmaxrateRL(i)=max(postdata.dcurve_RL(i,:)*1000);
        posthighbinsLR(i)=sum(postdata.dcurve_LR(i,:)>max(postdata.dcurve_LR(i,:))*.7);
        posthighbinsRL(i)=sum(postdata.dcurve_RL(i,:)>max(postdata.dcurve_RL(i,:))*.7);
        postdata.infoLR(i) = postdata.dcurve_LR_bps(i);
        postdata.infoRL(i) = postdata.dcurve_RL_bps(i);
    end
    
    
    % only the peak rate thresh is recalculated per 20 iteration shuffle; place score and
    % spikes per beeline come from the 200 iteration place field classifier
    preflags.placeLR=prehighbinsLR(:)<maxhighbins & predata.place_scoreLR(:)>minplacescore & predata.spikes_per_LRbee(:)>minspikes;
    preflags.placeRL=prehighbinsRL(:)<maxhighbins & predata.place_scoreRL(:)>minplacescore & predata.spikes_per_RLbee(:)>minspikes;
    postflags.placeLR=posthighbinsLR(:)<maxhighbins & postdata.place_scoreLR(:)>minplacescore & postdata.spikes_per_LRbee(:)>minspikes;
    postflags.placeRL=posthighbinsRL(:)<maxhighbins & postdata.place_scoreRL(:)>minplacescore & postdata.spikes_per_RLbee(:)>minspikes;
    
    predata.isplaceLR_N=preflags.placeLR;
    predata.isplaceRL_N=preflags.placeRL;
    predata.isplace_N=predata.isplaceLR_N | predata.isplaceRL_N;
    postdata.isplaceLR_N=postflags.placeLR;
    postdata.isplaceRL_N=postflags.placeRL;
    postdata.isplace_N=postdata.isplaceLR_N | postdata.isplaceRL_N;
    
    predata.isrecurringplace=false(1,length(predata.isplace_N));
    postdata.isrecurringplace=false(1,length(postdata.isplace_N));
    
    
    analysis_meaninfo(r,:) = [0 nanmean([predata.infoLR(predata.isplaceLR_N(:)') predata.infoRL(predata.isplaceRL_N(:)')]) nanmean([postdata.infoLR(postdata.isplaceLR_N(:)') postdata.infoRL(postdata.isplaceRL_N(:)')])];
    analysis_meanspeed(r,:) = [0 nanmean(predata.subspd_median) nanmean(postdata.subspd_median)];
    analysis_meanfsize(r,:) = 12*[0 nanmean([prehighbinsLR(predata.isplaceLR_N(:)') prehighbinsRL(predata.isplaceRL_N(:)')]) nanmean([posthighbinsLR(postdata.isplaceLR_N(:)') posthighbinsRL(postdata.isplaceRL_N(:)')])];
    analysis_meanmaxrate(r,:) = [0 nanmean([premaxrateLR(predata.isplaceLR_N(:)') premaxrateRL(predata.isplaceRL_N(:)')]) nanmean([postmaxrateLR(postdata.isplaceLR_N(:)') postmaxrateRL(postdata.isplaceRL_N(:)')])];
    analysis_Nplacecells(r,:) = [0 sum(predata.isplace_N) sum(postdata.isplace_N)];
    analysis_meanplaceperc(r,:) = [0 mean(sum(predata.isplace_N))/size(predata.dcurve_LR,1) mean(sum(postdata.isplace_N))/size(postdata.dcurve_LR,1)];
        
    clear *matrow *matrowLR *matrowRL;
    preplacematrowLR=[];
    preplacematrowRL=[];
    shkplacematrowLR=[];
    shkplacematrowRL=[];
    postplacematrowLR=[];
    postplacematrowRL=[];
    predexLR=find(predata.isplaceLR_N);
    predexRL=find(predata.isplaceRL_N);
    for i=predexLR(:)'
        if ~isempty(find(cmap(:,pre_column)==i))
            preplacematrowLR(i)=find(cmap(:,pre_column)==i);
        end
    end
    for i=predexRL(:)'
        if ~isempty(find(cmap(:,pre_column)==i))
            preplacematrowRL(i)=find(cmap(:,pre_column)==i);
        end
    end
    postdexLR=find(postdata.isplaceLR_N);
    postdexRL=find(postdata.isplaceRL_N);
    for i=postdexLR(:)'
        if ~isempty(find(cmap(:,post_column)==i))
            postplacematrowLR(i)=find(cmap(:,post_column)==i);
        end
    end
    for i=postdexRL(:)'
        if ~isempty(find(cmap(:,post_column)==i))
            postplacematrowRL(i)=find(cmap(:,post_column)==i);
        end
    end
    
    analysis_field_Nplace(r,:)=[0 sum(predata.isplace_N) sum(postdata.isplace_N)];
    analysis_field_Nnonplace(r,:)=[0 sum(~predata.isplace_N) sum(~postdata.isplace_N)];
    analysis_field_Nonedir(r,:)=[0 sum( (~predata.isplaceLR_N & predata.isplaceRL_N) | (predata.isplaceLR_N & ~predata.isplaceRL_N) )  sum( (~postdata.isplaceLR_N & postdata.isplaceRL_N) | (postdata.isplaceLR_N & ~postdata.isplaceRL_N) )];
    analysis_field_Ntwodir(r,:)=[0 sum( predata.isplaceLR_N & predata.isplaceRL_N  )  sum( postdata.isplaceLR_N & postdata.isplaceRL_N )];
    
    preplacematrowLR=preplacematrowLR(preplacematrowLR>0);
    preplacematrowRL=preplacematrowRL(preplacematrowRL>0);
    postplacematrowLR=postplacematrowLR(postplacematrowLR>0);
    postplacematrowRL=postplacematrowRL(postplacematrowRL>0);
    
    recur_shockpost_matdex=find(cmap(:,post_column)>0 & cmap(:,pre_column)>0); %rows of cell2index map that contain cells recurring in shock & pre sessions
    recur_shockpost_postflag=false(1,size(postframe.deconv,1)); %initialize recurrence flags
    recur_shockpost_postflag(cmap(recur_shockpost_matdex,post_column))=true; %logical vector of cells recorded in SHOCK session, flagged if they recur in the other session
    recur_shockpost_shockflag=false(1,size(preframe.deconv,1)); %initialize recurrence flags
    recur_shockpost_shockflag(cmap(recur_shockpost_matdex,pre_column))=true; %logical vector of cells recorded in SHOCK session, flagged if they recur in the other session
    recur_shockpost_matdex=recur_shockpost_matdex(find(recur_shockpost_matdex>0));
    
    for d=1:length(recur_shockpost_matdex)
        hmap_post(d,:)=[postdata.dcurve_LR(cmap(recur_shockpost_matdex(d),post_column),:) postdata.dcurve_RL(cmap(recur_shockpost_matdex(d),post_column),:)];
        hmap_shockpost(d,:)=[predata.dcurve_LR(cmap(recur_shockpost_matdex(d),pre_column),:) predata.dcurve_RL(cmap(recur_shockpost_matdex(d),pre_column),:)];
        posthighbinsLR=sum(postdata.dcurve_LR(cmap(recur_shockpost_matdex(d),post_column),:)>max(postdata.dcurve_LR(cmap(recur_shockpost_matdex(d),post_column),:))*.7);
        preposthighbinsLR=sum(predata.dcurve_LR(cmap(recur_shockpost_matdex(d),pre_column),:)>max(predata.dcurve_LR(cmap(recur_shockpost_matdex(d),pre_column),:))*.7);
        posthighbinsRL=sum(postdata.dcurve_RL(cmap(recur_shockpost_matdex(d),post_column),:)>max(postdata.dcurve_RL(cmap(recur_shockpost_matdex(d),post_column),:))*.7);
        preposthighbinsRL=sum(predata.dcurve_RL(cmap(recur_shockpost_matdex(d),pre_column),:)>max(predata.dcurve_RL(cmap(recur_shockpost_matdex(d),pre_column),:))*.7);
        isplace_postLR(d) = (posthighbinsLR<maxhighbins & postdata.place_scoreLR(cmap(recur_shockpost_matdex(d),post_column))>minplacescore & postdata.spikes_per_LRbee(cmap(recur_shockpost_matdex(d),post_column))>minspikes )...
            | (preposthighbinsLR<maxhighbins & predata.place_scoreLR(cmap(recur_shockpost_matdex(d),pre_column))>minplacescore & predata.spikes_per_LRbee(cmap(recur_shockpost_matdex(d),pre_column))>minspikes );
        isplace_postRL(d) = (posthighbinsRL<maxhighbins & postdata.place_scoreRL(cmap(recur_shockpost_matdex(d),post_column))>minplacescore & postdata.spikes_per_RLbee(cmap(recur_shockpost_matdex(d),post_column))>minspikes )...
            | (preposthighbinsRL<maxhighbins & predata.place_scoreRL(cmap(recur_shockpost_matdex(d),pre_column))>minplacescore & predata.spikes_per_RLbee(cmap(recur_shockpost_matdex(d),pre_column))>minspikes );
        %                    isplace_preRL(d) = predata.isplaceRL(cmap(recur_shockpre_matdex(d),pre_column)) | data.isplaceRL(cmap(recur_shockpre_matdex(d),pre_column));
        isplace_post(d) = (isplace_postLR(d) | isplace_postRL(d));
        nmax=max([hmap_post(d,:) hmap_shockpost(d,:)]);
        if isplace_post(d)
            postdata.isrecurringplace(cmap(recur_shockpost_matdex(d),post_column))=true;
            predata.isrecurringplace(cmap(recur_shockpost_matdex(d),pre_column))=true;
        end
        
        [temp, peakbin_postLR(d)]=max(postdata.dcurve_LR(d,:));
        [temp, peakbin_postRL(d)]=max(postdata.dcurve_RL(d,:));
        [temp, peakbin_shockpostLR(d)]=max(predata.dcurve_LR(d,:));
        [temp, peakbin_shockpostRL(d)]=max(predata.dcurve_RL(d,:));
        
    end
    
    analysis_onlyboth_recur(r,:)=[sum(isplace_postLR(:) & isplace_postRL(:)) sum(isplace_postLR(:) & ~isplace_postRL(:)) sum(~isplace_postLR(:) & isplace_postRL(:))];
    
    analysis_field_Nrecur(r,:)=sum(predata.isrecurringplace(:));
    
    analysis_field_Nnonrecur(r,:)=sum(postdata.isplace_N(:) & ~postdata.isrecurringplace(:)) + sum(predata.isplace_N(:) & ~predata.isrecurringplace(:));
    
    
    if (r<7 & ~(r==5)) | (r>12 & r<17) | (r==26) | (r==27)
        hm_shockpost_shk=[hm_shockpost_shk; hmap_shockpost];
        hm_post_shk_group=[hm_post_shk_group; hmap_post(:,12)*0+r];
        hm_post_shk=[hm_post_shk; hmap_post];
        %               hm_pflags_pre_shk=[hm_pflags_pre_shk; [isplace_preLR(:) isplace_preRL(:)] ];
        hm_pflags_post_shk=[hm_pflags_post_shk; [isplace_postLR(:) isplace_postRL(:)] ];
    elseif r>6 & r<13
        hm_shockpost_bar=[hm_shockpost_bar; hmap_shockpost];
        hm_post_bar_group=[hm_post_bar_group; hmap_post(:,12)*0+r];
        hm_post_bar=[hm_post_bar; hmap_post];
        %       hm_pflags_pre_bar=[hm_pflags_pre_bar; [isplace_preLR(:) isplace_preRL(:)] ];
        hm_pflags_post_bar=[hm_pflags_post_bar; [isplace_postLR(:) isplace_postRL(:)] ];
    elseif r>16 & r<26
        hm_shockpost_scpshk=[hm_shockpost_scpshk; hmap_shockpost];
        hm_post_scpshk_group=[hm_post_scpshk_group; hmap_post(:,12)*0+r];
        hm_post_scpshk=[hm_post_scpshk; hmap_post];
        %      hm_pflags_pre_scpshk=[hm_pflags_pre_scpshk; [isplace_preLR(:) isplace_preRL(:)] ];
        hm_pflags_post_scpshk=[hm_pflags_post_scpshk; [isplace_postLR(:) isplace_postRL(:)] ];
    end
    
    [pmat_postR, pmat_postP] = pv_heatmap_cormatrix([hmap_post(isplace_postLR,1:23); hmap_post(isplace_postRL,24:end)],[hmap_shockpost(isplace_postLR,1:23); hmap_shockpost(isplace_postRL,24:end)]);
    
    [pmatLR_postR, pmatLR_postP] = pv_heatmap_cormatrix([hmap_post(isplace_postLR,1:23)],[hmap_shockpost(isplace_postLR,1:23)]);
    
    [pmatRL_postR, pmatRL_postP] = pv_heatmap_cormatrix([hmap_post(isplace_postRL,24:end)],[hmap_shockpost(isplace_postRL,24:end)]);
    
    clear pmat_postR;
    pmat_postR(1,:,:)=pmatLR_postR; pmat_postR(2,:,:)=pmatRL_postR; %pmat_postR=squeeze(nanmean(pmat_postR));
    
    shortpath_postR=[];
    shortpath_postP=[];
    
    middle_postR=[];
    middle_postP=[];
    left_postR=[];
    left_postP=[];
    right_postR=[];
    right_postP=[];
    
    for i=1:23
        
        PC_diag_postR(i)=nanmean(pmat_postR(:,i,i));
        
    end
    
    adjM_zoneLR=M_zone;
    adjM_zoneRL=M_zone;
    
    analysis_adjzone(r,:)=[adjM_zoneLR(1) adjM_zoneRL(1)];
    
    adj_zone=unique([adjM_zoneLR adjM_zoneRL]);
    
    LR_offdiagpost=0;
    RL_offdiagpost=0;
   % N_offdiag=0;
    diagthresh=5;
    for i=1:23
        for j=i:23
            if abs(i-j)>diagthresh
                LR_offdiagpost=[LR_offdiagpost; pmatLR_postR(i,j); pmatLR_postR(j,i)]; 
                RL_offdiagpost=[RL_offdiagpost; pmatRL_postR(i,j); pmatRL_postR(j,i)];
    %            N_offdiag=N_offdiag+2;
            end
        end
    end
    LR_offdiagpost=nanmean(LR_offdiagpost);  
    RL_offdiagpost=nanmean(RL_offdiagpost);     
    
    for i=1:23
        
        if ~(r==25)%if analyzing LR journeys
            if i>=2 & i<22 %not at an end bin
                shortpath_postR=[shortpath_postR pmatLR_postR(i,i) pmatLR_postR(i+1,i) pmatLR_postR(i,i+1)];
            end
            if i==22 %last non end bin
                shortpath_postR=[shortpath_postR pmatLR_postR(i,i)];
            end
            if i>=adjM_zoneLR(1) & i<=adjM_zoneLR(end) %in the roaming middle zone
                middle_postR=[middle_postR pmatLR_postR(i,i) pmatLR_postR(i-1,i) pmatLR_postR(i,i-1)];
                middle_postP=[middle_postP pmatLR_postP(i,i) pmatLR_postP(i-1,i) pmatLR_postP(i,i-1)];
                if i==adjM_zoneLR(end) %last bin of the roaming middle zone
                    middle_postR=[middle_postR pmatLR_postR(i+1,i) pmatLR_postR(i,i+1)];
                    middle_postP=[middle_postP pmatLR_postP(i+1,i) pmatLR_postP(i,i+1)];
                end
            elseif i<=L_zone(end) & i>=2%in the static left zone
                left_postR=[left_postR pmatLR_postR(i,i)];
                left_postP=[left_postP pmatLR_postP(i,i)];
                if i<L_zone(end) %not at last bin of the static left zone
                    left_postR=[left_postR pmatLR_postR(i+1,i) pmatLR_postR(i,i+1)];
                    left_postP=[left_postP pmatLR_postP(i+1,i) pmatLR_postP(i,i+1)];
                end
            elseif i>=R_zone(1) & i<=22 %int the static right zone
                right_postR=[right_postR pmatLR_postR(i,i)];
                right_postP=[right_postP pmatLR_postP(i,i)];
                if i>R_zone(1) %not at first bin of the static right zone
                    right_postR=[right_postR pmatLR_postR(i-1,i) pmatLR_postR(i,i-1)];
                    right_postP=[right_postP pmatLR_postP(i-1,i) pmatLR_postP(i,i-1)];
                end
            end
        end
        
        if i>=2 & i<22 %not at an end bin
            shortpath_postR=[shortpath_postR pmatRL_postR(i,i) pmatRL_postR(i+1,i) pmatRL_postR(i,i+1)];
        end
        if i==22 %last non end bin
            shortpath_postR=[shortpath_postR pmatRL_postR(i,i)];
        end
        if i>=adjM_zoneRL(1) & i<=adjM_zoneRL(end) %in the roaming middle zone
            middle_postR=[middle_postR pmatRL_postR(i,i) pmatRL_postR(i-1,i) pmatRL_postR(i,i-1)];
            middle_postP=[middle_postP pmatRL_postP(i,i) pmatRL_postP(i-1,i) pmatRL_postP(i,i-1)];
            if i==adjM_zoneRL(end) %last bin iof the roaming middle zone
                middle_postR=[middle_postR pmatRL_postR(i+1,i) pmatRL_postR(i,i+1)];
                middle_postP=[middle_postP pmatRL_postP(i+1,i) pmatRL_postP(i,i+1)];
            end
        elseif i<=L_zone(end) & i>=2%in the static left zone
            left_postR=[left_postR pmatRL_postR(i,i)];
            left_postP=[left_postP pmatRL_postP(i,i)];
            if i<L_zone(end) %not at last bin of the static left zone
                left_postR=[left_postR pmatRL_postR(i+1,i) pmatRL_postR(i,i+1)];
                left_postP=[left_postP pmatRL_postP(i+1,i) pmatRL_postP(i,i+1)];
            end
        elseif i>=R_zone(1) & i<=22%int the static right zone
            right_postR=[right_postR pmatRL_postR(i,i)];
            right_postP=[right_postP pmatRL_postP(i,i)];
            if i>R_zone(1) %not at first bin of the static right zone
                right_postR=[right_postR pmatRL_postR(i-1,i) pmatRL_postR(i,i-1)];
                right_postP=[right_postP pmatRL_postP(i-1,i) pmatRL_postP(i,i-1)];
            end
        end
    end
    
    safeL_postR=[left_postR ];
    
    safeR_postR=[right_postR ];
    
    centerz_postR=middle_postR;
    
    if r==25
        safeR_postR=[safeR_postR*NaN safeR_postR];
        
        safeL_postR=[safeL_postR*NaN safeL_postR];
        
        centerz_postR=[centerz_postR*NaN centerz_postR];
    end
    
    safeL_postP=[left_postP ];
    
    safeR_postP=[right_postP ];
    
    centerz_postP=middle_postP;
    
    
    
    clear stackR stackP;
    
    if r<5 %shocked first
        Rmatrix_place_shock1_post((r*2-1):(r*2),:,:)=squeeze(pmat_postR);
    elseif r==6 %shocked first
        Rmatrix_place_shock1_post(9:10,:,:)=squeeze(pmat_postR);
    elseif r>=7 & r<=12 %barrier
        Rmatrix_place_barrier_post(((r-6)*2-1):((r-6)*2),:,:)=squeeze(pmat_postR);
    elseif r<=16 | r>=26%first shock after barrier
        if r<=16
            rf=7;
        else
            rf=16;
        end
        Rmatrix_place_shock1_post(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR);
    elseif r<=25 %scop + first shock
        rf=16;
        Rmatrix_place_scopshock_post(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR);
    else %scop alone
        Rmatrix_place_scop_post(r-24,:,:)=squeeze(pmat_postR);
    end
    
    %-------------------- plot sorted heatmaps
    hmapLR_post=hmap_post(isplace_postLR,1:23);
    hmapLR_shockpost=hmap_shockpost(isplace_postLR,1:23);
    [temp, peak]=max(hmapLR_shockpost');
    [temp, peak2]=max(hmapLR_post');
    
    clear hmap_250 hmap_shock250;
    for iii=1:sum(isplace_postLR)
        hmap_250(iii,:)=interp1(1:23,hmapLR_post(iii,:),linspace(1,23,250));
        hmap_shock250(iii,:)=interp1(1:23,hmapLR_shockpost(iii,:),linspace(1,23,250));
    end
    [temp, peak_250]=max(hmap_shock250');
    [temp, peak2_250]=max(hmap_250');
    shiftdistLR_post=abs(peak_250-peak2_250);
    
    peak=peak2;%(~moreplaceyLR_shockpost(isplace_postLR))=peak2(~moreplaceyLR_shockpost(isplace_postLR))+100;
    hmapLR_post=sortrows([peak(:) hmapLR_post],1);
    hmapLR_shockpost=sortrows([peak(:) hmapLR_shockpost],1);
%     subplot_tight(4,7,[1 8 ]+3); imagesc(1000*hmapLR_post(:,2:end)); title(['post ' num2str(cellmat_shockcols(r,postcol))]); set(gca,'XTick',[]); caxis([0 10]);;
%     subplot_tight(4,7,[1 8 ]+2); imagesc(1000*hmapLR_shockpost(:,2:end)); title(['shkpost ' num2str(cellmat_shockcols(r,precol))]); set(gca,'XTick',[]); caxis([0 10]);;
    
    hmapRL_post=hmap_post(isplace_postRL,24:end);
    hmapRL_shockpost=hmap_shockpost(isplace_postRL,24:end);
    [temp, peak]=max(hmapRL_shockpost');
    [temp, peak2]=max(hmapRL_post');
    
    clear hmap_250 hmap_shock250;
    for iii=1:sum(isplace_postRL)
        hmap_250(iii,:)=interp1(1:23,hmapRL_post(iii,:),linspace(1,23,250));
        hmap_shock250(iii,:)=interp1(1:23,hmapRL_shockpost(iii,:),linspace(1,23,250));
    end
    [temp, peak_250]=max(hmap_shock250');
    [temp, peak2_250]=max(hmap_250');
    shiftdistRL_post=abs(peak_250-peak2_250);
    
    peak=peak2;%(~moreplaceyRL_shockpost(isplace_postRL))=peak2(~moreplaceyRL_shockpost(isplace_postRL))+100;
    hmapRL_post=sortrows([peak(:) hmapRL_post],1);
    hmapRL_shockpost=sortrows([peak(:) hmapRL_shockpost],1);
%     subplot_tight(4,7,[1 8 ]+3+14); imagesc(1000*hmapRL_post(:,2:end)); title(['post ' num2str(cellmat_shockcols(r,postcol))]); set(gca,'XTick',[]); caxis([0 10]);;
%     subplot_tight(4,7,[1 8 ]+2+14); imagesc(1000*hmapRL_shockpost(:,2:end)); title(['shkpost ' num2str(cellmat_shockcols(r,precol))]); set(gca,'XTick',[]); caxis([0 10]);;
    
    [post_v_shockLR_R, post_v_shockLR_P] = pv_heatmap_cormatrix(hmapLR_post(:,2:end)',hmapLR_shockpost(:,2:end)');
    [post_v_shockRL_R, post_v_shockRL_P] = pv_heatmap_cormatrix(hmapRL_post(:,2:end)',hmapRL_shockpost(:,2:end)');
    
    prefs=nanmean([sum([predata.dcurve_LR(predata.isplaceLR_N,:)>repmat(max(predata.dcurve_LR(predata.isplaceLR_N,:)')',1,23)*.7]') sum([predata.dcurve_RL(predata.isplaceRL_N,:)>repmat(max(predata.dcurve_RL(predata.isplaceRL_N,:)')',1,23)*.7]')]*12);
    postfs=nanmean([sum([postdata.dcurve_LR(postdata.isplaceLR_N,:)>repmat(max(postdata.dcurve_LR(postdata.isplaceLR_N,:)')',1,23)*.7]') sum([postdata.dcurve_RL(postdata.isplaceRL_N,:)>repmat(max(postdata.dcurve_RL(postdata.isplaceRL_N,:)')',1,23)*.7]')]*12);
    
    prespd=250/nanmedian([predata.LRint(:,2)-predata.LRint(:,1); predata.RLint(:,2)-predata.RLint(:,1)]/11.4);
    postspd=250/nanmedian([postdata.LRint(:,2)-postdata.LRint(:,1); postdata.RLint(:,2)-postdata.RLint(:,1)]/11.4);
    
    preocc=median([preframe.vmap_LR preframe.vmap_RL])/1000;
    postocc=median([postframe.vmap_LR postframe.vmap_RL])/1000;
    
    predata.isplace=predata.isplaceLR_N | predata.isplaceRL_N;
    postdata.isplace=postdata.isplaceLR_N | postdata.isplaceRL_N;
    
    
    if ismember(r,[1:4 6 13:16 26:27])
        
        analysis_shift_dfpost=[analysis_shift_dfpost; [shiftdistLR_post'; shiftdistRL_post']];
        
        dfpre_placecellmap=[dfpre_placecellmap; [predata.dcurve_LR(predata.isplace,:) predata.dcurve_RL(predata.isplace,:)]];
        dfpost_placecellmap=[dfpost_placecellmap; [postdata.dcurve_LR(postdata.isplace,:) postdata.dcurve_RL(postdata.isplace,:)]];
        
        hm_post_shk_safeR=[hm_post_shk_safeR; safeR_postR];
        
        hm_post_shk_safeL=[hm_post_shk_safeL; safeL_postR];
        
        hm_post_shk_centerz=[hm_post_shk_centerz; centerz_postR];
        
    end
    
    if ismember(r,[17:25])
        analysis_shift_scpost=[analysis_shift_scpost; [shiftdistLR_post'; shiftdistRL_post']];
        
        scpre_placecellmap=[scpre_placecellmap; [predata.dcurve_LR(predata.isplace,:) predata.dcurve_RL(predata.isplace,:)]];
        scpost_placecellmap=[scpost_placecellmap; [postdata.dcurve_LR(postdata.isplace,:) postdata.dcurve_RL(postdata.isplace,:)]];
        
        hm_post_scp_safeR=[hm_post_scp_safeR; safeR_postR];
        
        hm_post_scp_safeL=[hm_post_scp_safeL; safeL_postR];
        
        hm_post_scp_centerz=[hm_post_scp_centerz; centerz_postR];
    end
    
    if ismember(r,[7:12])
        analysis_shift_barpost=[analysis_shift_barpost; [shiftdistLR_post'; shiftdistRL_post']];
        
        barpre_placecellmap=[barpre_placecellmap; [predata.dcurve_LR(predata.isplace,:) predata.dcurve_RL(predata.isplace,:)]];
        barpost_placecellmap=[barpost_placecellmap; [postdata.dcurve_LR(postdata.isplace,:) postdata.dcurve_RL(postdata.isplace,:)]];
        
        hm_post_bar_safeR=[hm_post_bar_safeR; safeR_postR];
        
        hm_post_bar_safeL=[hm_post_bar_safeL; safeL_postR];
        
        hm_post_bar_centerz=[hm_post_bar_centerz; centerz_postR];
    end
    
    analysis_results(r,:)=[r 0 sum(predata.isplaceLR_N | predata.isplaceRL_N)/length(predata.isplaceLR_N) sum(postdata.isplaceLR_N | postdata.isplaceRL_N)/length(postdata.isplaceLR_N) ... % 1 - 4
        0 sum(predata.isplaceLR_N | predata.isplaceRL_N)  sum(postdata.isplaceLR_N | postdata.isplaceRL_N) ... %5 6 7
        0 nanmedian([predata.dcurve_LR_bps(predata.isplaceLR_N) predata.dcurve_RL_bps(predata.isplaceRL_N)]) nanmedian([postdata.dcurve_LR_bps(postdata.isplaceLR_N) postdata.dcurve_RL_bps(postdata.isplaceRL_N)])... %8 9 10
        0 ... % 11
        0 length(predata.LRint)+length(predata.RLint) length(postdata.LRint)+length(postdata.RLint) ...  % 12 13 14
        0 nanmedian(max([predata.dcurve_LR(predata.isplace,:) predata.dcurve_RL(predata.isplace,:)]')) nanmedian(max([postdata.dcurve_LR(postdata.isplace,:) postdata.dcurve_RL(postdata.isplace,:)]')) ... %15 16 17
        0  nanmedian(safeL_postR)  0 ... %18 19 20
        0  0  nanmedian(centerz_postR) ... %21 22 23
        0 0  ... %24 25
        0 prefs postfs ... %26 27 28
        0 prespd postspd ...  %29 30 31
        0 nanmedian(safeR_postR) ...
        0 nanmedian([shiftdistLR_post'; shiftdistRL_post'])]; %32 33
    
    [h, p_post_sp] = ttest(safeR_postR);
    
    [h, p_post_safeL] = ttest(safeL_postR);
    
    [h, p_post_centerz] = ttest(centerz_postR);
    
    analysis_sigR(r,:) = [0 0 0 p_post_sp 0 p_post_safeL 0 p_post_centerz];
    
    clear temp*;
    for d=1:length(recur_shockpost_matdex)
        temp1(d)=0;%nanmean(data.dcurve_LR(cmap(recur_shockpre_matdex(d),shk_column),adjM_zoneLR))-nanmean(predata.dcurve_LR(cmap(recur_shockpre_matdex(d),pre_column),adjM_zoneLR));
        temp2(d)=0;%nanmean(data.dcurve_RL(cmap(recur_shockpre_matdex(d),shk_column),adjM_zoneRL))-nanmean(predata.dcurve_RL(cmap(recur_shockpre_matdex(d),pre_column),adjM_zoneRL));
    end
    
    for d=1:length(recur_shockpost_matdex)
        temp3(d)=nanmean(postdata.dcurve_LR(cmap(recur_shockpost_matdex(d),post_column),adjM_zoneLR))-nanmean(predata.dcurve_LR(cmap(recur_shockpost_matdex(d),pre_column),adjM_zoneLR));
        temp4(d)=nanmean(postdata.dcurve_RL(cmap(recur_shockpost_matdex(d),post_column),adjM_zoneRL))-nanmean(predata.dcurve_RL(cmap(recur_shockpost_matdex(d),pre_column),adjM_zoneRL));
    end
    
    analysis_deltamax_pre_centerz(r)=median([temp1 temp2]);
    analysis_deltamax_post_centerz(r)=median([temp3 temp4]);
    
    analysis_deltamax_pre_centerzLR(r)=median(temp1);
    analysis_deltamax_pre_centerzRL(r)=median(temp2);
    analysis_deltamax_post_centerzLR(r)=median(temp3);
    analysis_deltamax_post_centerzRL(r)=median(temp4);
    
    analysis_deltamean_pre_centerz(r)=mean(abs([temp1 temp2]));
    analysis_deltamean_post_centerz(r)=mean(abs([temp3 temp4]));
    
    analysis_deltamean_pre_centerzLR(r)=mean(abs(temp1));
    analysis_deltamean_pre_centerzRL(r)=mean(abs(temp2));
    analysis_deltamean_post_centerzLR(r)=mean(abs(temp3));
    analysis_deltamean_post_centerzRL(r)=mean(abs(temp4));
    
    
    
    temp=unique([cmap(recur_shockpost_matdex(isplace_post),pre_column)]);
    event_place_cells(r)=length(temp);
    
    recur_shockpost_shockflag=find(recur_shockpost_shockflag);
    recur_shockpost_shockmap=[recur_shockpost_matdex(isplace_post) hmap_shockpost(isplace_post,:)];
    recurplace_shockmap=[ recur_shockpost_shockmap];
    [keepercells, ia, ic]=unique(recurplace_shockmap(:,1));
    keepercells=unique(ia);
    recurplace_shockmap=recurplace_shockmap(keepercells,2:end);
    
    postdata.dcurve_LR_bps=max(postdata.dcurve_LR');
    postdata.dcurve_RL_bps=max(postdata.dcurve_RL');
    
    
end %rat loop -----------------------------------------------------------------------------------------------------------------



