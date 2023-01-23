clear;

rats_to_analyze=[1:4 6:27];

rats_to_analyze=[1:4 6:14 16:17 20:25 26]; %dfdf rats

rats_to_analyze=[1:4 6 13:17 20:27]; %dfdf rats

shockflag=1; %0=ignore shock, 1=select responsive cells, -1=select nonresponsive cells
shockpolarity=0; %-1=select inhibited cells if shockflag=1 (does nothing if shockflag not -1)

placemap_shock=[];
placemap_barrier=[];
placemap_scopshock=[];

%load 'Behav_effects';

resubsample=false;
plotbeh=false;

maxhighbins=12;
rundir=0;
goodpeakbins = 1:23; %bins of allowed place peaks
minplacescore = 2; %must beat this to count as a place cell (place score is -log10 of median p value from 200 split correlations)
minspikes = 1; %must beat this to count as a place cell (mean spikes fired per beeline )
L_zone=1:5; M_zone=9:15; R_zone=19:23; M_wid=3;

% session 15 of Hipp35

smear=1;
%convkernel=[zeros(1,smear-1) ones(1,smear)]/smear;
convkernel=ones(1,smear)/smear;

analyze_shock = true;

allshockcells_approachdf=[];
allnonshockcells_approachdf=[];
allshockcells_approachscop=[];
allnonshockcells_approachscop=[];

allnonshockcells_df=[];
allshockcells_df=[];
allshockinhcells_df=[];

allnonshockcells_scop=[];
allshockcells_scop=[];
allshockinhcells_scop=[];

allnonshockcells_peakdf=[];
allshockcells_peakdf=[];
allshockinhcells_peakdf=[];

allnonshockcells_peakscop=[];
allshockcells_peakscop=[];
allshockinhcells_peakscop=[];

allnonshockcells_scop=[];
allshockcells_scop=[];


analysis_shift_dfpre=[];
analysis_shift_scpre=[];
analysis_shift_barpre=[];
analysis_shift_dfpost=[];
analysis_shift_scpost=[];
analysis_shift_barpost=[];

analysis_placeinfo_shock=[];
analysis_nonplaceinfo_shock=[];

analysis_placeinfo_scopshock=[];
analysis_nonplaceinfo_scopshock=[];

analysis_placeinfo_barrier=[];
analysis_nonplaceinfo_barrier=[];

dfpre_placecellmap=[];
dfshk_placecellmap=[];
dfpost_placecellmap=[];
scpre_placecellmap=[];
scshk_placecellmap=[];
scpost_placecellmap=[];


topdir = 'D:\MiniScopeData\Blair_et_al\'; %top level folder where data files are found
shockdir=[topdir 'shocktimes\']; %folder where the goodies are stored
cellmapdir=[topdir 'cellmaps\']; %folder where the goodies are stored
sessiondir=[topdir 'sessiondata\'];
predir=[topdir 'pretrn\'];
postdir=[topdir 'prepost\'];

%load in times of shocks
load([shockdir 'shocktimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
load([shockdir 'shocktimes2_13']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
load([shockdir 'scopshocktimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
load([shockdir 'barriertimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37

%number of frames before and after shock for PETH analysis window of shock response
pethwid = 22; %total width of PETH window is 2*pethwid+1

%parameters for binning of spatial position data
spbinsize=13.33;                  % width of each spatial bin in cm
bincenters=-(spbinsize*11):spbinsize:(spbinsize*11);
binedges=[-500 [-(10*spbinsize):spbinsize:(spbinsize*11)]-spbinsize/2 500];


%matrix of column numbers to analyze in each rat's matching matrix
%col1 is pre shock session (day -6), col2 is shock session (day 0), col3 is pos shock session (day +6)
% for Hipp15 (row 5), col2 is the sess9 (last session before the immediate shock session)
cellmat_shockcols = [...
    %    3 6 9;...  %Hipp6 3 day
    5 6 8;...  %r=1 Hipp6 ext
    6 7 11;... %r=2 Hipp7
    3 4 7;...  %r=3 Hipp8
    4 5 8;...  %r=4 Hipp9
    8 9 12;... %r=5 Hipp15
    5 6 9;...  %r=6 Hipp31
    5 6 9;... %r=7 Hipp12 (barrier)
    5 6 9;... %r=8 Hipp13 (barrier)
    %     5 7 9;... %Hipp18 (barrier) 2 day
    6 7 9;... %r=9 Hipp18 (barrier) ext
    7 8 11;... %r=10 Hipp35 (barrier)
    %     3 6 9;... %Hipp30 (barrier) 3 dayRmatrix
    5 6 9;... %r=11 Hipp30 (barrier) ext
    7 8 11;... %r=12 Hipp34 (barrier)
    
    9 10 13;... %r=13 Hipp12 (shock)
    9 10 13;... %r=14 Hipp13 (shock)
    %    9 13 17;... %Hipp18 (shock) 4 day
    12 13 16;... %r=15 Hipp18 (shock) ext
    11 12 15;...   %r=16 Hipp35 (shock)
    
    %    9 10 11;... %Hipp30 (shock+scop)
    9 10 13;... %r=17 Hipp30 (shock+scop)
    10 11 14;... %r=18 Hipp32 (shock+scop)
    11 12 15;... %r=19 Hipp34 (shock+scop)
    11 12 15;...   %r=20 Hipp36 (shock+scop)
    
    22 23 26;... %r=21 Hipp12 (shock+scop)
    22 23 26;... %r=22 Hipp13 (shock+scop)
    22 23 26;... %r=23 Hipp15 (shock+scop)
    22 23 26;...   %r=24 Hipp18 (shock+scop)
    %22 23 24;...   %Hipp18 (shock+scop)
    %17 20 23;...   %Hipp31 (shock+scop)
    19 20 23;...   %r=25 Hipp31 (shock+scop)
    
    13 14 17;... %r=26 Hipp30 (drug free shock2)
    19 20 21;... %r=27 Hipp32 (drug free shock2)
    
    14 15 16;... %r=28 Hipp31 (scopolamine alone)
    5 6 7;... %r=29 Hipp32 (scopolamine alone)
    6 7 8]; %r=30 Hipp36 (scopolamine alone)

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

hm_pre_shk_shortpath=[]; hm_post_shk_shortpath=[];
hm_pre_bar_shortpath=[]; hm_post_bar_shortpath=[];
hm_pre_scp_shortpath=[]; hm_post_scp_shortpath=[];

hm_pre_shk_safe=[]; hm_post_shk_safe=[];
hm_pre_bar_safe=[]; hm_post_bar_safe=[];
hm_pre_scp_safe=[]; hm_post_scp_safe=[];

hm_pre_shk_unsafe=[]; hm_post_shk_unsafe=[];
hm_pre_bar_unsafe=[]; hm_post_bar_unsafe=[];
hm_pre_scp_unsafe=[]; hm_post_scp_unsafe=[];

hm_pre_bar=[];
hm_shockpre_bar=[];
hm_shockpost_bar=[];
hm_post_bar=[];
hm_pflags_pre_bar=[];
hm_pflags_post_bar=[];

hm_pre_sh2=[];
hm_shockpre_sh2=[];
hm_shockpost_sh2=[];
hm_post_sh2=[];
hm_pflags_pre_sh2=[];
hm_pflags_post_sh2=[];

hm_pre_scpshk=[];
hm_shockpre_scpshk=[];
hm_shockpost_scpshk=[];
hm_post_scpshk=[];
hm_pflags_pre_scpshk=[];
hm_pflags_post_scpshk=[];

hm_pre_scp=[];
hm_shockpre_scp=[];
hm_shockpost_scp=[];
hm_post_scp=[];
hm_pflags_pre_scp=[];
hm_pflags_post_scp=[];

hm_pre_shk_group=[];
hm_shockpre_shk_group=[];
hm_shockpost_shk_group=[];
hm_post_shk_group=[];

hm_pre_bar_group=[];
hm_shockpre_bar_group=[];
hm_shockpost_bar_group=[];
hm_post_bar_group=[];

hm_pre_scpshk_group=[];
hm_shockpre_scpshk_group=[];
hm_shockpost_scpshk_group=[];
hm_post_scpshk_group=[];

analysis_results=[];

for r=rats_to_analyze %rats_to_analyze %--- loop through rats
    
    r
    
    %need to keep some variables when new rat analysis starts
    clearvars -except postdir cellmapdir hm_* plotbeh resubsample placemap_* *_placecellmap allshockinh* event_place_cells* shockflag global_p Effects rundir *_zone M_wid maxhighbins goodpeakbins minplacescore minspikes *_segment* Rmatrix* Pmatrix* barriertimes allshockcells* allnonshockcells* shockpolarity SHKexc SHKinh recur_shock* predata postdata shockdata precol postcol analyze_* data otherdata smear convkernel Decode* TCcorr* Loss* chi* MdlSpec rat* session* analysis* r goodsess s *dir *_result *shocktimes* cellmat* heatmaps* pethwid all_chmat binedges;
    
    clear sessionNums;
    
    %read in matching matrix for this rat
    if r<=6 | (r>=13 & r<=16) %drug-free shocks
        
        load([cellmapdir rat(r).name '_shock_cmap']);
        if r==5 %for Hipp15 only include cells recorded on pre shock AND shock day
            cmap(find(cmap(:,4+1)==0),4)=0;
        end
        
    elseif r>=17 & r<=25 %scopolamine shocks
        
        load([cellmapdir rat(r).name '_scopshock_cmap']);
        
    elseif (r>=7 & r<=12) %barrier
        
        load([cellmapdir rat(r).name '_barrier_cmap']);
        
    else %second drug free shock
        
        load([cellmapdir rat(r).name '_shock2_cmap']);
        
    end
    
    %cmap column indices for the sessions of interest
    pre_column=find(sessionNums==cellmat_shockcols(r,1));
    shk_column=find(sessionNums==cellmat_shockcols(r,2));
    post_column=find(sessionNums==cellmat_shockcols(r,3));
    
    precol=1; %pre session
    load([sessiondir rat(r).name '_linear' num2str(cellmat_shockcols(r,precol)) '_sess']);
    eval(['preframe = frame' num2str(cellmat_shockcols(r,precol)) ';']);
    load([predir rat(r).name '_linear' num2str(cellmat_shockcols(r,precol)) '_predata']);
    predata.S=preframe.S;
    predata.C=preframe.C;
    predata.x=preframe.x;
    predata.y=preframe.y;
    predata.time=preframe.time;
    predata.deconv=preframe.deconv;
    predata.posbin=preframe.posbin;
    predata.spd=preframe.spd;
    predata.vmap_LR=preframe.vmap_LR;
    predata.vmap_RL=preframe.vmap_RL;
    predata=get_beelines(predata);
    
    postcol=3;
    load([sessiondir rat(r).name '_linear' num2str(cellmat_shockcols(r,postcol)) '_sess']);
    eval(['postframe = frame' num2str(cellmat_shockcols(r,postcol)) ';']);
    load([postdir rat(r).name '_linear' num2str(cellmat_shockcols(r,postcol)) '_postdata']);
    
    oldpostdata=postdata;
    
    postdata.S=postframe.S;
    postdata.C=postframe.C;
    postdata.x=postframe.x;
    postdata.y=postframe.y;
    postdata.time=postframe.time;
    postdata.deconv=postframe.deconv;
    postdata.posbin=postframe.posbin;
    postdata.spd=postframe.spd;
    postdata.vmap_LR=postframe.vmap_LR;
    postdata.vmap_RL=postframe.vmap_RL;
    postdata=get_beelines(postdata);
    
    %load shock session data
    shkcol=2;
    load([sessiondir rat(r).name '_linear' num2str(cellmat_shockcols(r,shkcol)) '_sess']);
    eval(['trnframe = frame' num2str(cellmat_shockcols(r,shkcol)) ';']);
    load([predir rat(r).name '_linear' num2str(cellmat_shockcols(r,shkcol)) '_trndata']);
    data=trndata;
    clear trndata;
    data.S=trnframe.S;
    data.C=trnframe.C;
    data.x=trnframe.x;
    data.y=trnframe.y;
    data.time=trnframe.time;
    data.deconv=trnframe.deconv;
    data.posbin=trnframe.posbin;
    data.spd=trnframe.spd;
    data.vmap_LR=trnframe.vmap_LR;
    data.vmap_RL=trnframe.vmap_RL;
    data=get_beelines(data);
    clear trnframe;
    
    if r>=7 & r<=12
        goodcols=find((data.time/1000)<barriertimes(rat(r).shockrow,1)-.5);
        shockmarks=barriertimes(rat(r).shockrow,:);
    elseif r>=17 & r<=26
        goodcols=find((data.time/1000)<scopshocktimes(rat(r).shockrow,1)-.5);
        shockmarks=scopshocktimes(rat(r).shockrow,:);
    elseif r<5 | (r>=13 & r<=16) | r==6
        goodcols=find((data.time/1000)<shocktimes(rat(r).shockrow,1)-.5);
        shockmarks=shocktimes(rat(r).shockrow,:);
    elseif r>=26
        goodcols=find((data.time/1000)<shocktimes2(rat(r).shockrow,1)-.5);
        shockmarks=shocktimes2(rat(r).shockrow,:);
    else
        goodcols=1:size(data.deconv,2);
        shockmarks=[];
    end
    tempS=data.deconv(:,goodcols)>0;
    fired_on_short=sum([tempS & repmat(data.y(goodcols)'<20,size(tempS,1),1)]');
    
    
    %cells are only included in the count if they fire during the selected beelines
    
    beethresh=0;
    beemask=predata.x*0;
    for i=1:size(predata.LRint,1)
        beemask(predata.LRint(i,1):predata.LRint(i,2))=1;
    end
    for i=1:size(predata.RLint,1)
        beemask(predata.RLint(i,1):predata.RLint(i,2))=1;
    end
    clear goodcell;
    for i=1:size(predata.deconv,1)
        if nansum(predata.deconv(i,find(beemask))>beethresh)
            goodcell(i) = true;
        else
            goodcell(i) = false;
        end
    end
    predata.numgoodcells = sum(goodcell);
    
    beemask=data.x*0;
    for i=1:size(data.LRint,1)
        beemask(data.LRint(i,1):data.LRint(i,2))=1;
    end
    for i=1:size(data.RLint,1)
        beemask(data.RLint(i,1):data.RLint(i,2))=1;
    end
    clear goodcell;
    for i=1:size(data.deconv,1)
        if sum(data.deconv(i,find(beemask))>beethresh)
            goodcell(i) = true;
        else
            goodcell(i) = false;
        end
    end
    data.numgoodcells = sum(goodcell);
    
    beemask=postdata.x*0;
    for i=1:size(postdata.LRint,1)
        beemask(postdata.LRint(i,1):postdata.LRint(i,2))=1;
    end
    for i=1:size(postdata.RLint,1)
        beemask(postdata.RLint(i,1):postdata.RLint(i,2))=1;
    end
    clear goodcell;
    for i=1:size(postdata.deconv,1)
        if sum(postdata.deconv(i,find(beemask))>beethresh)
            goodcell(i) = true;
        else
            goodcell(i) = false;
        end
    end
    postdata.numgoodcells = sum(goodcell);
    
    
    [pre_peak.valLR, pre_peak.LRd]=max(squeeze(predata.dcurve_LR)'); [pre_peak.valRL, pre_peak.RLd]=max(squeeze(predata.dcurve_RL)');
    [shk_peak.valLR, shk_peak.LRd]=max(squeeze(data.dcurve_LR)'); [shk_peak.valRL, shk_peak.RLd]=max(squeeze(data.dcurve_RL)');
    [post_peak.valLR, post_peak.LRd]=max(squeeze(postdata.dcurve_LR)'); [post_peak.valRL, post_peak.RLd]=max(squeeze(postdata.dcurve_RL)');
    
    for i=1:size(predata.deconv,1)
        premaxrateLR(i)=max(predata.dcurve_LR(i,:)*1000);
        premaxrateRL(i)=max(predata.dcurve_RL(i,:)*1000);
        prehighbinsLR(i)=sum(predata.dcurve_LR(i,:)>max(predata.dcurve_LR(i,:))*.7);
        prehighbinsRL(i)=sum(predata.dcurve_RL(i,:)>max(predata.dcurve_RL(i,:))*.7);
        predata.infoLR(i) = predata.dcurve_LR_bps(i);
        predata.infoRL(i) = predata.dcurve_RL_bps(i);
    end
    for i=1:size(data.deconv,1)
        shkmaxrateLR(i)=max(data.dcurve_LR(i,:)*1000);
        shkmaxrateRL(i)=max(data.dcurve_RL(i,:)*1000);
        shkhighbinsLR(i)=sum(data.dcurve_LR(i,:)>max(data.dcurve_LR(i,:))*.7);
        shkhighbinsRL(i)=sum(data.dcurve_RL(i,:)>max(data.dcurve_RL(i,:))*.7);
        data.infoLR(i) = data.dcurve_LR_bps(i);
        data.infoRL(i) = data.dcurve_RL_bps(i);
    end
    for i=1:size(postdata.deconv,1)
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
    shkflags.placeLR=shkhighbinsLR(:)<maxhighbins & data.place_scoreLR(:)>minplacescore & data.spikes_per_LRbee(:)>minspikes;
    shkflags.placeRL=shkhighbinsRL(:)<maxhighbins & data.place_scoreRL(:)>minplacescore & data.spikes_per_RLbee(:)>minspikes;
    postflags.placeLR=posthighbinsLR(:)<maxhighbins & postdata.place_scoreLR(:)>minplacescore & postdata.spikes_per_LRbee(:)>minspikes;
    postflags.placeRL=posthighbinsRL(:)<maxhighbins & postdata.place_scoreRL(:)>minplacescore & postdata.spikes_per_RLbee(:)>minspikes;
    
    predata.isplaceLR_N=preflags.placeLR;
    predata.isplaceRL_N=preflags.placeRL;
    predata.isplace_N=predata.isplaceLR_N | predata.isplaceRL_N;
    data.isplaceLR_N=shkflags.placeLR;
    data.isplaceRL_N=shkflags.placeRL;
    data.isplace_N=data.isplaceLR_N | data.isplaceRL_N;
    postdata.isplaceLR_N=postflags.placeLR;
    postdata.isplaceRL_N=postflags.placeRL;
    postdata.isplace_N=postdata.isplaceLR_N | postdata.isplaceRL_N;
    
    if ~(shockflag==0)
        if r<17 | r>=26
            if ~(r>6 & r<13) & ~(r==5)
                analyze_shock_responses_drugfree;
            end
        else
            analyze_shock_responses_scopo;
        end
    end
    
    
    
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
    shkdexLR=find(data.isplaceLR_N);
    shkdexRL=find(data.isplaceRL_N);
    for i=shkdexLR(:)'
        if ~isempty(find(cmap(:,shk_column)==i))
            shkplacematrowLR(i)=find(cmap(:,shk_column)==i);
        end
    end
    for i=shkdexRL(:)'
        if ~isempty(find(cmap(:,shk_column)==i))
            shkplacematrowRL(i)=find(cmap(:,shk_column)==i);
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
    
    analysis_field_Nplace(r,:)=[sum(predata.isplace_N) sum(data.isplace_N) sum(postdata.isplace_N)];
    analysis_field_Nnonplace(r,:)=[sum(~predata.isplace_N) sum(~data.isplace_N) sum(~postdata.isplace_N)];
    analysis_field_Nonedir(r,:)=[sum( (~predata.isplaceLR_N & predata.isplaceRL_N) | (predata.isplaceLR_N & ~predata.isplaceRL_N) ) sum( (~data.isplaceLR_N & data.isplaceRL_N) | (data.isplaceLR_N & ~data.isplaceRL_N) ) sum( (~postdata.isplaceLR_N & postdata.isplaceRL_N) | (postdata.isplaceLR_N & ~postdata.isplaceRL_N) )];
    analysis_field_Ntwodir(r,:)=[sum( predata.isplaceLR_N & predata.isplaceRL_N  ) sum( data.isplaceLR_N & data.isplaceRL_N ) sum( postdata.isplaceLR_N & postdata.isplaceRL_N )];
    
    preplacematrowLR=preplacematrowLR(preplacematrowLR>0);
    preplacematrowRL=preplacematrowRL(preplacematrowRL>0);
    shkplacematrowLR=shkplacematrowLR(shkplacematrowLR>0);
    shkplacematrowRL=shkplacematrowRL(shkplacematrowRL>0);
    postplacematrowLR=postplacematrowLR(postplacematrowLR>0);
    postplacematrowRL=postplacematrowRL(postplacematrowRL>0);
    analysis_field_Nrecur(r,:)=[length([unique([preplacematrowLR shkplacematrowLR]) unique([preplacematrowRL shkplacematrowRL])]) length([unique([postplacematrowLR shkplacematrowLR]) unique([postplacematrowRL shkplacematrowRL])])];
    analysis_field_Nnonrecur(r,:)=[length(setdiff(preplacematrowLR, shkplacematrowLR))+length(setdiff(preplacematrowLR, shkplacematrowLR)) length(setdiff(postplacematrowLR, shkplacematrowLR))+length(setdiff(postplacematrowRL, shkplacematrowRL))];
    
    %determine which cells recorded in the shock session also recurred in the pre session
    
    recur_shockpre_matdex=find(cmap(:,pre_column)>0 & cmap(:,shk_column)>0); %rows of cmap map that contain cells recurring in shock & pre sessions
    recur_shockpre_preflag=false(1,size(predata.deconv,1)); %initialize recurrence flags
    recur_shockpre_preflag(cmap(recur_shockpre_matdex,pre_column))=true; %logical vector of cells recorded in SHOCK session, flagged if they recur in the PRE session
    recur_shockpre_shockflag=false(1,size(data.deconv,1)); %initialize recurrence flags
    recur_shockpre_shockflag(cmap(recur_shockpre_matdex,shk_column))=true; %logical vector of cells recorded in PRE session, flagged if they recur in the SHOCK session
    
    if ~(shockflag==0)
        if ~(r>6 & r<13) & ~(r==5) %if not Hipp15 or barrier
            if shockflag==1
                if shockpolarity<0
                    recur_shockpre_shockflag(~SHKinh)=false;
                    for i=1:length(recur_shockpre_matdex)
                        if ~SHKinh(cmap(recur_shockpre_matdex(i),shk_column))
                            recur_shockpre_matdex(i)=0;
                        end
                    end
                else
                    recur_shockpre_shockflag(~SHKexc)=false;
                    for i=1:length(recur_shockpre_matdex)
                        if ~SHKexc(cmap(recur_shockpre_matdex(i),shk_column))
                            recur_shockpre_matdex(i)=0;
                        end
                    end
                end
            end
            if shockflag==-1
                recur_shockpre_shockflag(SHKexc)=false;
                for i=1:length(recur_shockpre_matdex)
                    if SHKexc(cmap(recur_shockpre_matdex(i),shk_column)) %| SHKinh(cmap(recur_shockpre_matdex(i),shk_column))
                        recur_shockpre_matdex(i)=0;
                    end
                end
            end
        end
    end
    recur_shockpre_matdex=recur_shockpre_matdex(find(recur_shockpre_matdex>0)); %drop non recurring cells from the list
    
    if length(recur_shockpre_matdex)>0
        
        for d=1:length(recur_shockpre_matdex) %loop through recurring cells
            hmap_pre(d,:)=[predata.dcurve_LR(cmap(recur_shockpre_matdex(d),pre_column),:) predata.dcurve_RL(cmap(recur_shockpre_matdex(d),pre_column),:)];
            hmap_shockpre(d,:)=[data.dcurve_LR(cmap(recur_shockpre_matdex(d),shk_column),:) data.dcurve_RL(cmap(recur_shockpre_matdex(d),shk_column),:)];
            prehighbinsLR=sum(predata.dcurve_LR(cmap(recur_shockpre_matdex(d),pre_column),:)>max(predata.dcurve_LR(cmap(recur_shockpre_matdex(d),pre_column),:))*.7);
            shkprehighbinsLR=sum(data.dcurve_LR(cmap(recur_shockpre_matdex(d),shk_column),:)>max(data.dcurve_LR(cmap(recur_shockpre_matdex(d),shk_column),:))*.7);
            prehighbinsRL=sum(predata.dcurve_RL(cmap(recur_shockpre_matdex(d),pre_column),:)>max(predata.dcurve_RL(cmap(recur_shockpre_matdex(d),pre_column),:))*.7);
            shkprehighbinsRL=sum(data.dcurve_RL(cmap(recur_shockpre_matdex(d),shk_column),:)>max(data.dcurve_RL(cmap(recur_shockpre_matdex(d),shk_column),:))*.7);
            isplace_preLR(d) = (prehighbinsLR<maxhighbins & predata.place_scoreLR(cmap(recur_shockpre_matdex(d),pre_column))>minplacescore & predata.spikes_per_LRbee(cmap(recur_shockpre_matdex(d),pre_column))>minspikes )...
                | (shkprehighbinsLR<maxhighbins & data.place_scoreLR(cmap(recur_shockpre_matdex(d),shk_column))>minplacescore & data.spikes_per_LRbee(cmap(recur_shockpre_matdex(d),shk_column))>minspikes );
            isplace_preRL(d) = (prehighbinsRL<maxhighbins & predata.place_scoreRL(cmap(recur_shockpre_matdex(d),pre_column))>minplacescore & predata.spikes_per_RLbee(cmap(recur_shockpre_matdex(d),pre_column))>minspikes )...
                | (shkprehighbinsRL<maxhighbins & data.place_scoreRL(cmap(recur_shockpre_matdex(d),shk_column))>minplacescore & data.spikes_per_RLbee(cmap(recur_shockpre_matdex(d),shk_column))>minspikes );
            isplace_pre(d) = (isplace_preLR(d) | isplace_preRL(d));
            nmax=max([hmap_pre(d,:) hmap_shockpre(d,:)]);
            
            [temp, peakbin_preLR(d)]=max(predata.dcurve_LR(d,:));
            [temp, peakbin_preRL(d)]=max(predata.dcurve_RL(d,:));
            [temp, peakbin_shockpreLR(d)]=max(data.dcurve_LR(d,:));
            [temp, peakbin_shockpreRL(d)]=max(data.dcurve_RL(d,:));
            
        end
        
    else
        
        hmap_pre=[];
        hmap_shockpre=[];
        isplace_preLR=[];
        isplace_preRL=[];
        isplace_pre=[];
        nmax=0;
        peakbin_preLR=[];
        peakbin_preRL=[];
        peakbin_shockpreLR=[];
        peakbin_shockpreRL=[];
        
    end
    
    %determine which cells recorded in the shock session also recurred in the post session
    recur_shockpost_matdex=find(cmap(:,post_column)>0 & cmap(:,shk_column)>0); %rows of cell2index map that contain cells recurring in shock & pre sessions
    recur_shockpost_postflag=false(1,size(postdata.deconv,1)); %initialize recurrence flags
    recur_shockpost_postflag(cmap(recur_shockpost_matdex,post_column))=true; %logical vector of cells recorded in SHOCK session, flagged if they recur in the other session
    recur_shockpost_shockflag=false(1,size(data.deconv,1)); %initialize recurrence flags
    recur_shockpost_shockflag(cmap(recur_shockpost_matdex,shk_column))=true; %logical vector of cells recorded in SHOCK session, flagged if they recur in the other session
    
    oldpostdata.shock_code=oldpostdata.place_scoreLR*0;
    for i=1:size(data.deconv,1)
        matdex=find(cmap(:,shk_column)==i);
        if ~isempty(matdex)
            postmatdex=cmap(matdex,post_column);
            if postmatdex>0
                if SHKexc(i)
                    oldpostdata.shock_code(postmatdex)=1;
                else
                    oldpostdata.shock_code(postmatdex)=-1;
                end
            end
        end
    end
    
    if ~(shockflag==0)
        if ~(r>6 & r<13) & ~(r==5)
            if shockflag==1
                if shockpolarity<0
                    recur_shockpost_shockflag(~SHKinh)=false;
                    for i=1:length(recur_shockpost_matdex)
                        if ~SHKinh(cmap(recur_shockpost_matdex(i),shk_column))
                            recur_shockpost_matdex(i)=0;
                        end
                    end
                else
                    recur_shockpost_shockflag(~SHKexc)=false;
                    for i=1:length(recur_shockpost_matdex)
                        if ~SHKexc(cmap(recur_shockpost_matdex(i),shk_column))
                            recur_shockpost_matdex(i)=0;
                        end
                    end
                end
            end
            if shockflag==-1
                recur_shockpost_shockflag(SHKexc)=false;
                for i=1:length(recur_shockpost_matdex)
                    if SHKexc(cmap(recur_shockpost_matdex(i),shk_column)) %| SHKinh(cmap(recur_shockpost_matdex(i),shk_column))
                        recur_shockpost_matdex(i)=0;
                    end
                end
            end
        end
    end
    recur_shockpost_matdex=recur_shockpost_matdex(find(recur_shockpost_matdex>0));
    recur_shockpost_matdex=recur_shockpost_matdex(find(cmap(recur_shockpost_matdex,post_column)<size(postdata.dcurve_LR,1)));
    
    if length(recur_shockpost_matdex)>0
        
        for d=1:length(recur_shockpost_matdex)
            hmap_post(d,:)=[postdata.dcurve_LR(cmap(recur_shockpost_matdex(d),post_column),:) postdata.dcurve_RL(cmap(recur_shockpost_matdex(d),post_column),:)];
            hmap_shockpost(d,:)=[data.dcurve_LR(cmap(recur_shockpost_matdex(d),shk_column),:) data.dcurve_RL(cmap(recur_shockpost_matdex(d),shk_column),:)];
            posthighbinsLR=sum(postdata.dcurve_LR(cmap(recur_shockpost_matdex(d),post_column),:)>max(postdata.dcurve_LR(cmap(recur_shockpost_matdex(d),post_column),:))*.7);
            shkposthighbinsLR=sum(data.dcurve_LR(cmap(recur_shockpost_matdex(d),shk_column),:)>max(data.dcurve_LR(cmap(recur_shockpost_matdex(d),shk_column),:))*.7);
            posthighbinsRL=sum(postdata.dcurve_RL(cmap(recur_shockpost_matdex(d),post_column),:)>max(postdata.dcurve_RL(cmap(recur_shockpost_matdex(d),post_column),:))*.7);
            shkposthighbinsRL=sum(data.dcurve_RL(cmap(recur_shockpost_matdex(d),shk_column),:)>max(data.dcurve_RL(cmap(recur_shockpost_matdex(d),shk_column),:))*.7);
            isplace_postLR(d) = (posthighbinsLR<maxhighbins & postdata.place_scoreLR(cmap(recur_shockpost_matdex(d),post_column))>minplacescore & postdata.spikes_per_LRbee(cmap(recur_shockpost_matdex(d),post_column))>minspikes )...
                | (shkposthighbinsLR<maxhighbins & data.place_scoreLR(cmap(recur_shockpost_matdex(d),shk_column))>minplacescore & data.spikes_per_LRbee(cmap(recur_shockpost_matdex(d),shk_column))>minspikes );
            isplace_postRL(d) = (posthighbinsRL<maxhighbins & postdata.place_scoreRL(cmap(recur_shockpost_matdex(d),post_column))>minplacescore & postdata.spikes_per_RLbee(cmap(recur_shockpost_matdex(d),post_column))>minspikes )...
                | (shkposthighbinsRL<maxhighbins & data.place_scoreRL(cmap(recur_shockpost_matdex(d),shk_column))>minplacescore & data.spikes_per_RLbee(cmap(recur_shockpost_matdex(d),shk_column))>minspikes );
            isplace_post(d) = (isplace_postLR(d) | isplace_postRL(d));
            nmax=max([hmap_post(d,:) hmap_shockpost(d,:)]);
            
            
            [temp, peakbin_postLR(d)]=max(postdata.dcurve_LR(d,:));
            [temp, peakbin_postRL(d)]=max(postdata.dcurve_RL(d,:));
            [temp, peakbin_shockpostLR(d)]=max(data.dcurve_LR(d,:));
            [temp, peakbin_shockpostRL(d)]=max(data.dcurve_RL(d,:));
            
        end
        
    else
        
        hmap_post=[];
        hmap_shockpost=[];
        isplace_postLR=[];
        isplace_postRL=[];
        isplace_post=[];
        nmax=0;
        peakbin_postLR=[];
        peakbin_postRL=[];
        peakbin_shockpostLR=[];
        peakbin_shockpostRL=[];
        
    end
    
    [sum(isplace_preLR) sum(isplace_preRL); sum(isplace_postLR) sum(isplace_postRL)]
    
    if (r<7 & ~(r==5)) | (r>12 & r<17) | (r==26) | (r==27)
        if ~isempty(hmap_pre)
            hm_pre_shk_group=[hm_pre_shk_group; zeros(size(hmap_pre,1),1)+r];
        end
        hm_pre_shk=[hm_pre_shk; hmap_pre];
        hm_shockpre_shk=[hm_shockpre_shk; hmap_shockpre];
        hm_shockpost_shk=[hm_shockpost_shk; hmap_shockpost];
        if ~isempty(hmap_post)
            hm_post_shk_group=[hm_post_shk_group; zeros(size(hmap_post,1),1)+r];
        end
        hm_post_shk=[hm_post_shk; hmap_post];
        hm_pflags_pre_shk=[hm_pflags_pre_shk; [isplace_preLR(:) isplace_preRL(:)] ];
        hm_pflags_post_shk=[hm_pflags_post_shk; [isplace_postLR(:) isplace_postRL(:)] ];
    elseif r>6 & r<13
        hm_pre_bar_group=[hm_pre_bar_group; zeros(size(hmap_pre,1),1)+r];
        hm_pre_bar=[hm_pre_bar; hmap_pre];
        hm_shockpre_bar=[hm_shockpre_bar; hmap_shockpre];
        hm_shockpost_bar=[hm_shockpost_bar; hmap_shockpost];
        hm_post_bar_group=[hm_post_bar_group; zeros(size(hmap_post,1),1)+r];
        hm_post_bar=[hm_post_bar; hmap_post];
        hm_pflags_pre_bar=[hm_pflags_pre_bar; [isplace_preLR(:) isplace_preRL(:)] ];
        hm_pflags_post_bar=[hm_pflags_post_bar; [isplace_postLR(:) isplace_postRL(:)] ];
    elseif r>16 & r<26
        if ~isempty(hmap_pre)
            hm_pre_scpshk_group=[hm_pre_scpshk_group; zeros(size(hmap_pre,1),1)+r];
        end
        hm_pre_scpshk=[hm_pre_scpshk; hmap_pre];
        hm_shockpre_scpshk=[hm_shockpre_scpshk; hmap_shockpre];
        hm_shockpost_scpshk=[hm_shockpost_scpshk; hmap_shockpost];
        if ~isempty(hmap_post)
            hm_post_scpshk_group=[hm_post_scpshk_group; zeros(size(hmap_post,1),1)+r];
        end
        hm_post_scpshk=[hm_post_scpshk; hmap_post];
        hm_pflags_pre_scpshk=[hm_pflags_pre_scpshk; [isplace_preLR(:) isplace_preRL(:)] ];
        hm_pflags_post_scpshk=[hm_pflags_post_scpshk; [isplace_postLR(:) isplace_postRL(:)] ];
    end
    
    %-------- correlation matrices (do this presort to avoid confusion)
    [allmat_preR, allmat_preP] = pv_heatmap_cormatrix(hmap_pre(:,:),hmap_shockpre(:,:));
    [allmat_postR, allmat_postP] = pv_heatmap_cormatrix(hmap_shockpost(:,:),hmap_post(:,:));
    
    try
        [pcorrLR_preR, pcorrLR_preP] = pv_heatmap_tuningcorr([hmap_pre(isplace_preLR,1:23)],[hmap_shockpre(isplace_preLR,1:23)]);
    catch
        pcorrLR_preR=[];
        pcorrLR_preP=[];
    end
    
    try
        [pcorrLR_postR, pcorrLR_postP] = pv_heatmap_tuningcorr([hmap_post(isplace_postLR,1:23)],[hmap_shockpost(isplace_postLR,1:23)]);
    catch
        pcorrLR_postR=[];
        pcorrLR_postP=[];
    end
    
    try
        [pcorrRL_preR, pcorrRL_preP] = pv_heatmap_tuningcorr([hmap_pre(isplace_preRL,24:end)],[hmap_shockpre(isplace_preRL,24:end)]);
    catch
        pcorrRL_preR=[];
        pcorrRL_preP=[];
    end
    
    try
        [pcorrRL_postR, pcorrRL_postP] = pv_heatmap_tuningcorr([hmap_post(isplace_postRL,24:end)],[hmap_shockpost(isplace_postRL,24:end)]);
    catch
        pcorrRL_postR=[];
        pcorrRL_postP=[];
    end
    
    try
        [global_p(r), h]=ranksum([pcorrLR_preR pcorrRL_preR],[pcorrLR_postR pcorrRL_postR]);
    end
    
    analysis_median_tcR(r,:)=[nanmedian([pcorrLR_preR pcorrRL_preR]) nanmedian([pcorrLR_postR pcorrRL_postR])];
    
    globalR_pre=[pcorrLR_preR'; pcorrRL_preR'];
    globalR_post=[pcorrLR_postR'; pcorrRL_postR'];
    
        [pmatLR_preR, pmatLR_preP] = pv_heatmap_cormatrix([hmap_pre(isplace_preLR,1:23)],[hmap_shockpre(isplace_preLR,1:23)]);
        [pmatLR_postR, pmatLR_postP] = pv_heatmap_cormatrix([hmap_post(isplace_postLR,1:23)],[hmap_shockpost(isplace_postLR,1:23)]);
        
        [pmatRL_preR, pmatRL_preP] = pv_heatmap_cormatrix([hmap_pre(isplace_preRL,24:end)],[hmap_shockpre(isplace_preRL,24:end)]);
        [pmatRL_postR, pmatRL_postP] = pv_heatmap_cormatrix([hmap_post(isplace_postRL,24:end)],[hmap_shockpost(isplace_postRL,24:end)]);
        
        
        pmat_preR(1,:,:)=pmatLR_preR; pmat_preR(2,:,:)=pmatRL_preR; %pmat_preR=squeeze(nanmean(pmat_preR));
        pmat_postR(1,:,:)=pmatLR_postR; pmat_postR(2,:,:)=pmatRL_postR; %pmat_postR=squeeze(nanmean(pmat_postR));
        pmat_preP(1,:,:)=pmatLR_preP; pmat_preP(2,:,:)=pmatRL_preP; %pmat_preP=squeeze(nanmean(pmat_preP));
        pmat_postP(1,:,:)=pmatLR_postP; pmat_postP(2,:,:)=pmatRL_postP; %pmat_postP=squeeze(nanmean(pmat_postP));
        
        PC_diffR=squeeze(nanmean(pmat_postR))-squeeze(nanmean(pmat_preR));
        PC_diffP=squeeze(nanmean(pmat_postP))-squeeze(nanmean(pmat_preP));
        
        PC_LRdiffR=pmatLR_postR-pmatLR_preR;
        PC_LRdiffP=pmatLR_postP-pmatLR_preP;
        PC_RLdiffR=pmatRL_postR-pmatRL_preR;
        PC_RLdiffP=pmatRL_postP-pmatRL_preP;
        
        
        shortpath_preR=[]; shortpath_postR=[];
        shortpath_preP=[]; shortpath_postP=[];
        
        middle_preR=[]; middle_postR=[];
        middle_preP=[]; middle_postP=[];
        left_preR=[]; left_postR=[];
        left_preP=[]; left_postP=[];
        right_preR=[]; right_postR=[];
        right_preP=[]; right_postP=[];
        
        for i=1:23
            
            PC_diagR(i)=PC_diffR(i,i);
            PC_diagP(i)=PC_diffP(i,i);
            PC_LRdiagR(i)=PC_LRdiffR(i,i);
            PC_LRdiagP(i)=PC_LRdiffP(i,i);
            PC_RLdiagR(i)=PC_RLdiffR(i,i);
            PC_RLdiagP(i)=PC_RLdiffP(i,i);
            
            PC_diag_preR(i)=nanmean(pmat_preR(:,i,i));
            PC_diag_preP(i)=nanmean(pmat_preP(:,i,i));
            PC_diag_postR(i)=nanmean(pmat_postR(:,i,i));
            PC_diag_postP(i)=nanmean(pmat_postP(:,i,i));
            
        end
        
        PC_LRdiagR=conv(PC_LRdiagR,gausswin(9)/sum(gausswin(9)),'same');
        if r==25
            PC_RLdiagR(6)=.46;
            PC_RLdiagR(11)=.14;
            pmatRL_preR(:,6)=mean(pmatRL_preR(:,[5 7])');
            pmatRL_preR(:,11)=mean(pmatRL_preR(:,[10 12])');
        end
        
        PC_RLdiagR=conv(PC_RLdiagR,gausswin(9)/sum(gausswin(9)),'same');
        PC_diag_postR=conv(PC_diag_postR,[.5 1 .5]/2,'same');
        
        dipcentLR=find(PC_LRdiagR(M_zone)==min(PC_LRdiagR(M_zone)));
        dipcentRL=find(PC_RLdiagR(M_zone)==min(PC_RLdiagR(M_zone)));
        
        if ~isempty(dipcentLR)
            adjM_zoneLR=M_zone(1)-1+dipcentLR-M_wid;
            adjM_zoneLR=adjM_zoneLR:(adjM_zoneLR+M_wid*2);
        else
            adjM_zoneLR=M_zone;
        end
        
        if ~isempty(dipcentRL)
            adjM_zoneRL=M_zone(1)-1+dipcentRL-M_wid;
            adjM_zoneRL=adjM_zoneRL:(adjM_zoneRL+M_wid*2);
        else
            adjM_zoneRL=M_zone;
        end
        
        
        analysis_adjzone(r,:)=[adjM_zoneLR(1) adjM_zoneRL(1)];
        
        adj_zone=unique([adjM_zoneLR adjM_zoneRL]);
               
        for i=1:23
            
            if ~(r==25)%if analyzing LR journeys
                if i>=2 & i<22 %not at an end bin
                    shortpath_preR=[shortpath_preR pmatLR_preR(i,i) pmatLR_preR(i+1,i) pmatLR_preR(i,i+1)];
                    shortpath_postR=[shortpath_postR pmatLR_postR(i,i) pmatLR_postR(i+1,i) pmatLR_postR(i,i+1)];
                end
                if i==22 %last non end bin
                    shortpath_preR=[shortpath_preR pmatLR_preR(i,i)];
                    shortpath_postR=[shortpath_postR pmatLR_postR(i,i)];
                end
                if i>=adjM_zoneLR(1) & i<=adjM_zoneLR(end) %in the roaming middle zone
                    middle_preR=[middle_preR pmatLR_preR(i,i) pmatLR_preR(i-1,i) pmatLR_preR(i,i-1)];
                    middle_postR=[middle_postR pmatLR_postR(i,i) pmatLR_postR(i-1,i) pmatLR_postR(i,i-1)];
                    middle_preP=[middle_preP pmatLR_preP(i,i) pmatLR_preP(i-1,i) pmatLR_preP(i,i-1)];
                    middle_postP=[middle_postP pmatLR_postP(i,i) pmatLR_postP(i-1,i) pmatLR_postP(i,i-1)];
                    if i==adjM_zoneLR(end) %last bin of the roaming middle zone
                        middle_preR=[middle_preR pmatLR_preR(i+1,i) pmatLR_preR(i,i+1)];
                        middle_postR=[middle_postR pmatLR_postR(i+1,i) pmatLR_postR(i,i+1)];
                        middle_preP=[middle_preP pmatLR_preP(i+1,i) pmatLR_preP(i,i+1)];
                        middle_postP=[middle_postP pmatLR_postP(i+1,i) pmatLR_postP(i,i+1)];
                    end
                elseif i<=L_zone(end) %in the static left zone
                    left_preR=[left_preR pmatLR_preR(i,i)];
                    left_postR=[left_postR pmatLR_postR(i,i)];
                    left_preP=[left_preP pmatLR_preP(i,i)];
                    left_postP=[left_postP pmatLR_postP(i,i)];
                    if i<L_zone(end) %not at last bin of the static left zone
                        left_preR=[left_preR pmatLR_preR(i+1,i) pmatLR_preR(i,i+1)];
                        left_postR=[left_postR pmatLR_postR(i+1,i) pmatLR_postR(i,i+1)];
                        left_preP=[left_preP pmatLR_preP(i+1,i) pmatLR_preP(i,i+1)];
                        left_postP=[left_postP pmatLR_postP(i+1,i) pmatLR_postP(i,i+1)];
                    end
                elseif i>=R_zone(1) %int the static right zone
                    right_preR=[right_preR pmatLR_preR(i,i)];
                    right_postR=[right_postR pmatLR_postR(i,i)];
                    right_preP=[right_preP pmatLR_preP(i,i)];
                    right_postP=[right_postP pmatLR_postP(i,i)];
                    if i>R_zone(1) %not at first bin of the static right zone
                        right_preR=[right_preR pmatLR_preR(i-1,i) pmatLR_preR(i,i-1)];
                        right_postR=[right_postR pmatLR_postR(i-1,i) pmatLR_postR(i,i-1)];
                        right_preP=[right_preP pmatLR_preP(i-1,i) pmatLR_preP(i,i-1)];
                        right_postP=[right_postP pmatLR_postP(i-1,i) pmatLR_postP(i,i-1)];
                    end
                end
            end
            
            if i>=2 & i<22 %not at an end bin
                shortpath_preR=[shortpath_preR pmatRL_preR(i,i) pmatRL_preR(i+1,i) pmatRL_preR(i,i+1)];
                shortpath_postR=[shortpath_postR pmatRL_postR(i,i) pmatRL_postR(i+1,i) pmatRL_postR(i,i+1)];
            end
            if i==22 %last non end bin
                shortpath_preR=[shortpath_preR pmatRL_preR(i,i)];
                shortpath_postR=[shortpath_postR pmatRL_postR(i,i)];
            end
            if i>=adjM_zoneRL(1) & i<=adjM_zoneRL(end) %in the roaming middle zone
                middle_preR=[middle_preR pmatRL_preR(i,i) pmatRL_preR(i-1,i) pmatRL_preR(i,i-1)];
                middle_postR=[middle_postR pmatRL_postR(i,i) pmatRL_postR(i-1,i) pmatRL_postR(i,i-1)];
                middle_preP=[middle_preP pmatRL_preP(i,i) pmatRL_preP(i-1,i) pmatRL_preP(i,i-1)];
                middle_postP=[middle_postP pmatRL_postP(i,i) pmatRL_postP(i-1,i) pmatRL_postP(i,i-1)];
                if i==adjM_zoneRL(end) %last bin iof the roaming middle zone
                    middle_preR=[middle_preR pmatRL_preR(i+1,i) pmatRL_preR(i,i+1)];
                    middle_postR=[middle_postR pmatRL_postR(i+1,i) pmatRL_postR(i,i+1)];
                    middle_preP=[middle_preP pmatRL_preP(i+1,i) pmatRL_preP(i,i+1)];
                    middle_postP=[middle_postP pmatRL_postP(i+1,i) pmatRL_postP(i,i+1)];
                end
            elseif i<=L_zone(end) %in the static left zone
                left_preR=[left_preR pmatRL_preR(i,i)];
                left_postR=[left_postR pmatRL_postR(i,i)];
                left_preP=[left_preP pmatRL_preP(i,i)];
                left_postP=[left_postP pmatRL_postP(i,i)];
                if i<L_zone(end) %not at last bin of the static left zone
                    left_preR=[left_preR pmatRL_preR(i+1,i) pmatRL_preR(i,i+1)];
                    left_postR=[left_postR pmatRL_postR(i+1,i) pmatRL_postR(i,i+1)];
                    left_preP=[left_preP pmatRL_preP(i+1,i) pmatRL_preP(i,i+1)];
                    left_postP=[left_postP pmatRL_postP(i+1,i) pmatRL_postP(i,i+1)];
                end
            elseif i>=R_zone(1) %int the static right zone
                right_preR=[right_preR pmatRL_preR(i,i)];
                right_postR=[right_postR pmatRL_postR(i,i)];
                right_preP=[right_preP pmatRL_preP(i,i)];
                right_postP=[right_postP pmatRL_postP(i,i)];
                if i>R_zone(1) %not at first bin of the static right zone
                    right_preR=[right_preR pmatRL_preR(i-1,i) pmatRL_preR(i,i-1)];
                    right_postR=[right_postR pmatRL_postR(i-1,i) pmatRL_postR(i,i-1)];
                    right_preP=[right_preP pmatRL_preP(i-1,i) pmatRL_preP(i,i-1)];
                    right_postP=[right_postP pmatRL_postP(i-1,i) pmatRL_postP(i,i-1)];
                end
            end
        end
        
        safe_preR=[left_preR right_preR];
        safe_postR=[left_postR right_postR];
        
        unsafe_preR=middle_preR;
        unsafe_postR=middle_postR;
        
        if r==25
            shortpath_preR=[shortpath_preR*NaN shortpath_preR];
            shortpath_postR=[shortpath_postR*NaN shortpath_postR];
            
            safe_preR=[safe_preR*NaN safe_preR];
            safe_postR=[safe_postR*NaN safe_postR];
            
            unsafe_preR=[unsafe_preR*NaN unsafe_preR];
            unsafe_postR=[unsafe_postR*NaN unsafe_postR];
        end
        
        safe_preP=[left_preP right_preP];
        safe_postP=[left_postP right_postP];
        
        unsafe_preP=middle_preP;
        unsafe_postP=middle_postP;
        
        goodbins=find(~isnan(safe_preR) & ~isnan(safe_postR));
        [hRy, pRt]=ttest(safe_preR(goodbins),safe_postR(goodbins)); tt_safe_slogp=-log10(pRt)*sign(mean(safe_postR(goodbins)-safe_preR(goodbins)));
        [pR, hR]=signrank(safe_preR(goodbins),safe_postR(goodbins)); sr_safe_slogp=-log10(pR)*sign(mean(safe_postR(goodbins)-safe_preR(goodbins)));
        
        goodbins=find(~isnan(unsafe_preR) & ~isnan(unsafe_postR));
        [hRy, pRt]=ttest(unsafe_preR(goodbins),unsafe_postR(goodbins)); tt_unsafe_slogp=-log10(pRt)*sign(mean(unsafe_postR(goodbins)-unsafe_preR(goodbins)));
        [pR, hR]=signrank(unsafe_preR(goodbins),unsafe_postR(goodbins)); sr_unsafe_slogp=-log10(pR)*sign(mean(unsafe_postR(goodbins)-unsafe_preR(goodbins)));
        
        goodbins=find(~isnan(right_preR) & ~isnan(right_postR));
        [hRy, pRt]=ttest(right_preR(goodbins),right_postR(goodbins)); tt_right_slogp=-log10(pRt)*sign(mean(right_postR(goodbins)-right_preR(goodbins)));
        [pR, hR]=signrank(right_preR(goodbins),right_postR(goodbins)); sr_right_slogp=-log10(pR)*sign(mean(right_postR(goodbins)-right_preR(goodbins)));
        
        
        
        
        
        clear stackR stackP;

    
    if r<5 %shocked first
        Rmatrix_place_shock1((r*2-1):(r*2),:,:)=squeeze(pmat_postR)-squeeze(pmat_preR);
        Rmatrix_place_shock1_pre((r*2-1):(r*2),:,:)=squeeze(pmat_preR);
        Rmatrix_place_shock1_post((r*2-1):(r*2),:,:)=squeeze(pmat_postR);
        Pmatrix_place_shock1(r,:,:)=squeeze(PC_diffP);
    elseif r==6 %shocked first
        Rmatrix_place_shock1(9:10,:,:)=squeeze(pmat_postR)-squeeze(pmat_preR);
        Rmatrix_place_shock1_pre(9:10,:,:)=squeeze(pmat_preR);
        Rmatrix_place_shock1_post(9:10,:,:)=squeeze(pmat_postR);
        Pmatrix_place_shock1(5,:,:)=squeeze(PC_diffP);
    elseif r>=7 & r<=12 %barrier
        Rmatrix_place_barrier(((r-6)*2-1):((r-6)*2),:,:)=squeeze(pmat_postR)-squeeze(pmat_preR);
        Rmatrix_place_barrier_pre(((r-6)*2-1):((r-6)*2),:,:)=squeeze(pmat_preR);
        Rmatrix_place_barrier_post(((r-6)*2-1):((r-6)*2),:,:)=squeeze(pmat_postR);
        Pmatrix_place_barrier(r-6,:,:)=squeeze(PC_diffP);
    elseif r<=16 | r>=26%first shock after barrier
        if r<=16
            rf=7;
        else
            rf=16;
        end
        Rmatrix_place_shock1(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR)-squeeze(pmat_preR);
        Rmatrix_place_shock1_pre(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_preR);
        Rmatrix_place_shock1_post(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR);
        Pmatrix_place_shock1(r-rf,:,:)=squeeze(PC_diffP);
    elseif r<=25 %scop + first shock
        rf=16;
        Rmatrix_place_scopshock(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR)-squeeze(pmat_preR);
        Rmatrix_place_scopshock_pre(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_preR);
        Rmatrix_place_scopshock_post(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR);
        Pmatrix_place_scopshock(r-16,:,:)=squeeze(PC_diffP);
    else %scop alone
        Rmatrix_place_scop(r-24,:,:)=squeeze(PC_diffR);
        Rmatrix_place_scop_pre(r-24,:,:)=squeeze(pmat_preR);
        Rmatrix_place_scop_post(r-24,:,:)=squeeze(pmat_postR);
        Pmatrix_place_scop(r-24,:,:)=squeeze(PC_diffP);
    end
    
    %-------------------- plot sorted heatmaps
    try
        hmapLR_pre=hmap_pre(isplace_preLR,1:23);
        hmapLR_shockpre=hmap_shockpre(isplace_preLR,1:23);
        [temp, peak]=max(hmapLR_shockpre');
        [temp, peak2]=max(hmapLR_pre');
        
        clear hmap_250 hmap_shock250;
        for iii=1:sum(isplace_preLR)
            hmap_250(iii,:)=interp1(1:23,hmapLR_pre(iii,:),linspace(1,23,250));
            hmap_shock250(iii,:)=interp1(1:23,hmapLR_shockpre(iii,:),linspace(1,23,250));
        end
        [temp, peak_250]=max(hmap_shock250');
        [temp, peak2_250]=max(hmap_250');
        shiftdistLR_pre=abs(peak_250-peak2_250);
        
        peak=peak2;
        hmapLR_pre=sortrows([peak(:) hmapLR_pre],1);
        hmapLR_shockpre=sortrows([peak(:) hmapLR_shockpre],1);
        
        hmapRL_pre=hmap_pre(isplace_preRL,24:end);
        hmapRL_shockpre=hmap_shockpre(isplace_preRL,24:end);
        [temp, peak]=max(hmapRL_shockpre');
        [temp, peak2]=max(hmapRL_pre');
        
        clear hmap_250 hmap_shock250;
        for iii=1:sum(isplace_preRL)
            hmap_250(iii,:)=interp1(1:23,hmapRL_pre(iii,:),linspace(1,23,250));
            hmap_shock250(iii,:)=interp1(1:23,hmapRL_shockpre(iii,:),linspace(1,23,250));
        end
        [temp, peak_250]=max(hmap_shock250');
        [temp, peak2_250]=max(hmap_250');
        shiftdistRL_pre=abs(peak_250-peak2_250);
        
        peak=peak2;
        hmapRL_pre=sortrows([peak(:) hmapRL_pre],1);
        hmapRL_shockpre=sortrows([peak(:) hmapRL_shockpre],1);
        
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
        
        peak=peak2;
        hmapLR_post=sortrows([peak(:) hmapLR_post],1);
        hmapLR_shockpost=sortrows([peak(:) hmapLR_shockpost],1);
        
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
        
        peak=peak2;
        hmapRL_post=sortrows([peak(:) hmapRL_post],1);
        hmapRL_shockpost=sortrows([peak(:) hmapRL_shockpost],1);
        
        [pre_v_shockLR_R, pre_v_shockLR_P] = pv_heatmap_cormatrix(hmapLR_pre(:,2:end)',hmapLR_shockpre(:,2:end)');
        [pre_v_shockRL_R, pre_v_shockRL_P] = pv_heatmap_cormatrix(hmapRL_pre(:,2:end)',hmapRL_shockpre(:,2:end)');
        [post_v_shockLR_R, post_v_shockLR_P] = pv_heatmap_cormatrix(hmapLR_post(:,2:end)',hmapLR_shockpost(:,2:end)');
        [post_v_shockRL_R, post_v_shockRL_P] = pv_heatmap_cormatrix(hmapRL_post(:,2:end)',hmapRL_shockpost(:,2:end)');
        
        
        
        isplace_allpre = (predata.place_scoreLR>minplacescore & predata.spikes_per_LRbee>minspikes) | (predata.place_scoreRL>minplacescore & predata.spikes_per_RLbee>minspikes);
        isplace_allshk = (data.place_scoreLR>minplacescore & data.spikes_per_LRbee>minspikes) | (data.place_scoreRL>minplacescore & data.spikes_per_RLbee>minspikes);
        isplace_allpost = (postdata.place_scoreLR>minplacescore & postdata.spikes_per_LRbee>minspikes) | (postdata.place_scoreRL>minplacescore & postdata.spikes_per_RLbee>minspikes);
        
        preplace_recur=length(isplace_allshk & recur_shockpre_shockflag)/(sum(isplace_allpre)+sum(isplace_allshk));
        postplace_recur=length(isplace_allshk & recur_shockpost_shockflag)/(sum(isplace_allpost)+sum(isplace_allshk));
        
        
        prefs=nanmean([sum([predata.dcurve_LR(predata.isplaceLR,:)>repmat(max(predata.dcurve_LR(predata.isplaceLR,:)')',1,23)*.7]') sum([predata.dcurve_RL(predata.isplaceRL,:)>repmat(max(predata.dcurve_RL(predata.isplaceRL,:)')',1,23)*.7]')]*12);
        shkfs=nanmean([sum([data.dcurve_LR(data.isplaceLR,:)>repmat(max(data.dcurve_LR(data.isplaceLR,:)')',1,23)*.7]') sum([data.dcurve_RL(data.isplaceRL,:)>repmat(max(data.dcurve_RL(data.isplaceRL,:)')',1,23)*.7]')]*12);
        postfs=nanmean([sum([postdata.dcurve_LR(postdata.isplaceLR,:)>repmat(max(postdata.dcurve_LR(postdata.isplaceLR,:)')',1,23)*.7]') sum([postdata.dcurve_RL(postdata.isplaceRL,:)>repmat(max(postdata.dcurve_RL(postdata.isplaceRL,:)')',1,23)*.7]')]*12);
        
        prespd=250/nanmedian([predata.LRint(:,2)-predata.LRint(:,1); predata.RLint(:,2)-predata.RLint(:,1)]/11.4);
        shkspd=250/nanmedian([data.LRint(:,2)-data.LRint(:,1); data.RLint(:,2)-data.RLint(:,1)]/11.4);
        postspd=250/nanmedian([postdata.LRint(:,2)-postdata.LRint(:,1); postdata.RLint(:,2)-postdata.RLint(:,1)]/11.4);
        
        preocc=median([predata.vmap_LR predata.vmap_RL])/1000;
        shkocc=median([data.vmap_LR data.vmap_RL])/1000;
        postocc=median([postdata.vmap_LR postdata.vmap_RL])/1000;
        
        predata.isplace=predata.isplaceLR | predata.isplaceRL;
        data.isplace=data.isplaceLR | data.isplaceRL;
        postdata.isplace=postdata.isplaceLR | postdata.isplaceRL;
        
    catch
        shiftdistLR_pre=[]; shiftdistLR_post=[]; globalR_pre=[];
        shortpath_preR=[]; safe_preR=[]; unsafe_preR=[];
        shiftdistRL_pre=[]; shiftdistRL_post=[]; globalR_post=[];
        shortpath_postR=[]; safe_postR=[]; unsafe_postR=[];
        sr_safe_slogp=NaN;  sr_unsafe_slogp=NaN;
        prefs=NaN; shkfs=NaN; postfs=NaN;
        prespd=NaN; shkspd=NaN; postspd=NaN;
        
    end
    
    if ismember(r,[1:4 6 13:16 26:27])
        
        
        analysis_shift_dfpre=[analysis_shift_dfpre; [shiftdistLR_pre'; shiftdistRL_pre']];
        analysis_shift_dfpost=[analysis_shift_dfpost; [shiftdistLR_post'; shiftdistRL_post']];

        
        dfpre_placecellmap=[dfpre_placecellmap; [predata.dcurve_LR(predata.isplace_N,:) predata.dcurve_RL(predata.isplace_N,:)]];
        dfshk_placecellmap=[dfshk_placecellmap; [data.dcurve_LR(data.isplace_N,:) data.dcurve_RL(data.isplace_N,:)]];
        dfpost_placecellmap=[dfpost_placecellmap; [postdata.dcurve_LR(postdata.isplace_N,:) postdata.dcurve_RL(postdata.isplace_N,:)]];
        
        hm_pre_shk_shortpath=[hm_pre_shk_shortpath; shortpath_preR];
        hm_post_shk_shortpath=[hm_post_shk_shortpath; shortpath_postR];
        
        hm_pre_shk_safe=[hm_pre_shk_safe; safe_preR];
        hm_post_shk_safe=[hm_post_shk_safe; safe_postR];
        
        hm_pre_shk_unsafe=[hm_pre_shk_unsafe; unsafe_preR];
        hm_post_shk_unsafe=[hm_post_shk_unsafe; unsafe_postR];
        
    end
    if ismember(r,[17:25])
        analysis_shift_scpre=[analysis_shift_scpre; [shiftdistLR_pre'; shiftdistRL_pre']];
        analysis_shift_scpost=[analysis_shift_scpost; [shiftdistLR_post'; shiftdistRL_post']];

        
        scpre_placecellmap=[scpre_placecellmap; [predata.dcurve_LR(predata.isplace_N,:) predata.dcurve_RL(predata.isplace_N,:)]];
        scshk_placecellmap=[scshk_placecellmap; [data.dcurve_LR(data.isplace_N,:) data.dcurve_RL(data.isplace_N,:)]];
        scpost_placecellmap=[scpost_placecellmap; [postdata.dcurve_LR(postdata.isplace_N,:) postdata.dcurve_RL(postdata.isplace_N,:)]];
        
        hm_pre_scp_shortpath=[hm_pre_scp_shortpath; shortpath_preR];
        hm_post_scp_shortpath=[hm_post_scp_shortpath; shortpath_postR];
        
        hm_pre_scp_safe=[hm_pre_scp_safe; safe_preR];
        hm_post_scp_safe=[hm_post_scp_safe; safe_postR];
        
        hm_pre_scp_unsafe=[hm_pre_scp_unsafe; unsafe_preR];
        hm_post_scp_unsafe=[hm_post_scp_unsafe; unsafe_postR];
    end
    if ismember(r,[7:12])
        analysis_shift_barpre=[analysis_shift_barpre; [shiftdistLR_pre'; shiftdistRL_pre']];
        analysis_shift_barpost=[analysis_shift_barpost; [shiftdistLR_post'; shiftdistRL_post']];
        
        analysis_globalR_barpre=[analysis_globalR_barpre; globalR_pre];
        analysis_globalR_barpost=[analysis_globalR_barpost; globalR_post];
        
        barpre_placecellmap=[barpre_placecellmap; [predata.dcurve_LR(predata.isplace_N,:) predata.dcurve_RL(predata.isplace_N,:)]];
        barshk_placecellmap=[barshk_placecellmap; [data.dcurve_LR(data.isplace_N,:) data.dcurve_RL(data.isplace_N,:)]];
        barpost_placecellmap=[barpost_placecellmap; [postdata.dcurve_LR(postdata.isplace_N,:) postdata.dcurve_RL(postdata.isplace_N,:)]];
        
        hm_pre_bar_shortpath=[hm_pre_bar_shortpath; shortpath_preR];
        hm_post_bar_shortpath=[hm_post_bar_shortpath; shortpath_postR];
        
        hm_pre_bar_safe=[hm_pre_bar_safe; safe_preR];
        hm_post_bar_safe=[hm_post_bar_safe; safe_postR];
        
        hm_pre_bar_unsafe=[hm_pre_bar_unsafe; unsafe_preR];
        hm_post_bar_unsafe=[hm_post_bar_unsafe; unsafe_postR];
    end
    
    try
        analysis_results(r,:)=[r sum(predata.isplaceLR_N | predata.isplaceRL_N)/length(predata.isplaceLR_N) sum(data.isplaceLR_N | data.isplaceRL_N)/length(data.isplaceLR_N) sum(postdata.isplaceLR_N | postdata.isplaceRL_N)/length(postdata.isplaceLR_N) ...
            sum(predata.isplaceLR_N | predata.isplaceRL_N)  sum(data.isplaceLR_N | data.isplaceRL_N)  sum(postdata.isplaceLR_N | postdata.isplaceRL_N) ... %5 6 7
            nanmedian([predata.dcurve_LR_bps(predata.isplaceLR_N) predata.dcurve_RL_bps(predata.isplaceRL_N)]) nanmedian([data.dcurve_LR_bps(data.isplaceLR_N) data.dcurve_RL_bps(data.isplaceRL_N)]) nanmedian([postdata.dcurve_LR_bps(postdata.isplaceLR_N) postdata.dcurve_RL_bps(postdata.isplaceRL_N)])... %8 9 10
            0 ... % 11
            length(predata.LRint)+length(predata.RLint) length(data.LRint)+length(data.RLint) length(postdata.LRint)+length(postdata.RLint) ...  % 12 13 14
            nanmedian(max([predata.dcurve_LR(predata.isplace,:) predata.dcurve_RL(predata.isplace,:)]')) nanmedian(max([data.dcurve_LR(data.isplace,:) data.dcurve_RL(data.isplace,:)]')) nanmedian(max([postdata.dcurve_LR(postdata.isplace,:) postdata.dcurve_RL(postdata.isplace,:)]')) ... %15 16 17
            nanmedian(safe_preR)  nanmedian(safe_postR)  nanmedian(unsafe_preR) ... %18 19 20
            sr_safe_slogp  sr_unsafe_slogp  nanmedian(unsafe_postR) ... %21 22 23
            nanmedian(globalR_pre) nanmedian(globalR_post)  ... %24 25
            prefs shkfs postfs ... %26 27 28
            prespd shkspd postspd ...  %29 30 31
            nanmedian(shortpath_preR) nanmedian(shortpath_postR) ...
            nanmedian([shiftdistLR_pre'; shiftdistRL_pre']) nanmedian([shiftdistLR_post'; shiftdistRL_post'])]; %32 33
    catch
        analysis_results(r,:)=[r sum(predata.isplaceLR_N | predata.isplaceRL_N)/length(predata.isplaceLR_N) sum(data.isplaceLR_N | data.isplaceRL_N)/length(data.isplaceLR_N) sum(postdata.isplaceLR_N | postdata.isplaceRL_N)/length(postdata.isplaceLR_N) ...
            sum(predata.isplaceLR_N | predata.isplaceRL_N)  sum(data.isplaceLR_N | data.isplaceRL_N)  sum(postdata.isplaceLR_N | postdata.isplaceRL_N) ... %5 6 7
            nanmedian([predata.dcurve_LR_bps(predata.isplaceLR_N) predata.dcurve_RL_bps(predata.isplaceRL_N)]) nanmedian([data.dcurve_LR_bps(data.isplaceLR_N) data.dcurve_RL_bps(data.isplaceRL_N)]) nanmedian([postdata.dcurve_LR_bps(postdata.isplaceLR_N) postdata.dcurve_RL_bps(postdata.isplaceRL_N)])... %8 9 10
            0 ... % 11
            length(predata.LRint)+length(predata.RLint) length(data.LRint)+length(data.RLint) length(postdata.LRint)+length(postdata.RLint) ...  % 12 13 14
            nanmedian(max([predata.dcurve_LR(predata.isplace_N,:) predata.dcurve_RL(predata.isplace_N,:)]')) nanmedian(max([data.dcurve_LR(data.isplace_N,:) data.dcurve_RL(data.isplace_N,:)]')) nanmedian(max([postdata.dcurve_LR(postdata.isplace_N,:) postdata.dcurve_RL(postdata.isplace_N,:)]')) ... %15 16 17
            nanmedian(safe_preR)  nanmedian(safe_postR)  nanmedian(unsafe_preR) ... %18 19 20
            NaN NaN  nanmedian(unsafe_postR) ... %21 22 23
            nanmedian(globalR_pre) nanmedian(globalR_post)  ... %24 25
            prefs shkfs postfs ... %26 27 28
            prespd shkspd postspd ...  %29 30 31
            nanmedian(shortpath_preR) nanmedian(shortpath_postR) ...
            nanmedian([shiftdistLR_pre'; shiftdistRL_pre']) nanmedian([shiftdistLR_post'; shiftdistRL_post'])]; %32 33
    end
    try
        [p_pre, h] = signrank(globalR_pre);
        [p_post, h] = signrank(globalR_post);
        
        [h, p_pre_sp] = ttest(shortpath_preR);
        [h, p_post_sp] = ttest(shortpath_postR);
        
        [h, p_pre_safe] = ttest(safe_preR);
        [h, p_post_safe] = ttest(safe_postR);
        
        [h, p_pre_unsafe] = ttest(unsafe_preR);
        [h, p_post_unsafe] = ttest(unsafe_postR);
    catch
        p_pre_sp = [];
        p_post_sp = [];
        p_pre_safe = [];
        p_post_safe = [];
        p_pre_unsafe = [];
        p_post_unsafe = [];
    end
    
    try
        analysis_sigR(r,:) = [p_pre p_post p_pre_sp p_post_sp p_pre_safe p_post_safe p_pre_unsafe p_post_unsafe];
    catch
        analysis_sigR(r,:) = [NaN NaN NaN NaN NaN NaN NaN NaN];
    end
    clear temp*;
    
    try
        
        for d=1:length(recur_shockpre_matdex)
            temp1(d)=nanmean(data.dcurve_LR(cmap(recur_shockpre_matdex(d),shk_column),adjM_zoneLR))-nanmean(predata.dcurve_LR(cmap(recur_shockpre_matdex(d),pre_column),adjM_zoneLR));
            temp2(d)=nanmean(data.dcurve_RL(cmap(recur_shockpre_matdex(d),shk_column),adjM_zoneRL))-nanmean(predata.dcurve_RL(cmap(recur_shockpre_matdex(d),pre_column),adjM_zoneRL));
        end
        
        for d=1:length(recur_shockpost_matdex)
            temp3(d)=nanmean(postdata.dcurve_LR(cmap(recur_shockpost_matdex(d),post_column),adjM_zoneLR))-nanmean(data.dcurve_LR(cmap(recur_shockpost_matdex(d),shk_column),adjM_zoneLR));
            temp4(d)=nanmean(postdata.dcurve_RL(cmap(recur_shockpost_matdex(d),post_column),adjM_zoneRL))-nanmean(data.dcurve_RL(cmap(recur_shockpost_matdex(d),shk_column),adjM_zoneRL));
        end
        
        analysis_deltamax_pre_unsafe(r)=median([temp1 temp2]);
        analysis_deltamax_post_unsafe(r)=median([temp3 temp4]);
        
        analysis_deltamax_pre_unsafeLR(r)=median(temp1);
        analysis_deltamax_pre_unsafeRL(r)=median(temp2);
        analysis_deltamax_post_unsafeLR(r)=median(temp3);
        analysis_deltamax_post_unsafeRL(r)=median(temp4);
        
        analysis_deltamean_pre_unsafe(r)=mean(abs([temp1 temp2]));
        analysis_deltamean_post_unsafe(r)=mean(abs([temp3 temp4]));
        
        analysis_deltamean_pre_unsafeLR(r)=mean(abs(temp1));
        analysis_deltamean_pre_unsafeRL(r)=mean(abs(temp2));
        analysis_deltamean_post_unsafeLR(r)=mean(abs(temp3));
        analysis_deltamean_post_unsafeRL(r)=mean(abs(temp4));
        
        [r size(predata.deconv,1) size(data.deconv,1) size(postdata.deconv,1)]
        
        
        
        temp=unique([cmap(recur_shockpre_matdex(isplace_pre),shk_column); cmap(recur_shockpost_matdex(isplace_post),shk_column)]);
        event_place_cells(r)=length(temp);
        
        recur_shockpost_shockflag=find(recur_shockpost_shockflag);
        recur_shockpre_shockflag=find(recur_shockpre_shockflag);
        
        recur_shockpre_shockmap=[recur_shockpre_matdex(isplace_pre) hmap_shockpre(isplace_pre,:)];
        recur_shockpost_shockmap=[recur_shockpost_matdex(isplace_post) hmap_shockpost(isplace_post,:)];
        recurplace_shockmap=[recur_shockpre_shockmap; recur_shockpost_shockmap];
        [keepercells, ia, ic]=unique(recurplace_shockmap(:,1));
        keepercells=unique(ia);
        recurplace_shockmap=recurplace_shockmap(keepercells,2:end);
        
        predata.dcurve_LR_bps=max(predata.dcurve_LR');
        predata.dcurve_RL_bps=max(predata.dcurve_RL');
        
        data.dcurve_LR_bps=max(data.dcurve_LR');
        data.dcurve_RL_bps=max(data.dcurve_RL');
        
        postdata.dcurve_LR_bps=max(postdata.dcurve_LR');
        postdata.dcurve_RL_bps=max(postdata.dcurve_RL');
        
        if ismember(r,[1:4 6 13:16])
            placemap_shock=[placemap_shock; recurplace_shockmap];
            analysis_placeinfo_shock=[analysis_placeinfo_shock predata.dcurve_LR_bps(predata.isplaceLR_N) predata.dcurve_RL_bps(predata.isplaceRL_N) data.dcurve_LR_bps(data.isplaceLR_N) data.dcurve_RL_bps(data.isplaceRL_N) postdata.dcurve_LR_bps(postdata.isplaceLR_N) postdata.dcurve_RL_bps(postdata.isplaceRL_N)];
            analysis_nonplaceinfo_shock=[analysis_nonplaceinfo_shock predata.dcurve_LR_bps(~predata.isplaceLR_N) predata.dcurve_RL_bps(~predata.isplaceRL_N) data.dcurve_LR_bps(~data.isplaceLR_N) data.dcurve_RL_bps(~data.isplaceRL_N) postdata.dcurve_LR_bps(~postdata.isplaceLR_N) postdata.dcurve_RL_bps(~postdata.isplaceRL_N)];
        elseif ismember(r,[7:12])
            placemap_barrier=[placemap_barrier; recurplace_shockmap];
            analysis_placeinfo_scopshock=[analysis_placeinfo_scopshock predata.dcurve_LR_bps(predata.isplaceLR_N) predata.dcurve_RL_bps(predata.isplaceRL_N) data.dcurve_LR_bps(data.isplaceLR_N) data.dcurve_RL_bps(data.isplaceRL_N) postdata.dcurve_LR_bps(postdata.isplaceLR_N) postdata.dcurve_RL_bps(postdata.isplaceRL_N)];
            analysis_nonplaceinfo_scopshock=[analysis_nonplaceinfo_scopshock predata.dcurve_LR_bps(~predata.isplaceLR_N) predata.dcurve_RL_bps(~predata.isplaceRL_N) data.dcurve_LR_bps(~data.isplaceLR_N) data.dcurve_RL_bps(~data.isplaceRL_N) postdata.dcurve_LR_bps(~postdata.isplaceLR_N) postdata.dcurve_RL_bps(~postdata.isplaceRL_N)];
        else
            placemap_scopshock=[placemap_scopshock; recurplace_shockmap];
            analysis_placeinfo_barrier=[analysis_placeinfo_barrier predata.dcurve_LR_bps(predata.isplaceLR_N) predata.dcurve_RL_bps(predata.isplaceRL_N) data.dcurve_LR_bps(data.isplaceLR_N) data.dcurve_RL_bps(data.isplaceRL_N) postdata.dcurve_LR_bps(postdata.isplaceLR_N) postdata.dcurve_RL_bps(postdata.isplaceRL_N)];
            analysis_nonplaceinfo_barrier=[analysis_nonplaceinfo_barrier predata.dcurve_LR_bps(~predata.isplaceLR_N) predata.dcurve_RL_bps(~predata.isplaceRL_N) data.dcurve_LR_bps(~data.isplaceLR_N) data.dcurve_RL_bps(~data.isplaceRL_N) postdata.dcurve_LR_bps(~postdata.isplaceLR_N) postdata.dcurve_RL_bps(~postdata.isplaceRL_N)];
        end
        
        
    end
    
    postdata=oldpostdata;
    save([postdir rat(r).name '_linear' num2str(cellmat_shockcols(r,postcol)) '_postdata'],'postdata');
    
end %rat loop -----------------------------------------------------------------------------------------------------------------

    


