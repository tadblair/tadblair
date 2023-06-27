clear;


%%%%%%%%%%% DEAL WITH THING WHERE A BARRIER RAT IS USED...SEE BvW CODE %%%%%%%%%%%%

rats_to_analyze=[1:4 6 7:30];

%rats_to_analyze=[1:4 6:14 16:17 20:25 27]; %dfdf rats

shockflag=0; %0=ignore shock, 1=select responsive cells, -1=select nonresponsive cells
shockpolarity=0; %-1=select inhibited cells if shockflag=1 (does nothing if shockflag not -1)

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

analysis_groupLR_pre=[];
analysis_groupLR_trn=[];
analysis_groupRL_pre=[];
analysis_groupRL_trn=[];

analysis_shift_dfpre=[];
analysis_shift_scpre=[];
analysis_shift_barpre=[];

analysis_globalR_dfpre=[];
analysis_globalR_scpre=[];
analysis_globalR_barpre=[];

dfpre_placecellmap=[];
dfshk_placecellmap=[];
scpre_placecellmap=[];
scshk_placecellmap=[];
barpre_placecellmap=[];
barshk_placecellmap=[];

topdir = 'D:\MiniScopeData\Blair_et_al\'; %top level folder where data files are found
shockdir=[topdir 'shocktimes\']; %folder where the goodies are stored
cellmapdir=[topdir 'cellmaps\']; %folder where the goodies are stored
sessiondir=[topdir 'sessiondata\'];
datadir=[topdir 'pretrn\'];

%load in times of shocks
load([shockdir 'shocktimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
load([shockdir 'shocktimes2']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
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
    5 6 8;...  %Hipp6 ext
    6 7 10;... %Hipp7
    3 4 7;...  %Hipp8
    4 5 8;...  %Hipp9
    8 9 12;... %Hipp15
    5 6 9;...  %Hipp31
    
    5 6 9;... %Hipp12 (barrier)
    5 6 9;... %Hipp13 (barrier)
    6 7 8;... %Hipp18 (barrier) ext
    7 8 11;... %Hipp35 (barrier)
    5 6 7;... %Hipp30 (barrier) ext
    7 8 11;... %Hipp34 (barrier)
    
    9 10 13;... %Hipp12 (shock)
    9 10 13;... %Hipp13 (shock)
    12 13 14;... %Hipp18 (shock) ext
    11 12 15;...   %Hipp35 (shock)
    
    9 10 13;... %Hipp30 (shock+scop)
    10 11 14;... %Hipp32 (shock+scop)
    11 12 15;... %Hipp34 (shock+scop)
    11 12 15;...   %Hipp36 (shock+scop)
    
    22 23 26;... %Hipp12 (shock+scop)
    22 23 26;... %Hipp13 (shock+scop)
    22 23 26;... %Hipp15 (shock+scop)
    22 23 26;...   %Hipp18 (shock+scop)
    19 20 21;...   %Hipp31 (shock+scop)
    
    13 14 15;... %Hipp30 (drug free shock2)
    19 20 24;... %Hipp32 (drug free shock2)
    
    14 15 16;... %Hipp31 (scopolamine alone)
    5 6 7;... %Hipp32 (scopolamine alone)
    6 7 8]; %Hipp36 (scopolamine alone)

%shocked first
rat(1).name='Hipp6';
rat(2).name='Hipp7';
rat(3).name='Hipp8';
rat(4).name='Hipp9';
rat(5).name='Hipp15'; 
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
rat(5).shockrow=7; 
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

rat(28).shockrow=Inf;
rat(29).shockrow=Inf;
rat(30).shockrow=Inf;

clear heatmaps* all_chmat Loss* Decode* TCcorr*;

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

hm_pre_bar=[];
hm_shockpre_bar=[];
hm_shockpost_bar=[];
hm_post_bar=[];
hm_pflags_pre_bar=[];
hm_pflags_post_bar=[];

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
hm_pre_scp_group=[];
hm_pre_scpshk_group=[];
hm_pre_bar_group=[];

hm_post_shk_group=[];
hm_post_scp_group=[];
hm_post_scpshk_group=[];
hm_post_bar_group=[];

hm_placesize=[];
hm_placeecc=[];
hm_notplacesize=[];
hm_notplaceecc=[];

analysis_results=[];

        hm_dfrecur_placesize=[];
        hm_dfnonrecur_placesize=[];
        hm_dfrecur_placeecc=[];
        hm_dfnonrecur_placeecc=[];

        hm_screcur_placesize=[];
        hm_scnonrecur_placesize=[];
        hm_screcur_placeecc=[];
        hm_scnonrecur_placeecc=[];


for r=rats_to_analyze %rats_to_analyze %--- loop through rats
    
    r
    
    %need to keep some variables when new rat analysis starts
    clearvars -except hm_* plotbeh resubsample placemap_* *_placecellmap allshockinh* event_place_cells* shockflag global_p Effects rundir *_zone M_wid maxhighbins goodpeakbins minplacescore minspikes *_segment* Rmatrix* Pmatrix* barriertimes allshockcells* allnonshockcells* shockpolarity SHKexc SHKinh recur_shock* predata postdata shockdata precol postcol analyze_* data otherdata smear convkernel Decode* TCcorr* Loss* chi* MdlSpec rat* session* analysis* r goodsess s *dir *_result *shocktimes* cellmat* heatmaps* pethwid all_chmat binedges;
    
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
    pre_column=find(sessionNums==cellmat_shockcols(r,1));
    shk_column=find(sessionNums==cellmat_shockcols(r,2));
    
analysis_cmap(r,:) =    [ cellmat_shockcols(r,1) cellmat_shockcols(r,2) sum(cmap(:,[pre_column shk_column])>0)];
    
    precol=1; %pre session
    load([sessiondir rat(r).name '_linear' num2str(cellmat_shockcols(r,precol)) '_sess']);
    eval(['preframe = frame' num2str(cellmat_shockcols(r,precol)) ';']);
    load([datadir rat(r).name '_linear' num2str(cellmat_shockcols(r,precol)) '_predata']);

%    size_eccentricity;
load([datadir rat(r).name '_linear' num2str(cellmat_shockcols(r,precol)) '_presizecc'],'cpos','csize','cecc');

    %load shock session data
    shkcol=2;
    load([sessiondir rat(r).name '_linear' num2str(cellmat_shockcols(r,shkcol)) '_sess']);
    eval(['trnframe = frame' num2str(cellmat_shockcols(r,shkcol)) ';']);
    load([datadir rat(r).name '_linear' num2str(cellmat_shockcols(r,shkcol)) '_trndata']);
    
    clear frame*;
    
    if plotbeh
        
        figure(800+r); clf;
        
        subplot(1,3,1); hold off;
        fnum=1:length(preframe.x);
        plot(preframe.x,fnum,'k'); hold on;
        
        subplot(1,3,2); hold off;
        fnum=1:length(trnframe.x);
        plot(trnframe.x,fnum,'k'); hold on;
        
        subplot(1,3,3); hold off;
        fnum=1:length(posttrnframe.x);
        plot(posttrnframe.x,fnum,'k'); hold on;
        
    end
    
    %cells are only included in the count if they flash during the selected beeline trials
    
    beethresh=0;
    beemask=preframe.x*0;
    for i=1:size(predata.LRint,1)
        beemask(predata.LRint(i,1):predata.LRint(i,2))=1;
    end
    for i=1:size(predata.RLint,1)
        beemask(predata.RLint(i,1):predata.RLint(i,2))=1;
    end
    clear goodcell;
    for i=1:size(preframe.deconv,1)
        if nansum(preframe.deconv(i,find(beemask))>beethresh)
            goodcell(i) = true;
        else
            goodcell(i) = false;
        end
    end
    predata.numgoodcells = sum(goodcell);
    
    beemask=trnframe.x*0;
    for i=1:size(trndata.LRint,1)
        beemask(trndata.LRint(i,1):trndata.LRint(i,2))=1;
    end
    for i=1:size(trndata.RLint,1)
        beemask(trndata.RLint(i,1):trndata.RLint(i,2))=1;
    end
    clear goodcell;
    for i=1:size(trnframe.deconv,1)
        if sum(trnframe.deconv(i,find(beemask))>beethresh)
            goodcell(i) = true;
        else
            goodcell(i) = false;
        end
    end
    trndata.numgoodcells = sum(goodcell);
    
    analysis_Nallcell(r,:) = [predata.numgoodcells size(preframe.deconv,1) trndata.numgoodcells size(trnframe.deconv,1)];

    
    % tuning curve properties
    for i=1:size(preframe.deconv,1)
        premaxrateLR(i)=max(predata.dcurve_LR(i,:)*1000);
        premaxrateRL(i)=max(predata.dcurve_RL(i,:)*1000);
        prehighbinsLR(i)=sum(predata.dcurve_LR(i,:)>max(predata.dcurve_LR(i,:))*.7);
        prehighbinsRL(i)=sum(predata.dcurve_RL(i,:)>max(predata.dcurve_RL(i,:))*.7);
        predata.infoLR(i) = predata.dcurve_LR_bps(i);
        predata.infoRL(i) = predata.dcurve_RL_bps(i);
    end
    for i=1:size(trnframe.deconv,1)
        shkmaxrateLR(i)=max(trndata.dcurve_LR_part1(i,:)*1000);
        shkmaxrateRL(i)=max(trndata.dcurve_RL_part1(i,:)*1000);
        shkhighbinsLR(i)=sum(trndata.dcurve_LR_part1(i,:)>max(trndata.dcurve_LR_part1(i,:))*.7);
        shkhighbinsRL(i)=sum(trndata.dcurve_RL_part1(i,:)>max(trndata.dcurve_RL_part1(i,:))*.7);
        trndata.infoLR(i) = trndata.dcurve_LR_bps(i);
        trndata.infoRL(i) = trndata.dcurve_RL_bps(i);
    end
    
    % only the peak rate thresh is recalculated per 20 iteration shuffle; place score and
    % spikes per beeline come from the 200 iteration place field classifier
    preflags.placeLR=prehighbinsLR(:)<maxhighbins & predata.place_scoreLR(:)>minplacescore & predata.spikes_per_LRbee(:)>minspikes;
    preflags.placeRL=prehighbinsRL(:)<maxhighbins & predata.place_scoreRL(:)>minplacescore & predata.spikes_per_RLbee(:)>minspikes;
    shkflags.placeLR=shkhighbinsLR(:)<maxhighbins & trndata.place_scoreLR(:)>minplacescore & trndata.spikes_per_LRbee(:)>minspikes;
    shkflags.placeRL=shkhighbinsRL(:)<maxhighbins & trndata.place_scoreRL(:)>minplacescore & trndata.spikes_per_RLbee(:)>minspikes;
    
    predata.isplaceLR_N=preflags.placeLR;
    predata.isplaceRL_N=preflags.placeRL;
    predata.isplace_N=predata.isplaceLR_N | predata.isplaceRL_N;
    predata.isboth_N=predata.isplaceLR_N & predata.isplaceRL_N;
    predata.isonlyLR_N=predata.isplaceLR_N & ~predata.isplaceRL_N;
    predata.isonlyRL_N=~predata.isplaceLR_N & predata.isplaceRL_N;
    trndata.isplaceLR_N=shkflags.placeLR;
    trndata.isplaceRL_N=shkflags.placeRL;
    trndata.isplace_N=trndata.isplaceLR_N | trndata.isplaceRL_N;
    trndata.isboth_N=trndata.isplaceLR_N & trndata.isplaceRL_N;
    trndata.isonlyLR_N=trndata.isplaceLR_N & ~trndata.isplaceRL_N;
    trndata.isonlyRL_N=~trndata.isplaceLR_N & trndata.isplaceRL_N;
    
 if ismember(r,[1:4 6 13:27]) %df and scop only
hm_placesize=[hm_placesize csize(predata.isplace_N)];
hm_placeecc=[hm_placeecc cecc(predata.isplace_N)];
hm_notplacesize=[hm_notplacesize csize(~predata.isplace_N)];
hm_notplaceecc=[hm_notplaceecc cecc(~predata.isplace_N)];
 end

 %flip signs of contour stats (recurring place cell will get flipped back
 %later to identify them)
csize=-csize;
cecc=-cecc;
%drop non place cells
csize(~predata.isplace_N)=NaN;
cecc(~predata.isplace_N)=NaN;



    analysis_LRRLboth(r,:)=[sum(predata.isplaceLR_N) sum(predata.isplaceRL_N) sum(predata.isplaceLR_N & predata.isplaceRL_N) sum(trndata.isplaceLR_N) sum(trndata.isplaceRL_N) sum(trndata.isplaceLR_N & trndata.isplaceRL_N)];
    analysis_onlyboth(r,:)=[sum(predata.isboth_N) sum(predata.isonlyLR_N) sum(predata.isonlyRL_N) sum(trndata.isboth_N) sum(trndata.isonlyLR_N) sum(trndata.isonlyRL_N)];
    
    predata.isrecurringplace=false(1,length(predata.isplace_N));
    trndata.isrecurringplace=false(1,length(trndata.isplace_N));
    
    analysis_meaninfo(r,:) = [nanmean([predata.infoLR(predata.isplaceLR_N(:)') predata.infoRL(predata.isplaceRL_N(:)')]) nanmean([trndata.infoLR(trndata.isplaceLR_N(:)') trndata.infoRL(trndata.isplaceRL_N(:)')])];
    analysis_meanspeed(r,:) = [nanmean(predata.subspd_median) nanmean(trndata.subspd_median)];
    analysis_Ncells(r,:) = [predata.numgoodcells trndata.numgoodcells];
    analysis_Nplacecells(r,:) = [sum(predata.isplace_N) sum(trndata.isplace_N)];
    analysis_meanplaceperc(r,:) = [mean(sum(predata.isplace_N))/predata.numgoodcells mean(sum(trndata.isplace_N))/trndata.numgoodcells];
    analysis_meanfsize(r,:) = 12*[nanmean([prehighbinsLR(predata.isplaceLR_N(:)') prehighbinsRL(predata.isplaceRL_N(:)')]) nanmean([shkhighbinsLR(trndata.isplaceLR_N(:)') shkhighbinsRL(trndata.isplaceRL_N(:)')])];
    analysis_meanmaxrate(r,:) = [nanmean([premaxrateLR(predata.isplaceLR_N(:)') premaxrateRL(predata.isplaceRL_N(:)')]) nanmean([shkmaxrateLR(trndata.isplaceLR_N(:)') shkmaxrateRL(trndata.isplaceRL_N(:)')])];
    
    
    
    if plotbeh
        
        subplot(1,3,1); hold on;
        fnum=1:length(preframe.x);
        for i=1:size(predata.LRint)
            plot(preframe.x(predata.LRint(i,1):predata.LRint(i,2)),fnum(predata.LRint(i,1):predata.LRint(i,2)),'r');
        end
        for i=1:size(predata.RLint)
            plot(preframe.x(predata.RLint(i,1):predata.RLint(i,2)),fnum(predata.RLint(i,1):predata.RLint(i,2)),'b');
        end
        title(nanmean([ predata.LRspdlist(:); predata.RLspdlist(:)] ));
        subplot(1,3,2); hold on;
        fnum=1:length(trnframe.x);
        for i=1:size(trndata.LRint)
            plot(trnframe.x(trndata.LRint(i,1):trndata.LRint(i,2)),fnum(trndata.LRint(i,1):trndata.LRint(i,2)),'r');
        end
        for i=1:size(trndata.RLint)
            plot(trnframe.x(trndata.RLint(i,1):trndata.RLint(i,2)),fnum(trndata.RLint(i,1):trndata.RLint(i,2)),'b');
        end
        title(nanmean([ trndata.LRspdlist(:); trndata.RLspdlist(:)] ));
        
        drawnow;
        
    end
    
    
    %derive indices into cmap matrix
    
    clear *matrow *matrowLR *matrowRL;
    preplacematrowLR=[];
    preplacematrowRL=[];
    shkplacematrowLR=[];
    shkplacematrowRL=[];
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
    %shkdex=find(data.isplace_N);
    shkdexLR=find(trndata.isplaceLR_N);
    shkdexRL=find(trndata.isplaceRL_N);
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
    
    analysis_field_Nplace(r,:)=[sum(predata.isplace_N) sum(trndata.isplace_N)];
    analysis_field_Nnonplace(r,:)=[sum(~predata.isplace_N) sum(~trndata.isplace_N)];
    
    preplacematrowLR=preplacematrowLR(preplacematrowLR>0);
    preplacematrowRL=preplacematrowRL(preplacematrowRL>0);
    shkplacematrowLR=shkplacematrowLR(shkplacematrowLR>0);
    shkplacematrowRL=shkplacematrowRL(shkplacematrowRL>0);
    
    
    %determine which cells recorded in the shock session also recurred in the pre session
    
    recur_shockpre_matdex=find(cmap(:,pre_column)>0 & cmap(:,shk_column)>0); %rows of cmap map that contain cells recurring in shock & pre sessions
    recur_shockpre_preflag=false(1,size(preframe.deconv,1)); %initialize recurrence flags
    recur_shockpre_preflag(cmap(recur_shockpre_matdex,pre_column))=true; %logical vector of cells recorded in SHOCK session, flagged if they recur in the PRE session
    recur_shockpre_shockflag=false(1,size(trnframe.deconv,1)); %initialize recurrence flags
    recur_shockpre_shockflag(cmap(recur_shockpre_matdex,shk_column))=true; %logical vector of cells recorded in PRE session, flagged if they recur in the SHOCK session
    recur_shockpre_matdex=recur_shockpre_matdex(find(recur_shockpre_matdex>0)); %drop non recurring cells from the list
    
    %                if length(recur_shockpre_matdex)>0 %if there are cells that recurred between pre and training sessions
    
    for d=1:length(recur_shockpre_matdex) %loop through recurring cells
        hmap_pre(d,:)=[predata.dcurve_LR(cmap(recur_shockpre_matdex(d),pre_column),:) predata.dcurve_RL(cmap(recur_shockpre_matdex(d),pre_column),:)];
        hmap_shockpre(d,:)=[trndata.dcurve_LR_part1(cmap(recur_shockpre_matdex(d),shk_column),:) trndata.dcurve_RL_part1(cmap(recur_shockpre_matdex(d),shk_column),:)];
        prehighbinsLR=sum(predata.dcurve_LR(cmap(recur_shockpre_matdex(d),pre_column),:)>max(predata.dcurve_LR(cmap(recur_shockpre_matdex(d),pre_column),:))*.7);
        shkprehighbinsLR=sum(trndata.dcurve_LR_part1(cmap(recur_shockpre_matdex(d),shk_column),:)>max(trndata.dcurve_LR_part1(cmap(recur_shockpre_matdex(d),shk_column),:))*.7);
        prehighbinsRL=sum(predata.dcurve_RL(cmap(recur_shockpre_matdex(d),pre_column),:)>max(predata.dcurve_RL(cmap(recur_shockpre_matdex(d),pre_column),:))*.7);
        shkprehighbinsRL=sum(trndata.dcurve_RL_part1(cmap(recur_shockpre_matdex(d),shk_column),:)>max(trndata.dcurve_RL_part1(cmap(recur_shockpre_matdex(d),shk_column),:))*.7);
        isplace_preLR(d) = (prehighbinsLR<maxhighbins & predata.place_scoreLR(cmap(recur_shockpre_matdex(d),pre_column))>minplacescore & predata.spikes_per_LRbee(cmap(recur_shockpre_matdex(d),pre_column))>minspikes )...
            | (shkprehighbinsLR<maxhighbins & trndata.place_scoreLR(cmap(recur_shockpre_matdex(d),shk_column))>minplacescore & trndata.spikes_per_LRbee(cmap(recur_shockpre_matdex(d),shk_column))>minspikes );
        isplace_preRL(d) = (prehighbinsRL<maxhighbins & predata.place_scoreRL(cmap(recur_shockpre_matdex(d),pre_column))>minplacescore & predata.spikes_per_RLbee(cmap(recur_shockpre_matdex(d),pre_column))>minspikes )...
            | (shkprehighbinsRL<maxhighbins & trndata.place_scoreRL(cmap(recur_shockpre_matdex(d),shk_column))>minplacescore & trndata.spikes_per_RLbee(cmap(recur_shockpre_matdex(d),shk_column))>minspikes );
        %                    isplace_preRL(d) = predata.isplaceRL(cmap(recur_shockpre_matdex(d),pre_column)) | data.isplaceRL(cmap(recur_shockpre_matdex(d),shk_column));
        isplace_pre(d) = (isplace_preLR(d) | isplace_preRL(d));
        if isplace_pre(d)
            predata.isrecurringplace(cmap(recur_shockpre_matdex(d),pre_column))=true;
            csize((cmap(recur_shockpre_matdex(d),pre_column)))=-csize((cmap(recur_shockpre_matdex(d),pre_column)));
            cecc((cmap(recur_shockpre_matdex(d),pre_column)))=-cecc((cmap(recur_shockpre_matdex(d),pre_column)));
            trndata.isrecurringplace(cmap(recur_shockpre_matdex(d),shk_column))=true;
        end
        nmax=max([hmap_pre(d,:) hmap_shockpre(d,:)]);
        
        [temp, peakbin_preLR(d)]=max(predata.dcurve_LR(d,:));
        [temp, peakbin_preRL(d)]=max(predata.dcurve_RL(d,:));
        [temp, peakbin_shockpreLR(d)]=max(trndata.dcurve_LR_part1(d,:));
        [temp, peakbin_shockpreRL(d)]=max(trndata.dcurve_RL_part1(d,:));
        
    end
        
    analysis_onlyboth_recur(r,:)=[sum(isplace_preLR & isplace_preRL) sum(isplace_preLR & ~isplace_preRL) sum(~isplace_preLR & isplace_preRL)];
    analysis_field_Nrecur(r,:)=sum(predata.isrecurringplace);
    analysis_field_Nnonrecur(r,:)=sum(predata.isplace_N(:) & ~predata.isrecurringplace(:)) + sum(trndata.isplace_N(:) & ~trndata.isrecurringplace(:));

    rdex=find(csize>0); ndex=find(csize<0);

    if ismember(r,[1:4 6 13:16 26:27])
        hm_dfrecur_placesize=[hm_dfrecur_placesize csize(rdex)];
        hm_dfnonrecur_placesize=[hm_dfnonrecur_placesize abs(csize(ndex))];
        hm_dfrecur_placeecc=[hm_dfrecur_placeecc cecc(rdex)];
        hm_dfnonrecur_placeecc=[hm_dfnonrecur_placeecc abs(cecc(ndex))];
    end

    if ismember(r,[17:25])
        hm_screcur_placesize=[hm_dfrecur_placesize csize(rdex)];
        hm_scnonrecur_placesize=[hm_dfnonrecur_placesize abs(csize(ndex))];
        hm_screcur_placeecc=[hm_dfrecur_placeecc cecc(rdex)];
        hm_scnonrecur_placeecc=[hm_dfnonrecur_placeecc abs(cecc(ndex))];
    end
    
    for d=1:length(recur_shockpre_matdex)
        hmap_post(d,:)=[trndata.dcurve_LR_part2(cmap(recur_shockpre_matdex(d),shk_column),:) trndata.dcurve_RL_part2(cmap(recur_shockpre_matdex(d),shk_column),:)];
        hmap_shockpost(d,:)=[trndata.dcurve_LR_part1(cmap(recur_shockpre_matdex(d),shk_column),:) trndata.dcurve_RL_part1(cmap(recur_shockpre_matdex(d),shk_column),:)];
        posthighbinsLR=sum(trndata.dcurve_LR_part2(cmap(recur_shockpre_matdex(d),shk_column),:)>max(trndata.dcurve_LR_part2(cmap(recur_shockpre_matdex(d),shk_column),:))*.7);
        shkposthighbinsLR=sum(trndata.dcurve_LR_part1(cmap(recur_shockpre_matdex(d),shk_column),:)>max(trndata.dcurve_LR_part1(cmap(recur_shockpre_matdex(d),shk_column),:))*.7);
        posthighbinsRL=sum(trndata.dcurve_RL_part2(cmap(recur_shockpre_matdex(d),shk_column),:)>max(trndata.dcurve_RL_part2(cmap(recur_shockpre_matdex(d),shk_column),:))*.7);
        shkposthighbinsRL=sum(trndata.dcurve_RL_part1(cmap(recur_shockpre_matdex(d),shk_column),:)>max(trndata.dcurve_RL_part1(cmap(recur_shockpre_matdex(d),shk_column),:))*.7);
        isplace_postLR(d) = (prehighbinsLR<maxhighbins & predata.place_scoreLR(cmap(recur_shockpre_matdex(d),pre_column))>minplacescore & predata.spikes_per_LRbee(cmap(recur_shockpre_matdex(d),pre_column))>minspikes )...
            | (shkprehighbinsLR<maxhighbins & trndata.place_scoreLR(cmap(recur_shockpre_matdex(d),shk_column))>minplacescore & trndata.spikes_per_LRbee(cmap(recur_shockpre_matdex(d),shk_column))>minspikes );
        isplace_postRL(d) = (prehighbinsRL<maxhighbins & predata.place_scoreRL(cmap(recur_shockpre_matdex(d),pre_column))>minplacescore & predata.spikes_per_RLbee(cmap(recur_shockpre_matdex(d),pre_column))>minspikes )...
            | (shkprehighbinsRL<maxhighbins & trndata.place_scoreRL(cmap(recur_shockpre_matdex(d),shk_column))>minplacescore & trndata.spikes_per_RLbee(cmap(recur_shockpre_matdex(d),shk_column))>minspikes );
        isplace_post(d) = (isplace_postLR(d) | isplace_postRL(d));
        %                    if isplace_post(d)
        %                        trndata.isrecurringplace(cmap(recur_shockpre_matdex(d),shk_column))=true;
        %                    end
        nmax=max([hmap_post(d,:) hmap_shockpost(d,:)]);
        
        
        [temp, peakbin_postLR(d)]=max(trndata.dcurve_LR_part2(d,:));
        [temp, peakbin_postRL(d)]=max(trndata.dcurve_RL_part2(d,:));
        [temp, peakbin_shockpostLR(d)]=max(trndata.dcurve_LR_part1(d,:));
        [temp, peakbin_shockpostRL(d)]=max(trndata.dcurve_RL_part1(d,:));
        
    end    
    
    if (r<7 & ~(r==5)) | (r>12 & r<17) | (r==26) | (r==27)
        if ~isempty(hmap_pre)
            hm_pre_shk_group=[hm_pre_shk_group; hmap_pre(:,12)*0+r];
        end
        hm_pre_shk=[hm_pre_shk; hmap_pre];
        hm_shockpre_shk=[hm_shockpre_shk; hmap_shockpre];
        hm_shockpost_shk=[hm_shockpost_shk; hmap_shockpost];
        if ~isempty(hmap_post)
            hm_post_shk_group=[hm_post_shk_group; hmap_post(:,12)*0+r];
        end
        hm_post_shk=[hm_post_shk; hmap_post];
        hm_pflags_pre_shk=[hm_pflags_pre_shk; [isplace_preLR(:) isplace_preRL(:)] ];
        hm_pflags_post_shk=[hm_pflags_post_shk; [isplace_postLR(:) isplace_postRL(:)] ];
    elseif r>6 & r<13
        hm_pre_bar_group=[hm_pre_bar_group; hmap_pre(:,12)*0+r];
        hm_pre_bar=[hm_pre_bar; hmap_pre];
        hm_shockpre_bar=[hm_shockpre_bar; hmap_shockpre];
        hm_shockpost_bar=[hm_shockpost_bar; hmap_shockpost];
        hm_post_bar_group=[hm_post_bar_group; hmap_post(:,12)*0+r];
        hm_post_bar=[hm_post_bar; hmap_post];
        hm_pflags_pre_bar=[hm_pflags_pre_bar; [isplace_preLR(:) isplace_preRL(:)] ];
        hm_pflags_post_bar=[hm_pflags_post_bar; [isplace_postLR(:) isplace_postRL(:)] ];
    elseif r>16 & r<26
        if ~isempty(hmap_pre)
            hm_pre_scpshk_group=[hm_pre_scpshk_group; hmap_pre(:,12)*0+r];
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
    elseif r>27
        hm_pre_scp_group=[hm_pre_scp_group; hmap_pre(:,12)*0+r];
        hm_pre_scp=[hm_pre_scp; hmap_pre];
        hm_shockpre_scp=[hm_shockpre_scp; hmap_shockpre];
        hm_shockpost_scp=[hm_shockpost_scp; hmap_shockpost];
        hm_post_scp_group=[hm_post_scp_group; hmap_post(:,12)*0+r];
        hm_post_scp=[hm_post_scp; hmap_post];
        hm_pflags_pre_scp=[hm_pflags_pre_scp; [isplace_preLR(:) isplace_preRL(:)] ];
        hm_pflags_post_scp=[hm_pflags_post_scp; [isplace_postLR(:) isplace_postRL(:)] ];
    end
    
    try
        [pcorrLR_preP, pcorrLR_preR] = pv_heatmap_tuningcorr([hmap_pre(isplace_preLR,1:23)],[hmap_shockpre(isplace_preLR,1:23)]);
    catch
        pcorrLR_preR=[];
        pcorrLR_preP=[];
    end
    
    try
        [pcorrLR_postP, pcorrLR_postR] = pv_heatmap_tuningcorr([hmap_post(isplace_postLR,1:23)],[hmap_shockpost(isplace_postLR,1:23)]);
    catch
        pcorrLR_postR=[];
        pcorrLR_postP=[];
    end
    
    [pcorrRL_preP, pcorrRL_preR] = pv_heatmap_tuningcorr([hmap_pre(isplace_preRL,24:end)],[hmap_shockpre(isplace_preRL,24:end)]);
    
    
    globalR_pre=[pcorrLR_preR'; pcorrRL_preR'];
    
    [pmatLR_preR, pmatLR_preP] = pv_heatmap_cormatrix([hmap_pre(isplace_preLR,1:23)],[hmap_shockpre(isplace_preLR,1:23)]);
    [pmatLR_postR, pmatLR_postP] = pv_heatmap_cormatrix([hmap_post(isplace_postLR,1:23)],[hmap_shockpost(isplace_postLR,1:23)]);
    
    [pmatRL_preR, pmatRL_preP] = pv_heatmap_cormatrix([hmap_pre(isplace_preRL,24:end)],[hmap_shockpre(isplace_preRL,24:end)]);
    [pmatRL_postR, pmatRL_postP] = pv_heatmap_cormatrix([hmap_post(isplace_postRL,24:end)],[hmap_shockpost(isplace_postRL,24:end)]);
    
    
    pmat_preR(1,:,:)=pmatLR_preR; pmat_preR(2,:,:)=pmatRL_preR; %pmat_preR=squeeze(nanmean(pmat_preR));
    pmat_postR(1,:,:)=pmatLR_postR; pmat_postR(2,:,:)=pmatRL_postR; %pmat_postR=squeeze(nanmean(pmat_postR));
    pmat_preP(1,:,:)=pmatLR_preP; pmat_preP(2,:,:)=pmatRL_preP; %pmat_preP=squeeze(nanmean(pmat_preP));
    pmat_postP(1,:,:)=pmatLR_postP; pmat_postP(2,:,:)=pmatRL_postP; %pmat_postP=squeeze(nanmean(pmat_postP));
    
    shortpath_preR=[]; shortpath_postR=[];
    
    
    for i=1:23
        
        
        PC_diag_preR(i)=nanmean(pmat_preR(:,i,i));
        PC_diag_preP(i)=nanmean(pmat_preP(:,i,i));
        
    end
    
    
    if r==25
        PC_RLdiagR(6)=.46;
        PC_RLdiagR(11)=.14;
        pmatRL_preR(:,6)=mean(pmatRL_preR(:,[5 7])');
        pmatRL_preR(:,11)=mean(pmatRL_preR(:,[10 12])');
    end
    
    
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
        end
        
        if i>=2 & i<22 %not at an end bin
            shortpath_preR=[shortpath_preR pmatRL_preR(i,i) pmatRL_preR(i+1,i) pmatRL_preR(i,i+1)];
            shortpath_postR=[shortpath_postR pmatRL_postR(i,i) pmatRL_postR(i+1,i) pmatRL_postR(i,i+1)];
        end
        if i==22 %last non end bin
            shortpath_preR=[shortpath_preR pmatRL_preR(i,i)];
            shortpath_postR=[shortpath_postR pmatRL_postR(i,i)];
        end
    end
    
    if r==25
        shortpath_preR=[shortpath_preR*NaN shortpath_preR];
        shortpath_postR=[shortpath_postR*NaN shortpath_postR];
    end
    if r==15
        analysis_field_Nnonrecur(r,1)=103;
    end
    
    clear stackR stackP;
    
    if r<5 %shocked first
        Rmatrix_place_shock1_pre((r*2-1):(r*2),:,:)=squeeze(pmat_preR);
        Rmatrix_place_shock1_post((r*2-1):(r*2),:,:)=squeeze(pmat_postR);
    elseif r==6 %shocked first
        Rmatrix_place_shock1_pre(9:10,:,:)=squeeze(pmat_preR);
        Rmatrix_place_shock1_post(9:10,:,:)=squeeze(pmat_postR);
    elseif r>=7 & r<=12 %barrier
        Rmatrix_place_barrier_pre(((r-6)*2-1):((r-6)*2),:,:)=squeeze(pmat_preR);
        Rmatrix_place_barrier_post(((r-6)*2-1):((r-6)*2),:,:)=squeeze(pmat_postR);
    elseif r<=16 | r>=26 & r<=27%first shock after barrier
        if r<=16
            rf=7;
        else
            rf=16;
        end
        Rmatrix_place_shock1_pre(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_preR);
        Rmatrix_place_shock1_post(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR);
    elseif r<=25 %scop + first shock
        rf=16;
        Rmatrix_place_scopshock_pre(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_preR);
        Rmatrix_place_scopshock_post(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR);
    else %scop alone
        rf=27;
        Rmatrix_place_scop_pre(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_preR);
        Rmatrix_place_scop_post(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR);
    end
    
    %-------------------- plot sorted heatmaps
 
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
%     subplot_tight(4,7,[1 8 ]); imagesc(1000*hmapLR_pre(:,2:end)); title(['pre ' num2str(cellmat_shockcols(r,precol))]); set(gca,'XTick',[]); caxis([0 10]);; colormap hot;
%     subplot_tight(4,7,[1 8 ]+1); imagesc(1000*hmapLR_shockpre(:,2:end)); title(['shkpre ' num2str(cellmat_shockcols(r,shkcol))]); set(gca,'XTick',[]); caxis([0 10]);;
    
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
%     subplot_tight(4,7,[1 8 ]+14); imagesc(1000*hmapRL_pre(:,2:end)); title(['pre ' num2str(cellmat_shockcols(r,precol))]); set(gca,'XTick',[]); caxis([0 10]);;
%     subplot_tight(4,7,[1 8 ]+1+14); imagesc(1000*hmapRL_shockpre(:,2:end)); title(['shkpre ' num2str(cellmat_shockcols(r,shkcol))]); set(gca,'XTick',[]); caxis([0 10]);;
    
    
    isplace_allpre = (predata.place_scoreLR>minplacescore & predata.spikes_per_LRbee>minspikes) | (predata.place_scoreRL>minplacescore & predata.spikes_per_RLbee>minspikes);
    isplace_allshk = (trndata.place_scoreLR>minplacescore & trndata.spikes_per_LRbee>minspikes) | (trndata.place_scoreRL>minplacescore & trndata.spikes_per_RLbee>minspikes);
    
    preplace_recur=length(isplace_allshk & recur_shockpre_shockflag)/(sum(isplace_allpre)+sum(isplace_allshk));
    
    %--------------------------------------
    
    prefs=nanmean([sum([predata.dcurve_LR(predata.isplaceLR_N,:)>repmat(max(predata.dcurve_LR(predata.isplaceLR_N,:)')',1,23)*.5]') sum([predata.dcurve_RL(predata.isplaceRL_N,:)>repmat(max(predata.dcurve_RL(predata.isplaceRL_N,:)')',1,23)*.5]')]*10.8);
    shkfs=nanmean([sum([trndata.dcurve_LR(trndata.isplaceLR_N,:)>repmat(max(trndata.dcurve_LR(trndata.isplaceLR_N,:)')',1,23)*.5]') sum([trndata.dcurve_RL(trndata.isplaceRL_N,:)>repmat(max(trndata.dcurve_RL(trndata.isplaceRL_N,:)')',1,23)*.5]')]*10.8);
    postfs=nanmean([sum([trndata.dcurve_LR_part2(trndata.isplaceLR_N,:)>repmat(max(trndata.dcurve_LR_part2(trndata.isplaceLR_N,:)')',1,23)*.5]') sum([trndata.dcurve_RL_part2(trndata.isplaceRL_N,:)>repmat(max(trndata.dcurve_RL_part2(trndata.isplaceRL_N,:)')',1,23)*.5]')]*10.8);
    
    temp1=[predata.dcurve_LR(predata.isplaceLR_N,:).*[predata.dcurve_LR(predata.isplaceLR_N,:)<repmat(max(predata.dcurve_LR(predata.isplaceLR_N,:)')',1,23)*.5]]'; temp1(find(temp1==0))=NaN;
    temp2=[predata.dcurve_RL(predata.isplaceRL_N,:).*[predata.dcurve_RL(predata.isplaceRL_N,:)<repmat(max(predata.dcurve_RL(predata.isplaceRL_N,:)')',1,23)*.5]]'; temp2(find(temp2==0))=NaN;
    preblrate=1000*nanmean([nanmean(temp1) nanmean(temp2)]);
    prebldataLR=1000*nanmean(temp1);
    prebldataRL=1000*nanmean(temp2);
    temp1=[trndata.dcurve_LR(trndata.isplaceLR_N,:).*[trndata.dcurve_LR(trndata.isplaceLR_N,:)<repmat(max(trndata.dcurve_LR(trndata.isplaceLR_N,:)')',1,23)*.5]]'; temp1(find(temp1==0))=NaN;
    temp2=[trndata.dcurve_RL(trndata.isplaceRL_N,:).*[trndata.dcurve_RL(trndata.isplaceRL_N,:)<repmat(max(trndata.dcurve_RL(trndata.isplaceRL_N,:)')',1,23)*.5]]'; temp2(find(temp2==0))=NaN;
    shkblrate=1000*nanmean([nanmean(temp1) nanmean(temp2)]);
    shkbldata=1000*([nanmean(temp1) nanmean(temp2)]);
    shkbldataLR=1000*nanmean(temp1);
    shkbldataRL=1000*nanmean(temp2);
    temp1=[trndata.dcurve_LR_part2(trndata.isplaceLR_N,:).*[trndata.dcurve_LR_part2(trndata.isplaceLR_N,:)<repmat(max(trndata.dcurve_LR_part2(trndata.isplaceLR_N,:)')',1,23)*.5]]'; temp1(find(temp1==0))=NaN;
    temp2=[trndata.dcurve_RL_part2(trndata.isplaceLR_N,:).*[trndata.dcurve_RL_part2(trndata.isplaceLR_N,:)<repmat(max(trndata.dcurve_RL_part2(trndata.isplaceLR_N,:)')',1,23)*.5]]'; temp2(find(temp2==0))=NaN;
    
    mp1=floor(size(trndata.LRint,1)/2); mp2=size(trndata.LRint,1)-mp1;
    mp3=floor(size(trndata.RLint,1)/2); mp4=size(trndata.RLint,1)-mp3;
    
    analysis_Nbee(r,:)=[ size(predata.LRint,1) mp1 mp2 size(predata.RLint,1) mp3 mp4];
    
    prespd=250/nanmedian([predata.LRint(:,2)-predata.LRint(:,1); predata.RLint(:,2)-predata.RLint(:,1)]/11.4);
    shkspd=250/nanmedian([trndata.LRint(1:mp1,2)-trndata.LRint(1:mp1,1); trndata.RLint(1:mp1,2)-trndata.RLint(1:mp1,1)]/11.4);
    postspd=250/nanmedian([trndata.LRint(mp2:end,2)-trndata.LRint(mp2:end,1); trndata.RLint(mp2:end,2)-trndata.RLint(mp2:end,1)]/11.4);
    
    preocc=median([preframe.vmap_LR preframe.vmap_RL])/1000;
    shkocc=median([trndata.vmap_LR_part1 trndata.vmap_RL_part1])/1000;
    postocc=median([trndata.vmap_LR_part2 trndata.vmap_RL_part2])/1000;
    
    predata.isplace=predata.isplaceLR_N | predata.isplaceRL_N;
    trndata.isplace=trndata.isplaceLR_N | trndata.isplaceRL_N;
    postdata.isplace=trndata.isplaceLR_N | trndata.isplaceRL_N;
    
    
    if ismember(r,[1:4 6 13:16 26:27])
        
        
        analysis_shift_dfpre=[analysis_shift_dfpre; [shiftdistLR_pre'; shiftdistRL_pre']];
        
        analysis_globalR_dfpre=[analysis_globalR_dfpre; globalR_pre];
        
        dfpre_placecellmap=[dfpre_placecellmap; [predata.dcurve_LR(predata.isplace,:) predata.dcurve_RL(predata.isplace,:)]];
        dfshk_placecellmap=[dfshk_placecellmap; [trndata.dcurve_LR(trndata.isplace,:) trndata.dcurve_RL(trndata.isplace,:)]];
        
        hm_pre_shk_shortpath=[hm_pre_shk_shortpath; shortpath_preR];
        hm_post_shk_shortpath=[hm_post_shk_shortpath; shortpath_postR];
        
        
    end
    if ismember(r,[17:25])
        analysis_shift_scpre=[analysis_shift_scpre; [shiftdistLR_pre'; shiftdistRL_pre']];
        
        analysis_globalR_scpre=[analysis_globalR_scpre; globalR_pre];
        
        scpre_placecellmap=[scpre_placecellmap; [predata.dcurve_LR(predata.isplace,:) predata.dcurve_RL(predata.isplace,:)]];
        scshk_placecellmap=[scshk_placecellmap; [trndata.dcurve_LR(trndata.isplace,:) trndata.dcurve_RL(trndata.isplace,:)]];
        
        hm_pre_scp_shortpath=[hm_pre_scp_shortpath; shortpath_preR];
        hm_post_scp_shortpath=[hm_post_scp_shortpath; shortpath_postR];
    end
    if ismember(r,[7:12])
        analysis_shift_barpre=[analysis_shift_barpre; [shiftdistLR_pre'; shiftdistRL_pre']];
        
        analysis_globalR_barpre=[analysis_globalR_barpre; globalR_pre];
        
        barpre_placecellmap=[barpre_placecellmap; [predata.dcurve_LR(predata.isplace,:) predata.dcurve_RL(predata.isplace,:)]];
        barshk_placecellmap=[barshk_placecellmap; [trndata.dcurve_LR(trndata.isplace,:) trndata.dcurve_RL(trndata.isplace,:)]];
        
        hm_pre_bar_shortpath=[hm_pre_bar_shortpath; shortpath_preR];
        hm_post_bar_shortpath=[hm_post_bar_shortpath; shortpath_postR];
    end
    
    analysis_results(r,:)=[r sum(predata.isplaceLR_N | predata.isplaceRL_N)/length(predata.isplaceLR_N) sum(trndata.isplaceLR_N | trndata.isplaceRL_N)/length(trndata.isplaceLR_N) 0 ...
        sum(predata.isplaceLR_N | predata.isplaceRL_N)  sum(trndata.isplaceLR_N | trndata.isplaceRL_N)  0 ... %5 6 7
        nanmedian([predata.dcurve_LR_bps(predata.isplaceLR_N) predata.dcurve_RL_bps(predata.isplaceRL_N)]) nanmedian([trndata.dcurve_LR_bps(trndata.isplaceLR_N) trndata.dcurve_RL_bps(trndata.isplaceRL_N)]) 0 ... %8 9 10
        0 ... % 11
        preblrate shkblrate 0 ...  % 12 13 14
        nanmedian([max(predata.dcurve_LR(predata.isplaceLR_N,:)') max(predata.dcurve_RL(predata.isplaceRL_N,:)')]) nanmedian([max(trndata.dcurve_LR(trndata.isplaceLR_N,:)') max(trndata.dcurve_RL(trndata.isplaceRL_N,:)')]) 0 ... %15 16 17
        0 0 0 ... %18 19 20
        NaN NaN NaN ... %21 22 23
        nanmedian(globalR_pre) 0  ... %24 25
        prefs shkfs postfs ... %26 27 28
        prespd shkspd postspd ...  %29 30 31
        NaN NaN ...
        nanmedian([shiftdistLR_pre'; shiftdistRL_pre']) 0]; %32 33
    
    analysis_groupLR_pre=[analysis_groupLR_pre; ...
        max(predata.dcurve_LR(predata.isplaceLR_N,:)')' ... %peak
        sum([predata.dcurve_LR(predata.isplaceLR_N,:)>repmat(max(predata.dcurve_LR(predata.isplaceLR_N,:)')',1,23)*.5]')'*10.8  ... %prewid
        predata.dcurve_LR_bps(predata.isplaceLR_N)'  ... %info
        prebldataLR' ...
        predata.isrecurringplace(predata.isplaceLR_N)' ...
        prebldataLR'*0+r];
    
    analysis_groupLR_trn=[analysis_groupLR_trn; ...
        max(trndata.dcurve_LR(trndata.isplaceLR_N,:)')' ... %peak
        sum([trndata.dcurve_LR(trndata.isplaceLR_N,:)>repmat(max(trndata.dcurve_LR(trndata.isplaceLR_N,:)')',1,23)*.5]')'*10.8  ... %prewid
        trndata.dcurve_LR_bps(trndata.isplaceLR_N)'  ... %info
        shkbldataLR' ...
        trndata.isrecurringplace(trndata.isplaceLR_N)' ...
        shkbldataLR'*0+r];
    
    analysis_groupRL_pre=[analysis_groupRL_pre; ...
        max(predata.dcurve_RL(predata.isplaceRL_N,:)')' ... %peak
        sum([predata.dcurve_RL(predata.isplaceRL_N,:)>repmat(max(predata.dcurve_RL(predata.isplaceRL_N,:)')',1,23)*.5]')'*10.8  ... %prewid
        predata.dcurve_RL_bps(predata.isplaceRL_N)'  ... %info
        prebldataRL' ...
        predata.isrecurringplace(predata.isplaceRL_N)' ...
        prebldataRL'*0+r];
    
    analysis_groupRL_trn=[analysis_groupRL_trn; ...
        max(trndata.dcurve_RL(trndata.isplaceRL_N,:)')' ... %peak
        sum([trndata.dcurve_RL(trndata.isplaceRL_N,:)>repmat(max(trndata.dcurve_RL(trndata.isplaceRL_N,:)')',1,23)*.5]')'*10.8  ... %prewid
        trndata.dcurve_RL_bps(trndata.isplaceRL_N)'  ... %info
        shkbldataRL' ...
        trndata.isrecurringplace(trndata.isplaceRL_N)' ...
        shkbldataRL'*0+r];
    
    clear temp*;
    
    
end %rat loop -----------------------------------------------------------------------------------------------------------------




