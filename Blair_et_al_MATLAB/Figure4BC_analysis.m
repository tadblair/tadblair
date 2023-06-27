%pre48 vs post, filters for shock responses and flips RL tuning curves

clear;

rats_to_analyze=[1:4 6:27];;%[1:4 6 13:14 16:17 20:25 27];
shockflag=-1; %set to 1 for fig 4B, -1 for fig 4C

placemap_shock=[];
placemap_barrier=[];
placemap_scopshock=[];

load 'Behav_effects';
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

analysis_globalR_dfpre=[];
analysis_globalR_scpre=[];
analysis_globalR_barpre=[];
analysis_globalR_dfpost=[];
analysis_globalR_scpost=[];
analysis_globalR_barpost=[];

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
barpre_placecellmap=[];
barshk_placecellmap=[];
barpost_placecellmap=[];

finaldir='D:\MiniScopeData\Hipp_LFOV - Copy\'; %folder where the goodies are stored

rootdir='D:\MiniScopeData\pvAnalysis\'; %folder where the goodies are stored
lineardir='D:\MiniScopeData\Shock_LFOV\'; %folder where the goodies are stored

%load in times of shocks
load([rootdir 'shocktimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
load([rootdir 'shocktimes2_13']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
load([rootdir 'scopshocktimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
load([rootdir 'barriertimes']); %row order: Hipp6, Hipp7, Hipp8, Hipp9, Hipp12, Hipp13, Hipp15, Hipp18, Hipp30, Hipp31, Hipp32, Hipp34-37
% shocktimes=[shocktimes(1:4,:); shocktimes([7 10],:); shocktimes([5:6 13],:)]; %reorder to match other structures
% shocktimes(5,:)=shocktimes(5,:)+.6;
% shocktimes(7,1)=shocktimes(7,1)+.7;

%number of frames before and after shock for PETH analysis window of shock response
pethwid = 22; %total width of PETH window is 2*pethwid+1

%parameters for binning of spatial position data
spbinsize=13.33;                  % width of each spatial bin in cm
bincenters=-(spbinsize*11):spbinsize:(spbinsize*11);
binedges=[-500 [-(10*spbinsize):spbinsize:(spbinsize*11)]-spbinsize/2 500];

%select output representation scheme for the decoder
%MdlSpec.CodingMatrix='onevsall'; %decoder output vector is sparsecode of each position bin
MdlSpec.CodingMatrix='ordinal'; %decoder output vector is ordinal line representatin of position

%select learning algorithm for the decoder
%MdlSpec.LearnerType = 'SVM'; %SVM network
MdlSpec.LearnerType = 'linear'; %linear regression
MdlSpec.Normalize_Weights = false;

%matrix of column numbers to analyze in each rat's matching matrix
%col1 is pre shock session (day -6), col2 is shock session (day 0), col3 is pos shock session (day +6)
% for Hipp15 (row 5), col2 is the sess9 (last session before the immediate shock session)

 cellmat_shockcols = [...
    %    3 6 9;...  %Hipp6 3 day
    6 5 8;...  %Hipp6 ext
    7 6 10;... %Hipp7
    4 3 7;...  %Hipp8
    5 4 8;...  %Hipp9
    9 8 12;... %Hipp15
    6 5 9;...  %Hipp31
    
    6 5 9;... %Hipp12 (barrier)
    6 5 9;... %Hipp13 (barrier)
    %     5 7 9;... %Hipp18 (barrier) 2 day
    7 6 9;... %Hipp18 (barrier) ext
    8 7 11;... %Hipp35 (barrier)
    %     3 6 9;... %Hipp30 (barrier) 3 dayRmatrix
    6 5 9;... %Hipp30 (barrier) ext
    8 7 11;... %Hipp34 (barrier)
    
    10 9 13;... %Hipp12 (shock)
    10 9 13;... %Hipp13 (shock)
    %    9 13 17;... %Hipp18 (shock) 4 day
    13 12 14;... %Hipp18 (shock) ext
    12 11 15;...   %Hipp35 (shock)
    
    %    9 10 11;... %Hipp30 (shock+scop)
    10 9 13;... %Hipp30 (shock+scop)
    11 10 14;... %Hipp32 (shock+scop)
    12 11 15;... %Hipp34 (shock+scop)
    12 11 15;...   %Hipp36 (shock+scop)
    
    23 22 26;... %Hipp12 (shock+scop)
    23 22 26;... %Hipp13 (shock+scop)
    23 22 26;... %Hipp15 (shock+scop)
    23 22 26;...   %Hipp18 (shock+scop)
        %22 23 24;...   %Hipp18 (shock+scop)
    %17 20 23;...   %Hipp31 (shock+scop)
        20 19 23;...   %Hipp31 (shock+scop)
    
     14 13 15;... %Hipp30 (drug free shock2)
     20 19 24;... %Hipp32 (drug free shock2)
    
    15 14 16;... %Hipp31 (scopolamine alone)
    6 5 7;... %Hipp32 (scopolamine alone)
    7 6 8]; %Hipp36 (scopolamine alone)

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


PC_segmentR=[];
AC_segmentR=[];
PC_segmentP=[];
AC_segmentP=[];
PC_segment_preR=[];
AC_segment_preR=[];
PC_segment_preP=[];
AC_segment_preP=[];
PC_segment_postR=[];
AC_segment_postR=[];
PC_segment_postP=[];
AC_segment_postP=[];

analysis_results=[];

for r=rats_to_analyze %rats_to_analyze %--- loop through rats
    
    r
    
    %need to keep some variables when new rat analysis starts
    clearvars -except hm_* plotbeh resubsample placemap_* allshockinh* *_placecellmap event_place_cells* shockflag global_p Effects rundir *_zone M_wid maxhighbins goodpeakbins minplacescore minspikes *_segment* Rmatrix* Pmatrix* barriertimes allshockcells* allnonshockcells* SHKexc recur_shock* predata postdata shockdata precol postcol analyze_* data otherdata smear convkernel Decode* TCcorr* Loss* chi* MdlSpec rat* session* analysis* r goodsess s *dir *_result *shocktimes* cellmat* heatmaps* pethwid all_chmat binedges;
    
    clear sessionNums;
    
    %read in matching matrix for this rat
    if r<=6 | (r>=13 & r<=16) %drug-free shocks
        
        load([rootdir 'cellmaps\' rat(r).name '_shock_cmap']);
        if r==5 %for Hipp15 only include cells recorded on pre shock AND shock day
            cmap(find(cmap(:,4+1)==0),4)=0;
        end
        
    elseif r>=17 & r<=25 %scopolamine shocks
        
        load([rootdir 'cellmaps\' rat(r).name '_scopshock_cmap']);
        
    elseif (r>=7 & r<=12) %barrier
        
        load([rootdir 'cellmaps\' rat(r).name '_barrier_cmap']);
        cmap2=cmap; sessionNums2=sessionNums;
        load([rootdir 'cellmaps\' rat(r).name '_barrier_cmap_dfdf']);
                
    else %second drug free shock
        
        load([rootdir 'cellmaps\' rat(r).name '_shock2_cmap']);
        
    end
    
    %cmap column indices for the sessions of interest
    pre_column=find(sessionNums==cellmat_shockcols(r,1));
    shk_column=find(sessionNums==cellmat_shockcols(r,2));
     %  shk_column2=find(sessionNums2==cellmat_shockcols2(r,2));
    post_column=find(sessionNums==cellmat_shockcols(r,3));
    
    cmap_df=[cmap(:,shk_column) cmap(:,post_column)];
    save([rootdir 'cellmaps_df\' rat(r).name '_cmap_df'],'cmap_df');
    
    precol=1; %pre session
    if resubsample
        load([finaldir '\' rat(r).name '_linear' num2str(cellmat_shockcols(r,precol)) '_data'],'data');
    else
        load([finaldir '\' rat(r).name '_linear' num2str(cellmat_shockcols(r,precol)) '_normbee'],'data');
    end
    data.S=data.deconv>0;
    fired_on_short=sum([data.S & repmat(data.y'<20,size(data.S,1),1)]');
    predata=data;
    [temp, predata.posbin]=histc(data.x,binedges);
    
    postcol=3;
    if resubsample
        load([finaldir '\' rat(r).name '_linear' num2str(cellmat_shockcols(r,postcol)) '_data'],'data');
    else
        load([finaldir '\' rat(r).name '_linear' num2str(cellmat_shockcols(r,postcol)) '_normbee2'],'data');
        hm_file_post{r}=[rat(r).name '_linear' num2str(cellmat_shockcols(r,postcol)) '_normbee2'];
    end
    data.S=data.deconv>0;
    fired_on_short=sum([data.S & repmat(data.y'<20,size(data.S,1),1)]');
    postdata=data;
    [temp, postdata.posbin]=histc(data.x,binedges);
    
    %load shock session data
    shkcol=2;
    if resubsample
        load([finaldir '\' rat(r).name '_linear' num2str(cellmat_shockcols(r,shkcol)) '_data'],'data');
    else
        load([finaldir '\' rat(r).name '_linear' num2str(cellmat_shockcols(r,shkcol)) '_normbee2'],'data');
        hm_file_pre48{r}=[rat(r).name '_linear' num2str(cellmat_shockcols(r,shkcol)) '_normbee2'];
    end
    
    
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
    
    
    if plotbeh
        
        figure(800+r); clf;
        
        subplot(1,3,1); hold off;
        fnum=1:length(predata.x);
        plot(predata.x,fnum,'k'); hold on;
        
        subplot(1,3,2); hold off;
        fnum=1:length(data.x);
        plot(data.x,fnum,'k'); hold on;
        
        subplot(1,3,3); hold off;
        fnum=1:length(postdata.x);
        plot(postdata.x,fnum,'k'); hold on;
        
    end
            
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
            %                 for i=1:length(shk_peak_valLR) %remove cells that dont fire during beelines in shock session
            %                     if (shk_peak_valLR(i)==0 & shk_peak_valRL(i)==0) | 1000*shk_peak_valLR(i)>11.3 | 1000*shk_peak_valRL(i)>11.3
            %                         matdex=find(cmap(:,shk_column)==i);
            %                         cmap(matdex,shk_column)=0;
            %                     end
            %                 end

            
            [similarityR_preLR resultP] = pv_heatmap_tuningcorr(squeeze(predata.dcurve_LR_even),squeeze(predata.dcurve_LR_odd));
            [similarityR_preRL resultP] = pv_heatmap_tuningcorr(squeeze(predata.dcurve_RL_even),squeeze(predata.dcurve_RL_odd));
            [similarityR_shkLR resultP] = pv_heatmap_tuningcorr(squeeze(data.dcurve_LR_even),squeeze(data.dcurve_LR_odd));
            [similarityR_shkRL resultP] = pv_heatmap_tuningcorr(squeeze(data.dcurve_RL_even),squeeze(data.dcurve_RL_odd));
            [similarityR_postLR resultP] = pv_heatmap_tuningcorr(squeeze(postdata.dcurve_LR_even),squeeze(postdata.dcurve_LR_odd));
            [similarityR_postRL resultP] = pv_heatmap_tuningcorr(squeeze(postdata.dcurve_RL_even),squeeze(postdata.dcurve_RL_odd));
                       
            [pvcR_pre(1) R resultP] = pv_heatmap_cormatrix_pdR(squeeze(predata.dcurve_LR_even(predata.isplaceLR,:)),squeeze(predata.dcurve_LR_odd(predata.isplaceLR,:)));
            [pvcR_pre(2) R resultP] = pv_heatmap_cormatrix_pdR(squeeze(predata.dcurve_RL_even(predata.isplaceRL,:)),squeeze(predata.dcurve_RL_odd(predata.isplaceRL,:)));
            [pvcR_shk(1) R resultP] = pv_heatmap_cormatrix_pdR(squeeze(data.dcurve_LR_even(data.isplaceLR,:)),squeeze(data.dcurve_LR_odd(data.isplaceLR,:)));
            [pvcR_shk(2) R resultP] = pv_heatmap_cormatrix_pdR(squeeze(data.dcurve_RL_even(data.isplaceRL,:)),squeeze(data.dcurve_RL_odd(data.isplaceRL,:)));
            [pvcR_post(1) R resultP] = pv_heatmap_cormatrix_pdR(squeeze(postdata.dcurve_LR_even(postdata.isplaceLR,:)),squeeze(postdata.dcurve_LR_odd(postdata.isplaceLR,:)));
            [pvcR_post(2) R resultP] = pv_heatmap_cormatrix_pdR(squeeze(postdata.dcurve_RL_even(postdata.isplaceRL,:)),squeeze(postdata.dcurve_RL_odd(postdata.isplaceRL,:)));
            
            
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
            %                 preflags.includedLR=preflags.placeLR & 0;
            %                 preflags.includedRL=preflags.placeLR & 0;
            %                 shkflags.included_preLR=shkflags.placeLR & 0;
            %                 shkflags.included_preRL=shkflags.placeLR & 0;
            %                 postflags.includedLR=postflags.placeLR & 0;
            %                 postflags.includedRL=postflags.placeLR & 0;
            %                 shkflags.included_postLR=shkflags.placeLR & 0;
            %                 shkflags.included_postRL=shkflags.placeLR & 0;
            
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


            % NOT SURE THIS IS USEFUL ANYMORE
            
            % clear *matrow;
            % predex=find(predata.isplace_N);
            % for i=predex(:)'
            %     if ~isempty(find(cmap(:,pre_column)==i))
            %     preplacematrow(i)=find(cmap(:,pre_column)==i);
            %     end
            % end
            % shkdex=find(data.isplace_N);
            % for i=shkdex(:)'
            %     if ~isempty(find(cmap(:,shk_column)==i))
            %     shkplacematrow(i)=find(cmap(:,shk_column)==i);
            %     end
            % end
            % postdex=find(postdata.isplace_N);
            % for i=postdex(:)'
            %     if ~isempty(find(cmap(:,post_column)==i))
            %     postplacematrow(i)=find(cmap(:,post_column)==i);
            %     end
            % end
            
        analysis_meaninfo(r,:) = [nanmean([predata.infoLR(predata.isplaceLR(:)') predata.infoRL(predata.isplaceRL(:)')]) nanmean([data.infoLR(data.isplaceLR(:)') data.infoRL(data.isplaceRL(:)')]) nanmean([postdata.infoLR(postdata.isplaceLR(:)') postdata.infoRL(postdata.isplaceRL(:)')])];
        analysis_meanspeed(r,:) = [nanmean(predata.subspd_median) nanmean(data.subspd_median) nanmean(postdata.subspd_median)];
        analysis_meanplaceperc(r,:) = [mean(sum(predata.isplace_N))/predata.numgoodcells mean(sum(data.isplace_N))/data.numgoodcells mean(sum(postdata.isplace_N))/postdata.numgoodcells];
        analysis_meanfsize(r,:) = 12*[nanmean([prehighbinsLR(predata.isplaceLR(:)') prehighbinsRL(predata.isplaceRL(:)')]) nanmean([shkhighbinsLR(data.isplaceLR(:)') shkhighbinsRL(data.isplaceRL(:)')]) nanmean([posthighbinsLR(postdata.isplaceLR(:)') posthighbinsRL(postdata.isplaceRL(:)')])];
        analysis_meanmaxrate(r,:) = [nanmean([premaxrateLR(predata.isplaceLR(:)') premaxrateRL(predata.isplaceRL(:)')]) nanmean([shkmaxrateLR(data.isplaceLR(:)') shkmaxrateRL(data.isplaceRL(:)')]) nanmean([postmaxrateLR(postdata.isplaceLR(:)') postmaxrateRL(postdata.isplaceRL(:)')])];
        analysis_similarityR(r,:) = [nanmean([similarityR_preLR(predata.isplaceLR(:)') similarityR_preRL(predata.isplaceRL(:)')]) nanmean([similarityR_shkLR(data.isplaceLR(:)') similarityR_shkRL(data.isplaceRL(:)')]) nanmean([similarityR_postLR(postdata.isplaceLR(:)') similarityR_postRL(postdata.isplaceRL(:)')])];
        analysis_pvcR(r,:) = [nanmean(pvcR_pre(:)) nanmean(pvcR_shk(:)) nanmean(pvcR_post(:))];
        
%         analysis_results(r,:)=[0 analysis_meanplaceperc(r,:) ...
%             sum(predata.isplaceLR | predata.isplaceRL)  sum(data.isplaceLR | data.isplaceRL)  sum(postdata.isplaceLR | postdata.isplaceRL) ... %5 6 7
%             analysis_meaninfo(r,:)... %8 9 10
%             0 ... % 11
%             length(predata.LRint)+length(predata.RLint) length(data.LRint)+length(data.RLint) length(postdata.LRint)+length(postdata.RLint) ...  % 12 13 14
%             analysis_meanmaxrate(r,:) ... %15 16 17
%             0 0 0 ... %18 19 20
%             0 0 0 ... %21 22 23
%             0 0 ... %24 25
%             analysis_meanfsize(r,:) ... %26 27 28
%             analysis_meanspeed(r,:) ... %29-31
%             analysis_similarityR(r,:) ... %32-34
%             analysis_pvcR(r,:) ... %35-37
%             predata.numgoodcells data.numgoodcells postdata.numgoodcells]; %38-40
                
    
    if plotbeh
        
        subplot(1,3,1); hold on;
        fnum=1:length(predata.x);
        for i=1:size(predata.LRint)
            plot(predata.x(predata.LRint(i,1):predata.LRint(i,2)),fnum(predata.LRint(i,1):predata.LRint(i,2)),'r');
        end
        for i=1:size(predata.RLint)
            plot(predata.x(predata.RLint(i,1):predata.RLint(i,2)),fnum(predata.RLint(i,1):predata.RLint(i,2)),'b');
        end
        title(nanmean([ predata.LRspdlist(:); predata.RLspdlist(:)] ));
        subplot(1,3,2); hold on;
        fnum=1:length(data.x);
        for i=1:size(data.LRint)
            plot(data.x(data.LRint(i,1):data.LRint(i,2)),fnum(data.LRint(i,1):data.LRint(i,2)),'r');
        end
        for i=1:size(data.RLint)
            plot(data.x(data.RLint(i,1):data.RLint(i,2)),fnum(data.RLint(i,1):data.RLint(i,2)),'b');
        end
        title(nanmean([ data.LRspdlist(:); data.RLspdlist(:)] ));
        subplot(1,3,3); hold on;
        fnum=1:length(postdata.x);
        for i=1:size(postdata.LRint)
            plot(postdata.x(postdata.LRint(i,1):postdata.LRint(i,2)),fnum(postdata.LRint(i,1):postdata.LRint(i,2)),'r');
        end
        for i=1:size(postdata.RLint)
            plot(postdata.x(postdata.RLint(i,1):postdata.RLint(i,2)),fnum(postdata.RLint(i,1):postdata.RLint(i,2)),'b');
        end
        title(nanmean([ postdata.LRspdlist(:); postdata.RLspdlist(:)] ));
        
        drawnow;
        
    end
    
    
clear *matrow *matrowLR *matrowRL;
    preplacematrowLR=[];
    preplacematrowRL=[];
    shkplacematrowLR=[];
    shkplacematrowRL=[];
    postplacematrowLR=[];
    postplacematrowRL=[];
%predex=find(predata.isplace_N);
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
%postdex=find(postdata.isplace_N);
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
                        recur_shockpre_shockflag(~SHKexc)=false;
                        for i=1:length(recur_shockpre_matdex)
                            if ~SHKexc(cmap(recur_shockpre_matdex(i),shk_column))
                                recur_shockpre_matdex(i)=0;
                            end
                        end
                    end
                    if shockflag==-1
                        recur_shockpre_shockflag(SHKexc)=false;
                        for i=1:length(recur_shockpre_matdex)
                            if SHKexc(cmap(recur_shockpre_matdex(i),shk_column))
                                recur_shockpre_matdex(i)=0;
                            end
                        end
                    end
                end
end
                recur_shockpre_matdex=recur_shockpre_matdex(find(recur_shockpre_matdex>0)); %drop non recurring cells from the list
                
                for d=1:length(recur_shockpre_matdex) %loop through recurring cells
                    hmap_pre(d,:)=[predata.dcurve_LR(cmap(recur_shockpre_matdex(d),pre_column),:) predata.dcurve_RL(cmap(recur_shockpre_matdex(d),pre_column),23:-1:1)]; 
                    hmap_shockpre(d,:)=[data.dcurve_LR(cmap(recur_shockpre_matdex(d),shk_column),:) data.dcurve_RL(cmap(recur_shockpre_matdex(d),shk_column),23:-1:1)]; 
                        prehighbinsLR=sum(predata.dcurve_LR(cmap(recur_shockpre_matdex(d),pre_column),:)>max(predata.dcurve_LR(cmap(recur_shockpre_matdex(d),pre_column),:))*.7);
                        shkprehighbinsLR=sum(data.dcurve_LR(cmap(recur_shockpre_matdex(d),shk_column),:)>max(data.dcurve_LR(cmap(recur_shockpre_matdex(d),shk_column),:))*.7);
                        prehighbinsRL=sum(predata.dcurve_RL(cmap(recur_shockpre_matdex(d),pre_column),:)>max(predata.dcurve_RL(cmap(recur_shockpre_matdex(d),pre_column),:))*.7);
                        shkprehighbinsRL=sum(data.dcurve_RL(cmap(recur_shockpre_matdex(d),shk_column),:)>max(data.dcurve_RL(cmap(recur_shockpre_matdex(d),shk_column),:))*.7);
                    isplace_preLR(d) = (prehighbinsLR<maxhighbins & predata.place_scoreLR(cmap(recur_shockpre_matdex(d),pre_column))>minplacescore & predata.spikes_per_LRbee(cmap(recur_shockpre_matdex(d),pre_column))>minspikes )...
                                     | (shkprehighbinsLR<maxhighbins & data.place_scoreLR(cmap(recur_shockpre_matdex(d),shk_column))>minplacescore & data.spikes_per_LRbee(cmap(recur_shockpre_matdex(d),shk_column))>minspikes );
                    isplace_preRL(d) = (prehighbinsRL<maxhighbins & predata.place_scoreRL(cmap(recur_shockpre_matdex(d),pre_column))>minplacescore & predata.spikes_per_RLbee(cmap(recur_shockpre_matdex(d),pre_column))>minspikes )...
                                     | (shkprehighbinsRL<maxhighbins & data.place_scoreRL(cmap(recur_shockpre_matdex(d),shk_column))>minplacescore & data.spikes_per_RLbee(cmap(recur_shockpre_matdex(d),shk_column))>minspikes );
%                    isplace_preRL(d) = predata.isplaceRL(cmap(recur_shockpre_matdex(d),pre_column)) | data.isplaceRL(cmap(recur_shockpre_matdex(d),shk_column));
                    isplace_pre(d) = (isplace_preLR(d) | isplace_preRL(d));
                    nmax=max([hmap_pre(d,:) hmap_shockpre(d,:)]); 
                    
                    [temp, peakbin_preLR(d)]=max(predata.dcurve_LR(d,:));
                    [temp, peakbin_preRL(d)]=max(predata.dcurve_RL(d,:));
                    [temp, peakbin_shockpreLR(d)]=max(data.dcurve_LR(d,:));
                    [temp, peakbin_shockpreRL(d)]=max(data.dcurve_RL(d,:));

                end
                
                %determine which cells recorded in the shock session also recurred in the post session
                recur_shockpost_matdex=find(cmap(:,post_column)>0 & cmap(:,shk_column)>0  & cmap(:,pre_column)>0); %rows of cell2index map that contain cells recurring in shock & pre sessions
                recur_shockpost_postflag=false(1,size(postdata.deconv,1)); %initialize recurrence flags
                recur_shockpost_postflag(cmap(recur_shockpost_matdex,post_column))=true; %logical vector of cells recorded in SHOCK session, flagged if they recur in the other session
          %      postflags.recur=recur_shockpost_postflag;
                recur_shockpost_shockflag=false(1,size(data.deconv,1)); %initialize recurrence flags
                recur_shockpost_shockflag(cmap(recur_shockpost_matdex,shk_column))=true; %logical vector of cells recorded in SHOCK session, flagged if they recur in the other session
%                 shkflags.recur_post=recur_shockpost_shockflag;

if ~(shockflag==0)                
                if ~(r>6 & r<13) & ~(r==5)
                if shockflag==1
                    recur_shockpost_shockflag(~SHKexc)=false;
                    for i=1:length(recur_shockpost_matdex)
                        if ~SHKexc(cmap(recur_shockpost_matdex(i),shk_column))
                            recur_shockpost_matdex(i)=0;
                        end
                    end
                end
                if shockflag==-1
                    recur_shockpost_shockflag(SHKexc)=false;
                    for i=1:length(recur_shockpost_matdex)
                        if SHKexc(cmap(recur_shockpost_matdex(i),shk_column))
                            recur_shockpost_matdex(i)=0;
                        end
                    end
                end
                end
end
                recur_shockpost_matdex=recur_shockpost_matdex(find(recur_shockpost_matdex>0));
                
                for d=1:length(recur_shockpost_matdex)
                    hmap_post(d,:)=[postdata.dcurve_LR(cmap(recur_shockpost_matdex(d),post_column),:) postdata.dcurve_RL(cmap(recur_shockpost_matdex(d),post_column),23:-1:1)]; 
                    hmap_shockpost(d,:)=[data.dcurve_LR(cmap(recur_shockpost_matdex(d),shk_column),:) data.dcurve_RL(cmap(recur_shockpost_matdex(d),shk_column),23:-1:1)]; 
                        posthighbinsLR=sum(postdata.dcurve_LR(cmap(recur_shockpost_matdex(d),post_column),:)>max(postdata.dcurve_LR(cmap(recur_shockpost_matdex(d),post_column),:))*.7);
                        shkposthighbinsLR=sum(data.dcurve_LR(cmap(recur_shockpost_matdex(d),shk_column),:)>max(data.dcurve_LR(cmap(recur_shockpost_matdex(d),shk_column),:))*.7);
                        posthighbinsRL=sum(postdata.dcurve_RL(cmap(recur_shockpost_matdex(d),post_column),:)>max(postdata.dcurve_RL(cmap(recur_shockpost_matdex(d),post_column),:))*.7);
                        shkposthighbinsRL=sum(data.dcurve_RL(cmap(recur_shockpost_matdex(d),shk_column),:)>max(data.dcurve_RL(cmap(recur_shockpost_matdex(d),shk_column),:))*.7);
                    isplace_postLR(d) = (posthighbinsLR<maxhighbins & postdata.place_scoreLR(cmap(recur_shockpost_matdex(d),post_column))>minplacescore & postdata.spikes_per_LRbee(cmap(recur_shockpost_matdex(d),post_column))>minspikes )...
                                     | (shkposthighbinsLR<maxhighbins & data.place_scoreLR(cmap(recur_shockpost_matdex(d),shk_column))>minplacescore & data.spikes_per_LRbee(cmap(recur_shockpost_matdex(d),shk_column))>minspikes );
                    isplace_postRL(d) = (posthighbinsRL<maxhighbins & postdata.place_scoreRL(cmap(recur_shockpost_matdex(d),post_column))>minplacescore & postdata.spikes_per_RLbee(cmap(recur_shockpost_matdex(d),post_column))>minspikes )...
                                     | (shkposthighbinsRL<maxhighbins & data.place_scoreRL(cmap(recur_shockpost_matdex(d),shk_column))>minplacescore & data.spikes_per_RLbee(cmap(recur_shockpost_matdex(d),shk_column))>minspikes );
%                    isplace_preRL(d) = predata.isplaceRL(cmap(recur_shockpre_matdex(d),pre_column)) | data.isplaceRL(cmap(recur_shockpre_matdex(d),shk_column));
                    isplace_post(d) = (isplace_postLR(d) | isplace_postRL(d));
                    nmax=max([hmap_post(d,:) hmap_shockpost(d,:)]); 
                    
%                     postflags.includedLR(cmap(recur_shockpost_matdex(d),post_column))=isplace_postLR(d);
%                     postflags.includedRL(cmap(recur_shockpost_matdex(d),post_column))=isplace_postRL(d);
%                     shkflags.included_postLR(cmap(recur_shockpost_matdex(d),shk_column))=isplace_postLR(d);
%                     shkflags.included_postRL(cmap(recur_shockpost_matdex(d),shk_column))=isplace_postRL(d);
                    
                    [temp, peakbin_postLR(d)]=max(postdata.dcurve_LR(d,:));
                    [temp, peakbin_postRL(d)]=max(postdata.dcurve_RL(d,:));
                    [temp, peakbin_shockpostLR(d)]=max(data.dcurve_LR(d,:));
                    [temp, peakbin_shockpostRL(d)]=max(data.dcurve_RL(d,:));
                    
                    %session specific place flags
%                     isplace_postLRRL(d) = (posthighbinsLR<maxhighbins & postdata.place_scoreLR(cmap(recur_shockpost_matdex(d),post_column))>minplacescore & ismember(post_peak_LRd(cmap(recur_shockpost_matdex(d),post_column)),goodpeakbins) & postdata.spikes_per_LRbee(cmap(recur_shockpost_matdex(d),post_column))>minspikes )...
%                                      | (posthighbinsRL<maxhighbins & postdata.place_scoreRL(cmap(recur_shockpost_matdex(d),post_column))>minplacescore & ismember(post_peak_RLd(cmap(recur_shockpost_matdex(d),post_column)),goodpeakbins) & postdata.spikes_per_RLbee(cmap(recur_shockpost_matdex(d),post_column))>minspikes );
%                     isplace_shockpostLRRL(d) = (shkposthighbinsLR<maxhighbins & data.place_scoreLR(cmap(recur_shockpost_matdex(d),shk_column))>minplacescore & ismember(shk_peak_LRd(cmap(recur_shockpost_matdex(d),shk_column)),goodpeakbins) & data.spikes_per_LRbee(cmap(recur_shockpost_matdex(d),shk_column))>minspikes ) ...
%                                      | (shkposthighbinsRL<maxhighbins & data.place_scoreRL(cmap(recur_shockpost_matdex(d),shk_column))>minplacescore & ismember(shk_peak_RLd(cmap(recur_shockpost_matdex(d),shk_column)),goodpeakbins) & data.spikes_per_RLbee(cmap(recur_shockpost_matdex(d),shk_column))>minspikes );
                                 
%                     isplace_postLRRL_LR(d) = (posthighbinsLR<maxhighbins & postdata.place_scoreLR(cmap(recur_shockpost_matdex(d),post_column))>minplacescore & ismember(post_peak_LRd(cmap(recur_shockpost_matdex(d),post_column)),goodpeakbins) & postdata.spikes_per_LRbee(cmap(recur_shockpost_matdex(d),post_column))>minspikes );
%                     isplace_postLRRL_RL(d) = (posthighbinsRL<maxhighbins & postdata.place_scoreRL(cmap(recur_shockpost_matdex(d),post_column))>minplacescore & ismember(post_peak_RLd(cmap(recur_shockpost_matdex(d),post_column)),goodpeakbins) & postdata.spikes_per_RLbee(cmap(recur_shockpost_matdex(d),post_column))>minspikes );
%                     isplace_shockpostLRRL_LR(d) = (shkposthighbinsLR<maxhighbins & data.place_scoreLR(cmap(recur_shockpost_matdex(d),shk_column))>minplacescore & ismember(shk_peak_LRd(cmap(recur_shockpost_matdex(d),shk_column)),goodpeakbins) & data.spikes_per_LRbee(cmap(recur_shockpost_matdex(d),shk_column))>minspikes );
%                     isplace_shockpostLRRL_RL(d) = (shkposthighbinsRL<maxhighbins & data.place_scoreRL(cmap(recur_shockpost_matdex(d),shk_column))>minplacescore & ismember(shk_peak_RLd(cmap(recur_shockpost_matdex(d),shk_column)),goodpeakbins) & data.spikes_per_RLbee(cmap(recur_shockpost_matdex(d),shk_column))>minspikes );
%                    hmap_post(d,:)=hmap_post(d,:)/nmax; 
%                    hmap_shockpost(d,:)=hmap_shockpost(d,:)/nmax; 
%                     hmap_post(d,:)=-log10(1+hmap_post(d,:)); 
%                     hmap_shockpost(d,:)=-log10(1+hmap_shockpost(d,:)); 
%                     if isplace_postLR(d)
%                         posthighbinsLR=[posthighbinsLR sum(postdata.deconv(cmap(recur_shockpost_matdex(d),post_column),:)>0)];
%                         shkposthighbinsLR=[shkposthighbinsLR sum(data.deconv(cmap(recur_shockpost_matdex(d),shk_column),:)>0)];
%                     end 
%                     if isplace_postRL(d)
%                         posthighbinsRL=[posthighbinsRL sum(postdata.deconv(cmap(recur_shockpost_matdex(d),post_column),:)>0)];
%                         shkposthighbinsRL=[shkposthighbinsRL sum(data.deconv(cmap(recur_shockpost_matdex(d),shk_column),:)>0)];
%                     end
                end
                
%[sum(isplace_preLR) sum(isplace_preRL); sum(isplace_postLR) sum(isplace_postRL)]
% figure(300); clf; 
% subplot(2,2,1); scatter(prehighbinsLR/size(predata.deconv,2),shkprehighbinsLR/size(data.deconv,2),'.');
% subplot(2,2,2); scatter(prehighbinsRL/size(predata.deconv,2),shkprehighbinsRL/size(data.deconv,2),'.');
% subplot(2,2,3); scatter(posthighbinsLR/size(postdata.deconv,2),shkposthighbinsLR/size(data.deconv,2),'.');
% subplot(2,2,4); scatter(posthighbinsRL/size(postdata.deconv,2),shkposthighbinsRL/size(data.deconv,2),'.');
                if (r<7 & ~(r==5)) | (r>12 & r<17) | (r==26) | (r==27)                     
                    hm_pre_shk_group=[hm_pre_shk_group; hmap_pre(:,12)*0+r];
                    hm_pre_shk=[hm_pre_shk; hmap_pre];
                    hm_shockpre_shk=[hm_shockpre_shk; hmap_shockpre];
                    hm_shockpost_shk=[hm_shockpost_shk; hmap_shockpost];
                    hm_post_shk_group=[hm_post_shk_group; hmap_post(:,12)*0+r];
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
                    hm_pre_scpshk_group=[hm_pre_scpshk_group; hmap_pre(:,12)*0+r];
                    hm_pre_scpshk=[hm_pre_scpshk; hmap_pre];
                    hm_shockpre_scpshk=[hm_shockpre_scpshk; hmap_shockpre];
                    hm_shockpost_scpshk=[hm_shockpost_scpshk; hmap_shockpost];
                    hm_post_scpshk_group=[hm_post_scpshk_group; hmap_post(:,12)*0+r];
                    hm_post_scpshk=[hm_post_scpshk; hmap_post];
                    hm_pflags_pre_scpshk=[hm_pflags_pre_scpshk; [isplace_preLR(:) isplace_preRL(:)] ];
                    hm_pflags_post_scpshk=[hm_pflags_post_scpshk; [isplace_postLR(:) isplace_postRL(:)] ];
                end

                %-------- correlation matrices (do this presort to avoid confusion)
                [allmat_preR, allmat_preP] = pv_heatmap_cormatrix(hmap_pre(:,:),hmap_shockpre(:,:));
                [allmat_postR, allmat_postP] = pv_heatmap_cormatrix(hmap_shockpost(:,:),hmap_post(:,:));

%                 [pmat_preR, pmat_preP] = pv_heatmap_cormatrix_plot([hmap_pre(isplace_preLR,1:23); hmap_pre(isplace_preRL,24:end)],[hmap_shockpre(isplace_preLR,1:23); hmap_shockpre(isplace_preRL,24:end)],500);
%                 [pmat_postR, pmat_postP] = pv_heatmap_cormatrix_plot([hmap_post(isplace_postLR,1:23); hmap_post(isplace_postRL,24:end)],[hmap_shockpost(isplace_postLR,1:23); hmap_shockpost(isplace_postRL,24:end)],501);

                [pcorrLR_preR, pcorrLR_preP] = pv_heatmap_tuningcorr([hmap_pre(isplace_preLR,1:23)],[hmap_shockpre(isplace_preLR,1:23)]);
                [pcorrLR_postR, pcorrLR_postP] = pv_heatmap_tuningcorr([hmap_post(isplace_postLR,1:23)],[hmap_shockpost(isplace_postLR,1:23)]);
                
                [pcorrRL_preR, pcorrRL_preP] = pv_heatmap_tuningcorr([hmap_pre(isplace_preRL,24:end)],[hmap_shockpre(isplace_preRL,24:end)]);
                [pcorrRL_postR, pcorrRL_postP] = pv_heatmap_tuningcorr([hmap_post(isplace_postRL,24:end)],[hmap_shockpost(isplace_postRL,24:end)]);

                figure(1); clf;

                subplot(4,7,28); hold off;
                histogram([pcorrLR_preR pcorrRL_preR],-1:.1:1); hold on;
                histogram([pcorrLR_postR pcorrRL_postR],-1:.1:1); hold on;
                [global_p(r), h]=ranksum([pcorrLR_preR pcorrRL_preR],[pcorrLR_postR pcorrRL_postR]);
                
                analysis_median_tcR(r,:)=[nanmedian([pcorrLR_preR pcorrRL_preR]) nanmedian([pcorrLR_postR pcorrRL_postR])];
                                
                globalR_pre=[pcorrLR_preR'; pcorrRL_preR'];
%                 globalR_pre(find(globalR_pre>.99))=.99;
%                 globalR_pre(find(globalR_pre<-.99))=-.99;
%                 globalR_pre=fisherz(globalR_pre);
                globalR_post=[pcorrLR_postR'; pcorrRL_postR'];
%                 globalR_post(find(globalR_post>.99))=.99;
%                 globalR_post(find(globalR_post<-.99))=-.99;
%                 globalR_post=fisherz(globalR_post);
                
                [pmatLR_preR, pmatLR_preP] = pv_heatmap_cormatrix([hmap_pre(isplace_preLR,1:23)],[hmap_shockpre(isplace_preLR,1:23)]);
                [pmatLR_postR, pmatLR_postP] = pv_heatmap_cormatrix([hmap_post(isplace_postLR,1:23)],[hmap_shockpost(isplace_postLR,1:23)]);

                [pmatRL_preR, pmatRL_preP] = pv_heatmap_cormatrix([hmap_pre(isplace_preRL,24:end)],[hmap_shockpre(isplace_preRL,24:end)]);
                [pmatRL_postR, pmatRL_postP] = pv_heatmap_cormatrix([hmap_post(isplace_postRL,24:end)],[hmap_shockpost(isplace_postRL,24:end)]);
 [sum(isplace_postLR) sum(isplace_postRL)]

                pmat_preR(1,:,:)=pmatLR_preR; pmat_preR(2,:,:)=pmatRL_preR; %pmat_preR=squeeze(nanmean(pmat_preR)); 
                pmat_postR(1,:,:)=pmatLR_postR; pmat_postR(2,:,:)=pmatRL_postR; %pmat_postR=squeeze(nanmean(pmat_postR)); 
                pmat_preP(1,:,:)=pmatLR_preP; pmat_preP(2,:,:)=pmatRL_preP; %pmat_preP=squeeze(nanmean(pmat_preP)); 
                pmat_postP(1,:,:)=pmatLR_postP; pmat_postP(2,:,:)=pmatRL_postP; %pmat_postP=squeeze(nanmean(pmat_postP)); 
                
%                 [allcorrLR_preR, allcorrLR_preP] = pv_heatmap_tuningcorr(hmap_pre(:,1:23),hmap_shockpre(:,1:23));
%                 [allcorrRL_preR, allcorrRL_preP] = pv_heatmap_tuningcorr(hmap_pre(:,24:end),hmap_shockpre(:,24:end));
%                 [allcorrLR_postR, allcorrLR_postP] = pv_heatmap_tuningcorr(hmap_post(:,1:23),hmap_shockpost(:,1:23));
%                 [allcorrRL_postR, allcorrRL_postP] = pv_heatmap_tuningcorr(hmap_post(:,24:end),hmap_shockpost(:,24:end));

%                 all_diffR=allmat_postR-allmat_preR;
%                 all_diffP=allmat_postP-allmat_preP;
                
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
subplot(4,7,7); plot(PC_diagR); hold on; plot(PC_diag_preR); plot(PC_diag_postR); title(rat(r).name); set(gca,'XTick',[]);

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

subplot(4,7,21); plot(PC_diagP); hold on; plot(PC_diag_preP); plot(PC_diag_postP); set(gca,'XTick',[]);

% shkflags.adjM_zoneLR=adjM_zoneLR;
% shkflags.adjM_zoneRL=adjM_zoneRL;


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

%tack on last R values
%             shortpath_preR=[shortpath_preR pmatRL_preR(22,22) pmatLR_preR(22,22)];
%             shortpath_postR=[shortpath_postR pmatRL_postR(22,22) pmatLR_postR(22,22)];            

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

try
    
goodbins=find(~isnan(safe_preR) & ~isnan(safe_postR));
[hRy, pRt]=ttest(safe_preR(goodbins),safe_postR(goodbins)); tt_safe_slogp=-log10(pRt)*sign(mean(safe_postR(goodbins)-safe_preR(goodbins)));               
[pR, hR]=signrank(safe_preR(goodbins),safe_postR(goodbins)); sr_safe_slogp=-log10(pR)*sign(mean(safe_postR(goodbins)-safe_preR(goodbins)));                

goodbins=find(~isnan(unsafe_preR) & ~isnan(unsafe_postR));
[hRy, pRt]=ttest(unsafe_preR(goodbins),unsafe_postR(goodbins)); tt_unsafe_slogp=-log10(pRt)*sign(mean(unsafe_postR(goodbins)-unsafe_preR(goodbins)));               
[pR, hR]=signrank(unsafe_preR(goodbins),unsafe_postR(goodbins)); sr_unsafe_slogp=-log10(pR)*sign(mean(unsafe_postR(goodbins)-unsafe_preR(goodbins)));                

goodbins=find(~isnan(right_preR) & ~isnan(right_postR));
[hRy, pRt]=ttest(right_preR(goodbins),right_postR(goodbins)); tt_right_slogp=-log10(pRt)*sign(mean(right_postR(goodbins)-right_preR(goodbins)));               
[pR, hR]=signrank(right_preR(goodbins),right_postR(goodbins)); sr_right_slogp=-log10(pR)*sign(mean(right_postR(goodbins)-right_preR(goodbins)));                

catch
    
tt_safe_slogp=NaN;               
sr_safe_slogp=NaN;                

tt_unsafe_slogp=NaN;               
sr_unsafe_slogp=NaN;                

tt_right_slogp=NaN;               
sr_right_slogp=NaN;                

end

figure(400); clf;
Rmax=.45;%max([pmatLR_preR(:)' pmatLR_postR(:)' pmatRL_preR(:)' pmatRL_postR(:)']);
Rmin=-.45%min([pmatLR_preR(:)' pmatLR_postR(:)' pmatRL_preR(:)' pmatRL_postR(:)']);
subplot(2,2,1); imagesc(pmatLR_preR); caxis([Rmin Rmax]); axis square; set(gca,'XTick',[],'YTick',[]); title('pre'); ylabel('LR'); colormap jet;%imagesc(PC_diffR);
subplot(2,2,2); imagesc(pmatLR_postR); caxis([Rmin Rmax]); axis square; set(gca,'XTick',[],'YTick',[]); title('post'); %imagesc(PC_diffR); title('place cells');
subplot(2,2,3); imagesc(pmatRL_preR); caxis([Rmin Rmax]); axis square; set(gca,'XTick',[],'YTick',[]); title('pre'); ylabel('LR'); colormap jet;%imagesc(PC_diffR);
subplot(2,2,4); imagesc(pmatRL_postR); caxis([Rmin Rmax]); axis square; set(gca,'XTick',[],'YTick',[]); title('post'); %imagesc(PC_diffR); title('place cells');

figure(1); clf;

Rmax=.7;%max([pmatLR_preR(:)' pmatLR_postR(:)' pmatRL_preR(:)' pmatRL_postR(:)']);
Rmin=-.3;%min([pmatLR_preR(:)' pmatLR_postR(:)' pmatRL_preR(:)' pmatRL_postR(:)']);
subplot(4,7,5); imagesc(pmatLR_preR); caxis([Rmin Rmax]); set(gca,'XTick',[],'YTick',[]); title('pre'); ylabel('LR'); %imagesc(PC_diffR);
xlabel(num2str(tt_unsafe_slogp));
subplot(4,7,6); imagesc(pmatLR_postR); caxis([Rmin Rmax]); set(gca,'XTick',[],'YTick',[]); title('post'); %imagesc(PC_diffR); title('place cells');
xlabel(num2str(sr_unsafe_slogp));

try
    
[hRt, pRt]=ttest(unsafe_preP,unsafe_postP);
[pR, hR]=signrank(unsafe_preP,unsafe_postP);
subplot(4,7,19); imagesc(pmatLR_preP); set(gca,'XTick',[],'YTick',[]); title('pre');  ylabel('LR'); %imagesc(PC_diffR);
xlabel(num2str(pR));
subplot(4,7,20); imagesc(pmatLR_postP); set(gca,'XTick',[],'YTick',[]); title('post');  %imagesc(PC_diffR); title('place cells');
xlabel(num2str(pRt));
subplot(4,7,12); imagesc(pmatRL_preR); caxis([Rmin Rmax]); set(gca,'XTick',[],'YTick',[]); title('pre'); ylabel('RL'); %imagesc(PC_diffR);
subplot(4,7,13); imagesc(pmatRL_postR); caxis([Rmin Rmax]); set(gca,'XTick',[],'YTick',[]); title('post'); %imagesc(PC_diffR); title('place cells');
xlabel(mean(unsafe_postR-unsafe_preR));
subplot(4,7,26); imagesc(pmatRL_preP); set(gca,'XTick',[],'YTick',[]); title('pre');  ylabel('RL'); %imagesc(PC_diffR);
subplot(4,7,27); imagesc(pmatRL_postP); set(gca,'XTick',[],'YTick',[]); title('post'); %imagesc(PC_diffR); title('place cells');
xlabel(mean(unsafe_postP-unsafe_preP));

end

PC_segmentR=[PC_segmentR; [mean(PC_diagR(L_zone)) mean([PC_LRdiagR(adjM_zoneLR) PC_RLdiagR(adjM_zoneRL)]) mean(PC_diagR(R_zone))]];

PC_segmentP=[PC_segmentP; [mean(PC_diagP(L_zone)) mean([PC_LRdiagP(adjM_zoneLR) PC_RLdiagP(adjM_zoneRL)]) mean(PC_diagP(R_zone))]];

                
                clear stackR stackP;
%                 stackR(1,:,:)=squeeze(all_diffR(1:23,1:23));
%                 stackR(2,:,:)=squeeze(all_diffR([1:23]+23,[1:23]+23));
%                 stackP(1,:,:)=squeeze(all_diffP(1:23,1:23));
%                 stackP(2,:,:)=squeeze(all_diffP([1:23]+23,[1:23]+23));
% 
%                 subplot(4,7,12); imagesc(squeeze(nanmean(stackR))); set(gca,'XTick',[],'YTick',[]);%subplot(4,7,10); imagesc(allmat_postR);
%                 subplot(4,7,26); imagesc(squeeze(nanmean(stackP))); set(gca,'XTick',[],'YTick',[]);%subplot(4,7,20); imagesc(allmat_postP);
                
                if r<5 %shocked first
%                     Rmatrix_all_shock1(r,:,:)=squeeze(nanmean(stackR));
%                     Pmatrix_all_shock1(r,:,:)=squeeze(nanmean(stackP));
                    Rmatrix_place_shock1((r*2-1):(r*2),:,:)=squeeze(pmat_postR)-squeeze(pmat_preR);
                    Rmatrix_place_shock1_pre((r*2-1):(r*2),:,:)=squeeze(pmat_preR);
                    Rmatrix_place_shock1_post((r*2-1):(r*2),:,:)=squeeze(pmat_postR);
                    Pmatrix_place_shock1(r,:,:)=squeeze(PC_diffP);
                elseif r==6 %shocked first
%                     Rmatrix_all_shock1(5,:,:)=squeeze(nanmean(stackR));
%                     Pmatrix_all_shock1(5,:,:)=squeeze(nanmean(stackP));
                    Rmatrix_place_shock1(9:10,:,:)=squeeze(pmat_postR)-squeeze(pmat_preR);
                    Rmatrix_place_shock1_pre(9:10,:,:)=squeeze(pmat_preR);
                    Rmatrix_place_shock1_post(9:10,:,:)=squeeze(pmat_postR);
                    Pmatrix_place_shock1(5,:,:)=squeeze(PC_diffP);
                elseif r>=7 & r<=12 %barrier
%                     Rmatrix_all_barrier(r-6,:,:)=squeeze(nanmean(stackR));
%                     Pmatrix_all_barrier(r-6,:,:)=squeeze(nanmean(stackP));
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
%                     Rmatrix_all_shock1(r-rf,:,:)=squeeze(nanmean(stackR));
%                     Pmatrix_all_shock1(r-rf,:,:)=squeeze(nanmean(stackP));
                    Rmatrix_place_shock1(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR)-squeeze(pmat_preR);
                    Rmatrix_place_shock1_pre(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_preR);
                    Rmatrix_place_shock1_post(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR);
                    Pmatrix_place_shock1(r-rf,:,:)=squeeze(PC_diffP);
                elseif r<=25 %scop + first shock
                    rf=16;
%                     Rmatrix_all_scopshock(r-16,:,:)=squeeze(nanmean(stackR));
%                     Pmatrix_all_scopshock(r-16,:,:)=squeeze(nanmean(stackP));
                    Rmatrix_place_scopshock(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR)-squeeze(pmat_preR);
                    Rmatrix_place_scopshock_pre(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_preR);
                    Rmatrix_place_scopshock_post(((r-rf)*2-1):((r-rf)*2),:,:)=squeeze(pmat_postR);
                    Pmatrix_place_scopshock(r-16,:,:)=squeeze(PC_diffP);
                else %scop alone
%                     Rmatrix_all_scop(r-24,:,:)=squeeze(nanmean(stackR));
%                     Pmatrix_all_scop(r-24,:,:)=squeeze(nanmean(stackP));
                    Rmatrix_place_scop(r-24,:,:)=squeeze(PC_diffR);
                    Rmatrix_place_scop_pre(r-24,:,:)=squeeze(pmat_preR);
                    Rmatrix_place_scop_post(r-24,:,:)=squeeze(pmat_postR);
                    Pmatrix_place_scop(r-24,:,:)=squeeze(PC_diffP);
                end
                
                %-------------------- plot sorted heatmaps
figure(2); clf;            
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
                
                peak=peak2;%(~moreplaceyLR_shockpre(isplace_preLR))=peak2(~moreplaceyLR_shockpre(isplace_preLR))+100;                 
                hmapLR_pre=sortrows([peak(:) hmapLR_pre],1);
                hmapLR_shockpre=sortrows([peak(:) hmapLR_shockpre],1);
                subplot_tight(4,7,[1 8 ]); imagesc(1000*hmapLR_pre(:,2:end)); title(['pre ' num2str(cellmat_shockcols(r,precol))]); set(gca,'XTick',[]); caxis([0 10]);; colormap hot;
                subplot_tight(4,7,[1 8 ]+1); imagesc(1000*hmapLR_shockpre(:,2:end)); title(['shkpre ' num2str(cellmat_shockcols(r,shkcol))]); set(gca,'XTick',[]); caxis([0 10]);; 
                
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
                
                peak=peak2;%(~moreplaceyRL_shockpre(isplace_preRL))=peak2(~moreplaceyRL_shockpre(isplace_preRL))+100;                 
                hmapRL_pre=sortrows([peak(:) hmapRL_pre],1);
                hmapRL_shockpre=sortrows([peak(:) hmapRL_shockpre],1);
                subplot_tight(4,7,[1 8 ]+14); imagesc(1000*hmapRL_pre(:,2:end)); title(['pre ' num2str(cellmat_shockcols(r,precol))]); set(gca,'XTick',[]); caxis([0 10]);; 
                subplot_tight(4,7,[1 8 ]+1+14); imagesc(1000*hmapRL_shockpre(:,2:end)); title(['shkpre ' num2str(cellmat_shockcols(r,shkcol))]); set(gca,'XTick',[]); caxis([0 10]);; 
                
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
                subplot_tight(4,7,[1 8 ]+3); imagesc(1000*hmapLR_post(:,2:end)); title(['post ' num2str(cellmat_shockcols(r,postcol))]); set(gca,'XTick',[]); caxis([0 10]);; 
                subplot_tight(4,7,[1 8 ]+2); imagesc(1000*hmapLR_shockpost(:,2:end)); title(['shkpost ' num2str(cellmat_shockcols(r,shkcol))]); set(gca,'XTick',[]); caxis([0 10]);; 

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
                subplot_tight(4,7,[1 8 ]+3+14); imagesc(1000*hmapRL_post(:,2:end)); title(['post ' num2str(cellmat_shockcols(r,postcol))]); set(gca,'XTick',[]); caxis([0 10]);; 
                subplot_tight(4,7,[1 8 ]+2+14); imagesc(1000*hmapRL_shockpost(:,2:end)); title(['shkpost ' num2str(cellmat_shockcols(r,shkcol))]); set(gca,'XTick',[]); caxis([0 10]);; 
 
figure(3); clf; 
                [pre_v_shockLR_R, pre_v_shockLR_P] = pv_heatmap_cormatrix(hmapLR_pre(:,2:end)',hmapLR_shockpre(:,2:end)');
                subplot_tight(2,2,1); imagesc(pre_v_shockLR_R); colormap jet;
                [pre_v_shockRL_R, pre_v_shockRL_P] = pv_heatmap_cormatrix(hmapRL_pre(:,2:end)',hmapRL_shockpre(:,2:end)');
                subplot_tight(2,2,2); imagesc(pre_v_shockRL_R); 
                [post_v_shockLR_R, post_v_shockLR_P] = pv_heatmap_cormatrix(hmapLR_post(:,2:end)',hmapLR_shockpost(:,2:end)');
                subplot_tight(2,2,3); imagesc(post_v_shockLR_R); 
                [post_v_shockRL_R, post_v_shockRL_P] = pv_heatmap_cormatrix(hmapRL_post(:,2:end)',hmapRL_shockpost(:,2:end)');
                subplot_tight(2,2,4); imagesc(post_v_shockRL_R); 

%                 figure(300); clf;
% %                 subplot_tight(4,7,[ 22]); plot3(predata.x,predata.y,predata.time/1000); set(gca,'XTick',[],'YTick',[],'ZTick',[],'ZLim',[-50 900]);
%                 %subplot_tight(4,7,[ 22]+1);
%                 plot3(data.x,data.y,data.time/1000); set(gca,'XTick',[],'YTick',[],'ZTick',[],'ZLim',[-50 900]); hold on;
%                          for i=1:size(data.LRint,1)
%                              plot3(data.x(data.LRint(i,1):data.LRint(i,2)),data.y(data.LRint(i,1):data.LRint(i,2)),data.time(data.LRint(i,1):data.LRint(i,2))/1000,'c');
%                          end                
%                          for i=1:size(data.RLint,1)
%                              plot3(data.x(data.RLint(i,1):data.RLint(i,2)),data.y(data.RLint(i,1):data.RLint(i,2)),data.time(data.RLint(i,1):data.RLint(i,2))/1000,'m');
%                          end                
%                          for i=1:length(shockmarks)
%                              shkframe=find((data.time/1000)>shockmarks(i),1,'first');
%                              scatter3(data.x(shkframe),data.y(shkframe),data.time(shkframe)/1000,'.r');
%                          end                
% %                 subplot_tight(4,7,[ 22]+2); plot3(postdata.x,postdata.y,postdata.time/1000); set(gca,'XTick',[],'YTick',[],'ZTick',[],'ZLim',[-50 900]);
                
                figure(600); clf;
                %subplot_tight(4,7,[ 22]+1); 
                plot(data.x,1:length(data.time)); set(gca,'XTick',[],'YTick',[]); hold on;
                         for i=1:size(data.LRint,1)
                             plot(data.x(data.LRint(i,1):data.LRint(i,2)),(data.LRint(i,1):data.LRint(i,2)),'c');
                         end                
                         for i=1:size(data.RLint,1)
                             plot(data.x(data.RLint(i,1):data.RLint(i,2)),(data.RLint(i,1):data.RLint(i,2)),'m');
                         end                
                         for i=1:length(shockmarks)
                             shkframe=find((1:length(data.time))>shockmarks(i),1,'first');
                             scatter(data.x(shkframe),(shkframe),'.r');
                         end                

                figure(1);

                subplot_tight(4,7,[ 22]); plot(predata.x,predata.time/1000); set(gca,'XTick',[],'YTick',[],'ZTick',[],'ZLim',[-50 900]);
                subplot_tight(4,7,[ 22]+2); plot(postdata.x,postdata.time/1000); set(gca,'XTick',[],'YTick',[],'ZTick',[],'ZLim',[-50 900]);
                subplot(4,7,25); hold off; bar([length(predata.LRint) 0 0 length(data.LRint) 0 0 length(postdata.LRint) 0 0],1);
                hold on; bar([0 length(predata.RLint) 0 0 length(data.RLint) 0 0 length(postdata.RLint) 0 0],1); set(gca,'XTick',[]);

isplace_allpre = (predata.place_scoreLR>minplacescore & predata.spikes_per_LRbee>minspikes) | (predata.place_scoreRL>minplacescore & predata.spikes_per_RLbee>minspikes);
isplace_allshk = (data.place_scoreLR>minplacescore & data.spikes_per_LRbee>minspikes) | (data.place_scoreRL>minplacescore & data.spikes_per_RLbee>minspikes);
isplace_allpost = (postdata.place_scoreLR>minplacescore & postdata.spikes_per_LRbee>minspikes) | (postdata.place_scoreRL>minplacescore & postdata.spikes_per_RLbee>minspikes);

preplace_recur=length(isplace_allshk & recur_shockpre_shockflag)/(sum(isplace_allpre)+sum(isplace_allshk));
postplace_recur=length(isplace_allshk & recur_shockpost_shockflag)/(sum(isplace_allpost)+sum(isplace_allshk));

% for i=1:23
%     minimumLR(i)=min([predata.vmap_LR(i) data.vmap_LR(i) postdata.vmap_LR(i)]);
%     minimumRL(i)=min([predata.vmap_RL(i) data.vmap_RL(i) postdata.vmap_RL(i)]);
% end

% try
% temp=data.dcurve_LR_bps;
% catch
% data=pv_tuning_curves_isplace_deconv_nomocc(data,2,0,0,minimumLR,minimumRL);
% % end
% % try
% % temp=predata.dcurve_LR_bps;
% % catch
% predata=pv_tuning_curves_isplace_deconv_nomocc(predata,2,0,0,minimumLR,minimumRL);
% % end
% % try
% % temp=postdata.dcurve_LR_bps;
% % catch
% postdata=pv_tuning_curves_isplace_deconv_nomocc(postdata,2,0,0,minimumLR,minimumRL);
%end

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


if ismember(r,[1:4 6 13:16 26:27])
    
    
      analysis_shift_dfpre=[analysis_shift_dfpre; [shiftdistLR_pre'; shiftdistRL_pre']];
      analysis_shift_dfpost=[analysis_shift_dfpost; [shiftdistLR_post'; shiftdistRL_post']];

      analysis_globalR_dfpre=[analysis_globalR_dfpre; globalR_pre];
      analysis_globalR_dfpost=[analysis_globalR_dfpost; globalR_post];

      dfpre_placecellmap=[dfpre_placecellmap; [predata.dcurve_LR(predata.isplace,:) predata.dcurve_RL(predata.isplace,:)]];
      dfshk_placecellmap=[dfshk_placecellmap; [data.dcurve_LR(data.isplace,:) data.dcurve_RL(data.isplace,:)]];
      dfpost_placecellmap=[dfpost_placecellmap; [postdata.dcurve_LR(postdata.isplace,:) postdata.dcurve_RL(postdata.isplace,:)]];
      
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
      
      analysis_globalR_scpre=[analysis_globalR_scpre; globalR_pre];
      analysis_globalR_scpost=[analysis_globalR_scpost; globalR_post];

      scpre_placecellmap=[scpre_placecellmap; [predata.dcurve_LR(predata.isplace,:) predata.dcurve_RL(predata.isplace,:)]];
      scshk_placecellmap=[scshk_placecellmap; [data.dcurve_LR(data.isplace,:) data.dcurve_RL(data.isplace,:)]];
      scpost_placecellmap=[scpost_placecellmap; [postdata.dcurve_LR(postdata.isplace,:) postdata.dcurve_RL(postdata.isplace,:)]];

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

      barpre_placecellmap=[barpre_placecellmap; [predata.dcurve_LR(predata.isplace,:) predata.dcurve_RL(predata.isplace,:)]];
      barshk_placecellmap=[barshk_placecellmap; [data.dcurve_LR(data.isplace,:) data.dcurve_RL(data.isplace,:)]];
      barpost_placecellmap=[barpost_placecellmap; [postdata.dcurve_LR(postdata.isplace,:) postdata.dcurve_RL(postdata.isplace,:)]];

                    hm_pre_bar_shortpath=[hm_pre_bar_shortpath; shortpath_preR];
                    hm_post_bar_shortpath=[hm_post_bar_shortpath; shortpath_postR];

                    hm_pre_bar_safe=[hm_pre_bar_safe; safe_preR];
                    hm_post_bar_safe=[hm_post_bar_safe; safe_postR];
                    
                    hm_pre_bar_unsafe=[hm_pre_bar_unsafe; unsafe_preR];
                    hm_post_bar_unsafe=[hm_post_bar_unsafe; unsafe_postR];      
end
               
                analysis_results(r,:)=[r sum(predata.isplaceLR | predata.isplaceRL)/length(predata.isplaceLR) sum(data.isplaceLR | data.isplaceRL)/length(data.isplaceLR) sum(postdata.isplaceLR | postdata.isplaceRL)/length(postdata.isplaceLR) ...
                    sum(predata.isplaceLR | predata.isplaceRL)  sum(data.isplaceLR | data.isplaceRL)  sum(postdata.isplaceLR | postdata.isplaceRL) ... %5 6 7
                    nanmedian([predata.dcurve_LR_bps(predata.isplaceLR) predata.dcurve_RL_bps(predata.isplaceRL)]) nanmedian([data.dcurve_LR_bps(data.isplaceLR) data.dcurve_RL_bps(data.isplaceRL)]) nanmedian([postdata.dcurve_LR_bps(postdata.isplaceLR) postdata.dcurve_RL_bps(postdata.isplaceRL)])... %8 9 10
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

[p_pre, h] = signrank(globalR_pre);                
[p_post, h] = signrank(globalR_post);

[h, p_pre_sp] = ttest(shortpath_preR);                
[h, p_post_sp] = ttest(shortpath_postR);

[h, p_pre_safe] = ttest(safe_preR);                
[h, p_post_safe] = ttest(safe_postR);

[h, p_pre_unsafe] = ttest(unsafe_preR);                
[h, p_post_unsafe] = ttest(unsafe_postR);

analysis_sigR(r,:) = [p_pre p_post p_pre_sp p_post_sp p_pre_safe p_post_safe p_pre_unsafe p_post_unsafe];

clear temp*;
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
analysis_placeinfo_shock=[analysis_placeinfo_shock predata.dcurve_LR_bps(predata.isplaceLR) predata.dcurve_RL_bps(predata.isplaceRL) data.dcurve_LR_bps(data.isplaceLR) data.dcurve_RL_bps(data.isplaceRL) postdata.dcurve_LR_bps(postdata.isplaceLR) postdata.dcurve_RL_bps(postdata.isplaceRL)];
analysis_nonplaceinfo_shock=[analysis_nonplaceinfo_shock predata.dcurve_LR_bps(~predata.isplaceLR) predata.dcurve_RL_bps(~predata.isplaceRL) data.dcurve_LR_bps(~data.isplaceLR) data.dcurve_RL_bps(~data.isplaceRL) postdata.dcurve_LR_bps(~postdata.isplaceLR) postdata.dcurve_RL_bps(~postdata.isplaceRL)];
      elseif ismember(r,[7:12])
          placemap_barrier=[placemap_barrier; recurplace_shockmap];
analysis_placeinfo_scopshock=[analysis_placeinfo_scopshock predata.dcurve_LR_bps(predata.isplaceLR) predata.dcurve_RL_bps(predata.isplaceRL) data.dcurve_LR_bps(data.isplaceLR) data.dcurve_RL_bps(data.isplaceRL) postdata.dcurve_LR_bps(postdata.isplaceLR) postdata.dcurve_RL_bps(postdata.isplaceRL)];
analysis_nonplaceinfo_scopshock=[analysis_nonplaceinfo_scopshock predata.dcurve_LR_bps(~predata.isplaceLR) predata.dcurve_RL_bps(~predata.isplaceRL) data.dcurve_LR_bps(~data.isplaceLR) data.dcurve_RL_bps(~data.isplaceRL) postdata.dcurve_LR_bps(~postdata.isplaceLR) postdata.dcurve_RL_bps(~postdata.isplaceRL)];
      else
          placemap_scopshock=[placemap_scopshock; recurplace_shockmap];
analysis_placeinfo_barrier=[analysis_placeinfo_barrier predata.dcurve_LR_bps(predata.isplaceLR) predata.dcurve_RL_bps(predata.isplaceRL) data.dcurve_LR_bps(data.isplaceLR) data.dcurve_RL_bps(data.isplaceRL) postdata.dcurve_LR_bps(postdata.isplaceLR) postdata.dcurve_RL_bps(postdata.isplaceRL)];
analysis_nonplaceinfo_barrier=[analysis_nonplaceinfo_barrier predata.dcurve_LR_bps(~predata.isplaceLR) predata.dcurve_RL_bps(~predata.isplaceRL) data.dcurve_LR_bps(~data.isplaceLR) data.dcurve_RL_bps(~data.isplaceRL) postdata.dcurve_LR_bps(~postdata.isplaceLR) postdata.dcurve_RL_bps(~postdata.isplaceRL)];
      end
      


               
    end %rat loop -----------------------------------------------------------------------------------------------------------------

%Figure2;
    


