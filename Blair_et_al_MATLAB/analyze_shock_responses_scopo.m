shk_isplace = data.isplace_N;
shk_isplace1 = data.isplace_N & 0;
shk_isplace2 = data.isplace_N & 0;

alphathresh=.05; 

data.S=data.deconv>0;

    eventtimes=scopshocktimes;
    figure(102); %graph window for shock responses
    if r==17; clf; end
    splot=r-16;

analyze_peakshk1=NaN(size(data.S,1),1);
analyze_peakshk2=NaN(size(data.S,1),1);


%--- measure shock1 and shock2 responses in ROI for each cell
ROIstart=pethwid+2; ROIend=pethwid*2+1; %ROI surrounding shock during which shock responsiveness will be analyzed
shockframe1=find(data.time/1000>eventtimes(rat(r).shockrow,1),1,'first')-3; %frame number of first shock in session
shock1_PETH=false(size(data.S,1),pethwid*2+1); %initialize binary shock1 response PETH to all false (so memory does not need to be dynamically allocated)
cellmap_index=zeros(size(data.S,1),1); %on first loop through cells, we will gather index numbers into cellmap for later use
for j=1:size(data.S,1) %loop through cells
    try
    temp=find(cmap(:,shk_column)==j); %get cellmap row index of this cell
    if ~isempty(temp) %store row index if it exists
        cellmap_index(j)=temp;
    end
    end
    shock1_PETH(j,:)=data.S(j,(shockframe1-pethwid):(shockframe1+pethwid)); %binary response PETH to shock1
    ROIcount1(j)=sum(1*shock1_PETH(j,ROIstart:ROIend));%count spikes in ROI after shock1
end
shockframe2=find(data.time/1000>eventtimes(rat(r).shockrow,2),1,'first')-3; %frame number of second shock in session
shock2_PETH=false(size(data.S,1),pethwid*2+1); %initialize binary shock2 response PETH to all false (so memory does not need to be dynamically allocated)
for j=1:size(data.S,1)
    shock2_PETH(j,:)=data.S(j,(shockframe2-pethwid):(shockframe2+pethwid)); %binary response PETH to shock2
    ROIcount2(j)=sum(1*shock2_PETH(j,ROIstart:ROIend));%count spikes in ROI after shock2
end

%            figure(102); clf; plot3(data.x,data.y,data.time/1000); %plot path for shock session


%--- measure responses in baseline period for each cell in each runnning direction (LR and RL)
LR_BL=false(size(data.S,1),size(data.LRint,1),pethwid*2+1); %initialize LR journey baseline period spike counts to 0
for i=1:size(data.LRint,1) %loop through LR journeys
    %frame0 = frame number when rat first entered shock zone on this LR journey
    frame0=data.LRint(i,1)+find(data.x(data.LRint(i,1):data.LRint(i,2))<-30 & data.x(1+[data.LRint(i,1):data.LRint(i,2)])>=-30)-1;
    if frame0<pethwid %corrective measures if PETH window starts before beginning of session
        temp=frame0-1;
        sdex=1+pethwid-temp;
    else
        temp=pethwid;
        sdex=1;
    end
    for j=1:size(data.S,1) %loop through cells
        LR_BL(j,i,sdex:end)=data.S(j,(frame0-temp):(frame0+pethwid)); %compute BL PETH triggered on this journey's shock zone entry for this cell
        LR_ROI(j).count(i)=sum(1*LR_BL(j,i,ROIstart:ROIend)); %compute spike count during ROI
    end
end
RL_BL=false(size(data.S,1),size(data.RLint,1),pethwid*2+1); %initialize RL journey baseline period spike counts to 0
for i=1:size(data.RLint,1) %loop through RL journeys
    %frame0 = frame number when rat first entered shock zone on this LR journey
    frame0=data.RLint(i,1)+find(data.x(data.RLint(i,1):data.RLint(i,2))>15 & data.x(1+[data.RLint(i,1):data.RLint(i,2)])<=15)-1;
    if frame0<pethwid %corrective measures if PETH window starts before beginning of session
        temp=frame0-1;
        sdex=1+pethwid-temp;
    else
        temp=pethwid;
        sdex=1;
    end
    for j=1:size(data.S,1)
        RL_BL(j,i,sdex:end)=data.S(j,(frame0-temp):(frame0+pethwid)); %compute BL PETH triggered on this journey's shock zone entry for this cell
        RL_ROI(j).count(i)=sum(1*RL_BL(j,i,ROIstart:ROIend)); %compute spike count during ROI
    end
end

%-- classify shock-excited cells as those that beat 99% confidence against a Poisson distribution fitted to BL shock counts
SHKexc1=false(1,size(data.S,1)); %initialize cell flags for excitation by shock1
SHKexc2=false(1,size(data.S,1)); %initialize cell flags for excitation by shock2

SHKinh1=false(1,size(data.S,1)); %initialize cell flags for excitation by shock1
SHKinh2=false(1,size(data.S,1)); %initialize cell flags for excitation by shock2

    labelshk1='X';
    labelshk2='X';

approach_1=false(1,size(data.S,1)); %initialize cell flags for excitation by shock1
approach_2=false(1,size(data.S,1)); %initialize cell flags for excitation by shock2    
    
if data.x(shockframe1-15)<data.x(shockframe1) %if shock1 occurred during an LR pass, then we use LR baseline for testing response to shock1
    shkdir1=-1;
    labelshk1=['LR' num2str(size(data.LRint,1))];
        data.shockLR=true;
    for j=1:size(data.S,1) %loop through cells
[temp analyze_peakshk1(j)]=max(data.dcurve_LR(j,:));
if analyze_peakshk1(j)<13 | (analyze_peakshk1(j)==13 & sum(data.dcurve_LR(j,1:12))>sum(data.dcurve_LR(j,14:end)))
    approach_1(j)=true;
end
shk_isplace1(j)=data.isplaceLR_N(j);
        if ROIcount1(j)>0 %if the cell spiked during the ROI, then see if it responded significantly above baseline
            fitpoiss = fitdist(LR_ROI(j).count','Poisson'); %fit a Poisson distribution to baseline spike counts
%            fitpoiss = fitdist([0 0 0 0 0 0]','Poisson'); %six BL passes with zero spikes each are sufficient for ci(2)<1, making one shock-evoked spike a significant response
            ci = paramci(fitpoiss,'Alpha',alphathresh); %find .01 confidence limits on fitted Poisson distribution
            ci_list1(j)=ci(2);
            if ROIcount1(j)>=ci(2) %if observed spike count exceeds 99% confidence,
                SHKexc1(j)=true;   %then this cell responded to the shock
            end
            if ROIcount1(j)<ci(1) %if observed spike count exceeds 99% confidence,
                SHKinh1(j)=true;   %then this cell responded to the shock
            end
        else
            ci_list1(j)=-6;
        end
    end
end
if data.x(shockframe2-15)<data.x(shockframe2) %if shock2 occurred during an LR pass, then we use LR baseline for testing response to shock2
    shkdir2=-1;
    labelshk2=['LR' num2str(size(data.LRint,1))];
        data.shockLR=true;
    for j=1:size(data.S,1)
[temp analyze_peakshk2(j)]=max(data.dcurve_LR(j,:));
if analyze_peakshk2(j)<13 | (analyze_peakshk2(j)==13 & sum(data.dcurve_LR(j,1:12))>sum(data.dcurve_LR(j,14:end)))
    approach_2(j)=true;
end
shk_isplace2(j)=data.isplaceLR_N(j);
        if ROIcount2(j)>0
            fitpoiss = fitdist(LR_ROI(j).count','Poisson');
            ci = paramci(fitpoiss,'Alpha',alphathresh);
            ci_list2(j)=ci(2);
            if ROIcount2(j)>=ci(2)
                SHKexc2(j)=true;
            end
            if ROIcount2(j)<ci(1) %if observed spike count exceeds 99% confidence,
                SHKinh2(j)=true;   %then this cell responded to the shock
            end
        else
            ci_list2(j)=-6;
        end
    end
end
if data.x(shockframe1-15)>data.x(shockframe1) %if shock1 occurred during an RL pass, then we use RL baseline for testing response to shock1
    shkdir1=1;
    labelshk1=['RL' num2str(size(data.RLint,1))];
        data.shockRL=true;
    for j=1:size(data.S,1)
[temp analyze_peakshk1(j)]=max(data.dcurve_RL(j,23:-1:1));
if analyze_peakshk1(j)<13 | (analyze_peakshk1(j)==13 & sum(data.dcurve_RL(j,1:12))>sum(data.dcurve_RL(j,14:end)))
    approach_1(j)=true;
end
shk_isplace1(j)=data.isplaceRL_N(j);
        if ROIcount1(j)>0
            fitpoiss = fitdist(RL_ROI(j).count','Poisson');
            ci = paramci(fitpoiss,'Alpha',alphathresh);
            ci_list1(j)=ci(2);
            if ROIcount1(j)>=ci(2) %if observed spike count exceeds 99% confidence,
                SHKexc1(j)=true;   %then this cell responded to the shock
            end
            if ROIcount1(j)<ci(1) %if observed spike count exceeds 99% confidence,
                SHKinh1(j)=true;   %then this cell responded to the shock
            end
        else
            ci_list1(j)=-6;
        end
    end
end
if data.x(shockframe2-15)>data.x(shockframe2) %if shock2 occurred during an RL pass, then we use RL baseline for testing response to shock2
    shkdir2=1;
    labelshk2=['RL' num2str(size(data.RLint,1))];
        data.shockRL=true;
    for j=1:size(data.S,1)
[temp analyze_peakshk2(j)]=max(data.dcurve_RL(j,23:-1:1));
if analyze_peakshk2(j)<13 | (analyze_peakshk2(j)==13 & sum(data.dcurve_RL(j,1:12))>sum(data.dcurve_RL(j,14:end)))
    approach_2(j)=true;
end
shk_isplace2(j)=data.isplaceRL_N(j);
        if ROIcount2(j)>0
            fitpoiss = fitdist(RL_ROI(j).count','Poisson');
            ci = paramci(fitpoiss,'Alpha',alphathresh);
            ci_list2(j)=ci(2);
            if ROIcount2(j)>=ci(2)
                SHKexc2(j)=true;
            end
            if ROIcount2(j)<ci(1) %if observed spike count exceeds 99% confidence,
                SHKinh2(j)=true;   %then this cell responded to the shock
            end
        else
            ci_list2(j)=-6;
        end
    end
end

SHKexc=[SHKexc1 | SHKexc2]; %flag cell as shock responsive if it responded to either shock1 or shock2

SHKinh=[SHKinh1 | SHKinh2]; %flag cell as shock responsive if it responded to either shock1 or shock2
    
if SHKexc & SHKinh %cell is classified as non responsive if it has conflicting responses
    SHKexc=false;
    SHKinh=false;
end
    
    %plot population averaged responses of shock responsive and non responsive cells to shock1 and shock2
    %NOTE: if only 1 shock was received, same response appears in both shock1 and shock2 windows
    %NOT Z SCORE NORMALIZED BEFORE TAKING POP AVERAGE!!!! (need to fix this if figure is to be published)


    subplot(11,2,(splot-1)*2+1); hold off;%pop avg of shock responsive cells is plotted in blue, non responsive cells in black
    shadedErrorBar(-pethwid:pethwid,mean(1*shock1_PETH(SHKexc1,:)),std(1*shock1_PETH(SHKexc1,:))/sqrt(sum(1*SHKexc1)),'lineProps','b'); hold on;
    shadedErrorBar(-pethwid:pethwid,mean(1*shock1_PETH(~SHKexc1 & ~SHKinh1,:)),std(1*shock1_PETH(~SHKexc1 & ~SHKinh1,:))/sqrt(sum(1*(~SHKexc1 & ~SHKinh1))),'lineProps','k'); hold on;
    shadedErrorBar(-pethwid:pethwid,mean(1*shock1_PETH(~SHKinh1,:)),std(1*shock1_PETH(~SHKinh1,:))/sqrt(sum(1*~SHKinh1)),'lineProps','r'); hold on;
    ylabel(sum(SHKexc)/length(SHKexc));
%    ylabel(['Rat ' num2str(r)]); set(gca,'XTick',[]);
    title(labelshk1);
    subplot(11,2,(splot-1)*2+2); hold off;
    shadedErrorBar(-pethwid:pethwid,mean(1*shock2_PETH(SHKexc2,:)),std(1*shock2_PETH(SHKexc2,:))/sqrt(sum(1*SHKexc2)),'lineProps','b'); hold on;
    shadedErrorBar(-pethwid:pethwid,mean(1*shock2_PETH(~SHKexc2 & ~SHKinh2,:)),std(1*shock2_PETH(~SHKexc2 & ~SHKinh2,:))/sqrt(sum(1*(~SHKexc2 & ~SHKinh2))),'lineProps','k'); hold on;
    shadedErrorBar(-pethwid:pethwid,mean(1*shock2_PETH(~SHKinh2,:)),std(1*shock2_PETH(~SHKinh2,:))/sqrt(sum(1*~SHKinh2)),'lineProps','r'); hold on;
    set(gca,'XTick',[]);
    title(labelshk2);    
% for i=1:size(data.deconv,1)
%                 highbinsLR(i)=sum(data.dcurve_LR(i,:)>max(data.dcurve_LR(i,:))*.7);
%                 highbinsRL(i)=sum(data.dcurve_RL(i,:)>max(data.dcurve_RL(i,:))*.7);
% end
% 
% minplacescore=2;
% minspikes=1;
% maxhighbins=12;
% 
% shk_isplaceLR=highbinsLR<maxhighbins & data.place_scoreLR>minplacescore & data.spikes_per_LRbee>minspikes
% shk_isplaceRL=highbinsRL<maxhighbins & data.place_scoreRL>minplacescore & data.spikes_per_RLbee>minspikes
% shk_isplace=shk_isplaceLR | shk_isplaceRL;    

% if r==5 | r==6
% allshockcells = [allshockcells; 1*[shock2_PETH(SHKexc2,:)]];
% allnonshockcells = [allnonshockcells; 1*[shock2_PETH(~SHKexc2,:)]];
% else

% allshockcells_approachscop=[allshockcells_approachscop; approach_1(SHKexc1(:) & shk_isplace1(:))']; %initialize cell flags for excitation by shock1
% allnonshockcells_approachscop=[allnonshockcells_approachscop; approach_1(~SHKexc1(:) & shk_isplace1(:))']; %initialize cell flags for excitation by shock1
% allshockcells_scop = [allshockcells_scop; 1*shock1_PETH(SHKexc1(:) & shk_isplace1(:),:)];
% allnonshockcells_scop = [allnonshockcells_scop; 1*shock1_PETH(~SHKexc1(:) & shk_isplace1(:),:)];
% %allshockcells_peakscop = [allshockcells_peakscop; analyze_peakshk1(SHKexc1(:) & shk_isplace1(:))]; 
% %allnonshockcells_peakscop = [allnonshockcells_peakscop; analyze_peakshk1(~SHKexc1(:) & shk_isplace1(:))];
% 
% allshockcells_approachscop=[allshockcells_approachscop; approach_1(SHKexc1(:) & shk_isplace1(:))']; %initialize cell flags for excitation by shock1
% allnonshockcells_approachscop=[allnonshockcells_approachscop; approach_1(~SHKexc1(:) & shk_isplace1(:))']; %initialize cell flags for excitation by shock1
% allshockcells_scop = [allshockcells_scop; 1*shock1_PETH(SHKexc1(:) & shk_isplace(:),:)];
% allnonshockcells_scop = [allnonshockcells_scop; 1*shock1_PETH(~SHKexc1(:) & shk_isplace(:),:)];
% allshockcells_peakscop = [allshockcells_peakscop; analyze_peakshk1(SHKexc1(:) & shk_isplace1(:))]; 
% allnonshockcells_peakscop = [allnonshockcells_peakscop; analyze_peakshk1(~SHKexc1(:) & shk_isplace1(:))];

analyze_percplaceindir_shockresp(r) = sum( (SHKexc1(:) & shk_isplace1(:))  ); %/sum(SHKexc(:) & shk_isplace(:));
analyze_percnotplaceindir_shockresp(r) = sum( (SHKexc1(:) & ~shk_isplace1(:))  ); %/sum(SHKexc(:) & shk_isplace(:));
ncurves=sum(SHKexc1(:)); 

allnonshockcells_peakscop = [allnonshockcells_peakscop; analyze_peakshk1(~SHKexc(:) & shk_isplace(:))];
allnonshockcells_approachscop=[allnonshockcells_approachscop; approach_1(~SHKexc(:) & shk_isplace(:))']; %initialize cell flags for excitation by shock1

if ~(shockframe1 == shockframe2) & ~(shkdir1 == shkdir2)
   % allshockcells_approachscop=[allshockcells_approachscop; approach_2(SHKexc2(:) & shk_isplace2(:))']; %initialize cell flags for excitation by shock1
   % allnonshockcells_approachscop=[allnonshockcells_approachscop; approach_2(~SHKexc2(:) & shk_isplace2(:))']; %initialize cell flags for excitation by shock1
    %allshockcells_peakscop = [allshockcells_peakscop; analyze_peakshk2(SHKexc2(:) & shk_isplace2(:))]; 
    %allnonshockcells_peakscop = [allnonshockcells_peakscop; analyze_peakshk2(~SHKexc2(:) & shk_isplace2(:))];
    
    analyze_percplaceindir_shockresp(r) = analyze_percplaceindir_shockresp(r)+sum( (SHKexc2(:) & shk_isplace2(:))  ); %/sum(SHKexc(:) & shk_isplace(:));
    analyze_percnotplaceindir_shockresp(r) = analyze_percnotplaceindir_shockresp(r)+sum( (SHKexc2(:) & ~shk_isplace2(:))  ); %/sum(SHKexc(:) & shk_isplace(:));

    allnonshockcells_peakscop = [allnonshockcells_peakscop; analyze_peakshk2(~SHKexc(:) & shk_isplace(:))];
    allnonshockcells_approachscop=[allnonshockcells_approachscop; approach_2(~SHKexc(:) & shk_isplace(:))']; %initialize cell flags for excitation by shock1
    %PETH is average of responses to first and second shocks
    %allshockcells_scop = [allshockcells_scop; (1*[shock1_PETH(SHKexc(:) & shk_isplace(:),:)] + 1*[shock2_PETH(SHKexc(:) & shk_isplace(:),:)])/2 ];
    %allnonshockcells_scop = [allnonshockcells_scop; (1*[shock1_PETH(~SHKexc(:) & shk_isplace(:),:)] + 1*[shock2_PETH(~SHKexc(:) & shk_isplace(:),:)])/2];

    ncurves=ncurves+sum(SHKexc2(:)); 
end
    %PETH is just response to first shock
    allshockcells_scop = [allshockcells_scop; 1*shock1_PETH(SHKexc1(:) & shk_isplace(:),:); 1*shock2_PETH(~SHKexc1(:) & SHKexc2(:) & shk_isplace(:),:)];
    allnonshockcells_scop = [allnonshockcells_scop; 1*shock1_PETH(~SHKexc(:) & shk_isplace(:),:)];
    
    allshockcells_peakscop = [allshockcells_peakscop; analyze_peakshk1(SHKexc1(:) & shk_isplace(:)); analyze_peakshk2(~SHKexc1(:) & SHKexc2(:) & shk_isplace(:))]; 
    allshockcells_approachscop=[allshockcells_approachscop; approach_1(SHKexc1(:) & shk_isplace(:))'; approach_2(~SHKexc1(:) & SHKexc2(:) & shk_isplace(:))']; %initialize cell flags for excitation by shock1
%end

analyze_percplaceindir_shockresp(r) = analyze_percplaceindir_shockresp(r)/ncurves;
analyze_percnotplaceindir_shockresp(r) = analyze_percnotplaceindir_shockresp(r)/ncurves;

analyze_percplace_shockresp(r) = sum( (SHKexc(:) & shk_isplace(:))  )/sum(shk_isplace(:));

analyze_percplaceindir_shockresp(r) = sum( (SHKexc1(:) & shk_isplace(:) & shk_isplace1(:)) | (SHKexc2(:) & shk_isplace(:) & shk_isplace2(:)) )/sum(shk_isplace(:));
analyze_percnotplaceindir_shockresp(r) = sum( (SHKexc1(:) & shk_isplace(:) & ~shk_isplace1(:)) | (SHKexc2(:) & shk_isplace(:) & ~shk_isplace2(:)) )/sum(SHKexc(:) & shk_isplace(:));

analyze_numplace_shockresp(r) = sum( (SHKexc(:) & shk_isplace(:)) );

analyze_percplace_shockinh(r) = sum(SHKinh(:) & shk_isplace(:))/sum(shk_isplace(:));
analyze_numplace_shockinh(r) = sum(SHKinh(:) & shk_isplace(:));

analyze_percplace_shocknon(r) = sum(~SHKexc(:) & ~SHKinh(:) & shk_isplace(:))/sum(shk_isplace(:));
analyze_numplace_shocknon(r) = sum(~SHKexc(:) & ~SHKinh(:) & shk_isplace(:));
% %fetch cellmap row indices for cells that recurred between pre-shock and shock session
% recurdex1=find(cmap(:,cellmat_shockcols(r,1))>0 & cmap(:,cellmat_shockcols(r,2))>0);
% recur1=false(1,size(data.S,1));  %initialize recurrence flags
% recur1(cmap(recurdex1,cellmat_shockcols(r,2)))=true; %logical vector of cells recorded in shock session, flagged if they recur in pre-shock session
% %the code below will drop cells if they are not place cells in at least
% %one of the two sessions in the pair being analyzed
% if analysis_code==4 | analysis_code==2  | analysis_code==3
%     is_place_in_other=false(size(data.S,1),1);  %initialize recurrence flags
%     for j=1:length(otherdata.isplace) %find which of this session's cells are place cells in the partner session
%         if otherdata.isplace(j) %if this was a place cell in the other session
%             temp=find(cmap(:,cellmat_shockcols(r,otherdex))==j); %get cellmap row index of this cell
%             if ~isempty(temp) %store row index if it exists
%                 temp2=cmap(temp,cellmat_shockcols(r,ddex)); %cell number in the current session
%                 if temp2>0
%                     is_place_in_other(temp2)=true;
%                 end
%             end
%         end
%     end
%     recur1=recur1 & (data.isplace(:)' | is_place_in_other(:)'); %only keep recurring cells that are place cells in the current OR other session session
% end

% % % % % % % % %fetch cellmap row indices for cells that recurred between post-shock and shock session
% % % % % % % % recurdex2=find(cmap(:,cellmat_shockcols(r,3))>0 & cmap(:,cellmat_shockcols(r,2))>0);
% % % % % % % % recur2=false(1,size(data.S,1)); %initialize recurrence flags
% % % % % % % % recur2(cmap(recurdex2,cellmat_shockcols(r,2)))=true; %logical vector of cells recorded in shock session, flagged if they recur in post-shock session
% % % % % % % % %the code below will drop cells if they are not place cells in at least
% % % % % % % % %one of the tweo sessions in the pair being analyzed
% % % % % % % % if analysis_code==5 | analysis_code==2  | analysis_code==3
% % % % % % % %     is_place_in_other=false(size(data.S,1),1);  %initialize recurrence flags
% % % % % % % %     for j=1:length(otherdata.isplace) %find which of this session's cells are place cells in the partner session
% % % % % % % %         if otherdata.isplace(j) %if this was a place cell in the other session
% % % % % % % %             temp=find(cmap(:,cellmat_shockcols(r,otherdex))==j); %get cellmap row index of this cell
% % % % % % % %             if ~isempty(temp) %store row index if it exists
% % % % % % % %                 temp2=cmap(temp,cellmat_shockcols(r,ddex)); %cell number in the current session
% % % % % % % %                 if temp2>0
% % % % % % % %                     is_place_in_other(temp2)=true;
% % % % % % % %                 end
% % % % % % % %             end
% % % % % % % %         end
% % % % % % % %     end
% % % % % % % %     recur2=recur2(:)' & (data.isplace(:)' | is_place_in_other(:)'); %only keep recurring cells that are place cells in the current OR other session session
% % % % % % % % end



% if analysis_code==2
%     recur1=recur1 & SHKexc;
%     recur2=recur2 & SHKexc;
% elseif analysis_code==3
%     recur1=recur1 & ~SHKexc;
%     recur2=recur2 & ~SHKexc;
% end

% if ~isfield(data,'SHKexc')
%     data.SHKexc=SHKexc;
%     save([rootdir rat(r).name '_linear' num2str(cellmat_shockcols(r,ddex)) '_data'],'data');
% end
% 
% %chi square tests on whether shock responsive cells were more likely to recur in pre vs post shock sessions
% chisqmat=[sum(1*(recur1 & SHKexc)) sum(1*(recur2 & SHKexc)); sum(1*(recur1 & ~SHKexc)) sum(1*(recur2 & ~SHKexc))];
% if ~exist('all_chmat')
%     all_chmat=chisqmat;
% else
%     all_chmat=all_chmat+chisqmat;
% end
% [ chi_sq(r), chi_DF(r), chi_p(r) ] = chiSquareTestTable(chisqmat);
% 
% ALLrows=false(size(cmap,1),1);
% goodex=cellmap_index(cellmap_index>0);
% ALLrows(goodex)=true;
% 
% %flag rows for shock responsive cells
% SHKrows=false(size(cmap,1),1);
% SHKrows(goodex)=SHKexc(cmap(goodex,cellmat_shockcols(r,2)));