load 'Behav_effects.mat'; % gjb

%index numbers for rat groups to be analyzed
df_shock1=[1 2 3 4 10   5 6 8 13 9   11 12 14 15];
df_female=[1 4 8 12:15];
scop_shock1=[10 11 12 14 15];
scop_shock2=[5 6 7 8 10 ];
scop_female=[7 8 12 14 15];
barrier=[5 6 8 13 9 12];
bar_female=[8 12 13];
scop_only1=[11 14 10];
scop_only2=[9];

% convert total number of short path rewards during first 10 m of each session into rewards per minute 

df_s1=Effects.num_short.shock(df_shock1(1:10),[3 6])/10;
df_s2=Effects.num_short.shock(df_shock1(11:14),[3 6])/10;
df_f=Effects.num_short.shock(df_female,[3 6])/10;

sc_s1=Effects.num_short.scoposhk(scop_shock1,[3 6])/10;
sc_s2=Effects.num_short.scoposhk(scop_shock2,[3 6])/10;
sc_f=Effects.num_short.scoposhk(scop_female,[3 6])/10;

scspath_s1=Effects.num_short.scoposhk(scop_shock1,[3 4])/10;
scspath_s2=Effects.num_short.scoposhk(scop_shock2,[3 4])/10;

barr=Effects.num_short.barrier(barrier,[3 6])/10;
bar_f=Effects.num_short.barrier(bar_female,[3 6])/10;

so1=Effects.num_short.scopo(scop_only1,[3 6])/10;
so2=Effects.num_short.scopo(scop_only2,[3 6])/10;

figure(300); clf;

subplot(3,4,[3 4]); % 48h retention plots

%drug free shock condition
plot(df_s1','o-k'); %shk1 rats have circle symbols
hold on;
plot(df_s2','s-k'); %shk2 rats have square symbols

plot(df_f','-c'); %highlight females in cyan
[p_df_48, h]=signrank([df_s1(:,1); df_s2(:,1)],[df_s1(:,2); df_s2(:,2)])

%scopolamine shock condition
plot(3:4,sc_s1','o-r'); hold on;
plot(3:4,sc_s2','s-r');
plot(3:4,sc_f','-c');
[p_sc_48, h]=signrank([sc_s1(:,1); sc_s2(:,1)],[sc_s1(:,2); sc_s2(:,2)])

%barrier training condition
plot(5:6,barr','o-b');
plot(5:6,bar_f','-c');
[p_barr_48, h]=signrank(barr(:,1),barr(:,2))

[p_scspath_48, h]=signrank([scspath_s1(:,1); scspath_s2(:,1)],[scspath_s1(:,2); scspath_s2(:,2)])

set(gca,'XLim',[0 9]);

subplot(3,4,10); %shk1 vs shk2 compairson for 48 h retention scatter plot
 plot(1:2,[median(df_s1(:,2)./df_s1(:,1)) median(df_s2(:,2)./df_s2(:,1))],'o-k')%,[std(df_s1(:,2)./df_s1(:,1))/sqrt(size(df_s1,1)-1) std(df_s2(:,2)./df_s2(:,1))/sqrt(size(df_s2,1)-1)],'k');
 hold on;
 scatter(.85+0*(df_s1(:,2)./df_s1(:,1)),(df_s1(:,2)./df_s1(:,1)),'ok')%,[std(df_s1(:,2)./df_s1(:,1))/sqrt(size(df_s1,1)-1) std(df_s2(:,2)./df_s2(:,1))/sqrt(size(df_s2,1)-1)],'k');
 scatter(1.85+0*(df_s2(:,2)./df_s2(:,1)),(df_s2(:,2)./df_s2(:,1)),'sk')%,[std(df_s1(:,2)./df_s1(:,1))/sqrt(size(df_s1,1)-1) std(df_s2(:,2)./df_s2(:,1))/sqrt(size(df_s2,1)-1)],'k');
 plot(1:2,[median(sc_s1(:,2)./sc_s1(:,1)) median(sc_s2(:,2)./sc_s2(:,1))],'s-r');%,[std(sc_s1(:,2)./sc_s1(:,1))/sqrt(size(sc_s1,1)-1) std(sc_s2(:,2)./sc_s2(:,1))/sqrt(size(sc_s2,1)-1)],'r');
 scatter(1.15+0*(sc_s1(:,2)./sc_s1(:,1)),(sc_s1(:,2)./sc_s1(:,1)),'or')%,[std(sc_s1(:,2)./sc_s1(:,1))/sqrt(size(sc_s1,1)-1) std(sc_s2(:,2)./sc_s2(:,1))/sqrt(size(sc_s2,1)-1)],'k');
 scatter(2.15+0*(sc_s2(:,2)./sc_s2(:,1)),(sc_s2(:,2)./sc_s2(:,1)),'sr')%,[std(df_s1(:,2)./df_s1(:,1))/sqrt(size(df_s1,1)-1) std(df_s2(:,2)./df_s2(:,1))/sqrt(size(df_s2,1)-1)],'k');
set(gca,'XLim',[0 3],'YLim',[0 1.5]);

[p_scop_s1v2, h]=ranksum((sc_s1(:,2)./sc_s1(:,1)), (sc_s2(:,2)./sc_s2(:,1)))
[p_df_s1v2, h]=ranksum(df_s1(:,2)./df_s1(:,1), df_s2(:,2)./df_s2(:,1))


subplot(3,4,12); %shk1 vs shk2 compairson for 48 h retention box plot (bars not in same order as paper figure)
x=[(df_s1(:,2)./df_s1(:,1)); (sc_s1(:,2)./sc_s1(:,1)); (df_s2(:,2)./df_s2(:,1)); (sc_s2(:,2)./sc_s2(:,1))];
g=[(1+0*df_s1(:,2)./df_s1(:,1)); 2+0*(df_s2(:,2)./df_s2(:,1)); 3+0*(sc_s1(:,2)./sc_s1(:,1)); 4+0*(sc_s2(:,2)./sc_s2(:,1))];
boxplot(x,g);
set(gca,'XLim',[0 5],'YLim',[0 1.5]);

[p_drug, h]=ranksum([(df_s1(:,2)./df_s1(:,1)); (df_s2(:,2)./df_s2(:,1))],[(sc_s1(:,2)./sc_s1(:,1)); (sc_s2(:,2)./sc_s2(:,1))]);
[p_shnum, h]=ranksum([(df_s1(:,2)./df_s1(:,1)); (sc_s1(:,2)./sc_s1(:,1))],[(df_s2(:,2)./df_s2(:,1)); (sc_s2(:,2)./sc_s2(:,1))]);

[p_drugshock1, h]=ranksum(df_s1(:,2)./df_s1(:,1), (sc_s1(:,2)./sc_s1(:,1)) );
[p_drugshock2, h]=ranksum(df_s2(:,2)./df_s2(:,1), (sc_s2(:,2)./sc_s2(:,1)) );

[p_drugshock, h]=ranksum(df_s1(:,2)./df_s1(:,1), (sc_s1(:,2)./sc_s1(:,1)) );


df_s1=Effects.num_short.shock(df_shock1(1:10),[4 5]);
df_s1(:,1)=df_s1(:,1)/10;
df_s1(:,2)=df_s1(:,2)/5;
df_s2=Effects.num_short.shock(df_shock1(11:14),[4 5]);
df_s2(:,1)=df_s2(:,1)/10;
df_s2(:,2)=df_s2(:,2)/5;
df_f=Effects.num_short.shock(df_female,[4 5]);
df_f(:,1)=df_f(:,1)/10;
df_f(:,2)=df_f(:,2)/5;

sc_s1=Effects.num_short.scoposhk(scop_shock1,[4 5]);
sc_s1(:,1)=sc_s1(:,1)/10;
sc_s1(:,2)=sc_s1(:,2)/5;
sc_s2=Effects.num_short.scoposhk(scop_shock2,[4 5]);
sc_s2(:,1)=sc_s2(:,1)/10;
sc_s2(:,2)=sc_s2(:,2)/5;
sc_f=Effects.num_short.scoposhk(scop_female,[4 5]);
sc_f(:,1)=sc_f(:,1)/10;
sc_f(:,2)=sc_f(:,2)/5;

barr=Effects.num_short.barrier(barrier,[4 5]);
barr(:,1)=barr(:,1)/10;
barr(:,2)=barr(:,2)/5;
bar_f=Effects.num_short.barrier(bar_female,[4 5]);
bar_f(:,1)=bar_f(:,1)/10;
bar_f(:,2)=bar_f(:,2)/5;

so1=Effects.num_short.scopo(scop_only1,[4 5]);
so1(:,1)=so1(:,1)/10;
so1(:,2)=so1(:,2)/5;

subplot(3,4,[1 2]); %immediate retention plots
plot(df_s1','o-k'); hold on;
plot(df_s2','s-k');
plot(df_f','-c');
[p_df_imm, h]=signrank([df_s1(:,1); df_s2(:,1)],[df_s1(:,2); df_s2(:,2)])

plot(3:4,sc_s1','o-r'); hold on;
plot(3:4,sc_s2','s-r');
plot(3:4,sc_f','-c');
[p_sc_imm, h]=signrank([sc_s1(:,1); sc_s2(:,1)],[sc_s1(:,2); sc_s2(:,2)])

plot(5:6,barr','o-b');
plot(5:6,bar_f','-c');
[p_barr_imm, h]=signrank(barr(:,1));%,barr(:,2))

set(gca,'XLim',[0 9]); title(['df:' num2str(p_df_imm) ' sc:' num2str(p_df_imm)]);

subplot(3,4,9); %shk1 vs shk2 compairson for immediate retention scatter plot
 plot(1:2,[median(df_s1(:,2)./df_s1(:,1)) median(df_s2(:,2)./df_s2(:,1))],'o-k')%,[std(df_s1(:,2)./df_s1(:,1))/sqrt(size(df_s1,1)-1) std(df_s2(:,2)./df_s2(:,1))/sqrt(size(df_s2,1)-1)],'k');
 hold on;
 scatter(.85+0*(df_s1(:,2)./df_s1(:,1)),(df_s1(:,2)./df_s1(:,1)),'ok')%,[std(df_s1(:,2)./df_s1(:,1))/sqrt(size(df_s1,1)-1) std(df_s2(:,2)./df_s2(:,1))/sqrt(size(df_s2,1)-1)],'k');
 scatter(1.85+0*(df_s2(:,2)./df_s2(:,1)),(df_s2(:,2)./df_s2(:,1)),'sk')%,[std(df_s1(:,2)./df_s1(:,1))/sqrt(size(df_s1,1)-1) std(df_s2(:,2)./df_s2(:,1))/sqrt(size(df_s2,1)-1)],'k');
 plot(1:2,[median(sc_s1(:,2)./sc_s1(:,1)) median(sc_s2(:,2)./sc_s2(:,1))],'s-r');%,[std(sc_s1(:,2)./sc_s1(:,1))/sqrt(size(sc_s1,1)-1) std(sc_s2(:,2)./sc_s2(:,1))/sqrt(size(sc_s2,1)-1)],'r');
 scatter(1.15+0*(sc_s1(:,2)./sc_s1(:,1)),(sc_s1(:,2)./sc_s1(:,1)),'or')%,[std(sc_s1(:,2)./sc_s1(:,1))/sqrt(size(sc_s1,1)-1) std(sc_s2(:,2)./sc_s2(:,1))/sqrt(size(sc_s2,1)-1)],'k');
 scatter(2.15+0*(sc_s2(:,2)./sc_s2(:,1)),(sc_s2(:,2)./sc_s2(:,1)),'sr')%,[std(df_s1(:,2)./df_s1(:,1))/sqrt(size(df_s1,1)-1) std(df_s2(:,2)./df_s2(:,1))/sqrt(size(df_s2,1)-1)],'k');
set(gca,'XLim',[0 3],'YLim',[0 1.5]);

subplot(3,4,11); %shk1 vs shk2 compairson for immediate retention box plot
x=[(df_s1(:,2)./df_s1(:,1)); (sc_s1(:,2)./sc_s1(:,1)); (df_s2(:,2)./df_s2(:,1)); (sc_s2(:,2)./sc_s2(:,1))];
g=[(1+0*df_s1(:,2)./df_s1(:,1)); 2+0*(df_s2(:,2)./df_s2(:,1)); 3+0*(sc_s1(:,2)./sc_s1(:,1)); 4+0*(sc_s2(:,2)./sc_s2(:,1))];
boxplot(x,g);
set(gca,'XLim',[0 5],'YLim',[0 1.5]);

[p_sh1_v_sh2_imm, h]=ranksum((df_s1(:,2)./df_s1(:,1)),(df_s2(:,2)./df_s2(:,1)));

[p_drug_sday, h]=ranksum([(df_s1(:,2)./df_s1(:,1)); (df_s2(:,2)./df_s2(:,1))],[(sc_s1(:,2)./sc_s1(:,1)); (sc_s2(:,2)./sc_s2(:,1))]);
[p_shnum_sday, h]=ranksum([(df_s1(:,2)./df_s1(:,1)); (sc_s1(:,2)./sc_s1(:,1))],[(df_s2(:,2)./df_s2(:,1)); (sc_s2(:,2)./sc_s2(:,1))]);

[p_drug_sday p_shnum_sday]

subplot(3,4,[7 8]); hold off; %3-day extinction curves

df_s1=Effects.num_short.shock(df_shock1(1:10),[6:8])./Effects.num_short.shock(df_shock1(1:10),[3 3 3]);%/10
df_s2=Effects.num_short.shock(df_shock1(11:14),[6:8])./Effects.num_short.shock(df_shock1(11:14),[3 3 3]);
df_f=Effects.num_short.shock(df_female,[6:8])./Effects.num_short.shock(df_female,[3 3 3]);

sc_s1=Effects.num_short.scoposhk(scop_shock1,[6:8])./Effects.num_short.scoposhk(scop_shock1,[3 3 3]);
sc_s2=Effects.num_short.scoposhk(scop_shock2,[6:8])./Effects.num_short.scoposhk(scop_shock2,[3 3 3]);
sc_f=Effects.num_short.scoposhk(scop_female,[6:8])./Effects.num_short.scoposhk(scop_female,[3 3 3]);

barr=Effects.num_short.barrier(barrier,[6:8])./Effects.num_short.barrier(barrier,[3 3 3]);
bar_f=Effects.num_short.barrier(bar_female,[6:8])./Effects.num_short.barrier(bar_female,[3 3 3]);

plot(df_s1','c'); hold on;
plot(df_s2','c');
plot(df_f','m');
plot(nanmedian([df_s1; df_s2]),'o-k');

plot(5:7,sc_s1','c'); hold on;
plot(5:7,sc_s2','c');
plot(5:7,sc_f','m');
plot(5:7,nanmedian([sc_s1; sc_s2]),'o-r');

plot(9:11,barr','c'); hold on;
plot(9:11,bar_f','m'); hold on;
plot(9:11,nanmedian(barr),'o-b');

set(gca,'XLim',[0 12]);

subplot(3,4,5); hold off; %drug free cdf

extshk=[2 3 3 3 2 3 3 1 3 4 4 1 1 5 6]; 
extscopshk=[1 1 1 1 1 1 3 2 1 1]; 
extbarrier=[1 1 1 1 1 1];

cdfplot(extshk); set(gca,'XLim',[1 5]);

subplot(3,4,6); hold off; %scopolamine cdf
cdfplot(extscopshk); set(gca,'XLim',[1 5]);

p_df_imm
p_sc_imm
p_barr_imm

p_df_48
p_sc_48
p_barr_48

topdir

