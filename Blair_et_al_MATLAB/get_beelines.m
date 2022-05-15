function data = get_beelines(data)

% RL_middleSpeed = [-118.81 -105.50 -107.16 -76.79 -98.62];
% LR_middleSpeed = [121.74 106.01 105.36 90.94 92.61];



%ms.time=ms.time(frameNum);

%stdth=2.5;                          % standard deviation parameter for spike detection
% frameint=mean(diff(ms.time));  % mean interframe interval for miniscope data in seconds
% temp = diff(ms.x'); allSpeed = [temp NaN]./frameint;

% compute spatial occupancy distributions (visit maps) on the track

spbinsize=13.33;    %size of each spatial bin in cm

%numcells=size(ms.C,1);

% figure(100); clf; 
% plot(ms.x,ms.time); hold on;

dx=diff(data.x(:));
SL=data.x(:)*0; %extract right-to-left trips
temp=find(data.x(1:end-1)>100 & data.x(2:end)<=100 & abs(data.y(1:end-1))<60);
temp=[1; temp(:); length(data.x(:))-1];
j = 0;
for i=2:(length(temp)-1)
   d1=temp(i)-find(dx(temp(i):-1:temp(i-1))>0,1,'first');
   if isempty(d1)
       d1=temp(i-1);
   end
   d2=temp(i)+find(dx(temp(i):temp(i+1))>0,1,'first')-1;
   if isempty(d2)
       d2=temp(i+1);
   end
   if data.x(d1)>100 && data.x(d2)<-100 % beeline
%        fit_x = [-140 0 140];
%        fit_y = [-40 abs(RL_middleSpeed(animalNum))/2 -40]; %10/9 changed -20 to -40 for both ends
%        f = fit(fit_x',fit_y','poly2');
%        check_spatial = data.x(d1:d2)';
%        check_speed = abs(allSpeed(d1:d2));
%        
%        checkRecorder=0;
%        for s=1:length(check_spatial)
%            if f(check_spatial(s)) > check_speed(s)
%                checkRecorder = checkRecorder +1;
%            end
%        end
%        
%        if checkRecorder == 0
           j=j+1;
           SL(d1:d2)=1;
           RL_index(j,1)=d1;
           RL_index(j,2)=d2;
           data.RLint(j,1)=d1;
           data.RLint(j,2)=d2;
%         end
%plot(data.x(d1:d2),data.time(d1:d2),'r');
   end
end

SL(isnan(SL))=0;
RLx2=data.x(:); RLx2(~SL)=NaN;
if exist('RL_index')
    percent_result.RL_beeline = size(RL_index,1);
    percent_result.RL_index = RL_index;
else
    percent_result.RL_beeline = 0;
    percent_result.RL_index = NaN;
end

SR=data.x(:)*0; %extract left-to-right trips
temp=find(data.x(1:end-1)<-100 & data.x(2:end)>=-100 & abs(data.y(1:end-1))<60);
temp=[1; temp(:); length(data.x(:))-1];
j = 0;
for i=2:(length(temp)-1)
   d1=temp(i)-find(dx(temp(i):-1:temp(i-1))<0,1,'first');
   if isempty(d1)
       d1=temp(i-1);
   end
   d2=temp(i)+find(dx(temp(i):temp(i+1))<0,1,'first')-1;
   if isempty(d2)
       d2=temp(i+1);
   end
   if data.x(d1)<-100 && data.x(d2)>100 % beeline
%        fit_x = [-140 0 140];
%        fit_y = [-40 abs(LR_middleSpeed(animalNum))/2 -40]; %10/9 changed -20 to -40 for both ends
%        f = fit(fit_x',fit_y','poly2');
%        check_spatial = data.x(d1:d2)';
%        check_speed = abs(allSpeed(d1:d2));
%        
%        checkRecorder=0;
%        for s=1:length(check_spatial)
%            if f(check_spatial(s)) > check_speed(s)
%                checkRecorder = checkRecorder +1;
%            end
%        end
%        
%        if checkRecorder == 0
           j=j+1;
           SR(d1:d2)=1;
           LR_index(j,1)=d1;
           LR_index(j,2)=d2;
           data.LRint(j,1)=d1;
           data.LRint(j,2)=d2;
%        end
   end
end


