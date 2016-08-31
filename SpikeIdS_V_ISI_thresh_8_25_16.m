%SpikeIdS.m
%Detects spikes by finding maximum of potentials above thresholds 

function [newspikeV,newspikeTime,newminSpikeV,newminSpikeTime,VBurstThresh,ISIThresh,ISIdiff] = SpikeIdS_V_ISI_thresh_8_25_16(yy1,Time)

V=yy1;%X(1,:);
% fig=figure;
% plot(Time,V)
% xlabel('Time(s)')
% ylabel('Voltage')
% title('Voltage vs Time')
    
[spikeV spikeId]=findpeaks(V);
spikeTime=Time(spikeId);
% hold on
% plot(spikeTime,spikeV,'r+');   

% (minimum V) Voltage and time realated to each spike
for j=1:length(spikeId)-1
[minSpikeV(j) minSpikeNaux(j)]=min(V(spikeId(j):spikeId(j+1)));
end
minSpikeID=minSpikeNaux(:)+spikeId(1:length(spikeId)-1)-1;
minSpikeTime=Time(minSpikeID);
minSpikeV=minSpikeV';
% hold on
% plot(minSpikeTime,minSpikeV,'g+');


% %Time difference
lenSpikeT=length(spikeTime);
TimeDiff1=spikeTime(2:lenSpikeT-1)-minSpikeTime(1:end-1);
TimeDiff2=minSpikeTime(2:end)-spikeTime(2:end-1);

% % % %Amplitude of spike
lenspikeV=length(spikeV);
lenminSpikeV=length(minSpikeV);

spikeV;
minSpikeV;

Amplitude1=spikeV(2:lenspikeV-1)-minSpikeV(1:lenminSpikeV-1);
Amplitude2=spikeV(2:lenspikeV-1)-minSpikeV(2:end);

lenTimeDiff1=length(TimeDiff1);
lenTimeDiff2=length(TimeDiff2);
lenAmp1=length(Amplitude1);
lenAmp2=length(Amplitude2);

% figure
% % scatter(VDiff,TimeDiff)
% % scatter(spikeTime(2:lenSpikeT),Amplitude)
% scatter(TimeDiff2,Amplitude2);
% xlabel('TimeDiff(s)')
% ylabel('Amplitude(V)')
% title('Amlitude vs TimeDiff')
% 
% 
badSpike=find(Amplitude1<0.0035| Amplitude2<0.02| TimeDiff1<0.02| TimeDiff1<0.02)+1;
% % Amplitude2<0.004|TimeDiff1<0.07
% % % badSpike=find(Amplitude<0.005)+1;
% % 
%     badSpikeT=spikeTime(badSpike);
%     for i=1:length(badSpike)
%     hold on
%     plot(badSpikeT(i,:)*ones(size(V)),V,'LineWidth', 1)
%     end


    newspikeV=spikeV;
%     newspikeV(badSpike)=[];
    newspikeTime=spikeTime;
%     newspikeTime(badSpike)=[];
%     newlenSpikeV=length(newspikeV);
%     
    newminSpikeV=minSpikeV;
%     newminSpikeV(badSpike)=[];  
    newminSpikeTime=minSpikeTime;
%     newminSpikeTime(badSpike)=[];
    
%     Amplitude1(badSpike(1:end-1))=[]
%     Amplitude2(badSpike(1:end-1))=[]
%      
%     TimeDiff1(badSpike(1:end-1))=[]
%     TimeDiff2(badSpike(1:end-1))=[]
%     
%     hold on
%     plot(newspikeTime,newspikeV,'bo');    
%     hold on
%     plot(newminSpikeTime,newminSpikeV,'ko');

for i = badSpike(end:-1:1)'
    if newspikeV(i-1) > newspikeV(i) && newspikeV(i+1) > newspikeV(i) 
        newspikeV(i)=[];
        newspikeTime(i)=[];
    elseif newspikeV(i-1) < newspikeV(i+1)
        newspikeV(i-1)=[];
        newspikeTime(i-1)=[];
    else
        newspikeV(i+1)=[];
        newspikeTime(i+1)=[];        
    end
end% % % %
for i = badSpike(end:-1:1)'
    if newminSpikeV(i-1) > newminSpikeV(i)
        newminSpikeV(i-1)=[];
        newminSpikeTime(i-1)=[];
    else
        newminSpikeV(i)=[];
        newminSpikeTime(i)=[];
    end
end
% % when you plot V vs t, you can uncomment this part
%     hold on
%     plot(newspikeTime,newspikeV,'bo');    
%     hold on
%     plot(newminSpikeTime,newminSpikeV,'ko');
    
%     newlenSpikeV=length(newspikeV);
    newlennewminSpikeV=length(newminSpikeV);

    % now we have correct spikes' voltages and minimum voltages. Then we
    % want to find the lowest voltages which are the minimum voltage of
    % the first spike of the each burst
    
    % finding new inter spike intervals
    for j= 1: length(newspikeTime)-1
    newISI(j)=newspikeTime(j+1)-newspikeTime(j);
    end
    
    newlenISI=length(newISI);
    
%      fig2=figure;
%     % scatter(newminSpikeV(2:length(newminSpikeV)-1),newISI(1:length(newISI)))
%     xlabel('minimum Voltage(V)')
%     ylabel('Inter spike interval(s)')
%     title('newISI vs minium Voltages') 

% %     % Identifying the threshold Voltage to detect the (lowest) voltage of the first spike in each burst. So we need to concentrate only on
% %     % minimum voltages
% %     sortminV=sort(newminSpikeV); % sort all minimum voltages in ascending order
% %     diffSortminV=diff(sortminV); % take the difference between adjacent elements
% %     maxDiffminV=max(diffSortminV); % choose the maximum difference which is the considerable gap
% %     minVl= sortminV(find(maxDiffminV==diffSortminV));% left voltage of the gap
% %     minVr=sortminV(find(maxDiffminV==diffSortminV)+1);  % right voltage of the gap
% %     VBurstThresh=(minVl+minVr)/2;
% % %     %plotting a line goes through middle of the gap, parallel to y-axis, on ISI vs minV plot
% %     hold on
% %     plot(VBurstThresh*ones(size(newISI)), newISI, 'LineWidth', 1)
% % 
% %     % threshold newISI 
% %     minVindexL=find(newminSpikeV<VBurstThresh); % finding the indicies of the minimum voltages which are voltages on left side of meanV
% %     newISIsL=newISI(minVindexL);          % finding the corresponding newISI
% %     lowestnewISIL=min(newISIsL);          % find the lowest newISI of that particular newISIs
% % 
% %     minVindexR=find(newminSpikeV>VBurstThresh) ;% finding the indicies of the minimum voltages which are voltages on right side of meanV
% %     newISIsR=newISI(minVindexR);          % finding the corresponding newISI
% %     maxnewISIR=max(newISIsR);           % find the maximum of that particular newISIs
% %     ISIThresh=(lowestnewISIL+maxnewISIR)/2;
% % 
% % %     %plotting a line goes through middle of the gap, parallel to x-axis, on ISI vs minV plot
% %     hold on
% %     plot(newminSpikeV,ISIThresh*ones(size(newminSpikeV)), 'LineWidth', 1)
% % 
% % 
% % N = length(newminSpikeV);
% % x = newminSpikeV;
% % y = newISI;
% % % Initialize a blue map
% % colorMap = [zeros(N, 1), zeros(N, 1), ones(N,1)];
% % % If y > 0.5, make the markers red.
% % for k = 1 : length(y)
% % 	if y(k) >= ISIThresh && x(k) < VBurstThresh
% % % if  x(k) < meanV
% % 		colorMap(k, :) = [1,0,0]; % Red
% %     else
% % 		colorMap(k, :) = [0,0,1]; % Blue
% % 	end
% % end
% % scatter(x,y,30* ones(length(y), 1), colorMap);
% % legend('Inter Spike Intervals(within bursts)','Inter Burst Intervals')
% % xlabel('minimum Voltage(V)')
% % ylabel('Inter spike interval(s)')
% % title('newISI vs minium Voltages')     
% % 
% % ISIdiff=-1;

    
 [cidx2,cmeans2] = kmeans(newminSpikeV,2,'dist','sqeuclidean');
% [silh2,h] = silhouette(newminSpikeV,cidx2,'sqeuclidean');
    
% ptsymb = {'bs','r^','md','go','c+'};
% for i = 1:2
%     clust = find(cidx2==i);
%     plot(newminSpikeV(clust),newISI(clust),ptsymb{i});
%     hold on
% end

cmeans2;
    clust1 = find(cidx2==1);
    clust2 = find(cidx2==2);
    
    if cmeans2(1)<cmeans2(2)
        clustISI=clust1;
        clustV=clust2;
    else
        clustISI=clust2;
        clustV=clust1;
    end   

     figure;
    plot(newminSpikeV(clustISI),newISI(clustISI),'bs');
    hold on   
    plot(newminSpikeV(clustV),newISI(clustV),'r^');
       
%     mean(clust1)
%     mean(clust2)
% Identifying the threshold Voltage to detect the (lowest) voltage of the first spike in each burst.   
minVl=max(newminSpikeV(clustISI));
minVr=min(newminSpikeV(clustV));
VBurstThresh=(minVl+minVr)/2;

hold on
plot(VBurstThresh*ones(size(newISI)), newISI, 'LineWidth', 1)

% Identifying thershold ISI
ISIup=min(newISI(clustISI));
ISIdown=max(newISI(clustV));
ISIdiff=(ISIup-ISIdown)/2;
ISIThresh=ISIdown+ISIdiff;

hold on
plot(newminSpikeV,ISIThresh*ones(size(newminSpikeV)), 'LineWidth', 1)

% if VThresh<0
% VBurstThresh=VThresh;
% ISIThresh=ISITh;
% else
%     VBurstThresh=ISITh;
%     ISIThresh=VThresh;
% end
 
% hold off
% figure
% xlabel('minimum Voltage(V)')
% ylabel('Inter spike interval(s)')
% title('newISI vs minium Voltages') 
% scatter(newminSpikeV,newISI)

% toc

clear yy1 V Time Amplitude1 Amplitude2 TimeDiff1 TimeDiff2 cidx2 cmeans2 clust1 clust2 clustV clustISI minVl minVr badSpike newISI ISIdown ISIup minSpikeID minSpikeNaux minSpikeTime minSpikeV spikeId spikeTime spikeV
% whos
% newspikeV,newspikeTime,newminSpikeV,newminSpikeTime,VBurstThresh,ISIThresh,ISIdiff
end
