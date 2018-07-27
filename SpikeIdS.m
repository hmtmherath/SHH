%SpikeIdS.m
%Detects spikes by finding maximum of potentials above thresholds 

function [spikeV,spikeTime,minSpikeV,minSpikeTime,VBurstThresh] = SpikeIdS(yy1,Time)

V=yy1;%X(1,:);
t=Time;
% figure;
% plot(Time,V)
% xlabel('Time(s)')
% ylabel('Voltage')
% title('Voltage vs Time')
  
% finding peaks
[spikeV spikeId]=findpeaks(V);

% finding peaks' times
spikeTime=Time(spikeId);
% hold on
% plot(spikeTime,spikeV,'r+');   

% (minimum V) Voltage and time realated to each minima
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

% %Amplitude of spike
lenspikeV=length(spikeV);
lenminSpikeV=length(minSpikeV);

A1=spikeV(2:lenspikeV-1)-minSpikeV(1:lenminSpikeV-1); %Amplitude1
A2=spikeV(2:lenspikeV-1)-minSpikeV(2:end); % Amplitude2

% want to remove very small peaks in between two bursts. so identifying the
% threshold amplitude to remove those peaks

% figure;
% scatter(Amplitude1,Amplitude2)
% xlabel('Amplitude-rising(V)')
% ylabel('Amplitude-lowering(V)')
% title('Amlitude-lowering vs rising')

% scatter(Amplitude2,Amplitude1)
% xlabel('Amplitude-lowering(V)')
% ylabel('Amplitude-rising(V)')
% title('Amlitude- rising vs lowering')

% badPeaks=find(A1<0.01 | A2<0.01)+1;
% 
%     badPeaksT=spikeTime(badPeaks);
%     for i=1:length(badPeaks)
%     hold on
%     plot(badPeaksT(i,:)*ones(size(V)),V,'LineWidth', 1)
%     end

badPeaks=find(A1<0.01 | A2<0.01)+1;
    LenBadP=length(badPeaks);
if (LenBadP==0 )
%       hold on
%       plot(spikeTime,spikeV,'bo');

elseif LenBadP==1
    spikeTime(badPeaks)=[];
    spikeV(badPeaks)=[];
    spikeId(badPeaks)=[];

%     hold on
%     plot(spikeTime,spikeV,'bo');

else
    while(LenBadP>1)
        LenBadP;
        lenBadPks=LenBadP;
         i=1;
         badPeaks;
    %      clear badp
        while(lenBadPks>1)
           b=badPeaks(2)-badPeaks(1);
            if b==1
        %         disp(b);
        %         spikeV(badPeaks(1));
        %         spikeTime(badPeaks(1));
               % % removing double near by bad spikes/peaks
                  if spikeV(badPeaks(1))>spikeV(badPeaks(2))
        %                 spikeV(badPeaks(2))=0;
        %                 spikeTime(badPeaks(2))=0;  
                   badp(i)=badPeaks(2);
                    else
        %                 spikeV(badPeaks(1))=0;
        %                 spikeTime(badPeaks(1))=0;
                  badp(i)=badPeaks(1);
                  end
        %           badPeaks(1)
        %           badPeaks(2)
                    badPeaks([1,2])=[];
        %             badPeaks
                    lenBadPks=length(badPeaks);

        %             
                else
                % removing single bad spikes/peaks
        %         spikeV(badPeaks(1))
        %         spikeTime(badPeaks(1))

        %             spikeV(badPeaks(1))=0;
        %             spikeTime(badPeaks(1))=0;
                    badp(i)=badPeaks(1);
                    badPeaks(1)=[];
        %             badPeaks
                    lenBadPks=length(badPeaks);

            end
            i=i+1;
            badp;
        end


    spikeTime(badp)=[];
    spikeV(badp)=[];
    spikeId(badp)=[];
%     hold on
%     plot(spikeTime,spikeV,'bo');
    % spikeId=spikeId'
    clear minSpikeV minSpikeNaux
    % Finding new minima
    for j=1:length(spikeId)-1
    [minSpikeV(j) minSpikeNaux(j)]=min(V(spikeId(j):spikeId(j+1)));
    end
    % size(minSpikeNaux(:))
    % size(spikeId(1:length(spikeId)-1)-1)
    minSpikeID=minSpikeNaux(:)+spikeId(1:length(spikeId)-1)-1;
    minSpikeTime=Time(minSpikeID);
    minSpikeV=minSpikeV';
    % hold on
    % plot(minSpikeTime,minSpikeV,'g+');

    clear badp badPeaks A1 A2 lenspikeV lenminSpikeV
    % %Finiding new amplitudes of spikes
    lenspikeV=length(spikeV);
    lenminSpikeV=length(minSpikeV);

    A1=spikeV(2:lenspikeV-1)-minSpikeV(1:lenminSpikeV-1); %Amplitude1
    A2=spikeV(2:lenspikeV-1)-minSpikeV(2:end); % Amplitude2

    badPeaks=find(A1<0.01 | A2<0.01)+1;
    LenBadP=length(badPeaks);
    end

    LenBadP;
    badPeaks;

    spikeTime(badPeaks)=[];
    spikeV(badPeaks)=[];
    spikeId(badPeaks)=[];
%     minSpikeTime(badPeaks-1)=[];
%     minSpikeV(badPeaks-1)=[];
    
%     hold on
%     plot(spikeTime,spikeV,'bo');

end

% hold on
% plot(spikeTime,spikeV,'bo');   

 clear minSpikeV minSpikeNaux
    % Finding new minima
    for j=1:length(spikeId)-1
    [minSpikeV(j) minSpikeNaux(j)]=min(V(spikeId(j):spikeId(j+1)));
    end
    % size(minSpikeNaux(:))
    % size(spikeId(1:length(spikeId)-1)-1)
    minSpikeID=minSpikeNaux(:)+spikeId(1:length(spikeId)-1)-1;
    minSpikeTime=Time(minSpikeID);
    minSpikeV=minSpikeV';
    
% hold on
% plot(minSpikeTime,minSpikeV,'ko');

   
%     now we have correct spikes' voltages and minimum voltages. Then we
%     want to find the lowest voltages which are the minimum voltage of
%     the first spike of the each burst
    
%   finding new inter spike intervals
    for j= 1: length(spikeTime)-1
    ISI(j)=spikeTime(j+1)-spikeTime(j);
    end
%     
%     newlenISI=length(newISI);
    
%     length(newminSpikeV)
%     length(newminSpikeV(1:end-3))
%     length(newISI)
%     length(newISI(1:(end-3)))
         


% legend('Selected region','Inter Burst Intervals','Inter Spike Intervals(within bursts)','Location','east')
% xlabel('minimum Voltage(V)')
% ylabel('Inter spike interval(s)')
% title('newISI vs minium Voltages') 

 % clustering 'minSpikeV' data   
%  [cidx2,cmeans2] = kmeans(minSpikeV,2);
% minSpikeV
% ISI

% X=[minSpikeV ISI'];
% X=[minSpikeV(1:end-1) (ISI(1:end-1))'];

% figure;
% plot(X(:,1),X(:,2),'o')
% xlabel('minimum Voltage(V)')
% ylabel('Inter spike interval(s)')
% title('ISI vs minium Voltages') 
% 
% SortminV=sort(X(:,1));
% diffminV=diff(SortminV);
% maxDiff=max(diffminV);
% mxID=find(diffminV==maxDiff);
% SortminV(mxID);
% VBurstThresh=(SortminV(mxID)+SortminV(mxID+1))/2;
% % VBurstThresh=-0.0345;
% % if minSpikeV==VBurstThresh
% %     exit
% % end
% hold on
% % % plot(SortminV(mxID)*ones(size(X(:,2))), X(:,2), 'LineWidth', 1)
% % % plot(SortminV(mxID+1)*ones(size(X(:,2))), X(:,2), 'LineWidth', 1)
% plot(VBurstThresh*ones(size(X(:,2))), X(:,2), 'LineWidth', 1)
% %  plot(Time,VBurstThresh*ones(size(Time)), 'LineWidth', 1)

 

% Threshold should lie below zero. So we can choose all minV s less than
% zero for both clusters
% minSpikeV=minSpikeV(minSpikeV<0);
X=[minSpikeV(1:end-1) (ISI(1:end-1))'];
XV=X(:,1); XI=X(:,2);
[cidx,cmeans] = kmeans(X,2, 'Start', [min(XV),mean(XI); mean(XV),min(XI)]);

% Identifying the corresponding indicies for cluster 1 and cluster 2
    clust1 = find(cidx==1);
    clust2 = find(cidx==2);
    
% Finding minV gap for each cluster
    clust1minVdiff=max(X(clust1,1))-min(X(clust1,1));
    clust2minVdiff=max(X(clust2,1))-min(X(clust2,1));
    
    max(X(clust1,1));
    min(X(clust1,1));
    cmeans(1);
    cmeans(2);

% Checking the centroids and gap of minV. Initially assuming, clust1= clustV and clust2=clustISI    
    if cmeans(1)<cmeans(2) & clust1minVdiff<clust2minVdiff
        clustISI=clust1;
        clustV=clust2;
    else cmeans(1)>cmeans(2) & clust1minVdiff>clust2minVdiff
        disp('Clust1= Clust ISI');
        clustISI=clust2;
        clustV=clust1;
    end   

%     figure;
% %     hold on  
% %     plot(cmeans(:,1),cmeans(:,2),'kx','MarkerSize',15,'LineWidth',3)
% % %     plot(min(minSpikeV),max(ISI),'kx', 'MarkerSize',15,'LineWidth',3)

%     plot(X(clustISI,1),X(clustISI,2),'bs');
%     hold on   
%     plot(X(clustV,1),X(clustV,2),'r^');
       
% %     mean(clust1);
% %     mean(clust2);

% Identifying the threshold Voltage to detect the (lowest) voltage of the first spike in each burst. 
minVISI=XV(clustISI); % minimum Voltages corresponding to Clust ISI
maxVcI=max(minVISI); % maximum voltage of clust ISI

minVV=XV(clustV); % minimum Voltages corresponding to Clust V
minVcV=min(minVV); % minimum voltage of clust V

% hold on
% plot(maxVcI*ones(size(X(:,2))), X(:,2), 'LineWidth', 1)
% plot(minVcV*ones(size(X(:,2))), X(:,2), 'LineWidth', 1)

cmeans(1);
cmeans(2);
meanC=(cmeans(1)+cmeans(2))/2;

if (maxVcI>minVcV)  
    if (maxVcI>meanC && minVcV<meanC) % minVcV and maxVcI have falsely detected.
        VcI=minVISI(minVISI<meanC);
        VcV=minVV(minVV>meanC);
        maxVcI=max(VcI);
        minVcV=min(VcV);
        fprintf('maxVcI>meanC && minVcV<meanC');
        
    elseif (maxVcI>meanC) % maxVcI has falsely detected and lie inside clust V
        VcI=minVISI(minVISI<meanC);
        maxVcI=max(VcI);
        fprintf('maxVcI>meanC');
    
    else (minVcV<meanC) % minVcV has falsely detected and lie inside clust ISI
        maxVcI;
        VcV=minVV(minVV>meanC);
        minVcV=min(VcV);
        fprintf('minVcV<meanC');
         
    end
         
    VBurstTh=(maxVcI+minVcV)/2;
    if (VBurstTh>minVcV || VBurstTh<maxVcI)
        fprintf('Threshold voltage is not deteced correctly\n');
        return;
    else
        VBurstThresh=VBurstTh;
    end

       
else
VBurstThresh=(maxVcI+minVcV)/2;
end

% % hold on
% plot(maxVcI*ones(size(X(:,2))), X(:,2), 'LineWidth', 1)
% plot(minVcV*ones(size(X(:,2))), X(:,2), 'LineWidth', 1)

% hold on
% plot(VBurstThresh*ones(size(X(:,2))), X(:,2), 'LineWidth', 1)
% % % plot(Time,VBurstThresh*ones(size(Time)), 'LineWidth', 1)
% 
%  % figure
% xlabel('minimum Voltage(V)')
% ylabel('Inter spike interval(s)')
% title('ISI vs minium Voltages') 
% legend('Inter Burst Intervals','Inter Spike Intervals(within bursts)','Location','east')


end

