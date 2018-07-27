%BurstIdS.m
% Edited by Thakshila M Herath Nov 20 2016

%spikeTime - time when spike appeared
%spikeNN - number in the series of decsrete time flow
%spikeV - maximum membrane potential during spike

function [BNSp,BDur,IntB,BPer,BISI,lenMinSp,tMinSp,nBF,nBL,spikeTime] = BurstIdS(spikeV,spikeTime,minSpikeV,minSpikeTime,VBurstThresh)

% BNSp=Number of Spikes
% BDur=Burst Duration
% IntB=Interburst Interval

% BPerAv,BFreqAv,BDCAv,BNSpAv,BDurAv,IntBAv

% ISIburst=1; %Determin minimum interburst interval 
MinSpInBurst=2;%Minimum spikes which burst should have

NSpikes=length(spikeTime);

nB=find(minSpikeV <= VBurstThresh );%Gives the index of last spike in a burst. index of each inter burst

lennB=length(nB)-1;

% for i=1:lennB
% hold on 
% plot(spikeTime(nB(i)),spikeV(nB(i)),'-r+')
% end

if(lennB<2)
    fprintf('Less than 2 bursts detected.\n');
     return;
end
nBI=nB(2:lennB)-nB(1:lennB-1);%Number of spikes in each burst

NotBurstN=find(nBI<MinSpInBurst);
tMinSp=spikeTime(nB(NotBurstN)+1)
lenMinSp=length(NotBurstN);
nBF=nB(1:lennB)+1; %Number for the first spike in a burst
nBF(NotBurstN)=[]; %nBF excluding bursts which have two spikes or less
nBI(NotBurstN)=[] ; %Number of spikes in the burst
lennBI=length(nBI);

nBL=nBF(1:lennBI)+nBI-1;  %Number for the last spike in a burst
nBFln=length(nBF);
nBLln=length(nBL);
% if nBFln>nBLln
% %     nBF(nBFln)
%     nBF(nBFln)=[];
% end

for i=1:nBLln
BISI{i}=spikeTime(nBF(i)+1:nBL(i))-spikeTime(nBF(i):nBL(i)-1);
% BISI=spikeTime(nBF(i)+1:nBL(i))-spikeTime(nBF(i):nBL(i)-1)
% BFreq(i)=mean(1./BISI);
end
% BISI
BPer=spikeTime(nBF(2:nBFln))-spikeTime(nBF(1:nBFln-1));
BDur=spikeTime(nBL)-spikeTime(nBF(1:nBLln));
% BDC=100.*BDur./BPer;
% BNSp=nBF(2:end)-nBF(1:end-1);
BNSp=nBI;
IntB=BPer-BDur;


% BFreqAv=mean(BFreq);
% BFreqStd=std(BFreq);
% BPerAv=mean(BPer);
% BPerStd=std(BPer);
% BDurAv=mean(BDur);
% BDurStd=std(BDur);
% IntBAv=mean(IntB);
% IntBStd=std(IntB);
% BDCAv=mean(BDC);
% BDCStd=std(BDC);
% BNSpAv=mean(BNSp);
% BNSpStd=std(BNSp);

% hold on
% plot(spikeTime,spikeV,'m+');
% hold on;
% plot(spikeTime(nBF),spikeV(nBF),'m*');
% % 
% hold on;
% plot(spikeTime(nBL),spikeV(nBL),'c*');
% % 
% hold on
% for i=1:nBFln-1
% plot([spikeTime(nBF(i)) spikeTime(nBL(i))],[spikeV(nBF(i)) spikeV(nBF(i))],'-r+');
% text(spikeTime(nBF(i)),spikeV(nBF(i)),['NS=',num2str(nBI(i))],'HorizontalAlignment','right','VerticalAlignment','top','FontSize',8)
% text(spikeTime(nBF(i)),spikeV(nBF(i)),['BD=',num2str(BDur(i))],'HorizontalAlignment','left','VerticalAlignment','top','FontSize',8)
% hold on
% end


% fprintf('Period= %6.3f std %6.3f\n',BPerAv,BPerStd);
% fprintf('Frequency= %6.3f std %6.3f\n',BFreqAv,BFreqStd);
% fprintf('DutyCycle= %6.3f std %6.3f\n',BDCAv,BDCStd);
% fprintf('Number of Spikes= %6.3f std %6.3f\n',BNSpAv,BNSpStd);
% fprintf('Burst Duration= %6.3f std %6.3f\n',BDurAv,BDurStd);
% fprintf('Interburst Interval = %6.3f std %6.3f\n',IntBAv,IntBStd);
% 
% fprintf('%d bursts had less than %d spikes and were skipped.\n',length(NotBurstN),MinSpInBurst);

clear nB nBI newspikeV newspikeTime newminSpikeV 
clear newminSpikeTime VBurstThresh ISIThresh ISIdiff NSpikes 
clear lennB ISI nBFln nBLln 
% whos
