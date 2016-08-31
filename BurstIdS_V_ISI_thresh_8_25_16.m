%BurstIdS.m
%Based on trigwIPSC6.m;
%Identify Bursts and their major parameters
% G.S.Cymbalyuk Dec 18 2002
%clear all;
%file_nameInp=input('Input filename:','s'); 

%spikeTime - time when spike appeared
%spikeNN - number in the series of decsrete time flow
%spikeV - maximum membrane potential during spike

function [BNSp,BDur,IntB] = BurstIdS_V_ISI_thresh_8_25_16(newspikeV,newspikeTime,newminSpikeV,newminSpikeTime,VBurstThresh,ISIThresh,ISIdiff)

% BNSp=Number of Spikes
% BDur=Burst Duration
% IntB=Interburst Interval

% BPerAv,BFreqAv,BDCAv,BNSpAv,BDurAv,IntBAv

% ISIburst=1; %Determin minimum interburst interval 
MinSpInBurst=2;%Minimum spikes which burst should have

NSpikes=length(newspikeTime);
ISI=newspikeTime(2:NSpikes)-newspikeTime(1:NSpikes-1);
length(ISI);
% % minspV = minspikeV(2:length(minspikeV));
% % ISpV = ISI(1:length(minspV));
% length(spikeV)
% length(minSpikeV)

if ISIdiff<0
    nB=find(newminSpikeV <= VBurstThresh )+1;%index of each inter burst
%     disp('ISI diff is negative')
else
nB=find(newminSpikeV <= VBurstThresh & ISI >= ISIThresh)+1;%index of each inter burst
end

% nB=find(newminSpikeV <= VBurstThresh)
% length(nB);
% nB1=find(ISI >= ISIThresh)+1
% length(nB1);
lennB=length(nB);

% NSpikes=length(spikeTime);
% ISI=spikeTime(2:NSpikes)-spikeTime(1:NSpikes-1);
% nB=find(ISI>=ISIburst);
% lennB=length(nB);

if(lennB<2)
    fprintf('Less than 2 bursts detected.\n');
     return;
end
nBI=nB(2:lennB)-nB(1:lennB-1);%Number of spikes in each burst
% lennB=length(nB);
NotBurstN=find(nBI<MinSpInBurst);
nB(NotBurstN)=[];
% lennB=length(nB);
nBF=nB+1;%Number for the first spike in a burst
nBL=nB(2:lennB)-1;%Number for the last spike in a burst
nBFln=length(nBF);
nBLln=length(nBL);
if nBFln>nBLln
    nBF(nBFln)=[];
end
% for i=1:nBLln
% BISI=spikeTime(nBF(i)+1:nBL(i))-spikeTime(nBF(i):nBL(i)-1);
% BFreq(i)=mean(1./BISI);
% end

BPer=newspikeTime(nB(2:lennB))-newspikeTime(nB(1:lennB-1));
BDur=newspikeTime(nBL)-newspikeTime(nBF);
% BDC=100.*BDur./BPer;
BNSp=nB(2:lennB)-nB(1:lennB-1);
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

% plot(spikeTime,spikeV,'m+');
% hold on;
% plot(spikeTime(nBF),spikeV(nBF),'m*');
% 
% hold on;
% plot(spikeTime(nBL),spikeV(nBL),'c*');
% 
% for i=1:lennB-1
% plot([spikeTime(nBF(i)) spikeTime(nBL(i))],[spikeV(nBF(i)) spikeV(nBF(i))],'-r+');
% end


% fprintf('Period= %6.3f std %6.3f\n',BPerAv,BPerStd);
% fprintf('Frequency= %6.3f std %6.3f\n',BFreqAv,BFreqStd);
% fprintf('DutyCycle= %6.3f std %6.3f\n',BDCAv,BDCStd);
% fprintf('Number of Spikes= %6.3f std %6.3f\n',BNSpAv,BNSpStd);
% fprintf('Burst Duration= %6.3f std %6.3f\n',BDurAv,BDurStd);
% fprintf('Interburst Interval = %6.3f std %6.3f\n',IntBAv,IntBStd);
% 
% fprintf('%d bursts had less than %d spikes and were skipped.\n',length(NotBurstN),MinSpInBurst);

clear nB nBI nBF nBL newspikeV newspikeTime newminSpikeV newminSpikeTime VBurstThresh ISIThresh ISIdiff NSpikes lennB BPer ISI 
% whos