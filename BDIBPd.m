function [RelNo,RunNo] = BDIBPd(vFile,k)

data=load(vFile);

    
    t1=data(:,1); V1=data(:,2); hNa1=data(:,3); mH1=data(:,4);mK21=data(:,5); 
    Noise=data(:,6);
    
   figure;plot(t1,V1);

%     delete 'InitialConditions.txt'
    fid=fopen('InitialConditions.txt','wt');
    fprintf(fid,'V=%.15f',V1(end));  
    fprintf(fid,'hNa=%.15f',hNa1(end));
    fprintf(fid,'mH=%.15f',mH1(end));
    fprintf(fid,'mK2=%.15f',mK21(end));
    fclose(fid);
    
[spikeV,spikeTime,minSpikeV,minSpikeTime,VBurstThresh]=SpikeIdS(V1,t1);
[NS,BD,IB,BPer,BISI,NMinSp,TMinSp,nBF,nBL,spikeTime]=BurstIdS(spikeV,spikeTime,minSpikeV,minSpikeTime,VBurstThresh);

NSt{k}=NS;  %Number of Spikes
BDt{k}=BD;   %Burst Duration
IBt{k}=IB ;  %Interburst Interval
Pdt{k}=BPer ; % Period
Bisi{k}=BISI; % Inter spike intervals within a burst
lnMinSp{k}=NMinSp; % No. of bursts which has less than 2 spikes
tMinSp{k}=TMinSp;

spiketf=spikeTime(nBF(1:end-1));
spiketl=spikeTime(nBL);
%size(spiketf)
%size(spiketl)

    for i=1:length(nBL)
    tidf(i)=find(t1==spiketf(i));
    tidl(i)=find(t1==spiketl(i));
    end
% hold on
% plot(spiketf,hNa(tidf),'r*') % first
% plot(spiketl,hNa(tidl),'b*') % last

    for j=1:length(tidf)
    hNabd(j)=mean(hNa1(tidf(j):tidl(j)));
    mHbd(j)=mean(mH1(tidf(j):tidl(j)));
    mK2bd(j)=mean(mK21(tidf(j):tidl(j)));
    Noisebd(j)=mean(Noise(tidf(j):tidl(j)));
    end
%     hNabd
%     mean(hNabd)
    hNaBD{k}=hNabd';
    mHBD{k}=mHbd';
    mK2BD{k}=mK2bd';
    NoiseBD{k}=Noisebd';
    
    for m=1:length(tidf)-1
    hNaib(m)=mean(hNa1(tidl(m):tidf(m+1)));
    hNaPd(m)=mean(hNa1(tidf(m):tidf(m+1)));
    mHib(m)=mean(mH1(tidl(m):tidf(m+1)));
    mHPd(m)=mean(mH1(tidf(m):tidf(m+1)));
    mK2ib(m)=mean(mK21(tidl(m):tidf(m+1)));
    mK2Pd(m)=mean(mK21(tidf(m):tidf(m+1)));
    Noiseib(m)=mean(Noise(tidl(m):tidf(m+1)));
    NoisePd(m)=mean(Noise(tidf(m):tidf(m+1)));
    end
    hNaIB{k}=hNaib';
    hNaPeriod{k}=hNaPd';
    mHIB{k}=mHib';
    mHPeriod{k}=mHPd';
    mK2IB{k}=mK2ib';
    mK2Period{k}=mK2Pd';
    NoiseIB{k}=Noiseib';
    NoisePeriod{k}=NoisePd';
    
%     mhNaBD{k}=mean(hNaBD);
%     mhNaIB{k}=mean(hNaIB);
%     mhNaPd{k}=mean(hNaPeriod);

%     mmHBD{k}=mean(mHBD);
%     mmHIB{k}=mean(mHIB);
%     mmHPd{k}=mean(mHPeriod);
% 
%     mmK2BD{k}=mean(mK2BD);
%     mmK2IB{k}=mean(mK2IB);
%     mmK2Pd{k}=mean(mK2Period);
% 
%     mNoiseBD{k}=mean(NoiseBD);
%     mNoiseIB{k}=mean(NoiseIB);
%     mNoisePd{k}=mean(NoisePeriod);
   
clear BISI BPer Files NMinSp TMinSp data filenames hNa 
clear mH  mK2 minSpikeV spiketf minSpikeTime  spikeV 
clear spikeTime tidl nBL tidf nBF spiketl Noise

% k=k+1;
RelNo=k;

% storing all BDt,IBt cell arrays into one array
totBD=cat(1,BDt{:});
totIB=cat(1,IBt{:});
totNS=cat(1,NSt{:});
totPd=cat(1,Pdt{:});
totBISI=cat(2,Bisi{:});
totMinSp=cat(1,lnMinSp{:});

%save totBD, totIB
save(['totBD_',datestr(now,'mm-dd-yy'),'.mat'],'totBD')
save(['totIB_',datestr(now,'mm-dd-yy'),'.mat'],'totIB')
save(['totNS_',datestr(now,'mm-dd-yy'),'.mat'],'totNS')
save(['totPd_',datestr(now,'mm-dd-yy'),'.mat'],'totPd')
save(['totBISI_',datestr(now,'mm-dd-yy'),'.mat'],'totBISI')
save(['totMinSp_',datestr(now,'mm-dd-yy'),'.mat'],'totMinSp')

% hNa
tothNaBD=cat(1,hNaBD{:});
tothNaIB=cat(1,hNaIB{:});
tothNaPeriod=cat(1,hNaPeriod{:});

% mBDhNa=mean(hNaBD)
% mIBhNa=mean(hNaIB)
% mPdhNa=mean(hNaPeriod)

save(['AVGhNaBDIBPd_',datestr(now,'mm-dd-yy'),'.mat'],'tothNaBD','tothNaIB','tothNaPeriod')

%mH
totmHBD=cat(1,mHBD{:});
totmHIB=cat(1,mHIB{:});
totmHPeriod=cat(1,mHPeriod{:});

% mBDmH=mean(mHBD)
% mIBmH=mean(mHIB)
% mPdmH=mean(mHPeriod)

save(['AVGmHBDIBPd_',datestr(now,'mm-dd-yy'),'.mat'],'totmHBD','totmHIB','totmHPeriod')

% mK2
totmK2BD=cat(1,mK2BD{:});
totmK2IB=cat(1,mK2IB{:});
totmK2Period=cat(1,mK2Period{:});

% mBDmK2=mean(mK2BD)
% mIBmK2=mean(mK2IB)
% mPdmK2=mean(mK2Period)

save(['AVGmK2BDIBPd_',datestr(now,'mm-dd-yy'),'.mat'],'totmK2BD','totmK2IB','totmK2Period')

% Noise
totNoiseBD=cat(1,NoiseBD{:});
totNoiseIB=cat(1,NoiseIB{:});
totNoisePeriod=cat(1,NoisePeriod{:});

% mBDNoise=mean(NoiseBD)
% mIBNoise=mean(NoiseIB)
% mPdNoise=mean(NoisePeriod)

save(['AVGNoiseBDIBPd_',datestr(now,'mm-dd-yy'),'.mat'],'totNoiseBD','totNoiseIB','totNoisePeriod')

RunNo=1;
RunNo=RunNo+1;

end