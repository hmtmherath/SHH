% yy1 = zeros(9, 200000);  % Not necessary, but it's faster to preallocate space
% tt= zeros(9, 200000); 
clear all
% for k = 1 % hNa_5_10-3
% for k=2 % mH_5_10-2
%  for k=9 % mH_5_10-2
for k=1:10; % mH_5_10-2
    tic;
%    A=load([Voltage 'num2str(k,%d)'.mat]);
vFile = sprintf('Results-%d.txt',k);

    if exist(vFile, 'file')
        dat = importdata(vFile);
%         yy1(k,:) = V.yy1(1,:); %Voltage
           
    t=dat(:,1)*10^-3;
    V=dat(:,2);
    mH=dat(:,4);
    Noise=dat(:,6);
    
%     fig=figure;
%     plot(t,V)
%     xlabel('Time(s)')
%     ylabel('Voltage')
%     title('Voltage vs Time')
%     
%      fig=figure;
%     plot(t,mH)
%     xlabel('Time(s)')
%     ylabel('mH')
%     title('mH vs Time')
%     
%     fig=figure;
%     plot(t,Noise)
%     xlabel('Time(s)')
%     ylabel('Noise')
%     title('Noise vs Time')
%     
%     hist(Noise,100)

%     t=dat((433000:450000),1)*10^-3;
%     V=dat((433000:450000),2);

%     t=dat((853000:865000),1)*10^-3;
%     V=dat((853000:865000),2);
    
%     t=dat((306000:316000),1)*10^-3;
%     V=dat((306000:316000),2);
    
%     t=dat((320000:350000),1)*10^-3;
%     V=dat((320000:350000),2);
    
%     t=dat((210000:226000),1)*10^-3;
%     V=dat((210000:226000),2);

%     t=dat((55000:69000),1)*10^-3;
%     V=dat((55000:69000),2);
    
%      t=dat((1:20000),1)*10^-3;
%     V=dat((1:20000),2);

%     t=dat((12000:14500),1)*10^-3;
%     V=dat((12000:14500),2);
    
%      t=dat((300000:320000),1)*10^-3;
%     V=dat((300000:320000),2);
      
%     t=dat((1:20000),1)*10^-3;
%     V=dat((1:20000),2);
    
%     t=dat((735000:750000),1)*10^-3;
%     V=dat((735000:750000),2);
% 
%     hNa=dat(:,3);
%     mH=dat(:,4);
%     mK2=dat(:,5);
         
    

[newspikeV,newspikeTime,newminSpikeV,newminSpikeTime,VBurstThresh,ISIThresh,ISIdiff]=SpikeIdS_V_ISI_thresh_8_25_16(V,t);
[NS,BD,IB]=BurstIdS_V_ISI_thresh_8_25_16(newspikeV,newspikeTime,newminSpikeV,newminSpikeTime,VBurstThresh,ISIThresh,ISIdiff);
NSt{k}=NS;  %Number of Spikes
BDt{k}=BD;   %Burst Duration
IBt{k}=IB ;  %Interburst Interval

% % Properties of no.of spikes
NSAv(k)=mean(NSt{k}); % mean no. of spikes in each burst
NSstd(k)=std(NSt{k}); % strandard deviation of no. of spikes in each burst
NSCoefVar=NSstd./NSAv; % coefficient of variation=std/mean*100

BDAv(k)=mean(BDt{k}); % mean burst duration of bursts
BDstd(k)=std(BDt{k}); % strandard deviation of burst duration
BDCoefVar=BDstd./BDAv; % coefficient of variation=std/mean*100

% Properties of interburst interval
IBAv(k)=mean(IBt{k}); % mean interburst interval
IBstd(k)=std(IBt{k}); % strandard deviation of interburst interval
IBCoefVar=IBstd./IBAv; % coefficient of variation=std/mean*100


else
        fprintf('File %s does not exist.\n', vFile);
    end

end
% % % % storing all BDt,IBt cell arrays into one array
% % % totBD=[BDt{:}];
% % % totIB=[IBt{:}];
% % % totNS=[NSt{:}];
% % % 
% % % %save totBD, totIB
% % % save('totBDn.mat','totBD')
% % % save('totIBn.mat','totIB')
% % % save('totNSn.mat','totNS')
% % % 
% % %  
% remove zeros from BDAv array
BDAv=BDAv(BDAv~=0);
IBAv=IBAv(IBAv~=0);
%         
% % % %save BDAv, IBAv
% % % save('BDAvN.mat','BDAv')
% % % save('IBAvN.mat','IBAv')
% % % 
mBD=mean(BDAv)
mIB=mean(IBAv)

% remove NaN s from BDCoefVar array
BDCoefVar(isnan(BDCoefVar)) = [];
IBCoefVar(isnan(IBCoefVar)) = [];

mBDCoef=mean(BDCoefVar)*100
mIBCoef=mean(IBCoefVar)*100


clear newspikeV newspikeTime newminSpikeV newminSpikeTime VBurstThresh ISIThresh ISIdiff NS BD IB dat
