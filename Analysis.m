clear all

for k=1:199
      
    clear NS BD IB lennB newspikeV newspikeTime newminSpikeV newminSpikeTime VBurstThresh ISIThresh ISIdiff dat V t
    tic;

vFile = sprintf('Results-%d.txt',k);
 
    if exist(vFile, 'file')

      [RelNo,RunNo] = BDIBPd(vFile,k)

      fprintf('%s is done.\n', vFile);
    else
      fprintf('File %s does not exist.\n', vFile);
    end
end