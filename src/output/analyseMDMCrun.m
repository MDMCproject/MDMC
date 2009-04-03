function analyseMDMCrun(filename)

s = xmlread(filename);

EleAccept = s.getElementsByTagName('accept');
EleRejected = s.getElementsByTagName('rejected');

nAccept = EleAccept.getLength();
acceptXaxis = zeros([1 nAccept]);
acceptYaxis = zeros([1 nAccept]);
acceptYaxisSigma = zeros([1 nAccept]);
acceptYaxisEpsilon = zeros([1 nAccept]);

nRejected = EleRejected.getLength();
rejectedXaxis = zeros([1 nRejected]);
rejectedYaxis = zeros([1 nRejected]);
rejectedYaxisSigma = zeros([1 nRejected]);
rejectedYaxisEpsilon = zeros([1 nRejected]);

docElement = s.getDocumentElement();
allNodes = docElement.getChildNodes();

iAcc = 1; iRej = 1; iMC = 1;
for i = 0 : allNodes.getLength()-1
  if ( allNodes.item(i).getNodeType() == 1 )
    if strcmp(allNodes.item(i).getNodeName(), 'accept')
      acceptXaxis(iAcc) = iMC;
      acceptYaxis(iAcc) = str2double(allNodes.item(i).getAttribute('val'));
      acceptYaxisSigma(iAcc) = str2double(allNodes.item(i).getAttribute('sigma'));
      acceptYaxisEpsilon(iAcc) = str2double(allNodes.item(i).getAttribute('epsilon'));
      iAcc = iAcc + 1; iMC = iMC + 1;
    end       
    if strcmp(allNodes.item(i).getNodeName(), 'rejected')
      rejectedXaxis(iRej) = iMC;
      rejectedYaxis(iRej) = str2double(allNodes.item(i).getAttribute('val'));
      rejectedYaxisSigma(iRej) = str2double(allNodes.item(i).getAttribute('sigma'));
      rejectedYaxisEpsilon(iRej) = str2double(allNodes.item(i).getAttribute('epsilon'));
      iRej = iRej + 1; iMC = iMC + 1;       
    end
  end
end
    
subplot(1,3,1)
plot(acceptXaxis, log10(acceptYaxis), 'bx')
hold on
plot(rejectedXaxis, log10(rejectedYaxis), 'rx')
title(['percentage accepted = ' num2str(100*nAccept/(nAccept+nRejected)) '%'])
ylabel('log10(FOM)')
xlabel('MC step number')
hold off
subplot(1,3,2)
plot(acceptXaxis, acceptYaxisSigma, 'bx')
hold on
plot(rejectedXaxis, rejectedYaxisSigma, 'rx')
ylabel('\sigma')
xlabel('MC step number')
hold off
subplot(1,3,3)
plot(acceptXaxis, acceptYaxisEpsilon, 'bx')
hold on
plot(rejectedXaxis, rejectedYaxisEpsilon, 'rx')
ylabel('\epsilon')
xlabel('MC step number')
hold off