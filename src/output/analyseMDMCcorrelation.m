function analyseMDMCcorrelation(filename)

[acceptX, acceptFOM, acceptSigma, acceptEpsilon] = analyseMDMCrun(filename);

start = 20;

subplot(3,2,1)
plot(acceptSigma(start:end), acceptEpsilon(start:end),'x');
xlabel('\sigma')
ylabel('\epsilon')
title(['starting from the ' int2str(start) ' accepted move'])
subplot(3,2,2)
plot(acceptX(start:end), acceptSigma(start:end),'x')
xlabel('Nmc')
ylabel('\sigma')
subplot(3,2,3)
plot(acceptEpsilon(start:end),acceptFOM(start:end),'x');
xlabel('\epsilon')
ylabel('FOM')
subplot(3,2,4)
plot(acceptSigma(start:end),acceptFOM(start:end),'x');
xlabel('\sigma')
ylabel('FOM')
subplot(3,2,5)
plot(acceptX(start:end), acceptEpsilon(start:end),'x')
xlabel('Nmc')
ylabel('\epsilon')
subplot(3,2,6)
plot(acceptX(start:end), acceptFOM(start:end),'x')
xlabel('Nmc')
ylabel('FOM')

[Y,I] = min(acceptFOM);
minAtNmc = acceptX(I)

covariance = [acceptFOM(start:end); acceptSigma(start:end); acceptEpsilon(start:end)]';
%cov(covariance)
corrcoef(covariance)
