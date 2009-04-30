[acceptX1, acceptFOM1, acceptSigma1, acceptEpsilon1, rejectX1, rejectFOM1, rejectSigma1, rejectEpsilon1] = ...
    analyseMDMCrun('mdmc_results_new_long.xml');
[acceptX2, acceptFOM2, acceptSigma2, acceptEpsilon2, rejectX2, rejectFOM2, rejectSigma2, rejectEpsilon2] = ...
    analyseMDMCrun('mdmc_results_start_from_good_sol.xml');

acceptX2 = acceptX2 + length(acceptX1) + length(rejectX1);
rejectX2 = rejectX2 + length(acceptX1) + length(rejectX1);

acceptX = [acceptX1  acceptX2];
rejectX = [rejectX1  rejectX2];
acceptFOM = [acceptFOM1  acceptFOM2];
rejectFOM = [rejectFOM1  rejectFOM2];

    plot(acceptX, log10(acceptFOM), 'bx')
    hold on
    plot(rejectX, log10(rejectFOM), 'ro')
    ylabel('log10(FOM)')
    xlabel('MC step')
    hold off