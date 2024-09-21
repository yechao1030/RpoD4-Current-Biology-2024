load("Ptrc-rpoD4 dose response curve data.mat")

[hillCoeff ec50]=doseResponse(d,r);
xlabel('[IPTG] \muM','FontWeight','bold')
ylabel('Cell length (\mum)','FontWeight','bold')
grid on
set(gca,'FontSize',16)

xlim([0.1 10000])