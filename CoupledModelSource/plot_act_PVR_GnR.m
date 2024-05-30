function fig = plot_act_PVR_GnR(PVR_clinical_act, fig_num)

 fig = figure(fig_num);
 gamma_act_iter = 0: 100 / (length(PVR_clinical_act) - 1) : 100;
 plot(gamma_act_iter, PVR_clinical_act, 'Color', 'k', 'LineWidth',2.0)
 xlabel('Reduction in Active Tone (%)'); ylabel('Change in PVR (%)');
 set(gca, 'FontSize', 14, 'LineWidth', 2.0, 'FontWeight', 'bold')
 set(gcf,'Units','inches','Position',[.25 .25 12 8],'PaperUnits','inches','PaperPosition',[.25 .25 6 4])

end