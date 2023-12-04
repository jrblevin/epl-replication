nbins = 100;
hist_edges = linspace(0, maxiter, nbins);


figure
t = tiledlayout(2,1);

% EPL
ax1 = nexttile;
h1 = histogram(iter_cepl, hist_edges, 'Normalization', 'pdf');
set(h1, 'FaceColor',[0.4 0.4 0.4]);
if selexper == 1
    legend('EPL', 'location', 'northeast')
else
    legend('EPL', 'location', 'northwest')
end

% NPL
ax2 = nexttile;
h2 = histogram(iter_cnpl, hist_edges, 'Normalization', 'pdf');
set(gca, 'Ydir', 'reverse')
set(h2, 'FaceColor',[0.9 0.9 0.9]);
if selexper == 1
    legend('NPL', 'location', 'southeast')
else
    legend('NPL', 'location', 'southwest')
end

% Link and move plots closer together
linkaxes([ax1,ax2],'x');
xticklabels(ax1,{})
t.TileSpacing = 'none';
set(ax1, 'XTick', [], 'YTick', []);
set(ax2, 'YTick', []);

% Fix y axis
ylim1 = ylim(ax1);
ylim2 = ylim(ax2);
ymax = max(ylim1(2), ylim2(2));
ylim(ax1, [0, ymax]);
ylim(ax2, [0, ymax]);

saveas(gcf, sprintf('mc_epl_exper_%d_%d_obs_iter_hist.pdf', selexper, nobs));

close all
