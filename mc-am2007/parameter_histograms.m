nbins = 100;
hist_edges = {}
hist_edges{1} = linspace(-2.5, -1, nbins);
hist_edges{2} = linspace(-2.5, -1, nbins);
hist_edges{3} = linspace(-2.5, -1, nbins);
hist_edges{4} = linspace(-2.5, -1, nbins);
hist_edges{5} = linspace(-2.5, -1, nbins);
hist_edges{6} = linspace(0.5, 1.5, nbins);
if selexper == 1
    hist_edges{7} = linspace(0, 2, nbins);
elseif selexper == 2
    hist_edges{7} = linspace(1, 4, nbins);
elseif selexper == 3
    hist_edges{7} = linspace(2, 5, nbins);
end
hist_edges{8} = linspace(0.7, 1.4, nbins);

for k = 1:kparam
    figure
    t = tiledlayout(2,1);

    % EPL
    ax1 = nexttile;
    h1 = histogram(bmat_cepl(:,k), hist_edges{k}, 'Normalization', 'pdf');
    xl = xline(theta_true(k),'-','LineWidth',2);
    set(h1, 'FaceColor',[0.4 0.4 0.4]);
    legend('EPL', 'location', 'northwest')

    % NPL
    ax2 = nexttile;
    h2 = histogram(bmat_cnpl(:,k), hist_edges{k}, 'Normalization', 'pdf');
    xl = xline(theta_true(k),'-','LineWidth',2);
    set(gca, 'Ydir', 'reverse')
    set(h2, 'FaceColor',[0.9 0.9 0.9]);
    legend('NPL', 'location', 'southwest')

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

    saveas(gcf, sprintf('mc_epl_exper_%d_%d_obs_param_%d_hist.pdf', selexper, nobs, k));

end
close all

