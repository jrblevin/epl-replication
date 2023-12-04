% mc_main.m -- Control program for running Monte Carlo experiments.

clear
global selexper
global nobs

for selexper = 1:3
    for k = 0:1
        nobs = 1600*4^k;
        diary(sprintf('mc_epl_exper_%d_%d_obs.log', selexper, nobs));
        mc_epl
        diary off;
        outfile = sprintf('mc_epl_exper_%d_%d_obs.mat', selexper, nobs);
        save(outfile, 'nrepli', 'selexper', 'nobs', 'theta_true', 'kparam', 'namesb', 'bmat_1epl', 'bmat_1npl', 'bmat_cepl', 'bmat_cnpl', 'btrue_epl', 'btrue_npl', 'iter_cepl', 'iter_cnpl', 'time_cnpl','time_cepl', 'time_cnpl_iter','time_cepl_iter', 'fail_cnpl', 'fail_cepl', 'maxiter');
        parameter_histograms;
        time_histograms;
        iter_histograms;
    end
end
