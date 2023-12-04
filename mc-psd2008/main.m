% Estimates the Pesendorfer and Schmidt-Dengler (2008) model with
% several sample sizes and values of c.

nsims = 1000;

for eqm_dgp = [ 1, 2, 3 ]
    diary(sprintf('psd2008_estimate_eqm%d.log', eqm_dgp));
    for nobs = [ 250, 1000 ]

        % Generate data from the model
        generate_data

        % Estimate the model with consistent CCPs
        c = 0.0;
        nstart = 1;
        psd2008_estimate

        % Estimate the model with noisy CCPs
        c = 0.5;
        nstart = 1;
        psd2008_estimate

        % Estimate the model with estimated CCPs and 4 random starting values
        c = 0.0;
        nstart = 5;
        psd2008_estimate

        % Estimate the model with all 5 random starting values
        c = 1.0;
        nstart = 5;
        psd2008_estimate
    end
    diary off
end
