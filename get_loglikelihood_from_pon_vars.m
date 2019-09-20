function pon_lml = get_loglikelihood_from_pon_vars(pon_vars, alt_count, ref_count)

sz = size(pon_vars, 1);
nsamp = sum(pon_vars(find(sum(pon_vars, 2) > 0, 1), :));

af_fracs=[0 .001 .003 .03 .2 1];
pon_var_like = NaN(sz, 5);

for a = 1:length(af_fracs),
  pon_var_like(:, a) = betacdf(af_fracs(a), alt_count + 1, ref_count + 1); 
end
pon_var_like = diff(pon_var_like, 1, 2);

pon_vars_cumprior = pon_vars(:, 3:7);
pon_vars_cumprior(:, 5) = pon_vars_cumprior(:, 5) + pon_vars(:, 8);
pon_vars_cumprior = fliplr(cumsum(fliplr(pon_vars_cumprior), 2))./nsamp;

pon_posterior = pon_var_like.*pon_vars_cumprior;

pon_lml = log10(sum(pon_posterior, 2) + 1e-20);
