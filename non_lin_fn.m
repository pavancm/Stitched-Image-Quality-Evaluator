function wt = non_lin_fn(wt,sig)

wt = 1 - exp(-(wt/sig).^2);
