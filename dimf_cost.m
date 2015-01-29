function retval = dimf_cost( gamma_coeff, signal )
  t = ( (1:length(signal)) - 0.5 ) / length(signal);
  gamma = dimf_gamma( gamma_coeff, t );
  samples = dimf_samples( gamma );
  retval = sumsq( samples - signal ) + punishing_term( gamma_coeffs );
endfunction