function retval = dimf_cost_off( gamma_coeffs, signal )
  t = ( (1:length(signal)) - 0.5 ) / length(signal);
  gamma = dimf_gamma( gamma_coeffs, t );
  offset = dimf_offset( gamma, signal );
  samples = dimf_samples( gamma + offset );
  retval = sumsqr( samples - signal ) + dimf_punishing_term( gamma_coeffs );
end