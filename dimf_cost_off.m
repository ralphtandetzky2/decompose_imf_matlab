function retval = dimf_cost_off( gamma_coeff, signal )
  t = ( (1:length(signal)) - 0.5 ) / length(signal);
  gamma = dimf_gamma( gamma_coeff, t );
  offset = dimf_offset( gamma, signal );
  samples = dimf_samples( gamma + offset );
  retval = sumsq( samples - signal ) + punishing_term( gamma_coeffs );
endfunction