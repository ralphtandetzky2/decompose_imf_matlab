function retval = dimf_samples( gamma )
  retval = exp( real(gamma) ) .* cos( imag(gamma) );
end
 