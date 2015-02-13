% finds an optimal value `offset` such that dimf_samples(gamma+offset)
% is closest to `signal` in an l2 sense. This means
%     sumsq( dimf_samples(gamma+offset) - signal )
% is minimized.
function retval = dimf_offset( gamma, signal )
  samples_complex = exp( gamma );
  if ( any(size(gamma) ~= size(signal)) )
    error( 'The sizes of the input parameters gamma and signal must be equal.' );
  end
  AB = signal / [real(samples_complex); imag(samples_complex)];
  A = AB(1);
  B = AB(2);
  retval = log(A-1i*B);
end