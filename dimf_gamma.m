% calculates
%     sum_{i=1}^N gamma_coeffs(i) * b_spline( (N-2)*t - i + 3 )
function retval = dimf_gamma( gamma_coeffs, t )
  N = length(gamma_coeffs);
  if ( N < 3 )
    error('gamma_coeffs must contain at least three items.');
  end
  x = (N-2)*t+3;
  retval = 0*t;
  for iter = 1:N
    retval = retval + gamma_coeffs(iter) * dimf_b_spline( x - iter );
  end
end
