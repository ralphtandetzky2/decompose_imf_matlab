function retval = dimf_punishing_term( gamma_coeffs )
  tau = gamma_coeffs(2:end) - gamma_coeffs(1:end-1);
  a = sum( max( 0, -imag(tau) ) );
  tau_heads = tau(1:end-1);
  tau_tails = tau(2:end  );
  tau_diff = abs( tau_tails - tau_heads );
  b = sum( max( 0, tau_diff - (imag(tau_heads)).^2 ) );
  c = sum( max( 0, tau_diff - (imag(tau_tails)).^2 ) );
  d = a+b+c;
  if ( d > 0 )
      d = d + 1;
  end
  retval = 1.e20*d;
  end