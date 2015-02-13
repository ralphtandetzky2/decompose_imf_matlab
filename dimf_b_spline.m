function retval = dimf_b_spline( t )
  t  = min( t, 3-t );
  t1 = max( t,0 );
  t2 = max( t-1,0 );
  retval = 0.5*t1.*t1 - 1.5*t2.*t2;
end
