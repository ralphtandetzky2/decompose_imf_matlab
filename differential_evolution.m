% The differential evolution algorithm as described on Wikipedia
% on http://en.wikipedia.org/wiki/Differential_evolution, 29 Jan 2015. 
function retval = differential_evolution( 
    swarm, cost_func, n_iters, differential_weight, crossover_probability )
  N = rows(swarm);
  if ( N < 4 )
    error( "The swarm is too small for differential evolution." );
  endif
  costs = 1:N;
  for row = 1:N
    costs(row) = cost_func(swarm(row,:));
  endfor
  for iter = 1:n_iters
    xrow = mod(iter-1, N)+1;
    do arow = unidrnd(N); until (  arow != xrow );
    do brow = unidrnd(N); until ( (brow != xrow) & (brow != arow) );
    do crow = unidrnd(N); until ( (crow != xrow) & (crow != brow) & (crow != arow) );
    x = swarm(xrow,:);
    a = swarm(arow,:);
    b = swarm(brow,:);
    c = swarm(crow,:);
    R = unidrnd(length(x));
    z = a + differential_weight * (b-c);
    coin = discrete_rnd( [0,1], [1-crossover_probability,crossover_probability], 1, length(x) );
    y = x + coin.*(z-x);
    y(R) = z(R);
    ycost = cost_func(y);
    if ( ycost < costs(xrow) )
      swarm(xrow,:) = y;
      costs(xrow) = ycost;
    endif
  endfor
  [min_cost,min_idx] = min(costs);
  retval = swarm(min_idx,:);
endfunction