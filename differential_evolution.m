% A global optimization algorithm. 
% The differential evolution algorithm is implemented as described in Wikipedia
% on http://en.wikipedia.org/wiki/Differential_evolution, 29 Jan 2015. 
%
% INPUT PARAMETERS
%
% swarm ... A distribution of candidate solutions that will be improved.
%           It is given in form of a matrix where each row contains a 
%           candidate vector. The number of candidate vectors must be at 
%           least 4. 
% cost_func ... a function that can be applied to row vectors of swarm. 
% n_iters ... The number of iterations for improving the swarm. 
%             The total number of evaluations of the cost_function 
%             is O(rows(swarm)+n_iters).
% differential_weight ... a value in the interval (0,2].
% crossover_probability ... a value in the interval [0,1]. 
%
% OUTPUT PARAMETERS
%
% best_fit ... The swarm row with the lower cost after n_iters iterations.
% best_fit_cost ... The value of cost_func(best_fit_cost).   
function [best_fit,best_fit_cost] = differential_evolution( ...
    swarm, cost_func, n_iters, differential_weight, crossover_probability )
  N = size(swarm,1);
  if ( N < 4 )
    error( 'The swarm is too small for differential evolution.' );
  end
  costs = 1:N;
  for row = 1:N
    costs(row) = cost_func(swarm(row,:));
  end
  for iter = 1:n_iters
    xrow = mod(iter-1, N)+1;
    arow = xrow;
    while arow == xrow
      arow = unidrnd(N);
    end
    brow = xrow;
    while brow == xrow || brow == arow
      brow = unidrnd(N);
    end
    crow = xrow;
    while crow == xrow || crow == arow || crow == brow
      crow = unidrnd(N);
    end
    x = swarm(xrow,:);
    a = swarm(arow,:);
    b = swarm(brow,:);
    c = swarm(crow,:);
    R = unidrnd(length(x));
    z = a + differential_weight * (b-c);
%    coin = discrete_rnd( [0,1], ...
%      [1-crossover_probability,crossover_probability], 1, length(x) );
    y = x + (rand(size(x)) <= crossover_probability) .* (z-x);
    y(R) = z(R);
    ycost = cost_func(y);
    if ( ycost < costs(xrow) )
      swarm(xrow,:) = y;
      costs(xrow) = ycost;
    end
  end
  [best_fit_cost,min_idx] = min(costs);
  best_fit = swarm(min_idx,:);
end