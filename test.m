function [max,loc] = test( input )
  loc = 1;
  max = input(loc);
  for iter = 2:length(input)
    if ( input(iter) > max )
      max = input(iter);
      loc = iter;
    end
  end
end