function [max,loc] = test( input )
  loc = 1;
  max = input(loc);
  for i = 2:length(input)
    if ( input(i) > max )
      max = input(i);
      loc = i;
    endif
  endfor
endfunction