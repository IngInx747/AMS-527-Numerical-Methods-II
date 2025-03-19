# Track position during an optimal process
# Null input: reset the tracker
function tracker(x = [])

  global __xs = [];

  if isempty(x)
    __xs = [];
    return
  endif

  if isempty(__xs)
    __xs = x;
  else
    __xs = [__xs(:,:), x];
  endif

endfunction

