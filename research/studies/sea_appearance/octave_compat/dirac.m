function y = dirac(varargin)
% Numeric compatibility for HDR-VDP's CIE99 Fourier expression in Octave.
% Only zero and second derivatives occur. At nonzero frequencies both are
% zero; the singular zero-frequency bin is explicitly replaced by 1 in the
% official hdrvdp_otf_cie99.m. This function is not a symbolic distribution.
  if nargin == 1
    x = varargin{1};
  elseif nargin == 2 && varargin{1} == 2
    x = varargin{2};
  else
    error('Unsupported numeric Dirac use');
  end
  y = zeros(size(x));
  y(x==0) = Inf;
end
