function g = dftfilt2(f, H, classout)
%DFTFILT Performs frequency domain filtering.
%   G = DFTFILT(F, H) filters F in the frequency domain using the
%   filter transfer function H. The output, G, is the filtered
%   image, which has the same size as F.
%
%   Valid values of CLASSOUT are:
%
%   'original'  The output is of the same class as the input.
%               This is the default if CLASSOUT is not included
%               in the call.
%   'fltpoint'  The output is floating point of class single, unless
%               both f and H are of class double, in which case the
%               output also is of class double.
%
%   DFTFILT automatically pads f to be the same size as H. Both f
%   and H must be real. In addition, H must be an uncentered,
%   circularly-symmetric filter function. 

%   Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
%   Digital Image Processing Using MATLAB, Prentice-Hall, 2004
%   $Revision: 1.5 $  $Date: 2003/08/25 14:28:22 $

% Convert the input to floating point.
[f, revertclass] = tofloat(f);

% Obtain the FFT of the padded input.
F = fft2(f, size(H, 1), size(H, 2));

% Perform filtering. 
g = ifft2(H.*F);

% Crop to original size.
g = g(1:size(f, 1), 1:size(f, 2)); % g is of class single here.

% Convert the output to the same class as the input if so specified.
if nargin == 2 || strcmp(classout, 'original')
    g = revertclass(g);
elseif strcmp(classout, 'fltpoint')
    return
else
    error('Undefined class for the output image.')
end