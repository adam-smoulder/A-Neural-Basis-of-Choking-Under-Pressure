function rtrans = fisherTransform(r)
% perform a fisher transform of input data
% Syntax is
%     rtrans = fisherTransform(r)
% -------------------------------------------------------
% O.Codol 9th Sept 2019
% email codol.olivier@gmail.com
% -------------------------------------------------------
rtrans = .5*log((1+r)./(1-r));
end