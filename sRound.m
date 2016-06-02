function out = sRound(val, prec)
% sRound - smartRound
% Use this small helper function to round scalars, vectors or matrices at a  
% defined place after the decimal point. 
%
% Usage:    out = sRound(val, prec)
% 
% Example1: a = sRound(0.46,1);
%           a = 0.5000
% 
% Example2: b = sRound(-5.478,2);
%           b = -5.4800
% 
% Example3: c = sRound([0.8236 2.6958 3.3171 0.9502 0.0364], 3);
%           c = 0.8240    2.6960    3.3170    0.9500    0.0360    0.34330
% 
% Example4: d = sRound([-0.1694 -1.4351 0.3925; -0.0522 -2.6 1.3433], 3);
%           d = -0.1690   -1.4350    0.3930
%               -0.0520   -2.6000    1.3430
% Example5: e = sRound(3.5);
%           e = 4
% 
% See also: round, ceil, floor, fix

% Tobias Otto
% 1.2
% 28.04.2011

% 21.04.2011, Tobias: first draft
% 28.04.2011, Tobias: cosmetics, more user info, more error handling

%% Check input arguments
if(nargin ~= 2)
    prec = 0;
end

if(~isinteger(prec))
    disp('The second argument has to be an integer.')
    disp(['Your second argument: ' num2str(prec) ' is rounded to ' ...
        num2str(round(prec))]);
    prec = round(prec);
end

if(nargin == 0)
    error('sRound needs two arguments: out = sRound(val, precision)');
end

%% Prepare variables
res     = val - floor(val);
prec    = 10^prec;

%% Start computing
out     = round(res*prec)/prec;
out     = out + floor(val);
