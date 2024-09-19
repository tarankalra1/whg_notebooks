function [f,Snn,A,ff,Af] = comp_spec(ts,sfreq)
% Based on xb_get_spectrum.m (XBeach toolbox), by Bas Hoonhout & Robert Mc Call
%
% % sampling freq in Hz!!
% 
% AvR, feb 2013
%%
% initialize spectrum
[n m] = size(ts);
 
% required number of samples
nr = 2^(nextpow2(n));

% number of Welch repetitions
nw = ceil((n-nr)/(0.5*nr))+1;

if nr > n
    nr = n;
end
% allocate matrices
Snn     = zeros(floor(nr/2),m);

% compute spectrum
idxe    = round(linspace(nr,n,nw));
idxb    = idxe-nr+1;

T       = nr/sfreq;
df      = 1/T;
ff      = df*[0:1:round(nr/2) -1*floor(nr/2)+1:1:-1];
f       = ff(1:floor(nr/2));

for i = 1:m
    P   = squeeze(ts(:,i));

    for j = 1:nw
        Pj  = P(idxb(j):idxe(j));
        
        Q   = fft(Pj,[],1)/nr;
        V   = 2/df*abs(Q).^2;
        Snn(:,i) = Snn(:,i) + squeeze(V(1:floor(nr/2)))/nw; % variance density spectrum
        
        % compute amplitude spectrum (following HOlthuijsen 2007)
        A(:,i) = sqrt(Snn(:,i)*df*2);
        Af(:,i) = sqrt(squeeze(V)/nw*df*2);
    end

end
