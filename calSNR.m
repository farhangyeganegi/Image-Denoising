% %% This function calculates the snr of a signal with reference to original 
% signal. SNR can be calculated for 1-D/2-D/3-D signals.
%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%
% %% output: snr-> snr in dB
%            
% %% input:  orgSig-> original 1-D/2-D/3-D signal (or reference signal)
%            recSig-> reconstructed (1-D/2-D/3-D) signal/ signal obtained 
%            from the experiment/ signal, of which snr is to be calculated 
%            with reference to original signal.
%            boun-> boun is the boundary left at the corners for the 
%            snr calculation.  default value = 0
%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%
function snr = calSNR(orgSig,recSig,varargin)

if isempty(varargin)
    boun = 0;
else boun = varargin{1};
end

if size(orgSig,2)==1       % if signal is 1-D
    orgSig = orgSig(boun+1:end-boun,:);
    recSig = recSig(boun+1:end-boun,:);
else                       % if signal is 2-D or 3-D
    orgSig = orgSig(boun+1:end-boun,boun+1:end-boun,:);
    recSig = recSig(boun+1:end-boun,boun+1:end-boun,:);
end
sigEner = norm(orgSig(:))^2;
errEner = norm(orgSig(:)-recSig(:))^2;
snr = 10*log10(sigEner/errEner);
end

