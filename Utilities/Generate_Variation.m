function [RR] = Generate_Variation(RMSSD, Constraint, TraceLength, BaseRR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates a fetal heart rate with bounded heart rate
% variability
% Based off of Jansz MS, Seed M, Van Amerom JFP, Wong D, Grosse-Wortmann L, Yoo SJ, MacGowan CK. Metric optimized gating for fetal cardiac MRI. Magnetic Resonance in Medicine. 2010;64:1304–1314.

% RMSSD         : RMS beat-to-beat differences (in ms) (6 is good)
% Constraint    : Multiplicative factor to bound HRV (0.14 is good)
% TraceLength   : Time of measurement (in seconds)
% BaseRR        : Base RR length
% RR            : Simulated RR intervals (CTG Trace)

% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear RR
RR(1) = BaseRR + RMSSD*randn(1,1);
while sum(RR) < TraceLength
   RR(end+1) = RR(end) + RMSSD*randn(1,1) - Constraint*(RR(end)-BaseRR);
end

end