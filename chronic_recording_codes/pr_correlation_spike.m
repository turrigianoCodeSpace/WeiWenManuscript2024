function [pSpike] = pr_correlation_spike(spikeV1,spikeV2)

%Equation from ATP 2019 PNAS
%two spike trains should have the same inner dimesion (to perform matrix
%multiplication)

pSpike = abs(mean((spikeV1-mean(spikeV1)) .* (spikeV2-mean(spikeV2))))/(std(spikeV1)*std(spikeV2));