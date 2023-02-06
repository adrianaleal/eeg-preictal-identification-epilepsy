function attractor = getAttractor(input_signal, tau, eDim, plotFigure)

% Input:
%   - input_signal: signal data
%   - tau: time delay
%   - eDim: embedding dimension
%   - plotFigure: flag to plot reconstructed phase space

% Output:
%   - attractor: reconstructed phase space


if nargin<4
    plotFigure = 0;
end

N = length(input_signal);

% Total points on the phase space:
T = N-(eDim-1)*tau;

% Initialize the phase space:
attractor = zeros(T,eDim);

for ii = 1:T
    attractor(ii,:) = input_signal(ii+(eDim-1)*tau-sort((0:eDim-1),'descend')*tau)';
end

if plotFigure
    if size(attractor,2)<4
        figure()
        if size(attractor,2)==3
            plot3(attractor(:,1),attractor(:,2),attractor(:,3),'-');
            zlabel('x(t+2$\tau$)');
        elseif size(attractor,2)==2
            plot(attractor(:,1),attractor(:,2),'-');
        end
        title('Signal time-delay embedding - state space plot');
        grid on;
        xlabel('x(t)');
        ylabel('x(t+$\tau$)');
    end
end

end
