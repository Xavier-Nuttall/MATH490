maxDepth = 20;
barrierDepth = 5;
period = 5;
positiveIncidentAmplitude = 1;
negativeIncidentAmplitude = 0;


modeTruncation = 10;
truncation = 1000;
colocationPointCount = 50;


colocationRange = 5:5:400;
truncationRange = 25:25:2000;

Cs = zeros(length(colocationRange), length(truncationRange));
Bs = zeros(length(colocationRange), length(truncationRange));
ijs = strings(length(colocationRange), length(truncationRange));

for i = 1:length(colocationRange)
    for j = 1:length(truncationRange)
        colocationPointCount = colocationRange(i);
        truncation = truncationRange(j);

        [~, B, C, ~] = verticalBarrierSingleTruncation(maxDepth, barrierDepth, period, colocationPointCount, truncation, modeTruncation, positiveIncidentAmplitude, negativeIncidentAmplitude);

        Cs(i,j) = C(1);
        Bs(i,j) = B(1);
        ijs(i,j) = sprintf('i: %d, j: %d', colocationPointCount, truncation);
    end
end


CsError = abs(Cs(2:end,2:end)- Cs(1:end-1,2:end)/2 - Cs(2:end,1:end-1)/2);%./abs(Cs(2:end,2:end));
% Plot heatmap of CsError
figure(2)
clf
heatmap(truncationRange(2:end), colocationRange(2:end), log(CsError))


figure(1)
clf
hold on
Cdiff = abs(Cs(end,1:end-1) - Cs(end,2:end))./abs(Cs(end,1:end-1));
semilogy(truncationRange(1:end-1), Cdiff, 'DisplayName', 'C(1)')

% Cs(1,:)
xlabel('Truncation')
ylabel('Coefficient')
title('Coefficients for colocation 200')
hold off