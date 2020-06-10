function C = InterpMtrx(A, IntRat)


%% calculate interpolation
B = interp2(A, max(IntRat), 'cubic');
SizeB = size(B);
%figure;
%surf(B);
StepX =2^(max(IntRat) - IntRat(1));
StepY =2^(max(IntRat) - IntRat(2));
C = B(1:StepX:SizeB(1), 1:StepY:SizeB(2));

% A = rand(5,6);
% D = InterpMtrx(A, [5, 5]);
