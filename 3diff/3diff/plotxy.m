exp = dlmread('explicitPoints.txt');
imp = dlmread('implicitPoints.txt');
%number of funcs
n = size(imp, 1) - 1;
%1 func
if n == 1
plot(exp(1,:), exp(2,:), imp(1,:), imp(2,:));
legend('explicit','implicit');
% 2 funcs
elseif n == 2
plot(exp(1,:), exp(2,:), exp(1,:), exp(3,:), imp(1,:), imp(2,:), imp(1,:), imp(3,:));
legend('explicit 1','explicit 2','implicit 1','implicit 2');
else
% 3 funcs
plot(exp(1,:), exp(2,:), exp(1,:), exp(3,:), exp(1,:), exp(4,:), imp(1,:), imp(2,:), imp(1,:), imp(3,:), imp(1,:), imp(4,:));
legend('explicit 1','explicit 2','explicit 3','implicit 1','implicit 2','implicit 3');
end
xlabel('x');
ylabel('y');