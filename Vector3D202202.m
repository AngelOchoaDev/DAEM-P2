function [outputArg1,outputArg2] = Vector3D202202(inputArg1, inputArg2)
%VECTOR3D202202 Summary of this function goes here
%   Detailed explanation goes here
Q=(0:10:360)';
x=inputArg1{1,2}*cosd(Q)/2;
y=inputArg1{1,2}*sind(Q)/2;
z=ones(size(Q));
hold on;
if strcmp(inputArg1{1,3},'z')==1
    surf(inputArg1{1,1}*[x x 2*x 0*x] + inputArg2(1,1),...
        inputArg1{1,1}*[y y 2*y 0*y] + inputArg2(1,2),...
        inputArg1{1,1}*[0*z 2*z/3 2*z/3 z] + inputArg2(1,3));
elseif strcmp(inputArg1{1,3},'x')==1
    surf(inputArg1{1,1}*[0*z 2*z/3 2*z/3 z] + inputArg2(1,1),...
        inputArg1{1,1}*[y y 2*y 0*y] + inputArg2(1,2),...
        inputArg1{1,1}*[x x 2*x 0*x] + inputArg2(1,3));
elseif strcmp(inputArg1{1,3},'y')==1
    surf(inputArg1{1,1}*[x x 2*x 0*x] + inputArg2(1,1),...
        inputArg1{1,1}*[0*z 2*z/3 2*z/3 z] + inputArg2(1,2),...
        inputArg1{1,1}*[y y 2*y 0*y] + inputArg2(1,3));
end
axis xy; grid on; view(130,30);
shading interp; colormap('jet');

