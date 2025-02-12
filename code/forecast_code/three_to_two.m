function y=three_to_two(x)
y=[];
for z=1:size(x,3)
    y=[y;x(:,:,z)];
end