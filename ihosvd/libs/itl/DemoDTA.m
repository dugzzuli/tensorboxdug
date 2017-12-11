clear all

A = full(ttensor(tenrand([2,3,4]),{rand(10,2),rand(30,3), rand(40,4)}));
B = full(ttensor(tenrand([2,3,4]),{rand(10,2),rand(30,3), rand(40,4)}));
D = {A,A,A,A,A,B,B,B,B,B};
[T,C] = DTA(D{1},[2,3,4]);
for i = 2:10
    [T,C] = DTA(D{i},[2,3,4],C);
    err = norm(full(T)-D{i})/norm(D{i});
    fprintf('tensor #%d has error %f\n',i,err);
end
