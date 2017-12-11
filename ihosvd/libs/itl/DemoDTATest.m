% 稠密张量分解函数 案例

clear A;
clear AC;
clear core_return;
clear U_return;
clear A_return;
A(: , : ,1)=[1,0,0;1,0,0;0,0,0];
A(: , : ,2)=[0,0,0;0,1,0;0,0,0];
A(: , : ,3)=[0,0,0;0,0,0;0,0,1];
A(: , : ,4)=[1,0,0;1,0,0;0,0,0];
A(: , : ,5)=[0,0,0;0,1,0;0,0,0];
A(: , : ,6)=[0,0,0;0,0,0;0,0,1];
B=A;
B=tensor(B)
[T,C]=DTA(B,[3,3,6]);


% 增加的张量
AC(: , : ,1)=[1,0,0;1,0,0;0,0,0];
AC(: , : ,2)=[0,0,0;0,1,0;0,0,0];
AC(: , : ,3)=[0,0,0;0,0,0;0,0,1];
AC(: , : ,4)=[1,0,0;1,0,0;0,0,0];
AC(: , : ,5)=[0,0,0;0,1,0;0,0,0];
AC(: , : ,6)=[0,0,0;0,0,0;0,0,1];
AC=tensor(AC);
[T,C]=DTA(AC,[3,3,6],C);
Ho=full(T);

A(: , : ,7)=[1,0,0;1,0,0;0,0,0];
A(: , : ,8)=[0,0,0;0,1,0;0,0,0];
A(: , : ,9)=[0,0,0;0,0,0;0,0,1];
A(: , : ,7)=[1,0,0;1,0,0;0,0,0];
A(: , : ,8)=[0,0,0;0,1,0;0,0,0];
A(: , : ,9)=[0,0,0;0,0,0;0,0,1];
A=tensor(A);
[core_return,U_return,A_return]=HOSVD_matlab_dense_modify(A,1);
A_return=tensor(A_return);
norm(A_return-Ho)


