function mat = antisym_mat(matIn)  %�������󷴶Գ���
    mat = [ 0,      -matIn(3),  matIn(2); 
          matIn(3),   0,      -matIn(1); 
          -matIn(2),  matIn(1),   0     ];   
end