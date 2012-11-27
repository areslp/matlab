function S = construct_similarity(C)
   n = length(C);
   S = zeros(n);
   classes = unique(C);
   row = 1;
   for i=1:length(classes)
       points(i) = sum(C==classes(i));
       S(row:row+points(i)-1,row:row+points(i)-1) = 1;
       row = row+points(i);
   end

       