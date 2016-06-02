function confusion_matrix = confusionmatrix(class_matrix)

% class_matrix: An CxM matrix, where M is the number of test samples in the largest class, and C is the number of classes 
%               Each entry is a number to which the test sample belongs, (e.g. 1 for class1, 2 for class 2, etc).  
%               This matrix is often the "memory" output variable of the leaveoneout classifier m-files.

[C,unused] = size(class_matrix);
confusion_matrix = zeros(C+1);

%count number of classifications and fill in entries of confusion matrix
for k = 1 : C
    for m = 1 : C
        confusion_matrix(k,m) = length(find(class_matrix(k,:)==m));
    end
end

%compute producer accuracies (across rows) [sensitivity]
for k = 1 : C
    confusion_matrix(k,C+1) = confusion_matrix(k,k)/sum(confusion_matrix(k,1:C));
end

%compute user accuracies (down columns)
for k = 1 : C
    confusion_matrix(C+1,k) = confusion_matrix(k,k)/sum(confusion_matrix(1:C,k));
end

%compute overall accuracy
confusion_matrix(C+1,C+1) = sum(diag(confusion_matrix(1:C,1:C)))/sum(sum(confusion_matrix(1:C,1:C)));

%compute specificity
for k = 1 : C
    numerator=zeros(C,1);
    for m = 1 : C
        if k==m;
            continue
        else
        numerator(m) = sum(confusion_matrix([1:k-1,k+1:C],m))-confusion_matrix(k,m);
        end
    end

    confusion_matrix(k,C+2) = sum(numerator)/(sum(sum(confusion_matrix(1:C,1:C)))-sum(confusion_matrix(k,1:C)));
end




