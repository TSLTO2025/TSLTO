function [preci, recall, f1] = anomaly_test(trueR,preR)
real_units = nnz(trueR);
detected_units = nnz(preR);
rd_units = nnz(trueR.*preR);
preci = rd_units/detected_units;
recall = rd_units/real_units;
f1 = 2*preci*recall/(preci+recall);
fprintf('Precision = %f\n Recall = %f\n F1 = %f\n',preci,recall,f1);