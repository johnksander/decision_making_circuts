clear
clc
format compact

%do paul's model
PMweights = [.5:.02:1.5];
for idx = 1:numel(PMweights)
    behavior_test_PM(PMweights(idx))
end

%do my model
JKweights = [.5:.02:1.5];
for idx = 1:numel(JKweights)
    behavior_test_JK_noD(JKweights(idx))
end


behavior_test_JK_noD(1)

