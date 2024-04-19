CoM = readtable('com_forward_log.csv');
CoM = table2array(CoM);
delRow = any(CoM==0,2);
CoM(delRow,:) = [];
clear delRow

rightFoot = readtable('forward_right_foot_log.csv');
rightFoot = table2array(rightFoot);

leftFoot = readtable('forward_left_foot_log.csv');
leftFoot = table2array(leftFoot);
