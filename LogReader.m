CoM = readtable('com_log.csv');
CoM = table2array(CoM);
delRow = any(CoM==0,2);
CoM(delRow,:) = [];
clear delRow

rightFoot = readtable('right_foot_log.csv');
rightFoot = table2array(rightFoot);

leftFoot = readtable('left_foot_log.csv');
leftFoot = table2array(leftFoot);
