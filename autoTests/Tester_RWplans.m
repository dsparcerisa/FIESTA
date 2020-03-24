clear all;
planFile = 'plan.txt';
planFile2 = 'plan2.txt';
planFile3 = 'plan3.txt';

plan = readPlan(planFile);
writePlan(plan, planFile2);
plan2 = readPlan(planFile2);
writePlan(plan2, planFile3);