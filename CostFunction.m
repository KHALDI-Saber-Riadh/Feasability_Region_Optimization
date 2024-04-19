function P = CostFunction(x)
Cost = @(x) FeasibilityRegion(x(1), x(2));

P = Cost(x)^2;
end