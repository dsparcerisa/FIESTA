function isInvalid = planIsInvalid(plan, stageLimits)
isInvalid = false;

if numel(plan.X) ~= plan.numSpots
    isInvalid = true;
    return;
end

if numel(plan.Y) ~= plan.numSpots
    isInvalid = true;
    return;
end

if numel(plan.Z) ~= plan.numSpots
    isInvalid = true;
    return;
end

if strcmp(plan.mode, 'FLASH')
    if ~isfield(plan, 'tRendija') || plan.tRendija <= 0
        isInvalid = true;
        return;
    end
    if ~isfield(plan, 'Nshots') || numel(plan.Nshots) ~= plan.numSpots
        isInvalid = true;
        return;
    end    
elseif strcmp(plan.mode, 'CONV')
    if ~isfield(plan, 't_s') || numel(plan.t_s) ~= plan.numSpots
        isInvalid = true;
        return;
    end
else
    isInvalid = true;
    return;
end
   
if any(plan.X < stageLimits(1))
    isInvalid = true;
    return;
end

if any(plan.X > stageLimits(2))
    isInvalid = true;
    return;
end

if any(plan.Y < stageLimits(3))
    isInvalid = true;
    return;
end

if any(plan.Y > stageLimits(4))
    isInvalid = true;
    return;
end

if any(plan.Z < stageLimits(5))
    isInvalid = true;
    return;
end

if any(plan.Z > stageLimits(6))
    isInvalid = true;
    return;
end

end

