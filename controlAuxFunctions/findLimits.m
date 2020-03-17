function limits = findLimits(COM)
limits = nan(6,1);

fprintf(COM,'C,S1M6000,I1M-0,I1M400,R'); % Move to beginning of X and return 400 steps

while(true)
   [status, posX, posY, posZ] = monitorStatus(COM);
   if status==1
       limits(1) = posX;
       break
   end
end

fprintf(COM,'C,I1M0,I1M-400,R') % Move to maximum position of X

while(true)
   [status, posX, posY, posZ] = monitorStatus(COM);
   if status==1
       limits(2) = posX;
       break
   end
end

fprintf(COM,'C,S2M6000,I2M-0,I2M400,R'); % Move to beginning of Y and return 400 steps

while(true)
   [status, posX, posY, posZ] = monitorStatus(COM);
   if status==1
       limits(3) = posY;
       break
   end
end

fprintf(COM,'C,I2M0,I2M-400,R') % Move to maximum position of Y

while(true)
   [status, posX, posY, posZ] = monitorStatus(COM);
   if status==1
       limits(4) = posY;
       break
   end
end

fprintf(COM,'C,S3M6000,I3M-0,I3M400,R'); % Move to beginning of Z and return 400 steps

while(true)
   [status, posX, posY, posZ] = monitorStatus(COM);
   if status==1
       limits(5) = posZ;
       break
   end
end

fprintf(COM,'C,I3M0,I3M-400,R') % Move to maximum position of Z

while(true)
   [status, posX, posY, posZ] = monitorStatus(COM);
   if status==1
       limits(6) = posZ;
       break
   end
end

end

