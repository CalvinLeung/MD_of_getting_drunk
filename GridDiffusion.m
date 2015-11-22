gridSize = 10000; % make this even
particles = 100;
wallGap = 1;
xPos = zeros(1,particles);
yPos = zeros(1,particles);
wallProfile = zeros(1,gridSize);
wallPos = gridSize/2;
timesteps = 5000;

% initialize the wall
for i = 1:gridSize
    if mod(i,2*wallGap)>=wallGap
        wallProfile(i) = 1; % 1 is a wall, 0 is a hole
    end
end

% initialize particle positions (indexed as 1 through "particles"):
for i = 1:particles
    loop = true;
    while loop
        randx = randi(wallPos-1);
        randy = randi(gridSize);
        % make sure no particle is there by checking them all
        loop = false; % default is that no particle is there
        for j = 1:particles
            if and((xPos(j) == randx), (yPos(j) == randy))
                loop = true; % try again
            end
        end
        xPos(i) = randx;
        yPos(i) = randy;
    end
end

% run timesteps
for t = 1:timesteps
    leftparticles = 0;
    rightparticles = 0;
    for p = 1:particles
        % only allow certain directions
        allowedDirs = [1,1,1,1]; % up, down, right, left
        if yPos(p) == gridSize
            allowedDirs(1) = 0; % can't go up
        elseif and((xPos(p) == wallPos), (wallProfile(yPos(p)+1) == 1))
            allowedDirs(1) = 0; % can't go up
        end
        if yPos(p) == 1
            allowedDirs(2) = 0; % can't go down
        elseif and((xPos(p) == wallPos), (wallProfile(yPos(p)-1) == 1))
            allowedDirs(2) = 0; % can't go down
        end
        if xPos(p) == gridSize
            allowedDirs(3) = 0; % can't go right
        elseif and((xPos(p) == wallPos-1), (wallProfile(yPos(p)) == 1))
            allowedDirs(3) = 0; % can't go right
        end
        if xPos(p) == 1
            allowedDirs(4) = 0; % can't go left
        elseif and((xPos(p) == wallPos+1), (wallProfile(yPos(p)) == 1))
            allowedDirs(4) = 0; % can't go left
        end
        
        % now update position
        loop = true;
        while loop
            dir = randi(4);
            if allowedDirs(dir) == 1
                loop = false; % we can go in this direction!
            end
        end
        if dir == 1
            yPos(p)=yPos(p)+1; % go up
        end
        if dir == 2
            yPos(p)=yPos(p)-1; % go down
        end
        if dir == 3
            xPos(p)=xPos(p)+1; % go right
        end
        if dir == 4
            xPos(p)=xPos(p)-1; % go left
        end
    end 
    %disp(t);
    %disp(leftparticles);
    %disp(rightparticles);
    end

% calculate grid matrix 
gridMatrix = zeros(gridSize);
for k = 1:gridSize
     gridMatrix(wallPos, k) = 2*wallProfile(k);
end
for p = 1:particles
    gridMatrix(xPos(p), yPos(p)) = 1;
    if xPos(p) >= wallPos
        rightparticles = rightparticles + 1;
    else
        leftparticles = leftparticles + 1;
    end
end
    
% output particle distribution
disp(leftparticles);
disp(rightparticles);
%pcolor(gridMatrix);
