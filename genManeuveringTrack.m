function [xTrue,T,numSamples]=genManeuveringTrack()
T=5;
TTotal=125+90+125+30+125;
%The extra sample is the very first step in the trajectory.
numSamples=TTotal/T+1;
xTrue=zeros(4,numSamples);

%An indicator of the maneuvering mode 1 for not maneuvering, and 2 for
%maneuvering.
mode=zeros(numSamples,1);

%The plane first flies westward at 120 m/s for 125 seconds.
xInit=[25e3;10e3; -120; 0];
xTrue(:,1)=xInit;
mode(1)=1;
numStepSeg=125/T;%The number of steps in this segment.

%The state transition matrix for constant velocity motion.
F=FPolyKal(T,4,1);

%Propagate the constant velocity model for 125 seconds.
baseStep=1;
for curStep=(baseStep+1):(baseStep+numStepSeg)
    xTrue(:,curStep)=F*xTrue(:,curStep-1);
    mode(curStep)=1;
end
baseStep=baseStep+numStepSeg;

%Next, the plane performs a 1 degree per second coordinated turn for 90
%seconds.
turnRate=1*(pi/180);%Convert from degrees per second to radians per second.

numStepSeg=90/T;
for curStep=(baseStep+1):(baseStep+numStepSeg)
    %Get the discrete-time propagation matrix for the state. The
    %FCoordTurn2D function assumes that the turn rate is part of the state.
    %In this instance, we just add it to the state, and then remove the
    %extra elements from the returned matrix.
    F=FCoordTurn2D(T,[xTrue(:,curStep-1);turnRate]);
    F=F(1:4,1:4);
    
    xTrue(:,curStep)=F*xTrue(:,curStep-1);
    mode(curStep)=2;
end
baseStep=baseStep+numStepSeg;

%Next, go straight (South after the previous turn) for another 125 seconds.
numStepSeg=125/T;%The number of steps in this segment.

%The state transition matrix for constant velocity motion.
F=FPolyKal(T,4,1);

for curStep=(baseStep+1):(baseStep+numStepSeg)
    xTrue(:,curStep)=F*xTrue(:,curStep-1);
    mode(curStep)=1;
end
baseStep=baseStep+numStepSeg;

%Make another turn. This time, it is at a rate of -three degrees per second
%for 30 seconds.
turnRate=-3*(pi/180);%Degrees per second to radians per second.

numStepSeg=30/T;%The number of steps in this segment.
for curStep=(baseStep+1):(baseStep+numStepSeg)
    F=FCoordTurn2D(T,[xTrue(:,curStep-1);turnRate]);
    F=F(1:4,1:4);
    
    xTrue(:,curStep)=F*xTrue(:,curStep-1);
    mode(curStep)=2;
end
baseStep=baseStep+numStepSeg;

%Finally, go straight again for another 125 seconds.
numStepSeg=125/T;%The number of steps in this segment.

%The state transition matrix for constant velocity motion.
F=FPolyKal(T,4,1);

for curStep=(baseStep+1):(baseStep+numStepSeg)
    xTrue(:,curStep)=F*xTrue(:,curStep-1);
    mode(curStep)=1;
end
end
