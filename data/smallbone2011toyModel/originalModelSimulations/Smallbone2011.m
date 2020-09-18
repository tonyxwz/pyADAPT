function xdot = Smallbone2011(time, x_values,p)
% function vanHeerden1 takes
%
% either	1) no arguments
%       	    and returns a vector of the initial values
%
% or    	2) time - the elapsed time since the beginning of the reactions
%       	   x_values    - vector of the current values of the variables
%       	    and returns a vector of the rate of change of value of each of the variables
%
% vanHeerden1 can be used with MATLABs odeN functions as 
%
%	[t,x] = ode23(@vanHeerden1, [0, t_end], vanHeerden1)
%
%			where  t_end is the end time of the simulation
%
%The variables in this model are related to the output vectors with the following indices
%	Index	Variable name
%	  1  	  Glc (glucose_in)
%	  2  	  G1P
%	  3  	  G6P
%	  4  	  Trehalose
%	  5  	  T6P
%	  6  	  UDP_glucose
%   Index   Fixed variables
%	  7  	  ADP
%	  8  	  ATP
%	  9   	  diphosphate
%     10      F6P
%	  11  	  H+
%	  12  	  phosphate
%	  13  	  UDP
%	  14  	  UTP
%	  15  	  water
%	  16  	  glucose_ext
%
%--------------------------------------------------------
% output vector

xdot = zeros(6, 1);

%--------------------------------------------------------
% compartment values

compartment = 1;

%--------------------------------------------------------
% parameter values

% inside 'setParameterStructure_smallbone2011.m'

%--------------------------------------------------------
% initial values of variables - these may be overridden by assignment rules
% NOTE: any use of initialAssignments has been considered in calculating the initial values

if (nargin == 0)

	% initial time
	time = 0;
    
    % initial values
    Glc         = 0.09765; %mM
    G1P         = 0.1; %mM
    G6P         = 2.675; %mM
    Trehalose   = 0.05; %mM
    T6P         = 0.02; %mM
    UDP_glucose = 0.7; %mM

else
	% floating variable values
    Glc         = x_values(1);
    G1P         = x_values(2);
    G6P         = x_values(3);
    Trehalose   = x_values(4);
    T6P         = x_values(5);
    UDP_glucose = x_values(6);

end

%--------------------------------------------------------
% assignment rules

%--------------------------------------------------------
% algebraic rules

%--------------------------------------------------------
% calculate concentration values

if (nargin == 0)

	% initial values
	xdot(1) = 0.09765;
	xdot(2) = 0.1;
	xdot(3) = 2.675;
	xdot(4) = 0.05;
	xdot(5) = 0.02;
	xdot(6) = 0.7;

else

	% rate equations
    rateEquations_smallbone2011;

    % ODE system
    xdot(1) = + v_GLT - v_HXK + 2 .* v_NTH1;        %Glc
    xdot(2) = + v_PGM1 - v_UDP;                     %G1P
    xdot(3) = + v_HXK - v_PGM1 - v_TPS1 - v_PGI;    %G6P
    xdot(4) = + v_TPS2 - v_NTH1;                    %Trehalose
    xdot(5) = + v_TPS1 - v_TPS2;                    %T6P
    xdot(6) = + v_UDP - v_TPS1;                     %UDP_glucose
    
end
