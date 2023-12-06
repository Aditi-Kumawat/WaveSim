classdef fns_Source
    methods (Static)
        %%
        function y = double_triangle_wave(t, T, A0)
            % double_triangle_wave: A function to compute the value of a double triangular wave
            % with period T and amplitude A0 at time t.

            t = mod(t, T);  % Ensure periodicity
            if t < T/4
                y = 4*A0*t/T;
            elseif t < T/2
                y = -4*A0*t/T + 2*A0;
            elseif t < 3*T/4
                y = -4*A0*t/T + 2*A0;
            else
                y = 4*A0*t/T - 4*A0;
            end
        end

        function [Aomega,Domega]=get_Domega(A0,T,omega)
            % Integration results from 0 to T/4
            int1 = (4*A0*(-1 + exp(-1/4*1i*T.*omega).*(1 + 1i*T.*omega/4)))./(T.*omega.^2);

            % Integration results from T/4 to T/2
            int2 = (A0*exp(-1/2*1i*T.*omega).*(-4 + exp(1i*T.*omega/4).*(4 - 1i*T.*omega)))./(T.*omega.^2);

            % Integration results from T/2 to 3T/4
            int3 = (A0*exp(-3/4*1i*T.*omega).*(-1i*T.*omega + 4*exp(1i*T.*omega/4) - 4))./(T.*omega.^2);

            % Integration results from 3T/4 to T
            int4 = (A0*exp(-1i*T.*omega).*(4 + exp(1i*T.*omega/4).*(-4 + 1i*T.*omega)))./(T.*omega.^2);

            % Summing up the results
            Aomega = int1 + int2 + int3 + int4;
            Domega=Aomega./(-omega.^2);
        end
    end
end