
%% Returns the coefficients of the acoustic paths of an active noise control system.
% In this function the acoustic paths are modelated as IIR filters. The
% coefficients was obtained from the disk of book 'Active Noise Control:
% Algorithms and DSP Implementations', Sen M. Kuo and Dennis R. Morgan.
%
% autor: Marcos Matsuo
% date: 04/04/2011
% -------------------------------------------------------------------------
% Inputs
% name: (string) acoustic path's name. Options: 'primary', 'secondary' and
%       'feedback'.
% -------------------------------------------------------------------------
% Outputs
% num: (vector) numerator of acoustic path
% den: (vector) denominator of acoustic path
% -------------------------------------------------------------------------

function [num, den] = acoustic_path_iir(name)

if strcmpi('primary', name)
    num = [ -0.237467, 0.806502, -1.882598, 3.601750, -5.938800,...
        8.734397, -11.601771, 13.989734, -15.218549, 14.661575,...
        -11.943175, 7.090693, -0.801205, -5.311028, 11.085567,...
        -15.154731, 16.927727, -17.206356, 15.141826, -12.420147,...
        8.971280, -5.755672, 3.278546, -1.155305, 0.303144 ];
    
    den = [ 1.000000, -2.616971, 5.868128, -10.254027, 16.472815,...
        -23.574261, 31.578068, -39.180908, 46.125938, -51.048965,...
        53.950325, -54.055084, 51.845238, -47.205044, 41.234455,...
        -34.143761, 27.083752, -20.176018, 14.395070, -9.473912,...
        5.902801, -3.258903, 1.655440, -0.656055, 0.208196 ];
    
elseif strcmpi('secondary', name)
    num = [ -0.073844, 0.280609, -0.778907, -0.321869, 2.563927,...
        -5.876712, 11.504314, -18.114344, 25.967108, -33.871861,...
        41.345085, -47.070049, 50.785450, -51.863720, 50.043129,...
        -45.990845, 39.786583, -33.078777, 25.692730, -18.861567,...
        12.781915, -8.174841, 4.567914, -2.076808, 0.728670 ];
    
    den = [ 1.000000, -2.126286, 4.237403, -6.964616, 10.683328,...
        -14.553812, 18.877407, -22.745882, 26.226549, -28.574720,...
        29.879738, -29.674311, 28.296484, -25.745836, 22.470392,...
        -18.625055, 14.694556, -10.964108, 7.748295, -4.980707,...
        2.969857, -1.507602, 0.719325, -0.266916, 0.055631 ];
    
elseif strcmpi('feedback', name)
    num = [ 0.028552, -0.155516, 0.454342, -0.981401, 1.857463,...
        -3.132221, 4.800774, -6.799313, 8.965320, -11.037522,...
        12.739935, -13.814743, 14.133910, -14.179461, 12.947546,...
        -11.282589, 9.336633, -6.914559, 5.270569, -3.164825,...
        1.657921, -0.479267, 0.027254, -0.066006, -0.092367 ];
    
    den = [ 1.000000, -2.335795, 4.744195, -8.090243, 12.767008,...
        -18.256697, 24.395975, -30.629076, 36.431721, -41.124641,...
        44.368660, -45.695057, 45.019215, -42.287392, 38.038147,...
        -32.600128, 26.678467, -20.568045, 15.046448, -10.192312,... 
        6.476933, -3.708652, 1.886758, -0.784532, 0.233771 ];
end

