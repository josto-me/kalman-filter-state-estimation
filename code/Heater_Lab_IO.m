classdef Heater_Lab_IO < matlab.System
    % Heater_Lab_IO  MATLAB System block: reads T1/T2 and writes heater Q1/Q2 of the TCLab.
    %
    % Johannes Stockhammer.
    % The Arduino connection and the sensor/heater functions in setupImpl
    % are adapted from tclab.m of APMonitor/arduino
    % (https://github.com/APMonitor/arduino, 0_Test_Device/MATLAB/tclab.m),
    % licensed under the Apache License 2.0. Changes: wrapped in a
    % matlab.System class, temperatures returned in K. See THIRD_PARTY.md.

    % Public, tunable properties
    properties

    end

    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)
        T1C
        T2C
        h1
        h2

    end

    methods(Access = protected)
        %% Function inputs
        % u_tclab = [Q1,Q2]
        % Q1 control input heater 1 %
        % Q2 control input heater 2 %
        % y_k = [T1,T2]
        % T1 sensor temperature heater 1 K
        % T2 sensor temperature heater 1 K
        %% Function outputs
        % y = [T1,T2]
        % T1 temperature heater 1 K
        % T2 temperature heater 2 K
        %% Set up inputs and outputs
        function num = getNumInputsImpl(~)
            num = 1;                                    % number of inputs
        end % end function
        function num = getNumOutputsImpl(~)
            num = 1;                                    % number of outputs
        end % end function
        function dt1 = getOutputDataTypeImpl(~)
        	dt1 = 'double';                             % output data y
        end % end function
        function dt1 = getInputDataTypeImpl(~)
        	dt1 = 'double';                             % input data type u
        end % end function
        function sz1 = getOutputSizeImpl(~)
        	sz1 = [1,2];                                % vector dimension output y
        end % end function
        function sz1 = getInputSizeImpl(~)
        	sz1 = [1,2];                                % vector dimension input u
        end % end function
        function cp1 = isInputComplexImpl(~)
        	cp1 = false;                                % input u not complex
        end % end function
        function cp1 = isOutputComplexImpl(~)
        	cp1 = false;                                % output y not complex
        end % end function
        function fz1 = isInputFixedSizeImpl(~)
        	fz1 = true;                                 % input u fixed size
        end % end function
        function fz1 = isOutputFixedSizeImpl(~)
        	fz1 = true;                                 % output y fixed size
        end % end function
        function setupImpl(obj,~)
            % Perform one-time calculations, such as computing constants
          % connect to Arduino
          try
              a = arduino;
              disp(a)
          catch
              warning('Unable to connect, user input required')
              disp('For Windows:')
              disp('  Open device manager, select "Ports (COM & LPT)"')
              disp('  Look for COM port of Arduino such as COM4')
              disp('For MacOS:')
              disp('  Open terminal and type: ls /dev/*.')
              disp('  Search for /dev/tty.usbmodem* or /dev/tty.usbserial*. The port number is *.')
              disp('For Linux')
              disp('  Open terminal and type: ls /dev/tty*')
              disp('  Search for /dev/ttyUSB* or /dev/ttyACM*. The port number is *.')
              disp('')
              com_port = input('Specify port (e.g. COM4 for Windows or /dev/ttyUSB0 for Linux): ','s');
              a = arduino(com_port,'Uno');
              disp(a)
          end
          
          % voltage read functions
          v1 = @() readVoltage(a, 'A0');
          v2 = @() readVoltage(a, 'A2');
          
          % temperature calculations as a function of voltage for TMP36
          TC = @(V) (V - 0.5)*100.0;          % Celsius
          TK = @(V) TC(V) + 273.15;           % Kelvin
          TF = @(V) TK(V) * 9.0/5.0 - 459.67; % Fahrenhiet
          
          % temperature read functions
          obj.T1C = @() TC(v1());
          obj.T2C = @() TC(v2());
          
          % LED function (0 <= level <= 1)
          led = @(level) writePWMDutyCycle(a,'D9',max(0,min(1,level)));  % ON
          
          % heater output (0 <= heater <= 100)
          % limit to 0-0.9 (0-100%)
          obj.h1 = @(level) writePWMDutyCycle(a,'D3',max(0,min(100,level))*0.9/100);
          % limit to 0-0.5 (0-100%)
          obj.h2 = @(level) writePWMDutyCycle(a,'D5',max(0,min(100,level))*0.5/100);
        end

        function y_sensor = stepImpl(obj,u_control)
            y_sensor = [obj.T1C()+273.15, obj.T2C()+273.15];
            obj.h1(u_control(1));
            obj.h2(u_control(2));
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            %y = u;
        end

        function resetImpl(obj)
            obj.h1(0);
            obj.h2(0);
            % Initialize / reset discrete-state properties
        end
    end
end
