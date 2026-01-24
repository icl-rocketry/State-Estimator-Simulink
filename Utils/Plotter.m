classdef Plotter < matlab.System
    % Live 3D Sensor + NED + Cuboid Plotter for Simulink
    % Inputs:
    % accel, h_accel, mag, gyro, quat, north, east, down

    properties(Access = private)
        fig
        ax
        q_acc
        q_hacc
        q_mag
        q_gyro
        q_north
        q_east
        q_down
        cuboid
        baseVerts
    end

    methods(Access = protected)

        function setupImpl(obj)
            coder.extrinsic('figure','clf','axes','quiver3','hold','grid','view','xlabel','ylabel','zlabel','legend','patch');

            obj.fig = figure(1001); clf(obj.fig);
            obj.fig.Name = "3D Sensor + NED Viewer";
            obj.fig.Color = [1 1 1];

            obj.ax = axes(obj.fig);
            hold(obj.ax,'on');
            grid(obj.ax,'on');

            xlabel(obj.ax,'X');
            ylabel(obj.ax,'Y');
            zlabel(obj.ax,'Z');

            view(obj.ax,3);
            axis(obj.ax,[-2 2 -2 2 -2 2]);

            % Sensor vectors
            obj.q_acc   = quiver3(obj.ax,0,0,0,0,0,0,'g','LineWidth',2);
            obj.q_hacc  = quiver3(obj.ax,0,0,0,0,0,0,'c','LineWidth',2);
            obj.q_mag   = quiver3(obj.ax,0,0,0,0,0,0,'b','LineWidth',2);
            obj.q_gyro  = quiver3(obj.ax,0,0,0,0,0,0,'r','LineWidth',2);

            % NED axes (external inputs)
            obj.q_north = quiver3(obj.ax,0,0,0,0,0,0,'y','LineWidth',3);
            obj.q_east  = quiver3(obj.ax,0,0,0,0,0,0,'m','LineWidth',3);
            obj.q_down  = quiver3(obj.ax,0,0,0,0,0,0,'k','LineWidth',3);

            % Cuboid
            L = 1; W = 0.3; H = 0.1;

            [X,Y,Z] = ndgrid([-L/2 L/2],[-W/2 W/2],[-H/2 H/2]);
            verts = [X(:) Y(:) Z(:)];
            obj.baseVerts = verts; % store original shape

            faces = [1 2 4 3;
                     5 6 8 7;
                     1 2 6 5;
                     3 4 8 7;
                     1 3 7 5;
                     2 4 8 6];

            obj.cuboid = patch(obj.ax,'Faces',faces,'Vertices',verts,...
                'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.4);

            legend(obj.ax,{'Accel','High Accel','Mag','Gyro','North','East','Down','Body'});
        end

        function stepImpl(obj,accel,h_accel,mag,gyro,quat,north,east,down)
            coder.extrinsic('set','drawnow','disp');

            % Update vectors
            set(obj.q_acc, 'UData',accel(1),'VData',accel(2),'WData',accel(3));
            set(obj.q_hacc,'UData',h_accel(1),'VData',h_accel(2),'WData',h_accel(3));
            set(obj.q_mag, 'UData',mag(1),'VData',mag(2),'WData',mag(3));
            set(obj.q_gyro,'UData',gyro(1),'VData',gyro(2),'WData',gyro(3));

            % ✅ Use externally supplied NED vectors
            set(obj.q_north,'UData',north(1),'VData',north(2),'WData',north(3));
            set(obj.q_east, 'UData',east(1),'VData',east(2),'WData',east(3));
            set(obj.q_down, 'UData',down(1),'VData',down(2),'WData',down(3));

            % ✅ Rotate cuboid using quaternion
            if numel(quat)==4
                R = quat2rotm_local(quat);
                verts = (R * obj.baseVerts')';
                set(obj.cuboid,'Vertices',verts);
            end

            drawnow limitrate nocallbacks
        end

        function releaseImpl(obj)
            coder.extrinsic('close');
            try close(obj.fig); end
        end

        function icon = getIconImpl(~)
            icon = "3D + NED";
        end
    end

    methods(Static, Access = protected)

        function simMode = getSimulateUsingImpl
            simMode = "Interpreted execution";
        end

        function flag = isSupportedContext(context)
            flag = context.isSimulating;
        end
    end
end

function R = quat2rotm_local(q)
q = q(:)'/norm(q);
w=q(1); x=q(2); y=q(3); z=q(4);

R=[1-2*(y^2+z^2) 2*(x*y-w*z) 2*(x*z+w*y);
   2*(x*y+w*z) 1-2*(x^2+z^2) 2*(y*z-w*x);
   2*(x*z-w*y) 2*(y*z+w*x) 1-2*(x^2+y^2)];
end
