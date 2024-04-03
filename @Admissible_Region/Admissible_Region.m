classdef Admissible_Region < handle

    methods (Access = public)

        function obj = Admissible_Region(center, size)  
            obj.compute_region(center(1:2), center(3), size);
        end
        
        function compute_region(obj, center, angle, size)
            obj.center_ = center(1:2);
            obj.angle_ = angle;
            obj.size_ = size;
            obj.size_(3) = 0;
            obj.center_(3) = 0;
            
            obj.R = zeros(3);
            obj.R(1, 1) = cos(obj.angle_);
            obj.R(1, 2) = -sin(obj.angle_);
            obj.R(2, 1) = sin(obj.angle_);
            obj.R(2, 2) = cos(obj.angle_);
            obj.R(3, 3) = 1;
            
            obj.upper_left_corner = obj.center_' + obj.R * [-obj.size_(1) / 2; obj.size_(2) / 2; 0];
            obj.upper_right_corner = obj.center_' + obj.R * [obj.size_(1) / 2; obj.size_(2) / 2; 0];
            obj.lower_left_corner = obj.center_' + obj.R * [-obj.size_(1) / 2; -obj.size_(2) / 2; 0];
            obj.lower_right_corner = obj.center_' + obj.R * [obj.size_(1) / 2; -obj.size_(2) / 2; 0];
            
            polygone_corners = [obj.upper_left_corner, obj.upper_right_corner, obj.lower_right_corner,obj.lower_left_corner];
            obj.polygone_Normals = zeros(2, 4);
            obj.polygone_edges_center = zeros(2,4);
            obj.polygone_vertices = zeros(2,4);
            obj.offset = zeros(4,1);
            
            for c = 1 : 4
               iter = mod(c+1,4);
               if iter == 0 
                   iter = 4;
               end
               point_1 = polygone_corners(:,c);
               point_2 = polygone_corners(:,iter);
               vertice = (point_2 - point_1) / norm(point_2 - point_1);
               normal = cross([0, 0, 1], vertice);
               normal = normal / norm(normal);
               obj.polygone_Normals(1,c) = normal(1);
               obj.polygone_Normals(2,c) = normal(2);
               obj.polygone_vertices(1,c) = vertice(1);
               obj.polygone_vertices(2,c) = vertice(2);
               
               R_Vertices_0 = zeros(2);
               R_Vertices_0(1, 1) = obj.polygone_Normals(1, c);
               R_Vertices_0(2, 1) = obj.polygone_vertices(1, c);
               R_Vertices_0(2, 1) = obj.polygone_Normals(2, c);
               R_Vertices_0(2, 2) = obj.polygone_vertices(2, c);
               
               temp = R_Vertices_0' * [point_1(1) ; point_1(2)];
               obj.offset(c) = temp(1);
            end
        end
        function corners = Get_corners(obi)
            corners = [obj.upper_left_corner'; obj.upper_right_corner'; obj.lower_right_corner'; obj.lower_left_corner'];
        end
            
    end
    properties (Access = public)
    polygone_Normals;
    offset;
    end

    properties (Access = private)
       center_;
       size_;
       angle_;
       R;
       upper_left_corner;
       upper_right_corner;
       lower_right_corner;
       lower_left_corner;
       polygone_vertices;
       polygone_edges_center;
    end
    
end
