function plate_bending_analysis()
    % Plate Bending Analysis - Cleaned up version without reordering
    clear; clc; close all;

    x = sym("x");
    y = sym("y");
    
    fprintf('Plate Bending Analysis\n');
    fprintf('=====================\n');
    
    % User inputs with validation
    plate_length = input('Enter the length of the plate: ');
    plate_breadth = input('Enter the breadth of the plate: ');
    divisions = input('Enter the number of divisions on length and breadth: ');
    
    % Validate divisions is a positive integer
    if divisions <= 0 || floor(divisions) ~= divisions
        error('Number of divisions must be a positive integer.');
    end
    
    % Calculate dimensions
    half_length = double(plate_length / (divisions * 2));
    half_breadth = double(plate_breadth / (divisions * 2));
    
    % Calculate total DOFs
    total_nodes = (divisions + 1)^2;
    total_dofs = total_nodes * 3;
    
    % Ensure integer for matrix creation
    total_dofs = double(total_dofs);
    
    fprintf('Mesh info: %d nodes, %d degrees of freedom\n', total_nodes, total_dofs);
    
    % Polynomial basis for plate bending (12 terms)
    poly_basis = [1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3, x^3*y, x*y^3];
    
    % Node matrix: [w; theta_y; theta_x]
    node_matrix = sym(zeros(3, 12));
    node_matrix(1,:) = poly_basis;                    % w displacement
    node_matrix(2,:) = diff(poly_basis, y);           % theta_y = ∂w/∂y
    node_matrix(3,:) = -diff(poly_basis, x);          % theta_x = -∂w/∂x
    
    % =========================================================================
    % ELEMENT MATRIX ASSEMBLY
    % =========================================================================
    
    % Define node coordinates in natural counter-clockwise order
    node_coords = [
        -half_length,  half_breadth;   % Node 1: bottom-left
         half_length,  half_breadth;   % Node 2: bottom-right  
         half_length, -half_breadth;   % Node 3: top-right
        -half_length, -half_breadth    % Node 4: top-left
    ];
    
    % Build element matrix by evaluating at each node
    element_matrix = sym(zeros(12, 12));
    for node_idx = 1:4
        x_val = node_coords(node_idx, 1);
        y_val = node_coords(node_idx, 2);
        temp = subs(node_matrix, {x, y}, {x_val, y_val});
        element_matrix((node_idx-1)*3+1:node_idx*3, :) = temp;
    end
    
    % Curvature matrix (strain-displacement relationship in polynomial space)
    curvature_matrix = sym(zeros(3, 12));
    curvature_matrix(1,:) = -diff(diff(poly_basis, x), x);        % ∂²w/∂x²
    curvature_matrix(2,:) = -diff(diff(poly_basis, y), y);        % ∂²w/∂y²
    curvature_matrix(3,:) = -2 * diff(diff(poly_basis, x), y);    % 2∂²w/∂x∂y
    
    % B-matrix calculation (strain-displacement matrix in physical coordinates)
    b_matrix = curvature_matrix / element_matrix;
    
    % Material properties
    poisson_ratio = input('Enter Poisson ratio of the plate: ');
    plate_thickness = input('Enter the thickness of the plate: ');
    youngs_modulus = input('Enter Young''s modulus of the plate: ');
    
    % Load function
    fprintf('Enter load function (e.g., 100 for uniform load): ');
    load_func_str = input('', 's');
    load_function = str2sym(load_func_str);
    
    % Boundary condition selection
    fprintf('\nBoundary Conditions:\n');
    fprintf('1: Fixed plate from all sides (clamped)\n');
    fprintf('2: Simply supported from all sides\n');
    fprintf('3: Simply supported from 2 sides, fixed from other two\n');
    fprintf('4: Clamped at opposite edges\n');
    boundary_choice = input('Enter your choice: ');
    
    % Flexural rigidity and constitutive matrix
    D = youngs_modulus * plate_thickness^3 / (12 * (1 - poisson_ratio^2));
    constitutive_matrix = D * [1, poisson_ratio, 0;
                              poisson_ratio, 1, 0;
                              0, 0, (1 - poisson_ratio)/2];
    
    % Element stiffness matrix
    fprintf('Computing element stiffness matrix...\n');
    integrand = b_matrix.' * constitutive_matrix * b_matrix;
    element_stiffness = int(int(integrand, x, -half_length, half_length), ...
                          y, -half_breadth, half_breadth);
    
    % Element load vector
    fprintf('Computing element load vector...\n');
    shape_functions = poly_basis / element_matrix;
    element_load = int(int(shape_functions.' * load_function, ...
                          x, -half_length, half_length), ...
                          y, -half_breadth, half_breadth);
    
    % =========================================================================
    % GLOBAL MATRIX ASSEMBLY 
    % =========================================================================
    
    fprintf('Assembling global matrices...\n');
    global_stiffness = zeros(total_dofs, total_dofs);
    global_load = zeros(total_dofs, 1);
    
    for elem = 1:divisions^2
        % Calculate row and column indices for this element
        row = floor((elem-1)/divisions);
        col = mod(elem-1, divisions);
        
        % Node connectivity in natural order (matches element_matrix ordering)
        node1 = row * (divisions + 1) + col + 1;           % bottom-left
        node2 = node1 + 1;                                 % bottom-right
        node3 = node1 + divisions + 2;                     % top-right
        node4 = node1 + divisions + 1;                     % top-left
        
        % DOF indices for each node (3 DOFs per node)
        dofs1 = (node1-1)*3 + (1:3);
        dofs2 = (node2-1)*3 + (1:3);
        dofs3 = (node3-1)*3 + (1:3);
        dofs4 = (node4-1)*3 + (1:3);
        
        all_dofs = [dofs1, dofs2, dofs3, dofs4];
        
        % Assemble stiffness matrix
        for i = 1:12
            for j = 1:12
                global_stiffness(all_dofs(i), all_dofs(j)) = ...
                    global_stiffness(all_dofs(i), all_dofs(j)) + double(element_stiffness(i,j));
            end
        end
        
        % Assemble load vector
        for i = 1:12
            global_load(all_dofs(i)) = global_load(all_dofs(i)) + double(element_load(i));
        end
    end
    
    % =========================================================================
    % BOUNDARY CONDITIONS
    % =========================================================================
    
    penalty = 1e20;
    
    % Get boundary nodes
    bottom_nodes = 1:(divisions+1);
    top_nodes = (divisions*(divisions+1)+1):total_nodes;
    left_nodes = 1:(divisions+1):total_nodes;
    right_nodes = (divisions+1):(divisions+1):total_nodes;
    
    switch boundary_choice
        case 1 % Fixed all sides
            fixed_nodes = unique([bottom_nodes, top_nodes, left_nodes, right_nodes]);
            apply_fixed_constraints(fixed_nodes);
            
        case 2 % Simply supported all sides
            % Bottom and top: constrain w and θ_x
            for nodes = [bottom_nodes, top_nodes]
                for i = 1:length(nodes)
                    node = nodes(i);
                    global_stiffness((node-1)*3+1, (node-1)*3+1) = global_stiffness((node-1)*3+1, (node-1)*3+1) * penalty;
                    global_stiffness((node-1)*3+3, (node-1)*3+3) = global_stiffness((node-1)*3+3, (node-1)*3+3) * penalty;
                    global_load((node-1)*3+1) = 0;
                    global_load((node-1)*3+3) = 0;
                end
            end
            % Left and right: constrain w and θ_y
            for nodes = [left_nodes, right_nodes]
                for i = 1:length(nodes)
                    node = nodes(i);
                    global_stiffness((node-1)*3+1, (node-1)*3+1) = global_stiffness((node-1)*3+1, (node-1)*3+1) * penalty;
                    global_stiffness((node-1)*3+2, (node-1)*3+2) = global_stiffness((node-1)*3+2, (node-1)*3+2) * penalty;
                    global_load((node-1)*3+1) = 0;
                    global_load((node-1)*3+2) = 0;
                end
            end
            
        case 3 % Mixed: fixed in x, simply supported in y
            % Left and right edges fixed
            fixed_nodes = unique([left_nodes, right_nodes]);
            apply_fixed_constraints(fixed_nodes);
            % Top and bottom simply supported (w and θ_x constrained)
            for nodes = [bottom_nodes, top_nodes]
                for i = 1:length(nodes)
                    node = nodes(i);
                    global_stiffness((node-1)*3+1, (node-1)*3+1) = global_stiffness((node-1)*3+1, (node-1)*3+1) * penalty;
                    global_stiffness((node-1)*3+3, (node-1)*3+3) = global_stiffness((node-1)*3+3, (node-1)*3+3) * penalty;
                    global_load((node-1)*3+1) = 0;
                    global_load((node-1)*3+3) = 0;
                end
            end
            
        case 4 % Clamped at opposite edges (left and right)
            fixed_nodes = unique([left_nodes, right_nodes]);
            apply_fixed_constraints(fixed_nodes);
    end
    
    % =========================================================================
    % SOLUTION
    % =========================================================================
    
    fprintf('Solving system...\n');
    nodal_displacements = global_stiffness \ global_load;
    
    % Extract results
    displacements = nodal_displacements(1:3:end);
    slopes_y = nodal_displacements(2:3:end);
    slopes_x = -nodal_displacements(3:3:end); % Convert to conventional slope
    
    % Display center point results
    center_node = (divisions + 1) * floor(divisions/2) + floor(divisions/2) + 1;
    fprintf('\nCenter Point Results:\n');
    fprintf('Displacement: %.6e\n', displacements(center_node));
    fprintf('Slope in x-direction: %.6e\n', slopes_x(center_node));
    fprintf('Slope in y-direction: %.6e\n', slopes_y(center_node));
    
    % Plot results
    plot_results();
    
    % =========================================================================
    % NESTED FUNCTIONS
    % =========================================================================
    
    function apply_fixed_constraints(nodes)
        % Apply fixed boundary conditions (w = θ_x = θ_y = 0)
        for i = 1:length(nodes)
            node = nodes(i);
            for dof = 1:3
                dof_index = (node-1)*3 + dof;
                global_stiffness(dof_index, dof_index) = global_stiffness(dof_index, dof_index) * penalty;
                global_load(dof_index) = 0;
            end
        end
    end
    
    function plot_results()
        % Create mesh grid for plotting
        [X, Y] = meshgrid(0:plate_length/divisions:plate_length, ...
                          0:plate_breadth/divisions:plate_breadth);
        
        % Reshape results to grid format
        displacement_grid = reshape(displacements, divisions+1, divisions+1);
        slope_x_grid = reshape(slopes_x, divisions+1, divisions+1);
        slope_y_grid = reshape(slopes_y, divisions+1, divisions+1);
        
        figure('Position', [100, 100, 1200, 800]);
        
        % Plot displacement
        subplot(2,2,1);
        surf(X, Y, displacement_grid);
        title('Plate Deflection (w)');
        xlabel('X'); ylabel('Y'); zlabel('Deflection');
        colorbar;
        
        % Plot slope in x-direction
        subplot(2,2,2);
        surf(X, Y, slope_x_grid);
        title('Slope in X-direction (θ_x)');
        xlabel('X'); ylabel('Y'); zlabel('Slope');
        colorbar;
        
        % Plot slope in y-direction
        subplot(2,2,3);
        surf(X, Y, slope_y_grid);
        title('Slope in Y-direction (θ_y)');
        xlabel('X'); ylabel('Y'); zlabel('Slope');
        colorbar;
        
        % Plot contour of displacement
        subplot(2,2,4);
        contourf(X, Y, displacement_grid, 20);
        title('Deflection Contour');
        xlabel('X'); ylabel('Y');
        colorbar;
        
        fprintf('\nResults plotted successfully.\n');
    end
end