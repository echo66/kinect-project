function [ pcloud, transforms ] = pivwitharuco( image_names, visible_aruco, aruco_cam, depth_cam, rgb_cam, Rdtrgb, Tdtrgb, Rcalib, Tcalib, mapa_RT, optional )
%-----PARAMETERS-----
%image_names - an array of structures with the names of the images. Each element of the array is a structure 
%
%  image_name[k].depth - a string with the path of the .mat file with depth data of image k. 
%  image_name[k].rgb - a string with the path of a jpeg file with rgb image k
%
%visible_aruco - an array of structures with the 3D coordinates of the corners of visible markers, its image projection and the aruco labels visible in frame k.
%
% visible_aruco[k].xymrks=[ x_nk  y_nk] - an array with nk rows and 2 columns, with the image coordinates of the 3D aruco points.
% visible_aruco[k].xyzbase=[X_nk Y_nk Z_nk] - an array with nk rows and 3 columns with the 3D coordinates of the aruco points
% visible_aruco[k].ind=[ a b  c ...] - an array with the indices of the visible arucos in image k
%
%aruco_cam - A structure with the intrinsic parameters of the static camera that detected the aruco tags.
%
% aruco_cam.K - a 3x3 matrix with the intrinsic parameters
% aruco_cam.DistCoef - a 1x5 array with the lens distortion coeficients 
%
%depth_cam - A structure with the intrinsic parameters of the depth camera
%
% depth_cam.K - a 3x3 matrix with the intrinsic parameters
% depth_cam.DistCoef - a 1x5 array with the lens distortion coeficients 
%
%rgb_cam - A structure with the intrinsic parameters of the rgb camera
%
% rgb_cam.K - a 3x3 matrix with the intrinsic parameters
% rgb_cam.DistCoef - a 1x5 array with the lens distortion coeficients 
%
%Rdtrgb - a 3x3 rotation matrix 
%
% Tdrgb - a 3x1 vector
%
% Rdtrgb and Tdrgb allow  transforming 3D coordinates represented in the depth camera  reference frame to the RGB camera frame
%
% Rcalib, Tcalib - The rotation matrix and translation vector that transforms 3D coordinates in the depth camera reference frame to the "box" reference frame (center of aruco #2).
%
% mapa_RT{ i } - a 12 x 1 cell array where each element contains a structure with the following fields:
%
% mapa_RT{ i }.R  - 3x3 rotation matrix
% mapa_RT{ i }.T  - 3x1 translation vector
%
% R and T transform the 3D coordinates expressed in the reference frame of aruco number  i+1 to the box reference frame (aruco number 2)
%
% optional - an optional argument (see nargin in matlab) that you may use to show some extra feature of your program (for example if you use voxels, it may be the size of the voxel). You can pass anything using this variable, but the program must also run without this argument (by default). In other words, you must test if this variable was introduced or not and if not you should return the whole point cloud.
%
%-----RETURNS-----
%
%pcloud ->  [Xi Yi Zi Ri Gi Bi] - One  N x 6 matrix with the 3D points and RGB data of each point,  represented in the world reference frame (Aruco camera).
%
%transforms -> an array of structures with the same size of image_name where each element contains the transformation between the depth camera reference frame and the world reference frame for image k.
%
% transforms[k].R - The rotation matrix between depth camera and world for image k
% transforms[k].T - The translation vector between depth camera and world for image k

    do_voxelization = true;
    voxel_size = 0.025;
    voxel_min_count = 4;
        
    if nargin > 11

        if isfield(optional,'voxel_size'),
            voxel_size = optional.voxel_size;
        end
        
        if isfield(optional,'voxel_min_count'),
            voxel_min_count = optional.voxel_min_count;
        end
        
        if isfield(optional,'do_voxelization'),
            do_voxelization = optional.do_voxelization;
        end
    
    end

    nimages = size(image_names,2);
    transforms(1:nimages) = struct('R', [], 'T', []);
    pcloud = zeros(size(image_names,2),6);
    
    for imageIdx=1:nimages
      
        K   = aruco_cam.K;
        ids = visible_aruco(imageIdx).ind;
        uv  = visible_aruco(imageIdx).xymrks;
        xyz = visible_aruco(imageIdx).xyzbase;
        
        % Check if there's no empty cell
        if (size(ids) == 0)
            continue;
        end
        
        % Additional improvement - throw away visible arucos that are too skew.
        [ids, uv, xyz] = check_arucos( ids, uv, xyz );
        if (size(ids) == 0)
            disp(['For frame number ' imageIdx ' there is no visible arucos.']);  
            disp('-----------------------------------');
            continue;
        end
        
        % Check visible faces and decide about algorithm type
        ret = choose_algorithm(ids');
        facesIds{imageIdx} = ret.faces;
        
        % Transform to visible face
        j = ids(1);
        xyz_ = (mapa_RT{j-1}.R)' * xyz' - (mapa_RT{j-1}.R)' * repmat(mapa_RT{j-1}.T, 1, size(xyz,1));
        xyz = xyz_';
        
        % Obtain [R,T] transformation using appropriate algorithm
        if strcmp(ret.alg,'non-coplanar')==1 
            [R, T] = non_coplanar_estimation( xyz, uv, K);
        else
            [R, T] = coplanar_estimation( xyz, uv, K);
        end
        
        %variable depth in the workspace has the depth image;
        load(image_names(imageIdx).depth); 
        imrgb=imread(image_names(imageIdx).rgb);
        figure(1);
        subplot(131),imagesc(repmat( double(depth_array) / max( double(depth_array(:)) ), [1 1 3] ));
        subplot(132),imagesc(imrgb);
        
        % create cloud of points from depth image (26th September lecture)
        xyzd = get_xyz_asus(depth_array(:),[480 640], depth_array(:)~=0, depth_cam.K, 1, 0); %Compute 3D coordinates
        % connect rgb image to cloud of points
        imrgbd = get_rgbd(xyzd,imrgb,Rdtrgb,Tdtrgb,rgb_cam.K);
        subplot(133);imagesc(imrgbd);
        drawnow;
        
        % Build the return values
        transforms(imageIdx).R = R * (mapa_RT{j-1}.R)' * Rcalib;
        transforms(imageIdx).T = T + R * (mapa_RT{j-1}.R)' * ( Tcalib - mapa_RT{j-1}.T );
        
        raux=double(imrgbd(:,1)); % CHANGED rgbd TO imrgbd
        gaux=double(imrgbd(:,2)); % CHANGED rgbd TO imrgbd 
        baux=double(imrgbd(:,3)); % CHANGED rgbd TO imrgbd 
        
        % Transform cloud of points
        xyzcamera = transforms(imageIdx).R * xyzd' + repmat( transforms(imageIdx).T, 1, size(xyzd,1) );
        
        pcloud = [pcloud; [xyzcamera' raux(:) gaux(:) baux(:)]]; 
    end
    
    if do_voxelization,%voxelize ...option has the voxel size
        xyzmin=min(pcloud(:,1:3));
        xyzmax=max(pcloud(:,1:3));%the total size of the box
        pcloud(:,1:3)=fix(pcloud(:,1:3)/voxel_size);
        [~,ia,ic]=unique(pcloud(:,1:3),'rows');%get unique coordinates
        disp(['Number of initial voxels: ' num2str(size(ia,1))]);
        counts_ = zeros(size(ic'));
%        colors_ = zeros(size(ic,1),3);
        disp('Counting how many pixels are in each voxel');
        for i=ic',
            counts_(i) = counts_(i) + 1;
%            colors_(i,:) = colors_(i,:) + pcloud(i,4:6);
        end
%        disp('Averaging the color using each voxel pixels.');
        A = ismember(ic',find(counts_ >= voxel_min_count));
        A = find(A==1);
%        pcloud(A,4) = colors_(A,1)' ./ counts_(A);
%        pcloud(A,5) = colors_(A,2)' ./ counts_(A);
%        pcloud(A,6) = colors_(A,3)' ./ counts_(A);
%        pcloud(A,4:6) = mean(pcloud(pix_idxs{i},4:6));
        pcloud = pcloud(A,:);

        [paux,ia,~]=unique(pcloud(:,1:3),'rows');%get unique coordinates
        pcloud=[paux pcloud(ia,4:6)];%get the rgb of those unique coordinates - 
        
    end
    
    for imageIdx=1:nimages
        if size(transforms(imageIdx).R,1)>0,
            pcloud = [pcloud; [drawcam(transforms(imageIdx).R,transforms(imageIdx).T) zeros(6,1) 255*ones(6,1) zeros(6,1)]];
        end
    end 

end


function [new_ids, new_uv, new_xyz]  = check_arucos( ids, uv, xyz )
% Additional improvement to check which arucos are properly vivible 
% (they are not skew too much)
    new_ids = [];
    new_uv  = [];
    new_xyz = [];
    
    for k=1:length(ids)
        BA = [ uv(k*4-3,1) - uv(k*4-2,1) ; uv(k*4-3,2) - uv(k*4-2,2) ];
        BC = [ uv(k*4-1,1) - uv(k*4-2,1) ; uv(k*4-1,2) - uv(k*4-2,2) ];
        angle = atan2( BA(2)*BC(1) - BA(1)*BC(2) ,  BA(1)*BC(1) + BA(2)*BC(2) ) * 180 / pi;
        if( abs(angle) > 40 && abs(angle) < 140 )
            new_ids = [new_ids; ids(k)];
            new_uv  = [new_uv; uv(k*4-3:k*4,:)];
            new_xyz = [new_xyz; xyz(k*4-3:k*4,:)];
        end
    end
    
end

function ret = choose_algorithm( ids )
% Small function to decide which solution (coplanar or non-coplanar)
% should be applied to estimate. 
% Parameters:
%  ids: [1 x N] vector, 2 <= N <= 12
% Returns:
%  { 
%    alg:"coplanar"/"non-coplanar", 
%    faces: 1xN Array, with the visiblefaces IDs
%  }
    faces = zeros(1,4);
    
    sideIndex = [-1];
    sideIndex = [sideIndex [1 1 1]];
    sideIndex = [sideIndex [2 2 2]];
    sideIndex = [sideIndex [3 3 3]];
    sideIndex = [sideIndex [4 4 4]];
    
    for i=1:size(ids,2) 
        faces(sideIndex(ids(i))) = faces(sideIndex(ids(i))) + 1;
    end
    
    faces = find(faces>0); 
    
    if size(faces,2) > 1 
        alg = 'non-coplanar';
    else 
        alg = 'coplanar';
    end   
    
    ret = {};
    ret.alg = alg;
    ret.faces = faces;
end

function [ R, T ] = coplanar_estimation( xyz, uv, K, transfms, index)
% Function calculating [R,T] matrix in coplanar mode, using homography.
% Parameters:
%  xyz - 3D enviroment points measured in the box frame
%  uv  - 2D image points representing the projection of xyz in the image plane 
%  K   - World frame camera intrinsics parameters
% Returns:
% [R,T] - Rotation and translation 

    npairs = size(uv,1);
	
    h_uv = [uv ones(npairs,1)];
    h_xy = [xyz(:,1) xyz(:,3) ones(npairs,1)];
    
    M = zeros(npairs*2,9);
    for i = 1:npairs
        % [ x y 1   0 0 0   -u*x -u*y -u*1]
		M(2*i-1,:) = [ h_xy(i,:), zeros(1,3), -h_uv(i,1)*h_xy(i,:) ];
        % [ 0 0 0   x y 1   -v*x -v*y -v*1]
		M(2*i,:)   = [ zeros(1,3), h_xy(i,:), -h_uv(i,2)*h_xy(i,:) ];
    end;
    
    % SVD M'*M matrix -> vectors from V is eigen vectors of M'M
    [~,~,V] = svd(M'*M); 

    % the last one for the smallest eigen value
    H_matrix = reshape(V(:,9),3,3)'; 
    
    % normalization 
    H_matrix = inv(K) * H_matrix;
    if (H_matrix(3,3)<0)
        H_matrix = - H_matrix;
    end
    
    h1 = H_matrix(:,1);
    h2 = H_matrix(:,2);
    h3 = H_matrix(:,3);
    
    alpha_mag = norm([h1 h2],'fro') / sqrt(2);
    rc1 = h1 / alpha_mag;
    h2_ = h2 / alpha_mag;
    aux2 = h2_ - dot(h2_, rc1) * rc1;
    rc2 = aux2 / norm(aux2);
    
    %cross product positive sign: c1 -> c2 -> c3
    rc3 = -cross(rc1,rc2); 
    
    R = [rc1 rc3 rc2];
    T = h3 / alpha_mag;
end

function [ R, T ] = non_coplanar_estimation(xyz, uv, K, transfms, index)
% Function calculating [R,T] matrix in non-coplanar mode.
% Parameters:
%  xyz - 3D enviroment points measured in the box frame
%  uv  - 2D image points representing the projection of xyz in the image plane 
%  K   - World frame camera intrinsics parameters
% Returns:
% [R,T] - Rotation and translation 

    npairs = size(uv,1);
    
	% homogeneous coordinates
    h_uv = [uv ones(npairs,1)];
    h_xyz = [xyz ones(npairs,1)];
    
	% DLT 
	M = zeros(npairs*2,12);
    for i = 1:npairs
		M(2*i-1,:) = [ h_xyz(i,:), zeros(1,4), -h_uv(i,1)*h_xyz(i,:) ]; 
		M(2*i,:)   = [ zeros(1,4), h_xyz(i,:), -h_uv(i,2)*h_xyz(i,:) ];
    end;
    [~,~,V] = svd(M'*M); % SVD M'*M matrix -> vectors from V is eigen vectors of M'M
    
    P_matrix = reshape(V(:,12),4,3)'; % the last one for the smallest eigen value

    % normalization 
    P_matrix = inv(K) * P_matrix;      
    % alpha*[R,t] = [Q,q] = P_matrix: conditions:   R'R=I && det(R)=1 &&
    %                                               Square of Frobenius norm of R equal 3
    Q = P_matrix(:,1:3);
    q = P_matrix(:,4); 
    alpha_abs = norm(Q,'fro')/sqrt(3); alpha_sign = sign(det(Q)); 
    alpha = alpha_abs * alpha_sign;
    R = alpha * Q;
    
    % Proscrutres Approximation 
    [u,~,v] = svd(R);
    
    R = u*v';
    T = q / alpha;
end

function rgbd = get_rgbd(xyz, rgb, R, T, K_rgb)

    Kx = K_rgb(1,1);
    Cx = K_rgb(1,3);
    Ky = K_rgb(2,2);
    Cy = K_rgb(2,3);

    xyz_rgb = R * xyz';
    xyz_rgb = [xyz_rgb(1,:) + T(1); xyz_rgb(2,:) + T(2); xyz_rgb(3,:) + T(3)];

    x = xyz_rgb(1,:);
    y = xyz_rgb(2,:);
    z = xyz_rgb(3,:);

    u = round(Kx * x./z + Cx);
    v = round(Ky * y./z + Cy);

    rgb_size = size(rgb);
    n_pixels = numel(rgb(:,:,1));

    v(v > rgb_size(1)) = 1;
    v(v < 1) = 1;
    u(u > rgb_size(2)) = 1;
    u(u < 1) = 1;

    rgb_inds = sub2ind(rgb_size, v, u);
    rgb_aux = reshape(rgb,480*640,3);
    rgbd = rgb_aux(rgb_inds,:);
    rgbd( xyz(:,3) == 0,:) = 0;

end

function xyz_ret = get_xyz_asus(im_vec, im_orig_size, good_inds, K, alpha, beta)
    % im_vec - depth image vectorized (Nx1)
    % im_orig_size - original image size (HxW) : [H, W]
    % goot_inds - indexes of the image that are valid, i.e., different from 0.

    persistent u;
    persistent v;
    persistent im_size;
    persistent xyz;
    persistent z;

    Kx = K(1,1);
    Cx = K(1,3);
    Ky = K(2,2);
    Cy = K(2,3);

    if isempty(im_size)
        %     im_size = size(im);
        im_size = im_orig_size;

        u = repmat(1:im_size(2),im_size(1),1);
        u = u(:)-Cx;
        v = repmat((1:im_size(1))',im_size(2),1);
        v=v(:)-Cy;
        xyz=zeros(length(u),3);
    end

    % tmp = im(:);
    xyz(:,3) = double(im_vec)*0.001; % Convert to meters
    xyz(good_inds,3) = alpha*xyz(good_inds,3) + beta;
    xyz(:,1) = (xyz(:,3)/Kx) .* u ;
    xyz(:,2) = (xyz(:,3)/Ky) .* v;

    %plot3(x,y,z,'.');axis equal
    xyz_ret = xyz;

end

function eixos = drawcam(R,T)
    %function drawcam(R,T,opt)
    % R - Camera Pose
    % T - Camera Position
    % opt - if no arugment is passed draw using traditional camera axis convention
    %       if opt=1, draw XYZ 
    f=.3;
    w=.2;
      eixos=[0 0 0;2*f 0 0;0 0 0;0 2*f 0;0 0 0;0 0 2*f];
      eixos(:,3)=eixos(:,3)-f;
      tri1=[0 0 0;w  -w f;w w f];
      tri2=[0 0 0;w -w f ;-w -w f];
      tri3=[0 0 0;-w -w f;-w w f];
      tri4=[0 0 0;-w w f;w w f];
    tri1=tri1*R'+repmat(T',3,1);
    tri2=tri2*R'+repmat(T',3,1);
    tri3=tri3*R'+repmat(T',3,1);
    tri4=tri4*R'+repmat(T',3,1);
    eixos=eixos*R'+repmat(T',6,1);
    if nargin==2,
        tri1=[tri1(:,1) tri1(:,3) -tri1(:,2) ];
        tri2=[tri2(:,1) tri2(:,3) -tri2(:,2) ];
        tri3=[tri3(:,1) tri3(:,3) -tri3(:,2) ];
        tri4=[tri4(:,1) tri4(:,3) -tri4(:,2) ];
        eixos=[eixos(:,1) eixos(:,3) -eixos(:,2)];
    end
end