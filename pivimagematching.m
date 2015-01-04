function [pcloud, transforms, matches] = pivimagematching( image_names, depth_cam, rgb_cam, Rdtrgb, Tdtrgb, optional)
%image_names - an array of structures with the names of the images. Each element of the array is a structure 
%  image_name[k].depth - a string with the path of the .mat file with depth data of image k. 
%  image_name[k].rgb - a string with the path of a jpeg file with rgb image k
%
%depth_cam - A structure with the intrinsic parameters of the depth camera
% depth_cam.K - a 3x3 matrix with the intrinsic parameters
% depth_cam.DistCoef - a 1x5 array with the lens distortion coeficients 
%
%rgb_cam - A structure with the intrinsic parameters of the rgb camera
% rgb_cam.K - a 3x3 matrix with the intrinsic parameters
% rgb_cam.DistCoef - a 1x5 array with the lens distortion coeficients 
%
%Rdtrgb - a 3x3 rotation matrix 
%Tdrgb - a 3x1 vector
% Rdtrgb and Tdrgb allow  transforming 3D coordinates represented in the depth camera  reference frame to the RGB camera frame
%
% RETURN VALUES
%
%pcloud -  [Xi Yi Zi Ri Gi Bi] - One  N x 6 matrix with the 3D points and RGB data of each point,  represented in the world reference frame (you choose the reference frame! For example, could be  the depth camera coordinate system of the first image).
%
%transforms - an array of structures with the same size of image_name where each element contains the transformation between the depth camera reference frame and the world reference frame (as selected above) for image k.
% the fields are the following:
% transforms[k].R - The rotation matrix between depth camera and world for image k
% transforms[k].T - The translation vector between depth camera and world for image k

    do_voxelization = true;
    noise_factor = 0;
    voxel_size = 0.025;
    ransac_iterations = 100;
    voxel_min_count = 4;
    hop_size = 1;
    descriptor_matching_threshold = 1.5;
    hop_in_both_directions = false;
    depth_loading_method = 1;
    
    if nargin > 5
        
        if isfield(optional,'depth_loading_method')
            depth_loading_method = optional.depth_loading_method;
        end
        
        if isfield(optional,'hop_in_both_directions'),
            hop_in_both_directions = optional.hop_in_both_directions;
        end
        
        if isfield(optional,'noise_factor'),
            noise_factor = optional.noise_factor;
        end

        if isfield(optional,'voxel_size'),
            voxel_size = optional.voxel_size;
        end
        
        if isfield(optional,'ransac_iterations'),
            ransac_iterations = optional.ransac_iterations;
        else
            ransac_iterations = 100;
        end
        
        if isfield(optional,'voxel_min_count'),
            voxel_min_count = optional.voxel_min_count;
        end
        
        if isfield(optional,'descriptor_matching_threshold'),
            descriptor_matching_threshold = optional.descriptor_matching_threshold;
        end
        
        if isfield(optional,'do_voxelization'),
            do_voxelization = optional.do_voxelization;
        end
        
        if isfield(optional,'hop_size'),
            hop_size = optional.hop_size;
        end
    
    end

    nimages = size(image_names,2);
    
    transforms(1:nimages) = struct('R', [], 'T', []);
    pcloud = zeros(size(image_names,2),6);
    
    matches(1:nimages) = struct('from', [], 'to', []);
    
    AdjencyMatrix = hops_list( nimages, hop_size, hop_in_both_directions);
    I = []; J = []; W = []; TT = cell(nimages, nimages);
    
    remaining_els = size(find(AdjencyMatrix==1),1);
    
    for i=1:size(AdjencyMatrix,1),
        for j=1:size(AdjencyMatrix,2),
            if AdjencyMatrix(i,j)==1,
                if ~isfield(TT{i,j},'R')
                    [R, T, P1, P2, error] = point_cloud_transform(image_names(i),image_names(j), depth_cam, Rdtrgb, Tdtrgb, rgb_cam, noise_factor, ransac_iterations, descriptor_matching_threshold, depth_loading_method);
                    
                    remaining_els = remaining_els - 2;
                    disp(remaining_els);
                    
                    if error==Inf,
                        continue;
                    end
                    
                    % X0 = R*X1+T
                    TT{i,j}.R = R; 
                    TT{i,j}.T = T; 
                    TT{i,j}.x0 = i;
                    TT{i,j}.x0_rgb = image_names(i).rgb;
                    TT{i,j}.x0_depth = image_names(i).depth;
                    TT{i,j}.x1 = j;
                    TT{i,j}.x1_rgb = image_names(j).rgb;
                    TT{i,j}.x1_depth = image_names(j).depth;
                    TT{i,j}.P1 = P1;
                    TT{i,j}.P2 = P2;
                    
                    % X1 = R'*X0 - R'*T
                    TT{j,i}.R = R';
                    TT{j,i}.T = -R'*T;
                    TT{j,i}.x0_rgb = image_names(j).rgb;
                    TT{j,i}.x0_depth = image_names(j).depth;
                    TT{j,i}.x0 = j;
                    TT{j,i}.x1 = i;
                    TT{j,i}.x1_rgb = image_names(i).rgb;
                    TT{j,i}.x1_depth = image_names(i).depth;
                    TT{j,i}.P1 = P2;
                    TT{j,i}.P2 = P1;                    
                    
                    % graph indexes and weights
                    I = [I; i; j];
                    J = [J; j; i];
                    W = [W; error; error];
                end
            end
        end
    end
    
    DG = sparse(I,J,W);
    save('all.mat','AdjencyMatrix','depth_cam','DG','I','image_names','J','nimages','Rdtrgb','rgb_cam','Tdtrgb','TT','W');
    [ST,~] = graphminspantree(DG);
    mean_path_error = zeros(nimages,1);
    for j=1:nimages, % nó rotas
        [ST,~] = graphminspantree(DG, 1);
        [dist, ~, ~] = graphshortestpath(ST, j, 'Directed', false);
        mean_path_error(j) = mean(dist);
    end
    [~,best_path_source] = max(mean_path_error);
    [~,paths,~]=graphshortestpath(ST,best_path_source,'Directed',false);
    
    transforms(best_path_source).R = eye(3);
    transforms(best_path_source).T = [0;0;0];
    transforms(best_path_source).to_rgb_file = image_names(best_path_source).rgb;
    transforms(best_path_source).to_depth_file = image_names(best_path_source).depth;
    matches(best_path_source).P1 = [];
    matches(best_path_source).P2 = [];

    im1 = imread(image_names(best_path_source).rgb);
    if depth_loading_method == 1
        depth1 = load(image_names(best_path_source).depth); depth1 = depth1.depth_array;
    else
        depth1 = imread(image_names(best_path_source).depth); 
    end
    xyzd1 = get_xyz_asus(depth1(:),[480 640], depth1(:)~=0, depth_cam.K, 1, 0); %Compute 3D coordinates
    imrgbd1 = get_rgbd(xyzd1, im1, Rdtrgb, Tdtrgb, rgb_cam.K);
    raux=double(imrgbd1(:,1)); 
    gaux=double(imrgbd1(:,2)); 
    baux=double(imrgbd1(:,3)); 
    xyzcamera = transforms(best_path_source).R * xyzd1' + repmat( transforms(best_path_source).T, 1, size(xyzd1,1) );
    pcloud = []; 
    pcloud = [pcloud; [xyzcamera' raux(:) gaux(:) baux(:)]]; 

    for i=1:size(paths,2),
        if i == best_path_source,
            continue;
        end
        %disp(i)
        %transforms(i).rgb_files_sequence = cell(size(paths,2),1);
        path_ = paths(i); path_ = path_{1};
        cRot = transforms(best_path_source).R;
        cTrans = transforms(best_path_source).T;
        for j=2:size(path_,2),
            
            transforms(i).from_rgb_file = TT{path_(j-1),path_(j)}.x1_rgb;
            transforms(i).from_depth_file = TT{path_(j-1),path_(j)}.x1_depth;
            transforms(i).to_rgb_file = TT{path_(j-1),path_(j)}.x0_rgb;
            transforms(i).to_depth_file = TT{path_(j-1),path_(j)}.x0_depth;
            
            transforms(i).T = cTrans + cRot * TT{path_(j-1),path_(j)}.T;
            transforms(i).R = cRot * TT{path_(j-1),path_(j)}.R;
            
            matches(i).P1 = TT{path_(j-1),path_(j)}.P1;
            matches(i).P2 = TT{path_(j-1),path_(j)}.P2;

            cRot = cRot * TT{path_(j-1),path_(j)}.R; 
            cTrans = transforms(i).T;
        end

        im2 = imread(image_names(i).rgb);

        if depth_loading_method==1,
            depth2 = load(image_names(i).depth); 
            depth2 = depth2.depth_array;
        else
            depth2 = imread(image_names(i).depth); 
        end

        xyzd2 = get_xyz_asus(depth2(:),[480 640], depth2(:)~=0, depth_cam.K, 1, 0); %Compute 3D coordinates

        imrgbd2 = get_rgbd(xyzd2, im2, Rdtrgb, Tdtrgb, rgb_cam.K);

        raux=double(imrgbd2(:,1)); % CHANGED rgbd TO imrgbd
        gaux=double(imrgbd2(:,2)); % CHANGED rgbd TO imrgbd 
        baux=double(imrgbd2(:,3)); % CHANGED rgbd TO imrgbd 
        
        %disp('---');

        %size(xyzd2)
        %size(transforms(i).R)
        %size(repmat( transforms(i).T, 1, size(xyzd2,1) ))
        xyzcamera = transforms(i).R * xyzd2' + repmat( transforms(i).T, 1, size(xyzd2,1) );

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

function pcloud = voxelize(pcloud_, resolution) 
    xyzmin=min(pcloud_(:,1:3));
    xyzmax=max(pcloud_(:,1:3));%the total size of the box
    pcloud(:,1:3)=fix(pcloud_(:,1:3)/resolution);
    [paux,ia,ic]=unique(pcloud(:,1:3),'rows');%get unique coordinates
    pcloud=[paux pcloud(ia,4:6)];%get the rgb of those unique coordinates - 
end

function ret = createPs( subset, imd1, f1, imd2, f2)

    P1 = []; P2 = [];
    
    for k=1:size(subset,2)
        i1 = subset(1,k);
        x1 = imd1(round(f1(2,i1)), round(f1(1,i1)), 1);
        y1 = imd1(round(f1(2,i1)), round(f1(1,i1)), 2);
        z1 = imd1(round(f1(2,i1)), round(f1(1,i1)), 3);

        i2 = subset(2,k);
        x2 = imd2(round(f2(2,i2)), round(f2(1,i2)), 1);
        y2 = imd2(round(f2(2,i2)), round(f2(1,i2)), 2);
        z2 = imd2(round(f2(2,i2)), round(f2(1,i2)), 3);

        P1 = [P1 [x1;y1;z1]];
        P2 = [P2 [x2;y2;z2]];
        
    end
    
    ret.P1 = P1; ret.P2 = P2;
    
end

function [imd1, imd2, im1, im2] = setup_image_and_depth_structs(image_files_1, image_files_2, depth_cam, Rdtrgb, Tdtrgb, rgb_cam, depth_loading_method)
    im1 = imread(image_files_1.rgb);
    im2 = imread(image_files_2.rgb);

    if depth_loading_method == 1,
        depth1 = load(image_files_1.depth); 
        depth1 = depth1.depth_array;
        depth2 = load(image_files_2.depth); 
        depth2 = depth2.depth_array;
    else
        depth1 = imread(image_files_1.depth);
        depth2 = imread(image_files_2.depth);
    end

    xyzd1 = get_xyz_asus(depth1(:),[480 640], depth1(:)~=0, depth_cam.K, 1, 0); %Compute 3D coordinates
    xyzd2 = get_xyz_asus(depth2(:),[480 640], depth2(:)~=0, depth_cam.K, 1, 0); %Compute 3D coordinates

    imrgbd1 = get_rgbd(xyzd1, im1, Rdtrgb, Tdtrgb, rgb_cam.K);
    imrgbd2 = get_rgbd(xyzd2, im2, Rdtrgb, Tdtrgb, rgb_cam.K);

    imd1 = reshape(xyzd1,480,640,3);
    imd2 = reshape(xyzd2,480,640,3);
end

function [nmatches, matches, PP1, PP2, f1, d1, f2, d2, scores, valid] = generate_point_clouds(im1, imd1, im2, imd2, noise_factor, descriptor_matching_threshold)
    %% SIFT descriptors acquisition and removal matches that have pz=0 (that is, invalid depth values).
    
    % make grayscale
    if size(im1,3) > 1, im1g = rgb2gray(im1) ; else im1g = im1 ; end
    if size(im2,3) > 1, im2g = rgb2gray(im2) ; else im2g = im2 ; end
    
    im1g = im2single(im1g) ;
    im2g = im2single(im2g) ;

    [f1,d1] = vl_sift(im1g) ;
    [f2,d2] = vl_sift(im2g) ;

    thres = 7;
    %disp('gen cloud thres');
    while thres > 1.5
        %disp('-1-');
        [matches, scores] = vl_ubcmatch(d1,d2,thres) ;

        matches_ = [];
        scores_ = [];
        for j=1:size(matches,2)
            z1 = imd1(round(f1(2,matches(1,j))), round(f1(1,matches(1,j))), 3);
            z2 = imd2(round(f2(2,matches(2,j))), round(f2(1,matches(2,j))), 3);
            if z1~=0 && z2~=0
                matches_ = [matches_ matches(:,j)];
                scores_ = [scores_ scores(j)];
            end
        end
        matches = matches_; 
        scores = scores_;
        nmatches = size(matches,2);
        %disp('-2-');
        
        %nmatches
        %thres
        
        if nmatches>=10,
            break;
        end
        
        thres = thres - 0.5;
    end
    
    if nmatches < 10 && thres == 1.5
        valid = false;
        nmatches = []; matches= []; PP1 = []; PP2 = []; f1 = []; d1 = []; f2 = [];  d2 = []; scores =[];
        disp('there is no valid matching between these two images');
        return;
    end

    ret = createPs( matches, imd1, f1, imd2, f2); 

    PP1 = ret.P1 + noise_factor * rand(3,nmatches); 
    PP2 = ret.P2 + noise_factor * rand(3,nmatches);
    
    valid = true;
end

function [R, T, P1, P2, subset, error] = ransac_iteration(matches, imd1, f1, PP1, imd2, f2, PP2, nmatches, noise_factor)

    [subset, Y] = vl_colsubset(matches,4);

    ret = createPs( subset, imd1, f1, imd2, f2);

    P1 = ret.P1 + noise_factor * rand(3,size(subset,2));
    P2 = ret.P2 + noise_factor * rand(3,size(subset,2));

    [R,T] = rigid_transform_3D(P2',P1');

    PP1_ = R * PP2 + repmat(T,1,nmatches);
    dx = PP1_(1,:) - PP1(1,:) ;
    dy = PP1_(2,:) - PP1(2,:) ;
    dz = PP1_(3,:) - PP1(3,:) ;
    error = (1/nmatches)*sum(sqrt(dx.*dx + dy.*dy + dz.*dz));
    
end

% Transform from image 2 to image 1.
function [R, T, P1, P2, error] = point_cloud_transform(image_files_1, image_files_2, depth_cam, Rdtrgb, Tdtrgb, rgb_cam, noise_factor, ransac_iterations, descriptor_matching_threshold, depth_loading_method) 

    disp(['from: ' image_files_1.rgb]);
    
    disp(['to: ' image_files_2.rgb]);
    
    [imd1, imd2, im1, im2] = setup_image_and_depth_structs(image_files_1, image_files_2, depth_cam, Rdtrgb, Tdtrgb, rgb_cam, depth_loading_method);

    [nmatches, matches, PP1, PP2, f1, d1, f2, d2, scores, valid] = generate_point_clouds(im1, imd1, im2, imd2, noise_factor, descriptor_matching_threshold);
    
    if ~valid,
        error = Inf;
        R = []; T = []; P1 = []; P2 = [];
        return;
    end

    %% RANSAC
    for j=1:ransac_iterations
        [R, T, P1, P2, subset, e] = ransac_iteration(matches, imd1, f1, PP1, imd2, f2, PP2, nmatches, noise_factor);

        error(j) = e;

        aux(j).R = R;
        aux(j).T = T;
        aux(j).P1 = P1;
        aux(j).P2 = P2;
        aux(j).subset = subset;

    end

    [min_score,min_idx] = min(error);

    T = aux(min_idx).T;
    R = aux(min_idx).R;
    P1 = aux(j).P1;
    P2 = aux(j).P2;
    error = min_score;
        
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


% expects row data
function [R,t] = rigid_transform_3D(A, B)
    if nargin ~= 2
	    error('Missing parameters');
    end

    %assert(size(A) == size(B))

    centroid_A = mean(A);
    centroid_B = mean(B);

    N = size(A,1);

    H = (A - repmat(centroid_A, N, 1))' * (B - repmat(centroid_B, N, 1));

    [U,~,V] = svd(H);

    R = V*U';

    if det(R) < 0
        %disp('Reflection detected\n');
        V(:,3) = -1 * V(:,3);
        R = V*U';
    end

    t = -R*centroid_A' + centroid_B';
end

function AdjencyMatrix = hops_list( size, hop, hop_in_both_directions )
    
    AdjencyMatrix = zeros(size,size);
    
    for i=1:size
        if hop_in_both_directions,
            lower_bound = max(1,(i-hop));
        else
            lower_bound = i+1;
        end
        for j=lower_bound:min((i+hop),size)
            if j~=i,
                AdjencyMatrix(i,j) = 1;
            end
        end
    end

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