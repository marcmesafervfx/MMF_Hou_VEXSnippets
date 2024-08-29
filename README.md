# Houdini VEX Snippets
## Initial Provisions
This repository is designated to be a place where I put some of the VEX snippets I've been using to fix, check, create, and manipulate information in different contexts. If something needs to be revisited, let me know so I can check for it and commit any of the requested modifications.

## Vector Along Curve
*Reference Code*: 72854126

> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to the curve.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Create tangent based on neighbours in a line. """;

// Get neighbours of current point and capture it's position.
int neigh[] = neighbours(0, @ptnum);
vector pos = point(0, "P", neigh[-1]);

// Get direction vector by subtracting the current position to
// the neighbour one and normalize the vector to get the proper
// length to work with.
vector tan = normalize(v@P-pos);

// Check if the current point is equal to the maximum points
// minus one (ptnum starts from 0) and negate tangent to obtain
// the opposite direction.
if(@ptnum==@numpt-1) tan*=-1;

// Set attribute.
v@tan = tan;
```

## Blur Point Positions
*Reference Code*: 26692791
> [!TIP]
> Create a for-loop and loop the attribute wrangle with feedback.

### blur_position
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a scatter node.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Blue based on nearpoints. """;

// Get maximum distance and maximum points.
float maxdist = chf("maxdist");
int maxpts = chi("maxpts");

// Get near points based on distance and max points.
int nearpts[] = nearpoints(0, v@P, maxdist, maxpts);

// Initialize pos variable with current position value.
vector pos = v@P;

// Iterate for each of the near points.
foreach(int i; nearpts){
    
    // Add value to the pos variable.
    pos += point(0, "P", i); 
}

// Divide position by the amount of near points that you iterated.
// We add an additional value because of the initial value of the
// current position.
pos/=len(nearpts)+1;

// Set attribute value.
v@P = pos;
```

## Cluster By Point Proximity
*Reference Code*: 45043176
> [!NOTE]
> The following snippet contains two variables: *primpoints* and *nearpoint*. Both output similar results, but the methodology is different.

### primpoints
> [!IMPORTANT]
> **Mode:** Primitives.
> - **Input 0:** connected to a geometry.
> - **Input 1:** connected to a scatter node.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Compute clusters avoiding promoting parameter. """;

// Get primitive points.
int pts[] = primpoints(0, @primnum);

// Initialize cluster as -1.
int cluster=-1;

// Iterate for each primitive points.
foreach(int i; pts){

    // Get position of the current point.
    vector pos = point(0, "P", i);
    
    // Get near point from second input.
    int new_cluster = nearpoint(1, pos);
    
    //If the previous cluster is bigger keep it (emulates the max method).
    if(new_cluster>cluster){       
        cluster = new_cluster;
    }
}

// Set the cluster.
i@cluster = cluster;
```
### nearpoint
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** connected to a scatter node.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Set cluster by proximity. """;

// Set cluster value.
i@cluster = nearpoint(1, v@P);
```

## NGon Detector
*Reference Code*: 50655883

### ngon_detector
> [!IMPORTANT]
> **Mode:** Primitives.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Group ngons to fix them later. """;

// Detect amount of points composing each primitive.
int pts[] = primpoints(0, i@primnum);

// Check if length of the array is bigger than 4.
if(len(pts)>4){
    
    // Set point group.
    setprimgroup(0, "ngons", i@primnum, 1);
}
```

## Normalize Distance
*Reference Code*: 89906276
> [!NOTE]
> The following snippet contains two variables: *get_distance + normalize_distance* and *normalize_distance_detail*. Both output similar results, but the methodology is different.

### get_distance
> [!IMPORTANT]
> **Mode:** Primitives.
> - **Input 0:** connected to a geometry.
> - **Input 1:** connected to a reference point.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Get the distance for each of the points. """;

// Get position from the second input.
vector pos = point(1, "P", 0);

// Get distance between current point and input 2 positon.
float dist = distance(pos, v@P);

// Set distance attribute.
f@dist = dist;
```

> [!NOTE]
> Use a promote attribute parameter to create a maximum distance value in Detail mode without removing the previous values to follow the next step.

### normalize_distance
> [!IMPORTANT]
> **Mode:** Primitives.
> - **Input 0:** connected to a geometry.
> - **Input 1:** connected to a reference point.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Normalize distance using the computed max distance. """;

// Get max distance from detail.
float max_dist = detail(0, "max_dist"); 

// Normalize distance.
float norm_dist = f@dist/max_dist;

// Set color attrivute to show the normalized distance.
v@Cd = chramp("color", norm_dist);
```

### normalize_distance_detail
> [!IMPORTANT]
> **Mode:** Detail.
> - **Input 0:** connected to a geometry.
> - **Input 1:** connected to a reference point.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Normalize distance attribute. """;

// Get amount of point from first input.
int pts = npoints(0);

// Get position from the second input.
vector ref_pos = point(1, "P", 0);

// Create handle using the pcopen.
int handle = pcopen(0, "P", ref_pos, 1e09, int(1e09));

// Get farthest distance of the point cloud.
float max_dist = pcfarthest(handle);

// Iterate for each point.
for(int pt=0; pt<pts; pt++){

    // Get current point position.
    vector curr_pos = point(0, "P", pt);
    
    // Compute current distance.
    float dist = distance(curr_pos, ref_pos);
    
    // Normalize distance and remap color.
    float norm_dist = dist/max_dist;
    vector color = chramp("ramp_color", norm_dist);
    
    // Set color attrivute to show the normalized distance.
    setpointattrib(0, "Cd", pt, color);
}
```

## Frustum Camera
*Reference Code*: 38002708

### frustum_camera
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a default box.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Create camera frustum and make available the expansions. """;

// Get camera path to create the frustum from.
string cam = chs("camera");

// Initialize expansion values (x,y).
float expand_top = chf("expand_top");
float expand_bottom = chf("expand_bottom");
float expand_right = chf("expand_right");
float expand_left = chf("expand_left");

// Initialize clipping values (z).
float near_clip = chf("near_clip");
float far_clip = chf("far_clip");

// Offset position to "convert" box position into normalized
// coordinates.
vector offset_pos = set(0.5, 0.5, -0.5);
vector pos = v@P+offset_pos;

// Apply expansions based on normalized positions.
if(pos.y==1) pos.y+=expand_top;
if(pos.y==0) pos.y-=expand_bottom;
if(pos.x==1) pos.x+=expand_right;
if(pos.x==0) pos.x-=expand_left;

// Apply clipping based on normalized positions.
if(pos.z==-1) pos.z-=far_clip;
if(pos.z==0) pos.z-=near_clip;

// Set position converting from NDC coordinates to world space.
v@P = fromNDC(cam, pos);

/* In this case we inverted the process. We created an "NDC" and
we are transforming back to world space. */
```

## Remove by threshold
*Reference Code*: 9067034
> [!TIP]
> Remember that you can use the @id or @prim_id attributes instead of the @ptnum or @primnum to be consistent, but it needs to be precomputed.

### remove_points_by_threshold
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Remove by threshold. """;

// Create a random value betweem 0 and 1 for each point.
float rand_value = rand(@ptnum);

// Check if the value is smaller than the threshold.
if(rand_value<chf("threshold")){

    // Remove point.
    removepoint(0, @ptnum);
}
```
### remove_prims_by_threshold
> [!IMPORTANT]
> **Mode:** Primitives.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Remove by threshold. """;

// Create a random value betweem 0 and 1 for each prim.
float rand_value = rand(@primnum);

// Check if the value is smaller than the threshold.
if(rand_value<chf("threshold")){

    // Remove prim.
    removeprim(0, @primnum, 1);
}
```

## Primitive Centroid
*Reference Code*: 39725183

### primitive_centroid
> [!IMPORTANT]
> **Mode:** Primitive.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Compute bbox center. """;

// Get bounding box center.
vector pos = v@P;

// Create point using the computed position.
addpoint(0, pos);

// Remove unused primitive.
removeprim(0, @primnum, 1);
```

## Compute Curveu From Line
*Reference Code*: 49138898

### curveu
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a polyline.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Compute curveu. """;

// Compute curveu based on amount of points and current ptnum.
float curveu = float(@ptnum)/float(@numpt-1); 

// Set value.
f@curveu = curveu;
```

## Angle Between Two Vectors
*Reference Code*: 89221217
> [!NOTE]
> In the example code, the v@up and the v@axis are computed already. If you need to compute the angle between two other vectors or other attributes, you can susbtitute the value of the up and axis variables.

### angle_vectors
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Compute angle between two vectors. """;

// Initialize values.
vector up = v@up;
vector axis = v@axis;

// Get angle between two vectors.
float angle = degrees(acos(dot(up, axis)));

// Compute stable axis to check for values over 180.
vector stable_axis = normalize(cross(up, cross({0,1,0}, up)));

// Compute full 360 angle. 
float full_angle = (int(sign(dot(axis, stable_axis)))==-1)? angle:360-angle;

// Set angle value.
f@angle = full_angle;
```

## Unshared Points
*Reference Code*: 82391305
### unshared_points
> [!IMPORTANT]
> **Mode:** Primitives.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Group unshared points. """;

// Get amount of points composing the current primitive.
int pts = len(primpoints(0, @primnum));

// Get initial half edge of the current primitive.
int init_hedge = primhedge(0, @primnum);

// Iterate as many times as points the primitive has.
for(int i=0; i<pts; i++){

    // Get next equivalent (opposite half edge) of the current half edge.
    int equiv = hedge_nextequiv(0, init_hedge);
    
    // Get primitive number that contains the equivalent half edge.
    int prim = hedge_prim(0, equiv);

    // If the opposite primitive is the same as current, it means
    // that it doesn't have another primitive connected.
    if(prim==@primnum){
    
        // Check for the destiny and source point for the half edge.
        int pt_dst = hedge_dstpoint(0, equiv);
        int pt_src = hedge_srcpoint(0, equiv);
        
        // Set unshared point group.
        setpointgroup(0, "unshared", pt_dst, 1);
        setpointgroup(0, "unshared", pt_src, 1);
    }
    
    // Once iteration is finished, move to the next half edge of the 
    // same primitive.
    init_hedge = hedge_next(0, init_hedge);
}
```

## Basic Transform With Matrix
*Reference Code*: 32956689
### transform
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Transform object based on input values.""";

// Initialize world space matrix.
matrix orig_matrix = ident();

// Create scale vector and transform matrix.
vector scale = chv("scale");
scale(orig_matrix, scale);

// Create rotation vector.
vector rot = chv("rotation");

// Rotate function uses radians, so we convert degrees to radians.
// Value 1 stands for X. Value 2 stands for Y. Value 4 stands for Z.
rotate(orig_matrix, radians(rot.x), 1);
rotate(orig_matrix, radians(rot.y), 2);
rotate(orig_matrix, radians(rot.z), 4);

// Create translate vector and transform matrix.
vector trans = chv("translation");
translate(orig_matrix, trans);

// Set transformations.
v@P*=orig_matrix;
```

## Extract Transform
*Reference Code*: 30376309
> [!NOTE]
> This method pretends to output a similar result as a extract transform node would do.

### get_offset_matrix
> [!IMPORTANT]
> **Mode:** Detail.
> - **Input 0:** no-connected.
> - **Input 1:** connected to a rest geometry.
> - **Input 2:** connected to a transformed geometry.
> - **Input 3:** no-connected.

``` c
"""Retrieve offset comparing original rest and moving poses. """;

// Get point neighbour point.
int nb_pts = neighbour(1, 0, 0);

// Get the position of transformed point.
vector nb_pos = point(2, 'P', 0);

// Compute zaxis based on current and neighbour point direction. 
vector zaxis = normalize(nb_pos - point(2, 'P', nb_pts));

// Use normal to create y axis.
vector yaxis = point(2, "N", 0);

// Store transformation matrix from transformed axis and position point.
matrix trans_xform = maketransform(zaxis, yaxis, nb_pos);

// Get the position of rest point.
nb_pos = point(1, "P", 0);

// Compute zaxis based on current and neighbour point direction. 
zaxis = normalize(point(1, "P", 0) - point(1, 'P', nb_pts));

// Use normal to create y axis.
yaxis = point(1, "N", 0);

// Store transformation matrix from rest axis and position point.
matrix rest_xform = maketransform(zaxis, yaxis, nb_pos);

// Compute offset matrix by inverting rest matrix.
matrix totalxform = invert(rest_xform)*trans_xform;

// Create point and set the attribute.
int pt = addpoint(0, {0,0,0});
setpointattrib(0, "xform", pt, totalxform);
```
### set_transformations
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a rest geometry.
> - **Input 1:** connected to the get_offset_matrix node.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Reset transformations. """;

// Retrieve xform attribute.
matrix xform = point(1, "xform", 0);

// Apply transformations to the position.
v@P*=xform;

// Apply transformations to the normal.
v@N*=matrix3(xform);
```

## Create Checkboard
*Reference Code*: 43837465

### checkerboard
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Create checker using ternary conditions. """;

// Set the frequency for the scene.
float freq = chf("frequency");

// Compute vertical sections.
v@Cd = (sin(v@P.x*freq)<0)?0:1;

// Add horizontal sections.
v@Cd += (sin(v@P.z*freq)<0)?0:1;

// Check if there's coincidence and multiply by 0.
v@Cd *= (v@Cd.r==2)?0:1;
```

## Basic Point Deform
*Reference Code*: 39569619
> [!NOTE]
> The simplicity of the method makes the process quite limitated. The deformed rest geometry should be quite similar as the geometry that you want to deform in terms of shape.

### point_deform
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry to deform.
> - **Input 1:** connected to a deformed rest geometry.
> - **Input 2:** connected to a deformed animated geometry.
> - **Input 3:** no-connected.

``` c
""" Deform poitn position based on rest and animated geometry. """; 

// Initialize primitive and uv intrinsic coordinates.
int prim; vector uvw;
xyzdist(1, v@P, prim, uvw);

// Retrieve the position value from the animated geometry.
vector pos = primuv(2, "P", prim, uvw);

// Set position.
v@P = pos;
```

## Advanced Point Deform
*Reference Code*: 209363

### capture
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry to deform.
> - **Input 1:** connected to a deformed rest geometry.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
"""Capture the close points and weights""";

//Capture close points between lowres and highres rest pose.
float maxdist = chf('maxdist');
int maxpts = chi('maxpts');
int npts[] = nearpoints(1, v@P, maxdist, maxpts);

// Store captured points.
i[]@npts = npts;

//Iterate for each captured point and compute the distance to generate the weights.
float weights[];

foreach(int val; npts){
    vector npos = point(1, 'P', val);
    float ndist = distance(v@P, npos);
    //Invert values to have higher values being closer to its reference point
    ndist = fit(ndist, 0, maxdist, 1, 0);
    push(weights, ndist);
}

// Store captured weights.
f[]@weights = weights;
```
### get_offset_matrix
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a deformed rest geometry.
> - **Input 1:** connected to a deformed animated geometry.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Retrieve offset comparing original rest and anim geometry. """;

//Get the basic translation offset
vector npos = point(1, 'P', @ptnum);
v@offset = npos-v@P;

//Get neighbours to create the animated axis for the matrix
int npts[] = neighbours(0, @ptnum);
vector ndir = point(1, 'P', npts[0]);
vector tan1 = npos - ndir;
ndir = point(1, 'P', npts[1]);
vector tan2 = npos - ndir;
vector up = cross(tan1, tan2);

//Create the animated matrix for each of the points
matrix3 xformnew = maketransform(normalize(tan1), normalize(up));

//Create the reference axis for the matrix
ndir = point(0, 'P', npts[0]);
tan1 = v@P - ndir;
ndir = point(0, 'P', npts[1]);
tan2 = v@P - ndir;
up = cross(tan1, tan2);

matrix3 xformold = maketransform(normalize(tan1), normalize(up));

//Create the offset using addition method = invert(reference matrix)* new matrix 
matrix3 totalxform = invert(xformold)*xformnew;

3@xform = totalxform;
```
### set_deform
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to capture node.
> - **Input 1:** connected to get_offset_matrix node.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
"""Set deformation for the new geometry""";

//Initialize values and store attributes inside statements
float weights[] = f[]@weights;
int npts[] = i[]@npts;
float sumweights = 0;
vector sumoffsets = {0,0,0};
int val = 0;

//Offset based on its weigth and capture point
foreach(int npt; npts){
    //Retrieve xform, position and offset from anim
    vector opos = point(1, 'P', npt);
    matrix3 xform = point(1, 'xform', npt);
    vector offset = point(1, 'offset', npt);
    
    //Transform to the center, apply xform transformations and bring back transforms
    vector pos = v@P;
    pos -= opos;
    pos *= xform;
    pos += opos;
    pos -= v@P;
    
    //Add basic displacement to the point
    offset += pos;
    
    //Multiply offset by weights to transform based on relative position, add all influenced offsets and add all the weights
    offset *= weights[val];
    sumweights += weights[val];
    sumoffsets += offset;
    val++;
}

//Get final offset based on the influence of the weights
vector finaloffset = sumoffsets / sumweights;
v@P += finaloffset;
```

## Noise Edge Mask
*Reference Code*: 72728404

### noise_edge_mask
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Noise mask edge based on x axis. """;

// Displace position and store it in a variable.
vector pos = v@P+noise(v@P*chf("frequency"))*chf("amplitude");

// Normalize the xaxis position.
float xaxis = fit(pos.x, getbbox_min(0).x, getbbox_max(0).x, 0, 1);

// Remap normalize x axis position values and contrast them.
xaxis = chramp("axis", xaxis);

// Set color.
v@Cd = xaxis;
```

## Dihedral Offset
*Reference Code*: 29263487
> [!NOTE]
> In the example, the normal is used to compute the offset matrix. You can eventually use your custom attribute just by changing the corresponding variable values.

### compute_offset_matrix
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a rest geometry.
> - **Input 1:** connected to a deformed geometry.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Compute rotation offset matrix. """;

// Get the up vector from input 1.
vector N = v@opinput1_N;

// Compute offset rotation with dihedral.
matrix3 rot_offset = dihedral(v@N, N);

// Store offset matrix.
3@rot = rot_offset;
```

## Primitive Dimensions
*Reference Code*: 9604585
> [!TIP]
> You can eventually accumulate the values to get a similar output as the measure node would return.

### measure
> [!IMPORTANT]
> **Mode:** Primitives.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Get dimension values from intrinsic attributes. """;

// Set volume, area and perimeter.
f@volume = primintrinsic(0, "measuredvolume", @primnum);
f@area = primintrinsic(0, "measuredarea", @primnum);
f@perimeter = primintrinsic(0, "measuredarea", @primnum);
```

## Primitive Type Checker
*Reference Code*: 914613
> [!NOTE]
> Checking the primitive type can allow you to treat the input geometry in a different way. For example, if the input geometry is packed, you'll probably require to apply a different treatment compared to an unpacked one. 

### prim_type
> [!IMPORTANT]
> **Mode:** Primitives.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Get primitive type from intrinsic attributes. """;

// Set primitive type.
s@prim_type = primintrinsic(0, "typename", @primnum);
```

## Jitter Points
*Reference Code*: 55905505

### jitter_pts
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a scatter node.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Jitter point position. """;

// Create function to rand values between -1 and 1.
function float randVal(int id; float seed){
    
    // Randomize values between -1 and 1 and return the value.
    float rand = fit01(rand(id+seed), -1, 1);
    return rand;
}

// Initialize jitter parameters.
int is_id = chi("id");
string id_attr = chs("id_attribute");
float scale = chf("scale");
float seed = chf("seed");
vector axis_scale = chv("axis_scale");

// Check if user wants to input id attribute.
int id = (is_id==1)?point(0, id_attr, @ptnum):@ptnum;

// Create directional vector.
vector dir = set(randVal(id, seed),
                 randVal(id, seed+1),
                 randVal(id, seed+2))*axis_scale;

// Set new position.
v@P+=(dir*scale);
```

## Check Point Inside Geometry
*Reference Code*: 69877006
> [!NOTE]
> Note that there are two different approaches to getting a similar result. The check_inside_pts_geo is limited because it requires a clean geometry without self-intersections, but it doesn't require an additional step converting the geometry into sdf. Meantime, the check_inside_pts_vol is more stable for self-intersecting geometries, but it requires the mentioned additional step.

### check_inside_pts_geo
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** connected to a bounding geometry.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Check if point is inside of geometry without volumes. """;

// Initialize primitive and parametric coords.
vector uvw; int prim;

// Capture primitive and parametric coords.
xyzdist(1, v@P, prim, uvw);

// Get position and normal ray intersection information.
vector pos = primuv(1, "P", prim, uvw);
vector norm = primuv(1, "N", prim, uvw);

// Compute direction vector between two points.
vector dir = normalize(v@P-pos);

// Check dot product between direction and normal.
float dot = dot(dir, norm);

// Create a group to contain the inside points.
if(dot<0) setpointgroup(0, "inside", @ptnum, 1);
```
### check_inside_pts_vol
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** connected to a bounding sdf volume.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Check if point is inside of geometry with volumes. """;

// Sample signed distance field volume.
float dist = volumesample(1, 0, v@P);

// Create a group to contain the inside points.
if(dist<0) setpointgroup(0, "inside", @ptnum, 1);
```

## Ray Minimum Distance
*Reference Code*: 77468763

### ray_min_dist
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** connected to a snapping geometry.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Ray minimum distance without volumes. """;

// Get closes position value.
vector minpos = minpos(1, v@P);

// Set position.
v@P = minpos; 
```

## Voxel Index And Rest
*Reference Code*: 8712550
> [!TIP]
> This code exports the necessary values to work with the volumerasterizelattice. You can potentially deform the exported points and the volumerasterizelattice will deform the volumes using the i@ix, i@iy, i@iz and v@rest attributes.

> [!NOTE]
> In order to create geometry with a volumewrangle, you'll need to turn on the Only Output Created Geometry in the Bindings tab.

### voxel_index
> [!IMPORTANT]
> **Mode:** Volume.
> - **Input 0:** connected to a density volume.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Compute rest and index values for volumes. """;

// Create point for each voxel.
int pt = addpoint(0, v@P);

// Convert position to rest position.
vector rest = volumepostoindex(0, 0, v@P);

// Store rest position.
setpointattrib(0, "rest", pt, rest);

// Create index for each axis and export it as integer.
setpointattrib(0, "ix", pt, int(rest.x));
setpointattrib(0, "iy", pt, int(rest.y));
setpointattrib(0, "iz", pt, int(rest.z));

// Set density to 0 to be able to generate geometry.
@density=0;
```
