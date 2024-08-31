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
vector index = volumepostoindex(0, 0, v@P);

// Store rest position.
setpointattrib(0, "rest", pt, v@P);

// Create index for each axis and export it as integer.
setpointattrib(0, "ix", pt, int(index.x));
setpointattrib(0, "iy", pt, int(index.y));
setpointattrib(0, "iz", pt, int(index.z));

// Set density to 0 to be able to generate geometry.
@density=0;
```

## Remove Array Duplicates
*Reference Code*: 31145437
> [!NOTE]
> This example is being created in order to show how we can remove duplicates from an array, not to output a specific example using actual geometry. It can be used in different contexts and situations where the input geometry might change.

### remove_duplicates
> [!IMPORTANT]
> **Mode:** Detail.
> - **Input 0:** no-connected.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Remove array duplicates. """;

// Create or input an array.
float smp_array[] = array(0,1,2,2,3,3,4,5);

// Initialize new array and loop through the vales of the old array. 
float new_smp_array[];
foreach(float val; smp_array){

    // Add value if it is not found already in the new array.
    if(find(new_smp_array, val)<0) append(new_smp_array, val);
}

// Store array attribute.
f[]@array = new_smp_array;
```

## Most Repeated Value Array
*Reference Code*: 30011921
> [!NOTE]
> This example is being created in order to show how we can get the most repeated element from an array, not to output a specific example using actual geometry. It can be used in different contexts and situations where the input geometry might change.

### most_repeated_value
> [!IMPORTANT]
> **Mode:** Detail.
> - **Input 0:** no-connected.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Find the most repeated element in the array. """;

// Create or input an array.
float smp_array[] = array(0,1,2,2,3,3,4,5);

// Initialize unique values.
float unique_vals[];

// Compute unique values from the array.
foreach(float item; smp_array){
    int checker = find(unique_vals, item);
    if(checker<0) append(unique_vals, item);
}

// Initialize counting.
float max_repeated = 0.0;
int max_items = 0;

// Check for the most repeated element in the array.
foreach(float val; unique_vals){
    
    // Get values from original array that matches current value.
    int items[] = find(smp_array, val);
    
    // Count how many items are in the array.
    int item_count = len(items);
    
    // If the count is bigger, it is the most repeated one until now.
    if(item_count>max_items){
        max_items = item_count;
        max_repeated = val;
    }
}

// Store max repeated attribute.
f@max_repeated = max_repeated;
```

## STMap Lens Shader
*Reference Code*: 77666090
> [!NOTE]
> This example needs to be tested inside a CVEX Shader Builder with a Inline Code node. It requires the creation of some external parameters to link the STMap file, aperture and focal length of the camera.

> [!NOTE]
> To be able to use it in Karma, you'll need to create an HDA and embbed the code into it. For mantra works without having to explicitly create that HDA.

> [!TIP]
> You can link the parameters from the CVEX Shader Builder to the actual camera values. 

### lens_shader
> [!IMPORTANT]
> **Mode:** VEX Shader.
> - **Input 0:** bind x axis.
> - **Input 1:** bind y axis.
> - **Input 2:** bind aspect.
> - **Input 3:** aperture value parameter node.
> - **Input 4:** focal length value parameter node.
> - **Input 5:** STMap file value parameter node.

> - **Output 0:** P as a Vector value.
> - **Output 1:** I as a Vector value.

``` c
""" Camera lens shader based on STMap. """;

// x and y is on a -1 to 1 space, we need to switch to ndc.
float ox = fit(x, -1, 1, 0, 1);
float oy = fit(y, -1, 1, 0, 1);

// Get color from the STMap and move distortion to center.
vector c = colormap(file, ox, oy)-0.5;

// Set postion.
$P = set(0, 0, 0);

// Set ray direction and length based on focus length and aperture.
$I = set(c.x, c.y / aspect, (fo/ap));
```

## Two Vector Intersect
*Reference Code*: 70358649
> [!NOTE]
> This example is being created in order to show how we can get the position information from two vectors, not to output a specific example using actual geometry. It can be used in different contexts and situations where the input geometry might change.

### vector_intersect
> [!IMPORTANT]
> **Mode:** Detail.
> - **Input 0:** no-connected.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Get intersection position between two vectors in a position. """;

// Initialize point positions.
vector origin_pt1 = set(0.5, 0, 0);
vector origin_pt2 = set(-0.5, 0, 0);

// Initialize direction vector.
vector dir_pt1 = set(-0.7, 0.7, 0);
vector dir_pt2 = set(0.7, 0.7, 0);

// Compute final position using the direction and origin position.
vector final_pt1 = origin_pt1 + dir_pt1 * 1.0e08;
vector final_pt2 = origin_pt2 + dir_pt2 * 1.0e08;

// Compute intersection position using the line-line equation by Eric Wolfgang Weisstein.
float den = (origin_pt1.x - final_pt1.x) *
            (origin_pt2.y - final_pt2.y) -
            (origin_pt1.y - final_pt1.y) *
            (origin_pt2.x - final_pt2.x);

float mag_vec = (origin_pt1.x * final_pt1.y - origin_pt1.y * final_pt1.x);
            
vector pos = set(mag_vec * (origin_pt2.x - final_pt2.x) + (origin_pt1.x - final_pt1.x) * mag_vec,
                 mag_vec * (origin_pt2.y - final_pt2.y) + (origin_pt1.y - final_pt1.y) * mag_vec,
                 0) / den;

// Store intersection position.
v@intersect_pos = pos;
```

## Create Circle
*Reference Code*: 32305122

### create_circle
> [!IMPORTANT]
> **Mode:** Details.
> - **Input 0:** no-connected.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Create circle. """;

// Initialize parameters.
int div = clamp(chi("divisions"), 2, int(1e09));
float uniform_scale = chf("uniform_scale");
vector2 radius = chu("radius");
int open = chi("open_arc");

// Initialize point array to store created points.
int pts[];

// Iterate for each division points.
for(int pt=0; pt<div; pt++){
    
    // Compute radians for each iteration.
    float rad = $PI*2/div*pt;
    
    // Compute position using radians and scaling values.
    vector pos = set(cos(rad)*radius.x,
                     sin(rad)*radius.y, 
                     0)*uniform_scale;
    
    // Create point and add number to point array.
    int p = addpoint(0, pos);
    append(pts, p);
}

// Append last number to close circle.
append(pts, pts[0]);

// If user wants open geo, create a polyline. Otherwise, create a closed poly.
(open==1)? addprim(0, "polyline", pts):addprim(0, "poly", pts);
```

## Point Attribute Transfer
*Reference Code*: 53451011
> [!NOTE]
> This code is able to handle the following datatypes: int, float, vector2 and vector. If you want to create transfer other datatypes, you can follow the pattern of the other ones.

> [!NOTE]
> There's another way to transfer attributes using point clouds and near points. The fastest approach with a large amount of points is this one as far as I tested.

> [!TIP]
> Note that this code has a falloff parameter that allows the user to modify how the data is tranferred using a ramp. In addition, there's a parameter to input the attributes that the user want to transfer.

### attribute_transfer
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry to copy attributes to.
> - **Input 1:** connected to a geometry to copy attributes from.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Point attribute transfer with falloff. """;

// Initialize transfer values.
int enable_falloff = chi("enable_falloff"); 
float max_dist = chf("max_distance");

// Get list of attributes and split based on spaces.
string list_attr = chs("point_attributes");
string attrs[] = split(list_attr, " ");

// Get closest point from the current point.
int nearpt = nearpoint(1, v@P);

// Get closest point position and compute distance.
vector nearpt_pos = point(1, "P", nearpt);
float dist = distance(v@P, nearpt_pos);

// If the point is not further than the maximum distance.
if(dist<max_dist){

    // Remap distance and store it as bias.
    float bias = fit(dist, 0, max_dist, 1, 0);
    bias = chramp("falloff_ramp", bias);
    
    // Loop for each of the input values.
    foreach(string att; attrs){
    
        // Get the type of the attribute and its size.
        int attrtype = attribtype(1, "point", att);
        int attrsize = attribsize(1, "point", att);
        
        // Check if the attribute is int and filter.
        if(attrtype==0) setpointattrib(0, att, @ptnum, int(point(1, att, nearpt)));

        // Check if the attribute is string and filter.
        if(attrtype==2) setpointattrib(0, att, @ptnum, string(point(1, att, nearpt)));
        
        // Check if the attribute is float.
        if(attrtype==1 && attrsize==1){
            
            // Filter float information.
            float filter = point(1, att, nearpt);
            
            // Check if user wants falloff and set the value.
            filter = (enable_falloff)? lerp(point(0, att, nearpt), filter, bias):filter;
            setpointattrib(0, att, @ptnum, filter);
        }
        
        // Check if the attribute is vector2.
        if(attrtype==1 && attrsize==2){
        
            // Filter vector2 information.
            vector2 filter = point(1, att, nearpt);
            
            // Check if user wants falloff and set the value.
            filter = (enable_falloff)? lerp(point(0, att, nearpt), filter, bias):filter;
            setpointattrib(0, att, @ptnum, filter);
        }
        
        // Check if the attribute is vector3.
        if(attrtype==1 && attrsize==3){
        
            // Filter vector3 information.
            vector filter = point(1, att, nearpt);
            
            // Check if user wants falloff and set the value.
            filter = (enable_falloff)? lerp(point(0, att, nearpt), filter, bias):filter;
            setpointattrib(0, att, @ptnum, filter);
            
        }
    }
}
```

## Find Equivalent Ptnum
*Reference Code*: 39128540
> [!NOTE]
> In this example I use an integer value, but you can use a string one if you'd like. In addition, this process can be done in different types of geometries.

### find_equiv_ptnum
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** connected to a geometry with at least a coinciding attribute with Input 0.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Find equivalent value and set the color using integer values. """;

// Find equivalent point using int attribute.
int equiv_pt = findattribval(1, "point", "id", i@id);

// Get attribute using equivalent point.
vector color = point(1, "Cd", equiv_pt);

// Set attribute.
v@Cd = color;
```

## Velocity Point Trail
*Reference Code*: 50993777

### point_trail
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry with v@v attribute.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Compute basic particle trail using velocity. """; 

// Compute velocity using fps.
vector comp_vel = v@v*f@TimeInc;

// Compute new position.
vector new_pos = v@P-comp_vel;

// Create point based on velocity.
int pt = addpoint(0, new_pos);

// Create polyline.
addprim(0, "polyline", @ptnum, pt);
```

## Convert Integer To String
*Reference Code*: 44470061
> [!NOTE]
> This example is being created in order to show how we can convert integers into strings, not to output a specific example using actual geometry. It can be used in different contexts and situations where the input geometry might change.

> [!TIP]
> You can use that snippet inputting the interations of a loop. For example, that would help you out to create proper naming for fracturing pieces.

### int_to_str
> [!IMPORTANT]
> **Mode:** Detail.
> - **Input 0:** no-connected.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Convert integer values into strings. """;

// Create or import the integer value.
int value = 2;

// Create a name including the converted integer value.
string name = "piece_" + itoa(value);

// Export string attribute.
s@name = name;
```

## Convert String To Float
*Reference Code*: 46533271
> [!NOTE]
> This example is being created in order to show how we can convert strings into floats, not to output a specific example using actual geometry. It can be used in different contexts and situations where the input geometry might change.

> [!NOTE]
> The str_to_flt_shash is unlikely to return the same values using different strings, but it can happen. The str_to_flt_utf won't return the same values if we use the raw utf value.

### str_to_flt_shash
> [!IMPORTANT]
> **Mode:** Detail.
> - **Input 0:** no-connected.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Convert string values into floating values (basic way). """;

// Create or import the integer value.
string name = "string";

// Create a random value using the string.
float value = rand(random_shash(name));

// Export float attribute.
f@name = value;
```
### str_to_flt_utf
> [!IMPORTANT]
> **Mode:** Detail.
> - **Input 0:** no-connected.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Convert string values into floating values (basic way). """;

// Create or import the integer value.
string name = "string";

// Decode string using utf8.
int values[] = decodeutf8(name);

// Initialize utfval.
string utfval="";

// Loop for each of the values.
foreach(int v; values){
    
    // Convert integer to string.
    string val = itoa(v);
    
    // Concat with previous iteration. 
    utfval+=val;
}

// Create random value converting string to integer.
float value = rand(atoi(utfval));

// Export float attribute.
f@value = value;
```

## Vector Between Positions
*Reference Code*: 4300925
> [!NOTE]
> Note that this example is a basic exemplification of how we can get the direction vector between two positions. This snippet usually comes with additional functionality that is not being implemented here.
 
### vector_between_pos
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** connected to a geometry.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Get direction vector between two points. """;

// Get point position to point to.
vector to_pos = point(1, "P", @ptnum);

// Get current point position to point from.
vector from_pos = v@P;

// Compute direction vector.
vector dir = to_pos-from_pos;

// Export direction attribute.
v@dir = dir;
```

## Hanging Catenary Wire
*Reference Code*: 63831594

### hanging_pts
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to reference points.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Create hanging catenary wires. """; 

// Get the parabola amplitude and resolution from lines.
float a = 1-ch('parabola_amplitude');                        
float res = chi('res')+1;                                  

// Initialize point array.
int pts[];

// Iterate for each input point. 
for(int i=1;i<npoints(0);i++){
    
    // Append current input point to the point array.
    append(pts, i-1);
    
    // Get current point position and next point position. 
    vector curr_pos = point(0,'P',i-1);
    vector next_pos = point(0,'P',i);
    
    // Iterate to create points for the line based on resolution factor.
    for(int b=1;b<res;b++){
        
        // Compute normalized line point value.
        float curveu = b/res;
        
        // Compute positions for the points.
        vector pos = curr_pos + ((next_pos-curr_pos)/res)*b;
        
        // Compute the amount of anchor displacement.
        float disp_corner = a * cosh((curveu-0.5)/a);
        
        // Compute the amount of center displacement. 
        float center_base = a * cosh((-0.5)/a);
        
        // Apply offsets/
        pos.y += disp_corner-center_base;
        
        // Create point.
        int pt = addpoint(0, pos);
        
        // Append point number to the point array.
        append(pts, pt);
    
    }                                   
}

// Append last point to the list and create the polyline.
append(pts, npoints(0)-1);
addprim(0,"polyline",pts);
```

## UDIM Connectivity
*Reference Code*: 93149288

### connectivity_udim
> [!IMPORTANT]
> **Mode:** Primitives.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Create attribute per UDIM. """;

// Get uv information using parametric uvs (just to not have to convert to points or vertices). 
vector2 uv = primuv(0, "uv", @primnum, {0.5,0.5});

// Round up the uv values.
uv = ceil(uv);

// Construct UDIM value convention. 
float uv_num = 1000 + (uv.x) + ((uv.y-1)*10);

// Export uv name attribute.
i@uv_name = int(uv_num);
```

## Convert Attribute To Group
*Reference Code*: 78120939
> [!NOTE]
> This example is being created in order to show how to convert a point string attribute into a point group. You can do the same process for the other geometry types, but you have to make sure that the attribute that you are using is a string.

> [!WARNING]
> Don't use attributes with a lot of different values... You don't want to get a lot of groups.  

### attr_to_grp
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

## Convert Group To Attribute
*Reference Code*: 29125957
> [!NOTE]
> This example is being created in order to show how to convert a point group into a point attribute. You can do the same process for the other geometry types.

> [!NOTE]
> Note that in this code you have two different methods: attr_creation and attr_name_array. Both output a really different result, so check what would be convinient for you.

> [!WARNING]
> If you want to use the attr_creation method is better to check for how many groups you have in your current geometry... You don't want to get a lot of attributes.  

### attr_name_array
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Convert attribute into groups. """;

// Get all point groups.
string grps[] = detailintrinsic(0, "pointgroups");

// Initialize point group array.
string pt_grps[];

// Iterate for each available point group.
foreach(string grp; grps){
    
    // Check if point in point group and append if so.
    if (inpointgroup(0, grp, @ptnum)) append(pt_grps, grp);
}

// Export array of point groups for current point.
s[]@grp_attr = pt_grps;
```
### attr_creation
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Convert attribute into groups. """;

// Get all point groups.
string grps[] = detailintrinsic(0, "pointgroups");

// Iterate for each available point group.
foreach(string grp; grps){
    
    // Check if point in point group and create attribute if so.
    if(inpointgroup(0, grp, @ptnum)) setpointattrib(0, grp, @ptnum, 1);
}
```

## Push Points Volume
*Reference Code*: 73825197

### push_points
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to points.
> - **Input 1:** connected to sdf volume.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Push Points Outisde Volumes. """;

// Compute signed distance field value.
float dist = volumesample(1, 0, v@P);

// Compute gradient from volume.
vector grad = volumegradient(1, 0, v@P);

// Push points.
if(dist<0.001) v@P -= normalize(grad)*dist;
```

## Peak Geometry
*Reference Code*: 92501258

### peak
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a geometry.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Peak object using normals. """;

// Get peak value.
float peak = chf("peak");

// Mult normals and peak value.
vector push = v@N*peak;

// Add peak.
v@P+=push;
```

## Color Normalized Age
*Reference Code*: 66781383
> [!NOTE]
> You can use the @nage built-in attribute, which contains the same value as the nage variable from the code.

> [!TIP]
> Change the color ramp user parameter to color in order to visualize the proper ramp.

### colorize
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to points with f@age and f@life attributes.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Normalize particle age. """;

// Compute normalized age.
float nage = f@age/f@life;

// Remap color using normalized age.
vector color = chramp("color", nage);

// Export color attribute.
v@Cd = color;
```

## Create Box
*Reference Code*: 93703926
> [!NOTE]
> Note that the code doesn't allow to do modifications to the geometry because it is intended to be a really default box. In case you want to translate, rotate or scale, you can do it in the code by applying some transformation matrix.

### create_box
> [!IMPORTANT]
> **Mode:** Detail.
> - **Input 0:** no-connected.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Create default box. """;

// Create reusable void function.
void createPrim(int pt0, pt1, pt2, pt3){
    
    // Create primitives using arrays.
    int prim_pts[] = array(pt0, pt1, pt2, pt3);
    addprim(0,'poly', prim_pts);
}

// Create points.
int pt0 = addpoint(0, {0.5,-0.5,0.5});
int pt1 = addpoint(0, {0.5,-0.5,-0.5});
int pt2 = addpoint(0, {-0.5,-0.5,-0.5});
int pt3 = addpoint(0, {-0.5,-0.5,0.5});
int pt4 = addpoint(0, {0.5,0.5,0.5});
int pt5 = addpoint(0, {0.5,0.5,-0.5});
int pt6 = addpoint(0, {-0.5,0.5,-0.5});
int pt7 = addpoint(0, {-0.5,0.5,0.5});

// Create primitives using void function.
createPrim(pt0, pt1, pt2, pt3);
createPrim(pt0, pt4, pt5, pt1);
createPrim(pt1, pt5, pt6, pt2);
createPrim(pt2, pt6, pt7, pt3);
createPrim(pt3, pt7, pt4, pt0);
createPrim(pt6, pt5, pt4, pt7);
```

## Create Tube
*Reference Code*: 61391039
> [!NOTE]
> Note that the code only allows you to do modifications to the columns because it is intended to be a really default tube. In case you want to translate, rotate or scale, you can do it in the code by applying some transformation matrix.

### create_tube
> [!IMPORTANT]
> **Mode:** Detail.
> - **Input 0:** no-connected.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Create default tube. """;

// Set columns.
int col = chi("columns");

// Initialize point number and point position arrays.
int all_pts[]; 
vector all_pos[];

// Loop 2 times to create reference winding faces.
for(int a=0; a<2;a++){
    
    // Loop col times to create the primitive points.
    for(int i=0; i<col; i++){
        
        // Compute degrees in radians.
        float rad = $PI*2/col;
        
        // Compute the x and y axes.
        float x = sin(rad*i);
        float z = cos(rad*i);
        
        // Compute positions using x and z. Use y to place points in height.
        vector pos = (a)? set(x,-1,z)/2 : set(x,1,z)/2;
        
        // Create points.
        int pt = addpoint(0, pos);
        
        // Append points and positions to arrays.
        append(all_pts, pt);
        append(all_pos, pos);
    } 
}

// Separate core primitive points into two groups.
int first_prim[] = all_pts[:col];
int second_prim[] = all_pts[col:];

// Create primitives. One with reversed windings.
addprim(0, "poly", reverse(first_prim));
addprim(0, "poly", second_prim);

// Iterate for each point of the first primitive.
for(int i=0; i<len(first_prim); i++){
    
    // Find next equivalent position. Set first position if current point 
    // is the last point in the first primitive.
    vector next_pos_equiv = (first_prim[i]==first_prim[-1])? all_pos[0]:all_pos[i+1];
    
    // Invert z axis.
    next_pos_equiv.y*=-1;
    
    // Find equivalent point index. all_pos index = all_pts value 
    int equiv_pt = find(all_pos, next_pos_equiv);
    
    // Canstruct winding order. 
    // Set first point of the first prim and last point of the second prim
    // if current point is the last point in the first primitive.
    int prim_pts[] = (first_prim[i]==first_prim[-1])? array(i, first_prim[0], equiv_pt, second_prim[-1]) : 
                                   array(first_prim[i], first_prim[i+1], equiv_pt, equiv_pt-1);
    
    // Create primitive using the winding order.
    addprim(0, "poly", prim_pts);
}

```

## Create Line
*Reference Code*: 58694772
> [!NOTE]
> Note that the code only allows you to do modifications to the length and points because it is intended to be a really default line. In case you want to translate, rotate or scale, you can do it in the code by applying some transformation matrix.

### create_line
> [!IMPORTANT]
> **Mode:** Detail.
> - **Input 0:** no-connected.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Create line. """;

// Get number of points and length for the line.
int points = chi("points");
float length = chf("length");

// Initialize point array.
int pts[];

// Iterate for each points.
for(int i=0; i<points; i++){
    
    // Compute height value.
    float height = length/(points-1);
    
    // Create position value.
    vector pos = set(0, height*i, 0);
    
    // Create point and append to point array list. 
    int pt = addpoint(0, pos);
    append(pts, pt);
}

// Create poly line using point array.
addprim(0, "polyline", pts);
```

## Carve Primitive
*Reference Code*: 46938032

### carve
> [!IMPORTANT]
> **Mode:** Primitives.
> - **Input 0:** connected to a polyline.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Carve curve. """;

// Include groom library.
#include <groom.h>

// Get the distance to keep.
float dist = chf('distance');

// Get intrinsic perimeter
float p = primintrinsic(0, 'measuredperimeter', @primnum);

// Clamp distance.
float trim = clamp(dist, 0, p);

// Adjust the primitive length.
adjustPrimLength(0, @primnum, p, trim);
```

## Remove By Speed Threshold
*Reference Code*: 9971578

### remove_speed
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to points with v@v attribute.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

``` c
""" Remove based on speed. """;

// Get speed threshold.
float speed_thr = chf("speed_threshold");

// Compute speed.
float speed = length(v@v);

// Remove based on speed.
if(speed<speed_thr)removepoint(0, @ptnum);
```
