# Houdini VEX Snippets
## Initial Provisions
This repository is designated to be a place where I put some of the VEX snippets I've been using to fix, check, create, and manipulate information in different contexts. If something needs to be revisited, let me know so I can check for it and commit any of the requested modifications.

## Vector along curve
*Reference Code*: 72854126

### create_tangents_line
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

### blur_position
> [!IMPORTANT]
> **Mode:** Points.
> - **Input 0:** connected to a scatter node.
> - **Input 1:** no-connected.
> - **Input 2:** no-connected.
> - **Input 3:** no-connected.

> [!TIP]
> Create a for-loop and loop the attribute wrangle with feedback.

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
    
    // Add value to the pos vairiable.
    pos += point(0, "P", i); 
}

// Divide position by the amount of near points that you iterated.
// We add an additional value because of the initial value of the
// current position.
pos/=len(nearpts)+1;

// Set attribute value.
v@P = pos;
```
