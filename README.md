# CircleFit
[![View CircleFit on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/112975-circlefit)

CircleFit is simply fitting noisy, weighted 2D data points into a single circle.
More information is written in the script "CircleFit.m"

# CircleFitV2

CircleFitV2 is an extended version of CircleFit, where we fit 2 parallel road lane by cocentric circles.
Unlike CircleFit, CircleFitV2 needs separate data points and weights for both lanes.
More deeper analysis is done in CircleFitV2 for my personal research(Vehicle Lane Mapping)

# CircleFitV3(Latest!)

CircleFitV3 is probably the final version of circle fit, providing 3 different types of modes for circle fitting

1. Free ends

2. One end fixed(To be accurate, the end points are constrained to a line)

3. Both ends fixed(Two ends have line constraints)

# CurveFit(Latest!)

CurveFit fits 2 parallel vehicle lane points with straight lines and arcs.

Works with LineFitV2, CircleFitV3, StringCreator, FindSegOrder files

StringCreator and FindSegOrder may need more modification in the future

## CircleFitV2 Example
![ArcParam1](https://user-images.githubusercontent.com/50237894/173098378-ec0e9892-1b34-4ad5-bb8e-5945d3240cc5.jpg)
> From data idx 1~2700

![ArcParam2](https://user-images.githubusercontent.com/50237894/173098612-d7ec5fd5-252f-4fd6-bdff-e42a19bce54b.jpg)
