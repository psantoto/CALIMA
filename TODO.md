## High priority
+ [] Correct license (LGPL instead of GPL)

## Error fixes
+ [] Make errors always positive
+ [] Check false positives (additional condition?)

## Improvements
+ [] Include parameters used in the run in the header of the output file
+ [] Change how CALIMA reads the input file so that the user does not have to
     write it on the command window
+ [] Check if any part of the code (especially the ellipse condition) can be
     further optimized.
+ [] Include an input for the fit box (right now it is the same as the source
     box, but I believe the source box wants to be smaller and the fit box,
     bigger). Make it so if no fit box is specified, the source box is used
     instead.
+ [] Highlight (somehow) the brightest peak.
+ [] Better way to report errors and warning to the user (?)
