Notes
----------------------------------------------------
May not work well with very large images.

May not work if the intersection point or coordinate is outside
the image. In other words, all corners of a paper or plane
should to be inside the image ideally.

May not work for all images, but I tweaked parameters such as thresholds to
work for most of the images I tested.

Some images may take a while to process if many local 
maxima are detected. For example, "a.png" took significantly more
time to process than "b.png". You can find these "a.png" and "b.png" in
the zip file.
-----------------------------------------------------
-----------------------------------------------------
Instructions
-----------------------------------------------------

1) Run main.m and 3 figures will pop up.
2) Figure 1 is the lines superimposed on the edge image.
3) Figure 2 is the 2D hough with the location of the lines
   marked with a red 'x'.
4) Figure 3 is the rectified image.

Others)
The following parameters can be changed for testing, but some changes
will probably break things more often than not. 
- line 5 and 6 are the images I showed on the PDF
- lines 9 is N or kernel size for edge detection
- line 10 is sigma value for edge detection
- line 11 and 12 is low and high threshold for edge detection
- line 52 is local maxima threshold
- line 79 is theta threshold
- line 87 is Rho threshold
- line 89 is length threshold
