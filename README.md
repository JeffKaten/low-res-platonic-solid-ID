# low-res-platonic-solid-ID
The following is part of my ongoing personal project in which I aim to implement and train a basic computer vision
model from scratch. This particular model will seek to categorize low resolution 2D images of 3D geometric objects.
In particular, it will sort images into categories based on which of the five Platonic solids (tetrahedron, cube, octahedron, icosahedron, and dodecahedron) the image contains.

The current repository represents the data set creation phase. 

There are currently two files:

#data-set-generation is a Python file which will create and label 50,000 PNG image files: 10,000 unique images of
each of the five Platonic solids. They are produced by applying pseudorandom transformations to a set of vertices 
of an initial object of each category and calculating the correct pairs of vertices to connect with edges. Next, we 
check for uniqueness before drawing to a low-resolution PNG using the PIL library.
