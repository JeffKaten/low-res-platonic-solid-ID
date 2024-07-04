import math
import os
import numpy as np
from PIL import Image, ImageDraw

#PNGs will be saved to this folder:
file_location = #### Insert desired destination folder path here ####

#important constants: Golden ratio and its reciprocal are used to calculate coordinates of dodecahedron and icosahedron
phi1 = (1+ math.sqrt(5))/2
phi2 = (math.sqrt(5)-1)/2

#small parameter used when calling generate_initial_edge_list
E = 0.05

#number of PNGs produced per shape
output_size = 10000

#Draw line segments connecting specified pairs of points to an image
def draw_lines_to_image(line_coords_list, image):

    # Get a drawing context
    draw = ImageDraw.Draw(image)

    for line_coords in line_coords_list:
        # Draw the line segment
        draw.line(line_coords, fill=(0, 0, 0), width=1)


#Create a 40 pixel by 40 pixel PNG of an image of a specified list of line segments
def create_image_with_shape(line_coords_list,size=(40, 40), color=(255,255,255)):
    # Create a blank image with the specified size and background color
    image = Image.new('RGB', size, color)

    draw_lines_to_image(line_coords_list,image)

    image = image.convert('L')
    return image

#save the image as PNG
def save_image(output_path,image):
    image.save(output_path, 'PNG', compress_level=0)


#Given a list of vertices of a polyhedron and a list specifying which pairs of points form edges, perform specified series of transformations.
#Output will be the list of edges of the transformed polyhedron
def generate_edge_list(initial_vertex_list, initial_edge_list, diameter, stretch, rotation1, rotation2, shift1, shift2):
    rot1 = np.array(
        [[math.cos(rotation1), -math.sin(rotation1), 0], [math.sin(rotation1), math.cos(rotation1), 0], [0, 0, 1]])

    rot2 = np.array(
        [[1, 0, 0], [0, math.cos(rotation2), -math.sin(rotation2)], [0, math.sin(rotation2), math.cos(rotation2)]])

    str = np.array([[stretch / diameter, 0, 0], [0, stretch / diameter, 0], [0, 0, stretch / diameter]])

    proj = np.array([[1, 0], [0, 1], [0, 0]])

    shi = np.array([shift1, shift2])

    u = []
    for v in initial_vertex_list:
        u.append(tuple((v.dot(rot1).dot(rot2).dot(str).dot(proj) + shi).astype(int)))

    edge_list = []
    for edge in initial_edge_list:
        edge_list.append([u[edge[0]], u[edge[1]]])
    return edge_list


#In any regular convex polyhedron, the two vertices are connected if and only if their Euclidean distance is equal to edge_length. Due to rounding,
#we will instead check if an edge is present by checking if the distance between two vertices are less than edge_length+E for some small parameter E. For the
#purpose of this project, E=0.05 will suffice. This function returns the edge list of a regular convex polyhedron given its vertices and edge length.
def generate_initial_edge_list(initial_vertex_list, edge_length, E):
    n = len(initial_vertex_list)
    initial_edge_list = []
    for i in range(n):
        for j in range(i + 1, n):
            v = initial_vertex_list[i]
            u = initial_vertex_list[j]
            if math.sqrt((v[0] - u[0]) ** 2 + (v[1] - u[1]) ** 2 + (v[2] - u[2]) ** 2) < edge_length+E:
                initial_edge_list.append([i, j])
    return initial_edge_list


#Given a single polyhedron, generate a variety of PNGS of 2D projections of that polyhedron and place them in a folder.
#Folder name will be the name of the polyhedron. Use a python set() data type to ensure there are no duplicate images.
def generate_individual_dataset(initial_vertex_list,initial_edge_list,diameter,shape_name,file_location):
    folder = file_location +'\\'+ shape_name
    if not os.path.isdir(folder):
        os.makedirs(folder)
    i = 0
    image_data = set()
    while i<output_size:
        stretch = np.random.uniform(15.0, 35.0)
        rotation1 = np.random.uniform(0, math.pi)
        rotation2 = np.random.uniform(0, 2 * math.pi)
        shift1 = np.random.uniform(20 - (40 - stretch) / 2.5, 20 + (40 - stretch) / 2.5)
        shift2 = np.random.uniform(20 - (40 - stretch) / 2.5, 20 + (40 - stretch) / 2.5)
        my_shape = generate_edge_list(initial_vertex_list, initial_edge_list, diameter, stretch, rotation1, rotation2, shift1, shift2)

        this_file_name = folder + '\\' + shape_name + str(i) + '.png'

        my_image = create_image_with_shape(my_shape)

        #Only save the image if it is distinct from all other previously generated images
        binary_image = my_image.point(lambda x: 1 if x == 0 else 0, '1')
        my_tuple = tuple(list(binary_image.getdata()))
        if my_tuple not in image_data:
            image_data.add(my_tuple)
            save_image(this_file_name,my_image)
            i+=1
    print(shape_name + " dataset is created.")

#Given a dictionary of regular convex polyhedra data, generate all the data for each polyhedron type, organized into folders by polyhedron name within the specified location.
def generate_full_dataset(shapes_dict,file_location):
    for shape in shapes_dict:
        shape_vertex_list = shapes_dict[shape][0]
        shape_edge_list = generate_initial_edge_list(shape_vertex_list, shapes_dict[shape][1], E)
        generate_individual_dataset(shape_vertex_list, shape_edge_list, shapes_dict[shape][2], shape, file_location)

######################################################################################################################

#First, calculate coordinates of vertices of 5 platonic solids centered at 0 and record the edge length and diameter of each.
#Store info in dictionary of the form {shape_name: (vertex_list, edge_length, diameter)}

PlatonicSolids = {}

# Tetrahedra

#Initial vertices of tetrahedron
v1 = np.array([1,0,-math.sqrt(2)/2])
v2 = np.array([-1,0,-math.sqrt(2)/2])
v3 = np.array([0,1,math.sqrt(2)/2])
v4 = np.array([0,-1,math.sqrt(2)/2])
tetrahedron_vertex_list = [v1,v2,v3,v4]

tetrahedron_edge_length = 2

tetrahedron_diameter = 2

PlatonicSolids["tetrahedron"] = (tetrahedron_vertex_list,tetrahedron_edge_length,tetrahedron_diameter)


#Cubes
#Initial vertices of cube
cube_vertex_list = []
for i in [1,-1]:
    for j in [1,-1]:
        for k in [1,-1]:
            cube_vertex_list.append(np.array([i,j,k]))

cube_edge_length = 2
cube_diameter = math.sqrt(3)*2

PlatonicSolids["cube"] = (cube_vertex_list,cube_edge_length,cube_diameter)

#octahedra

#Initial vertices of octahedron
octahedron_vertex_list = []
for i in [1,-1]:
    octahedron_vertex_list.append(np.array([i,0,0]))
    octahedron_vertex_list.append(np.array([0,i,0]))
    octahedron_vertex_list.append(np.array([0,0,i]))

octahedron_edge_length = math.sqrt(2)
octahedron_diameter = 2

PlatonicSolids["octahedron"] = (octahedron_vertex_list,octahedron_edge_length,octahedron_diameter)

#dodecahedra

dodecahedron_vertex_list = []
for i in [1,-1]:
    for j in [1,-1]:
        dodecahedron_vertex_list.append(np.array([i * phi1, j * phi2, 0]))
        dodecahedron_vertex_list.append(np.array([0, i * phi1, j * phi2]))
        dodecahedron_vertex_list.append(np.array([j * phi2, 0, i * phi1]))
        for k in [1,-1]:
            dodecahedron_vertex_list.append(np.array([i,j,k]))

dodecahedron_edge_length = 2*phi2
dodecahedron_diameter = 2*math.sqrt(3)

PlatonicSolids["dodecahedron"] = (dodecahedron_vertex_list,dodecahedron_edge_length,dodecahedron_diameter)

#icosahedra

icosahedron_vertex_list = []

for i in [1,-1]:
    for j in [1,-1]:
        icosahedron_vertex_list.append(np.array([i,j*phi1,0]))
        icosahedron_vertex_list.append(np.array([0,i, j * phi1]))
        icosahedron_vertex_list.append(np.array([j * phi1,0,i]))

icosahedron_edge_length = 2
icosahedron_diameter = 2*math.sqrt(phi1**2+1)

PlatonicSolids["icosahedron"] = (icosahedron_vertex_list,icosahedron_edge_length,icosahedron_diameter)

#Generate the data
generate_full_dataset(PlatonicSolids,file_location)
