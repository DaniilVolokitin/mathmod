import gmsh
import sys

gmsh.initialize()

gmsh.model.add("hollow_torus")

# R = distance from the center of the torus to the center of tube
# inner and outer radiuses = radiuses of the tube 
R = 1.0 
outer_radius = 0.35
inner_radius = 0.20

outer_torus = gmsh.model.occ.addTorus(0, 0, 0, R, outer_radius)
inner_torus = gmsh.model.occ.addTorus(0, 0, 0, R, inner_radius)

gmsh.model.occ.cut([(3, outer_torus)], [(3, inner_torus)])

gmsh.model.occ.synchronize()

# The wall thickness is 0.15. With mesh_size = 0.04 there are about
# 3-4 tetrahedra across the wall thickness.
mesh_size = 0.04
gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size)
gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)

gmsh.model.mesh.generate(3)

gmsh.write("hollow_torus.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
