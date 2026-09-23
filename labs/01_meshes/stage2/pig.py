import gmsh
import sys

gmsh.initialize()

gmsh.model.add("pig")

gmsh.merge("pig_lowpoly.stl")


surfaces = gmsh.model.getEntities(2)
surface_tags = [surface[1] for surface in surfaces]

surface_loop = gmsh.model.geo.addSurfaceLoop(surface_tags)
gmsh.model.geo.addVolume([surface_loop])

gmsh.model.geo.synchronize()

mesh_size = 2.0
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.4)
gmsh.option.setNumber("Mesh.MeshSizeMax", 2.0)

gmsh.model.mesh.generate(3)

gmsh.write("pig_mesh.msh")

if "-nopopup" not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()