import math
import os

import gmsh
import numpy as np
import vtk


class CalcMesh:
    def __init__(self, node_tags, node_coords, tetra_tags):
        # Координаты узлов: отдельная строка для каждой оси.
        self.nodes = np.array([
            node_coords[0::3],
            node_coords[1::3],
            node_coords[2::3],
        ], dtype=np.double)
        self.center = np.mean(self.nodes, axis=1)

        # Продольная ось этой STL-модели лежит по диагонали плоскости XY.
        diagonal = 1.0 / math.sqrt(2.0)
        x = self.nodes[0] - self.center[0]
        y = self.nodes[1] - self.center[1]
        self.along = (-x + y) * diagonal
        self.side = (x + y) * diagonal
        self.height = self.nodes[2] - self.center[2]
        self.twist_factor = self.along / np.max(np.abs(self.along))
        self.diagonal = diagonal

        tag_to_index = {int(tag): i for i, tag in enumerate(node_tags)}
        self.tetrahedra = np.array(
            [tag_to_index[int(tag)] for tag in tetra_tags], dtype=np.int64
        ).reshape(-1, 4)

        self.velocity = np.zeros_like(self.nodes)
        self.twist_angle = np.zeros(len(self.nodes[0]))
        self.time = 0.0
        self.period = 2.4
        self.amplitude = math.pi / 4.0 # 45 градусов
        self.update_fields()

    # концы поворачиваются в разные стороны.
    def update_fields(self):
        omega = 2.0 * math.pi / self.period
        phase = omega * self.time + math.pi / 4.0
        angle = self.amplitude * self.twist_factor * math.sin(phase)
        angle_rate = self.amplitude * self.twist_factor * omega * math.cos(phase)

        # Поворачиваем каждое поперечное сечение вокруг продольной оси.
        side = self.side * np.cos(angle) - self.height * np.sin(angle)
        height = self.side * np.sin(angle) + self.height * np.cos(angle)
        self.nodes[0] = self.center[0] + (-self.along + side) * self.diagonal
        self.nodes[1] = self.center[1] + (self.along + side) * self.diagonal
        self.nodes[2] = self.center[2] + height

        side_velocity = -height * angle_rate
        self.velocity[0] = side_velocity * self.diagonal
        self.velocity[1] = side_velocity * self.diagonal
        self.velocity[2] = side * angle_rate
        self.twist_angle = angle

    def move(self, tau):
        self.time += tau
        self.update_fields()

    # Записываем один кадр с тетраэдрами, скоростью и углом поворота.
    def snapshot(self, number, output_dir):
        grid = vtk.vtkUnstructuredGrid()
        points = vtk.vtkPoints()

        angle_array = vtk.vtkDoubleArray()
        angle_array.SetName("twist_angle")

        velocity_array = vtk.vtkDoubleArray()
        velocity_array.SetNumberOfComponents(3)
        velocity_array.SetName("velocity")

        for i in range(len(self.nodes[0])):
            points.InsertNextPoint(
                self.nodes[0, i], self.nodes[1, i], self.nodes[2, i]
            )
            angle_array.InsertNextValue(self.twist_angle[i])
            velocity_array.InsertNextTuple((
                self.velocity[0, i], self.velocity[1, i], self.velocity[2, i]
            ))

        grid.SetPoints(points)
        grid.GetPointData().AddArray(angle_array)
        grid.GetPointData().AddArray(velocity_array)

        for nodes in self.tetrahedra:
            tetra = vtk.vtkTetra()
            for j in range(4):
                tetra.GetPointIds().SetId(j, int(nodes[j]))
            grid.InsertNextCell(tetra.GetCellType(), tetra.GetPointIds())

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputDataObject(grid)
        writer.SetFileName(os.path.join(
            output_dir, "pig-spiral-step-%03d.vtu" % number
        ))
        writer.Write()


script_dir = os.path.dirname(os.path.abspath(__file__))
stl_path = os.path.join(
    script_dir, "..", "01_meshes", "stage2", "pig_lowpoly.stl"
)
output_dir = os.path.join(script_dir, "frames")

gmsh.initialize()
gmsh.model.add("pig_spiral")
gmsh.merge(stl_path)

# Замкнутая поверхность STL задаёт объём, который Gmsh заполнит тетраэдрами.
surfaces = gmsh.model.getEntities(2)
surface_tags = [surface[1] for surface in surfaces]
surface_loop = gmsh.model.geo.addSurfaceLoop(surface_tags)
gmsh.model.geo.addVolume([surface_loop])
gmsh.model.geo.synchronize()

gmsh.option.setNumber("Mesh.MeshSizeMin", 0.4)
gmsh.option.setNumber("Mesh.MeshSizeMax", 2.0)
gmsh.model.mesh.generate(3)

node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
element_types, _, element_node_tags = gmsh.model.mesh.getElements()

tetra_tags = None
for i in range(len(element_types)):
    if element_types[i] == 4:
        tetra_tags = element_node_tags[i]
        break

if tetra_tags is None:
    gmsh.finalize()
    raise RuntimeError("Gmsh did not create tetrahedra")

print("The model has %d nodes and %d tetrahedra" % (
    len(node_tags), len(tetra_tags) // 4
))
gmsh.finalize()

# 120 кадров с шагом 0.02.
os.makedirs(output_dir, exist_ok=True)
mesh = CalcMesh(node_tags, node_coords, tetra_tags)
steps = 120
tau = 0.02

for i in range(steps):
    mesh.snapshot(i, output_dir)
    mesh.move(tau)

print("Saved %d VTK snapshots to %s" % (steps, output_dir))
