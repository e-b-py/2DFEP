import vtk

# Reader
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("../out/mesh.vtk")
#reader.SetFileName("gmsh_output.vtk")
reader.Update()
# Extract edges
edges = vtk.vtkExtractEdges()
edges.SetInput(reader.GetOutput())
edges.Update()
# Mapper
element_mapper = vtk.vtkDataSetMapper()
element_mapper.SetInput(reader.GetOutput())
edge_mapper = vtk.vtkDataSetMapper()
edge_mapper.SetInput(edges.GetOutput())
# Actor
actor = vtk.vtkActor()
actor.SetMapper(element_mapper)
actor.GetProperty().SetColor(0.4, 0.4, 0.4)
edge_actor = vtk.vtkActor()
edge_actor.SetMapper(edge_mapper)
edge_actor.GetProperty().SetLineWidth(2)
edge_actor.GetProperty().SetColor(1.0, 1.0, 1.0)
# Renderer
renderer = vtk.vtkRenderer()
renderer.SetBackground(0.0, 0.0, 0.0)
renderer.AddActor(actor)
renderer.AddActor(edge_actor)
# Renderer Window
render_window = vtk.vtkRenderWindow()
render_window.SetWindowName("Undeformed Mesh")
render_window.SetSize(600, 600)
render_window.AddRenderer(renderer)
# Interactor
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(render_window)
# Initialize thei nteractor and start the rendering loop
interactor.Initialize()
render_window.Render()
interactor.Start()
