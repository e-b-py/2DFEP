import vtk
#
# Source / Reader (in this case Source)
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("cantilever.vtk")
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()
reader.Update()
# 2nd reader for the connected moments
reader2 = vtk.vtkUnstructuredGridReader()
reader2.SetFileName("cantilever_moments.vtk")
reader2.Update()
# Try Hedgehod
hedge = vtk.vtkHedgeHog()
hedge.SetInput(reader.GetOutput())
hedge_mapper = vtk.vtkPolyDataMapper()
hedge_mapper.SetInputConnection(hedge.GetOutputPort())
hedge_actor = vtk.vtkActor()
hedge_actor.SetMapper(hedge_mapper)
hedge_actor.GetProperty().SetLineWidth(1.0)
hedge_actor.GetProperty().SetColor(0.0,0.0,1.0)
# Creat the mapper which maps the 3D cube data to geometrical primitives
data_mapper = vtk.vtkDataSetMapper()
data_mapper.SetInput(reader.GetOutput())
# Moment mapper
moment_mapper = vtk.vtkDataSetMapper()
moment_mapper.SetInput(reader2.GetOutput())
# Connect this to an actor
data_actor = vtk.vtkActor()
data_actor.SetMapper(data_mapper)
data_actor.GetProperty().SetLineWidth(2.0)
data_actor.GetProperty().SetColor(1.0,1.0,1.0)
# Moment actor
moment_actor = vtk.vtkActor()
moment_actor.SetMapper(moment_mapper)
moment_actor.GetProperty().SetLineWidth(1.0)
moment_actor.GetProperty().SetColor(0.0,0.0,1.0)
# Create a renderer and add the cube actor to it
renderer = vtk.vtkRenderer()
renderer.SetBackground(0.0, 0.0, 0.0) # Black
renderer.AddActor(data_actor)
renderer.AddActor(hedge_actor)
renderer.AddActor(moment_actor)
# Create the render window
render_window = vtk.vtkRenderWindow()
render_window.SetWindowName("VTK moment diagram")
render_window.SetSize(400, 400)
render_window.AddRenderer(renderer)
#1.0) Create an interactor
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(render_window)
# Initialize the interactor and start the rendering loop
interactor.Initialize()
render_window.Render()
interactor.Start()

