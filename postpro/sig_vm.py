import vtk

# Reader
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("../out/sigvm.vtk")
reader.Update()
# Mapper
element_mapper = vtk.vtkDataSetMapper()
element_mapper.SetInput(reader.GetOutput())
# Find minimum and maximum stresses to get a proper lookup table
scalar_data = reader.GetOutput()
min_value = scalar_data.GetScalarRange()[0]
max_value = scalar_data.GetScalarRange()[1]
lookup_table = vtk.vtkLookupTable()
lookup_table.Build()
lookup_table.SetRange(0, scalar_data.GetScalarRange()[1])
element_mapper.SetScalarRange(0, scalar_data.GetScalarRange()[1])
element_mapper.SetLookupTable(lookup_table)
print element_mapper.GetLookupTable().GetRange()
# Actor
actor = vtk.vtkActor()
actor.SetMapper(element_mapper)
scalar_bar = vtk.vtkScalarBarActor()
scalar_bar.SetLookupTable(lookup_table)
scalar_bar.SetHeight(0.5)
scalar_bar.SetPosition(0.02, 0.25)
scalar_bar.SetPosition2(0.10, 0.75)
# Renderer
renderer = vtk.vtkRenderer()
renderer.SetBackground(0.0, 0.0, 0.0)
renderer.AddActor(actor)
renderer.AddActor(scalar_bar)
# Renderer Window
render_window = vtk.vtkRenderWindow()
render_window.SetWindowName("sigma_vm")
render_window.SetSize(800, 800)
render_window.AddRenderer(renderer)
# Interactor
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(render_window)
# Initialize the interactor and start the rendering loop
interactor.Initialize()
render_window.Render()
interactor.Start()
