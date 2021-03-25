# Display custom box (no rotation)

# Create NGL Viewer widget
view = nglview.NGLWidget()

# Add protein
p1 = view.add_component(reference)
p1.clear()
p1.add_cartoon(color='grey')

# Add box
view.shape.add_cylinder(box['box']['p1'], box['box']['p2'], [ 1, 0, 0 ], 0.1 , 'x')
view.shape.add_cylinder(box['box']['p1'], box['box']['p3'], [ 0, 1, 0 ], 0.1 , 'y')
view.shape.add_cylinder(box['box']['p1'], box['box']['p4'], [ 0, 0, 1 ], 0.1 , 'z')
view.shape.add_cylinder(box['box']['p2'], [box['box']['p2'][0], box['box']['p3'][1], box['box']['p2'][2]], [ 0, 0, 0 ], 0.1, '0')
view.shape.add_cylinder(box['box']['p2'], [box['box']['p2'][0], box['box']['p2'][1], box['box']['p4'][2]], [ 0, 0, 0 ], 0.1, '1')
view.shape.add_cylinder(box['box']['p3'], [box['box']['p2'][0], box['box']['p3'][1], box['box']['p2'][2]], [ 0, 0, 0 ], 0.1, '2')
view.shape.add_cylinder(box['box']['p3'], [box['box']['p3'][0], box['box']['p3'][1], box['box']['p4'][2]], [ 0, 0, 0 ], 0.1, '3')
view.shape.add_cylinder(box['box']['p4'], [box['box']['p2'][0], box['box']['p2'][1], box['box']['p4'][2]], [ 0, 0, 0 ], 0.1, '4')
view.shape.add_cylinder(box['box']['p4'], [box['box']['p3'][0], box['box']['p3'][1], box['box']['p4'][2]], [ 0, 0, 0 ], 0.1, '5')
view.shape.add_cylinder([box['box']['p2'][0], box['box']['p3'][1], box['box']['p2'][2]], [box['box']['p2'][0], box['box']['p3'][1], box['box']['p4'][2]], [ 0, 0, 0 ], 0.1, '6')
view.shape.add_cylinder([box['box']['p3'][0], box['box']['p3'][1], box['box']['p4'][2]], [box['box']['p2'][0], box['box']['p3'][1], box['box']['p4'][2]], [ 0, 0, 0 ], 0.1, '7')
view.shape.add_cylinder([box['box']['p2'][0], box['box']['p2'][1], box['box']['p4'][2]], [box['box']['p2'][0], box['box']['p3'][1], box['box']['p4'][2]], [ 0, 0, 0 ], 0.1, '8')

view
