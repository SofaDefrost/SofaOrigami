import math
import json
import numpy

def readFold(filename):
    with open(filename, 'r') as myfile:
        data=myfile.read()


    # parse file
    foldPattern = json.loads(data)

    positions = foldPattern["vertices_coords"]
    scaledPositions = [[100*p[0],100*p[1],100*p[2]] for p in positions]
    positions = scaledPositions
    faces = foldPattern["faces_vertices"]
    edges = foldPattern["edges_vertices"]
    edgesAss = foldPattern["edges_assignment"]
    triangles = []
    for face in faces:
        print("face",face)
        if len(face)==4:
            triangles = triangles# + face[0:3] + face[2:4] + [face[0]]
        else:
            triangles = triangles + [x for x in face]
    print(triangles)    
    origamiEdges = []
    for i in range(len(edges)):
        if edgesAss[i] == 'B':
            print("Border !")
            for face in faces:                            
                if edges[i][0] in face and edges[i][1] in face:
                    for v in face:
                        if v not in edges[i]:
                            p1 = positions[edges[i][0]]
                            p2 = positions[edges[i][1]]
                            print('p1:', p1)
                            print('p2:', p2)
                            print('norm:',numpy.linalg.norm([p1[i]-p2[i] for i in range(3)]) )
                            origamiEdges = origamiEdges + edges[i] + [v,-1,0] + [500,5] + [numpy.linalg.norm([p1[i]-p2[i] for i in range(3)])]                                
                    break
        else:
            nodeDirect = -1
            nodeIndirect = -1
            for face in faces:
                if edges[i][0] in face and edges[i][1] in face:
                    posEdge0inFace = face.index(edges[i][0])
                    posEdge1inFace = face.index(edges[i][1])
                    if posEdge1inFace == ((posEdge0inFace + 1) % 3):
                        nodeDirect = face[((posEdge0inFace + 2) % 3)]
                    else:
                        print("face",face)
                        print("posEdge0inFace",posEdge0inFace)
                        print("posEdge0inFace",((posEdge0inFace + 1) % 3))
                        nodeIndirect = face[((posEdge0inFace + 1) % 3)]
                            
                    if nodeDirect >= 0 and nodeIndirect >= 0:
                        p1 = positions[edges[i][0]]
                        p2 = positions[edges[i][1]]
                        if edgesAss[i] == 'M':
                            origamiEdges = origamiEdges + edges[i] + [nodeDirect,nodeIndirect] + [1] + [500,5] + [numpy.linalg.norm([p1[i]-p2[i] for i in range(3)])] 
                        else:
                            origamiEdges = origamiEdges + edges[i] + [nodeDirect,nodeIndirect] + [-1]+ [500,5] + [numpy.linalg.norm([p1[i]-p2[i] for i in range(3)])] 
                        break
    return origamiEdges, triangles,edges,edgesAss, positions                                
