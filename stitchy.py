import math, sys, argparse
from pathlib import Path
import numpy as np
import gdspy
import json

#input handling
print("stitchy - ebl writefield stitching suppression for raith eline")
parser = argparse.ArgumentParser(description="creates gdsii files with gradient on writefield boundaries")
parser.add_argument('-i', "--inputFile",        required="true", type=Path,  help="gdsii file input path")
parser.add_argument("-w", "--writefieldSize",   required="true", type=int,   help="writefieldsize in um")
parser.add_argument("-s", "--scaleWritefield",                   type=float, help="writefield scale, defaults to 1.02")
parser.add_argument("-c", "--cornerWritefield", required="true", type=ascii, help="writefield offset of lower left corner in um")
parser.add_argument("-d", "--doseSteps",                         type=int,   help="how many discrete dose values should the generator use")
parser.add_argument("-o", "--outputFile",       required="true", type=Path,  help="gdsii file output path")
args = parser.parse_args()

writefieldSize     = args.writefieldSize
writefieldLLOffset = np.fromstring(args.cornerWritefield.strip("'[]'"), sep=',', dtype=float)
writefieldScaler   = 1.02
if args.scaleWritefield is not None:
    writefieldScaler = args.scaleWritefield
doseSteps = 10
if args.doseSteps is not None:
    doseSteps = args.doseSteps
    
#gdsii parsing
gdspy.library.use_current_library = False
gdsiiIn = gdspy.GdsLibrary(infile=args.inputFile)
gdsiiOut = gdspy.GdsLibrary()
outTopCell = gdsiiOut.new_cell("TOP")
inTopCell = gdsiiIn.top_level()
polygons = inTopCell[0].get_polygonsets()
paths = inTopCell[0].get_paths()

def seperateIntoLayerAndDtypes(polyOrPathList):
    sameLayerAndAtt = lambda x,y: x.layers[0] == y.layers[0] and x.datatypes[0] == y.datatypes[0]
    res = []
    for element in polyOrPathList:
        l = next((x for x in res if sameLayerAndAtt(element, x[0])), [])
        if l == []:
            res.append(l)
        l.append(element)
    return res
    
polyAndPaths = seperateIntoLayerAndDtypes(polygons)
#scanPaths do not work, since they are converted to rectangles
#also bug in https://github.com/heitzmann/gdspy/blob/5234778068e8d85dde955056517dd5f95379b6d4/gdspy/operation.py#L75C28-L75C28 , should be p instead of args

def getIsodoseCoords(mindose, maxdose):
    print(f"mindose, maxdose {mindose}, {maxdose}")
    cornerSteps = 5
    boundaryThickness = writefieldSize*(writefieldScaler-1)
    minDoseDist = boundaryThickness*np.arccos(2*mindose-1)/np.pi
    maxDoseDist = boundaryThickness*np.arccos(2*maxdose-1)/np.pi
    print(np.arccos(2*mindose-1)/np.pi)
    print(np.arccos(2*maxdose-1)/np.pi)
    cornerPoints = np.array([[1,1], [1,-1], [-1,-1], [-1,1]])
    coords = []
    for cornerIdx, corner in enumerate(cornerPoints):
        for cornerAngle in np.linspace(0,np.pi/2,cornerSteps)+cornerIdx*np.pi/2:
            point2d = minDoseDist * np.array([np.sin(cornerAngle), np.cos(cornerAngle)])
            point2d += corner*writefieldSize*0.5
            point2d -= corner*boundaryThickness*0.5
            coords.append(point2d)
    if maxDoseDist > 1e-5:
        #add line to connect inside with outside path
        coords.append(coords[0].copy())
        coords.append(coords[0].copy())
        coords[-1][1] = maxDoseDist+writefieldSize*0.5
        coords[-1][1] -= boundaryThickness*0.5
        
        for cornerIdx, corner in enumerate(cornerPoints[::-1]):
            for cornerAngle in np.linspace(0,np.pi/2,cornerSteps)[::-1]+(3-cornerIdx)*np.pi/2:
                point2d = maxDoseDist * np.array([np.sin(cornerAngle), np.cos(cornerAngle)])
                point2d += corner*writefieldSize*0.5
                point2d -= corner*boundaryThickness*0.5
                coords.append(point2d)
    return coords

for pOrP in polyAndPaths:
    layAndDtypeDict = {
        "layer": pOrP[0].layers[0],
        "datatype": pOrP[0].datatypes[0]
    }
    if layAndDtypeDict["datatype"] == 0:
        layAndDtypeDict["datatype"] = 1000
    newCell = gdspy.Cell("temp")
    newCell.add(pOrP)
    boundingBox =  newCell.get_bounding_box() - writefieldLLOffset[None,:]
    boundingBox /= writefieldSize
    wfCentersX =  writefieldSize * np.arange(math.floor(boundingBox[0][0]), math.ceil(boundingBox[1][0])).astype("float")
    wfCentersX += writefieldLLOffset[0] + 0.5*writefieldSize
    wfCentersY =  writefieldSize * np.arange(math.floor(boundingBox[0][1]), math.ceil(boundingBox[1][1])).astype("float")
    wfCentersY += writefieldLLOffset[1] + 0.5*writefieldSize
    
    for wfCenter in np.stack(np.meshgrid(wfCentersX, wfCentersY), -1).reshape(-1, 2):
        doseList = np.linspace(0, 1, doseSteps)
        for doseIdx in range(len(doseList)-1):
            maxDose = doseList[doseIdx+1]
            if maxDose == 1.0:
                minDose = 1.0
            else:
                minDose = doseList[doseIdx]
            avgDose = 0.5*(maxDose+minDose)
            cutoutMask   = gdspy.Polygon(getIsodoseCoords(doseList[doseIdx], doseList[doseIdx+1]) + wfCenter)
            cutoutResult = gdspy.boolean(cutoutMask, pOrP, "and", layer=layAndDtypeDict["layer"], datatype=int(layAndDtypeDict["datatype"]*avgDose))
            print()
            if cutoutResult is not None:
                cutoutResult = cutoutResult.scale(1/writefieldScaler, center=wfCenter)
                outTopCell.add(cutoutResult)
                print(wfCenter)
gdsiiOut.write_gds(args.outputFile)
