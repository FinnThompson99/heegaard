import time

#The isomorphism signature that produces the gadget (G)
global GADGET_ISOSIG
GADGET_ISOSIG = 'fHLMabddeaaaa'

"""Normal coordinates for the normal surfaces of G that form all triangles and quadrilaterals
gadT0, gadT1, gadT2, gadT3: triangles around vertices 0, 1, 2, 3
gadQ0123, gadQ0213, gadQ0312: quadrilaterals partitioning vertices 01/23, 02/13, 03/12
"""
global gadT0, gadT1, gadT2, gadT3, gadQ0123, gadQ0213, gadQ0312
gadT0=[1,0,0,0,0,0,0, 1,0,0,0,0,0,0, 1,0,0,0,0,0,0, 1,0,0,0,0,0,0, 1,0,0,0,0,0,0]
gadT1=[0,1,0,0,0,0,0, 0,1,0,0,0,0,0, 0,1,0,0,0,0,0, 0,1,0,0,0,0,0, 0,1,0,0,0,0,0]
gadT2=[0,0,1,0,0,0,0, 0,0,1,0,0,0,0, 0,0,1,0,0,0,0, 0,0,1,0,0,0,0, 0,0,1,0,0,0,0]
gadT3=[0,0,0,1,0,0,0, 0,0,0,1,0,0,0, 0,0,0,1,0,0,0, 0,0,0,1,0,0,0, 0,0,0,1,0,0,0]
gadQ0123=[0,0,0,0,1,0,0, 0,0,0,0,1,0,0, 0,0,0,0,1,0,0, 0,0,0,0,1,0,0, 0,0,0,0,1,0,0]
gadQ0213=[0,0,0,0,0,1,0, 1,0,1,0,0,0,0, 0,0,0,0,0,1,0, 1,0,1,0,0,0,0, 0,0,0,0,0,1,0]
gadQ0312=[0,1,1,0,0,0,0, 0,1,1,0,0,0,0, 0,1,1,0,0,0,0, 0,1,1,0,0,0,0, 0,0,0,0,0,0,1]

"""Normal coordinates for the normal surfaces of G that correspond to annuli pieces
gadTubeT0T2, gadTubeT0T3, gadTubeT1T2, gadTubeT1T3: annuli between triangles 0/2, 0/3, 1/2, 1/3
gadTubeT0Q0312, gadTubeT3Q0312: annuli between quadrilateral q0312, and triangles 0, 3.
"""
gadTubeT0T2=[1,0,1,0,0,0,0, 0,0,0,0,0,1,0, 1,0,1,0,0,0,0, 0,0,0,0,0,1,0, 1,0,1,0,0,0,0]
gadTubeT0T3=[1,0,0,1,0,0,0, 1,0,0,1,0,0,0, 0,0,0,0,0,0,1, 1,0,0,1,0,0,0, 0,0,0,0,0,0,1]
gadTubeT1T2=[0,0,0,0,0,0,1, 0,0,0,0,0,0,1, 0,0,0,0,0,0,1, 0,0,0,0,0,0,1, 0,1,1,0,0,0,0]
gadTubeT1T3=[0,0,0,0,0,1,0, 0,0,0,0,0,1,0, 0,1,0,1,0,0,0, 0,1,0,1,0,0,0, 0,1,0,1,0,0,0]
gadTubeT0Q0312=[1,1,1,0,0,0,0, 0,1,0,0,0,1,0, 1,1,1,0,0,0,0, 0,1,0,0,0,1,0, 1,0,0,0,0,0,1]
gadTubeT3Q0312=[0,0,1,0,0,1,0, 0,0,1,0,0,1,0, 0,1,1,1,0,0,0, 0,1,1,1,0,0,0, 0,0,0,1,0,0,1]

"""Normal coordinates for the normal surfaces of G that correspond to octagon pieces
gadOctagonO0213, gadOctagonO0312: octagons partitioning vertices 02/13, 03/12
"""
gadOctagonO0213=[1,0,1,0,0,0,0, 0,0,0,0,0,1,0, 1,0,0,1,1,0,0, 0,1,0,1,0,0,0, 1,0,0,1,1,0,0]
gadOctagonO0312=[1,0,0,1,0,0,0, 1,0,0,1,0,0,0, 0,0,0,0,0,0,1, 1,0,0,1,0,0,0, 0,1,1,0,0,0,0]

def removeOctsGadget(surf, tetNum):
	"""Given a normal surface surf with an octagon in tetNum, inserts the gadget into surf's triangulation, and expresses the octagon in terms of normal surfaces, in T+G.
	This surface MUST be encoded in NS_AN_STANDARD.
	In general, this can be used for converting a surface into coordinates in T+G where there is no octagon in the given tetrahedron.
	
	PARAMETERS
	NormalSurface surf: the normal surface in the base triangulation
	int tetNum: the index in surf.triangulation() of the tetrahedron where the gadget is inserted. i.e., where the octagon is.
	
	RETURNS
	NormalSurface surfOctNormalised: the normal surface that is isomorphic to surf, but embedded in T+G
	"""
	T = surf.triangulation()

	#Saving the coordinates of surf
	surfVecOrig = list(surf.vector())
	
	#Save the coordinates in the tetrahedron from tetNum (note that here, using NS_AN_STANDARD, there are 10 coords per tet)
	surfVecTetNum = surfVecOrig[tetNum*10:tetNum*10+10]
	#Save the remainder of the coordinates
	surfVecRemainder = surfVecOrig[0:tetNum*10].extend(surfVecOrig[tetNum*10+10:])
	#Remove the last three coordinates in each block of 10 (these are the octagon coordinates) so that we can encode our surface with 7*T coordinates in NS_STANDARD.
	#surfVecStart is the 'chunk' of the normal coordinates of the surf that will be unchanged
	surfVecStart = []
	for i in range(len(surfVecRemainder)):
		if i%10 in [7,8,9]:
			surfVecStart.append(vec[surfVecRemainder])

	#Saving number of triangles and quadrilaterals and octagons in tetNum of the original surface
	(sT0, sT1, sT2, sT3, sQ0123, sQ0213, sQ0312, sO0123, sO0213, sO0312) = tuple(surfVecTetNum)
	
	#If we have an octagon, check its type, and then glue in the gadget accordingly.
	#o0213 and o0312 exist naturally in the gadget. Also, quads and octagons cannot co-exist, so if there are a nonzero number of octagons, we don't bother including sums of quadrilateral pieces.
	if sO0213!=0:
		newTri = gadgetPermGluer(T, tetNum, Perm4([0,1,2,3]))
		surfGadgetPart = [sT0*gadT0[i] + sT1*gadT1[i] + sT2*gadT2[i] + sT3*gadT3[i] + sO213*gadOctagonO0213[i] for i in range(35)]
	elif sO0312!=0:
		newTri = gadgetPermGluer(T, tetNum, Perm4([0,1,2,3]))
		surfGadgetPart = [sT0*gadT0[i] + sT1*gadT1[i] + sT2*gadT2[i] + sT3*gadT3[i] + sO312*gadOctagonO0312[i] for i in range(35)]
	elif sO0123!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([2,0,1,3]))
		surfGadgetPart = [sT*gadT0[i] + sT2*gadT1[i] + sT0*gadT2[i] + sT3*gadT3[i] + sO0123*gadOctagonO0213[i] for i in range(35)]
	#No octagons at all, just use an identity permutation. This could be included in the O0213 case, but this makes no difference. In this case, note we could have quadrilaterals. The quads aren't split into cases because there can only be one type of them in a tetrahedron.
	else:
		newTri = gadgetPermGluer(T, tetNum, Perm4([0,1,2,3]))
		surfGadgetPart = [sT0*gadT0[i] + sT1*gadT1[i] + sT2*gadT2[i] + sT3*gadT3[i] + sQ0123*gadQ0123[i] + sQ0213*gadQ0213[i] + sQ0312*gadQ0312[i] for i in range(35)]
	surfVecNew = surfVecStart + surfGadgetPart
	surfOctNormalised = NormalSurface(newTri, NS_STANDARD, surfVecNew)
	return surfOctNormalised

def surfaceConverter(surf, PERM, tetNum):
	"""Given a normal surface surf in a triangulation T, a tetrahedron and a permutation for a gadget gluing, return the equivalent surface that exists in the triangulation T+G, as defined by the permutation and tetrahedron.
	This surface MUST be encoded in NS_STANDARD.
	(If this has octagons in tetNum, use removeOctsGadget. If there are octagons in other tets, it does not make sense to convert the surface in tetNum, as this function is intended for use with adding surf to a tubed surface.)
	
	PARAMETERS
	NormalSurface surf: the normal surface in the base triangulation
	Permutation<4> PERM: the permutation of the gadget gluing.
	int tetNum: the index in surf.triangulation() of the tetrahedron where the gadget is inserted.
	
	RETURNS
	NormalSurface surfInTplusG: the normal surface that is isomorphic to surf, but embedded in T+G
	"""
	T = surf.triangulation()
	
	#Saving the coordinates of surf
	surfVecOrig = list(surf.vector())
	
	#The triangulation resulting from gluing the gadget in
	TplusG = gadgetPermGluer(T, tetNum, PERM)
	
	#Save the coordinates in tetrahedra before, in and after tetNum. The surface is already encoded in NS_STANDARD so we don't need to remove 0s in the coordinates.
	(surfVecOrig1, surfVecOrig2, surfVecOrig3) = (surfVecOrig[0:tetNum*7], surfVecOrig[tetNum*7:tetNum*7+7], surfVecOrig[tetNum*7+7:])
	#The 'chunk' of the normal coordinates of surf that are unchanged
	surfVecStart = surfVecOrig1 + surfVecOrig3
	#Saving number of triangles and quadrilaterals in tetNum of the original surface
	(sT0, sT1, sT2, sT3, sQ0123, sQ0213, sQ0312) = tuple(surfVecOrig2)
	if PERM == Perm4([0,1,2,3]):
		surfPortionInG = [sT0*gadT0[i] + sT1*gadT1[i] + sT2*gadT2[i] + sT3*gadT3[i] + sQ0123*gadQ0123[i] + sQ0213*gadQ0213[i] + sQ0312*gadQ0312[i] for i in range(35)]
	elif PERM == Perm4([0,3,1,2]):
		surfPortionInG = [sT0*gadT0[i] + sT1*gadT3[i] + sT2*gadT1[i] + sT3*gadT2[i] + sQ0123*gadQ0312[i] + sQ0213*gadQ0123[i] + sQ0312*gadQ0213[i] for i in range(35)]
	elif PERM == Perm4([1,0,3,2]):
		surfPortionInG = [sT0*gadT1[i] + sT1*gadT0[i] + sT2*gadT3[i] + sT3*gadT2[i] + sQ0123*gadQ0123[i] + sQ0213*gadQ0213[i] + sQ0312*gadQ0312[i] for i in range(35)]
	elif PERM == Perm4([1,2,0,3]):
		surfPortionInG = [sT0*gadT1[i] + sT1*gadT2[i] + sT2*gadT0[i] + sT3*gadT3[i] + sQ0123*gadQ0312[i] + sQ0213*gadQ0123[i] + sQ0312*gadQ0213[i] for i in range(35)]
	elif PERM == Perm4([2,0,1,3]):
		surfPortionInG = [sT0*gadT2[i] + sT1*gadT0[i] + sT2*gadT1[i] + sT3*gadT3[i] + sQ0123*gadQ0213[i] + sQ0213*gadQ0312[i] + sQ0312*gadQ0123[i] for i in range(35)]
	elif PERM == Perm4([3,1,0,2]):
		surfPortionInG = [sT0*gadT3[i] + sT1*gadT1[i] + sT2*gadT0[i] + sT3*gadT2[i] + sQ0123*gadQ0213[i] + sQ0213*gadQ0312[i] + sQ0312*gadQ0123[i] for i in range(35)]
	surfPortionInG = [sT0*gadT0[i] + sT1*gadT1[i] + sT2*gadT2[i] + sT3*gadT3[i] + sQ0123*gadQ0123[i] + sQ0213*gadQ0213[i] + sQ0312*gadQ0312[i] for i in range(35)]
	surfVecNew = surfVecStart + surfPortionInG
	surfNew = NormalSurface(TplusG, NS_STANDARD, surfVecNew)
	return surfNew

def surfaceGadgetCombiner(triRaw, tetNum, surf):
	"""Given a normal surface surf in the triangulation triRaw, we return all possible surfaces with tubes placed in tetrahedron tetNum. This assumes that s.eulerChar()<0 is already checked.
	The surface must be encoded in NS_STANDARD.
	
	PARAMETERS
	Triangulation3 triRaw: the full triangulation
	int tetNum: the index in triRaw of the tetrahedron where tubes are to be inserted
	NormalSurface surf: the normal surface to add tubes to
	
	RETURNS
	List<NormalSurface> allTubings: the list of all normal surfaces with tubes inserted in tetNum
	"""
	
	allTubings = []
	
	#Local copy of the triangulation (for safety)
	tri = Triangulation3(triRaw)

	#Saving the coordinates of surf
	surfVecOrig = list(surf.vector())
	
	#Save the coordinates in tetrahedra before, in and after tetNum
	(surfVecOrig1, surfVecOrig2, surfVecOrig3) = (surfVecOrig[0:tetNum*7], surfVecOrig[tetNum*7:tetNum*7+7], surfVecOrig[tetNum*7+7:])
	
	#We will create surfGadgetPart in the for loop then append it to surfVecStart to make surfVecNew. This is just the pre-existing surface parts outside of tetNum.
	#When the gadget is added to the triangulation, its tetrahedra are the last five of the triangulation, so surfGadgetPart will be at the end of the tubed surfaces' coordinates
	surfVecStart = surfVecOrig1 + surfVecOrig3
	
	#Saving number of triangles and quadrilaterals in tetNum of the original surface
	(sT0, sT1, sT2, sT3, sQ0123, sQ0213, sQ0312) = tuple(surfVecOrig2)

	#In each case below, if there are sufficient triangle/quadrilateral pieces present to form a tube, then a pre-chosen permutation and surface pairing are used to create the normal coordinates of the surface in those five tetrahedra. Also, the new triangulation formed from this gluing is produced.
	#tri0 tri1 tube
	if sT0!=0 and sT1!=0 and sQ0213==0 and sQ0312==0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([2,0,1,3]))
		surfGadgetPart = [(sT1-1)*gadT0[i] + sT2*gadT1[i] + (sT0-1)*gadT2[i] + sT3*gadT3[i] + sQ0123*gadQ0213[i] + gadTubeT0T2[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri0 tri2 tube
	if sT0!=0 and sT2!=0 and sQ0123==0 and sQ0312==0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([0,1,2,3]))
		surfGadgetPart = [(sT0-1)*gadT0[i] + sT1*gadT1[i] + (sT2-1)*gadT2[i] + sT3*gadT3[i] + sQ0213*gadQ0213[i] + gadTubeT0T2[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri0 tri3 tube
	if sT0!=0 and sT3!=0 and sQ0123==0 and sQ0213==0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([0,1,2,3]))
		surfGadgetPart = [(sT0-1)*gadT0[i] + sT1*gadT1[i] + sT2*gadT2[i] + (sT3-1)*gadT3[i] + sQ0312*gadQ0312[i] + gadTubeT0T3[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri1 tri2 tube
	if sT1!=0 and sT2!=0 and sQ0123==0 and sQ0213==0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([0,1,2,3]))
		surfGadgetPart = [sT0*gadT0[i] + (sT1-1)*gadT1[i] + (sT2-1)*gadT2[i] + sT3*gadT3[i] + sQ0312*gadQ0312[i] + gadTubeT1T2[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri1 tri3 tube
	if sT1!=0 and sT3!=0 and sQ0123==0 and sQ0312==0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([0,1,2,3]))
		surfGadgetPart = [sT0*gadT0[i] + (sT1-1)*gadT1[i] + sT2*gadT2[i] + (sT3-1)*gadT3[i] + sQ0213*gadQ0213[i] + gadTubeT1T3[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri2 tri3 tube
	if sT2!=0 and sT3!=0 and sQ0213==0 and sQ0312==0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([2,0,1,3]))
		surfGadgetPart = [sT1*gadT0[i] + (sT2-1)*gadT1[i] + sT0*gadT2[i] + (sT3-1)*gadT3[i] + sQ0123*gadQ0213[i] + gadTubeT1T3[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri0 quad0123 tube
	if sT0!=0 and sQ0123!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([0,3,1,2]))
		surfGadgetPart = [(sT0-1)*gadT0[i] + sT2*gadT1[i] + sT3*gadT2[i] + sT1*gadT3[i] + (sQ0123-1)*gadQ0312[i] + gadTubeT0Q0312[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri0 quad0213 tube
	if sT0!=0 and sQ0213!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([3,1,0,2]))
		surfGadgetPart = [sT2*gadT0[i] + sT1*gadT1[i] + sT3*gadT2[i] + (sT0-1)*gadT3[i] + (sQ0213-1)*gadQ0312[i] + gadTubeT3Q0312[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri0 quad0312 tube
	if sT0!=0 and sQ0312!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([0,1,2,3]))
		surfGadgetPart = [(sT0-1)*gadT0[i] + sT1*gadT1[i] + sT2*gadT2[i] + sT3*gadT3[i] + (sQ0312-1)*gadQ0312[i] + gadTubeT0Q0312[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri1 quad0123 tube
	if sT1!=0 and sQ0123!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([0,3,1,2]))
		surfGadgetPart = [sT0*gadT0[i] + sT2*gadT1[i] + sT3*gadT2[i] + (sT1-1)*gadT3[i] + (sQ0123-1)*gadQ0312[i] + gadTubeT3Q0312[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri1 quad0213 tube
	if sT1!=0 and sQ0213!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([2,0,1,3]))
		surfGadgetPart = [(sT1-1)*gadT0[i] + sT2*gadT1[i] + sT0*gadT2[i] + sT3*gadT3[i] + (sQ0213-1)*gadQ0312[i] + gadTubeT0Q0312[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri1 quad0312 tube
	if sT1!=0 and sQ0312!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([1,0,3,2]))
		surfGadgetPart = [(sT1-1)*gadT0[i] + sT0*gadT1[i] + sT3*gadT2[i] + sT2*gadT3[i] + (sQ0312-1)*gadQ0312[i] + gadTubeT0Q0312[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri2 quad0123 tube
	if sT2!=0 and sQ0123!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([1,2,0,3]))
		surfGadgetPart = [(sT2-1)*gadT0[i] + sT0*gadT1[i] + sT1*gadT2[i] + sT3*gadT3[i] + (sQ0123-1)*gadQ0312[i] + gadTubeT0Q0312[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri2 quad0213 tube
	if sT2!=0 and sQ0213!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([3,1,0,2]))
		surfGadgetPart = [(sT2-1)*gadT0[i] + sT1*gadT1[i] + sT3*gadT2[i] + sT0*gadT3[i] + (sQ0213-1)*gadQ0312[i] + gadTubeT0Q0312[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri2 quad0312 tube
	if sT2!=0 and sQ0312!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([1,0,3,2]))
		surfGadgetPart = [sT1*gadT0[i] + sT0*gadT1[i] + sT3*gadT2[i] + (sT2-1)*gadT3[i] + (sQ0312-1)*gadQ0312[i] + gadTubeT3Q0312[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri3 quad0123 tube
	if sT3!=0 and sQ0123!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([1,2,0,3]))
		surfGadgetPart = [sT2*gadT0[i] + sT0*gadT1[i] + sT1*gadT2[i] + (sT3-1)*gadT3[i] + (sQ0123-1)*gadQ0312[i] + gadTubeT3Q0312[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri3 quad0213 tube
	if sT3!=0 and sQ0213!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([2,0,1,3]))
		surfGadgetPart = [sT1*gadT0[i] + sT2*gadT1[i] + sT0*gadT2[i] + (sT3-1)*gadT3[i] + (sQ0213-1)*gadQ0312[i] + gadTubeT3Q0312[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	#tri3 quad0312 tube
	if sT3!=0 and sQ0312!=0:
		newTri = gadgetPermGluer(tri, tetNum, Perm4([0,1,2,3]))
		surfGadgetPart = [sT0*gadT0[i] + sT1*gadT1[i] + sT2*gadT2[i] + (sT3-1)*gadT3[i] + (sQ0312-1)*gadQ0312[i] + gadTubeT3Q0312[i] for i in range(35)]
		surfVecNew = surfVecStart + surfGadgetPart
		surfNew = NormalSurface(newTri, NS_STANDARD, surfVecNew)
		allTubings.append(surfNew)
	return allTubings

#perm MUST be 0123, 0312, 1032, 1203, 2013 or 3102
#this will be a privately used function anyway.
def gadgetPermGluer(triRaw, tetNum, PERM):
	"""Given a triangulation triRaw and a permutation PERM, this returns the triangulation where tetrahedron tetNum is removed and replaced with the gadget, glued according to PERM.
	
	PARAMETERS
	Triangulation3 triRaw: the full triangulation
	int tetNum: the index in triRaw of the tetrahedron where the gadget is to be inserted
	Permutation<4> PERM: the permutation to glue the gadget with. This is PERM:tet->gad.
	
	RETURNS
	Triangulation3 tri: the tetrahedron with the gadget glued into tetNum, according to PERM
	"""
	tri = Triangulation3(triRaw)

	#Define the gadget using fromGluings for future-proofing
	gad = Triangulation3.fromGluings(5, [
	[0,0,1,Perm4(0,1,2,3)],[0,2,1,Perm4(0,1,2,3)],
	[ 0, 3, 2, Perm4(0,1,2,3) ], [ 1, 1, 3, Perm4(0,1,2,3) ],
	[ 1, 3, 3, Perm4(0,1,2,3) ], [ 2, 0, 3, Perm4(0,1,2,3) ],
	[ 2, 1, 4, Perm4(0,1,2,3) ], [ 2, 2, 4, Perm4(0,1,2,3) ]])
	tri.insertTriangulation(gad)
	tet = tri.tetrahedron(tetNum)

	#This grabs the info for the faces that will be glued to the gadget.
	#tABC refers to what is adjacent to face ABC of tetrahedron tetNum.
	#Perms here are P:tet->tABC. Take P.inverse() to get the perm tABC->tet
	(t012simp, t012face, t012perm) = (tet.adjacentSimplex(3), tet.adjacentFacet(3), tet.adjacentGluing(3))
	(t013simp, t013face, t013perm) = (tet.adjacentSimplex(2), tet.adjacentFacet(2), tet.adjacentGluing(2))
	(t023simp, t023face, t023perm) = (tet.adjacentSimplex(1), tet.adjacentFacet(1), tet.adjacentGluing(1))
	(t123simp, t123face, t123perm) = (tet.adjacentSimplex(0), tet.adjacentFacet(0), tet.adjacentGluing(0))
	
	#Saving these for iteration:
	tSimps = [t012simp, t013simp, t023simp, t123simp]
	tFaces = [t012face, t013face, t023face, t123face]
	tPerms = [t012perm, t013perm, t023perm, t123perm]

	#This grabs the gadget from tri+gad (should be the only thing with boundary)
	gadBC = tri.boundaryComponent(0)
	(g012simp, g012face) = (gadBC.triangle(3).front().tetrahedron(), gadBC.triangle(3).front())
	(g013simp, g013face) = (gadBC.triangle(2).front().tetrahedron(), gadBC.triangle(2).front())
	(g023simp, g023face) = (gadBC.triangle(0).front().tetrahedron(), gadBC.triangle(0).front())
	(g123simp, g123face) = (gadBC.triangle(1).front().tetrahedron(), gadBC.triangle(1).front())

	#selfGluer[PERM][(tri1, tri2)]=[gadSimp1, gadFace1, gadSimp2] describes how, in cases where tet.triangle(tri1)==tet.triangle(tri2), the gadget should be glued. We are gluing from gadSimp1 to gadSimp2, along triangle gadFace1 of gadSimp1.
	#Here, (tri1, tri2) always has tri1>tri2 - the code using this loops over triangle numbers accordingly.
	selfGluer = {
			Perm4([0,1,2,3]): {(3,2): (g012simp, 3, g013simp),
							   (3,1): (g012simp, 3, g023simp),
							   (3,0): (g012simp, 3, g123simp),
							   (2,1): (g013simp, 2, g023simp),
							   (2,0): (g013simp, 2, g123simp),
							   (1,0): (g023simp, 1, g123simp)},
			Perm4([0,3,1,2]): {(3,2): (g013simp, 2, g023simp),
							   (3,1): (g013simp, 2, g012simp),
							   (3,0): (g013simp, 2, g123simp),
							   (2,1): (g023simp, 1, g012simp),
							   (2,0): (g023simp, 1, g123simp),
							   (1,0): (g012simp, 3, g123simp)},
			Perm4([1,0,3,2]): {(3,2): (g013simp, 2, g012simp),
							   (3,1): (g013simp, 2, g123simp),
							   (3,0): (g013simp, 2, g023simp),
							   (2,1): (g012simp, 3, g123simp),
							   (2,0): (g012simp, 3, g023simp),
							   (1,0): (g123simp, 0, g023simp)},
			Perm4([1,2,0,3]): {(3,2): (g012simp, 3, g123simp),
							   (3,1): (g012simp, 3, g013simp),
							   (3,0): (g012simp, 3, g023simp),
							   (2,1): (g123simp, 0, g013simp),
							   (2,0): (g123simp, 0, g023simp),
							   (1,0): (g013simp, 2, g023simp)},
			Perm4([2,0,1,3]): {(3,2): (g012simp, 3, g023simp),
							   (3,1): (g012simp, 3, g123simp),
							   (3,0): (g012simp, 3, g013simp),
							   (2,1): (g023simp, 1, g123simp),
							   (2,0): (g023simp, 1, g013simp),
							   (1,0): (g123simp, 0, g013simp)},
			Perm4([3,1,0,2]): {(3,2): (g013simp, 2, g123simp),
							   (3,1): (g013simp, 2, g023simp),
							   (3,0): (g013simp, 2, g012simp),
							   (2,1): (g123simp, 0, g023simp),
							   (2,0): (g123simp, 0, g012simp),
							   (1,0): (g023simp, 1, g012simp)}
			}
	
	#Check for self-gluings
	selfJoin = ()
	#for i1 in [3,2,1,0]
	for i1 in range(3,-1,-1):
		#for i2 in [i1,...,0]
		for i2 in range(i1-1,-1,-1):
			if tet.triangle(i1) == tet.triangle(i2):
				#note this is going to be 3,2, 3,1, 3,0, 2,1, 2,0 or 1,0. If one is true, any other wont occur (if there are two face pairings, then we have a one-tet triangulation.)
				selfJoin = (i1, i2)
	
	#Unglue and reglue each face of tet, depending on the chosen PERM.
	#At each stage, if selfJoin is non-empty, then the gadget is glued to itself accordingly. If triAtoBperm:triA->triB and PERM:triA->gadA and PERM:triB->gadB, then the permutation mapping gadA to gadB is PERM^{-1}(triAtoBperm(PERM)), which is PERM * triAtoBperm * PERM^{-1}
	#If selfJoin is empty, the face adjacent to tet (the one we are removing) is glued to the gadget. If tABCperm:tet->tABC and PERM:tet->gad, then the permutation mapping tABC to gad is tABCperm^{-1}(PERM), which is PERM * tABCperm^{-1}
	if PERM == Perm4([0,1,2,3]):
		gSimps = [g012simp, g013simp, g023simp, g123simp]
	elif PERM == Perm4([0,3,1,2]):
		gSimps = [g013simp, g023simp, g012simp, g123simp]
	elif PERM == Perm4([1,0,3,2]):
		gSimps = [g013simp, g012simp, g123simp, g023simp]
	elif PERM == Perm4([1,2,0,3]):
		gSimps = [g012simp, g123simp, g013simp, g023simp]
	elif PERM == Perm4([2,0,1,3]):
		gSimps = [g012simp, g023simp, g123simp, g013simp]
	elif PERM == Perm4([3,1,0,2]):
		gSimps = [g013simp, g123simp, g023simp, g012simp]
	
	#If there is a self-joining,
	if selfJoin != ():
		#Save the triangles where A is glued to B
		(triA, triB) = selfJoin
		#Save the permutation from triA to triB
		triA_to_triB = tet.adjacentGluing(triA)
		#Look-up the gluing method from selfGluer
		(gadFROMsimp, gadFROMface, gadTOsimp) = selfGluer[PERM][selfJoin]
		#Join the faces accordingly
		gadFROMsimp.join(gadFROMface, gadTOsimp, PERM*triA_to_triB*(PERM.inverse()))
		
	#For each other triangle of tetNum that isn't part of selfJoin...
	for i in range(4):
		triNumI = 3 - i
		if not(triNumI in selfJoin):
			#Separate this triangle from the rest of the triangulation
			tet.unjoin(triNumI)
			#Glue this triangle to the gadget accordingly
			tSimps[i].join(tFaces[i], gSimps[i], PERM*(tPerms[i].inverse()))
	
	#Remove tet now that it's isolated
	tri.removeTetrahedron(tet)
	return tri

def filteredSurfaceGenerator(tri, eulerChars, vertex, almNorm):
	"""This generates the subset of (almost) normal surfaces for a triangulation with prescribed Euler characteristic(s).
	
	PARAMETERS
	Triangulation3 tri: the full triangulation
	List<int> eulerChars: list of desired Euler characteristics
	Boolean vertex: if True, generate vertex surfaces; if False, generate fundamental surfaces
	Boolean almNorm: if True, generate almost normal surfaces; if False, generate normal surfaces
	
	RETURNS
	List<NormalSurface> filteredSurfaces: the (almost) normal surfaces with prescribed Euler characteristic(s).
	"""
	#Generate surfaces according to vertex and almNorm booleans.
	if vertex:
		if almNorm:
			rawSurfaces = NormalSurfaces(tri, NS_AN_STANDARD, NS_VERTEX)
		else:
			rawSurfaces = NormalSurfaces(tri, NS_STANDARD, NS_VERTEX)
	else:
		if almNorm:
			rawSurfaces = NormalSurfaces(tri, NS_AN_STANDARD, NS_FUNDAMENTAL)
		else:
			rawSurfaces = NormalSurfaces(tri, NS_STANDARD, NS_FUNDAMENTAL)
	
	#Create a surface filter, restrict it to eulerChars
	filter = SurfaceFilterProperties()
	filter.setEulerChars(eulerChars)
	filteredSurfaces = NormalSurfaces(rawSurfaces, filter)
	#NTD: separate this into a surface generator and a filter, then save
	#	surfaces to each triangulation. Aaaaaaaaaaaaaa
	return filteredSurfaces

def checkSurfaces(tri, eulerChar, vertex, almNorm):
	"""This function generates normal surfaces with a given euler characteristic and cuts along each until a Heegaard splitting is found. This is currently only used in the initial Heegaard check, where cliques are not considered. NTD: should also just refer to checkOneSurface
	
	PARAMETERS
	Triangulation3 tri: the full triangulation
	int eulerChar: desired Euler characteristic
	Boolean vertex: if True, generate vertex surfaces; if False, generate fundamental surfaces
	Boolean almNorm: if True, generate almost normal surfaces; if False, generate normal surfaces
	
	RETURNS
	List<NormalSurface> filteredSurfaces: the (almost) normal surfaces with prescribed Euler characteristic(s).
	"""
	#Value to return if a splitting is found.
	success = False
	#Local value for genus
	g = (2 - eulerChar)/2
	
	#Generate potential surfaces
	filteredSurfaces = filteredSurfaceGenerator(tri, [eulerChar], vertex, almNorm)
	for surface in filteredSurfaces:
		if not surface.isConnected():
			#Surface isn't connected so will not split the manifold into handlebodies
			continue
		#Cut along the surface and simplify the two pieces
		cut = surface.cutAlong()
		cutTriangulated = cut.triangulateComponents()
		(cutTriangulated[0]).intelligentSimplify()
		handlebodyGenus0 = (cutTriangulated[0]).recogniseHandlebody()
		if handlebodyGenus0 == g:
			#These are checked separately to save on time spent simplifying
			(cutTriangulated[1]).intelligentSimplify()
			handlebodyGenus1 = (cutTriangulated[1]).recogniseHandlebody()
			if handlebodyGenus1 == g:
				success = True
				break
	return success

#This cuts along the given normal surface (first checking if it has the right euler characteristic and is orientable and connected) and checks if it is a splitting.
#aaaaaaaa NTD SHOULD MERGE WITH CHECKSURFACES()
def checkOneSurface(normSurf, eulerChar):
	"""This function cuts a triangulation along a given (almost) normal surface (first checking if it has the required Euler characteristic, and is orientable and connected), and checks if it is a splitting (both sides are handlebodies).
	
	PARAMETERS
	NormalSurface normSurf: the (almost) normal surface to cut along
	int eulerChar: desired Euler characteristic
	
	RETURNS
	Boolean success: True if normSurf was a Heegaard splitting, False otherwise.
	"""
	success = False
	#Local value for genus
	g = (2 - eulerChar)/2
	if normSurf.isConnected() and normSurf.eulerChar() == eulerChar and normSurf.isOrientable():
		cut = normSurf.cutAlong()
		cutTriangulated = cut.triangulateComponents()
		#Safety check (in case the surface is non-separating)
		if len(cutTriangulated) == 2:
			(cutTriangulated[0]).intelligentSimplify()
			if (cutTriangulated[0]).recogniseHandlebody() == g:
				(cutTriangulated[1]).intelligentSimplify()
				if (cutTriangulated[1]).recogniseHandlebody() == g:
					success = True
	return success

def surfaceGraphGenerator(tri, minEuler, maxEuler, vertex, almNorm):
	"""Generates a graph where surfaces are connected if they are locally compatible. Includes surfaces with Euler characterstic minEuler up to maxEuler. Adds double/triple/etc surfaces of those with nonzero Euler characteristic and sets their neighbours to be the neighbours of their original surface only (e.g. so 2*s and 3*s are not listed as compatible so that 5*s will not appear when testing 2*s+3*s).
	
	PARAMETERS
	Triangulation3 tri: triangulation
	int minEuler: minimum euler characteristic to be considered
	int maxEuler: maximum euler characteristic to be considered
	Boolean vertex: if True, generate vertex surfaces; if False, generate fundamental surfaces
	Boolean almNorm: if True, generate almost normal surfaces; if False, generate normal surfaces
	
	RETURNS
	List<NormalSurface> trueSurfaceList: list of surfaces that appear in the graph
	dict<int, List<int>> graph: if graph[i]=[a,b], then surface trueSurfaceList[i] is locally compatible with trueSurfaceList[a], and also trueSurfaceList[b].
	"""
	#Generate all surfaces with Euler characteristics in the prescribed range
	surfacesPossible = filteredSurfaceGenerator(tri, list(range(minEuler,maxEuler+1)), vertex, almNorm)
	trueSurfaceList = []
	graph = {}
	
	#For each surface, find all surfaces that are locally compatible with it, and insert into the graph.
	for i1 in range(surfacesPossible.size()):
		surf1 = surfacesPossible.surface(i1)
		trueSurfaceList.append(surf1)
		graph[i1]=[]
		for i2 in range(surfacesPossible.size()):
			surf2 = surfacesPossible.surface(i2)
			if surf1.locallyCompatible(surf2) and surf1!=surf2:
				graph[i1].append(i2)
	
	#If adding multiples of surf1 into the graph, as long as k*surf1 has Euler characteristic within the prescribed range. Local compatability is copied from graph[surf1] to save from interactions between (for example) 2*surf1 and 3*surf1.
	for i1 in range(surfacesPossible.size()):
		surf1 = surfacesPossible.surface(i1)
		#Work-around for Regina integer types:
		ecS1 = int((surf1.eulerChar()).stringValue())
		
		#We specifically avoid torii here by declaring that the Euler characteristic is nonzero. For a probabilisitcally appropriate treatment of torii, one could add muliples of each torus up to an appropriate number (such as 10) but ONLY for torii that are NOT thin edge links.
		#If the euler char is nonzero and (ecS1+ecS1) is still in our range:
		if ecS1 != 0 and minEuler <= ecS1*2 <= maxEuler:
			if ecS1<0:
				maxMult = minEuler//ecS1
			elif ecS1>0:
				maxMult = maxEuler//ecS1
			#Multiples of the surface from 2 up to (no. times ecs1 can go into ec)
			for i3 in range(2, 1 + maxMult):
				#Neighbours of i3*surf will be the same as neighbours of i3
				graph[len(graph)] = graph[i1]
				trueSurfaceList.append(surf1*i3)
	return trueSurfaceList, graph

def findWeightedCliques(graph, surfaceID, EC, gadgetFlag):
	"""Finds all cliques (groups of locally compatible surfaces) whose weights (euler characteristics) sum to EC, checking if they're splittings in the process. Breaks and returns True if successful at some point.
	
	PARAMETERS
	dict<int, List<int>> graph: according to surfaceGraphGenerator()
	List<NormalSurface> surfaceID: according to surfaceGraphGenerator()
	int EC: the prescribed weight of cliques to search for. This is the desired Euler characteristic of the sum of surfaces.
	Boolean gadgetFlag: if True, tubings from the gadget will decrease the Euler characteristic by 2, so we search for cliques with characteristic EC+2 instead. If False, business as usual.

	RETURNS
	Boolean success: if True, a splitting has been found from some clique.
	"""

	#Set of cliques found so far
	cliques = []
	
	#Setting n as the Euler characteristic, based on gadget flag.
	if gadgetFlag:
		n = EC + 2
	else:
		n = EC
	
	def cliqueWeight(clique, surfaceID):
		"""Getter method for summing up all weights in a given clique
	
		PARAMETERS
		List<int> clique: a list of integers, representing surfaces in the clique, numbered by surfaceID
		List<NormalSurface> surfaceID: a list of surfaces, where surfaceID[i] corresponds to the surface with integer ID i from clique

		RETURNS
		int: calculated weight of the clique
		"""
		return sum((surfaceID[node]).eulerChar() for node in clique)
	
	def bronKerboschWeightSum(R, P, X, weightSum):
		"""A version of the Bron-Kerbosch algorithm for finding maximal cliques - modified to target weights of cliques instead of restricting to maximal cliques.
	
		PARAMETERS
		set<int> R: represents surfaces considered in a particular clique so far.
		set<int> P: represents candidate surfaces which could be included in the particular clique.
		set<int> X: represents surfaces which cannot be included in the particular clique.
		int weightSum: calculate sum of the weights in the particular clique so far.

		RETURNS
		Boolean success: True if a valid clique Heegaard splitting has been
			found. Cancels the recursive clique search/generation.
		"""
		
		#e.g. we are searching with n=-4. If our weightSum is NOT in the range[-6,2], specifically below, then the only way for it to enter this range is by adding a surface to the current clique with Euler characteristic greater than 2 (so, a disconnected surface). If our weightSum is above this range, then it means that at some point, we added at least one surface with Euler characteristic greater than 2 (disconnected), or multiple with Euler characteristic 1 or 2 (not allowed, by Rubinstein's algorithm).
		#Alternatively, it's a dummy-check if the user has set up an invalid value of EC initially.
		#If we don't quit here, we may get stuck with a constantly increasing or decreasing Euler characteristic.
		if not(n-2 <= weightSum <= 2):
			return False
		
		#Default value:
		success = False
		
		#Valid weight sum, and R is not empty (because an empty set tautologically has Euler characteristic 0)
		if weightSum == n and R!=set():
			#Add our 'clique so far' R to the list of cliques.
			cliques.append(R)
			
			#We've found a potential surface. Now, to test it.
			
			#First, sum over the clique to get our surface
			cliqueList = list(R)
			surfSum = surfaceID[cliqueList[0]]
			for i in range(1, len(cliqueList)):
				surfSum += surfaceID[cliqueList[i]]
			
			#If we are using the gadget, construcy and test all tubings.
			if gadgetFlag:
				tri = surfSum.triangulation()
				allTubes = []
				#For each tetrahedron in tri, find all tubings.
				for tetNum in range(tri.countTetrahedra()):
					allTubesInTet = surfaceGadgetCombiner(tri, tetNum, surfSum)
					allTubes.extend(allTubesInTet)
				#Test each surface, breaking on success. Note we use EC as the Euler characterisitc instead of n. Regardless, tubedSurface.eulerChar()==EC anyway.
				for tubedSurface in allTubes:
					success = checkOneSurface(tubedSurface, EC)
					if success:
						break
						
			#Not using the gadget, test the surface normally.
			else:
				success = checkOneSurface(surfSum, EC)
			
			#We've found a valid splitting in the lines above! Return success to exit the recursive function.
			if success:
				return success
			
		#Either our previously tested clique was not a valid splitting, or our clique doesn't have the required weight yet. Test surfaces that are compatible with those in R.
		for node in P.copy():
			neighbours = set(graph[node])
			#Adding node (a common neighbour of R) to the clique, now considering neighbours of R that are also neighbours of node, adding the Euler characteristic of node to the current weightSum.
			bronSuccess=bronKerboschWeightSum(R.union([node]), P.intersection(neighbours), X.intersection(neighbours), weightSum + (surfaceID[node]).eulerChar())
			#If a successful splitting is ever found in the recursion, success=True will be returned. This is saved above as bronSuccess. This is another check to escape the recursion.
			if bronSuccess:
				return bronSuccess
				break
			#All cliques formed by including node have been considered, and are not splittings. Remove node from the candidates, and keep searching.
			P.remove(node)
			X.add(node)

	#Initial recursive call:
	success=bronKerboschWeightSum(set(), set(graph.keys()), set(), 0)
	#I believe if we get to this point it means we found no success and so it will return the default false value of it.
	return success

def genusNTest(tri, wantedG):
	"""Tests if tri has a Heegaard splitting of the specified genus. This tests vertex almost normal surfaces, then cliques of fundamental almost normal surfaces, then cliques of normal surfaces formed by tubings of fundamental normal surfaces.
	
	PARAMETERS
	Triangulation3 tri: the triangulation to determine Heegaard genus of. This should be zero-efficient AND it should be reasonably simplified.
	int wantedG: look for a Heegaard splitting of this genus.

	RETURNS
	Boolean attempt(1,2,3): True if a splitting was found (calculated using a particular method). If attempt1 is False, attempt2 is tested, etc.
	int wantedG: just spitting the genus back out.
	str: the method used. This is either "Vertex", "Cliques" or "Gadget with cliques"
	"""
	#Calculation for the Euler characteristic.
	wantedEc = 2 - 2*wantedG

	#This checks with vertex=True and almostnormal=True
	attempt1 = checkSurfaces(tri, wantedEc, True, True)
	if attempt1:
		return attempt1, wantedG, "Vertex"
	else:
		#Checks cliques of almost normal fundamental surfaces (so, vertex=False and almostnormal=True). Max Euler characteristic is set to -1 when generating surfaces so that tori are not used.
		#NTD: make surfacegraphgenerator use the existing set of surfaces
		surfaceID, graph = surfaceGraphGenerator(tri, wantedEc, -1, False, True)
		#Not using gadget so gadFlag=false
		attempt2 = findWeightedCliques(graph, surfaceID, wantedEc, False)
		if attempt2:
			return attempt2, wantedG, "Cliques"
		else:
			#Checks cliques of tubed normal fundamental surfaces created using the gadget (so, vertex=False and almostnormal=False)
			#I want max eulerchar 1 because tubing could make euler char go down, but I don't want to use 2 because vertex links do nothing.
			surfaceID, graph = surfaceGraphGenerator(tri, wantedEc, 1, False, False)
			#Using gadget so gadFlag=true
			attempt3 = findWeightedCliques(graph, surfaceID, wantedEc, True)
			if attempt3:
				return attempt3, wantedG, "Gadget and cliques"
			else:
				#The only time we return False.
				return attempt3, wantedG, "Fail"

def heegaardRunnerFunc(sig):
	"""This simplifies a triangulation, then checks for splittings of genus g=min(2, rank of homology), then for g+1, etc.
	
	PARAMETERS
	str sig: the isomorphism signature of the triangulation to determine Heegaard genus of.

	RETURNS
	str: a comma separated string including sig, the method used, the run time, the isomorphism signature of the simplified triangulation, the genus, the rank of the first homology group, and whether the triangulation is of a Haken manifold.
	"""
	startT = time.time()
	tri = Triangulation3(sig)
	tri.simplifyExhaustive(2)
	isHaken = tri.isHaken()
	hom1 = tri.homology().countInvariantFactors()
	minGenus = max(hom1, 2)
	success = False
	iter = 0
	while not success:
		(testSuccess, foundGenus, message) = genusNTest(tri, minGenus + iter)
		success = testSuccess
		iter += 1
	endT=time.time()
	print(sig+","+message+","+str(endT-startT)+","+tri.isoSig()+","+str(foundGenus)+","+str(hom1)+","+str(isHaken))
