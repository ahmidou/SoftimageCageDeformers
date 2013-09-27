#pragma once

#include <algorithm>
#include <xsi_application.h>

#include <xsi_math.h>
#include <math.h>

#include <xsi_primitive.h>
#include <xsi_geometry.h>
#include <xsi_kinematics.h>
#include <xsi_point.h>
#include <xsi_string.h>
#include <xsi_vector3.h>
#include <xsi_polygonface.h>
#include <vector>

using namespace std;
using namespace XSI::MATH;
using namespace XSI;

#define TOLERANCE 0.000000001
#define PI 3.1415926535897932

static Application app;

typedef enum
{
	ALL_VERTICES_AFFINITY,
	BOUNDARY_BARYCENTER_AFFINITY,
	BOUNDARY_AFFINITY,
}AFFINITY_TYPE;

typedef enum
{
	BOUNDARY_DISTANCE,
	CENTER_DISTANCE,
	HYBRID_DISTANCE,
	BASE_TRANSFORMATION_DISTANCE,
}DISTANCE_FUNCTION;

typedef enum
{
	MVC_ENGINE,
	GC_ENGINE,
	HC_ENGINE
}ENGINE_TYPE;



typedef enum
{
	VERTEX_CAGE,
	VERTEX_MESH
}VERTEX_TYPE;


typedef std::vector<CVector3> VECTOR3_VECTOR;
typedef std::vector<int> INTEGER_VECTOR;
typedef std::vector<float> FLOAT_VECTOR;
typedef std::vector<double> DOUBLE_VECTOR;
typedef std::vector<bool> BOOL_VECTOR;

typedef std::vector<INTEGER_VECTOR> INTEGER_MATRIX;


struct DEFORMATION_COORDINATES
{
	DOUBLE_VECTOR vertexCoordinates;
	DOUBLE_VECTOR faceCoordinates;
};

typedef std::vector<DEFORMATION_COORDINATES> DEFORMATION_COORDINATES_VECTOR;

inline CVector3 faceNormalGeneric( CVector3 a, CVector3 b, CVector3 c)
{
	CVector3 v1= b.SubInPlace(a);
	CVector3 v2=c.SubInPlace(a);

	return v1.Cross(v1, v2 );
}

