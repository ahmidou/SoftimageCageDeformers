#pragma once

#include <stdint.h>
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

/* ----------------------------------------------------------------------- */
/*
  Easy embeddable cross-platform high resolution timer function. For each 
  platform we select the high resolution timer. You can call the 'ns()' 
  function in your file after embedding this. 
*/

#if defined(__linux)
#  define HAVE_POSIX_TIMER
#  include <time.h>
#  ifdef CLOCK_MONOTONIC
#     define CLOCKID CLOCK_MONOTONIC
#  else
#     define CLOCKID CLOCK_REALTIME
#  endif
#elif defined(__APPLE__)
#  define HAVE_MACH_TIMER
#  include <mach/mach_time.h>
#elif defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#endif
static uint64_t ns() {
  static uint64_t is_init = 0;
#if defined(__APPLE__)
    static mach_timebase_info_data_t info;
    if (0 == is_init) {
      mach_timebase_info(&info);
      is_init = 1;
    }
    uint64_t now;
    now = mach_absolute_time();
    now *= info.numer;
    now /= info.denom;
    return now;
#elif defined(__linux)
    static struct timespec linux_rate;
    if (0 == is_init) {
      clock_getres(CLOCKID, &linux_rate);
      is_init = 1;
    }
    uint64_t now;
    struct timespec spec;
    clock_gettime(CLOCKID, &spec);
    now = spec.tv_sec * 1.0e9 + spec.tv_nsec;
    return now;
#elif defined(_WIN32)
    static LARGE_INTEGER win_frequency;
    if (0 == is_init) {
      QueryPerformanceFrequency(&win_frequency);
      is_init = 1;
    }
    LARGE_INTEGER now;
    QueryPerformanceCounter(&now);
    return (uint64_t) ((1e9 * now.QuadPart)  / win_frequency.QuadPart);
#endif
}
/* ----------------------------------------------------------------------- */

