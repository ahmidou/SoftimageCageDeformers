/******************************************************************/
/*                   GreenCoordinates Class H                     */
/******************************************************************/

#pragma once


#include "definitions.h"

#define EPSILON_GC 0.0001

using namespace std;
using namespace XSI::MATH;
using namespace XSI;


class GreenCoordinates
{
	private:
		DEFORMATION_COORDINATES coordinates;
		CVector3Array cageVerticesDeformed;
		CVector3Array cageVertices;

	protected:
		DEFORMATION_COORDINATES computeCoordinates ( Point v, CPointRefArray cageVertices, CFacetRefArray cageFaces );

		double GCTriInt( CVector3 p, CVector3 v1, CVector3 v2, CVector3 n );
		double GCTriInt2(CVector3 p, CVector3 v1, CVector3 v2, CVector3 e);
		double integralEvaluation( double theta, double c, double lambda );
		
	public:
		GreenCoordinates();
		~GreenCoordinates();

		void PrintVector( const CString & in_desc, const CVector3& in_v );
		DEFORMATION_COORDINATES computeCoordinates(  Point v, Geometry c );
		CVector3 computeDeformedVertexNode (DEFORMATION_COORDINATES* gc, Geometry * c, Geometry * dc );
} ;
