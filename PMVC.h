/******************************************************************/
/*                   MeanValueCoordinates Class H                     */
/******************************************************************/

#pragma once

#include "definitions.h"

#define EPSILON_MVC 0.00001

using namespace std;

class PMVCoordinates
{
	public:
		PMVCoordinates();
		~PMVCoordinates();

		double sum(const vector<double>& v);
		CVector3Array pointsOnSphere(int N);
		vector<float> computeCoordinates ( Point v, Geometry &cGeo, vector<vector<double> > points, vector<vector<double> > samples );
		CVector3 computeDeformedVertexNode ( Point  v, Geometry  c );

} ;
