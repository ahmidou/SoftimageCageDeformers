/******************************************************************/
/*                  MeanValueCoordinates Class CPP                    */
/******************************************************************/

#include "PMVC.h"

/*********************************************************/
/*   MeanValueCoordinates::MeanValueCoordinates(  )              */
/*Constructor per defecte de la classe MeanValueCoordinates  */
/*********************************************************/

PMVCoordinates::PMVCoordinates( )
{
}

/********************************************/
/*   ~MeanValueCoordinates::MeanValueCoordinates()  */
/*Destructor de la classe MeanValueCoordinates  */
/********************************************/

PMVCoordinates::~PMVCoordinates()
{
}


/*****************   FUNCION *****************/

//vector<float> PMVCoordinates::computeCoordinates ( Point v, Geometry &cGeo, vector<vector<double> > points, vector<vector<double> > samples )
//{
//	//cout<<"MVC Compute Coordinates"<<endl;
//	//VERTEX_VECTOR cageVertices = cage->getVertices();
//	//FACE_VECTOR cageFaces = cage->getFaces();
//	CVector3 pos = v.GetPosition();
//	double test[] = {0,0,1,1,2,2};
//	PointLocatorData closestPointLocator = cGeo.GetRaycastIntersections( 1, test, test, siSegmentIntersection );
//
//	//cout<<"MVC Computation --> "<<timeDef<<endl;
//	return coordinates;
//}

CVector3Array PMVCoordinates::pointsOnSphere(int N)
{
//	vector<vector<double> > samples(3,vector<double> ( N, 0 ));
	N = float(N);
	CVector3Array samples(N);
	double inc = PI * (3.0 - sqrt(5.0)) ;
	double off = 2.0 / N ;
	#pragma omp parallel for
	for (int k=0; k< int(N); k++)
	{
		double y = k * off - 1 + (off / 2);
		double r = sqrt(1 - y*y);
		double phi = k * inc ;
		CVector3 v(cos(phi)*r, y, sin(phi)*r);
		samples[k] = v;
	}
	return samples;
}

double PMVCoordinates::sum(const std::vector<double>& v)
{
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for(size_t ii=0; ii< v.size(); ++ii)
	{
		sum += v[ii];
	}
	return sum;
}
