/******************************************************************/
/*                  GreenCoordinates Class CPP                    */
/******************************************************************/

#include "GreenCoordinates.h"

/*********************************************************/
/*   GreenCoordinates::GreenCoordinates(  )              */
/*Constructor per defecte de la classe GreenCoordinates  */
/*********************************************************/

GreenCoordinates::GreenCoordinates( ) 
{
}

/********************************************/
/*   ~GreenCoordinates::GreenCoordinates()  */
/*Destructor de la classe GreenCoordinates  */
/********************************************/

GreenCoordinates::~GreenCoordinates() 
{
}

/*************************************************************************************************/
/*   DeformationCoordiantes GreenCoordinates::computeCoordinates ( Point n, double epsilon )      */
/*Metode que retorna les Green Coordinates del punt indicat per parametres                       */
/*************************************************************************************************/

void GreenCoordinates::PrintVector( const CString & in_desc, const CVector3& in_v )
{
	CString x = CValue(in_v.GetX()).GetAsText() ;
	CString y = CValue(in_v.GetY()).GetAsText() ;
	CString z = CValue(in_v.GetZ()).GetAsText() ;

	Application().LogMessage( in_desc + L" " + x + L"," + y + L"," + z ) ;
}

DEFORMATION_COORDINATES GreenCoordinates::computeCoordinates ( Point v, Geometry c )
{
	//cout<<"GC Compute Coordinates"<<endl;
	DEFORMATION_COORDINATES coordinates = computeCoordinates( v, c.GetPoints(), c.GetFacets() );
	return coordinates;
}

DEFORMATION_COORDINATES GreenCoordinates::computeCoordinates ( Point  v, CPointRefArray cageVertices, CFacetRefArray cageFaces )
{
	DEFORMATION_COORDINATES gc;
	//VERTEX_VECTOR cageVertices = c->getVertices();
	//FACE_VECTOR cageFaces = c->getFaces();
	CVector3 n = v.GetPosition();

	coordinates.vertexCoordinates.clear();
	coordinates.faceCoordinates.clear();

	coordinates.vertexCoordinates.reserve(cageVertices.GetCount());
	for ( int i=0;i<cageVertices.GetCount();i++ )
	{
		coordinates.vertexCoordinates.push_back ( 0.0 );
	}
	coordinates.faceCoordinates.clear();

	coordinates.faceCoordinates.reserve(cageFaces.GetCount());
	for ( int t=0;t<cageFaces.GetCount();t++ )
	{
		//Get vertices from the current cage face
		CPointRefArray pnts = Facet(cageFaces[t]).GetPoints();

		CLongArray ids;
		ids.Add(Point(pnts[0]).GetIndex());
		ids.Add(Point(pnts[1]).GetIndex());
		ids.Add(Point(pnts[2]).GetIndex());

		CVector3Array v;
		v.Add ( Point(pnts[0]).GetPosition());
		v.Add ( Point(pnts[1]).GetPosition());
		v.Add ( Point(pnts[2]).GetPosition());

		CVector3 normalT = faceNormalGeneric( v[0], v[1], v[2] );
		normalT.NormalizeInPlace();

		//Translate n to the origin and consequently, the vertices of the actual face
		for ( int l=0;l<3;l++ )
		{
			v[l] = v[l].SubInPlace(n);
		}

		//Project the point n onto the plane defined by the actual face
		CVector3 vect;
		vect.Scale(normalT.Dot( CVector3(v[0]) ),normalT);
		CVector3 p = vect ;
		CVector3 v0, v1;
		INTEGER_VECTOR s;
		DOUBLE_VECTOR I;
		DOUBLE_VECTOR II;
		CVector3Array N_;

		s.reserve(3);
		I.reserve(3);
		II.reserve(3);
		for ( int l=0;l<3;l++ )
		{
			v0 = v[l];
			if ( l==2 )
				v1 = v[0];
			else
				v1 = v[l+1];

			CVector3 pv0; pv0.Sub(v0 , p) ;
			CVector3 pv1; pv1.Sub(v1 , p) ;
			CVector3 prod = prod.Cross( pv0, pv1 );
			if ( prod.Dot ( normalT ) <0. )
				s.push_back ( -1 );
			else
				s.push_back ( 1 );

			I.push_back ( GCTriInt ( p, v0, v1, CVector3 ( 0,0,0 ) ) );
			II.push_back ( GCTriInt ( CVector3 ( 0,0,0 ), v1, v0, CVector3 ( 0,0,0 ) ) );
			CVector3 v0_ = v0;
			CVector3 v1_ = v1;
			CVector3 q = q.Cross (v1_, v0_ );
			q.NormalizeInPlace();
			N_.Add ( q );
//			PrintVector("",CVector3(q));

		}
		double Itotal = 0.0;
		CVector3 IItotal ( 0.0,0.0,0.0 );
		for ( int k = 0;k<3;k++ )
		{
			CVector3 tmp= tmp.Scale(double(II[k]),N_[k]);
			Itotal +=  s[k] * I[k] ;
			IItotal.AddInPlace( tmp );
//			app.LogMessage(CString(II[k]));
//			PrintVector("",N_[k]);
		}

		Itotal = -fabs ( Itotal );
		coordinates.faceCoordinates.push_back ( -Itotal );
		CVector3 w = w.Scale(Itotal,normalT).AddInPlace(IItotal);
//		PrintVector("",w);
//		for ( int l = 0;l<N_.GetCount();l++ )
//		{
//			PrintVector("",CVector3(N_[l]));
//		}
		if ( w.GetLength() > EPSILON_GC )
		{
			int nextV;
			for ( int l = 0;l<3;l++ )
			{
				if ( l==2 )
					nextV = 0;
				else
					nextV = l+1;

				CVector3 v_ = v[l] ;

				//Comprovem la distancia entre el vertex de la caixa i el vertex del model per tal de decidir si el descartem o no
// 				Vector3D vn( v[l], n );
// 				cout<<"Distancia: "<<vn.calcula_modul()<<endl;
// 				if( vn.calcula_modul() < 1.5 )
// 				{
//	 				cout<<"Suma: "<<( N_[nextV].producteEscalar(w) / N_[nextV].producteEscalar(v_))<<endl;
				coordinates.vertexCoordinates[ids[l]] += ( N_[nextV].Dot( w ) / N_[nextV].Dot( v_ ) );

//				app.LogMessage("GC Vertex: "+ CString(nextV));
//	 				cout<<"GC Vertex: "<<gc.gc_cage_vertexs[ids[l]]<<endl;
// 				}
			}
		}
	}

	//Comprovaci� de que el sumatori de coordenades de vertexs dona 1
	/*double suma = 0;
	for ( int i=0;i<cageVertices.size();i++ )
	{
		suma = suma + gc.vertexCoordinates[i];
	}*/
	//n->print();
	//cout<<"* Comprovaci� sumatori GC Vertexs: "<<suma<<endl;
	
	return coordinates;
}

/*****************************************************************************************/
/*   double GreenCoordinates::GCTriInt ( Vector3 p, Vector3 v1, Vector3 v2, Vector3 n )              */
/*Metode que calcula la integral del triangles format pels punts indicats per parametres */
/*****************************************************************************************/

double GreenCoordinates::GCTriInt ( CVector3 p, CVector3 v1, CVector3 v2, CVector3 n )
{
	CVector3 v1_v2 = v1_v2.Sub(v2,v1);
	CVector3 v1_p = v1_p.Sub(p ,v1);
	CVector3 p_v1 = p_v1.Sub(v1 , p);
	CVector3 p_v2 = p_v2.Sub(v2 , p);
	CVector3 np = np.Sub(p , n);

	double alfa = acos ( v1_v2.Dot ( v1_p ) / ( v1_v2.GetLength() *v1_p.GetLength() ) );
	double beta = acos ( p_v1.Dot ( p_v2 ) / ( p_v1.GetLength() *p_v2.GetLength() ) );
	double lambda = ( pow ( v1_p.GetLength(),2 ) ) * ( pow ( sin ( alfa ),2 ) );
	double c = pow ( np.GetLength(),2 );

	double integral;
	if ( fabs ( alfa-PI ) >0.0001 && fabs ( ( alfa+beta )-PI ) >0.0001 ) //if( alfa!=PI && (alfa+beta)!=PI )
	{
		double intPart1 = 1., intPart2 = 1.;
		intPart1 = integralEvaluation ( PI - alfa, c, lambda );
		intPart2 = integralEvaluation ( PI - alfa - beta, c, lambda );
		
		integral = ( -1./ ( 4.*PI ) ) * fabs ( intPart1 - intPart2 - ( sqrt ( c ) *beta ) );
	}
	else
	{
		integral = 0.0;
	}

	return integral;
}

//double GreenCoordinates::GCTriInt2(CVector3 p, CVector3 v1, CVector3 v2, CVector3 e)
//{
//	CVector3 v1_v2 = v1_v2.Sub(v2,v1).NormalizeInPlace();
//	CVector3 v1_p = v1_p.Sub(p ,v1);
//	CVector3 v1_pN = v1_pN.Normalize(v1_p);
//	CVector3 p_v1 = p_v1.Sub(v1 , p).NormalizeInPlace();
//	CVector3 p_v2 = p_v2.Sub(v2 , p).NormalizeInPlace();
//	CVector3 ep = ep.Sub(p , e);
//
//	const double alpha = acos( RANGED(-1.0, v1_v2.Dot ( v1_pN ), 1.0) );
//	const double beta  = acos( RANGED(-1.0, p_v1.Dot ( p_v2 ), 1.0) );
//	const double lambda = v1_p.GetLengthSquared() * sin(alpha) * sin(alpha);
//	const double c      = ep.GetLengthSquared();
//	const double theta[2] = { M_PI - alpha, M_PI - alpha - beta };
//
//	double I[2] = {0,0};
//
//	for (int i = 0; i < 2; ++i)
//	{
//			double S = sin(theta[i]);
//			double C = cos(theta[i]);
//
//			double sign = S < 0 ? -1.0 : 0 < S ? 1.0 : 0.0;
//
//			if (sign == 0.0)
//					I[i] = 0.0;
//			else
//			{
//					double M = (-sign / 2.0);
//					double N = 2 * sqrt(c) * atan((sqrt(c) * C) / sqrt(lambda + (S * S * c)));
//					double O = sqrt(lambda);
//					double P = (2 * sqrt(lambda) * S * S) / pow(1.0 - C, 2);
//					double denom = ( (c*(1+C) + lambda + sqrt((lambda * lambda) + (lambda * c * S * S)) ));
//					double Q = (2 * c * C) / denom;
//					double R = 1.0 - Q;
//
//					I[i] = M * (N + (O * log(P * R)));
//			}
//	}
//
//double ret = (-0.25 / M_PI) * fabs(I[0] - I[1] - sqrt(c) * beta);
//
//	return ret;
//}

/****************************************************************************/
/*   double GreenCoordinates::GCTriInt ( Point p, Point v1, Point v2, Point n ) */
/*Metode que avalua la integral a traves del parametres indicats            */
/****************************************************************************/

double GreenCoordinates::integralEvaluation ( double theta, double c, double lambda )
{
	double S = sin ( theta );
	double C = cos ( theta );

	double sub1 = -1./2.;
	if ( S<0. )
		sub1 = 1./2.;

	double sub2 = 2.*sqrt ( c ) *atan ( ( sqrt ( c ) *C ) / ( sqrt ( lambda + ( ( pow ( S,2 ) ) *c ) ) ) );
	double sub3 = 0. ;
	if( lambda != 0. )
		sub3 = sqrt ( lambda ) *log ( ( ( 2.*sqrt ( lambda ) * ( pow ( S,2 ) ) ) / ( pow ( ( 1.-C ),2 ) ) ) * ( 1. - ( ( 2.*c*C ) / ( ( c* ( 1.+C ) ) +lambda+sqrt ( ( pow ( lambda,2 ) ) + ( lambda*c* ( pow ( S,2 ) ) ) ) ) ) ) );

	double result = ( sub1 * ( sub2 + sub3 ) );
	return result;
}

CVector3 GreenCoordinates::computeDeformedVertexNode ( DEFORMATION_COORDINATES* gc, Geometry * c , Geometry * dc)
{
	//cout<<"GC Compute Deformed Vertex Node"<<endl;
	//cout<<"GREEN COORDINATES!!!"<<endl;
	CVector3 sum_vertexs ( 0.,0.,0. );
	CVector3 sum_faces ( 0.,0.,0. );

	cageVerticesDeformed = dc->GetPoints().GetPositionArray();
	cageVertices = c->GetPoints().GetPositionArray();
	CFacetRefArray cageFaces = c->GetFacets();


	for ( int i=0;i<cageVertices.GetCount();i++ )
	{
		CVector3 n = cageVerticesDeformed[i];
		sum_vertexs.AddInPlace( n.ScaleInPlace(gc->vertexCoordinates[i]) );
	}

	double scaling;
	int v0, v1, v2;
	CVector3 normalT, normalTD;
	for ( int j=0;j<cageFaces.GetCount();j++ )
	{
		CPointRefArray pnts = Facet(cageFaces[j]).GetPoints();
		v0 = Point(pnts[0]).GetIndex();
		v1 = Point(pnts[1]).GetIndex();
		v2 = Point(pnts[2]).GetIndex();
		
		normalT = faceNormalGeneric( cageVertices[v0], cageVertices[1], cageVertices[2] );
		normalTD = faceNormalGeneric( cageVerticesDeformed[v0], cageVerticesDeformed[1], cageVerticesDeformed[2] );

		if( normalT.GetLength()!= 0.)
			scaling = normalTD.GetLength() /normalT.GetLength();
		else
			scaling = 0.;
		normalTD.NormalizeInPlace();
		sum_faces.AddInPlace( normalTD.ScaleInPlace( gc->faceCoordinates[j]*scaling ) );
	}

	//Vector3 res = ( sum_vertexs + sum_faces );
	//if( _isnan( res.x) ||_isnan( res.y) ||_isnan( res.z) )
	//	cout<<"NAN!!!"<<this<<endl;
	//cout<<"Punt deformat: "<<res<<endl;
	CVector3 result = result.Add(sum_vertexs , sum_faces );
	return result;
}

