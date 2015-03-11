/****************************************************************************
 *  Copyright (C) 2012 Animal Logic Pty Ltd                                 *
 *                                                                          *
 *  These coded instructions, statements, and computer programs contain     *
 *  unpublished  proprietary  information of Animal Logic Pty Ltd,          *
 *  and are protected by copyright law.  They  may  not be disclosed        *
 *  to  third  parties  or copied or duplicated in any form, in whole or    *
 *  in part, without the prior written consent of  Animal Logic Pty Ltd.    *
 *                                                                          *
 ****************************************************************************/
#ifndef linux
#define WIN32_LEAN_AND_MEAN
	#include <windows.h> // Needed for OpenGL on windows
	typedef  UINT  uint;
#endif

#include <string.h>
#include <algorithm>
#include <xsi_application.h>
#include <xsi_context.h>

#include <xsi_math.h>
#include <math.h>
#include <xsi_argument.h>
#include <xsi_menu.h>
#include <xsi_menuitem.h>

#include <xsi_comapihandler.h>
#include <xsi_selection.h>
#include <xsi_x3dobject.h>
#include <xsi_primitive.h>
#include <xsi_polygonmesh.h>
#include <xsi_point.h>
#include <xsi_iceattribute.h>
#include <xsi_iceattributedataarray.h>
#include <xsi_iceattributedataarray2D.h>

#include "GC.h"
#include "PMVC.h"
#include <xsi_pluginregistrar.h>
#include <xsi_command.h>



XSIPLUGINCALLBACK CStatus XSILoadPlugin( PluginRegistrar& in_reg )
{
	in_reg.PutAuthor(L"ahmidoul");
	in_reg.PutName(L"Cages");
	in_reg.PutVersion(1,0);
//	in_reg.RegisterOperator(L"AL_mergeCurveOp");
	in_reg.RegisterCommand(L"CreateGreenCoordinates",L"CreateGreenCoordinates");
	in_reg.RegisterCommand(L"CreatePMVCoordinates",L"CreatePMVCoordinates");
	//RegistrationInsertionPoint - do not remove this line

	return CStatus::OK;
}

XSIPLUGINCALLBACK CStatus XSIUnloadPlugin( const PluginRegistrar& in_reg )
{
	CString strPluginName;
	strPluginName = in_reg.GetName();
	Application().LogMessage(strPluginName + L" has been unloaded.",siVerboseMsg);
	return CStatus::OK;
}

XSIPLUGINCALLBACK CStatus CreateGreenCoordinates_Init( CRef& in_ctxt )
{
	Context ctxt( in_ctxt );
	Command oCmd;
	oCmd = ctxt.GetSource();
	oCmd.PutDescription(L"Build the green coordinates");
	oCmd.SetFlag(siNoLogging,false);
	ArgumentArray oArgs = oCmd.GetArguments();
	oArgs.AddWithHandler( L"inputs","Collection" ,"");
	return CStatus::OK;
}

XSIPLUGINCALLBACK CStatus CreateGreenCoordinates_Execute( CRef& in_ctxt )
{
	Application().LogMessage(L"CreateGreenCoordinates_Execute called",siVerboseMsg);
	Context ctxt( in_ctxt );

	CValueArray args = ctxt.GetAttribute(L"Arguments");
	CValue inObj = args[0];
	CValueArray objArray = inObj;

	Geometry dGeo = X3DObject(objArray[0]).GetActivePrimitive().GetGeometry();
	CPointRefArray dPoints = dGeo.GetPoints();

	Geometry cage = X3DObject(objArray[1]).GetActivePrimitive().GetGeometry();

	CValueArray gcArray;
	CValueArray vArray;
	CValueArray fArray;
	GreenCoordinates gc;
//	#pragma omp parallel for
	for (LONG i=0; i<dPoints.GetCount(); i++ )
	{
		CValueArray v;
		CValueArray f;

		Point pnt(dPoints[i]);
		DEFORMATION_COORDINATES dc = gc.computeCoordinates ( pnt, cage);

		for (uint j=0; j<dc.vertexCoordinates.size(); j++ )
		{
			v.Add(dc.vertexCoordinates[j]);
		}
		for (uint j=0; j<dc.faceCoordinates.size(); j++ )
		{
			f.Add(dc.faceCoordinates[j]);
		}
		vArray.Add(v);
		fArray.Add(f);
	}
	gcArray.Add(vArray);
	gcArray.Add(fArray);

	ctxt.PutAttribute( L"ReturnValue", CValue(gcArray) );
	return CStatus::OK;
}


XSIPLUGINCALLBACK CStatus CreatePMVCoordinates_Init( CRef& in_ctxt )
{
	Context ctxt( in_ctxt );
	Command oCmd;
	oCmd = ctxt.GetSource();
	oCmd.PutDescription(L"Build the green coordinates");
	oCmd.SetFlag(siNoLogging,false);
	ArgumentArray oArgs = oCmd.GetArguments();
	oArgs.AddWithHandler( L"inputs","Collection" ,"");
	return CStatus::OK;
}




XSIPLUGINCALLBACK CStatus CreatePMVCoordinates_Execute( CRef& in_ctxt )
{
	Application().LogMessage(L"CreateGreenCoordinates_Execute called",siVerboseMsg);
	Context ctxt( in_ctxt );

	CValueArray args = ctxt.GetAttribute(L"Arguments");
	CValue inObj = args[0];
	CValueArray objArray = inObj;

	Geometry dGeo = X3DObject(objArray[0]).GetActivePrimitive().GetGeometry();
	CPointRefArray dPoints = dGeo.GetPoints();
	int nbSamples=64;
	PolygonMesh cage;
	cage = X3DObject(objArray[1]).GetActivePrimitive().GetGeometry();
	cage.SetupPointLocatorQueries(siClosestSurfaceRaycastIntersection,0,-1,0,nbSamples);

	PMVCoordinates pmv;
	CPointRefArray cPoints = cage.GetPoints();
	CValueArray v(dPoints.GetCount());
	CVector3Array samples =pmv.pointsOnSphere(nbSamples);

	double area = PI*4/nbSamples;

	//use a vector instead of a CVector3Array
	vector<vector<double> > points(3,vector<double>(dPoints.GetCount()));
	for (LONG i=0; i<dPoints.GetCount(); i++ )
	{
		CVector3 pos = Point(dPoints[i]).GetPosition();
		points[0][i]=pos[0];
		points[1][i]=pos[1];
		points[2][i]=pos[2];
	}

  vector<vector<float>> pmvArray2D;
  pmvArray2D.resize(dPoints.GetCount());
  vector<vector<LONG>> pmvID2D;
  pmvID2D.resize(dPoints.GetCount());
  vector<ULONG> subArraySizes;
  subArraySizes.resize(dPoints.GetCount());

  uint64_t start = ns();

	#pragma omp parallel for
	for (LONG i=0; i<dPoints.GetCount(); i++ )
	{	
		vector<float> pmvArray (cPoints.GetCount());
		vector<LONG> pmvID (cPoints.GetCount());
		vector<double> w(cPoints.GetCount(),0);
		CVector3 pntPos = Point(dPoints[i]).GetPosition();
		vector<double > posArray(nbSamples*3);
		for (LONG j=0; j<nbSamples; j++ )
		{
			posArray[j*3] = points[0][i];
			posArray[j*3+1] = points[1][i];
			posArray[j*3+2] = points[2][i];
		}
		PointLocatorData loc;
		#pragma omp critical
		{
			loc= cage.GetRaycastIntersections( nbSamples, (double*)&posArray[0], (double*)&samples[0], siSemiLineIntersection );
		}
		posArray.clear();
		const LONG lCount = loc.GetCount()*3;
        vector<double> pos(lCount);
        vector<LONG> triVtx(lCount);
        vector<float> triWei(lCount);
        cage.EvaluatePositions(loc, -1, 0, &pos[0]);
        cage.GetTriangleVertexIndexArray(loc, -1, 0, &triVtx[0]);
        cage.GetTriangleWeightArray(loc, -1, 0, &triWei[0]);
    	for (int s=0; s<nbSamples; s++ )
		{
			if(triVtx[s*3] != -1)
			{
				CVector3 dist(pos[s*3]-pntPos[0], pos[s*3+1]-pntPos[1], pos[s*3+2]-pntPos[2]) ;
    			double distScale = 1.0/dist.GetLength()*area;
    			w[triVtx[s*3]]+= triWei[s*3]*distScale;
    			w[triVtx[s*3+1]]+= triWei[s*3+1]*distScale;
    			w[triVtx[s*3+2]]+= triWei[s*3+2]*distScale;
			}
		}
		double tW = pmv.sum(w);
		uint k=0;
    for (uint j=0; j<w.size(); j++ )
    {
			if (w[j]/tW > 0.000001)
			{
				pmvArray[k] = w[j]/tW;
				pmvID[k] = LONG(j);
				k++;
			}
    }
		//pmvArray.Resize(0); //test
    pmvArray.resize(k);
		//pmvID.Resize(0);
		pmvID.resize(k);

    pmvArray2D[i] = pmvArray;
    pmvID2D[i] = pmvID;
    subArraySizes[i] = k;

    w.clear();
		pmvID.clear();
		pmvArray.clear();
	}
  uint64_t end = ns();
  app.LogMessage("coordinates generation: " + CString(1.0f/(end-start)) + " seconds");

  start = ns();

  // Set the ICE attributes
  CRefArray attr = dGeo.GetICEAttributes();
  ICEAttribute PMVWeight;
  ICEAttribute PMVID;
  if (attr.GetItem("PMVWeight").IsValid())
    PMVWeight = (dGeo.GetICEAttributeFromName( L"PMVWeight" ));
  else
	  PMVWeight = (CValue(dGeo.AddICEAttribute("PMVWeight",siICENodeDataFloat ,siICENodeStructureArray,siICENodeContextComponent0D)));

  if (attr.GetItem("PMVWeight").IsValid())
    PMVID = (dGeo.GetICEAttributeFromName( L"PMVID" ));
  else
	  PMVID = (CValue(dGeo.AddICEAttribute("PMVID",siICENodeDataLong ,siICENodeStructureArray,siICENodeContextComponent0D)));

  //ICEAttribute PMVWeight(CValue(dGeo.AddICEAttribute("PMVWeight",siICENodeDataFloat ,siICENodeStructureArray,siICENodeContextComponent0D)));
	//ICEAttribute PMVID(CValue(dGeo.AddICEAttribute("PMVID",siICENodeDataLong ,siICENodeStructureArray,siICENodeContextComponent0D)));

	CICEAttributeDataArray2D< float > w2D;
	PMVWeight.GetDataArray2D( w2D );
	CICEAttributeDataArray2D< LONG > id2D;
	PMVID.GetDataArray2D( id2D );
  
	for (LONG i=0; i<dPoints.GetCount(); i++ )
	{	

    CICEAttributeDataArray< float > wArray;
    CICEAttributeDataArray< LONG > idArray;

    w2D.GetSubArray( i, wArray );
    id2D.GetSubArray( i, idArray );

		wArray.SetArray(&pmvArray2D[i][0], pmvArray2D[i].size());
		idArray.SetArray(&pmvID2D[i][0], pmvID2D[i].size());
  }
  //w2D.SetArray2D((const float**)&pmvArray2D[0], pmvArray2D.size(), &subArraySizes[0]);
  //id2D.SetArray2D((const LONG**)&pmvID2D[0], pmvID2D.size(), &subArraySizes[0]);

  end = ns();
  app.LogMessage("ICE attr tranfer: " + CString(1.0f/(end-start)) + " seconds");

	ctxt.PutAttribute( L"ReturnValue", CValue("") );
  pmvID2D.clear();
	pmvArray2D.clear();
  subArraySizes.clear();
	samples.Clear();
	points.clear();
	
	return CStatus::OK;
}

//XSIPLUGINCALLBACK CStatus AL_mergeCurveOp_Define( CRef& in_ctxt )
//{
//	Context ctxt( in_ctxt );
//	CustomOperator oCustomOperator;
//	oCustomOperator.PutAlwaysEvaluate(false);
//	oCustomOperator.PutDebug(0);
//	return CStatus::OK;
//}
//
//XSIPLUGINCALLBACK CStatus AL_mergeCurveOp_Init( CRef& in_ctxt )
//{
//	OperatorContext ctxt( in_ctxt );
//	Application().LogMessage(L"mergeCurveOp_Init called",siVerboseMsg);
//	return CStatus::OK;
//}
//
//XSIPLUGINCALLBACK CStatus AL_mergeCurveOp_Term( CRef& in_ctxt )
//{
//	OperatorContext ctxt( in_ctxt );
//	Application().LogMessage(L"mergeCurveOp_Term called",siVerboseMsg);
//	return CStatus::OK;
//}
//
//XSIPLUGINCALLBACK CStatus AL_mergeCurveOp_Update( CRef& in_ctxt )
//{
//	OperatorContext ctxt( in_ctxt );
//	Operator op = ctxt.GetSource();
//	CRefArray iPorts = op.GetInputPorts();
//	Primitive out = ctxt.GetOutputTarget() ;
////	CTransformation srt = out.GetTransform();
//	CNurbsCurveDataArray outData;
//	double nb = iPorts.GetCount();
//	for (LONG i=0; i < nb; i++)
//	{
//	 	CNurbsCurveDataArray inData;
//		Primitive prim = (CRef)ctxt.GetInputValue(i);
//		NurbsCurveList curveList(prim.GetGeometry() );
//		curveList.Get(siSINurbs, inData);
//		for (LONG j=0; j < inData.GetCount(); j++)
//		{
//			outData.Add(CNurbsCurveData(inData[j]));
//		}
//	}
//	CStatus st;
//	if ( (st= NurbsCurveList(out.GetGeometry()).Set(outData, siSINurbs)) != CStatus::OK )
//	{
//		return st;
//	}
//	return CStatus::OK;
//}
//
//
