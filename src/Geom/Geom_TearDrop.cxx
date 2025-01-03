// Created on: 1993-03-10
// Created by: JCV
// Copyright (c) 1993-1999 Matra Datavision
// Copyright (c) 1999-2014 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.


#include <ElCLib.hxx>
#include <Geom_TearDrop.hxx>
#include <Geom_Geometry.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Teard.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>
#include <Standard_ConstructionError.hxx>
#include <Standard_RangeError.hxx>
#include <Standard_Type.hxx>

IMPLEMENT_STANDARD_RTTIEXT(Geom_TearDrop,Geom_Conic)

typedef gp_Ax1  Ax1;
typedef gp_Ax2  Ax2;
typedef gp_Pnt  Pnt;
typedef gp_Vec  Vec;
typedef gp_Trsf Trsf;
typedef gp_XYZ  XYZ;

//=======================================================================
//function : Copy
//purpose  : 
//=======================================================================

Handle(Geom_Geometry) Geom_TearDrop::Copy() const
{
  Handle(Geom_TearDrop) E;
  E = new Geom_TearDrop(pos, majorRadius, minorRadius);
  return E;
}




//=======================================================================
//function : Geom_TearDrop
//purpose  : 
//=======================================================================

Geom_TearDrop::Geom_TearDrop (const gp_Teard& E)
  : majorRadius (E.MajorRadius()), minorRadius (E.MinorRadius()) 
{
  pos = E.Position ();
}


//=======================================================================
//function : Geom_TearDrop
//purpose  : 
//=======================================================================

Geom_TearDrop::Geom_TearDrop ( const Ax2& A, 
                             const Standard_Real MajorRadius,
                             const Standard_Real MinorRadius) 
  : majorRadius (MajorRadius), minorRadius (MinorRadius) {

   if (MajorRadius < MinorRadius || MinorRadius < 0.0 ) {
     throw Standard_ConstructionError();
   }
   pos = A;
}


//=======================================================================
//function : IsClosed
//purpose  : 
//=======================================================================

Standard_Boolean Geom_TearDrop::IsClosed () const      { return Standard_True; }

//=======================================================================
//function : IsPeriodic
//purpose  : 
//=======================================================================

Standard_Boolean Geom_TearDrop::IsPeriodic () const    { return Standard_True; }

//=======================================================================
//function : FirstParameter
//purpose  : 
//=======================================================================

Standard_Real Geom_TearDrop::FirstParameter () const   { return 0.0; }

//=======================================================================
//function : LastParameter
//purpose  : 
//=======================================================================

Standard_Real Geom_TearDrop::LastParameter () const    { return 2.0 * M_PI; }

//=======================================================================
//function : MajorRadius
//purpose  : 
//=======================================================================

Standard_Real Geom_TearDrop::MajorRadius () const      { return majorRadius; }

//=======================================================================
//function : MinorRadius
//purpose  : 
//=======================================================================

Standard_Real Geom_TearDrop::MinorRadius () const      { return minorRadius; }

//=======================================================================
//function : SetElips
//purpose  : 
//=======================================================================

void Geom_TearDrop::SetElips (const gp_Teard& E) {

  majorRadius = E.MajorRadius();
  minorRadius = E.MinorRadius();
  pos = E.Position();
}


//=======================================================================
//function : SetMajorRadius
//purpose  : 
//=======================================================================

void Geom_TearDrop::SetMajorRadius (const Standard_Real MajorRadius) {

  if (MajorRadius < minorRadius)  throw Standard_ConstructionError();
  else                            majorRadius = MajorRadius; 
}


//=======================================================================
//function : SetMinorRadius
//purpose  : 
//=======================================================================

void Geom_TearDrop::SetMinorRadius (const Standard_Real MinorRadius) {

   if (MinorRadius < 0 || majorRadius < MinorRadius) {
     throw Standard_ConstructionError();
   }
   else { minorRadius = MinorRadius; }
}


//=======================================================================
//function : Elips
//purpose  : 
//=======================================================================

gp_Teard Geom_TearDrop::Elips () const {

  return gp_Teard(pos, majorRadius, minorRadius);
}


//=======================================================================
//function : ReversedParameter
//purpose  : 
//=======================================================================

Standard_Real Geom_TearDrop::ReversedParameter( const Standard_Real U) const 
{
  return ( 2. * M_PI - U);
}


//=======================================================================
//function : Directrix1
//purpose  : 
//=======================================================================

Ax1 Geom_TearDrop::Directrix1 () const {

    gp_Teard Ev (pos, majorRadius, minorRadius);
   return Ev.Directrix1();
}


//=======================================================================
//function : Directrix2
//purpose  : 
//=======================================================================

Ax1 Geom_TearDrop::Directrix2 () const {

    gp_Teard Ev (pos, majorRadius, minorRadius);
  return Ev.Directrix2();
}


//=======================================================================
//function : D0
//purpose  : 
//=======================================================================

void Geom_TearDrop::D0 (const Standard_Real U, gp_Pnt& P) const {

  P = ElCLib::TearDropValue(U, pos, majorRadius, minorRadius);
}


//=======================================================================
//function : D1
//purpose  : 
//=======================================================================

void Geom_TearDrop::D1 (const Standard_Real U, Pnt& P, Vec& V1) const {

  ElCLib::TearDropD1 (U, pos, majorRadius, minorRadius, P, V1);
}


//=======================================================================
//function : D2
//purpose  : 
//=======================================================================

void Geom_TearDrop::D2 (const Standard_Real U, Pnt& P, Vec& V1, Vec& V2) const {

  ElCLib::TearDropD2 (U, pos, majorRadius, minorRadius, P, V1, V2);
}


//=======================================================================
//function : D3
//purpose  : 
//=======================================================================

void Geom_TearDrop::D3 (const Standard_Real U, Pnt& P, Vec& V1, Vec& V2, Vec& V3) const {

  ElCLib::EllipseD3 (U, pos, majorRadius, minorRadius, P, V1, V2, V3);
}


//=======================================================================
//function : DN
//purpose  : 
//=======================================================================

Vec Geom_TearDrop::DN (const Standard_Real U, const Standard_Integer N) const {

   Standard_RangeError_Raise_if (N < 1, " ");
   return ElCLib::EllipseDN (U, pos, majorRadius, minorRadius, N);
}


//=======================================================================
//function : Eccentricity
//purpose  : 
//=======================================================================

Standard_Real Geom_TearDrop::Eccentricity () const {

  if (majorRadius == 0.0) { return 0.0; }
  else {
    return (Sqrt(majorRadius*majorRadius-minorRadius*minorRadius))/majorRadius;
  }
}


//=======================================================================
//function : Focal
//purpose  : 
//=======================================================================

Standard_Real Geom_TearDrop::Focal () const {

  return 2.0 * Sqrt(majorRadius * majorRadius - minorRadius * minorRadius);
}


//=======================================================================
//function : Focus1
//purpose  : 
//=======================================================================

Pnt Geom_TearDrop::Focus1 () const {

  Standard_Real C = Sqrt (majorRadius * majorRadius - minorRadius * minorRadius);
  Standard_Real Xp, Yp, Zp, Xd, Yd, Zd;
  pos.Location().Coord (Xp, Yp, Zp);
  pos.XDirection().Coord (Xd, Yd, Zd);
  return Pnt (Xp + C * Xd,  Yp + C * Yd,  Zp + C * Zd);
}


//=======================================================================
//function : Focus2
//purpose  : 
//=======================================================================

Pnt Geom_TearDrop::Focus2 () const {

  Standard_Real C = Sqrt (majorRadius * majorRadius - minorRadius * minorRadius);
  Standard_Real Xp, Yp, Zp, Xd, Yd, Zd;
  pos.Location().Coord (Xp, Yp, Zp);
  pos.XDirection().Coord (Xd, Yd, Zd);
  return Pnt (Xp - C * Xd,  Yp - C * Yd,  Zp - C * Zd);
}


//=======================================================================
//function : Parameter
//purpose  : 
//=======================================================================

Standard_Real Geom_TearDrop::Parameter () const {

  if (majorRadius == 0.0)  return 0.0;
  else                     return (minorRadius * minorRadius) / majorRadius;
}


//=======================================================================
//function : Transform
//purpose  : 
//=======================================================================

void Geom_TearDrop::Transform (const Trsf& T) {

  majorRadius = majorRadius * Abs(T.ScaleFactor());
  minorRadius = minorRadius * Abs(T.ScaleFactor());
  pos.Transform(T);
}

//=======================================================================
//function : DumpJson
//purpose  : 
//=======================================================================
void Geom_TearDrop::DumpJson (Standard_OStream& theOStream, Standard_Integer theDepth) const
{
  OCCT_DUMP_TRANSIENT_CLASS_BEGIN (theOStream)

  OCCT_DUMP_BASE_CLASS (theOStream, theDepth, Geom_Conic)

  OCCT_DUMP_FIELD_VALUE_NUMERICAL (theOStream, majorRadius)
  OCCT_DUMP_FIELD_VALUE_NUMERICAL (theOStream, minorRadius)
}
