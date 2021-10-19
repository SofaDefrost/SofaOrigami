/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once
#include <sofa/defaulttype/BaseMatrix.h>
#include <SofaOrigami/OrigamiForceField.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/type/Vec.h>
#include <sofa/helper/AdvancedTimer.h>

#include <sofa/core/visual/VisualParams.h>

namespace sofa::component::interactionforcefield
{

template<class DataTypes>
OrigamiForceField<DataTypes>::OrigamiForceField(double ks, double kd)
    : OrigamiForceField<DataTypes>(nullptr, nullptr, ks, kd)
{
}

template<class DataTypes>
OrigamiForceField<DataTypes>::OrigamiForceField(MechanicalState* object1, MechanicalState* object2, double ks, double kd)
    : SpringForceField<DataTypes>(object1, object2, ks, kd)
    , d_indices1(initData(&d_indices1, "indices1", "Indices of the source points on the first model"))
    , d_indices2(initData(&d_indices2, "indices2", "Indices of the fixed points on the second model"))
    , d_length(initData(&d_length, 0.0, "length", "uniform length of all springs"))
    , angleTarget(initData(&angleTarget, 0.0, "angleTarget", "fold angle to reach"))
    , origamiEdges(initData(&origamiEdges,"origamiSpring","pairs of indices + faces opposite indices, stiffness, damping, rest length"))
{
    this->addUpdateCallback("updateSprings", { &d_indices1, &d_indices2, &d_length, &this->ks, &this->kd}, [this](const core::DataTracker& t)
    {
        SOFA_UNUSED(t);
        createSpringsFromInputs();
        return sofa::core::objectmodel::ComponentState::Valid;
    }, {&this->springs});
}


template<class DataTypes>
void OrigamiForceField<DataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);

    if (d_indices1.isSet() && d_indices2.isSet())
    {
        createSpringsFromInputs();
    }
    this->SpringForceField<DataTypes>::init();

    const type::vector<OrigamiEdge>& springs= this->origamiEdges.getValue();
    Triangles.resize(0);
    for (unsigned int i=0; i<springs.size(); i++ )
    {
        int t1 = springs[i].m1;
        int t2 = springs[i].m2;
        int t3 = springs[i].m3;
        int t4 = springs[i].m4;
        std::vector<int> triangleLeft, triangleRight;
        triangleLeft.resize(3);
        triangleRight.resize(3);
        if (t1>t2)
        {
            int tmp=t1;
            t1 = t2;
            t2=tmp;
        }
        if (t3<t1)
        {
            triangleLeft[0] = t3;
            triangleLeft[1] = t1;
            triangleLeft[2] = t2;
        }
        else if (t3<t2)
        {
            triangleLeft[0] = t1;
            triangleLeft[1] = t3;
            triangleLeft[2] = t2;
        }
        else
        {
            triangleLeft[0] = t1;
            triangleLeft[1] = t2;
            triangleLeft[2] = t3;
        }
        if (t4<t1)
        {
            triangleRight[0] = t4;
            triangleRight[1] = t1;
            triangleRight[2] = t2;
        }
        else if (t4<t2)
        {
            triangleRight[0] = t1;
            triangleRight[1] = t4;
            triangleRight[2] = t2;
        }
        else
        {
            triangleRight[0] = t1;
            triangleRight[1] = t2;
            triangleRight[2] = t4;
        }
//        msg_warning() << triangleLeft[0] << " " << triangleLeft[1] << " " << triangleLeft[2] << " ---- " << triangleRight[0] << " " << triangleRight[1] << " " << triangleRight[2];
        int nbTriangles = Triangles.size();
        bool addLeft = true, addRight = true;
        for (unsigned int j=0; j<Triangles.size();j++)
        {
            if (triangleLeft[0] != -1)
            {
                if ( (Triangles[j][0] == triangleLeft[0]) && (Triangles[j][1] == triangleLeft[1]) && (Triangles[j][2] == triangleLeft[2]) )
                {
//                    msg_warning() << "Left triangle Already in!"<< triangleLeft[0] << " " << triangleLeft[1] << " " << triangleLeft[2];
                    addLeft = false;
                }
            }
            else
                addLeft = false;

            if (triangleRight[0] != -1)
            {
                if ( (Triangles[j][0] == triangleRight[0]) && (Triangles[j][1] == triangleRight[1]) && (Triangles[j][2] == triangleRight[2]) )
                {
//                    msg_warning() << "Right triangle Already in!"<< triangleRight[0] << " " << triangleRight[1] << " " << triangleRight[2];
                    addRight = false;
                }
            }
            else
                addRight = false;
        }
        if (triangleLeft[0] == -1)
            addLeft = false;
        if (triangleRight[0] == -1)
            addRight = false;
        if (addLeft)
        {
            nbTriangles++;
            Triangles.resize(nbTriangles);
            Triangles[nbTriangles-1] = triangleLeft;
//            msg_warning() << "Adding left triangle:" << triangleLeft[0] << " " << triangleLeft[1] << " " << triangleLeft[2];
        }
        if (addRight)
        {
            nbTriangles++;
            Triangles.resize(nbTriangles);
            Triangles[nbTriangles-1] = triangleRight;
//            msg_warning() << "Adding right triangle:" << triangleRight[0] << " " << triangleRight[1] << " " << triangleRight[2];

        }
//        Triangles[]
//        msg_warning() << "My springs at init" <<origamiSprings[i].m1;
    }





    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}

template<class DataTypes>
void OrigamiForceField<DataTypes>::addFaceConstraints(
        Real& potentialEnergy,
        VecDeriv& f1,
        const  VecCoord& p1,
        const VecDeriv& v1,
        VecDeriv& f2,
        const  VecCoord& p2,
        const  VecDeriv& v2,
        int i,
        const std::vector<int> triangleInOrigami)
{
//    msg_warning() << "Staaaaaaaaaaaaaaaaaaaart of face constraiiiiiiiiiiiiiiiiiiiiiint";
    int a = triangleInOrigami[0];
    int b = triangleInOrigami[1];
    int c = triangleInOrigami[2];
//    msg_warning() << "a is : " << a << " b is: " << b << " c is : " << c;
    Real kFace = 0.0;
    typename DataTypes::CPos ab = DataTypes::getCPos(p1[b])-DataTypes::getCPos(p1[a]);
    typename DataTypes::CPos ac = DataTypes::getCPos(p1[c])-DataTypes::getCPos(p1[a]);
    typename DataTypes::CPos bc = DataTypes::getCPos(p2[c])-DataTypes::getCPos(p1[b]);
    Real alphaA = acos( dot(ab,ac)/(ab.norm()*ac.norm()) );
    Real alphaB = acos( dot(-ab,bc)/(ab.norm()*bc.norm()) );
    Real alphaC = acos( dot(-ac,-bc)/(ac.norm()*bc.norm()) );
//    msg_warning() << "the alphas from Faces: " << alphaA << " " << alphaB << " " << alphaC;

    typename DataTypes::CPos normal=cross(ab,ac);
    normal = normal/normal.norm();
//    msg_warning() << "the normal: " << normal;
    bool direct;
    if (normal[2] > 0)
        direct = true;
    else
    {
        direct = false;
        normal = -normal;
    }

    typename DataTypes::CPos directionForA;
    typename DataTypes::CPos directionForB;
    typename DataTypes::CPos directionForC;

    // /////////////////////////////////
    // Forces due to angle A
    // /////////////////////////////////
    Real alphaA0;
    if (alphaA<1.1)
        alphaA0 = M_PI_4;
    else
        alphaA0 = M_PI_2;

    Real forceFaceVal;
    forceFaceVal = -kFace*(alphaA-alphaA0);
//    msg_warning() << "forceFaceVal: " << forceFaceVal;
    if (direct)
        directionForA = cross(normal,ab)/dot(ab,ab) - cross(normal,ac)/dot(ac,ac);
    else
        directionForA = -cross(normal,ab)/dot(ab,ab) + cross(normal,ac)/dot(ac,ac);
    directionForB = cross(normal,ab)/dot(ab,ab);
    directionForC = cross(normal,ac)/dot(ac,ac);
//    msg_warning() << "directions: " << directionForA;
//    msg_warning() << "directions: " << directionForB;
//    msg_warning() << "directions: " << directionForC;

    DataTypes::setDPos( f1[a], DataTypes::getDPos(f1[a]) + forceFaceVal*directionForA ) ;
    DataTypes::setDPos( f1[b], DataTypes::getDPos(f1[b]) + forceFaceVal*directionForB ) ;
    DataTypes::setDPos( f1[c], DataTypes::getDPos(f1[c]) + forceFaceVal*directionForC ) ;

    // /////////////////////////////////
    // Forces due to angle B
    // /////////////////////////////////
    Real alphaB0;
    if (alphaB<1.1)
        alphaB0 = M_PI_4;
    else
        alphaB0 = M_PI_2;

    forceFaceVal = -kFace*(alphaB-alphaB0);
//    msg_warning() << "forceFaceVal: " << forceFaceVal;
    directionForA = cross(normal,-ab)/dot(ab,ab) ;
    if (direct)
        directionForB = -cross(normal,-ab)/dot(ab,ab) + cross(normal,-bc)/dot(bc,bc) ;
    else
        directionForB = cross(normal,-ab)/dot(ab,ab) - cross(normal,-bc)/dot(bc,bc) ;
    directionForC = cross(normal,bc)/dot(bc,bc);
//    msg_warning() << "directions: " << directionForA;
//    msg_warning() << "directions: " << directionForB;
//    msg_warning() << "directions: " << directionForC;

    DataTypes::setDPos( f1[a], DataTypes::getDPos(f1[a]) + forceFaceVal*directionForA ) ;
    DataTypes::setDPos( f1[b], DataTypes::getDPos(f1[b]) + forceFaceVal*directionForB ) ;
    DataTypes::setDPos( f1[c], DataTypes::getDPos(f1[c]) + forceFaceVal*directionForC ) ;

    // /////////////////////////////////
    // Forces due to angle C
    // /////////////////////////////////
    Real alphaC0;
    if (alphaC<1.1)
        alphaC0 = M_PI_4;
    else
        alphaC0 = M_PI_2;

    forceFaceVal = -kFace*(alphaC-alphaC0);
//    msg_warning() << "++++++++++++++++++++++++++++ forceFaceVal: " << forceFaceVal;
    directionForA = cross(normal,-ac)/dot(ac,ac);
    directionForB = cross(normal,-bc)/dot(bc,bc);
    if (direct)
        directionForC = cross(normal,-ac)/dot(ac,ac) - cross(normal,-bc)/dot(bc,bc);
    else
        directionForC = -cross(normal,-ac)/dot(ac,ac) + cross(normal,-bc)/dot(bc,bc);
//    msg_warning() << "directions: " << directionForA;
//    msg_warning() << "directions: " << directionForB;
//    msg_warning() << "directions: " << directionForC;
    DataTypes::setDPos( f1[a], DataTypes::getDPos(f1[a]) + forceFaceVal*directionForA ) ;
    DataTypes::setDPos( f1[b], DataTypes::getDPos(f1[b]) + forceFaceVal*directionForB ) ;
    DataTypes::setDPos( f1[c], DataTypes::getDPos(f1[c]) + forceFaceVal*directionForC ) ;

}

template<class DataTypes>
void OrigamiForceField<DataTypes>::createSpringsFromInputs()
{
    if (d_indices1.getValue().size() != d_indices2.getValue().size())
    {
        msg_error() << "Inputs indices sets sizes are different: d_indices1: " << d_indices1.getValue().size()
            << " | d_indices2 " << d_indices2.getValue().size()
            << " . No springs will be created";
        return;
    }

    msg_warning() << "createSpringsFromInputs createSpringsFromInputs createSpringsFromInputs createSpringsFromInputs createSpringsFromInputs";

    type::vector<Spring>& _springs = *this->springs.beginEdit();
    _springs.clear();

    const SetIndexArray & indices1 = d_indices1.getValue();
    const SetIndexArray & indices2 = d_indices2.getValue();

    const SReal& _ks = this->ks.getValue();
    const SReal& _kd = this->kd.getValue();
    const SReal& _length = d_length.getValue();
    for (sofa::Index i = 0; i<indices1.size(); ++i)
        _springs.push_back(Spring(indices1[i], indices2[i], _ks, _kd, _length));

    this->springs.endEdit();
}


template<class DataTypes>
void OrigamiForceField<DataTypes>::addSpringForce(
        Real& potentialEnergy,
        VecDeriv& f1,
        const  VecCoord& p1,
        const VecDeriv& v1,
        VecDeriv& f2,
        const  VecCoord& p2,
        const  VecDeriv& v2,
        sofa::Index i,
        const OrigamiEdge& spring)
{
    //    this->cpt_addForce++;
    sofa::Index a = spring.m1;
    sofa::Index b = spring.m2;
    sofa::Index leftNode = spring.m3;
    sofa::Index rightNode = spring.m4;
    Real kcrease = 200;

    /// Get the positional part out of the dofs.
    typename DataTypes::CPos u = DataTypes::getCPos(p2[b])-DataTypes::getCPos(p1[a]);
    Real d = u.norm();
    if( spring.enabled && d>1.0e-9 && (!spring.elongationOnly || d>spring.initpos))
    {
        // F =   k_s.(l-l_0 ).U + k_d((V_b - V_a).U).U = f.U   where f is the intensity and U the direction
        Real inverseLength = 1.0f/d;
        u *= inverseLength;
        Real elongation = (Real)(d - spring.initpos);
        potentialEnergy += elongation * elongation * spring.ks / 2;
        typename DataTypes::DPos relativeVelocity = DataTypes::getDPos(v2[b])-DataTypes::getDPos(v1[a]);
        Real elongationVelocity = dot(u,relativeVelocity);
        Real forceIntensity = (Real)(spring.ks*elongation+spring.kd*elongationVelocity);
        typename DataTypes::DPos force = u*forceIntensity;
        if(leftNode != -1 && rightNode !=-1)
        {
            typename DataTypes::CPos toRightNode = DataTypes::getCPos(p1[rightNode])-DataTypes::getCPos(p1[a]);
            typename DataTypes::CPos toRightNode2 = DataTypes::getCPos(p1[rightNode])-DataTypes::getCPos(p1[b]);
            typename DataTypes::CPos normalRightTriangle=cross(toRightNode,toRightNode2);
            normalRightTriangle = normalRightTriangle/normalRightTriangle.norm();
            msg_warning() << "normalRightTriangle " << normalRightTriangle;
            toRightNode = toRightNode - u*dot(toRightNode,u);

            typename DataTypes::CPos toLeftNode = DataTypes::getCPos(p1[leftNode])-DataTypes::getCPos(p1[a]);
            typename DataTypes::CPos toLeftNode2 = DataTypes::getCPos(p1[leftNode])-DataTypes::getCPos(p1[b]);
            typename DataTypes::CPos normalLeftTriangle=cross(toLeftNode,toLeftNode2);
            msg_warning() << "normalLeftTriangle.norm() " << normalLeftTriangle.norm();
            normalLeftTriangle = -normalLeftTriangle/normalLeftTriangle.norm();
            msg_warning() << "normalLeftTriangle " << normalLeftTriangle;
            msg_warning() << "toLeftNode " << toLeftNode;
            msg_warning() << "toLeftNode2 " << toLeftNode2;
            toLeftNode = toLeftNode - u*dot(toLeftNode,u);

            msg_warning() << "distanceRight " << toRightNode.norm();
            msg_warning() << "distanceLeft " << toLeftNode.norm();
            Real theta = acos( std::min(dot(toRightNode,toLeftNode)/(toRightNode.norm()*toLeftNode.norm()),1.0 ));
            Real thetaFromNormals = acos( std::min(dot(normalRightTriangle,normalLeftTriangle)/(normalRightTriangle.norm()*normalLeftTriangle.norm()) ,1.0) );
            msg_warning() << "theta " << theta << " Theta from normals" << thetaFromNormals;
            msg_warning() << "dot(normalRightTriangle,normalLeftTriangle) " << dot(normalRightTriangle,normalLeftTriangle);
            msg_warning() << "(normalRightTriangle.norm()*normalLeftTriangle.norm()) " << (normalRightTriangle.norm()*normalLeftTriangle.norm());
            msg_warning() << " dot(normalRightTriangle,normalLeftTriangle)/(normalRightTriangle.norm()*normalLeftTriangle.norm())  " <<  dot(normalRightTriangle,normalLeftTriangle)/(normalRightTriangle.norm()*normalLeftTriangle.norm()) ;
            msg_warning() << "acos(1)" << std::min(1.0, 1.1);
//            msg_warning() << "normalLeftTriangle.norm() " << normalLeftTriangle.norm();
//            msg_warning() << "normalLeftTriangle.norm() " << normalLeftTriangle.norm();
            if (spring.foldType == -1)
            {
                if (dot(cross(normalRightTriangle,normalLeftTriangle),u) < 0)
                {
                    theta = M_PI + thetaFromNormals;
                }
                else
                {
                    theta = M_PI - thetaFromNormals;
                }
            }
            else if (spring.foldType == 1)
            {
                if (dot(cross(normalRightTriangle,normalLeftTriangle),u) < 0)
                {
                    theta = M_PI - thetaFromNormals;
                }
                else
                {
                    theta = M_PI + thetaFromNormals;
                }
            }
            else
                theta = angleTarget.getValue();
            typename DataTypes::CPos normalToNormals = cross(normalRightTriangle,normalLeftTriangle);

            msg_warning() << "theta " << theta << " Theta from normals" << thetaFromNormals;

            typename DataTypes::CPos v31 = DataTypes::getCPos(p1[rightNode])-DataTypes::getCPos(p1[a]);
            typename DataTypes::CPos v41 = DataTypes::getCPos(p1[rightNode])-DataTypes::getCPos(p1[b]);
            typename DataTypes::CPos v32 = DataTypes::getCPos(p1[leftNode])-DataTypes::getCPos(p1[a]);
            typename DataTypes::CPos v42 = DataTypes::getCPos(p1[leftNode])-DataTypes::getCPos(p1[b]);

            Real alphaRightDown = acos( dot(v31,u)/(v31.norm()*u.norm()) );
            Real alphaRightUp = acos( dot(-u,v41)/(v41.norm()*u.norm()) );
            Real alphaLeftDown = acos( dot(u,v32)/(v32.norm()*u.norm()) );
            Real alphaLeftUp = acos( dot(v42,-u)/(v42.norm()*u.norm()) );
            msg_warning() << "the alphas: " << alphaRightDown << " " << alphaRightUp << " " << alphaLeftDown << " " << alphaLeftUp;

            typename DataTypes::CPos ab = DataTypes::getCPos(p1[b])-DataTypes::getCPos(p1[a]);
            Real forceCreaseVal;
            if (spring.foldType == -1)
                forceCreaseVal = -kcrease*spring.initpos*(theta - angleTarget.getValue());
            else if (spring.foldType == 1)
                forceCreaseVal = kcrease*spring.initpos*(theta - angleTarget.getValue());
            else
                forceCreaseVal = 0.0;
            msg_warning() << "+++++++++++++++++++++++ crease force:" << forceCreaseVal << " Norm of ab: " << ab.norm();
            typename DataTypes::DPos forceRightCrease = forceCreaseVal* normalRightTriangle/toRightNode.norm();
            typename DataTypes::DPos forceLeftCrease = forceCreaseVal* normalLeftTriangle/toLeftNode.norm();
            typename DataTypes::DPos forceCreaseA = forceCreaseVal* ( -tan(M_PI_2 - alphaRightUp) / (tan(M_PI_2 -alphaRightDown) + tan(M_PI_2 -alphaRightUp)) ) * normalRightTriangle/toRightNode.norm() + forceCreaseVal* ( -tan(M_PI_2 - alphaLeftUp) / (tan(M_PI_2 -alphaLeftDown) + tan(M_PI_2 -alphaLeftUp)) ) * normalLeftTriangle/toLeftNode.norm();
            typename DataTypes::DPos forceCreaseB = forceCreaseVal* ( -tan(M_PI_2 - alphaRightDown) / (tan(M_PI_2 -alphaRightUp) + tan(M_PI_2 -alphaRightDown)) ) * normalRightTriangle/toRightNode.norm() + forceCreaseVal* ( -tan(M_PI_2 - alphaLeftDown) / (tan(M_PI_2 -alphaLeftUp) + tan(M_PI_2 -alphaLeftDown)) ) * normalLeftTriangle/toLeftNode.norm();

            DataTypes::setDPos( f1[leftNode], DataTypes::getDPos(f1[leftNode]) + forceLeftCrease ) ;
            DataTypes::setDPos( f1[rightNode], DataTypes::getDPos(f1[rightNode]) + forceRightCrease ) ;
            DataTypes::setDPos( f1[a], DataTypes::getDPos(f1[a]) + forceCreaseA ) ;
            DataTypes::setDPos( f1[b], DataTypes::getDPos(f2[b]) + forceCreaseB ) ;



        }













        DataTypes::setDPos( f1[a], DataTypes::getDPos(f1[a]) + force ) ;
        DataTypes::setDPos( f2[b], DataTypes::getDPos(f2[b]) - force ) ;

        // Compute stiffness dF/dX
        // The force change dF comes from length change dl and unit vector change dU:
        // dF = k_s.dl.U + f.dU
        // dU = 1/l.(I-U.U^T).dX   where dX = dX_1 - dX_0  and I is the identity matrix
        // dl = U^T.dX
        // dF = k_s.U.U^T.dX + f/l.(I-U.U^T).dX = ((k_s-f/l).U.U^T + f/l.I).dX
        Mat& m = this->dfdx[i];
        Real tgt = forceIntensity * inverseLength;
        for(sofa::Index j=0; j<N; ++j )
        {
            for(sofa::Index k=0; k<N; ++k )
            {
                m[j][k] = ((Real)spring.ks-tgt) * u[j] * u[k];
            }
            m[j][j] += tgt;
        }
    }
    else // null length, no force and no stiffness
    {
        Mat& m = this->dfdx[i];
        for(sofa::Index j=0; j<N; ++j )
        {
            for(sofa::Index k=0; k<N; ++k )
            {
                m[j][k] = 0;
            }
        }
    }
}

template<class DataTypes>
void OrigamiForceField<DataTypes>::addSpringDForce(VecDeriv& df1,const  VecDeriv& dx1, VecDeriv& df2,const  VecDeriv& dx2, sofa::Index i, const Spring& spring, double kFactor, double /*bFactor*/)
{
    const sofa::Index a = spring.m1;
    const sofa::Index b = spring.m2;
    const typename DataTypes::CPos d = DataTypes::getDPos(dx2[b]) - DataTypes::getDPos(dx1[a]);
    typename DataTypes::DPos dforce = this->dfdx[i]*d;

    dforce *= kFactor;

    DataTypes::setDPos( df1[a], DataTypes::getDPos(df1[a]) + dforce ) ;
    DataTypes::setDPos( df2[b], DataTypes::getDPos(df2[b]) - dforce ) ;
}

template<class DataTypes>
void OrigamiForceField<DataTypes>::addForce(const core::MechanicalParams* /*mparams*/, DataVecDeriv& data_f1, DataVecDeriv& data_f2, const DataVecCoord& data_x1, const DataVecCoord& data_x2, const DataVecDeriv& data_v1, const DataVecDeriv& data_v2 )
{
    VecDeriv&       f1 = *data_f1.beginEdit();
    const VecCoord& x1 =  data_x1.getValue();
    const VecDeriv& v1 =  data_v1.getValue();
    VecDeriv&       f2 = *data_f2.beginEdit();
    const VecCoord& x2 =  data_x2.getValue();
    const VecDeriv& v2 =  data_v2.getValue();

    const type::vector<OrigamiEdge>& springs= this->origamiEdges.getValue();
    this->dfdx.resize(springs.size());
    f1.resize(x1.size());
    f2.resize(x2.size());
    this->m_potentialEnergy = 0;
    for (sofa::Index i=0; i<springs.size(); i++)
    {
        this->addSpringForce(this->m_potentialEnergy,f1,x1,v1,f2,x2,v2, i, springs[i]);
    }
    for (unsigned int i=0; i<Triangles.size(); i++)
    {
//        msg_warning() << "adding force for triangle :" << i;
        this->addFaceConstraints(this->m_potentialEnergy,f1,x1,v1,f2,x2,v2, i, Triangles[i]);
//        msg_warning() << "Triangles:" << Triangles[i][0]<< Triangles[i][1]<< Triangles[i][2];
    }
    data_f1.endEdit();
    data_f2.endEdit();
}

template<class DataTypes>
void OrigamiForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams, DataVecDeriv& data_df1, DataVecDeriv& data_df2, const DataVecDeriv& data_dx1, const DataVecDeriv& data_dx2)
{
    VecDeriv&        df1 = *data_df1.beginEdit();
    VecDeriv&        df2 = *data_df2.beginEdit();
    const VecDeriv&  dx1 =  data_dx1.getValue();
    const VecDeriv&  dx2 =  data_dx2.getValue();
    Real kFactor       =  (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(mparams,this->rayleighStiffness.getValue());
    Real bFactor       =  (Real)sofa::core::mechanicalparams::bFactor(mparams);

    const type::vector<Spring>& springs = this->springs.getValue();
    df1.resize(dx1.size());
    df2.resize(dx2.size());

    for (sofa::Index i=0; i<springs.size(); i++)
    {
        this->addSpringDForce(df1,dx1,df2,dx2, i, springs[i], kFactor, bFactor);
    }

    data_df1.endEdit();
    data_df2.endEdit();
}




template<class DataTypes>
void OrigamiForceField<DataTypes>::addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    Real kFact = (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(mparams,this->rayleighStiffness.getValue());
    if (this->mstate1 == this->mstate2)
    {
        sofa::core::behavior::MultiMatrixAccessor::MatrixRef mat = matrix->getMatrix(this->mstate1);
        if (!mat) return;
        const sofa::type::vector<Spring >& ss = this->springs.getValue();
        const sofa::Size n = ss.size() < this->dfdx.size() ? sofa::Size(ss.size()) : sofa::Size(this->dfdx.size());
        for (sofa::Index e=0; e<n; e++)
        {
            const Spring& s = ss[e];
            sofa::Index p1 = mat.offset+Deriv::total_size*s.m1;
            sofa::Index p2 = mat.offset+Deriv::total_size*s.m2;
            const Mat& m = this->dfdx[e];
            for(sofa::Index i=0; i<N; i++)
            {
                for (sofa::Index j=0; j<N; j++)
                {
                    Real k = (Real)(m[i][j]*kFact);
                    mat.matrix->add(p1+i,p1+j, -k);
                    mat.matrix->add(p1+i,p2+j, k);
                    mat.matrix->add(p2+i,p1+j, k);//or mat->add(p1+j,p2+i, k);
                    mat.matrix->add(p2+i,p2+j, -k);
                }
            }
        }
    }
    else
    {
        sofa::core::behavior::MultiMatrixAccessor::MatrixRef mat11 = matrix->getMatrix(this->mstate1);
        sofa::core::behavior::MultiMatrixAccessor::MatrixRef mat22 = matrix->getMatrix(this->mstate2);
        sofa::core::behavior::MultiMatrixAccessor::InteractionMatrixRef mat12 = matrix->getMatrix(this->mstate1, this->mstate2);
        sofa::core::behavior::MultiMatrixAccessor::InteractionMatrixRef mat21 = matrix->getMatrix(this->mstate2, this->mstate1);

        if (!mat11 && !mat22 && !mat12 && !mat21) return;
        const sofa::type::vector<Spring >& ss = this->springs.getValue();
        const sofa::Size n = ss.size() < this->dfdx.size() ? sofa::Size(ss.size()) : sofa::Size(this->dfdx.size());
        for (sofa::Index e=0; e<n; e++)
        {
            const Spring& s = ss[e];
            unsigned p1 = /*mat.offset+*/Deriv::total_size*s.m1;
            unsigned p2 = /*mat.offset+*/Deriv::total_size*s.m2;
            Mat m = this->dfdx[e]* (Real) kFact;
            if (mat11)
            {
                for(sofa::Index i=0; i<N; i++)
                {
                    for (sofa::Index j=0; j<N; j++)
                    {
                        mat11.matrix->add(mat11.offset+p1+i,mat11.offset+p1+j, -(Real)m[i][j]);
                    }
                }
            }
            if (mat12)
            {
                for(sofa::Index i=0; i<N; i++)
                {
                    for (sofa::Index j=0; j<N; j++)
                    {
                        mat12.matrix->add(mat12.offRow+p1+i,mat12.offCol+p2+j,  (Real)m[i][j]);
                    }
                }
            }
            if (mat21)
            {
                for(sofa::Index i=0; i<N; i++)
                {
                    for (sofa::Index j=0; j<N; j++)
                    {
                        mat21.matrix->add(mat21.offRow+p2+i,mat21.offCol+p1+j,  (Real)m[i][j]);
                    }
                }
            }
            if (mat22)
            {
                for(sofa::Index i=0; i<N; i++)
                {
                    for (sofa::Index j=0; j<N; j++)
                    {
                        mat22.matrix->add(mat22.offset+p2+i,mat22.offset+p2+j, -(Real)m[i][j]);
                    }
                }
            }
        }
    }

}

template<class DataTypes>
void OrigamiForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    using namespace sofa::defaulttype;

    if (!((this->mstate1 == this->mstate2)?vparams->displayFlags().getShowForceFields():vparams->displayFlags().getShowInteractionForceFields())) return;
    const VecCoord& p1 =this->mstate1->read(core::ConstVecCoordId::position())->getValue();
    const VecCoord& p2 =this->mstate2->read(core::ConstVecCoordId::position())->getValue();

    std::vector< type::Vector3 > points[5];
    bool external = (this->mstate1!=this->mstate2);
    //const helper::vector<Spring>& springs = this->springs.getValue();
    const type::vector<OrigamiEdge>& springs = this->origamiEdges.getValue();
    std::vector<sofa::type::Vec4f> colorVector;
    std::vector<sofa::type::Vector3> vertices;

    for (unsigned int i=0; i<springs.size(); i++)
    {
        if (!springs[i].enabled) continue;
        Real d = (p2[springs[i].m2]-p1[springs[i].m1]).norm();
        type::Vector3 point2,point1,point3, point4;
        point1 = DataTypes::getCPos(p1[springs[i].m1]);
        point2 = DataTypes::getCPos(p2[springs[i].m2]);
        point3 = DataTypes::getCPos(p1[springs[i].m3]);
        point4 = DataTypes::getCPos(p2[springs[i].m4]);

        if (external)
        {
            if (springs[i].foldType == 1)
            {
                points[0].push_back(point1);
                points[0].push_back(point2);
            }
            else
            {
                points[1].push_back(point1);
                points[1].push_back(point2);
            }
        }
        else
        {
            if (springs[i].foldType == 1)
            {
                points[2].push_back(point1);
                points[2].push_back(point2);
            }
            else
            {
                points[3].push_back(point1);
                points[3].push_back(point2);
            }
        }
//        points[4].push_back(point3);
//        points[4].push_back(point4);


//        if (springs[i].m3 != -1){
//            colorVector.push_back(sofa::defaulttype::RGBAColor(0,1,0,1));
//            vertices.push_back(sofa::defaulttype::Vector3(point1));
//            colorVector.push_back(sofa::defaulttype::RGBAColor(0,0.5,0.5,1));
//            vertices.push_back(sofa::defaulttype::Vector3(point2));
//            colorVector.push_back(sofa::defaulttype::RGBAColor(0,0,1,1));
//            vertices.push_back(sofa::defaulttype::Vector3(point3));

//        }
        if (springs[i].m4 != -1){
            colorVector.push_back(sofa::type::RGBAColor(0,1,0,1));
            vertices.push_back(sofa::type::Vector3(point1));
            colorVector.push_back(sofa::type::RGBAColor(0,0.5,0.5,1));
            vertices.push_back(sofa::type::Vector3(point2));
            colorVector.push_back(sofa::type::RGBAColor(0,0,1,1));
            vertices.push_back(sofa::type::Vector3(point4));
        }


    }

    const type::RGBAColor c0 = type::RGBAColor::red();
    const type::RGBAColor c1 = type::RGBAColor::green();
    const type::RGBAColor c2 {1.0f, 0.5f, 0.0f, 1.0f };
    const type::RGBAColor c3{ 0.0f, 1.0f, 0.5f, 1.0f };

    if (this->showArrowSize.getValue()==0 || this->drawMode.getValue() == 0)
    {
        vparams->drawTool()->drawLines(points[0], 1, c0);
        vparams->drawTool()->drawLines(points[1], 1, c1);
        vparams->drawTool()->drawLines(points[2], 1, c2);
        vparams->drawTool()->drawLines(points[3], 1, c3);
    }
    else if (this->drawMode.getValue() == 1)
    {
        const auto numLines0=points[0].size()/2;
        const auto numLines1=points[1].size()/2;
        const auto numLines2=points[2].size()/2;
        const auto numLines3=points[3].size()/2;

        for (unsigned int i=0; i<numLines0; ++i) vparams->drawTool()->drawCylinder(points[0][2*i+1], points[0][2*i], this->showArrowSize.getValue(), c0);
        for (unsigned int i=0; i<numLines1; ++i) vparams->drawTool()->drawCylinder(points[1][2*i+1], points[1][2*i], this->showArrowSize.getValue(), c1);
        for (unsigned int i=0; i<numLines2; ++i) vparams->drawTool()->drawCylinder(points[2][2*i+1], points[2][2*i], this->showArrowSize.getValue(), c2);
        for (unsigned int i=0; i<numLines3; ++i) vparams->drawTool()->drawCylinder(points[3][2*i+1], points[3][2*i], this->showArrowSize.getValue(), c3);

    }
    else if (this->drawMode.getValue() == 2)
    {
        const auto numLines0=points[0].size()/2;
        const auto numLines1=points[1].size()/2;
        const auto numLines2=points[2].size()/2;
        const auto numLines3=points[3].size()/2;

        for (unsigned int i=0; i<numLines0; ++i) vparams->drawTool()->drawArrow(points[0][2*i+1], points[0][2*i], this->showArrowSize.getValue(), c0);
        for (unsigned int i=0; i<numLines1; ++i) vparams->drawTool()->drawArrow(points[1][2*i+1], points[1][2*i], this->showArrowSize.getValue(), c1);
        for (unsigned int i=0; i<numLines2; ++i) vparams->drawTool()->drawArrow(points[2][2*i+1], points[2][2*i], this->showArrowSize.getValue(), c2);
        for (unsigned int i=0; i<numLines3; ++i) vparams->drawTool()->drawArrow(points[3][2*i+1], points[3][2*i], this->showArrowSize.getValue(), c3);
    }
    else
    {
        msg_error()<< "No proper drawing mode found!";
    }



}

} // namespace sofa::component::interactionforcefield
