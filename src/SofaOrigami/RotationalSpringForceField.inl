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


#include "RotationalSpringForceField.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/type/RGBAColor.h>

namespace sofa::component::forcefield
{

template <class DataTypes>
RotationalSpringForceField<DataTypes>::
RotationalSpringForceField()
    : _indexedElements(nullptr)
    , _initialPoints(initData(&_initialPoints, "initialPoints", "Initial Position"))
    , m_topology(nullptr)
    , method(LARGE)
    , d_edgesAssignment(initData(&d_edgesAssignment,"edgesAssignment","Assignment of edges: Valley, Mountain, Flat or Border"))
    , d_angleTarget(initData(&d_angleTarget, 0.0, "angleTarget", "fold angle to reach"))
    , f_method(initData(&f_method,std::string("large"),"method","large: large displacements, small: small displacements"))
    , f_poisson(initData(&f_poisson,Real(0.3),"poissonRatio","Poisson ratio in Hooke's law"))
    , f_young(initData(&f_young,Real(1000.),"youngModulus","Young modulus in Hooke's law"))
    , f_thickness(initData(&f_thickness,Real(1.),"thickness","Thickness of the elements"))
//    , f_damping(initData(&f_damping,(Real)0.,"damping","Ratio damping/stiffness"))
    , f_planeStrain(initData(&f_planeStrain,false,"planeStrain","Plane strain or plane stress assumption"))
    , l_topology(initLink("topology", "link to the topology container"))
{
    f_poisson.setRequired(true);
    f_young.setRequired(true);
}

template <class DataTypes>
RotationalSpringForceField<DataTypes>::~RotationalSpringForceField()
{
}


template <class DataTypes>
void RotationalSpringForceField<DataTypes>::init()
{
    this->Inherited::init();

    // checking inputs using setter
    setMethod(f_method.getValue());
    setPoisson(f_poisson.getValue());
    setYoung(f_young.getValue());

    if (l_topology.empty())
    {
        msg_info() << "link to Topology container should be set to ensure right behavior. First Topology found in current context will be used.";
        l_topology.set(this->getContext()->getMeshTopologyLink());
    }

    m_topology = l_topology.get();
    msg_info() << "Topology path used: '" << l_topology.getLinkedPath() << "'";

    if (m_topology == nullptr)
    {
        msg_error() << "No topology component found at path: " << l_topology.getLinkedPath() << ", nor in current context: " << this->getContext()->name;
        sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    if (m_topology->getTriangles().empty() && m_topology->getQuads().empty())
    {
        msg_warning() << "No triangles found in linked Topology.";
        _indexedElements = &(m_topology->getTriangles());
    }
    else if (!m_topology->getTriangles().empty())
    {
        msg_info() << "Init using triangles mesh: " << m_topology->getTriangles().size() << " triangles.";
        _indexedElements = &(m_topology->getTriangles());
    }
    else if (m_topology->getNbQuads() > 0)
    {
        msg_info() << "Init using quads mesh: " << m_topology->getNbQuads() * 2 << " triangles.";
        sofa::core::topology::BaseMeshTopology::SeqTriangles* trias = new sofa::core::topology::BaseMeshTopology::SeqTriangles;
        int nbcubes = m_topology->getNbQuads();
        trias->reserve(nbcubes*2);
        for (int i=0; i<nbcubes; i++)
        {
            sofa::core::topology::BaseMeshTopology::Quad q = m_topology->getQuad(i);
            trias->push_back(Element(q[0],q[1],q[2]));
            trias->push_back(Element(q[0],q[2],q[3]));
        }
        _indexedElements = trias;
    }

    if (_initialPoints.getValue().size() == 0)
    {
        const VecCoord& p = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
        _initialPoints.setValue(p);
    }

    _strainDisplacements.resize(_indexedElements->size());
    _rotations.resize(_indexedElements->size());

    if (method == SMALL)
    {
        initSmall();
    }
    else
    {
        initLarge();
    }

    computeMaterialStiffnesses();
//    msg_warning() << "Get triangles tesssssssssssssssssssssssssssssssssssssssssssssst";
//    msg_warning() << "Nb edges: " << m_topology->getNbEdges();
    adjacentVertices.resize(m_topology->getNbEdges());
    edgeInitialLength.resize(m_topology->getNbEdges());
    edgeStiffness.resize(m_topology->getNbEdges());
    nodes.resize(m_topology->getNbEdges());
    forceCrease.resize(m_topology->getNbEdges());
    for (int i; i<m_topology->getNbEdges(); i++){
        int vertexDirect = -1, vertexIndirect = -1;
        auto trianglesAroundEdge = m_topology->getTrianglesAroundEdge(i);
        auto edge = m_topology->getEdge(i);
        if (trianglesAroundEdge.size() == 2){
            int posEdge0InTri, posEdge1InTri;

            auto triangle1 = m_topology->getTriangle(trianglesAroundEdge[0]);
            auto triangle2 = m_topology->getTriangle(trianglesAroundEdge[1]);
//            msg_warning() << "t1: " << triangle1 << " t2: " << triangle2;

            for (int index=0; index <3; index++)
            {
                if (triangle1[index] == edge[0])
                    posEdge0InTri = index;
                if (triangle1[index] == edge[1])
                    posEdge1InTri = index;
            }
            if (posEdge1InTri == ((posEdge0InTri + 1) % 3))
                vertexDirect = triangle1[((posEdge0InTri + 2) % 3)];
            else
                vertexIndirect = triangle1[((posEdge0InTri + 1) % 3)];
            for (int index=0; index <3; index++)
            {
                if (triangle2[index] == edge[0])
                    posEdge0InTri = index;
                if (triangle2[index] == edge[1])
                    posEdge1InTri = index;
            }
            if (posEdge1InTri == ((posEdge0InTri + 1) % 3))
                vertexDirect = triangle2[((posEdge0InTri + 2) % 3)];
            else
                vertexIndirect = triangle2[((posEdge0InTri + 1) % 3)];

        }
        adjacentVertices[i].resize(2);
        adjacentVertices[i][0] = vertexDirect;
        adjacentVertices[i][1] = vertexIndirect;
        nodes[i].resize(4);
        Coord p0,p1;
        p0[0] = m_topology->getPX(edge[0]);
        p0[1] = m_topology->getPY(edge[0]);
        p0[2] = m_topology->getPZ(edge[0]);
        p1[0] = m_topology->getPX(edge[1]);
        p1[1] = m_topology->getPY(edge[1]);
        p1[2] = m_topology->getPZ(edge[1]);
        edgeInitialLength[i] = norm(p1-p0);
    }

//for (int i; i<m_topology->getNbEdges(); i++)
//    msg_warning() << "edge:" << i << " adj vert: " << adjacentVertices[i][0] << " " << adjacentVertices[i][1];
}



template <class DataTypes>
void RotationalSpringForceField<DataTypes>::reinit()
{
    if (f_method.getValue() == "small")
        method = SMALL;
    else if (f_method.getValue() == "large")
        method = LARGE;

    if (method == SMALL)
    {
        //    initSmall();  // useful ? The rotations are recomputed later
    }
    else
    {
        initLarge(); // compute the per-element strain-displacement matrices
    }

    computeMaterialStiffnesses();
}
template <class DataTypes>
type::Mat3x3 RotationalSpringForceField<DataTypes>::diamond(Coord a, Coord b){
    return tensorProduct(a,b) + tensorProduct(b,a);
}

template <class DataTypes>
void RotationalSpringForceField<DataTypes>::addForce(const core::MechanicalParams* /* mparams */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& /* v */)
{
    VecDeriv& f1 = *f.beginEdit();
    const VecCoord& positions = x.getValue();
//    msg_warning() << "--------------------------------------------------------------------------";
//    msg_warning() << positions;
//    msg_warning() << "--------------------------------------------------------------------------";
    const VecCoord& x1 = x.getValue();
    f1.resize(positions.size());
    Real kcrease = 200;
    Real forceCreaseVal, dForcecrease;
    auto edgesAssignment = d_edgesAssignment.getValue();
    msg_warning() << "Number of edges !!!!!!!!!!!!!!!! " << m_topology->getNbEdges();
    for (int edgeNum=0; edgeNum < m_topology->getNbEdges(); edgeNum++)
    {

        auto edge = m_topology->getEdge(edgeNum);
        auto adjacent = adjacentVertices[edgeNum];
        Coord pj = positions[edge[0]];
        Coord pk = positions[edge[1]];
        Coord pi = positions[adjacent[0]];
        Coord pl = positions[adjacent[1]];

        nodes[edgeNum][0] = adjacent[0];
        nodes[edgeNum][1] = edge[0];
        nodes[edgeNum][2] = edge[1];
        nodes[edgeNum][3] = adjacent[1];
        msg_warning() << "edge Number : " <<edgeNum;
        msg_warning() << "Nodes : " << nodes[edgeNum];
        Coord rij = pi - pj;
        Coord rkj = pk - pj;
        Coord rkl = pk - pl;
        Coord m = cross(rij,rkj);
        Coord n = cross(rkj,rkl);
        Real m4 = dot(m,m)*dot(m,m);
        Real n4 = dot(n,n)*dot(n,n);







        Real cosTheta = dot(m,n)/(norm(m)*norm(n));
        if (cosTheta >= 1.0 )
            cosTheta = 1.0;
        else if (cosTheta <= -1.0 )
            cosTheta = -1.0;

        Real theta = acos( cosTheta );
        msg_warning() << "cos theta : !!!!!!!!!!!!!!!! " << dot(m,n)/(norm(m)*norm(n));
        int signMdotRkl = 1;
        if (dot(m,rkl) != 0)
            signMdotRkl = (0 < dot(m,rkl)) - (dot(m,rkl)<0);

        theta = signMdotRkl * theta;
        theta = theta - floor( theta/(2*M_PI))*2*M_PI; // Modulo 2PI
        Deriv dTheta_xi = (norm(rkj)/dot(m,m))*m;
        Deriv dTheta_xl = -(norm(rkj)/dot(n,n))*n;
        Deriv dTheta_xj = (dot(rij,rkj)/dot(rkj,rkj) - 1)*dTheta_xi - (dot(rkl,rkj)/dot(rkj,rkj))*dTheta_xl;
        Deriv dTheta_xk = (dot(rkl,rkj)/dot(rkj,rkj) - 1)*dTheta_xl - (dot(rij,rkj)/dot(rkj,rkj))*dTheta_xi;
        msg_warning() << "dot(rij,rkj): " << dot(rij,rkj);
        msg_warning() << "dot(rkl,rkj): " << dot(rkl,rkj);

//        msg_warning() << "forces : !!!!!!!!!!!!!!!! " << fxi;
//        msg_warning() << "forces : !!!!!!!!!!!!!!!! " << fxl;
//        msg_warning() << "forces : !!!!!!!!!!!!!!!! " << fxj;
//        msg_warning() << "forces : !!!!!!!!!!!!!!!! " << fxk;
        if (edgesAssignment[edgeNum] == "M")
            msg_warning() << "theta : " << theta << "angle target : " << d_angleTarget.getValue();
        else if (edgesAssignment[edgeNum] == "V")
            msg_warning() << "theta : " << theta << "angle target : " << 2*M_PI - d_angleTarget.getValue();

        if (edgesAssignment[edgeNum] == "V"){
            forceCreaseVal = -kcrease*edgeInitialLength[edgeNum]*(theta - d_angleTarget.getValue());
            msg_warning() << "force detail:" << kcrease << " " << edgeInitialLength[edgeNum] << " " << (theta - d_angleTarget.getValue());
            dForcecrease = -kcrease*edgeInitialLength[edgeNum];
        } else if (edgesAssignment[edgeNum] == "M"){
            forceCreaseVal = -kcrease*edgeInitialLength[edgeNum]*(theta - (2*M_PI - d_angleTarget.getValue()));
            dForcecrease = -kcrease*edgeInitialLength[edgeNum];
        } else {
            forceCreaseVal = 0.0;
            dForcecrease = 0.0;
        }
//        dForcecrease = 0.0;
        forceCrease[edgeNum] = forceCreaseVal;
        msg_warning() << "forceCreaseVal !!!!!!!!!!!!!!!!" << forceCreaseVal;
        if ((edgesAssignment[edgeNum] == "M") || (edgesAssignment[edgeNum] == "V")) {
            f1[adjacent[0]] += forceCreaseVal*dTheta_xi;
            f1[adjacent[1]] += forceCreaseVal*dTheta_xl;
            f1[edge[0]] += forceCreaseVal*dTheta_xj;
            f1[edge[1]] += forceCreaseVal*dTheta_xk;
            msg_warning() << "dTheta_xj !!!!!!!!!!!!!!!!" << dTheta_xj;
            msg_warning() << "dTheta_xk !!!!!!!!!!!!!!!!" << dTheta_xk;
            msg_warning() << "dTheta_xi !!!!!!!!!!!!!!!!" << dTheta_xi;
            msg_warning() << "dTheta_xl !!!!!!!!!!!!!!!!" << dTheta_xl;
        }
        msg_warning() << "Starting to compute edgeStiffness ";
        Real A = dot(rij,rkj)/dot(rkj,rkj);
        Real B = dot(rkl,rkj)/dot(rkj,rkj);
        Coord dAj = 1/dot(rkj,rkj) * ((2*A-1)*rkj - rij);
        Coord dBj = 1/dot(rkj,rkj) * (2*B*rkj - rkl);
        Coord dAk = 1/dot(rkj,rkj) * (-2*A*rkj + rij);
        Coord dBk = 1/dot(rkj,rkj) * ((1- 2*B)*rkj + rkl);
        msg_warning() << "middle compute edgeStiffness ";
        type::Mat3x3 block_ii = - norm(rkj)/m4 * diamond(m,cross(rkj,m));
        type::Mat3x3 block_ll = norm(rkj)/n4 * diamond(n,cross(rkj,n));
        type::Mat3x3 block_ik = 1/(dot(m,m)*norm(rkj)) * tensorProduct(m,rkj) + norm(rkj)/m4 * diamond(m,cross(rij,m));
        type::Mat3x3 block_lj = 1/(dot(n,n)*norm(rkj)) * tensorProduct(n,rkj) - norm(rkj)/n4 * diamond(n,cross(rkl,n));
        type::Mat3x3 block_ij = 1/(dot(m,m)*norm(rkj)) * tensorProduct(m,rkj) + norm(rkj)/m4 * diamond(m,cross(rkj-rij,m));
        type::Mat3x3 block_lk = 1/(dot(n,n)*norm(rkj)) * tensorProduct(n,rkj) - norm(rkj)/n4 * diamond(n,cross(rkj-rkl,n));
        type::Mat3x3 block_jj = tensorProduct(dTheta_xi,dAj) + (A-1)* block_ij - ( tensorProduct(dTheta_xl,dBj) + B*block_lj );
        type::Mat3x3 block_jk = tensorProduct(dTheta_xi,dAk) + (A-1)* block_ik - ( tensorProduct(dTheta_xl,dBk) + B*block_lk );
        type::Mat3x3 block_kk = tensorProduct(dTheta_xl,dBk) + (B-1)* block_lk - ( tensorProduct(dTheta_xi,dAk) + A*block_ik );

        type::Mat3x3 block_ki = block_ik.transposed();
        type::Mat3x3 block_jl = block_lj.transposed();
        type::Mat3x3 block_ji = block_ij.transposed();
        type::Mat3x3 block_kj = block_jk.transposed();
        type::Mat3x3 block_kl = block_lk.transposed();
        msg_warning() << "Blocks computed ! ";

        for (int i=0;i<12;i++){
            for (int j=0;j<12;j++){
//                edgeStiffness[edgeNum].clear();
                edgeStiffness[edgeNum][i][j] = 0;
            }
        }
//        forceCreaseVal = 0;
        if ((edgesAssignment[edgeNum] == "M") || (edgesAssignment[edgeNum] == "V")) {
            msg_warning() << "Edge Stiffness for edge: " << edgeNum;
            for (int i=0;i<3;i++){
                for (int j=0;j<3;j++){
                    edgeStiffness[edgeNum][i][j] += forceCreaseVal*block_ii[i][j] + dForcecrease*tensorProduct(dTheta_xi,dTheta_xi)[i][j];
                    edgeStiffness[edgeNum][3+i][j] += forceCreaseVal*block_ji[i][j] + dForcecrease*tensorProduct(dTheta_xj,dTheta_xi)[i][j];
                    edgeStiffness[edgeNum][i][3+j] += forceCreaseVal*block_ij[i][j] + dForcecrease*tensorProduct(dTheta_xi,dTheta_xj)[i][j];
                    edgeStiffness[edgeNum][i][6+j] += forceCreaseVal*block_ik[i][j] + dForcecrease*tensorProduct(dTheta_xi,dTheta_xk)[i][j];
                    edgeStiffness[edgeNum][6+i][j] += forceCreaseVal*block_ki[i][j] + dForcecrease*tensorProduct(dTheta_xk,dTheta_xi)[i][j];
                    edgeStiffness[edgeNum][3+i][3+j] += forceCreaseVal*block_jj[i][j] + dForcecrease*tensorProduct(dTheta_xj,dTheta_xj)[i][j];
                    edgeStiffness[edgeNum][3+i][6+j] += forceCreaseVal*block_jk[i][j] + dForcecrease*tensorProduct(dTheta_xj,dTheta_xk)[i][j];
                    edgeStiffness[edgeNum][6+i][3+j] += forceCreaseVal*block_kj[i][j] + dForcecrease*tensorProduct(dTheta_xk,dTheta_xj)[i][j];
                    edgeStiffness[edgeNum][3+i][9+j] += forceCreaseVal*block_jl[i][j] + dForcecrease*tensorProduct(dTheta_xj,dTheta_xl)[i][j];
                    edgeStiffness[edgeNum][9+i][3+j] += forceCreaseVal*block_lj[i][j] + dForcecrease*tensorProduct(dTheta_xl,dTheta_xj)[i][j];
                    edgeStiffness[edgeNum][6+i][6+j] += forceCreaseVal*block_kk[i][j] + dForcecrease*tensorProduct(dTheta_xk,dTheta_xk)[i][j];
                    edgeStiffness[edgeNum][6+i][9+j] += forceCreaseVal*block_kl[i][j] + dForcecrease*tensorProduct(dTheta_xk,dTheta_xl)[i][j];
                    edgeStiffness[edgeNum][9+i][6+j] += forceCreaseVal*block_lk[i][j] + dForcecrease*tensorProduct(dTheta_xl,dTheta_xk)[i][j];
                    edgeStiffness[edgeNum][9+i][9+j] += forceCreaseVal*block_ll[i][j] + dForcecrease*tensorProduct(dTheta_xl,dTheta_xl)[i][j];
                }

            }
            msg_warning() << edgeStiffness[edgeNum];
            msg_warning() << "----------------------------------";
        }
//        msg_warning() << "Edge Stiffness Set ! " << edgeNum;
//        msg_warning() <<  edgeStiffness[edgeNum];

    }
//    if(method==SMALL)
//    {
//        typename VecElement::const_iterator it;
//        unsigned int i(0);

//        for(it = _indexedElements->begin() ; it != _indexedElements->end() ; ++it, ++i)
//        {
//            accumulateForceSmall( f1, x1, i, true );
//        }
//    }
//    else
//    {
//        typename VecElement::const_iterator it;
//        unsigned int i(0);

//        for(it = _indexedElements->begin() ; it != _indexedElements->end() ; ++it, ++i)
//        {
//            accumulateForceLarge( f1, x1, i, true );
//        }
//    }

    f.endEdit();
}


template <class DataTypes>
void RotationalSpringForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx)
{
    msg_warning() << "In addDDDDDDDDDDDDforce !!!";
    VecDeriv& df1 = *df.beginEdit();
    const VecDeriv& dx1 = dx.getValue();
    Real kFactor = (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(mparams, this->rayleighStiffness.getValue());

//    Real h=1;
    df1.resize(dx1.size());
    for (int edgeNum=0; edgeNum < m_topology->getNbEdges(); edgeNum++)
    {
        auto edgesAssignment = d_edgesAssignment.getValue();
        if (edgesAssignment[edgeNum] == "M" || edgesAssignment[edgeNum] == "V"){
            for (int i=0;i<4;i++){
                for (int j=0;j<12;j++){
//                    msg_warning() << "i: " << i << " j: " << j;
//                    msg_warning() << "edgeNum: "<< edgeNum << " nodes[edgeNum]: " << nodes[edgeNum] << " floor(j/3): " << floor(j/3) << " j%3: " << j%3;
                    df1[nodes[edgeNum][i]][0] += edgeStiffness[edgeNum][3*i][j]*dx1[nodes[edgeNum][floor(j/3)]][j%3] * kFactor;
                    df1[nodes[edgeNum][i]][1] += edgeStiffness[edgeNum][3*i+1][j]*dx1[nodes[edgeNum][floor(j/3)]][j%3]* kFactor;
                    df1[nodes[edgeNum][i]][2] += edgeStiffness[edgeNum][3*i+2][j]*dx1[nodes[edgeNum][floor(j/3)]][j%3]* kFactor;
                }
            }
        }
    }
//    if (method == SMALL)
//    {
//        applyStiffnessSmall( df1, h, dx1, kFactor );
//    }
//    else
//    {
//        applyStiffnessLarge( df1, h, dx1, kFactor );
//    }

//    df.endEdit();
    df.endEdit();
}

template <class DataTypes>
void RotationalSpringForceField<DataTypes>::applyStiffness( VecCoord& v, Real h, const VecCoord& x, const SReal &kFactor )
{
    if (method == SMALL)
    {
        applyStiffnessSmall( v, h, x, kFactor );
    }
    else
    {
        applyStiffnessLarge( v, h, x, kFactor );
    }
}


template <class DataTypes>
void RotationalSpringForceField<DataTypes>::computeStrainDisplacement( StrainDisplacement &J, Coord /*a*/, Coord b, Coord c )
{
    //    //Coord ab_cross_ac = cross(b, c);
    Real determinant = b[0] * c[1]; // Surface * 2

    J[0][0] = J[1][2] = -c[1] / determinant;
    J[0][2] = J[1][1] = (c[0] - b[0]) / determinant;
    J[2][0] = J[3][2] = c[1] / determinant;
    J[2][2] = J[3][1] = -c[0] / determinant;
    J[4][0] = J[5][2] = 0;
    J[4][2] = J[5][1] = b[0] / determinant;
    J[1][0] = J[3][0] = J[5][0] = J[0][1] = J[2][1] = J[4][1] = 0;

    /* The following formulation is actually equivalent:
      Let
      | alpha1 alpha2 alpha3 |                      | 1 xa ya |
      | beta1  beta2  beta3  | = be the inverse of  | 1 xb yb |
      | gamma1 gamma2 gamma3 |                      | 1 xc yc |
      The strain-displacement matrix is:
      | beta1  0       beta2  0        beta3  0      |
      | 0      gamma1  0      gamma2   0      gamma3 | / (2*A)
      | gamma1 beta1   gamma2 beta2    gamma3 beta3  |
      where A is the area of the triangle and 2*A is the determinant of the matrix with the xa,ya,xb...
      Since a0=a1=b1=0, the matrix is triangular and its inverse is:
      |  1              0              0  |
      | -1/xb           1/xb           0  |
      | -(1-xc/xb)/yc  -xc/(xb*yc)   1/yc |
      our strain-displacement matrix is:
      | -1/xb           0             1/xb         0            0     0    |
      | 0              -(1-xc/xb)/yc  0            -xc/(xb*yc)  0     1/yc |
      | -(1-xc/xb)/yc  -1/xb          -xc/(xb*yc)  1/xb         1/yc  0    |
      */

    //    Real beta1  = -1/b[0];
    //    Real beta2  =  1/b[0];
    //    Real gamma1 = (c[0]/b[0]-1)/c[1];
    //    Real gamma2 = -c[0]/(b[0]*c[1]);
    //    Real gamma3 = 1/c[1];

    //    // The transpose of the strain-displacement matrix is thus:
    //    J[0][0] = J[1][2] = beta1;
    //    J[0][1] = J[1][0] = 0;
    //    J[0][2] = J[1][1] = gamma1;

    //    J[2][0] = J[3][2] = beta2;
    //    J[2][1] = J[3][0] = 0;
    //    J[2][2] = J[3][1] = gamma2;

    //    J[4][0] = J[5][2] = 0;
    //    J[4][1] = J[5][0] = 0;
    //    J[4][2] = J[5][1] = gamma3;



}


template <class DataTypes>
void RotationalSpringForceField<DataTypes>::computeMaterialStiffnesses()
{
    _materialsStiffnesses.resize(_indexedElements->size());
    const VecCoord& p= _initialPoints.getValue();


    for(unsigned i = 0; i < _indexedElements->size(); ++i)
    {
        Index a = (*_indexedElements)[i][0];
        Index b = (*_indexedElements)[i][1];
        Index c = (*_indexedElements)[i][2];
        Real triangleVolume = (Real)0.5 * f_thickness.getValue() * cross( p[b]-p[a], p[c]-p[a] ).norm();

        if( f_planeStrain.getValue() == true )
        {
            _materialsStiffnesses[i][0][0] = 1-f_poisson.getValue();
            _materialsStiffnesses[i][0][1] = f_poisson.getValue();
            _materialsStiffnesses[i][0][2] = 0;
            _materialsStiffnesses[i][1][0] = f_poisson.getValue();
            _materialsStiffnesses[i][1][1] = 1-f_poisson.getValue();
            _materialsStiffnesses[i][1][2] = 0;
            _materialsStiffnesses[i][2][0] = 0;
            _materialsStiffnesses[i][2][1] = 0;
            _materialsStiffnesses[i][2][2] = 0.5f - f_poisson.getValue();

            _materialsStiffnesses[i] *= f_young.getValue() / ( (1 + f_poisson.getValue()) * (1-2*f_poisson.getValue()) ) * triangleVolume;
        }
        else // plane stress
        {
            _materialsStiffnesses[i][0][0] = 1;
            _materialsStiffnesses[i][0][1] = f_poisson.getValue();
            _materialsStiffnesses[i][0][2] = 0;
            _materialsStiffnesses[i][1][0] = f_poisson.getValue();
            _materialsStiffnesses[i][1][1] = 1;
            _materialsStiffnesses[i][1][2] = 0;
            _materialsStiffnesses[i][2][0] = 0;
            _materialsStiffnesses[i][2][1] = 0;
            _materialsStiffnesses[i][2][2] = 0.5f * (1 - f_poisson.getValue());

            _materialsStiffnesses[i] *= f_young.getValue() / ( (1 - f_poisson.getValue() * f_poisson.getValue())) * triangleVolume;
        }
    }
}


template <class DataTypes>
void RotationalSpringForceField<DataTypes>::computeForce( Displacement &F, const Displacement &Depl, const MaterialStiffness &K, const StrainDisplacement &J )
{
    type::Mat<3,6,Real> Jt;
    Jt.transpose( J );

    type::Vec<3,Real> JtD;

    // Optimisations: The following values are 0 (per computeStrainDisplacement )

    // Jt[0][1]
    // Jt[0][3]
    // Jt[0][4]
    // Jt[0][5]
    // Jt[1][0]
    // Jt[1][2]
    // Jt[1][4]
    // Jt[2][5]


    //	JtD = Jt * Depl;

    JtD[0] = Jt[0][0] * Depl[0] + /* Jt[0][1] * Depl[1] + */ Jt[0][2] * Depl[2]
            /* + Jt[0][3] * Depl[3] + Jt[0][4] * Depl[4] + Jt[0][5] * Depl[5] */ ;

    JtD[1] = /* Jt[1][0] * Depl[0] + */ Jt[1][1] * Depl[1] + /* Jt[1][2] * Depl[2] + */
            Jt[1][3] * Depl[3] + /* Jt[1][4] * Depl[4] + */ Jt[1][5] * Depl[5];

    JtD[2] = Jt[2][0] * Depl[0] + Jt[2][1] * Depl[1] + Jt[2][2] * Depl[2] +
            Jt[2][3] * Depl[3] + Jt[2][4] * Depl[4] /* + Jt[2][5] * Depl[5] */ ;

    type::Vec<3,Real> KJtD;

    //	KJtD = K * JtD;

    // Optimisations: The following values are 0 (per computeMaterialStiffnesses )

    // K[0][2]
    // K[1][2]
    // K[2][0]
    // K[2][1]

    KJtD[0] = K[0][0] * JtD[0] + K[0][1] * JtD[1] /* + K[0][2] * JtD[2] */;

    KJtD[1] = K[1][0] * JtD[0] + K[1][1] * JtD[1] /* + K[1][2] * JtD[2] */;

    KJtD[2] = /* K[2][0] * JtD[0] + K[2][1] * JtD[1] */ + K[2][2] * JtD[2];

    //	F = J * KJtD;


    // Optimisations: The following values are 0 (per computeStrainDisplacement )

    // J[0][1]
    // J[1][0]
    // J[2][1]
    // J[3][0]
    // J[4][0]
    // J[4][1]
    // J[5][0]
    // J[5][2]

    F[0] = J[0][0] * KJtD[0] + /* J[0][1] * KJtD[1] + */ J[0][2] * KJtD[2];

    F[1] = /* J[1][0] * KJtD[0] + */ J[1][1] * KJtD[1] + J[1][2] * KJtD[2];

    F[2] = J[2][0] * KJtD[0] + /* J[2][1] * KJtD[1] + */ J[2][2] * KJtD[2];

    F[3] = /* J[3][0] * KJtD[0] + */ J[3][1] * KJtD[1] + J[3][2] * KJtD[2];

    F[4] = /* J[4][0] * KJtD[0] + J[4][1] * KJtD[1] + */ J[4][2] * KJtD[2];

    F[5] = /* J[5][0] * KJtD[0] + */ J[5][1] * KJtD[1] /* + J[5][2] * KJtD[2] */ ;
}


/*
** SMALL DEFORMATION METHODS
*/
template <class DataTypes>
void RotationalSpringForceField<DataTypes>::initSmall()
{
    Transformation identity;
    identity[0][0]=identity[1][1]=identity[2][2]=1;
    identity[0][1]=identity[0][2]=0;
    identity[1][0]=identity[1][2]=0;
    identity[2][0]=identity[2][1]=0;
    for(unsigned i=0; i< _indexedElements->size() ; ++i)
    {
        _rotations[i] = identity;
    }
}


template <class DataTypes>
void RotationalSpringForceField<DataTypes>::accumulateForceSmall( VecCoord &f, const VecCoord &p, Index elementIndex, bool implicit )
{
    Index a = (*_indexedElements)[elementIndex][0];
    Index b = (*_indexedElements)[elementIndex][1];
    Index c = (*_indexedElements)[elementIndex][2];

    Coord deforme_a, deforme_b, deforme_c;
    deforme_b = p[b]-p[a];
    deforme_c = p[c]-p[a];
    deforme_a = Coord(0,0,0);

    // displacements
    Displacement D;
    D[0] = 0;
    D[1] = 0;
    D[2] = (_initialPoints.getValue()[b][0]-_initialPoints.getValue()[a][0]) - deforme_b[0];
    D[3] = 0;
    D[4] = (_initialPoints.getValue()[c][0]-_initialPoints.getValue()[a][0]) - deforme_c[0];
    D[5] = (_initialPoints.getValue()[c][1]-_initialPoints.getValue()[a][1]) - deforme_c[1];


    StrainDisplacement J;
    computeStrainDisplacement(J,deforme_a,deforme_b,deforme_c);
    if (implicit)
        _strainDisplacements[elementIndex] = J;

    // compute force on element
    Displacement F;
    computeForce( F, D, _materialsStiffnesses[elementIndex], J );

    f[a] += Coord( F[0], F[1], 0);
    f[b] += Coord( F[2], F[3], 0);
    f[c] += Coord( F[4], F[5], 0);
}

template <class DataTypes>
void RotationalSpringForceField<DataTypes>::applyStiffnessSmall(VecCoord &v, Real h, const VecCoord &x, const SReal &kFactor)
{
    typename VecElement::const_iterator it;
    unsigned int i(0);

    for(it = _indexedElements->begin() ; it != _indexedElements->end() ; ++it, ++i)
    {
        Index a = (*it)[0];
        Index b = (*it)[1];
        Index c = (*it)[2];

        Displacement X;

        X[0] = x[a][0];
        X[1] = x[a][1];

        X[2] = x[b][0];
        X[3] = x[b][1];

        X[4] = x[c][0];
        X[5] = x[c][1];

        Displacement F;
        computeForce( F, X, _materialsStiffnesses[i], _strainDisplacements[i] );

        v[a] += Coord(-h*F[0], -h*F[1], 0)*kFactor;
        v[b] += Coord(-h*F[2], -h*F[3], 0)*kFactor;
        v[c] += Coord(-h*F[4], -h*F[5], 0)*kFactor;
    }
}


/*
** LARGE DEFORMATION METHODS
*/

template <class DataTypes>
void RotationalSpringForceField<DataTypes>::initLarge()
{
    _rotatedInitialElements.resize(_indexedElements->size());

    typename VecElement::const_iterator it;
    unsigned int i(0);

    for(it = _indexedElements->begin() ; it != _indexedElements->end() ; ++it, ++i)
    {
        Index a = (*it)[0];
        Index b = (*it)[1];
        Index c = (*it)[2];

        // Rotation matrix (transpose of initial triangle/world)
        // first vector on first edge
        // second vector in the plane of the two first edges
        // third vector orthogonal to first and second
        Transformation R_0_1;
        computeRotationLarge( R_0_1, _initialPoints.getValue(), a, b, c );
        _rotations[i].transpose(R_0_1);

        // coordinates of the triangle vertices in their local frames
        _rotatedInitialElements[i][0] = R_0_1 * _initialPoints.getValue()[a];
        _rotatedInitialElements[i][1] = R_0_1 * _initialPoints.getValue()[b];
        _rotatedInitialElements[i][2] = R_0_1 * _initialPoints.getValue()[c];
        // set the origin of the local frame at vertex a
        _rotatedInitialElements[i][1] -= _rotatedInitialElements[i][0];
        _rotatedInitialElements[i][2] -= _rotatedInitialElements[i][0];
        _rotatedInitialElements[i][0] = Coord(0,0,0);

        computeStrainDisplacement(_strainDisplacements[i], _initialPoints.getValue()[a], _initialPoints.getValue()[b], _initialPoints.getValue()[c] );
    }
}


template <class DataTypes>
void RotationalSpringForceField<DataTypes>::computeRotationLarge( Transformation &r, const VecCoord &p, const Index &a, const Index &b, const Index &c)
{
    // first vector on first edge
    // second vector in the plane of the two first edges
    // third vector orthogonal to first and second

    Coord edgex = p[b] - p[a];
    edgex.normalize();

    Coord edgey = p[c] - p[a];
    edgey.normalize();

    Coord edgez;
    edgez = cross(edgex, edgey);
    edgez.normalize();

    edgey = cross(edgez, edgex);
    edgey.normalize();

    r[0][0] = edgex[0];
    r[0][1] = edgex[1];
    r[0][2] = edgex[2];
    r[1][0] = edgey[0];
    r[1][1] = edgey[1];
    r[1][2] = edgey[2];
    r[2][0] = edgez[0];
    r[2][1] = edgez[1];
    r[2][2] = edgez[2];
}


template <class DataTypes>
void RotationalSpringForceField<DataTypes>::accumulateForceLarge(VecCoord &f, const VecCoord &p, Index elementIndex, bool implicit )
{
    // triangle vertex indices
    Index a = (*_indexedElements)[elementIndex][0];
    Index b = (*_indexedElements)[elementIndex][1];
    Index c = (*_indexedElements)[elementIndex][2];

    // Rotation matrix (deformed and displaced Triangle/world)
    Transformation R_2_0, R_0_2;
    computeRotationLarge( R_0_2, p, a, b, c);


    // positions of the deformed points in the local frame
    Coord deforme_a, deforme_b, deforme_c;
    deforme_b = R_0_2 * (p[b]-p[a]);
    deforme_c = R_0_2 * (p[c]-p[a]);

    // displacements in the local frame
    Displacement D;
    D[0] = 0;
    D[1] = 0;
    D[2] = _rotatedInitialElements[elementIndex][1][0] - deforme_b[0];
    D[3] = 0;
    D[4] = _rotatedInitialElements[elementIndex][2][0] - deforme_c[0];
    D[5] = _rotatedInitialElements[elementIndex][2][1] - deforme_c[1];

    // Strain-displacement matrix
    StrainDisplacement B;
    computeStrainDisplacement(B,deforme_a,deforme_b,deforme_c);

    // compute force on element, in local frame
    Displacement F;
    computeForce( F, D, _materialsStiffnesses[elementIndex], B ); // F = Bt.S.B.D

    // project forces to world frame
    R_2_0.transpose(R_0_2);
    f[a] += R_2_0 * Coord(F[0], F[1], 0);
    f[b] += R_2_0 * Coord(F[2], F[3], 0);
    f[c] += R_2_0 * Coord(F[4], F[5], 0);

    // store for re-use in matrix-vector products
    if(implicit)
    {
        _strainDisplacements[elementIndex] = B;
        _rotations[elementIndex] = R_2_0 ;
    }

}


template <class DataTypes>
void RotationalSpringForceField<DataTypes>::applyStiffnessLarge(VecCoord &v, Real h, const VecCoord &x, const SReal &kFactor)
{
    typename VecElement::const_iterator it;
    unsigned int i(0);

    for(it = _indexedElements->begin() ; it != _indexedElements->end() ; ++it, ++i)
    {
        Index a = (*it)[0];
        Index b = (*it)[1];
        Index c = (*it)[2];

        Transformation R_0_2;
        R_0_2.transpose(_rotations[i]);

        Displacement X;
        Coord x_2;

        x_2 = R_0_2 * x[a];
        X[0] = x_2[0];
        X[1] = x_2[1];

        x_2 = R_0_2 * x[b];
        X[2] = x_2[0];
        X[3] = x_2[1];

        x_2 = R_0_2 * x[c];
        X[4] = x_2[0];
        X[5] = x_2[1];

        Displacement F;
        computeForce( F, X, _materialsStiffnesses[i], _strainDisplacements[i] );

        v[a] += (_rotations[i] * Coord(-h*F[0], -h*F[1], 0))*kFactor;
        v[b] += (_rotations[i] * Coord(-h*F[2], -h*F[3], 0))*kFactor;
        v[c] += (_rotations[i] * Coord(-h*F[4], -h*F[5], 0))*kFactor;
    }
}


template<class DataTypes>
void RotationalSpringForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowForceFields())
        return;

    vparams->drawTool()->saveLastState();
    vparams->drawTool()->disableLighting();

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0, true);

    std::vector<sofa::type::RGBAColor> colorVector;
    std::vector<sofa::type::Vector3> vertices;

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();

    typename VecElement::const_iterator it;
    for(it = _indexedElements->begin() ; it != _indexedElements->end() ; ++it)
    {
        Index a = (*it)[0];
        Index b = (*it)[1];
        Index c = (*it)[2];

        colorVector.push_back(sofa::type::RGBAColor(0,1,0,1));
        vertices.push_back(sofa::type::Vector3(x[a]));
        colorVector.push_back(sofa::type::RGBAColor(0,0.5,0.5,1));
        vertices.push_back(sofa::type::Vector3(x[b]));
        colorVector.push_back(sofa::type::RGBAColor(0,0,1,1));
        vertices.push_back(sofa::type::Vector3(x[c]));
    }
    vparams->drawTool()->drawTriangles(vertices,colorVector);

    vparams->drawTool()->restoreLastState();
}


template<class DataTypes>
void RotationalSpringForceField<DataTypes>::computeElementStiffnessMatrix( StiffnessMatrix& S, StiffnessMatrix& SR, const MaterialStiffness &K, const StrainDisplacement &J, const Transformation& Rot )
{
    type::MatNoInit<3, 6, Real> Jt;
    Jt.transpose( J );

    type::MatNoInit<6, 6, Real> JKJt;
    JKJt = J*K*Jt;  // in-plane stiffness matrix, 6x6

    // stiffness JKJt expanded to 3 dimensions
    type::Mat<9, 9, Real> Ke; // initialized to 0
    // for each 2x2 block i,j
    for(unsigned i=0; i<3; i++)
    {
        for(unsigned j=0; j<3; j++)
        {
            // copy the block in the expanded matrix
            for(unsigned k=0; k<2; k++)
                for(unsigned l=0; l<2; l++)
                    Ke[3*i+k][3*j+l] = JKJt[2*i+k][2*j+l];
        }
    }

    // rotation matrices. TODO: use block-diagonal matrices, more efficient.
    type::Mat<9, 9, Real> RR,RRt; // initialized to 0
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
        {
            RR[i][j]=RR[i+3][j+3]=RR[i+6][j+6]=Rot[i][j];
            RRt[i][j]=RRt[i+3][j+3]=RRt[i+6][j+6]=Rot[j][i];
        }

    S = RR*Ke;
    SR = S*RRt;
}


template<class DataTypes>
void RotationalSpringForceField<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal k, unsigned int &offset)
{
//    for(unsigned i=0; i< _indexedElements->size() ; i++)
//    {
//        StiffnessMatrix JKJt,RJKJtRt;
//        computeElementStiffnessMatrix(JKJt, RJKJtRt, _materialsStiffnesses[i], _strainDisplacements[i], _rotations[i]);
//        this->addToMatrix(mat,offset,(*_indexedElements)[i],RJKJtRt,-k);
//    }
    for (int edgeNum=0; edgeNum < m_topology->getNbEdges(); edgeNum++)
    {

        this->addToMatrix(mat,offset,nodes[edgeNum],edgeStiffness[edgeNum],-k);

    }
    msg_warning() << "mat : " << mat;
}

template<class DataTypes>
void RotationalSpringForceField<DataTypes>::setPoisson(Real val)
{
    if (val < 0)
    {
        msg_warning() << "Input Poisson Coefficient is not possible: " << val << ", setting default value: 0.3";
        f_poisson.setValue(0.3);
    }
    else if (val != f_poisson.getValue())
    {
        f_poisson.setValue(val);
    }
}

template<class DataTypes>
void RotationalSpringForceField<DataTypes>::setYoung(Real val)
{
    if (val < 0)
    {
        msg_warning() << "Input Young Modulus is not possible: " << val << ", setting default value: 1000";
        f_young.setValue(1000);
    }
    else if (val != f_young.getValue())
    {
        f_young.setValue(val);
    }
}

template<class DataTypes>
void RotationalSpringForceField<DataTypes>::setMethod(int val)
{
    if (val != 0 && val != 1)
    {
        msg_warning() << "Input Method is not possible: " << val << ", should be 0 (Large) or 1 (Small). Setting default value: Large";
        method = LARGE;
    }
    else if (method != val)
    {
        method = val;
    }
}

template<class DataTypes>
void RotationalSpringForceField<DataTypes>::setMethod(std::string val)
{
    if (val == "small")
        method = SMALL;
    else if (val == "large")
        method = LARGE;
    else
    {
        msg_warning() << "Input Method is not possible: " << val << ", should be 0 (Large) or 1 (Small). Setting default value: Large";
        method = LARGE;
    }
}


template<class DataTypes>
const type::fixed_array <typename RotationalSpringForceField<DataTypes>::Coord, 3>& RotationalSpringForceField<DataTypes>::getRotatedInitialElement(Index elemId)
{
    if (elemId != sofa::InvalidID && elemId < _rotatedInitialElements.size())
        return _rotatedInitialElements[elemId];

    msg_warning() << "Method getRotatedInitialElement called with element index: " << elemId
        << " which is out of bounds: [0, " << _rotatedInitialElements.size() << "]. Returning default empty array of coordinates.";
    return InvalidCoords;
}

template<class DataTypes>
const typename RotationalSpringForceField<DataTypes>::Transformation& RotationalSpringForceField<DataTypes>::getRotationMatrix(Index elemId)
{
    if (elemId != sofa::InvalidID && elemId < _rotations.size())
        return _rotations[elemId];

    msg_warning() << "Method getRotationMatrix called with element index: "
        << elemId << " which is out of bounds: [0, " << _rotations.size() << "]. Returning default empty rotation.";
    return InvalidTransform;
}

template<class DataTypes>
const typename RotationalSpringForceField<DataTypes>::MaterialStiffness& RotationalSpringForceField<DataTypes>::getMaterialStiffness(Index elemId)
{
    if (elemId != sofa::InvalidID && elemId < _materialsStiffnesses.size())
        return _materialsStiffnesses[elemId];

    msg_warning() << "Method getMaterialStiffness called with element index: "
        << elemId << " which is out of bounds: [0, " << _materialsStiffnesses.size() << "]. Returning default empty matrix.";
    return InvalidTransform;
}

template<class DataTypes>
const typename RotationalSpringForceField<DataTypes>::StrainDisplacement& RotationalSpringForceField<DataTypes>::getStrainDisplacements(Index elemId)
{
    if (elemId != sofa::InvalidID && elemId < _strainDisplacements.size())
        return _strainDisplacements[elemId];

    msg_warning() << "Method getStrainDisplacements called with element index: "
        << elemId << " which is out of bounds: [0, " << _strainDisplacements.size() << "]. Returning default empty displacements.";
    return InvalidStrainDisplacement;
}



} // namespace sofa::component::forcefield
