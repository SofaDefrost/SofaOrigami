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
    , f_kcrease(initData(&f_kcrease,Real(1000.),"kcrease","Crease stiffness"))
    , f_thickness(initData(&f_thickness,Real(1.),"thickness","Thickness of the elements"))
//    , f_damping(initData(&f_damping,(Real)0.,"damping","Ratio damping/stiffness"))
    , f_planeStrain(initData(&f_planeStrain,false,"planeStrain","Plane strain or plane stress assumption"))
    , l_topology(initLink("topology", "link to the topology container"))
{
    f_poisson.setRequired(true);
    f_kcrease.setRequired(true);
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
    setYoung(f_kcrease.getValue());

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

//    msg_warning() << "Get triangles tesssssssssssssssssssssssssssssssssssssssssssssst";
//    msg_warning() << "Nb edges: " << m_topology->getNbEdges();
    adjacentVertices.resize(m_topology->getNbEdges());
    edgeInitialLength.resize(m_topology->getNbEdges());
    edgeStiffness.resize(m_topology->getNbEdges());
    nodes.resize(m_topology->getNbEdges());
    forceCrease.resize(m_topology->getNbEdges());
    roundCount.resize(m_topology->getNbEdges());
    signTemp.resize(m_topology->getNbEdges());
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
        signTemp[i] = 0;
        roundCount[i] = 0;
    }

//for (int i; i<m_topology->getNbEdges(); i++)
//    msg_warning() << "edge:" << i << " adj vert: " << adjacentVertices[i][0] << " " << adjacentVertices[i][1];
}

template <class DataTypes>
type::Mat3x3 RotationalSpringForceField<DataTypes>::diamond(Coord a, Coord b){
    return tensorProduct(a,b) + tensorProduct(b,a);
}

template <class DataTypes>
void RotationalSpringForceField<DataTypes>::getRotationalStiffness(const Real theta,const int edgeNum, Real& M, Real& k){

    auto edgesAssignment = d_edgesAssignment.getValue();
    Real kcrease = f_kcrease.getValue();
    bool linear = true;
    Real theta1 = M_PI/2.0 ,theta0 = d_angleTarget.getValue(),theta2 = 3.0*M_PI/2.0;

    if (edgesAssignment[edgeNum] == "V"){
        if (linear){
            M = -kcrease*edgeInitialLength[edgeNum]*(theta - theta0); // Resistance Moment
            k = -kcrease*edgeInitialLength[edgeNum]; // tangent rotational stiffness = dM/dtheta
        }
        else
        {
            if (theta <= theta1){
                Real phi = M_PI*(theta - theta1)/(2*theta1);
                M = -kcrease*edgeInitialLength[edgeNum]*(theta1 - theta0) + (2*kcrease*theta1/M_PI) * tan(phi); // Resistance Moment
                k = -kcrease*edgeInitialLength[edgeNum]/(cos(phi)*cos(phi));
            }
            else if (theta <= theta2){
                M = -kcrease*edgeInitialLength[edgeNum]*(theta - theta0);
                k = -kcrease*edgeInitialLength[edgeNum];
            }
            else  // theta >theta2
            {
                Real phi = M_PI*(theta - theta2)/(4*M_PI - 2*theta2);
                M = -kcrease*edgeInitialLength[edgeNum]*(theta2 - theta0) + (2*kcrease*(2*M_PI - theta2)/M_PI) * tan(phi); // Resistance Moment
                k = -kcrease*edgeInitialLength[edgeNum]/(cos(phi)*cos(phi));
            }

        }

    } else if (edgesAssignment[edgeNum] == "M"){
        if (linear){
            M = -kcrease*edgeInitialLength[edgeNum]*(theta - (2*M_PI - theta0)); // Resistance Moment
            k = -kcrease*edgeInitialLength[edgeNum]; // tangent rotational stiffness = dM/dtheta
        }
        else
        {
            if (theta <= theta1){
                Real phi = M_PI*(theta - (2*M_PI - theta1))/(2*(2*M_PI - theta1));
                M = -kcrease*edgeInitialLength[edgeNum]*(theta1 - theta0) + (2*kcrease*theta1/M_PI) * tan(phi); // Resistance Moment
                k = -kcrease*edgeInitialLength[edgeNum]/(cos(phi)*cos(phi));
            }
            else if (theta <= theta2){
                M = -kcrease*edgeInitialLength[edgeNum]*(theta - (2*M_PI - theta0));
                k = -kcrease*edgeInitialLength[edgeNum];
            }
            else  // theta >theta2
            {
                Real phi = M_PI*(theta - (2*M_PI - theta2))/(4*M_PI - 2*(2*M_PI - theta2));
                M = -kcrease*edgeInitialLength[edgeNum]*(theta2 - theta0) + (2*kcrease*(2*M_PI - theta2)/M_PI) * tan(phi); // Resistance Moment
                k = -kcrease*edgeInitialLength[edgeNum]/(cos(phi)*cos(phi));
            }

        }
    } else {
        M = 0.0;
        k = 0.0;
    }

}

template <class DataTypes>
void RotationalSpringForceField<DataTypes>::addForce(const core::MechanicalParams* /* mparams */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& /* v */)
{
    VecDeriv& f1 = *f.beginEdit();
    const VecCoord& positions = x.getValue();
//    msg_warning() << "--------------------------------------------------------------------------";
//    std::cout << "[" << positions[0][0] << ", " << positions[0][1] << ", " << positions[0][2] << "]";
//    std::cout << "[" << positions[13][0] << ", " << positions[13][1] << ", " << positions[13][2] << "]";
//    std::cout << "[" << positions[23][0] << ", " << positions[23][1] << ", " << positions[23][2] << "]";
//    std::cout << "[" << positions[27][0] << ", " << positions[27][1] << ", " << positions[27][2] << "]";
//    std::cout << "[" << positions[31][0] << ", " << positions[31][1] << ", " << positions[31][2] << "]";
//    std::cout << "[" << positions[35][0] << ", " << positions[35][1] << ", " << positions[35][2] << "]";
//    std::cout << "[" << positions[37][0] << ", " << positions[37][1] << ", " << positions[37][2] << "]";
//    std::cout << "[" << positions[41][0] << ", " << positions[41][1] << ", " << positions[41][2] << "]";
//    std::cout << "[" << positions[45][0] << ", " << positions[45][1] << ", " << positions[45][2] << "]";
//    std::cout << "[" << positions[49][0] << ", " << positions[49][1] << ", " << positions[49][2] << "]";
//    std::cout << "[" << positions[53][0] << ", " << positions[53][1] << ", " << positions[53][2] << "]\n" ;

//    msg_warning() << "--------------------------------------------------------------------------";
    const VecCoord& x1 = x.getValue();
    f1.resize(positions.size());

    Real kcrease = f_kcrease.getValue();
    Real forceCreaseVal, dForcecrease;
    auto edgesAssignment = d_edgesAssignment.getValue();
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
//        msg_warning() << "edge Number : " <<edgeNum;
//        msg_warning() << "Nodes : " << nodes[edgeNum];
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
//        msg_warning() << "¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨";
//        if (edgesAssignment[edgeNum] == "M")
//            msg_warning() << "theta Mountain original: " << theta << "angle target : " << 2*M_PI - d_angleTarget.getValue();
//        else if (edgesAssignment[edgeNum] == "V")
//            msg_warning() << "theta Valley original: " << theta << "angle target : " << d_angleTarget.getValue();

//        msg_warning() << "cos theta : !!!!!!!!!!!!!!!! " << dot(m,n)/(norm(m)*norm(n));
        int signMdotRkl = 1;
        if (dot(m,rkl) != 0)
            signMdotRkl = (0 < dot(m,rkl)) - (dot(m,rkl)<0);
        bool fromNeg2Pos = false,fromPos2Neg=false;
        if ((signMdotRkl*signTemp[edgeNum] < 0) && (abs(theta - M_PI) > 0.01))
            if (signMdotRkl > 0)
                fromNeg2Pos = true;
            else
                fromPos2Neg = true;

        signTemp[edgeNum] = signMdotRkl;
        theta = signMdotRkl * theta;
        if (fromNeg2Pos)
            roundCount[edgeNum] += 1;
        else if (fromPos2Neg)
            roundCount[edgeNum] -= 1;
//        msg_warning() << "Updated roundCount[" << edgeNum << "]: " << roundCount[edgeNum];
//        if (edgesAssignment[edgeNum] == "M")
//            msg_warning() << "theta Mountain sign change: " << theta << "angle target : " << 2*M_PI - d_angleTarget.getValue();
//        else if (edgesAssignment[edgeNum] == "V")
//            msg_warning() << "theta Valley sign change: " << theta << "angle target : " << d_angleTarget.getValue();


        theta = theta - floor( theta/(2*M_PI))*2*M_PI; // Modulo 2PI
        theta += roundCount[edgeNum]*2*M_PI;
        Deriv dTheta_xi = (norm(rkj)/dot(m,m))*m;
        Deriv dTheta_xl = -(norm(rkj)/dot(n,n))*n;
        Deriv dTheta_xj = (dot(rij,rkj)/dot(rkj,rkj) - 1)*dTheta_xi - (dot(rkl,rkj)/dot(rkj,rkj))*dTheta_xl;
        Deriv dTheta_xk = (dot(rkl,rkj)/dot(rkj,rkj) - 1)*dTheta_xl - (dot(rij,rkj)/dot(rkj,rkj))*dTheta_xi;
//        msg_warning() << "dot(rij,rkj): " << dot(rij,rkj);
//        msg_warning() << "dot(rkl,rkj): " << dot(rkl,rkj);

//        msg_warning() << "forces : !!!!!!!!!!!!!!!! " << fxi;
//        msg_warning() << "forces : !!!!!!!!!!!!!!!! " << fxl;
//        msg_warning() << "forces : !!!!!!!!!!!!!!!! " << fxj;
//        msg_warning() << "forces : !!!!!!!!!!!!!!!! " << fxk;
//        if (edgesAssignment[edgeNum] == "M")
//            msg_warning() << "theta Mountain: " << theta << "angle target : " << 2*M_PI - d_angleTarget.getValue();
//        else if (edgesAssignment[edgeNum] == "V")
//            msg_warning() << "theta Valley : " << theta << "angle target : " << d_angleTarget.getValue();

        getRotationalStiffness(theta, edgeNum, forceCreaseVal, dForcecrease);


//        if (edgesAssignment[edgeNum] == "V"){
//            forceCreaseVal = -kcrease*edgeInitialLength[edgeNum]*(theta - d_angleTarget.getValue());
////            msg_warning() << "force detail:" << kcrease << " " << edgeInitialLength[edgeNum] << " " << (theta - d_angleTarget.getValue());
//            dForcecrease = -kcrease*edgeInitialLength[edgeNum];
//        } else if (edgesAssignment[edgeNum] == "M"){
//            forceCreaseVal = -kcrease*edgeInitialLength[edgeNum]*(theta - (2*M_PI - d_angleTarget.getValue()));
//            dForcecrease = -kcrease*edgeInitialLength[edgeNum];
//        } else {
//            forceCreaseVal = 0.0;
//            dForcecrease = 0.0;
//        }
//        dForcecrease = 0.0;
        forceCrease[edgeNum] = forceCreaseVal;
//        msg_warning() << "forceCreaseVal !!!!!!!!!!!!!!!!" << forceCreaseVal;
        if ((edgesAssignment[edgeNum] == "M") || (edgesAssignment[edgeNum] == "V")) {
            f1[adjacent[0]] += forceCreaseVal*dTheta_xi;
            f1[adjacent[1]] += forceCreaseVal*dTheta_xl;
            f1[edge[0]] += forceCreaseVal*dTheta_xj;
            f1[edge[1]] += forceCreaseVal*dTheta_xk;
        }
        Real A = dot(rij,rkj)/dot(rkj,rkj);
        Real B = dot(rkl,rkj)/dot(rkj,rkj);
        Coord dAj = 1/dot(rkj,rkj) * ((2*A-1)*rkj - rij);
        Coord dBj = 1/dot(rkj,rkj) * (2*B*rkj - rkl);
        Coord dAk = 1/dot(rkj,rkj) * (-2*A*rkj + rij);
        Coord dBk = 1/dot(rkj,rkj) * ((1- 2*B)*rkj + rkl);

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

        for (int i=0;i<12;i++){
            for (int j=0;j<12;j++){
//                edgeStiffness[edgeNum].clear();
                edgeStiffness[edgeNum][i][j] = 0;
            }
        }
        forceCreaseVal = 0;
        dForcecrease=0;
        if ((edgesAssignment[edgeNum] == "M") || (edgesAssignment[edgeNum] == "V")) {
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
//            std::cout << "[ ";
//            for (int row=0;row<12;row++){
//                for (int col=0;col<12;col++){
//                    std::cout << edgeStiffness[edgeNum][row][col] << "\t";
//                }
//                std::cout << "\n";
//            }
//            std::cout << "] ";
//                    msg_warning() << << edgeStiffness[edgeNum][row][0] << " , " << " , " << edgeStiffness[edgeNum][row][1] << " , "<< edgeStiffness[edgeNum][row][2] << " , "<< edgeStiffness[edgeNum][row][3] << " , "<< edgeStiffness[edgeNum][row][4] << " , "<< edgeStiffness[edgeNum][row][5] << " , "<< edgeStiffness[edgeNum][row][6] << " , "<< edgeStiffness[edgeNum][row][7] << " , "<< edgeStiffness[edgeNum][row][8] << " , "<< edgeStiffness[edgeNum][row][9] << " , "<< edgeStiffness[edgeNum][row][10] << " , "<< edgeStiffness[edgeNum][row][11];


//            msg_warning() << edgeStiffness[edgeNum];
        }
//        msg_warning() << "Edge Stiffness Set ! " << edgeNum;
//        msg_warning() <<  edgeStiffness[edgeNum];

    }

    f.endEdit();
}


template <class DataTypes>
void RotationalSpringForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx)
{
    VecDeriv& df1 = *df.beginEdit();
    const VecDeriv& dx1 = dx.getValue();
    Real kFactor = (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(mparams, this->rayleighStiffness.getValue());
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
    df.endEdit();
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
    std::vector<sofa::type::Vector3> vertices, pointsValley, pointsMountain, pointsBorder;
    auto edgesAssignment = d_edgesAssignment.getValue();
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
    for (int edgeNum=0; edgeNum < m_topology->getNbEdges(); edgeNum++)
    {
        auto edge = m_topology->getEdge(edgeNum);
        if (edgesAssignment[edgeNum] == "M") {
        pointsMountain.push_back(x[edge[0]]);
        pointsMountain.push_back(x[edge[1]]);
        }
        else if (edgesAssignment[edgeNum] == "V") {
            pointsValley.push_back(x[edge[0]]);
            pointsValley.push_back(x[edge[1]]);
        }
        else
        {
            pointsBorder.push_back(x[edge[0]]);
            pointsBorder.push_back(x[edge[1]]);
        }
    }
    vparams->drawTool()->drawLines(pointsMountain,5.0,sofa::type::RGBAColor(0.5,0,0.5,1));
    vparams->drawTool()->drawLines(pointsValley,5.0,sofa::type::RGBAColor(0.5,0.5,0.5,1));
    vparams->drawTool()->drawLines(pointsBorder,5.0,sofa::type::RGBAColor(0,0.5,0,1));

    vparams->drawTool()->drawTriangles(vertices,colorVector);

    vparams->drawTool()->restoreLastState();
}


template<class DataTypes>
void RotationalSpringForceField<DataTypes>::addKToMatrix(sofa::linearalgebra::BaseMatrix *mat, SReal k, unsigned int &offset)
{
//    for(unsigned i=0; i< _indexedElements->size() ; i++)
//    {
//        StiffnessMatrix JKJt,RJKJtRt;
//        computeElementStiffnessMatrix(JKJt, RJKJtRt, _materialsStiffnesses[i], _strainDisplacements[i], _rotations[i]);
//        this->addToMatrix(mat,offset,(*_indexedElements)[i],RJKJtRt,-k);
//    }
    auto edgesAssignment = d_edgesAssignment.getValue();
    for (int edgeNum=0; edgeNum < m_topology->getNbEdges(); edgeNum++)
    {
        if ((edgesAssignment[edgeNum] == "M") || (edgesAssignment[edgeNum] == "V")) {
            this->addToMatrix(mat,offset,nodes[edgeNum],edgeStiffness[edgeNum],k);
        }
    }
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
        msg_warning() << "Input crease stiffness is not possible: " << val << ", setting default value: 1000";
        f_kcrease.setValue(1000);
    }
    else if (val != f_kcrease.getValue())
    {
        f_kcrease.setValue(val);
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
