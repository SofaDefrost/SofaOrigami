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

#include <SofaOrigami/config.h>

#include <sofa/gui/qt/DataWidget.h>

#include <QLabel>
#include <QVBoxLayout>
#include <QSlider>
#include <QString>


namespace sofa::gui::qt
{

/**
 * \brief Customization of the representation of Data<unsigned> types
 * in the gui. In the .cpp file this widget is registered to represent
 * myData from MyBehaviorModel in the gui.
 **/
class SOFA_SOFAORIGAMI_API MyDataWidgetUnsigned : public TDataWidget<unsigned>
{
    Q_OBJECT
public :
    // The class constructor takes a TData<unsigned> since it creates
    // a widget for a that particular data type.
    MyDataWidgetUnsigned(QWidget* parent, const char* name, core::objectmodel::Data<unsigned> *data):
        TDataWidget<unsigned>(parent, name,data) {};

    // In this method we  create the widgets and perform the signal / slots
    // connections.
    virtual bool createWidgets();
    virtual void setDataReadOnly(bool readOnly);
protected slots:
    void change();
protected:
    ///Implements how update the widgets knowing the data value.
    virtual void readFromData();
    ///Implements how to update the data, knowing the widget value.
    virtual void writeToData();
    QSlider *m_qslider;
    QLabel *m_label1;
    QLabel *m_label2;
};


} // namespace sofa::gui::qt
