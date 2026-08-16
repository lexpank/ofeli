/****************************************************************************
**
** Copyright (C) 2010-2012 Fabien Bessy.
** Contact: fabien.bessy@gmail.com
**
** This file is part of project Ofeli.
**
** http://www.cecill.info/licences/Licence_CeCILL_V2-en.html
** You may use this file under the terms of the CeCILL license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Fabien Bessy and its Subsidiary(-ies) nor the
**     names of its contributors may be used to endorse or promote products
**     derived from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
**
****************************************************************************/

//! \file main.cpp
//! \brief Ofeli
//! \author Fabien Bessy
//! \date August 2012

#include <QApplication>
#include <QCoreApplication>
#include <QLibraryInfo>

#include "imageviewer.hpp"

int main(int argc, char *argv[])
{
    // to have the windows theme
    //QApplication::setStyle(new QWindowsStyle);

    QApplication app(argc, argv);
    QCoreApplication::setApplicationName(QStringLiteral("Ofeli"));
    QCoreApplication::setOrganizationName(QStringLiteral("Bessy"));
    QCoreApplication::setApplicationVersion(QStringLiteral(OFELI_VERSION));

    QSettings settings("Bessy", "Ofeli");
    int language = settings.value("language",0).toInt();

    QTranslator translator_qt_widget;
    QTranslator translator_ofeli;
    QString locale;

    switch( language )
    {
    case 0 :
    {
        locale = QLocale::system().name().section('_', 0, 0);
        break;
    }
    case 1 :
    {
        locale = QString("en");
        break;
    }
    case 2 :
    {
        locale = QString("fr");
        break;
    }
    default :
    {
        locale = QString("en");
    }
    }

    if( translator_qt_widget.load( QString("qt_") + locale, QLibraryInfo::path(QLibraryInfo::TranslationsPath) ) )
    {
        app.installTranslator(&translator_qt_widget);
    }

    if( translator_ofeli.load( QString(":/i18n/Ofeli_") + locale + QStringLiteral(".qm") ) )
    {
        app.installTranslator(&translator_ofeli);
    }

    ofeli::ImageViewer imageViewer;
    imageViewer.show();

    return app.exec();
}
