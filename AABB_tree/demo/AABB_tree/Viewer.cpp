#include "Viewer.h"
#include "Scene.h"
#include <QMouseEvent>
#include <QOpenGLFramebufferObject>
#include <QOpenGLFunctions>
#include <QMessageBox>
#include <QDebug>
#include <manipulatedCameraFrame.h>
#include <cmath>


Viewer::Viewer(QWidget* parent)
    : QGLViewer(parent),
      m_pScene(NULL),
      m_custom_mouse(false),
      settingPivotPoint(false),
      areOpenGLBuffersInitialized(false)

{
}

void Viewer::setScene(Scene* pScene)
{
    this->m_pScene = pScene;
}
void Viewer::draw()
{
    QGLViewer::draw();
    QOpenGLFunctions *gl = QOpenGLContext::currentContext()->functions();

    if(m_pScene != NULL )
    {
        gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        m_pScene->draw(this);
    }
    else
    {
        QMessageBox box;
        box.setText("Error : No scene");
        box.exec();
    }

   if(settingPivotPoint)
   {


    }
}

void Viewer::initializeGL()
{
    QGLViewer::initializeGL();
    // Set up the rendering context, load shaders and other resources, etc.:
    gl = QOpenGLContext::currentContext()->functions();
    gl->initializeOpenGLFunctions();
    setBackgroundColor(::Qt::white);
    m_pScene->setGL(gl);
    m_pScene->compile_shaders();
}



void Viewer::mousePressEvent(QMouseEvent* e)
{
    if ( e->modifiers() == Qt::ControlModifier  || frame_manipulation)
    {
        m_pScene->set_fast_distance(true);
        // Refresh distance function
        m_pScene->cutting_plane();
        m_custom_mouse = true;
    }

    //setMouseBinding(Qt::Key(0),Qt::NoModifier, Qt::LeftButton, NO_CLICK_ACTION);
    if(settingPivotPoint)
    {
        makeCurrent();
        if(!setPivotPointFromPixelGLES(m_pScene->programs, camera(), e->pos()))
            camera()->setPivotPoint(sceneCenter());
       setVisualHintsMask(1,2000);


    }
    else if(!settingPivotPoint && frame_manipulation)
    {
        setMouseBinding(Qt::Key(0),Qt::NoModifier, Qt::LeftButton, FRAME, ROTATE);


    }
    else
    {
        setMouseBinding(Qt::NoModifier, Qt::LeftButton, CAMERA, ROTATE);

    }

    QGLViewer::mousePressEvent(e);
}

void Viewer::mouseReleaseEvent(QMouseEvent* e)
{
    if ( m_custom_mouse  )
    {
        m_pScene->set_fast_distance(false);
        // Recompute distance function
        QApplication::setOverrideCursor(Qt::WaitCursor);
        m_pScene->cutting_plane();
        QApplication::restoreOverrideCursor();

        m_custom_mouse = false;
    }

    QGLViewer::mouseReleaseEvent(e);
}
/*
void Viewer::drawAxis(qreal length)
{

qglviewer::AxisData data;
    pos_arrows.resize(0);
    normals_arrows.resize(0);
    color_arrows.resize(0);

    data.vertices=&pos_arrows;
    data.normals=&normals_arrows;
    data.colors=&color_arrows;
    //float color[4];
    //color[0] = 0.7f;  color[1] = 0.7f;  color[2] = 1.0f;  color[3] = 1.0f;

    drawArrowGLES(0.5*length, 10, qglviewer::Vec(0,0,0),qglviewer::Vec(length,0,0),qglviewer::Vec(1,0,0),data);

    //color[0] = 1.0f;  color[1] = 0.7f;  color[2] = 0.7f;  color[3] = 1.0f;
    //QGLViewer::drawArrow(length, 0.01*length);
    drawArrowGLES(0.5*length, 10, qglviewer::Vec(0,0,0),qglviewer::Vec(0,length,0),qglviewer::Vec(0,1,0),data);

    //color[0] = 0.7f;  color[1] = 1.0f;  color[2] = 0.7f;  color[3] = 1.0f;
    //QGLViewer::drawArrow(length, 0.01*length);
    drawArrowGLES(0.5*length, 10, qglviewer::Vec(0,0,0),qglviewer::Vec(0,0,length),qglviewer::Vec(0,0,1),data);

}
*/
