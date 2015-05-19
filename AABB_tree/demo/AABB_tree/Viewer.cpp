#include "Viewer.h"
#include "Scene.h"
#include <QMouseEvent>
#include <QOpenGLFunctions>
#include <QMessageBox>
#include <QDebug>



Viewer::Viewer(QWidget* parent)
  : QGLViewer(parent),
    m_pScene(NULL),
    m_custom_mouse(false),
    translation_mode(false),
    areOpenGLBuffersInitialized(false)

{
}

void Viewer::setScene(Scene* pScene)
{
  this->m_pScene = pScene;
}

void Viewer::draw()
{
 QOpenGLFunctions *gl = QOpenGLContext::currentContext()->functions();
 gl->glClear(GL_COLOR_BUFFER_BIT);

  if(m_pScene != NULL && areOpenGLBuffersInitialized
          )
  {
      m_pScene->draw(this);
  }
  else
  {
      QMessageBox box;
      box.setText("Error : No scene");
      box.exec();
  }


}

void Viewer::initializeGL()
{
    QOpenGLWidget::initializeGL();
    // Set up the rendering context, load shaders and other resources, etc.:
    QOpenGLFunctions *gl = QOpenGLContext::currentContext()->functions();
    gl->initializeOpenGLFunctions();
    setBackgroundColor(::Qt::white);
    m_pScene->setGL(gl);
    m_pScene->compile_shaders();
    areOpenGLBuffersInitialized=true;
}

void Viewer::mousePressEvent(QMouseEvent* e)
{
  if ( e->modifiers() == Qt::ControlModifier )
  {
    m_pScene->set_fast_distance(true);
    // Refresh distance function
    m_pScene->cutting_plane();
    m_custom_mouse = true;
  }
  else if(translation_mode)
  {
      qDebug()<< "translate ON";
       setMouseBinding(Qt::Key(0),Qt::NoModifier, Qt::LeftButton, CAMERA, TRANSLATE);
  }
  else if(!translation_mode)
  {
      qDebug()<< "translate OFF";
       setMouseBinding(Qt::Key(0),Qt::NoModifier, Qt::LeftButton, CAMERA, ROTATE);
  }
  
  QGLViewer::mousePressEvent(e);
}

void Viewer::mouseReleaseEvent(QMouseEvent* e)
{
  if ( m_custom_mouse )
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

