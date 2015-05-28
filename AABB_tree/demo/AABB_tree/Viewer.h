#ifndef VIEWER_H
#define VIEWER_H

#include <QGLViewer/qglviewer.h>
#include <QOpenGLFunctions>

// forward declarations
class QWidget;
class Scene;

class Viewer : public QGLViewer{

  Q_OBJECT

public:
  Viewer(QWidget * parent);

  // overload several QGLViewer virtual functions
  void draw();
  //void paintGL();
  virtual void initializeGL();
  void setScene(Scene* pScene);
  bool areOpenGLBuffersInitialized;
  bool translation_mode;

protected:
  virtual void mousePressEvent(QMouseEvent* e);
  virtual void mouseReleaseEvent(QMouseEvent* e);
  
private:
  Scene* m_pScene;
  bool m_custom_mouse;
  QOpenGLFunctions *gl;
}; // end class Viewer

#endif // VIEWER_H
