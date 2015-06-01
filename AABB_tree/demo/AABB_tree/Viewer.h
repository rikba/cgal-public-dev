#ifndef VIEWER_H
#define VIEWER_H

#include <QGLViewer/qglviewer.h>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>

// forward declarations
class QWidget;
class Scene;

class Viewer : public QGLViewer{

  Q_OBJECT

public:
  Viewer(QWidget * parent);

  // overload several QGLViewer virtual functions
  void draw();

  /*renders the scene in a frame buffer object, in a gray scale that represents the depth value of the rendered objects,
  so it can use glReadPixels and determine the z-value of a point clicked in the scene. To do so the function requires a vector containing
  all the programs used in the draw function of the scene, so it can use them with a special fragment shader.
    */
  void draw_depth_value(std::vector<QOpenGLShaderProgram *>, QMouseEvent* e);
  //void paintGL();
  virtual void initializeGL();
  void setScene(Scene* pScene);
  bool areOpenGLBuffersInitialized;
  bool settingPivotPoint;
  bool frame_manipulation;
  int count;
  bool isMultiTouch;

protected:
  virtual void mousePressEvent(QMouseEvent* e);
  virtual void mouseReleaseEvent(QMouseEvent* e);
  
private:
  Scene* m_pScene;
  bool m_custom_mouse;
  QOpenGLFunctions *gl;
  std::vector<float> pos_pivot;
  QOpenGLBuffer buffers;
  QOpenGLVertexArrayObject vao;
  QOpenGLShaderProgram rendering_program;
  void compute_pivot();

  int pivot_Location;
  int mvpLocation;
  bool event(QEvent *e);
  void rotateBy(qglviewer::Quaternion q_);
  void translateBy(qreal x, qreal y, qreal z);
  void scaleBy(qreal factor);


}; // end class Viewer

#endif // VIEWER_H
