#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtOpenGL/qgl.h>
#include <CGAL/Qt/DemosMainWindow.h>

class QDragEnterEvent;
class QDropEvent;
class Scene;
class Viewer;
namespace Ui {
  class MainWindow;
}


class MainWindow : 
  public CGAL::Qt::DemosMainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget* parent = 0);
  ~MainWindow();

public slots:
  void updateViewerBBox();
  void open(QString filename);

protected slots:
  // settings
  void quit();
  void readSettings();
  void writeSettings();

  // drag & drop
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);
  void closeEvent(QCloseEvent *event);

  // file menu
  void on_actionLoadPolyhedron_triggered();
  void on_actionSaveApproximation_triggered();

  // edit menu
  void on_actionSaveSnapshot_triggered();

  // set metric menu
  void on_actionL21_triggered();
  void on_actionL2_triggered();
  void on_actionCompact_triggered();

  // operations menu
  void on_actionApproximation_triggered();
  void on_actionSeeding_triggered();
  void on_actionFit_triggered();
  void on_actionMeshing_triggered();
  void on_actionAdd_triggered();
  void on_actionTeleport_triggered();
  void on_actionSplit_triggered();

  // view menu
  void on_actionViewPolyhedron_triggered();
  void on_actionViewWireframe_triggered();
  void on_actionViewBoundary_triggered();
  void on_actionViewProxies_triggered();
  void on_actionViewAnchors_triggered();
  void on_actionViewApproximation_triggered();

private:
  Scene *m_pScene;
  Viewer *m_pViewer;
  Ui::MainWindow *ui;
};

#endif // ifndef MAINWINDOW_H
