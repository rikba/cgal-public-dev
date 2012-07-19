#include <fstream>
#include <cmath>
#include <list>
// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
//My two headers
#include <include/Delaunay_mesh_hyperbolic_criteria_2.h>
#include <Mathieu_traits.h>
//other
#include <CGAL/Aff_transformation_2.h>

#include <CGAL/point_generators_2.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QGraphicsEllipseItem>

// GraphicsView items and event filters (input classes)
#include "TriangulationCircumcircle.h"
#include "TriangulationMovingPoint.h"
#include "TriangulationPointInputAndConflictZone.h"
#include "TriangulationGraphicsItem.h"

// for viewportsBbox
#include <CGAL/Qt/utility.h>

// the two base classes
#include "ui_Delaunay_triangulation_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel R;
typedef CGAL::Mathieu_traits<R> K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> Delaunay;
typedef CGAL::Delaunay_mesh_hyperbolic_criteria_2<Delaunay> Criteria;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Point Point;
typedef Delaunay::All_vertices_iterator All_vertices_iterator;
typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Construct_reflexion Construct_reflexion;
typedef CGAL::Aff_transformation_2<K>  Transformation;


class MainWindow :
public CGAL::Qt::DemosMainWindow,
public Ui::Delaunay_triangulation_2
{
	Q_OBJECT
	
private:  
	Delaunay dt;
	QGraphicsEllipseItem* disk;
	QGraphicsScene scene;  
	
	CGAL::Qt::TriangulationGraphicsItem<Delaunay> * dgi;
	CGAL::Qt::TriangulationPointInputAndConflictZone<Delaunay> * pi;
	CGAL::Qt::TriangulationCircumcircle<Delaunay> *tcc;
	Criteria criteria;
	CGAL::Delaunay_mesher_2<Delaunay,Criteria> mesher;
	
	
public:
	MainWindow();
	
	
	
	public slots:
	
	void processInput(CGAL::Object o);
	
	void Apply_symetry();
	void apply_rotation();
	void apply_reflexion();
	void on_actionCircumcenter_toggled(bool checked);
	
	void on_actionShowDelaunay_toggled(bool checked);
	
	void on_actionInsertPoint_toggled(bool checked);
	
	void on_actionInsertRandomPoints_triggered();
	
	void on_actionLoadPoints_triggered();
	
	void on_actionSavePoints_triggered();
	
	void on_actionClear_triggered();
	
	void on_actionRecenter_triggered();
	
	virtual void open(QString fileName);
	
signals:
	void changed();
};


MainWindow::MainWindow()
: DemosMainWindow(), criteria(0.001), mesher(dt,criteria)
{
	
	
	setupUi(this);
	
	this->graphicsView->setAcceptDrops(false);
	/*double a = 0.707106781;
	 double b = 90;
	 Vertex_handle va = dt.insert(Point(b,0));
	 Vertex_handle vb = dt.insert(Point(b*a,b*a));
	 Vertex_handle vc = dt.insert(Point(0,b));
	 Vertex_handle vd = dt.insert(Point(-b*a,b*a));
	 Vertex_handle ve = dt.insert(Point(-b,0));
	 Vertex_handle vf = dt.insert(Point(-b*a,-b*a));
	 Vertex_handle vg = dt.insert(Point(0,-b));
	 Vertex_handle vh = dt.insert(Point(b*a,-b*a));
	 
	 dt.insert_constraint(va, vb);
	 dt.insert_constraint(vb, vc);
	 dt.insert_constraint(vc, vd);
	 dt.insert_constraint(vd, ve);
	 dt.insert_constraint(ve, vf);
	 dt.insert_constraint(vf, vg);
	 dt.insert_constraint(vg, vh);
	 dt.insert_constraint(vh, va);*/
	qreal phi = 3.14159265 / 8;
	qreal xhi = 3.14159265 / 3;
	qreal rho = sqrt(cos(xhi)*cos(xhi)-sin(phi)*sin(phi));
	K traits;
	const Point_2 o(0.,0.);
	const Point_2 a(100*(cos(xhi)-sin(phi))/rho,0);
	const Point_2 b(100*(cos(phi)*cos(phi+xhi))/rho,100*(sin(phi)*cos(phi+xhi)/rho));
	/*const Point_2 c =  traits.construct_reflexion_object()(a,b,o);
	 const Point_2 d =  traits.construct_reflexion_object()(b,c,a);
	 const Point_2 e =  traits.construct_reflexion_object()(d,c,b);
	 const Point_2 f =  traits.construct_reflexion_object()(e,c,d);
	 const Point_2 g =  traits.construct_reflexion_object()(e,f,c);*/
	
	Vertex_handle vo = dt.insert(o);
	Vertex_handle va = dt.insert(a);
	Vertex_handle vb = dt.insert(b);
	/*Vertex_handle vc = dt.insert(c);
	 Vertex_handle vd = dt.insert(d);
	 Vertex_handle ve = dt.insert(e);
	 Vertex_handle vf = dt.insert(f);
	 Vertex_handle vg = dt.insert(g);*/
	
	
	dt.insert_constraint(vo, va);
	dt.insert_constraint(va, vb);
	dt.insert_constraint(vb, vo);
	/*dt.insert_constraint(vb, vc);
	 dt.insert_constraint(va, vc);
	 dt.insert_constraint(vc, vd);
	 dt.insert_constraint(vd, ve);
	 dt.insert_constraint(ve, vc);
	 dt.insert_constraint(vf, vc);
	 dt.insert_constraint(vf, ve);*/
	
	
	
	
	
	//mesher.init ();
	
	
	
	// Add Poincaré disk
	qreal origin_x = 0, origin_y = 0, radius = 100, diameter = 2*radius;
	qreal left_top_corner_x = origin_x - radius;
	qreal left_top_corner_y = origin_y - radius;
	qreal width = diameter, height = diameter;
	
	disk = new QGraphicsEllipseItem(left_top_corner_x, left_top_corner_y, width, height);
	scene.addItem(disk);
	
	// Add a GraphicItem for the Delaunay triangulation
	dgi = new CGAL::Qt::TriangulationGraphicsItem<Delaunay>(&dt);
	
	QObject::connect(this, SIGNAL(changed()),
					 dgi, SLOT(modelChanged()));
	
	dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
	scene.addItem(dgi);
	
	
	// Setup input handlers. They get events before the scene gets them
	// and the input they generate is passed to the triangulation with 
	// the signal/slot mechanism    
	pi = new CGAL::Qt::TriangulationPointInputAndConflictZone<Delaunay>(&scene, &dt, this );
	
	QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
					 this, SLOT(processInput(CGAL::Object)));
	
	
	tcc = new CGAL::Qt::TriangulationCircumcircle<Delaunay>(&scene, &dt, this);
	tcc->setPen(QPen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
	
	// 
	// Manual handling of actions
	//
	
	QObject::connect(this->actionQuit, SIGNAL(triggered()), 
					 this, SLOT(close()));
	
	// We put mutually exclusive actions in an QActionGroup
	QActionGroup* ag = new QActionGroup(this);
	ag->addAction(this->actionInsertPoint);
	ag->addAction(this->actionCircumcenter);
	
	// Check two actions 
	this->actionInsertPoint->setChecked(true);
	this->actionShowDelaunay->setChecked(true);
	
	//
	// Setup the scene and the view
	//
	scene.setItemIndexMethod(QGraphicsScene::NoIndex);
	scene.setSceneRect(left_top_corner_x, left_top_corner_y, width, height);
	this->graphicsView->setScene(&scene);
	this->graphicsView->setMouseTracking(true);
	
	// Turn the vertical axis upside down
	this->graphicsView->matrix().scale(1, -1);
	
	// The navigation adds zooming and translation functionality to the
	// QGraphicsView
	this->addNavigation(this->graphicsView);
	
	this->setupStatusBar();
	this->setupOptionsMenu();
	this->addAboutDemo(":/cgal/help/about_Delaunay_triangulation_2.html");
	this->addAboutCGAL();
	
	this->addRecentFiles(this->menuFile, this->actionQuit);
	connect(this, SIGNAL(openRecentFile(QString)),
			this, SLOT(open(QString)));
}

void MainWindow::Apply_symetry()
{
	K traits;
	
	qreal phi = 3.14159265 / 8;
	qreal xhi = 3.14159265 / 3;
	qreal rho = sqrt(cos(xhi)*cos(xhi)-sin(phi)*sin(phi));
	const Point_2 o(0,0);
	const Point_2 a(100*(cos(xhi)-sin(phi))/rho,0);
	const Point_2 b(100*(cos(phi)*cos(phi+xhi))/rho,100*(sin(phi)*cos(phi+xhi)/rho));
	const Point_2 c =  traits.construct_reflexion_object()(a,b,o);
	const Point_2 d =  traits.construct_reflexion_object()(b,c,a);
	const Point_2 e =  traits.construct_reflexion_object()(d,c,b);
	const Point_2 f =  traits.construct_reflexion_object()(e,c,d);
	const Point_2 g =  traits.construct_reflexion_object()(e,f,c);
	
	
	std::list<Point> temp;
	for (All_vertices_iterator vi = dt.all_vertices_begin();vi!= dt.all_vertices_end();vi++)
	{
		temp.push_back(vi->point());
	}
	
	std::list<Point>::iterator li;
	Point p;
	for (li = temp.begin(); li!=temp.end();li++)
	{
		p = (*li);
		p = traits.construct_reflexion_object()(a,b,p);
		dt.insert(p);
		(*li)= p;
	}
	for (li = temp.begin(); li!=temp.end();li++)
	{
		p = (*li);
		p = traits.construct_reflexion_object()(c,b,p);
		dt.insert(p);
		(*li)= p;
	}
	for (li = temp.begin(); li!=temp.end();li++)
	{
		p = (*li);
		p = traits.construct_reflexion_object()(d,c,p);
		dt.insert(p);
		(*li)= p;
	}
	for (li = temp.begin(); li!=temp.end();li++)
	{
		p = (*li);
		p = traits.construct_reflexion_object()(e,c,p);
		dt.insert(p);
		(*li)= p;
	}
	for (li = temp.begin(); li!=temp.end();li++)
	{
		p = (*li);
		p = traits.construct_reflexion_object()(e,f,p);
		dt.insert(p);
		(*li)= p;
	}
}
void MainWindow::apply_reflexion()
	{
		K traits;
		std::list<Point> temp;
		temp.clear();
		std::list<Point>::iterator li;
		Point p;

		for (All_vertices_iterator vi = dt.all_vertices_begin();vi!= dt.all_vertices_end();vi++)
	{
		temp.push_back(vi->point());
	}
	for (li = temp.begin(); li !=temp.end(); li++) {
		p = (*li);
		Point_2 q(p.x(),-p.y());
		dt.insert(q);
		(*li)=q;
	}
	}
void MainWindow::apply_rotation()
{
	K traits;
	std::list<Point> temp;

	std::list<Point>::iterator li;
	Point p;

	Transformation rotate(CGAL::ROTATION,sqrt(0.5),sqrt(0.5));
	temp.clear();
	for (All_vertices_iterator vi = dt.all_vertices_begin();vi!= dt.all_vertices_end();vi++)
	{
		temp.push_back(vi->point());
	}
	for (int i =0;i<8;i++)
	{
		for (li = temp.begin(); li !=temp.end(); li++) 
		{
			p = (*li);
			p = rotate(p);
			dt.insert(p);
			(*li)=p;
		}
	}
	
	
}
int count = 0;
void
MainWindow::processInput(CGAL::Object o)
{
	Point_2 p;
	if(CGAL::assign(p, o)){
		QPointF qp(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
		
		// note that if the point is on the boundary then the disk contains the point
		if(disk->contains(qp)){
			mesher.init ();
			//for (int i = 0; i< 1;i++){mesher.step_by_step_refine_mesh ();}
			mesher.refine_mesh ();
			//Apply_symetry();
			
			
		}
		else {
			if (count ==0)
			{Apply_symetry();count++;}
			else if (count ==1) {apply_reflexion(); count++;}
			else {apply_rotation();}
		}
		
	}
	emit(changed());
}


/* 
 *  Qt Automatic Connections
 *  http://doc.trolltech.com/4.4/designer-using-a-component.html#automatic-connections
 * 
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */
void
MainWindow::on_actionInsertPoint_toggled(bool checked)
{
	if(checked){
		scene.installEventFilter(pi);
		//scene.installEventFilter(trv);
	} else {
		scene.removeEventFilter(pi);
		//scene.removeEventFilter(trv);
	}
}


void
MainWindow::on_actionCircumcenter_toggled(bool checked)
{
	if(checked){
		scene.installEventFilter(tcc);
		tcc->show();
	} else {  
		scene.removeEventFilter(tcc);
		tcc->hide();
	}
}


void
MainWindow::on_actionShowDelaunay_toggled(bool checked)
{
	dgi->setVisibleEdges(checked);
}




void
MainWindow::on_actionClear_triggered()
{
	dt.clear();
	emit(changed());
}


void
MainWindow::on_actionInsertRandomPoints_triggered()
{
	QRectF rect = CGAL::Qt::viewportsBbox(&scene);
	CGAL::Qt::Converter<K> convert;  
	Iso_rectangle_2 isor = convert(rect);
	CGAL::Random_points_in_iso_rectangle_2<Point_2> pg(isor.min(), isor.max());
	bool ok = false;
	const int number_of_points = 
    QInputDialog::getInteger(this, 
                             tr("Number of random points"),
                             tr("Enter number of random points"),
							 100,
							 0,
							 std::numeric_limits<int>::max(),
							 1,
							 &ok);
	
	if(!ok) {
		return;
	}
	
	// wait cursor
	QApplication::setOverrideCursor(Qt::WaitCursor);
	std::vector<Point_2> points;
	points.reserve(number_of_points);
	for(int i = 0; i < number_of_points; ++i){
		points.push_back(*pg++);
	}
	dt.insert(points.begin(), points.end());
	// default cursor
	QApplication::restoreOverrideCursor();
	emit(changed());
}


void
MainWindow::on_actionLoadPoints_triggered()
{
	QString fileName = QFileDialog::getOpenFileName(this,
													tr("Open Points file"),
													".");
	if(! fileName.isEmpty()){
		open(fileName);
	}
}


void
MainWindow::open(QString fileName)
{
	// wait cursor
	QApplication::setOverrideCursor(Qt::WaitCursor);
	std::ifstream ifs(qPrintable(fileName));
	
	K::Point_2 p;
	std::vector<K::Point_2> points;
	while(ifs >> p) {
		points.push_back(p);
	}
	dt.insert(points.begin(), points.end());
	
	// default cursor
	QApplication::restoreOverrideCursor();
	this->addToRecentFiles(fileName);
	actionRecenter->trigger();
	emit(changed());
    
}

void
MainWindow::on_actionSavePoints_triggered()
{
	QString fileName = QFileDialog::getSaveFileName(this,
													tr("Save points"),
													".");
	if(! fileName.isEmpty()){
		std::ofstream ofs(qPrintable(fileName));
		for(Delaunay::Finite_vertices_iterator 
			vit = dt.finite_vertices_begin(),
			end = dt.finite_vertices_end();
			vit!= end; ++vit)
		{
			ofs << vit->point() << std::endl;
		}
	}
}


void
MainWindow::on_actionRecenter_triggered()
{
	this->graphicsView->setSceneRect(dgi->boundingRect());
	this->graphicsView->fitInView(dgi->boundingRect(), Qt::KeepAspectRatio);  
}


#include "Hyperbolic_Octagon_Mesh_2_demo.moc"

int main(int argc, char **argv)
{
	QApplication app(argc, argv);
	
	app.setOrganizationDomain("geometryfactory.com");
	app.setOrganizationName("GeometryFactory");
	app.setApplicationName("Delaunay_Octagon_Mesh_2 demo");
	
	// Import resources from libCGALQt4.
	// See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
	Q_INIT_RESOURCE(File);
	Q_INIT_RESOURCE(Triangulation_2);
	Q_INIT_RESOURCE(Input);
	Q_INIT_RESOURCE(CGAL);
	
	MainWindow mainWindow;
	mainWindow.show();
	
	QStringList args = app.arguments();
	args.removeAt(0);
	Q_FOREACH(QString filename, args) {
		mainWindow.open(filename);
	}
	
	return app.exec();
}
