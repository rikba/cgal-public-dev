#include "Viewer.h"
#include "Scene.h"
#include <QMouseEvent>
#include <QOpenGLFramebufferObject>
#include <QOpenGLFunctions>
#include <QMessageBox>
#include <QDebug>
#include"QGLViewer/manipulatedCameraFrame.h"


Viewer::Viewer(QWidget* parent)
    : QGLViewer(parent),
      m_pScene(NULL),
      m_custom_mouse(false),
      settingPivotPoint(false),
      frame_manipulation(false),
      areOpenGLBuffersInitialized(false)

{
    count = 0;
    isMultiTouch = false;
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

    //renders a point to show the pivot point
    if(settingPivotPoint)
    {
        vao.bind();
        QMatrix4x4 mvpMatrix;
        float mat[16];
        camera()->getModelViewProjectionMatrix(mat);
        for(int i=0; i < 16; i++)
        {
            mvpMatrix.data()[i] = (float)mat[i];
        }
        rendering_program.bind();
        mvpLocation = rendering_program.uniformLocation("mvp_matrix");
        rendering_program.setUniformValue(mvpLocation, mvpMatrix);
        rendering_program.release();
        rendering_program.bind();
        gl->glDrawArrays(GL_POINTS, 0, pos_pivot.size()/3);
        rendering_program.release();
        vao.release();
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

    if(!buffers.create())
    {
        std::cerr<<"VBO Creation FAILED"<<std::endl;

    }
    if(!vao.create())
    {
        std::cerr<<"VAO Creation FAILED"<<std::endl;
    }


    //Vertex source code
    const char vertex_source[] =
    {
        // "#version 330 \n"
        "attribute highp vec4 vertex;\n"
        "uniform highp mat4 mvp_matrix;\n"
        "void main(void)\n"
        "{\n"
        "   gl_PointSize = 10.0; \n"
        "   gl_Position = mvp_matrix * vertex;\n"
        "}"
    };
    //Fragment source code
    const char fragment_source[] =
    {
        //"#version 330 \n"
        "uniform highp vec4 color; \n"
        "void main(void) { \n"
        "gl_FragColor = vec4(0.6, 0.0, 0.0, 0.4); \n"
        "} \n"
        "\n"
    };
    QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader->compileSourceCode(vertex_source))
    {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    QOpenGLShader *fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader->compileSourceCode(fragment_source))
    {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }

    if(!rendering_program.addShader(vertex_shader))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program.addShader(fragment_shader))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program.link())
    {
        //std::cerr<<"linking Program FAILED"<<std::endl;
         qDebug() << rendering_program.log();
    }
}

void Viewer::compute_pivot()
{
 pos_pivot.resize(3);
 pos_pivot[0] = camera()->pivotPoint().x;
 pos_pivot[1] = camera()->pivotPoint().y;
 pos_pivot[2] = camera()->pivotPoint().z;

 vao.bind();
 buffers.bind();
 buffers.allocate(pos_pivot.data(), pos_pivot.size()*sizeof(float));
 pivot_Location = rendering_program.attributeLocation("vertex");
 rendering_program.bind();
 rendering_program.enableAttributeArray(pivot_Location);
 rendering_program.setAttributeBuffer(pivot_Location,GL_FLOAT,0,3);
 buffers.release();
 rendering_program.release();
 vao.release();

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
    if(!isMultiTouch)
    {
        if(settingPivotPoint)
        {
            makeCurrent();
            draw_depth_value((m_pScene->programs), e);
            compute_pivot();

        }
       else if(!settingPivotPoint && frame_manipulation)
        {
            setMouseBinding(Qt::Key(0),Qt::NoModifier, Qt::LeftButton, FRAME, ROTATE);

        }
        else
        {
            setMouseBinding(Qt::Key(0),Qt::NoModifier, Qt::LeftButton, CAMERA, ROTATE);
        }
    }
    else
    {
        setMouseBinding(Qt::Key(0),Qt::NoModifier, Qt::LeftButton, NO_CLICK_ACTION);
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


bool Viewer::event(QEvent *e)
{

    switch(e->type())
    {
    case QEvent::TouchBegin:
    {
        count = 0;
        break;
    }
    case QEvent::TouchUpdate:
    {
        count = (count +1)%2;
        QTouchEvent* te = (static_cast<QTouchEvent*>(e));
        if(te->touchPoints().count()==2 )
        {
            isMultiTouch = true;
            if(count == 0)
            {

                QTouchEvent::TouchPoint p0 = te->touchPoints().first();
                QTouchEvent::TouchPoint p1 = te->touchPoints().last();
                QLineF line1(p0.lastPos(), p1.lastPos());
                QLineF line2(p0.pos(), p1.pos());
                qreal scaling =line1.length()-line2.length();
                QPointF c1(line1.x1()/2.0+line1.x2()/2.0, line1.y1()/2.0+line1.y2()/2.0);
                QPointF c2(line2.x1()/2.0+line2.x2()/2.0, line2.y1()/2.0+line2.y2()/2.0);
                QVector2D trans(c2-c1);

                translateBy(trans.x(), trans.y(), 0);
                scaleBy(2*scaling);

                //The eye's position in the World's coordinate system
                //the pivot point position in the screen coordinate system
                qglviewer::Vec pivot = qglviewer::Vec(line1.x1()/2.0+line1.x2()/2.0,line1.y1()/2.0+line1.y2()/2.0,camera()->zNear());
                //the pivot point position in the World coordinate system
                pivot = camera()->unprojectedCoordinatesOf(pivot);
                //pivot.z = eye.z-camera()->zNear();


                qglviewer::Vec U= camera()->frame()->transformOf(camera()->frame()->inverseTransformOf(qglviewer::Vec(0.0, 0.0, -1.0)));
                //qglviewer::Vec c = qglviewer::Vec(line1.x1()/2.0+line1.x2()/2.0,line1.y1()/2.0+line1.y2()/2.0,0);

                //qglviewer::Vec centre = camera()->unprojectedCoordinatesOf(pivot);

                qreal angle = line2.angleTo(line1)*M_PI/180;//qglviewer::Vec(0,0,1)
                qglviewer::Quaternion quaternion;
                if(frame_manipulation)
                    quaternion = qglviewer::Quaternion(U, -2*angle);
                else
                    quaternion = qglviewer::Quaternion(U, -2*angle);
                rotateBy(quaternion);
                break;
            }
            else if (count == 2)
            {
                update();
            }
        }
    }
    case QEvent::TouchEnd:
    {
        isMultiTouch = false;
        break;
    }
    default:
    {
        QGLViewer::event(e);
        update();
    }
    }

    return true;
}

void Viewer::rotateBy(qglviewer::Quaternion q_)
{
    if(!frame_manipulation)
        camera()->frame()->rotateAroundPoint(q_, camera()->pivotPoint());
    else
        manipulatedFrame()->rotateAroundPoint(q_, camera()->pivotPoint());
}

void Viewer::scaleBy(qreal factor)
{
    if(!frame_manipulation)
        camera()->frame()->translate(camera()->frame()->inverseTransformOf(qglviewer::Vec(0.0, 0.0, 0.2*camera()->flySpeed()*factor)));
    else
        manipulatedFrame()->translate(manipulatedFrame()->inverseTransformOf(qglviewer::Vec(0.0, 0.0, 0.2*camera()->flySpeed()*factor)));
}

void Viewer::translateBy(qreal x, qreal y, qreal z )
{

    qglviewer::Vec trans(-x, y, 0);
    if(!frame_manipulation)
    {
        // Scale to fit the screen mouse displacement
        trans *= 2.0 * tan(camera()->fieldOfView()/2.0) *
                fabs((camera()->frame()->coordinatesOf(camera()->frame()->pivotPoint())).z) / camera()->screenHeight();
        // Transform to world coordinate system.
        trans = camera()->frame()->orientation().rotate(camera()->frame()->translationSensitivity()*trans);
        // And then down to frame
        if (camera()->frame()->referenceFrame())
            trans = camera()->frame()->referenceFrame()->transformOf(trans);
        camera()->frame()->translate(2*trans);
    }
    else
    {
        // Scale to fit the screen mouse displacement
        trans *= 2.0 * tan(camera()->fieldOfView()/2.0) * fabs((camera()->frame()->coordinatesOf(manipulatedFrame()->position())).z) / camera()->screenHeight();
        // Transform to world coordinate system.
        trans = camera()->frame()->orientation().rotate(manipulatedFrame()->translationSensitivity()*trans);
        // And then down to frame
        if (manipulatedFrame()->referenceFrame()) trans = manipulatedFrame()->referenceFrame()->transformOf(trans);
        manipulatedFrame()->translate(-trans);

    }
}

void Viewer::draw_depth_value(std::vector<QOpenGLShaderProgram*> programs, QMouseEvent* e)
{
    struct couple
    {
        QByteArray code;
        int program_index;
        int shader_index;
    };

    static const int size = programs.size();

    std::vector<couple> original_shaders;
    //The fragmentertex source code
    const char grayscale_fragment_source[] =
    {
        //"#version 330 \n"
        "void main(void) { \n"
        "gl_FragColor = vec4(vec3(gl_FragCoord.z), 1.0); \n"
        "} \n"
        "\n"
    };


    for(int i=0; i<size; i++)
    {
        for(int j=0; j<programs[i]->shaders().size(); j++)
        {
            if(programs[i]->shaders().at(j)->shaderType() == QOpenGLShader::Fragment)
            {
                //copies the original shaders of each program
                couple c;
                c.code = programs[i]->shaders().at(j)->sourceCode();
                c.program_index = i;
                c.shader_index = j;
                original_shaders.push_back(c);
                //replace their fragment shaders so they display in a grayscale
                programs[i]->shaders().at(j)->compileSourceCode(grayscale_fragment_source);
            }
            programs[i]->link();
        }
    }
    //the FBO in which the grayscale image will be rendered
    QOpenGLFramebufferObject *fbo = new QOpenGLFramebufferObject(camera()->screenWidth(), camera()->screenWidth());
    fbo->bind();
    //make the lines thicker so it is easier to click
    gl->glLineWidth(10.0);
    //draws the image in the fbo
    paintGL();
    gl->glLineWidth(1.0);

    //determines the size of the buffer
    int deviceWidth = camera()->screenWidth();
    int deviceHeight = camera()->screenHeight();
    int rowLength = deviceWidth * 4; // data asked in RGBA,so 4 bytes.

    const static int dataLength = rowLength * deviceHeight;
    GLubyte* buffer = new GLubyte[dataLength];

    // Qt uses upper corner for its origin while GL uses the lower corner.
    gl->glReadPixels(e->pos().x(), camera()->screenHeight()-1-e->pos().y(), 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
    //reset the fbo to the one rendered on-screen, now that we have our information
    fbo->release();
    //resets the originals programs
    for(int i=0; i<original_shaders.size(); i++)
    {
        programs[original_shaders[i].program_index]->shaders().at(original_shaders[i].shader_index)->compileSourceCode(original_shaders[i].code);
        programs[original_shaders[i].program_index]->link();
    }
    //depth value needs to be between 0 and 1.
    float depth = buffer[0]/255.0;
    qglviewer::Vec point(e->pos().x(), e->pos().y(), depth);
    point = camera()->unprojectedCoordinatesOf(point);

    //if depth is 1, then it is the zFar plane that is hit, so there is nothing rendered along the ray.
    bool found = depth<1;

    if (found)
    {
        camera()->setPivotPoint(point);
    }
    else
        camera()->setPivotPoint(sceneCenter());

    setVisualHintsMask(1);

    update();

}
