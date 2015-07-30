
#include "opengl_tools.h"
#include "Scene_edit_polyhedron_item.h"
#include <boost/foreach.hpp>
#include <algorithm>


#include <CGAL/gl_render.h>
struct light_info
{
    //position
    GLfloat position[4];

    //ambient
    GLfloat ambient[4];

    //diffuse
    GLfloat diffuse[4];

    //specular
    GLfloat specular[4];
};
Scene_edit_polyhedron_item::Scene_edit_polyhedron_item
(Scene_polyhedron_item* poly_item,
 Ui::DeformMesh* ui_widget,
 QMainWindow* mw)
    : ui_widget(ui_widget),
      poly_item(poly_item),
      deform_mesh(*(poly_item->polyhedron()), Deform_mesh::Vertex_index_map(), Deform_mesh::Hedge_index_map(), Array_based_vertex_point_map(&positions)),
      is_rot_free(true),
      own_poly_item(true),
      ROI_points(0),
      control_points(0),
      control_color(0),
      ROI_color(0),
      pos_sphere(0),
      normals_sphere(0),
      k_ring_selector(poly_item, mw, Scene_polyhedron_item_k_ring_selection::Active_handle::VERTEX, true),
      Scene_item(20,8)
{
    mw->installEventFilter(this);
    // bind vertex picking
    connect(&k_ring_selector, SIGNAL(selected(const std::set<Polyhedron::Vertex_handle>&)), this,
            SLOT(selected(const std::set<Polyhedron::Vertex_handle>&)));

    poly_item->set_color_vector_read_only(true); // to prevent recomputation of color vector in changed()
    poly_item->update_vertex_indices();

    length_of_axis = bbox().diagonal_length() / 15.0;

    // interleave events of viewer (there is only one viewer)
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    viewer->installEventFilter(this);

    // create an empty group of control vertices for starting
    create_ctrl_vertices_group();

    // start QObject's timer for continuous effects
    // (deforming mesh while mouse not moving)
    startTimer(0);

    // Required for drawing functionality
    positions.resize(num_vertices(*polyhedron())*3);
    normals.resize(positions.size());
    Polyhedron::Vertex_iterator vb, ve;
    std::size_t counter = 0;
    for(vb=polyhedron()->vertices_begin(), ve = polyhedron()->vertices_end();vb != ve; ++vb, ++counter) {
        positions[counter*3] = vb->point().x();
        positions[counter*3+1] = vb->point().y();
        positions[counter*3+2] = vb->point().z();

        const Polyhedron::Traits::Vector_3& n =
                compute_vertex_normal<Polyhedron::Vertex, Polyhedron::Traits>(*vb);

        normals[counter*3] = n.x();
        normals[counter*3+1] = n.y();
        normals[counter*3+2] = n.z();
    }
    tris.resize(polyhedron()->size_of_facets()*3);
    counter = 0;
    for(Polyhedron::Facet_handle fb = polyhedron()->facets_begin(); fb != polyhedron()->facets_end(); ++fb, ++counter) {
        tris[counter*3] =  static_cast<unsigned int>(fb->halfedge()->vertex()->id());
        tris[counter*3+1] = static_cast<unsigned int>(fb->halfedge()->next()->vertex()->id());
        tris[counter*3+2] = static_cast<unsigned int>(fb->halfedge()->prev()->vertex()->id());
    }

    edges.resize(polyhedron()->size_of_halfedges());
    counter = 0;
    for(Polyhedron::Edge_iterator eb = polyhedron()->edges_begin(); eb != polyhedron()->edges_end(); ++eb, ++counter) {
        edges[counter*2] = static_cast<unsigned int>(eb->vertex()->id());
        edges[counter*2+1] = static_cast<unsigned int>(eb->opposite()->vertex()->id());
    }
    qFunc.initializeOpenGLFunctions();
    //Generates an integer which will be used as ID for each buffer

    const char vertex_shader_source_bbox[] =
    {

        "attribute highp vec3 vertex; \n"
        "attribute highp vec3 colors; \n"

        "uniform highp mat4 mvp_matrix; \n"
        "uniform highp mat4 rotations; \n"
        "uniform highp vec3 translation; \n"
        "uniform highp vec3 translation_2; \n"
        "varying highp vec3 fColors; \n"
        " \n"

        "void main(void) \n"
        "{ \n"
        "   fColors = colors; \n"
        "   gl_Position = mvp_matrix * (rotations *(vec4(translation_2,0.0)+vec4(vertex,1.0) )+ vec4(translation,0.0)) ; \n"
        "} \n"
    };
    const char fragment_shader_source[]=
    {
        "varying highp vec3 fColors; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        " gl_FragColor = vec4(fColors, 1.0); \n"
        "} \n"
    };
    bbox_program.addShaderFromSourceCode(QOpenGLShader::Vertex,vertex_shader_source_bbox);
    bbox_program.addShaderFromSourceCode(QOpenGLShader::Fragment,fragment_shader_source);
    bbox_program.link();

    //the spheres :
    create_Sphere(length_of_axis/15.0);
    program_list_is_empty = true;
    changed();
}

Scene_edit_polyhedron_item::~Scene_edit_polyhedron_item()
{
    setVisible(false);
    while(is_there_any_ctrl_vertices_group())
    {
        delete_ctrl_vertices_group(false);
    }
 //   gluDeleteQuadric(quadric);
    if (own_poly_item) delete poly_item;

}
/////////////////////////////
/// For the Shader gestion///
void Scene_edit_polyhedron_item::initialize_buffers(Viewer_interface *viewer =0) const
{

    //vao for the facets
    {
        program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);

        program->bind();

        vaos[0]->bind();
        buffers[0].bind();
        buffers[0].allocate(positions.data(), positions.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[0].release();

        buffers[1].bind();
        buffers[1].allocate(normals.data(), normals.size()*sizeof(float));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        buffers[1].release();
        vaos[0]->release();
        program->release();
    }
    //vao for the ROI points
    {   program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[1]->bind();
        buffers[2].bind();
        buffers[2].allocate(ROI_points.data(), ROI_points.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[2].release();

        buffers[3].bind();
        buffers[3].allocate(ROI_color.data(), ROI_color.size()*sizeof(float));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[3].release();
        vaos[1]->release();
        program->release();
    }


   //vao for the edges
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[2]->bind();
        buffers[4].bind();
        buffers[4].allocate(positions.data(), positions.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[4].release();

        buffers[5].bind();
        buffers[5].allocate(color_edges.data(), color_edges.size()*sizeof(float));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[5].release();
        vaos[2]->release();
        program->release();
    }
    //vao for the ROI spheres
    {
        program = getShaderProgram(PROGRAM_INSTANCED, viewer);
        program->bind();
        vaos[3]->bind();
        buffers[6].bind();
        buffers[6].allocate(pos_sphere.data(), pos_sphere.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[6].release();

        buffers[7].bind();
        buffers[7].allocate(normals_sphere.data(), normals_sphere.size()*sizeof(float));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        buffers[7].release();

        buffers[8].bind();
        buffers[8].allocate(color_sphere_ROI.data(), color_sphere_ROI.size()*sizeof(float));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[8].release();

        buffers[9].bind();
        buffers[9].allocate(centers_ROI.data(), centers_ROI.size()*sizeof(float));
        program->enableAttributeArray("center");
        program->setAttributeBuffer("center",GL_FLOAT,0,3);
        buffers[9].release();

      //  qFunc.glVertexAttribDivisor(program->attributeLocation("center"), 1);
      //  qFunc.glVertexAttribDivisor(program->attributeLocation("colors"), 1);
        vaos[3]->release();
    }
    //vao for the BBOX
    {
        bbox_program.bind();
        vaos[4]->bind();
        buffers[10].bind();
        buffers[10].allocate(pos_bbox.data(), pos_bbox.size()*sizeof(float));
        bbox_program.enableAttributeArray("vertex");
        bbox_program.setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[10].release();

        buffers[11].bind();
        buffers[11].allocate(color_bbox.data(), color_bbox.size()*sizeof(float));
        bbox_program.enableAttributeArray("colors");
        bbox_program.setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[11].release();
        vaos[4]->release();
        bbox_program.release();
    }
    //vao for the control points
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[5]->bind();
        buffers[12].bind();
        buffers[12].allocate(control_points.data(), control_points.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[12].release();

        buffers[13].bind();
        buffers[13].allocate(control_color.data(), control_color.size()*sizeof(float));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[13].release();
        vaos[5]->release();
        program->release();
    }
    //vao for the control spheres
    {
        program = getShaderProgram(PROGRAM_INSTANCED, viewer);
        program->bind();
        vaos[6]->bind();
        buffers[14].bind();
        buffers[14].allocate(pos_sphere.data(), pos_sphere.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[14].release();

        buffers[15].bind();
        buffers[15].allocate(normals_sphere.data(), normals_sphere.size()*sizeof(float));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        buffers[15].release();

        buffers[16].bind();
        buffers[16].allocate(color_sphere_control.data(), color_sphere_control.size()*sizeof(float));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[16].release();

        buffers[17].bind();
        buffers[17].allocate(centers_control.data(), centers_control.size()*sizeof(float));
        program->enableAttributeArray("center");
        program->setAttributeBuffer("center",GL_FLOAT,0,3);
        buffers[17].release();

       // qFunc.glVertexAttribDivisor(program->attributeLocation("center"), 1);
       // qFunc.glVertexAttribDivisor(program->attributeLocation("colors"), 1);
        vaos[6]->release();
    }
    //vao for the axis
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[7]->bind();
        buffers[18].bind();
        buffers[18].allocate(pos_axis.data(), pos_axis.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[18].release();

        buffers[19].bind();
        buffers[19].allocate(color_lines.data(), color_lines.size()*sizeof(float));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[19].release();
        vaos[7]->release();
        program->release();
    }
    if(program_list_is_empty)
    {
    foreach(QOpenGLShaderProgram* prog, shader_programs)
        viewer->program_list.push_back(prog);
    k_ring_selector.edit_programs = viewer->program_list;
    program_list_is_empty = false;
    }
    are_buffers_filled = true;

}

void Scene_edit_polyhedron_item::compute_normals_and_vertices(void)
{

    ROI_points.clear();
    control_points.clear();
    BOOST_FOREACH(vertex_descriptor vd, deform_mesh.roi_vertices())
    {
        if(!deform_mesh.is_control_vertex(vd))
        {//gl_draw_point( vd->point() );
            ROI_points.push_back(vd->point().x());
            ROI_points.push_back(vd->point().y());
            ROI_points.push_back(vd->point().z());
        }
    }
    centers_ROI.resize(ROI_points.size());
    ROI_color.resize(ROI_points.size());
    color_sphere_ROI.resize(ROI_points.size());
    for(int i=0; i<centers_ROI.size(); i++)
    {
        centers_ROI[i] = ROI_points[i];
    }
    for(int i=0; i<ROI_color.size(); i++)
    {
        if(i%3==1)
        {
        ROI_color[i]=1.0;
        color_sphere_ROI[i]=1.0;

        }
        else
        {
        ROI_color[i]=0.0;
        color_sphere_ROI[i]=0.0;
        }
    }
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    for(Ctrl_vertices_group_data_list::const_iterator hgb_data = ctrl_vertex_frame_map.begin(); hgb_data != ctrl_vertex_frame_map.end(); ++hgb_data)
    {
        if(hgb_data->frame == viewer->manipulatedFrame())
        {
            // draw axis

            if(ui_widget->ActivatePivotingCheckBox->isChecked())
            {
                // draw bbox
                compute_bbox(hgb_data->bbox);
            }
        }
        // draw control vertices
        if(hgb_data == active_group)
        {
            //set color to red
            control_color.push_back(1.0);
            control_color.push_back(0.0);
            control_color.push_back(0.0);
        }
        else
        {
            //set color to blue
            control_color.push_back(0.0);
            control_color.push_back(0.0);
            control_color.push_back(1.0);
        }
        for(std::vector<vertex_descriptor>::const_iterator hb = hgb_data->ctrl_vertices_group.begin(); hb != hgb_data->ctrl_vertices_group.end(); ++hb)
        {
            control_points.push_back((*hb)->point().x());
            control_points.push_back((*hb)->point().y());
            control_points.push_back((*hb)->point().z());

        }
        centers_control.resize(control_points.size());
        for(int i=0; i<centers_control.size(); i++)
        {
            centers_control[i]=control_points[i];
        }
    }
    color_sphere_control.resize(control_color.size());
    for(int i=0; i<color_sphere_control.size(); i++)
    {
        color_sphere_control[i] = control_color[i];
    }

    //The edges color
    color_edges.resize(edges.size());
    for(int i =0; i< edges.size(); i++)
        color_edges[i]=0.0;

    //The box color
    color_bbox.resize(pos_bbox.size());
    for(int i =0; i< pos_bbox.size(); i++)
        color_bbox[i]=0.0;

    for(int i =0; i< pos_bbox.size(); i+=3)
        color_bbox[i]=1.0;

    //The axis

    pos_axis.resize(18);
    for(int i =0; i< 18; i++)
        pos_axis[i]=0.0;
    pos_axis[3] = length_of_axis; pos_axis[10] = length_of_axis; pos_axis[17] = length_of_axis;
    color_lines.resize(18);
    for(int i =0; i< 18; i++)
        color_lines[i]=0.0;

    color_lines[2] = 1.0; color_lines[5] = 1.0;
    color_lines[6] = 1.0; color_lines[9] = 1.0;
    color_lines[13] = 1.0; color_lines[16] = 1.0;

}

/////////////////////////////////////////////////////////
/////////// Most relevant functions lie here ///////////
void Scene_edit_polyhedron_item::deform()
{
    if(!is_there_any_ctrl_vertices()) { return; }

    for(Ctrl_vertices_group_data_list::iterator it = ctrl_vertex_frame_map.begin(); it != ctrl_vertex_frame_map.end(); ++it)
    { it->set_target_positions(); }
    deform_mesh.deform();

    poly_item->changed(); // now we need to call poly_item changed to delete AABB tree
    emit itemChanged();
}

void Scene_edit_polyhedron_item::timerEvent(QTimerEvent* /*event*/)
{ // just handle deformation - paint like selection is handled in eventFilter()
    if(state.ctrl_pressing && (state.left_button_pressing || state.right_button_pressing)) {
        if(!ui_widget->ActivatePivotingCheckBox->isChecked()) {
            deform();
        }
        else {
            emit itemChanged(); // for redraw while Pivoting (since we close signals of manipulatedFrames while pivoting,
            // for now redraw with timer)
        }
    }
}
bool Scene_edit_polyhedron_item::eventFilter(QObject* /*target*/, QEvent *event)
{
    // This filter is both filtering events from 'viewer' and 'main window'
    Mouse_keyboard_state_deformation old_state = state;
    ////////////////// TAKE EVENTS /////////////////////
    // key events

    if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)
    {
        QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
        Qt::KeyboardModifiers modifiers = keyEvent->modifiers();

        ctrl_pressing = modifiers.testFlag(Qt::ControlModifier);
        shift_pressing = modifiers.testFlag(Qt::ShiftModifier);
    }
   QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    // mouse events

    if(event->type() == QEvent::MouseButtonPress || event->type() == QEvent::MouseButtonRelease)
    {

        if(viewer->frame_manipulation || ctrl_pressing  )
        {
            state.ctrl_pressing = true;
        }
        else
            state.ctrl_pressing = false;
        if(viewer->selection_mode || shift_pressing )
            state.shift_pressing = true;
        else
            state.shift_pressing = false;
        QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
        if(mouse_event->button() == Qt::LeftButton) {
            state.left_button_pressing = event->type() == QEvent::MouseButtonPress;
        }
        if(mouse_event->button() == Qt::RightButton) {
            state.right_button_pressing = event->type() == QEvent::MouseButtonPress;
        }
    }
    ////////////////// //////////////// /////////////////////

    if(!poly_item->visible()) { return false; } // if not visible just update event state but don't do any action

    // check state changes between old and current state
    bool ctrl_pressed_now = ctrl_pressing && !old_state.ctrl_pressing;
    bool ctrl_released_now = !ctrl_pressing && old_state.ctrl_pressing;
    if(ctrl_pressed_now || ctrl_released_now || event->type() == QEvent::HoverMove)
    {// activate a handle manipulated frame
        QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
        const QPoint& p = viewer->mapFromGlobal(QCursor::pos());
        bool need_repaint = activate_closest_manipulated_frame(p.x(), p.y());

        if(need_repaint) { emit itemChanged(); }
    }

    return false;
}

void Scene_edit_polyhedron_item::draw_edges(Viewer_interface* viewer) const {

    Scene_item::draw();
    if(!are_buffers_filled)
        initialize_buffers(viewer);
    vaos[2]->bind();
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITHOUT_LIGHT);
    program->bind();
    qFunc.glDrawElements(GL_LINES, (GLsizei) edges.size(), GL_UNSIGNED_INT, edges.data());
    program->release();
    vaos[2]->release();

    if(rendering_mode == Wireframe) {
        draw_ROI_and_control_vertices(viewer);
    }

}
void Scene_edit_polyhedron_item::draw(Viewer_interface* viewer) const {
    Scene_item::draw();
    if(!are_buffers_filled)
    {
        initialize_buffers(viewer);

    }
    vaos[0]->bind();
    program = getShaderProgram(PROGRAM_WITH_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITH_LIGHT);
    program->bind();
    QColor color = this->color();
    program->setAttributeValue("colors", color);
    qFunc.glDrawElements(GL_TRIANGLES, (GLsizei) tris.size(), GL_UNSIGNED_INT, tris.data());
    program->release();
    vaos[0]->release();
    draw_edges(viewer);
    draw_ROI_and_control_vertices(viewer);


}

void Scene_edit_polyhedron_item::draw_ROI_and_control_vertices(Viewer_interface* viewer) const {

    Scene_item::draw();
    //CGAL::GL::Color color;
    //CGAL::GL::Point_size point_size; point_size.set_point_size(5);

   // color.set_rgb_color(0, 1.f, 0);
       // qFunc.glEnable(GL_POINT_SPRITE);
        //qFunc.glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
    if(ui_widget->ShowROICheckBox->isChecked()) {

        if(!ui_widget->ShowAsSphereCheckBox->isChecked()) {

            vaos[1]->bind();
            program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
            attrib_buffers(viewer,PROGRAM_WITHOUT_LIGHT);
            program->bind();
            qFunc.glDrawArrays(GL_POINTS, 0, ROI_points.size()/3);
            program->release();
            vaos[1]->release();
        }
        else{
            vaos[3]->bind();
            program = getShaderProgram(PROGRAM_INSTANCED);
            attrib_buffers(viewer,PROGRAM_INSTANCED);
            program->bind();
       //     qFunc.glDrawArraysInstanced(GL_TRIANGLES, 0, pos_sphere.size()/3, ROI_points.size()/3);
            program->release();
            vaos[3]->release();
        }
    }

    if(!ui_widget->ShowAsSphereCheckBox->isChecked()) {
        vaos[5]->bind();
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
        attrib_buffers(viewer,PROGRAM_WITHOUT_LIGHT);
        program->bind();
        qFunc.glDrawArrays(GL_POINTS, 0, control_points.size()/3);
        program->release();
        vaos[5]->release();
    }
    else{
        vaos[6]->bind();
        program = getShaderProgram(PROGRAM_INSTANCED);
        attrib_buffers(viewer,PROGRAM_INSTANCED);
        program->bind();
   //     qFunc.glDrawArraysInstanced(GL_TRIANGLES, 0, pos_sphere.size()/3, control_points.size()/3);
        program->release();
        vaos[6]->release();
    }

    QGLViewer* viewerB = *QGLViewer::QGLViewerPool().begin();
    for(Ctrl_vertices_group_data_list::const_iterator hgb_data = ctrl_vertex_frame_map.begin(); hgb_data != ctrl_vertex_frame_map.end(); ++hgb_data)
    {
        if(hgb_data->frame == viewerB->manipulatedFrame())
        {
            GLfloat f_matrix[16];
            for(int i =0; i<16; i++)
                f_matrix[i] = hgb_data->frame->matrix()[i];
            QMatrix4x4 f_mat;
                for(int i=0; i<16; i++)
                    f_mat.data()[i] = (float)f_matrix[i];
            vaos[7]->bind();
            program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
            attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
            program->bind();
            program->setUniformValue("f_matrix", f_mat);
            qFunc.glDrawArrays(GL_LINES, 0, pos_axis.size()/3);
            program->release();
            vaos[7]->release();

            //QGLViewer::drawAxis(length_of_axis);
            // draw bbox
            if(!ui_widget->ActivatePivotingCheckBox->isChecked())
            {
                GLfloat colors[3];
                GLfloat f_matrix[16];
                GLfloat trans[3];
                GLfloat trans2[3];
                colors[0]=1.0;
                colors[1]=0.0;
                colors[2]=0.0;

                trans[0] = hgb_data->frame->position().x;
                trans[1] = hgb_data->frame->position().y;
                trans[2] = hgb_data->frame->position().z;

                trans2[0] = -hgb_data->frame_initial_center.x;
                trans2[1] = -hgb_data->frame_initial_center.y;
                trans2[2] = -hgb_data->frame_initial_center.z;

                for(int i =0; i<16; i++)
                    f_matrix[i] = hgb_data->frame->orientation().matrix()[i];
                QMatrix4x4 f_mat;
                QMatrix4x4 mvp_mat;

                QVector3D vec(trans[0], trans[1], trans[2]);
                QVector3D vec2(trans2[0], trans2[1], trans2[2]);
                    for(int i=0; i<16; i++)
                        f_mat.data()[i] = (float)f_matrix[i];
                    GLfloat temp_mat[16];
                    viewer->camera()->getModelViewProjectionMatrix(temp_mat);
                    for(int i=0; i<16; i++)
                        mvp_mat.data()[i] = (float)temp_mat[i];
                vaos[4]->bind();
                bbox_program.bind();
                bbox_program.setUniformValue("rotations", f_mat);
                bbox_program.setUniformValue("translation", vec);
                bbox_program.setUniformValue("translation_2", vec2);
                bbox_program.setUniformValue("mvp_matrix", mvp_mat);
                qFunc.glDrawArrays(GL_LINES, 0, pos_bbox.size()/3);
                bbox_program.release();
                vaos[4]->release();
            }
        }
    }

}


void Scene_edit_polyhedron_item::compute_bbox(const Scene_interface::Bbox& bb){
    pos_bbox.resize(24*3);

    pos_bbox[0]=bb.xmin; pos_bbox[1]=bb.ymin; pos_bbox[2]=bb.zmin;
    pos_bbox[3]=bb.xmax; pos_bbox[4]=bb.ymin; pos_bbox[5]=bb.zmin;
    pos_bbox[6]=bb.xmin; pos_bbox[7]=bb.ymin; pos_bbox[8]=bb.zmin;
    pos_bbox[9]=bb.xmin; pos_bbox[10]=bb.ymax; pos_bbox[11]=bb.zmin;

    pos_bbox[12]=bb.xmin; pos_bbox[13]=bb.ymin; pos_bbox[14]=bb.zmin;
    pos_bbox[15]=bb.xmin; pos_bbox[16]=bb.ymin; pos_bbox[17]=bb.zmax;
    pos_bbox[18]= bb.xmax; pos_bbox[19]=bb.ymin; pos_bbox[20]=bb.zmin;
    pos_bbox[21]= bb.xmax; pos_bbox[22]=bb.ymax; pos_bbox[23]=bb.zmin;

    pos_bbox[24]= bb.xmax; pos_bbox[25]=bb.ymin; pos_bbox[26]=bb.zmin;
    pos_bbox[27]= bb.xmax; pos_bbox[28]=bb.ymin; pos_bbox[29]=bb.zmax;
    pos_bbox[30]=bb.xmin; pos_bbox[31]=bb.ymax; pos_bbox[32]=bb.zmin;
    pos_bbox[33]=bb.xmax; pos_bbox[34]=bb.ymax; pos_bbox[35]=bb.zmin;

    pos_bbox[36]=bb.xmin; pos_bbox[37]=bb.ymax; pos_bbox[38]=bb.zmin;
    pos_bbox[39]=bb.xmin; pos_bbox[40]=bb.ymax; pos_bbox[41]=bb.zmax;
    pos_bbox[42]=bb.xmin; pos_bbox[43]=bb.ymin; pos_bbox[44]=bb.zmax;
    pos_bbox[45]=bb.xmax; pos_bbox[46]=bb.ymin; pos_bbox[47]=bb.zmax;

    pos_bbox[48]=bb.xmin; pos_bbox[49]=bb.ymin; pos_bbox[50]=bb.zmax;
    pos_bbox[51]=bb.xmin; pos_bbox[52]=bb.ymax; pos_bbox[53]=bb.zmax;
    pos_bbox[54]=bb.xmax; pos_bbox[55]=bb.ymax; pos_bbox[56]=bb.zmax;
    pos_bbox[57]=bb.xmin; pos_bbox[58]=bb.ymax; pos_bbox[59]=bb.zmax;

    pos_bbox[60]=bb.xmax; pos_bbox[61]=bb.ymax; pos_bbox[62]=bb.zmax;
    pos_bbox[63]=bb.xmax; pos_bbox[64]=bb.ymin; pos_bbox[65]=bb.zmax;
    pos_bbox[66]=bb.xmax; pos_bbox[67]=bb.ymax; pos_bbox[68]=bb.zmax;
    pos_bbox[69]=bb.xmax; pos_bbox[70]=bb.ymax; pos_bbox[71]=bb.zmin;

}

void Scene_edit_polyhedron_item::changed()
{

    compute_normals_and_vertices();
    update_normals();
    are_buffers_filled = false;

}

Scene_polyhedron_item* Scene_edit_polyhedron_item::to_polyhedron_item() {
    Scene_polyhedron_item* poly_item_tmp = poly_item;
    poly_item->set_color_vector_read_only(false);
    own_poly_item=false;
    return poly_item_tmp;
}

Polyhedron* Scene_edit_polyhedron_item::polyhedron()
{ return poly_item->polyhedron(); }
const Polyhedron* Scene_edit_polyhedron_item::polyhedron() const
{ return poly_item->polyhedron(); }
QString Scene_edit_polyhedron_item::toolTip() const
{
    if(!poly_item->polyhedron())
        return QString();

    return QObject::tr("<p>Polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                       "<p>Number of vertices: %2<br />"
                       "Number of edges: %3<br />"
                       "Number of facets: %4</p>")
            .arg(this->name())
            .arg(poly_item->polyhedron()->size_of_vertices())
            .arg(poly_item->polyhedron()->size_of_halfedges()/2)
            .arg(poly_item->polyhedron()->size_of_facets())
            .arg(this->renderingModeName())
            .arg(this->color().name());
}
bool Scene_edit_polyhedron_item::isEmpty() const {
    return poly_item->isEmpty();
}
Scene_edit_polyhedron_item::Bbox Scene_edit_polyhedron_item::bbox() const {
    return poly_item->bbox();
}

void Scene_edit_polyhedron_item::setVisible(bool b) {
    poly_item->setVisible(b);
    Scene_item::setVisible(b);
    if(!b) {
        (*QGLViewer::QGLViewerPool().begin())->setManipulatedFrame(NULL);
    }
}
void Scene_edit_polyhedron_item::setColor(QColor c) {
    poly_item->setColor(c);
    Scene_item::setColor(c);
}
void Scene_edit_polyhedron_item::setName(QString n) {
    Scene_item::setName(n);
    n.replace(" (edit)", "");
    poly_item->setName(n);
}
void Scene_edit_polyhedron_item::setRenderingMode(RenderingMode m) {
    poly_item->setRenderingMode(m);
    Scene_item::setRenderingMode(m);
}
Scene_edit_polyhedron_item* Scene_edit_polyhedron_item::clone() const {
    return 0;
}
void Scene_edit_polyhedron_item::select(
        float orig_x,
        float orig_y,
        float orig_z,
        float dir_x,
        float dir_y,
        float dir_z)
{
    Scene_item::select(orig_x,
                       orig_y,
                       orig_z,
                       dir_x,
                       dir_y,
                       dir_z);
    poly_item->select(orig_x,
                      orig_y,
                      orig_z,
                      dir_x,
                      dir_y,
                      dir_z);
}

bool Scene_edit_polyhedron_item::keyPressEvent(QKeyEvent* e)
{
    //setting/unsetting rotation constraints
    if (e->key()==Qt::Key_R && !state.ctrl_pressing)
    {
        is_rot_free = !is_rot_free;
        rot_constraint.setRotationConstraintType( is_rot_free?
                                                      qglviewer::AxisPlaneConstraint::FREE:
                                                      qglviewer::AxisPlaneConstraint::AXIS);
        return true;
    }

    return false;
}

void Scene_edit_polyhedron_item::create_Sphere(float R)
{

    float T, P;
    float x[4],y[4],z[4];
    int rings = 22, sectors = 45;


    //Top of the sphere
    for(int t=0; t<360; t+=sectors)
    {

        pos_sphere.push_back(0);
        pos_sphere.push_back(0);
        pos_sphere.push_back(R);


        normals_sphere.push_back(0);
        normals_sphere.push_back(0);
        normals_sphere.push_back(1);



        P = rings*M_PI/180.0;
        T = t*M_PI/180.0;
        x[1] = sin(P) * cos(T) ;
        y[1] = sin(P) * sin(T) ;
        z[1] = cos(P);
        pos_sphere.push_back(R * x[1]);
        pos_sphere.push_back(R * y[1]);
        pos_sphere.push_back(R * z[1]);

        normals_sphere.push_back(x[1]);
        normals_sphere.push_back(y[1]);
        normals_sphere.push_back(z[1]);

        //
        P = rings*M_PI/180.0;
        T = (t+sectors)*M_PI/180.0;
        x[2] = sin(P) * cos(T) ;
        y[2] = sin(P) * sin(T) ;
        z[2] = cos(P);
        pos_sphere.push_back(R * x[2]);
        pos_sphere.push_back(R * y[2]);
        pos_sphere.push_back(R * z[2]);

        normals_sphere.push_back(x[2]);
        normals_sphere.push_back(y[2]);
        normals_sphere.push_back(z[2]);

    }

    //Body of the sphere
    for (int p=rings; p<180-rings; p+=rings)
        for(int t=0; t<360; t+=sectors)
        {
            //A
            P = p*M_PI/180.0;
            T = t*M_PI/180.0;
            x[0] = sin(P) * cos(T) ;
            y[0] = sin(P) * sin(T) ;
            z[0] = cos(P);

            pos_sphere.push_back(R * x[0]);
            pos_sphere.push_back(R * y[0]);
            pos_sphere.push_back(R * z[0]);

            normals_sphere.push_back(x[0]);
            normals_sphere.push_back(y[0]);
            normals_sphere.push_back(z[0]);

            //B
            P = (p+rings)*M_PI/180.0;
            T = t*M_PI/180.0;
            x[1] = sin(P) * cos(T) ;
            y[1] = sin(P) * sin(T) ;
            z[1] = cos(P);
            pos_sphere.push_back(R * x[1]);
            pos_sphere.push_back(R * y[1]);
            pos_sphere.push_back(R * z[1]);

            normals_sphere.push_back(x[1]);
            normals_sphere.push_back(y[1]);
            normals_sphere.push_back(z[1]);

            //C
            P = p*M_PI/180.0;
            T = (t+sectors)*M_PI/180.0;
            x[2] = sin(P) * cos(T) ;
            y[2] = sin(P) * sin(T) ;
            z[2] = cos(P);
            pos_sphere.push_back(R * x[2]);
            pos_sphere.push_back(R * y[2]);
            pos_sphere.push_back(R * z[2]);

            normals_sphere.push_back(x[2]);
            normals_sphere.push_back(y[2]);
            normals_sphere.push_back(z[2]);
            //D
            P = (p+rings)*M_PI/180.0;
            T = (t+sectors)*M_PI/180.0;
            x[3] = sin(P) * cos(T) ;
            y[3] = sin(P) * sin(T) ;
            z[3] = cos(P);
            pos_sphere.push_back(R * x[3]);
            pos_sphere.push_back(R * y[3]);
            pos_sphere.push_back(R * z[3]);

            normals_sphere.push_back(x[3]);
            normals_sphere.push_back(y[3]);
            normals_sphere.push_back(z[3]);



            pos_sphere.push_back(R * x[1]);
            pos_sphere.push_back(R * y[1]);
            pos_sphere.push_back(R * z[1]);

            normals_sphere.push_back(x[1]);
            normals_sphere.push_back(y[1]);
            normals_sphere.push_back(z[1]);

            pos_sphere.push_back(R * x[2]);
            pos_sphere.push_back(R * y[2]);
            pos_sphere.push_back(R * z[2]);

            normals_sphere.push_back(x[2]);
            normals_sphere.push_back(y[2]);
            normals_sphere.push_back(z[2]);

        }
    //Bottom of the sphere
    for(int t=0; t<360; t+=sectors)
    {


        pos_sphere.push_back(0);
        pos_sphere.push_back(0);
        pos_sphere.push_back(-R);

        normals_sphere.push_back(0);
        normals_sphere.push_back(0);
        normals_sphere.push_back(-1);


        P = (180-rings)*M_PI/180.0;
        T = t*M_PI/180.0;
        x[1] = sin(P) * cos(T) ;
        y[1] = sin(P) * sin(T) ;
        z[1] = cos(P);
        pos_sphere.push_back(R * x[1]);
        pos_sphere.push_back(R * y[1]);
        pos_sphere.push_back(R * z[1]);

        normals_sphere.push_back(x[1]);
        normals_sphere.push_back(y[1]);
        normals_sphere.push_back(z[1]);


        P = (180-rings)*M_PI/180.0;
        T = (t+sectors)*M_PI/180.0;
        x[2] = sin(P) * cos(T) ;
        y[2] = sin(P) * sin(T) ;
        z[2] = cos(P);
        pos_sphere.push_back(R * x[2]);
        pos_sphere.push_back(R * y[2]);
        pos_sphere.push_back(R * z[2]);

        normals_sphere.push_back(x[2]);
        normals_sphere.push_back(y[2]);
        normals_sphere.push_back(z[2]);

    }
}

//#include "Scene_edit_polyhedron_item.moc"
