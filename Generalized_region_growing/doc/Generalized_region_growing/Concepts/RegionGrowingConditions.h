/*!
\ingroup PkgGeneralizedRegionGrowingConcepts
\cgalConcept

Concept describing the set of types and methods required by the class `CGAL::Region_growing::Generalized_region_growing`.

\cgalHasModel `CGAL::Region_growing::Region_growing_with_points::Points_conditions_2` `CGAL::Region_growing::Region_growing_with_points::Points_conditions_3`
`CGAL::Region_growing::Region_growing_with_mesh::Mesh_conditions`
*/

class RegionGrowingConditions {

public:

    /// Local condition to check if two elements `unassigned_elem` and `assigned_elem` are similar to be put in a region
    template < class ElementWithProperties, class Region_ >
    bool is_in_same_region(ElementWithProperties unassigned_elem, ElementWithProperties assigned_elem, Region_ region) {
        
    }

    /// Refit the shape to the region function
    template < class Region_ >
    void update_shape(Region_ region) {

    }

    /// Global condition to check the validity of the region
    template < class Region_ >
    bool is_valid(Region_ region) {

    }

};