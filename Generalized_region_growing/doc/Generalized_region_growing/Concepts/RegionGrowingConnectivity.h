/*!
\ingroup PkgGeneralizedRegionGrowingConcepts
\cgalConcept

Concept describing the set of types and methods required by the class `CGAL::Region_growing::Generalized_region_growing`.

\cgalHasModel `CGAL::Region_growing::Region_growing_with_points::Points_connectivity_circular_query` `CGAL::Region_growing::Region_growing_with_points::Points_connectivity_nearest_neighbors` `CGAL::Region_growing::Region_growing_with_mesh::Mesh_connectivity`
*/

class RegionGrowingConnectivity {

public:

    /// Find all elements that have connection with `query_element` and push into `neighbors`
    template < class ElementWithProperties, class Neighbors_ >
    void get_neighbors(ElementWithProperties query_element, Neighbors_ neighbors) {
        
    }

};