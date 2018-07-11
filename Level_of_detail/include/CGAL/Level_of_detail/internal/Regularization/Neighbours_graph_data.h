#ifndef CGAL_LEVEL_OF_DETAIL_NEIGHBOURS_GRAPH_DATA_H
#define CGAL_LEVEL_OF_DETAIL_NEIGHBOURS_GRAPH_DATA_H

// STL includes.
#include <vector>

// Eigen includes.
#include <Eigen/SparseCore>

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel>
		class Neighbours_graph_data {

        public:
            using Kernel = InputKernel;
            using FT     = typename Kernel::FT;

            using FT_triplet  = Eigen::Triplet<FT>;
            using Int_triplet = Eigen::Triplet<int>;
            
            using Mus       = std::vector<FT_triplet>;
            using Targets   = std::vector<FT_triplet>;
            using Relations = std::vector<Int_triplet>;

            Neighbours_graph_data() : 
            m_mus(), 
            m_targets(), 
            m_relations() 
            { }

            inline Mus &get_mus() {
                return m_mus;
            }

            inline const Mus &get_mus() const {
                return m_mus;
            }

            inline Targets &get_targets() {
                return m_targets;
            }

            inline const Targets &get_targets() const {
                return m_targets;
            }

            inline Relations &get_relations() {
                return m_relations;
            }

            inline const Relations &get_relations() const {
                return m_relations;
            }

            inline size_t number_of_mus() const {
                return get_mus().size();
            }

            inline size_t number_of_targets() const {
                return get_targets().size();
            }

            inline size_t number_of_relations() const {
                return get_relations().size();
            }

            void clear() {

                m_mus.clear();
                m_targets.clear();
                m_relations.clear();
            }

            inline bool filled() const {
                
                const bool mat = (m_mus.size() != 0) && (m_targets.size() != 0);
                const bool rel = (m_relations.size() != 0);
                
                return (mat && rel) || mat;
            }

        private:
            Mus       m_mus;
            Targets   m_targets;
            Relations m_relations;
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_NEIGHBOURS_GRAPH_DATA_H