#include "csentenceSSI.h"
#include "csentence.h"
#include "common.h"
#include "globalVars.h"
#include "kbGraph.h"
#include "wdict.h"

#include <string>
#include <vector>
#include <iosfwd>

namespace ukb {
    using namespace std;
    using namespace boost;

    CWordSSI::CWordSSI(const std::string & w_, const std::string & id_, char pos_, cwtype type_, float wght_) {

        w = w_;
        this->m_id = id_;
        this->m_pos = pos_;
        this->m_weight = wght_;
        this->m_type = type_;

        // Debug
        // cout << "Word: " << w << " + " << m_type << endl;

        switch (m_type) {
            case cw_concept:

                Kb_vertex_t u;
                bool P;

                tie(u, P) = ukb::Kb::instance().get_vertex_by_name(w);
                if (!P) {
                    throw std::runtime_error("CWord concept " + w + " not in KB. SSI");
                }
                m_syns.push_back(w);
                m_V.push_back(make_pair(u, 1.0f));
                m_ranks.push_back(0.0f);
                m_linkw_factor = 1.0;
                break;
            case cw_ctxword:
            case cw_tgtword:
            case cw_tgtword_nopv:
                if (!link_dict_concepts(w, m_pos)) {
                    // empty CWord
                    std::vector<std::string > ().swap(m_syns);
                }
                m_disamb = (1 == m_syns.size()); // monosemous words are disambiguated
                break;
            default:
                break;
        }
    }

    CWordSSI::CWordSSI(const std::string & word, const std::string & id, char pos, cwtype type, float weight, std::string vertex_string, Kb_vertex_t vertex) {
        if (type == cw_concept) {
            w = word;
            m_id = id;
            m_pos = pos;
            m_weight = weight;
            m_syns.push_back(vertex_string);
            m_V.push_back(make_pair(vertex, 1.0f));
            m_type = type;
            m_disamb = false;
        } else {
            w = word;
            m_id = id;
            m_pos = pos;
            m_weight = weight;
            m_syns.push_back(vertex_string);
            m_V.push_back(make_pair(vertex, 1.0f));
            m_type = type;
            m_disamb = true;
        }
    }

    std::ostream & CSentenceSSI::print_csent_simple_all(std::ostream & o, bool no_cnt, bool cnt_word) const {

        vector<CWordSSI>::const_iterator cw_it = v.begin();
        vector<CWordSSI>::const_iterator cw_end = v.end();

        for (; cw_it != cw_end; ++cw_it) {
            if (cw_it->size() == 0) continue;
            if (!cnt_word) {
                if (cw_it->is_cw_ctxword()) continue;
            }
            if (no_cnt) {
                if (!cw_it->is_disambiguated() && !glVars::output::ties) continue;
            }
            if (cw_it->is_monosemous() && !glVars::output::monosemous) return o; // Don't output monosemous words
            o << cs_id << " ";
            cw_it->print_cword_simple(o);
        }
        return o;
    }

}
