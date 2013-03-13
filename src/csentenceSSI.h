#ifndef CSENTENCESSI_H
#define CSENTENCESSI_H

#include "csentence.h"
#include <string>

namespace ukb {

    class CSentenceSSI; // forward declaration

    class CWordSSI : public CWord {
        
    public:
        CWordSSI(){};
        CWordSSI(const std::string & word, const std::string & id, char pos, cwtype type, float weight);
        CWordSSI(const std::string & word, const std::string & id, char pos, cwtype type, float weight, std::string vertex_string, Kb_vertex_t vertex);

        bool is_tgtword_nopv() const {
            return (m_type == cw_tgtword_nopv);
        }
        
        CWordSSI & operator=(const CWordSSI & cs_);

    private:
        std::vector<Kb_vertex_t> m_V2; // Esto no debería estar aqui
    };

    class CSentenceSSI : public CSentence {
        
    public:
        
        CSentenceSSI() {};
        CSentenceSSI(const std::vector<CWord> & cw_vec, const std::string id) {
            v = cw_vec;
            cs_id = id;
        }

        std::ostream & print_csent_simple_all(std::ostream & o, bool no_cnt, bool cnt_word) const;

    };
}
#endif
