#ifndef CSENTENCESSI_H
#define CSENTENCESSI_H

#include "csentence.h"
#include "kbGraph.h"
#include <string>
#include <vector>

namespace ukb {

    class CSentenceSSI; // forward declaration

    class CWordSSI : public CWord {
    public:

        CWordSSI() {
        };
        CWordSSI(const std::string & word, const std::string & id, char pos, cwtype type, float weight);
        CWordSSI(const std::string & word, const std::string & id, char pos, cwtype type, float weight, std::string vertex_string, Kb_vertex_t vertex);

    };

    class CSentenceSSI : public CSentence {
    public:

        CSentenceSSI() {
        };

        CSentenceSSI(const std::vector<CWordSSI> & cw_vec, const std::string id) {
            this->v = cw_vec;
            this->cs_id = id;
        }

        std::ostream & print_csent_simple_all(std::ostream & o, bool no_cnt, bool cnt_word) const;

    private:
        std::vector<CWordSSI> v;
    };
}
#endif
