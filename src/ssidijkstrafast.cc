// UKB stuff and more

#include "wdict.h"
#include "common.h"
#include "fileElem.h"
#include "globalVars.h"
#include "kbGraph.h"
#include "disambGraph.h"
#include "csentenceSSI.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <ostream>

// BoostGraph stuff

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp> // Dijkstra's Shortest Paths Algorithm
#include <boost/program_options.hpp> // Program options

// Time control stuff

#include <time.h>
#include <ctime>
#include <unistd.h>
#include <sys/time.h>


// Usage of namespaces

using namespace std;
using namespace boost;
using namespace ukb;

// Global variables

float **dijkstra_matrix; // Matrix in which we will store the distance between different Kb_vertex_t
int rows; // Number of rows of the Dijkstra matrix

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Matrix manipulation methods (SSI-DijkstraFast)

void create_distance_matrix(vector<Kb_vertex_t> I, vector<CWordSSI> P) {

    int number_of_rows = 0;
    int number_of_columns = 0;

    unsigned int isize = I.size();
    unsigned int psize = P.size();

    //Getting the maximum size of the distance matrix (sum of P senses)
    number_of_columns = isize + psize - 1;
    for (unsigned int j = 0; j < psize; j++) {
        vector<string> synset_vector = P.at(j).get_syns_vector();
        number_of_rows = synset_vector.size() + number_of_rows;
    }

    rows = number_of_rows;

    dijkstra_matrix = new float * [number_of_rows];
    for (int i = 0; i < number_of_rows; i++) {
        dijkstra_matrix[i] = new float[number_of_columns];
    }

    Kb_vertex_t previous_vertex = 0;
    for (unsigned int i = 0; i < isize; i++) {
        Kb_vertex_t actual_isynset = I.at(i);
        int r = 0;
        for (unsigned int j = 0; j < psize; j++) {
            CWordSSI actual_word = P.at(j);
            vector<string> synset_vector = actual_word.get_syns_vector();
            for (unsigned int k = 0; k < synset_vector.size(); k++) {
                Kb_vertex_t actual_synset;
                bool x;
                tie(actual_synset, x) = Kb::instance().get_vertex_by_name(synset_vector.at(k));
                if (x) {
                    dijkstra_matrix[r][i] = Kb::instance().obtain_distance_dijsktra_faster(actual_isynset, actual_synset, previous_vertex);
                    // cout << "Create: [" << r << "," << i << "]" << endl;
                    r++;
                    previous_vertex = actual_isynset;
                }
            }
        }
    }
}

void update_distance_matrix(Kb_vertex_t src, vector<CWordSSI> P, int rowm, int columnm) {

    Kb_vertex_t previous_vertex = 0;
    int modifing_row = rowm;
    int modifing_column = columnm;

    for (unsigned int i = 0; i < P.size(); i++) {
        CWordSSI actual_word = P.at(i);
        vector<string> synset_vector = actual_word.get_syns_vector();
        for (unsigned int k = 0; k < synset_vector.size(); k++) {
            Kb_vertex_t actual_synset;
            bool x;
            tie(actual_synset, x) = Kb::instance().get_vertex_by_name(synset_vector.at(k));
            if (x) {
                dijkstra_matrix[modifing_row][modifing_column] = Kb::instance().obtain_distance_dijsktra_faster(src, actual_synset, previous_vertex);
                // cout << "Update: [" << modifing_row << "," << modifing_column << "]" << endl;
                modifing_row++;
            }
        }
    }
}

void delete_distance_matrix() {
    for (int i = 0; i < rows; i++) {
        delete [] dijkstra_matrix[i];
    }
    delete [] dijkstra_matrix;
}

void show_path(Kb_vertex_t v1, Kb_vertex_t v2, const KbGraph & g) {

    vector<Kb_vertex_t> parents = Kb::instance().get_parents();
    graph_traits<KbGraph>::vertex_iterator vi, vend;
    vector<Kb_vertex_t> vertex_vector;
    for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) {
        if (*vi == v2) //Chosen vertex is our destination vertex (v2)
        {
            Kb_vertex_t auxiliar_vertex = v2;
            unsigned int vertex_index = *vi;
            while (auxiliar_vertex != v1) {
                vertex_vector.push_back(auxiliar_vertex);
                auxiliar_vertex = parents[vertex_index];
                for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) {
                    if (*vi == auxiliar_vertex) {
                        vertex_index = *vi;
                        break;
                    }
                }
            }
            vertex_vector.push_back(v1);
            break;
        }
    }
    cout << "Path: ";
    //	Printing backwards the vertices vector is for getting the way from v1 to v2 (v1 -> v2) and NOT the way from v2 to v1 (v2-> v1)
    for (int i = vertex_vector.size() - 1; i >= 1; i--) {
        cout << Kb::instance().get_vertex_name(vertex_vector.at(i)) << " -> ";
    }
    cout << Kb::instance().get_vertex_name(vertex_vector.at(0));
    cout << "\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SSI-Dijkstra, SSI-Dijkstra+ & SSI-Dijkstra-Fast

CSentenceSSI ssi_dijkstra_fast(CSentenceSSI cs, int option) {

    vector<Kb_vertex_t> interpreted; // I
    vector<CWordSSI> pending; // P
    vector<CWordSSI> solution_vector; // Auxiliar CWord for solution


    int row_modifier = 0; // How many matrix rows we already used
    int column_modifier = 0; // In which column we need to update the matrix info next

    CSentence::iterator it;
    CSentence::iterator end = cs.end();
    //@Aritza: Variable que usaremos para saber que elemento del Csentences hemos añadido ya a la lista de interpreted y borrarlo de los pending 	
    int posicion;

    string id = cs.id();

    // First step - part one: Initializing main data structs
    for (it = cs.begin(); it != end; ++it) {
        vector<string> syns = it->get_syns_vector();
        int senses = syns.size();
        if (senses == 0) // No senses for this word, dude!
        {
            continue; // As this word doesn't is in KbGraph, we continue processing other words
        } else {
            if (senses == 1) // Monosemical & disambiguated
            {
                Kb_vertex_t aux_vertex;
                bool u;
                tie(aux_vertex, u) = Kb::instance().get_vertex_by_name(syns.at(0));
                if (u) {
                    CWordSSI aux_cword;
                    //@Aritza: Cambiar a la nueva version con el type()
                    aux_cword = CWordSSI(it->word(), it->id(), it->get_pos(), it->type(), it->get_weight(), syns.at(0), aux_vertex);
                    //@Aritza: Creamos un if para el nuevo tipo de palabra 3, ya que estas no tienen q entrar en el verctor de interpreted y puede ser monosemica
                    if (!it->is_tgtword_nopv()) {
                        interpreted.push_back(aux_vertex);
                        solution_vector.push_back(aux_cword);
                    } else {
                        pending.push_back(aux_cword);
                    }
                }
            } else // Polisemical... pending disambiguation 
            {
                CWordSSI aux_cword;
                //@Aritza: Cambiar a la nueva version con el type()
                aux_cword = CWordSSI(it->word(), it->id(), it->get_pos(), it->type(), it->get_weight());
                pending.push_back(aux_cword);
            }
        }
    }

    // First step - part two: If there are no monosemical words, pass one pair <word, sense> from P to I using FSI

    //@Aritza: Cambiamos if por for por si la 1º o las siguientes palabras llevan un 3, lo recorre asta q encuentra una palabra q no sea 3

    for (unsigned int o = 0; interpreted.empty() && o < pending.size(); o++) {
        CWordSSI word_for_disambiguate = pending.at(o);
        vector<string> the_synset_vector = word_for_disambiguate.get_syns_vector();
        Kb_vertex_t BestSense;
        float MaxValue = 0.0;
        int index = 0;
        for (unsigned int j = 0; j < the_synset_vector.size(); j++) {
            Kb_vertex_t actual_synset;
            string actual_synset_word = the_synset_vector.at(j);
            bool v;
            tie(actual_synset, v) = Kb::instance().get_vertex_by_name(actual_synset_word);
            if (v) {
                float total_weight = 0.0;
                float d = 0.0;
                int number_of_synsets = 0;
                Kb_vertex_t previous_vertex = 0;
                for (unsigned int k = 1; k < pending.size(); k++) {
                    CWordSSI pending_actual_word = pending.at(k);
                    vector<string> p_synset_vector = pending_actual_word.get_syns_vector();
                    string pending_string_synset = p_synset_vector.at(0);
                    Kb_vertex_t actual_psynset;
                    bool w;
                    tie(actual_psynset, w) = Kb::instance().get_vertex_by_name(pending_string_synset); // Pending word FIRST synset
                    if (w) {
                        float actual_weight = 0.0;

                        if (option) {
                            actual_weight = Kb::instance().obtain_distance_dijsktra_faster(actual_synset, actual_psynset, previous_vertex);
                            previous_vertex = actual_synset;
                        } else {
                            actual_weight = Kb::instance().obtain_distance_dijsktra(actual_synset, actual_psynset);
                        }

                        d = actual_weight + 1.0;
                        total_weight = total_weight + (1 / d);
                        number_of_synsets++;
                    }
                }
                //@Aritza: Añadimos la condicion para que solo las palabras que no sean de tipo 3 se tengan en cuenta
                if (number_of_synsets > 0 && !word_for_disambiguate.is_tgtword_nopv()) {
                    float NewValue = 0.0;
                    NewValue = total_weight / number_of_synsets;
                    if (NewValue > MaxValue) {
                        MaxValue = NewValue;
                        BestSense = actual_synset;
                        index = j;
                        posicion = o;
                    }
                }
            }
        }
        if (MaxValue > 0) {
            interpreted.push_back(BestSense);
            CWordSSI resolved;
            //@Aritza: Cambiar a la nueva version con el type()
            resolved = CWordSSI(word_for_disambiguate.word(), word_for_disambiguate.id(), word_for_disambiguate.get_pos(), word_for_disambiguate.type(), word_for_disambiguate.get_weight(), the_synset_vector.at(index), BestSense);
            solution_vector.push_back(resolved);
            //@Aritza: Ahora no tenemos que eliminar la primera, sino la que hemos metido en el interpreted
            pending.erase(pending.begin() + posicion);
        }
    }
    // End FSI

    // Second step: Disambiguation

    //@Aritza: Si todas las palabras son de tipo 3 es imposible desambiguar el contexto
    if (interpreted.empty()) cout << "Todas las palabras son de tipo: target 'nopv' word (3), por tanto no es posible desambiguar el contexto: " << cs.id() << endl;

    if (!interpreted.empty()) // If there is only one word, interpreted will be empty... 
    {
        //We'll try to avoid calling Dijkstra too much by storing distances on a matrix
        create_distance_matrix(interpreted, pending);
        column_modifier = interpreted.size();

        CWordSSI actual_word;
        while (!pending.empty()) {
            actual_word = pending.at(0);
            vector<string> the_synset_vector = actual_word.get_syns_vector();
            Kb_vertex_t BestSense;
            float MaxValue = 0.0;
            int index = 0;
            for (unsigned int j = 0; j < the_synset_vector.size(); j++) {
                unsigned int row = j + row_modifier;
                string actual_synset_word = the_synset_vector.at(j);
                Kb_vertex_t actual_synset;
                bool x;
                tie(actual_synset, x) = Kb::instance().get_vertex_by_name(actual_synset_word);
                if (x) {
                    float total_weight = 0.0;
                    float d = 0.0;
                    int number_of_synsets = 0;
                    Kb_vertex_t previous_vertex = 0;
                    for (unsigned int k = 0; k < interpreted.size(); k++) {
                        Kb_vertex_t actual_isynset; // Already interpreted synset
                        actual_isynset = interpreted.at(k);
                        float actual_weight = 0.0;
                        if (option) {
                            actual_weight = dijkstra_matrix[row][k];
                            // cout << "Operate: [" << row << "," << k << "]" << endl;
                        } else {
                            actual_weight = Kb::instance().obtain_distance_dijsktra(actual_synset, actual_isynset);
                        }

                        d = actual_weight + 1.0;
                        total_weight = total_weight + (1 / d);
                        number_of_synsets++;
                    }

                    // Condition: if PoS of actual processing word is 'v' then use FSP for optimize answer

                    if (actual_word.get_pos() == 'v') {
                        for (unsigned int l = 1; l < pending.size(); l++) {
                            CWordSSI pending_actual_word = pending.at(l);
                            vector<string> p_synset_vector = pending_actual_word.get_syns_vector();
                            string pending_string_synset = p_synset_vector.at(0);
                            Kb_vertex_t actual_psynset;
                            bool y;
                            tie(actual_psynset, y) = Kb::instance().get_vertex_by_name(pending_string_synset); // Pending word FIRST synset
                            if (y) {
                                float actual_weight = 0.0;
                                float d = 0.0;

                                if (option) {
                                    actual_weight = Kb::instance().obtain_distance_dijsktra_faster(actual_synset, actual_psynset, previous_vertex);
                                    previous_vertex = actual_synset;
                                } else {
                                    actual_weight = Kb::instance().obtain_distance_dijsktra(actual_synset, actual_psynset);
                                }

                                d = actual_weight + 1.0;
                                total_weight = total_weight + (1 / d);
                                number_of_synsets++;
                            }
                        }
                    }
                    // End FSP	

                    if (number_of_synsets > 0) {
                        float NewValue = 0.0;
                        NewValue = total_weight / number_of_synsets;
                        if (NewValue > MaxValue) {
                            MaxValue = NewValue;
                            BestSense = actual_synset;
                            index = j;
                        }
                    }
                }
            }
            if (MaxValue > 0) {
                //@Aritza: Creamos un if para el nuevo tipo de palabra 3, ya que estas no tienen q entrar en el vector de interpreted
                if (!actual_word.is_tgtword_nopv()) {
                    interpreted.push_back(BestSense);
                }
                CWordSSI resolved;
                //@Aritza: Cambiar a la nueva version con el type()
                resolved = CWordSSI(actual_word.word(), actual_word.id(), actual_word.get_pos(), actual_word.type(), actual_word.get_weight(), the_synset_vector.at(index), BestSense);
                solution_vector.push_back(resolved);
                pending.erase(pending.begin());

                // Udpate matrix info with new interpreted item (only if pending has more items)
                if (!pending.empty()) {
                    row_modifier = row_modifier + the_synset_vector.size();
                    update_distance_matrix(BestSense, pending, row_modifier, column_modifier);
                    column_modifier++;
                }
            }
        }
    }
    delete_distance_matrix();
    return CSentenceSSI(solution_vector, id);
}

CSentenceSSI ssi_dijkstra_plus(CSentenceSSI cs, int option) {

    vector<Kb_vertex_t> interpreted; // I
    vector<CWordSSI> pending; // P
    vector<CWordSSI> solution_vector; // Auxiliar CWord for solution

    CSentence::iterator it;
    CSentence::iterator end = cs.end();
    //@Aritza: Variable que usaremos para saber que elemento del Csentences hemos añadido ya a la lista de interpreted y borrarlo de los pending 	
    int posicion;

    string id = cs.id();

    // First step - part one: Initializing main data structs
    for (it = cs.begin(); it != end; ++it) {
        vector<string> syns = it->get_syns_vector();
        int senses = syns.size();
        if (senses == 0) // No senses for this word, dude!
        {
            continue; // As this word doesn't is in KbGraph, we continue processing other words
        } else {
            if (senses == 1) // Monosemical & disambiguated
            {
                Kb_vertex_t aux_vertex;
                bool u;
                tie(aux_vertex, u) = Kb::instance().get_vertex_by_name(syns.at(0));
                if (u) {
                    CWordSSI aux_cword;
                    //@Aritza: Cambiar a la nueva version con el type()
                    aux_cword = CWordSSI(it->word(), it->id(), it->get_pos(), it->type(), it->get_weight(), syns.at(0), aux_vertex);
                    //@Aritza: Creamos un if para el nuevo tipo de palabra 3, ya que estas no tienen q entrar en el verctor de interpreted y puede ser monosemica
                    if (!it->is_tgtword_nopv()) {
                        interpreted.push_back(aux_vertex);
                        solution_vector.push_back(aux_cword);
                    } else {
                        pending.push_back(aux_cword);
                    }
                }
            } else // Polisemical... pending disambiguation 
            {
                CWordSSI aux_cword;
                //@Aritza: Cambiar a la nueva version con el type()
                aux_cword = CWordSSI(it->word(), it->id(), it->get_pos(), it->type(), it->get_weight());
                pending.push_back(aux_cword);
            }
        }
    }

    // First step - part two: If there are no monosemical words, pass one pair <word, sense> from P to I using FSI

    //@Aritza: Cambiamos if por for por si la 1º o las siguientes palabras llevan un 3, lo recorre asta q encuentra una palabra q no sea 3

    for (unsigned int o = 0; interpreted.empty() && o < pending.size(); o++) {
        CWordSSI word_for_disambiguate = pending.at(o);
        vector<string> the_synset_vector = word_for_disambiguate.get_syns_vector();
        Kb_vertex_t BestSense;
        float MaxValue = 0.0;
        int index = 0;
        for (unsigned int j = 0; j < the_synset_vector.size(); j++) {
            Kb_vertex_t actual_synset;
            string actual_synset_word = the_synset_vector.at(j);
            bool v;
            tie(actual_synset, v) = Kb::instance().get_vertex_by_name(actual_synset_word);
            if (v) {
                float total_weight = 0.0;
                float d = 0.0;
                int number_of_synsets = 0;
                Kb_vertex_t previous_vertex = 0;
                for (unsigned int k = 1; k < pending.size(); k++) {
                    CWordSSI pending_actual_word = pending.at(k);
                    vector<string> p_synset_vector = pending_actual_word.get_syns_vector();
                    string pending_string_synset = p_synset_vector.at(0);
                    Kb_vertex_t actual_psynset;
                    bool w;
                    tie(actual_psynset, w) = Kb::instance().get_vertex_by_name(pending_string_synset); // Pending word FIRST synset
                    if (w) {
                        float actual_weight = 0.0;

                        if (option) {
                            actual_weight = Kb::instance().obtain_distance_dijsktra_faster(actual_synset, actual_psynset, previous_vertex);
                            previous_vertex = actual_synset;
                        } else {
                            actual_weight = Kb::instance().obtain_distance_dijsktra(actual_synset, actual_psynset);
                        }

                        d = actual_weight + 1.0;
                        total_weight = total_weight + (1 / d);
                        number_of_synsets++;
                    }
                }
                //@Aritza: Añadimos la condicion para que solo las palabras que no sean de tipo 3 se tengan en cuenta
                if (number_of_synsets > 0 && !word_for_disambiguate.is_tgtword_nopv()) {
                    float NewValue = 0.0;
                    NewValue = total_weight / number_of_synsets;
                    if (NewValue > MaxValue) {
                        MaxValue = NewValue;
                        BestSense = actual_synset;
                        index = j;
                        posicion = o;
                    }
                }
            }
        }
        if (MaxValue > 0) {
            interpreted.push_back(BestSense);
            CWordSSI resolved;
            //@Aritza: Cambiar a la nueva version con el type()
            resolved = CWordSSI(word_for_disambiguate.word(), word_for_disambiguate.id(), word_for_disambiguate.get_pos(), word_for_disambiguate.type(), word_for_disambiguate.get_weight(), the_synset_vector.at(index), BestSense);
            solution_vector.push_back(resolved);
            //@Aritza: Ahora no tenemos que eliminar la primera, sino la que hemos metido en el interpreted
            pending.erase(pending.begin() + posicion);
        }
    }
    // End FSI

    // Second step: Disambiguation

    //@Aritza: Si todas las palabras son de tipo 3 es imposible desambiguar el contexto
    if (interpreted.empty()) cout << "Todas las palabras son de tipo: target 'nopv' word (3), por tanto no es posible desambiguar el contexto: " << cs.id() << endl;

    if (!interpreted.empty()) // If there is only one word, interpreted will be empty... 
    {
        CWordSSI actual_word;
        while (!pending.empty()) {
            actual_word = pending.at(0);
            vector<string> the_synset_vector = actual_word.get_syns_vector();
            Kb_vertex_t BestSense;
            float MaxValue = 0.0;
            int index = 0;
            for (unsigned int j = 0; j < the_synset_vector.size(); j++) {
                string actual_synset_word = the_synset_vector.at(j);
                Kb_vertex_t actual_synset;
                bool x;
                tie(actual_synset, x) = Kb::instance().get_vertex_by_name(actual_synset_word);
                if (x) {
                    float total_weight = 0.0;
                    float d = 0.0;
                    int number_of_synsets = 0;
                    Kb_vertex_t previous_vertex = 0;
                    for (unsigned int k = 0; k < interpreted.size(); k++) {
                        Kb_vertex_t actual_isynset; // Already interpreted synset
                        actual_isynset = interpreted.at(k);
                        float actual_weight = 0.0;

                        if (option) {
                            actual_weight = Kb::instance().obtain_distance_dijsktra_faster(actual_synset, actual_isynset, previous_vertex);
                            previous_vertex = actual_synset;
                            /* //Debug
                            if (actual_word.word().compare("together") == 0){
                                show_path(actual_synset, actual_isynset, Kb::instance().graph());
                            } */
                        } else {
                            actual_weight = Kb::instance().obtain_distance_dijsktra(actual_synset, actual_isynset);
                        }

                        d = actual_weight + 1.0;
                        total_weight = total_weight + (1 / d);
                        number_of_synsets++;
                    }

                    // Condition: if PoS of actual processing word is 'v' then use FSP for optimize answer

                    if (actual_word.get_pos() == 'v') {
                        for (unsigned int l = 1; l < pending.size(); l++) {
                            CWordSSI pending_actual_word = pending.at(l);
                            vector<string> p_synset_vector = pending_actual_word.get_syns_vector();
                            string pending_string_synset = p_synset_vector.at(0);
                            Kb_vertex_t actual_psynset;
                            bool y;
                            tie(actual_psynset, y) = Kb::instance().get_vertex_by_name(pending_string_synset); // Pending word FIRST synset
                            if (y) {
                                float actual_weight = 0.0;
                                float d = 0.0;

                                if (option) {
                                    actual_weight = Kb::instance().obtain_distance_dijsktra_faster(actual_synset, actual_psynset, previous_vertex);
                                    previous_vertex = actual_synset;
                                } else {
                                    actual_weight = Kb::instance().obtain_distance_dijsktra(actual_synset, actual_psynset);
                                }

                                d = actual_weight + 1.0;
                                total_weight = total_weight + (1 / d);
                                number_of_synsets++;
                            }
                        }
                    }
                    // End FSP	

                    if (number_of_synsets > 0) {
                        float NewValue = 0.0;
                        NewValue = total_weight / number_of_synsets;
                        if (NewValue > MaxValue) {
                            MaxValue = NewValue;
                            BestSense = actual_synset;
                            index = j;
                        }
                    }
                }
            }
            if (MaxValue > 0) {
                //@Aritza: Creamos un if para el nuevo tipo de palabra 3, ya que estas no tienen q entrar en el vector de interpreted
                if (!actual_word.is_tgtword_nopv()) interpreted.push_back(BestSense);
                CWordSSI resolved;
                //@Aritza: Cambiar a la nueva version con el type()
                resolved = CWordSSI(actual_word.word(), actual_word.id(), actual_word.get_pos(), actual_word.type(), actual_word.get_weight(), the_synset_vector.at(index), BestSense);
                solution_vector.push_back(resolved);
                pending.erase(pending.begin());
            }
        }
    }
    return CSentenceSSI(solution_vector, id);
}

CSentenceSSI ssi_dijkstra(CSentenceSSI cs, int option) {

    vector<Kb_vertex_t> interpreted; // I
    vector<CWordSSI> pending; // P
    vector<CWordSSI> solution_vector; // Auxiliar CWord for solution

    CSentenceSSI::iterator it;
    CSentenceSSI::iterator end = cs.end();

    string id = cs.id();

    // First step: Initializing main data structs
    for (it = cs.begin(); it != end; ++it) {
        vector<string> syns = it->get_syns_vector();
        int senses = syns.size();
        if (senses == 0) // No senses for this word, dude!
        {
            continue; // As this word doesn't is in KbGraph, we continue processing other words
        } else {
            if (senses == 1) // Monosemical & disambiguated
            {
                Kb_vertex_t aux_vertex;
                bool a;
                tie(aux_vertex, a) = Kb::instance().get_vertex_by_name(syns.at(0));
                if (a) {
                    CWordSSI aux_cword;
                    //@Aritza: Cambiar el m_type al de la nueva version
                    aux_cword = CWordSSI(it->word(), it->id(), it->get_pos(), it->type(), it->get_weight(), syns.at(0), aux_vertex);
                    //@Aritza: Creamos un if para el nuevo tipo de palabra 3, ya que estas no tienen q entrar en el verctor de interpreted y puede ser monosemica. Si es monosemica la introducimos directamente en el vector pending

                    if (!it->is_tgtword_nopv()) {
                        interpreted.push_back(aux_vertex);
                        solution_vector.push_back(aux_cword);
                    } else {
                        pending.push_back(aux_cword);
                    }
                }
            } else // Polisemical... pending disambiguation
            {
                CWordSSI aux_cword;
                //@Aritza: Cambiar el m_type al de la nueva version, incluyendo el it->type()
                aux_cword = CWordSSI(it->word(), it->id(), it->get_pos(), it->type(), it->get_weight());
                pending.push_back(aux_cword);
            }
        }
    }

    //Aritza: Si el contexto no va a poder ser desambiguado mostramos un mensaje por pantalla indicandolo
    if (interpreted.empty()) cout << "No hay ninguna palabra monosémica ni ningún concepto, por tanto no es posible desambiguar el contexto: " << cs.id() << endl;

    // Second step: Disambiguation
    if (!interpreted.empty()) {
        CWordSSI actual_word;
        while (!pending.empty()) {
            actual_word = pending.at(0);
            vector<string> the_synset_vector = actual_word.get_syns_vector();
            Kb_vertex_t BestSense;
            float MaxValue = 0.0;
            int index = 0;
            for (unsigned int j = 0; j < the_synset_vector.size(); j++) {
                string actual_synset_word = the_synset_vector.at(j);
                Kb_vertex_t actual_synset;
                bool z;
                tie(actual_synset, z) = Kb::instance().get_vertex_by_name(actual_synset_word);
                if (z) {
                    float total_weight = 0.0;
                    float d = 0.0;
                    int number_of_synsets = 0;
                    Kb_vertex_t previous_vertex = 0;
                    for (unsigned int k = 0; k < interpreted.size(); k++) {
                        Kb_vertex_t actual_isynset; // Already interpreted synset
                        actual_isynset = interpreted.at(k);
                        float actual_weight = 0.0;

                        if (option) {
                            actual_weight = Kb::instance().obtain_distance_dijsktra_faster(actual_synset, actual_isynset, previous_vertex);
                            previous_vertex = actual_synset;
                        } else {
                            actual_weight = Kb::instance().obtain_distance_dijsktra(actual_synset, actual_isynset);
                        }

                        d = actual_weight + 1.0;
                        total_weight = total_weight + (1 / d);
                        number_of_synsets++;
                    }
                    if (number_of_synsets > 0) {
                        float NewValue = 0.0;
                        NewValue = total_weight / number_of_synsets;
                        if (NewValue > MaxValue) {
                            MaxValue = NewValue;
                            BestSense = actual_synset;
                            index = j;
                        }
                    }
                }
            }
            if (MaxValue > 0) {
                //@Aritza: Creamos un if para el nuevo tipo de palabra 3, ya que estas no tienen q entrar en el vector de interpreted
                if (!actual_word.is_tgtword_nopv()) interpreted.push_back(BestSense);
                CWordSSI resolved;
                //@Aritza: Cambiar el m_type al de la nueva version
                resolved = CWordSSI(actual_word.word(), actual_word.id(), actual_word.get_pos(), actual_word.type(), actual_word.get_weight(), the_synset_vector.at(index), BestSense);
                solution_vector.push_back(resolved);
                pending.erase(pending.begin());
            }
        }
    }
    return CSentenceSSI(solution_vector, id);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Sorting stuff

struct polysemy_sorting {

    bool operator() (CWord cw1, CWord cw2) {
        vector<string> syns1 = cw1.get_syns_vector();
        vector<string> syns2 = cw2.get_syns_vector();

        bool value_to_return = (syns1.size() < syns2.size());

        return value_to_return;
    }
} sort_by_polysemy_degree;

struct explicit_sorting {

    bool operator() (CWord cw1, CWord cw2) {
        float w1 = cw1.get_weight();
        float w2 = cw2.get_weight();

        bool value_to_return = (w1 < w2);

        return value_to_return;
    }
} sort_by_explicit_sorting;

struct sorting_by_id {

    bool operator() (CWord cw1, CWord cw2) {
        string id1 = cw1.id();
        string id2 = cw2.id();

        bool value_to_return = (id1 < id2);

        return value_to_return;
    }
} id_sorting;

void sort_words(CSentence & cs, int sort_option) {

    switch (sort_option) {
        case 0: // Without sorting
            break;

        case 1: std::sort(cs.begin(), cs.end(), sort_by_polysemy_degree); // Sorting by polysemy degree
            break;

        case 2: std::sort(cs.begin(), cs.end(), sort_by_explicit_sorting); // Explicit sorting
            break;

        default: // Impossible...
            break;
    }
}

void sort_by_id(CSentence & cs) {
    sort(cs.begin(), cs.end(), id_sorting);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {

    // Name of the binfile. Remember to put the binfile in the correct path (that is same of the executable file).

    const string kb_default_binfile = "kb_wnet_weights.bin"; // Graph with weights for default

    string kb_binfile(kb_default_binfile);

    // -- boost::program_options --

    // Sorting option

    enum sorting {
        no_explicit,
        polysemy_degree,
        explicit_sorting
    };

    // Dijkstra Shortest Path option

    enum dijkstra_option {
        normal,
        faster
    };

    // Dijkstra Algorithm

    enum dijkstra_method {
        basic,
        plus,
        fast
    };

    sorting sort_method = no_explicit;
    dijkstra_option dsp_option = faster;
    dijkstra_method dij_option = fast;

    bool pretty_mode = false;
    //@Aritza:añadinos la nueva variable para especificar que no queremos mostros las palabras ya desambiguadas	
    bool no_cnt = false;

    bool cnt_word = false;

    string cmdline("!! -v ");
    cmdline += glVars::ukb_version;
    for (int i = 0; i < argc; ++i) {
        cmdline += " ";
        cmdline += argv[i];
    }

    vector<string> input_files;
    string fullname_in;

    using namespace boost::program_options;

    const char desc_header[] = "SSI_Dijkstra: Here comes a description";

    options_description po_general("General");

    po_general.add_options()
            ("help,h", "This page")
            ("version", "Show version")
            ("kb_binfile,K", value<string > (), "Binary file of KB (see compile_kb). Default is kb_wnet_weights.bin")
            ("dict_file,D", value<string > (), "Dictionary text file. Default is dict.txt")
            ("pretty", "Use this option for sort the words by id before disambiguation")
            ;

    options_description po_sorting("Sorting");
    po_sorting.add_options()
            ("nexp", "No explicit sorting. This is default option")
            ("poly", "Sorting by polysemy degree")
            ("expl", "Sorting by explicit sorting (given by user)")
            ;

    options_description po_method("SSI Dijkstra method");
    po_method.add_options()
            ("basic", "Basic SSI Dijkstra method")
            ("plus", "SSI Dijkstra plus method")
            ("fast", "SSI Dijstra fast method")
            ;

    options_description po_ssid("SSI Dijkstra options");
    po_ssid.add_options()
            ("no-opt", "Use a normal invocation to dijkstra_shortest_paths instead of faster one")
            ("ssid", "Use SSI-Dijkstra algorithm instead of SSI-Dijkstra+")
            //@Aritza: Añadimos una nueva opcion por si no queremos que se muestren los conceptos
            ("no_cnt", "No muestra el cnt")
            //@Aritza: Añadimos una nueva opcion por si queremos que las palabras de tipo cw_ctxword se muestren
            ("cnt_word", "Mostrar los cnt word");
    ;

    options_description po_hidden("Hidden");
    po_hidden.add_options()
            ("input_file", value<string > (), "Input file")
            ;

    options_description po_visible(desc_header);
    po_visible.add(po_general).add(po_sorting).add(po_method).add(po_ssid);

    options_description po_all("All options");
    po_all.add(po_visible).add(po_hidden);

    positional_options_description po_optdesc;
    po_optdesc.add("input_file", 1);

    try {

        variables_map vm;
        store(command_line_parser(argc, argv).options(po_all).positional(po_optdesc).run(), vm);

        notify(vm);

        // If asked for help don't do anything more

        // Program options: General

        if (vm.count("help")) {
            cout << po_visible << endl;
            exit(0);
        }

        if (vm.count("version")) {
            cout << glVars::ukb_version << endl;
            exit(0);
        }

        if (vm.count("kb_binfile")) {
            kb_binfile = vm["kb_binfile"].as<string > ();
        }

        if (vm.count("dict_file")) {
            glVars::dict_filename = vm["dict_file"].as<string > ();
        }

        if (vm.count("pretty")) {
            pretty_mode = true;
        }

        // Program options: Sorting

        if (vm.count("nexp")) {
            sort_method = no_explicit;
        }

        if (vm.count("poly")) {
            sort_method = polysemy_degree;
        }

        if (vm.count("expl")) {
            sort_method = explicit_sorting;
        }

        // Program options: SSI Dijkstra method

        if (vm.count("basic")) {
            dij_option = basic;
        }

        if (vm.count("plus")) {
            dij_option = plus;
        }

        if (vm.count("fast")) {
            dij_option = fast;
        }

        // Program options: SSI Dijkstra 

        if (vm.count("no-opt")) {
            dsp_option = normal;
        }

        //@Aritza: Additional options

        if (vm.count("no_cnt")) no_cnt = true;
        if (vm.count("cnt_word")) cnt_word = true;

        // Input file

        if (vm.count("input_file")) {
            fullname_in = vm["input_file"].as<string > ();
        }
    } catch (std::exception& e) {
        cerr << e.what() << "\n";
        throw (e);
    }

    if (!fullname_in.size()) {
        cout << po_visible << endl;
        cout << "Error: No input" << endl;
        exit(-1);
    }

    // Reading graph from binfile

    Kb::create_from_binfile(kb_binfile);
    cout << cmdline << endl;

    // Reading CSentence (context...)

    ifstream fh_in(fullname_in.c_str());
    if (!fh_in) {
        cerr << "Can't open " << fullname_in << endl;
        exit(-1);
    }

    CSentenceSSI cs;
    // Adding dictionary 
    /* 
       if (insert_all_dict) {
        //@Aritza: Ya no es necesario pasarle el peso
        kb.add_dictionary(); //TODO: Corregir referencia sin definir
    } 
     */

    size_t l_n = 0;
    try {
        while (cs.read_aw(fh_in, l_n)) {
            /* 
               if (!insert_all_dict) {
                // Add CSentence words to graph
                CSentenceSSI::iterator it = cs.begin();
                CSentenceSSI::iterator end = cs.end();
                for (; it != end; ++it) {
                    //@Aritza: Ya no es necesario pasarle el peso
                    kb.add_token(it->word()); //TODO: Corregir referencia sin definir
                }
            } 
             */

            sort_words(cs, sort_method);
            switch (dij_option) {
                case basic:
                    cs = ssi_dijkstra(cs, dsp_option);
                    break;

                case plus:
                    cs = ssi_dijkstra_plus(cs, dsp_option);
                    break;

                case fast:
                    cs = ssi_dijkstra_fast(cs, dsp_option);
                    break;
                default:
                    break;
            }
            if (pretty_mode) {
                sort_by_id(cs); // Just a little pretty thing
            }
            //@Aritza: Cambiamos el metodo por el cual se imprimen las palabras
            cs.print_csent_simple_all(cout, no_cnt, cnt_word);

            cs = CSentenceSSI(); // Next iteration ...		
        }
    } catch (std::exception & e) {
        cerr << "Error reading: " << fullname_in << " : " << e.what() << endl;
        throw (e);
    }
    return 0;
}
