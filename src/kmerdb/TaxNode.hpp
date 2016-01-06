#ifndef __TAXNODE__
#define __TAXNODE__

#include <string>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <stdint.h>
#include "metag_typedefs.hpp"

#define DEBUGME

using namespace std;

namespace metag {

/**
Encapsulates a taxonomy node; supplies methods for
finding parent, path to root, KmerNodes associated with
the tax node, etc.
*/

template<class tid_T>
class TaxTree;

template<class tid_T>
class TaxNode {
public:

  TaxNode() { }


  //! returns the taxonomy node id (either ncbi or kpath mfe)
  tid_T id() const { return m_id; }

  //! returns the id of the parent node;
  //! if no parent, this will be the same as
  //! the node id
  tid_T parent() const { return m_parent; }

  //! returns the node's name 
  const std::string & name() const { return m_name; }


  bool isLeaf() const {
    if (m_children.size() == 0) {
      return true;
    }
    return false;
  }

  void makeRoot() {
    m_parent = id();
  }

  //! returns true if 'id' is an ancestor of this
  //! TaxNode
  bool isAncestor(tid_T id) {
    if (m_set_to_root.find(id) == m_set_to_root.end()) {
      return false;
    }
    return true;
  }

  const std::set<tid_T>& getChildren() const { return m_children; } 

  void write(std::ostream &out = std::cout) {
    out << "id: " << m_id << " parent: " << m_parent << " name: " << m_name << " children: "; 
    for (typename std::set<tid_T>::const_iterator t = m_children.begin(); t != m_children.end(); t++) {
      std::cout <<*t << " ";
    }
    std::cout << std::endl;
  }

  uint16_t distToRoot() const {
    return m_path_to_root.size();
  }

  const std::vector<tid_T> * getPathToRoot() const {
    return &m_path_to_root;
  }

  int depth() const {
    return m_path_to_root.size()+1;
  }

  const string &getRank() const {
    return m_rank;
  }

  const string &getParentRank() const {
    return m_parent_rank;
  }

  void setRank(string &rank, string &parent_rank) {
    m_rank = rank;
    m_parent_rank = rank;
  }

  void printPathToRoot() const {
    std::cout << id() << " :: ";
    for (size_t j=0; j<m_path_to_root.size(); j++) {
      std::cout << m_path_to_root[j] << " ";
    }
    std::cout << endl;
  }

  const std::set<tid_T> * getSetToRoot() const {
    return &m_set_to_root;
  }

//  void setPathToRoot(TaxTree<tid_T> * tree);
void setPathToRoot(TaxTree<tid_T> *tree) {
  m_path_to_root.clear();
  m_set_to_root.clear();
  tree->getPathToRoot(id(), m_path_to_root);
  for (size_t j=0; j<m_path_to_root.size(); j++) {
    m_set_to_root.insert(m_path_to_root[j]);
  }
}

  void reversePathToRoot();

  void eraseChild(tid_T child) {
    m_children.erase(child);
  }

  static TaxNode * read(std::ifstream &in) {
    TaxNode *nd = new TaxNode();
    in >> nd->m_id;
    tid_T ct, id;
    in >> ct; //number of children

    for (tid_T j=0; j<ct; j++) {
      in >> id; //id of child
      if (id != nd->m_id) { //correct for possible error in kpath taxonomy
        nd->m_children.insert(id);
      }
    }
    in >> nd->m_parent;
    getline(in, nd->m_name);
    getline(in, nd->m_name);
    return nd;
  }

//friend template<class tid_T> class TaxTree;
private:

  tid_T m_id;

  tid_T m_parent;

  std::string m_name;

  std::string m_rank;
  std::string m_parent_rank;

  std::set<tid_T> m_children;

  //! IDs of node in the path to root
  std::vector<tid_T> m_path_to_root;
  std::set<tid_T> m_set_to_root;
};

}

#endif
