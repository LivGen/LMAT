#include "TaxNode.hpp"
#include "TaxTree.hpp"
#include <iostream>

using namespace metag;
using namespace std;

/*
void TaxNode::printPathToParent() const {
  for (size_t j=0; j<m_path_to_root.size(); j++) {
    cout << m_path_to_root[j] << " ";
  }
  cout << endl;
}
*/

/*
template<class tid_T>
void TaxNode<tid_T>::reversePathToRoot() {
  if (m_path_to_root.size() > 1) {
    reverse(m_path_to_root.begin(), m_path_to_root.end());
  }
}
*/


/*
template<class tid_T>
void TaxNode<tid_T>::setPathToRoot(TaxTree<tid_T> *tree) {
  m_path_to_root.clear();
  m_set_to_root.clear();
  tree->getPathToRoot(id(), m_path_to_root);
  for (size_t j=0; j<m_path_to_root.size(); j++) {
    m_set_to_root.insert(m_path_to_root[j]);
  }
}
*/
